from __future__ import annotations

import copy
import json
import math
import os
import selectors
import shlex
import signal
import subprocess
import sys
import time
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Protocol

from morse.construct import VaccineConstruct
from morse.router import RoutedSequence, ScorerInputRouter
from morse.scorer_scope import ScorerScopeConfig

if TYPE_CHECKING:
    from morse.config import ScorerConfig


@dataclass(frozen=True)
class ScoreResult:
    candidate_id: str
    outputs: dict[str, Any]


@dataclass(frozen=True)
class ConstructCandidate:
    candidate_id: str
    construct: VaccineConstruct


class SequenceScorer(Protocol):
    def score(self, sequences: list[str]) -> list[dict[str, Any]]: ...


@dataclass(frozen=True)
class ScorerBinding:
    scorer: SequenceScorer
    scope: ScorerScopeConfig
    router: ScorerInputRouter = field(default_factory=ScorerInputRouter)
    output_prefix: str | None = None
    namespace: str | None = None

    @property
    def name(self) -> str:
        return _display_name(self.scorer)

    def sequences(self, candidates: list[ConstructCandidate]) -> list[str]:
        return [item.sequence for item in self.routed_sequences(candidates)]

    def routed_sequences(
        self, candidates: list[ConstructCandidate]
    ) -> list[RoutedSequence]:
        return [
            self.router.route_with_context(candidate.construct, self.scope)
            for candidate in candidates
        ]


@dataclass(frozen=True)
class ComponentTiming:
    scorer: str
    batch_size: int
    elapsed_seconds: float
    candidate_ids: list[str]

    @property
    def seconds_per_candidate(self) -> float:
        return self.elapsed_seconds / self.batch_size if self.batch_size else 0.0


BASH_EXECUTABLE = "/bin/bash"


class BashCommandScorer:
    """Invoke an already prepared scorer once per JSON batch, then send EOF.

    Commands may be scripts, pipelines, curl, or docker exec. Only the local
    invocation's process group belongs to Core; external services and containers
    are never created, stopped, or restarted by this executor.
    """

    def __init__(self, command: str, timeout_seconds: float = 120,
                 name: str | None = None, env: dict[str, str] | None = None,
                 workdir: str | None = None,
                 parameters: dict[str, Any] | None = None) -> None:
        if not isinstance(command, str) or not command.strip():
            raise ValueError("Bash scorer command must be a non-empty string")
        if not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
            raise ValueError("Bash scorer timeout_seconds must be positive and finite")
        self.command = command
        self.timeout_seconds = timeout_seconds
        self.name = name or "Bash command"
        self.env = env or {}
        self.workdir = workdir
        from morse.config import parse_scorer_parameters
        self.parameters = parse_scorer_parameters({} if parameters is None else parameters)
        self.resolved_parameters: dict[str, Any] | None = None
        self._request_id = 0
        self._process: subprocess.Popen[str] | None = None

    def score(self, sequences: list[str]) -> list[dict[str, Any]]:
        return self.score_with_context(sequences, None)

    def score_with_context(self, sequences: list[str],
                           sequence_contexts: list[dict[str, object]] | None) -> list[dict[str, Any]]:
        if not sequences:
            return []
        if sequence_contexts is not None and len(sequence_contexts) != len(sequences):
            raise ValueError("sequence_contexts must contain one context per input sequence")
        self._request_id += 1
        request = {"request_id": self._request_id, "sequences": sequences}
        if sequence_contexts is not None:
            request["sequence_contexts"] = sequence_contexts
        if self.parameters:
            request["parameters"] = copy.deepcopy(self.parameters)
        environment = _merged_env(self.env)
        environment.setdefault("MORSE_PYTHON", sys.executable)
        try:
            process = subprocess.Popen(
                [BASH_EXECUTABLE, "-o", "pipefail", "-c", self.command],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, encoding="utf-8", errors="replace", env=environment, cwd=self.workdir,
                start_new_session=True,
            )
        except OSError as exc:
            raise RuntimeError(f"Bash scorer {self.name!r} could not start: {exc}") from exc
        self._process = process
        try:
            try:
                stdout, stderr = process.communicate(json.dumps(request, allow_nan=False) + "\n", timeout=self.timeout_seconds)
            except subprocess.TimeoutExpired as exc:
                self._kill_owned_group(process)
                _, stderr = process.communicate()
                detail = f"; stderr: {stderr.strip()[-4000:]}" if stderr.strip() else ""
                raise TimeoutError(f"Bash scorer {self.name!r} timed out after {self.timeout_seconds:g} s{detail}") from exc
            if stderr:
                sys.stderr.write(stderr)
                sys.stderr.flush()
            if process.returncode:
                detail = f"; stderr: {stderr.strip()[-4000:]}" if stderr.strip() else ""
                raise RuntimeError(f"Bash scorer {self.name!r} exited with code {process.returncode}{detail}")
            try:
                response = json.loads(stdout)
            except json.JSONDecodeError as exc:
                raise RuntimeError(f"Bash scorer {self.name!r} emitted invalid JSON; stdout must contain one JSON response (send diagnostics to stderr)") from exc
            if not isinstance(response, dict):
                raise RuntimeError(f"Bash scorer {self.name!r} must return one JSON response object")
            if type(response.get("request_id")) is not int or response.get("request_id") != self._request_id:
                raise RuntimeError(f"Bash scorer {self.name!r} returned request_id {response.get('request_id')!r}, expected {self._request_id}")
            if "error" in response:
                raise RuntimeError(f"Bash scorer {self.name!r} failed: {response['error']}")
            self.resolved_parameters = _response_parameters(response)
            try:
                return _outputs_from_batch_response(response, len(sequences))
            except RuntimeError as exc:
                raise RuntimeError(f"Bash scorer {self.name!r}: {exc}") from exc
        finally:
            self._process = None
            if process.poll() is None:
                self._kill_owned_group(process)
                process.wait(timeout=5)
            for pipe in (process.stdin, process.stdout, process.stderr):
                if pipe is not None:
                    pipe.close()

    @staticmethod
    def _kill_owned_group(process: subprocess.Popen[str]) -> None:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass

    def close(self) -> None:
        process = self._process
        if process is not None:
            self._kill_owned_group(process)


class PersistentProcessScorer:
    def __init__(
        self,
        command: list[str],
        timeout_seconds: int = 120,
        name: str | None = None,
        env: dict[str, str] | None = None,
        workdir: str | None = None,
        parameters: dict[str, Any] | None = None,
    ) -> None:
        self.command = command
        self.timeout_seconds = timeout_seconds
        self.name = name or " ".join(shlex.quote(part) for part in command)
        self.env = env or {}
        self.workdir = workdir
        from morse.config import parse_scorer_parameters
        self.parameters = parse_scorer_parameters({} if parameters is None else parameters)
        self.resolved_parameters: dict[str, Any] | None = None
        self._process: subprocess.Popen[str] | None = None
        self._request_id = 0

    def score(self, sequences: list[str]) -> list[dict[str, Any]]:
        return self.score_with_context(sequences, None)

    def score_with_context(
        self,
        sequences: list[str],
        sequence_contexts: list[dict[str, object]] | None,
        *,
        parameters: dict[str, Any] | None = None,
    ) -> list[dict[str, Any]]:
        if not sequences:
            return []
        if sequence_contexts is not None and len(sequence_contexts) != len(sequences):
            raise ValueError(
                "sequence_contexts must contain one context per input sequence"
            )
        if parameters is None:
            request_parameters = self.parameters
        else:
            from morse.config import parse_scorer_parameters
            request_parameters = parse_scorer_parameters(parameters)
        process = self._ensure_process()
        if process.stdin is None or process.stdout is None:
            raise RuntimeError(f"Persistent scorer {self.name!r} has closed pipes")

        self._request_id += 1
        request = {
            "request_id": self._request_id,
            "sequences": sequences,
            **(
                {"sequence_contexts": sequence_contexts}
                if sequence_contexts is not None
                else {}
            ),
        }
        if request_parameters:
            request["parameters"] = copy.deepcopy(request_parameters)
        process.stdin.write(json.dumps(request, allow_nan=False) + "\n")
        process.stdin.flush()

        line = self._readline_with_timeout(process)
        try:
            raw = json.loads(line)
        except json.JSONDecodeError as exc:
            raise RuntimeError(
                f"Persistent scorer {self.name!r} emitted invalid JSON: {line}"
            ) from exc
        if raw.get("request_id") != self._request_id:
            raise RuntimeError(
                f"Persistent scorer {self.name!r} returned request_id "
                f"{raw.get('request_id')}, expected {self._request_id}"
            )
        if "error" in raw:
            raise RuntimeError(f"Persistent scorer {self.name!r} failed: {raw['error']}")
        self.resolved_parameters = _response_parameters(raw)

        return _outputs_from_batch_response(raw, len(sequences))

    def close(self) -> None:
        process = self._process
        self._process = None
        if process is None:
            return
        if process.stdin is not None:
            try:
                process.stdin.close()
            except OSError:
                pass
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.terminate()
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=10)
        finally:
            for pipe in (process.stdout, process.stderr):
                if pipe is not None:
                    try:
                        pipe.close()
                    except OSError:
                        pass

    def _ensure_process(self) -> subprocess.Popen[str]:
        if self._process is not None and self._process.poll() is None:
            return self._process
        self._process = subprocess.Popen(
            self.command,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=None,
            text=True,
            bufsize=1,
            env=_merged_env(self.env),
            cwd=self.workdir,
        )
        return self._process

    def _readline_with_timeout(self, process: subprocess.Popen[str]) -> str:
        if process.stdout is None:
            raise RuntimeError(f"Persistent scorer {self.name!r} stdout is closed")
        selector = selectors.DefaultSelector()
        selector.register(process.stdout, selectors.EVENT_READ)
        events = selector.select(self.timeout_seconds)
        selector.close()
        if not events:
            self.close()
            raise TimeoutError(
                f"Persistent scorer {self.name!r} timed out after {self.timeout_seconds} s"
            )
        line = process.stdout.readline()
        if line == "":
            returncode = process.poll()
            raise RuntimeError(
                f"Persistent scorer {self.name!r} stopped before returning a response "
                f"(exit code {returncode})"
            )
        return line


class CompositeScorer:
    def __init__(self, bindings: list[ScorerBinding]) -> None:
        if not bindings:
            raise ValueError("CompositeScorer requires at least one scorer")
        self.bindings = bindings
        self.name = "+".join(binding.name for binding in bindings)
        self.last_timings: list[ComponentTiming] = []

    @property
    def resolved_parameters(self) -> dict[str, Any]:
        return {
            binding.namespace or binding.name: copy.deepcopy(binding.scorer.resolved_parameters)
            for binding in self.bindings
            if getattr(binding.scorer, "resolved_parameters", None) is not None
        }

    def score(self, candidates: list[ConstructCandidate]) -> list[ScoreResult]:
        self.last_timings = []
        merged: dict[str, dict[str, Any]] = {
            candidate.candidate_id: {} for candidate in candidates
        }
        for binding in self.bindings:
            routed = binding.routed_sequences(candidates)
            sequences = [item.sequence for item in routed]
            contexts = [item.context() for item in routed]
            started = time.perf_counter()
            score_with_context = getattr(binding.scorer, "score_with_context", None)
            batch_outputs = (
                score_with_context(sequences, contexts)
                if callable(score_with_context)
                else binding.scorer.score(sequences)
            )
            elapsed = time.perf_counter() - started
            self.last_timings.append(
                ComponentTiming(
                    scorer=binding.name,
                    batch_size=len(candidates),
                    elapsed_seconds=elapsed,
                    candidate_ids=[candidate.candidate_id for candidate in candidates],
                )
            )
            _validate_batch_outputs(batch_outputs, len(candidates), binding.name)
            for candidate, outputs in zip(
                candidates,
                batch_outputs,
                strict=True,
            ):
                bound_outputs = (
                    {f"{binding.namespace}.{name}": value for name, value in outputs.items()}
                    if binding.namespace is not None
                    else _bind_output_names(outputs, binding.output_prefix)
                )
                overlap = set(merged[candidate.candidate_id]) & set(bound_outputs)
                if overlap:
                    names = ", ".join(sorted(overlap))
                    raise RuntimeError(
                        f"Scorer {binding.name!r} returned duplicate output names: {names}"
                    )
                merged[candidate.candidate_id].update(bound_outputs)
        return [
            ScoreResult(candidate_id=candidate.candidate_id, outputs=merged[candidate.candidate_id])
            for candidate in candidates
        ]

    def close(self) -> None:
        for binding in self.bindings:
            close = getattr(binding.scorer, "close", None)
            if callable(close):
                close()


def _response_parameters(response: dict[str, Any]) -> dict[str, Any] | None:
    if "parameters" not in response:
        return None
    from morse.config import parse_scorer_parameters
    try:
        return parse_scorer_parameters(response["parameters"])
    except ValueError as exc:
        raise RuntimeError(f"Invalid resolved scorer parameters: {exc}") from exc


def docker_run_command(
    image: str,
    env: dict[str, str] | None = None,
    docker_args: list[str] | None = None,
) -> list[str]:
    command = ["docker", "run", "--rm", "-i"]
    command.extend(docker_args or [])
    for key, value in (env or {}).items():
        command.extend(["-e", f"{key}={value}"])
    command.append(image)
    return command


def make_raw_scorer(config: ScorerConfig, *, legacy_python: str | None = None) -> Any:
    """Create one runtime; sequence projection belongs to the caller's binding.

    ``legacy_python`` preserves Core CLI's historical interpreter substitution.
    Warm Starter leaves it unset to retain its legacy PATH-based argv behavior.
    Bash commands always use their own MORSE_PYTHON environment contract.
    """
    if config.mode == "bash-command":
        if not isinstance(config.command, str):
            raise ValueError(f"Bash scorer {config.name!r} must define a command string")
        return BashCommandScorer(
            config.command, timeout_seconds=config.timeout_seconds, name=config.name,
            env=config.env, workdir=config.workdir, parameters=config.parameters,
        )
    if config.mode == "persistent-container":
        if config.image is None:
            raise ValueError(f"Persistent container scorer {config.name!r} must define image")
        return PersistentProcessScorer(
            docker_run_command(config.image, env=config.env, docker_args=config.docker_args),
            timeout_seconds=config.timeout_seconds, name=config.name,
            parameters=config.parameters,
        )
    if config.mode == "persistent-process":
        if not config.command:
            raise ValueError(f"Persistent process scorer {config.name!r} must define command")
        command = list(config.command)
        if legacy_python is not None and command[0] in {"python", "python3", "@python"}:
            command[0] = legacy_python
        return PersistentProcessScorer(
            command, timeout_seconds=config.timeout_seconds, name=config.name,
            env=config.env, workdir=config.workdir, parameters=config.parameters,
        )
    raise ValueError(f"Unsupported scorer mode: {config.mode}")


def _outputs_from_batch_response(
    raw: dict[str, Any],
    expected_count: int,
) -> list[dict[str, Any]]:
    outputs = raw.get("outputs")
    if not isinstance(outputs, list):
        raise RuntimeError("Scorer response must define `outputs` as a list")
    if len(outputs) != expected_count:
        raise RuntimeError(
            f"Scorer returned {len(outputs)} outputs for {expected_count} sequences"
        )
    if not all(isinstance(item, dict) for item in outputs):
        raise RuntimeError("Every scorer output must be a JSON object")
    return [dict(values) for values in outputs]


def _validate_batch_outputs(
    outputs: list[dict[str, Any]],
    expected_count: int,
    scorer_name: str,
) -> None:
    if len(outputs) != expected_count:
        raise RuntimeError(
            f"Scorer {scorer_name!r} returned {len(outputs)} outputs for "
            f"{expected_count} sequences"
        )
    if not all(isinstance(item, dict) for item in outputs):
        raise RuntimeError(
            f"Scorer {scorer_name!r} must return one output mapping per sequence"
        )
    for row in outputs:
        for name, value in row.items():
            if isinstance(value, float) and not math.isfinite(value):
                raise RuntimeError(f"Scorer {scorer_name!r} returned non-finite metric {name!r}")


def _bind_output_names(
    outputs: dict[str, Any],
    output_prefix: str | None,
) -> dict[str, Any]:
    if output_prefix is None:
        return outputs
    return {f"{output_prefix}_{name}": value for name, value in outputs.items()}


def _merged_env(extra: dict[str, str]) -> dict[str, str]:
    env = os.environ.copy()
    env.update(extra)
    return env


def _display_name(scorer: Any) -> str:
    name = getattr(scorer, "name", None)
    if name:
        return str(name)
    image = getattr(scorer, "image", None)
    if image:
        return f"{scorer.__class__.__name__}:{image}"
    command = getattr(scorer, "command", None)
    if command:
        return " ".join(shlex.quote(part) for part in command)
    return scorer.__class__.__name__
