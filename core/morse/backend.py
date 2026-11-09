from __future__ import annotations

import argparse
import hashlib
import hmac
import json
import os
import queue
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid
import zipfile
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime, timezone
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, unquote, urlsplit

import yaml

from morse import __version__
from morse.input_bundle import stage_input_bundle
from morse.run_spec import RunSpec, load_run_spec, write_run_spec
from morse.visualization import render_sequence_scope_svg

TERMINAL_STATUSES = {"succeeded", "failed", "canceled", "interrupted"}
RUN_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
MAX_REQUEST_BYTES = 12 * 1024 * 1024
BACKEND_API_VERSION = 6
SCORER_PROTOCOL_VERSION = 3


def runtime_metadata() -> dict[str, Any]:
    """Return deployment identity used by notebook compatibility checks."""

    return {
        "core_version": __version__,
        "core_revision": _core_revision(),
        "backend_api_version": BACKEND_API_VERSION,
        "scorer_protocol_version": SCORER_PROTOCOL_VERSION,
    }


def _core_revision() -> str:
    configured = os.environ.get("MORSE_CORE_REVISION", "").strip()
    if configured:
        return configured
    project_root = Path(__file__).resolve().parent.parent
    revision_file = project_root / ".morse-revision"
    if revision_file.is_file():
        revision = revision_file.read_text(errors="replace").strip()
        if revision:
            return revision
    completed = subprocess.run(
        ["git", "-C", str(project_root), "rev-parse", "--short=12", "HEAD"],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode == 0 and completed.stdout.strip():
        return completed.stdout.strip()
    return "unknown"


class BackendError(Exception):
    def __init__(self, status: int, message: str) -> None:
        super().__init__(message)
        self.status = status
        self.message = message


@dataclass(frozen=True)
class BackendSettings:
    config_root: Path
    results_root: Path
    workdir: Path
    token: str

    @classmethod
    def from_values(
        cls,
        *,
        config_root: Path,
        results_root: Path,
        workdir: Path | None = None,
        token: str | None = None,
    ) -> BackendSettings:
        resolved_config_root = config_root.expanduser().resolve()
        resolved_results_root = results_root.expanduser().resolve()
        resolved_workdir = (
            workdir.expanduser().resolve()
            if workdir is not None
            else resolved_config_root.parent
        )
        resolved_token = (
            os.environ.get("MORSE_BACKEND_TOKEN", "") if token is None else token
        )
        if not resolved_token:
            raise ValueError(
                "MORSE_BACKEND_TOKEN is required; generate one with "
                "`python -c \"import secrets; print(secrets.token_urlsafe(32))\"`"
            )
        if not resolved_config_root.is_dir():
            raise ValueError(
                f"Backend config root does not exist: {resolved_config_root}"
            )
        if not resolved_workdir.is_dir():
            raise ValueError(f"Backend workdir does not exist: {resolved_workdir}")
        resolved_results_root.mkdir(parents=True, exist_ok=True)
        return cls(
            config_root=resolved_config_root,
            results_root=resolved_results_root,
            workdir=resolved_workdir,
            token=resolved_token,
        )


@dataclass(frozen=True)
class PreparedSubmission:
    config_name: str
    config_path: Path
    custom_configuration: bool
    run_name: str
    validated: Any


class JobManager:
    def __init__(
        self,
        settings: BackendSettings,
        *,
        command_builder: Callable[[dict[str, Any]], list[str]] | None = None,
        config_validator: Callable[[Path], Any] = load_run_spec,
    ) -> None:
        self.settings = settings
        self._command_builder = command_builder or self._default_command
        self._config_validator = config_validator
        self._lock = threading.RLock()
        self._archive_lock = threading.Lock()
        self._jobs: dict[str, dict[str, Any]] = {}
        self._processes: dict[str, subprocess.Popen[bytes]] = {}
        self._queue: queue.Queue[str | None] = queue.Queue()
        self._stopping = threading.Event()
        self._restore_jobs()
        self._worker = threading.Thread(
            target=self._worker_loop,
            name="morse-backend-worker",
            daemon=True,
        )
        self._worker.start()

    def list_configs(self) -> list[str]:
        configs: list[str] = []
        paths = list(self.settings.config_root.rglob("*.yaml"))
        paths.extend(self.settings.config_root.rglob("*.yml"))
        for path in sorted(paths):
            if path.is_file() and path.resolve().is_relative_to(
                self.settings.config_root
            ):
                configs.append(path.relative_to(self.settings.config_root).as_posix())
        return configs

    def list_jobs(self) -> list[dict[str, Any]]:
        with self._lock:
            jobs = [self._public_record(record) for record in self._jobs.values()]
        return sorted(jobs, key=lambda item: item["created_at"], reverse=True)

    def submit(
        self,
        payload: dict[str, Any],
        *,
        kind: str = "optimization",
    ) -> dict[str, Any]:
        if kind not in {"optimization", "pilot"}:
            raise BackendError(HTTPStatus.BAD_REQUEST, "Unknown job kind")
        pilot_samples = self._parse_pilot_samples(
            payload.get("pilot_samples"),
            kind=kind,
        )
        job_id = uuid.uuid4().hex
        job_dir = self.settings.results_root / job_id
        prepared = self._prepare_submission(
            payload,
            job_dir=job_dir,
        )
        initialization = (
            prepared.validated.run_config.initialization.to_dict()
            if isinstance(prepared.validated, RunSpec)
            else {"mode": "random"}
        )
        output_dir = job_dir / "results"
        config_fingerprint = hashlib.sha256(
            prepared.config_path.read_bytes()
        ).hexdigest()
        record: dict[str, Any] = {
            "id": job_id,
            "kind": kind,
            "status": "queued",
            "config": prepared.config_name,
            "config_path": str(prepared.config_path),
            "custom_configuration": prepared.custom_configuration,
            "run_name": prepared.run_name,
            "configuration_fingerprint": config_fingerprint,
            "initialization": initialization,
            "resume_from": None,
            "pilot_samples": pilot_samples,
            "job_dir": str(job_dir),
            "output_dir": str(output_dir),
            "log_path": str(job_dir / "run.log"),
            "created_at": _utc_now(),
            "started_at": None,
            "finished_at": None,
            "exit_code": None,
            "error": None,
            "timings": {},
        }
        with self._lock:
            self._jobs[job_id] = record
            self._write_record(record)
        self._queue.put(job_id)
        return self._public_record(record)

    def submit_pilot(self, payload: dict[str, Any]) -> dict[str, Any]:
        return self.submit(payload, kind="pilot")

    def preflight(self, payload: dict[str, Any]) -> dict[str, Any]:
        """Validate a complete configuration without starting or queuing a run."""
        request_id = f"preflight-{uuid.uuid4().hex}"
        preflight_root = self.settings.results_root / ".preflight"
        preflight_root.mkdir(exist_ok=True)
        job_dir = preflight_root / request_id
        try:
            prepared = self._prepare_submission(payload, job_dir=job_dir)
            spec = prepared.validated
            if not isinstance(spec, RunSpec):
                raise BackendError(HTTPStatus.BAD_REQUEST, "Preflight requires a valid run configuration")
            outputs, optimizer = spec.run_config.outputs, spec.run_config.optimizer
            return {
                "valid": True, "runtime": runtime_metadata(),
                "config": prepared.config_name, "run_name": prepared.run_name,
                "algorithm": optimizer.type, "population_size": optimizer.population_size,
                "generations": optimizer.generations, "seed": optimizer.seed,
                "objective_count": sum(output.role == "objective" for output in outputs),
                "threshold_count": sum(output.role == "threshold" for output in outputs),
                "scorers": [{key: value for key, value in scorer.items()
                             if key in {"name", "command", "scope", "metrics"}}
                            for scorer in spec.to_dict()["scorers"]],
                "configuration_fingerprint": hashlib.sha256(prepared.config_path.read_bytes()).hexdigest(),
                "sequence_scope_svg": render_sequence_scope_svg(spec),
            }
        finally:
            shutil.rmtree(job_dir, ignore_errors=True)

    def _prepare_submission(self, payload: dict[str, Any], *, job_dir: Path) -> PreparedSubmission:
        # A bearer token is a trusted execution credential: the submitted YAML
        # contains Bash commands for an already prepared environment.
        unknown = sorted(set(payload) - {"config", "configuration", "input_files", "pilot_samples"})
        if unknown:
            raise BackendError(HTTPStatus.BAD_REQUEST, "Unknown job fields: " + ", ".join(unknown))
        inline = payload.get("configuration")
        config_name = str(payload.get("config", "")).strip()
        if inline is not None:
            if config_name:
                raise BackendError(HTTPStatus.BAD_REQUEST, "Provide either config or configuration, not both")
            if not isinstance(inline, dict):
                raise BackendError(HTTPStatus.BAD_REQUEST, "configuration must be a JSON object")
            config_name = "notebook-generated"
        else:
            if "input_files" in payload:
                raise BackendError(HTTPStatus.BAD_REQUEST, "input_files requires configuration")
            source = self._resolve_config(config_name)
        job_dir.mkdir(parents=True)
        config_path = job_dir / "submitted_config.yaml"
        try:
            if inline is not None:
                raw = json.loads(json.dumps(inline))
                raw["name"] = raw.get("name") or "notebook-run"
                if "input_files" in payload:
                    stage_input_bundle(raw, payload["input_files"], self.settings.results_root / ".inputs")
                if isinstance(raw.get("output_dir"), str):
                    requested_output = Path(raw["output_dir"]).expanduser()
                    if not requested_output.is_absolute():
                        raw["output_dir"] = str((self.settings.config_root / requested_output).resolve())
                config_path.write_text(yaml.safe_dump(raw, sort_keys=False, allow_unicode=True))
                validated = (load_run_spec(config_path, input_root=self.settings.config_root)
                             if self._config_validator is load_run_spec else self._config_validator(config_path))
            else:
                validated = self._config_validator(source)
                shutil.copyfile(source, config_path)
            if isinstance(validated, RunSpec):
                run_name = validated.name
                # Freeze both inline and selected configurations. Resolving paths
                # before the copy preserves the source directory's input meaning.
                write_run_spec(config_path, validated)
            else:
                run_name = str(inline.get("name") or "notebook-run") if inline is not None else source.stem
            if not RUN_NAME_PATTERN.fullmatch(run_name):
                raise ValueError("config name must contain only letters, numbers, dot, underscore, or hyphen and be at most 128 characters")
        except Exception as exc:
            shutil.rmtree(job_dir, ignore_errors=True)
            raise BackendError(HTTPStatus.BAD_REQUEST, f"Config validation failed for {config_name!r}: {exc}") from exc
        return PreparedSubmission(config_name=config_name, config_path=config_path,
                                  custom_configuration=inline is not None, run_name=run_name, validated=validated)

    @staticmethod
    def _parse_pilot_samples(value: Any, *, kind: str) -> int | None:
        if kind != "pilot":
            if value is not None:
                raise BackendError(
                    HTTPStatus.BAD_REQUEST,
                    "pilot_samples is accepted only for pilot jobs",
                )
            return None
        try:
            samples = 24 if value is None else int(value)
        except (TypeError, ValueError) as exc:
            raise BackendError(
                HTTPStatus.BAD_REQUEST, "pilot_samples must be an integer"
            ) from exc
        if not 2 <= samples <= 512:
            raise BackendError(
                HTTPStatus.BAD_REQUEST,
                "pilot_samples must be between 2 and 512",
            )
        return samples

    def get(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            record = self._jobs.get(job_id)
            if record is None:
                raise BackendError(HTTPStatus.NOT_FOUND, "Job not found")
            return self._public_record(record)

    def log_tail(self, job_id: str, lines: int = 200) -> str:
        with self._lock:
            record = self._jobs.get(job_id)
            if record is None:
                raise BackendError(HTTPStatus.NOT_FOUND, "Job not found")
            log_path = Path(record["log_path"])
        return _tail_text(log_path, max(1, min(lines, 2_000)))

    def metrics(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            record = self._jobs.get(job_id)
            if record is None:
                raise BackendError(HTTPStatus.NOT_FOUND, "Job not found")
            status = str(record["status"])
            progress_path = Path(record["output_dir"]) / "progress.json"
            try:
                total_generations = load_run_spec(
                    Path(record["config_path"])
                ).run_config.optimizer.generations
            except Exception:
                # Old records may point at removed configs. Progress data, when
                # present, remains the authoritative fallback.
                total_generations = record.get("overrides", {}).get("generations")

        payload: dict[str, Any] = {
            "schema_version": 1,
            "generation": 0,
            "total_generations": total_generations,
            "progress_fraction": 0.0,
            "hypervolume": None,
            "best_hypervolume": None,
            "history": [],
            "objectives": [],
            "thresholds": [],
            "scorer_timings": [],
        }
        if progress_path.is_file():
            try:
                parsed = json.loads(progress_path.read_text())
            except (OSError, json.JSONDecodeError) as exc:
                raise BackendError(
                    HTTPStatus.CONFLICT,
                    "Optimizer metrics are being updated; retry shortly",
                ) from exc
            if not isinstance(parsed, dict):
                raise BackendError(
                    HTTPStatus.INTERNAL_SERVER_ERROR,
                    "Optimizer metrics must contain a JSON object",
                )
            payload.update(parsed)
        payload["job_id"] = job_id
        payload["status"] = status
        return payload

    def pilot_result(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            record = self._jobs.get(job_id)
            if record is None:
                raise BackendError(HTTPStatus.NOT_FOUND, "Job not found")
            if record.get("kind") != "pilot":
                raise BackendError(HTTPStatus.CONFLICT, "Job is not a pilot")
            status = str(record["status"])
            path = Path(record["output_dir"]) / "pilot.json"
        if status not in TERMINAL_STATUSES:
            raise BackendError(HTTPStatus.CONFLICT, "Pilot is still running")
        if not path.is_file():
            raise BackendError(HTTPStatus.NOT_FOUND, "Pilot result was not produced")
        try:
            result = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise BackendError(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                "Pilot result is unreadable",
            ) from exc
        if not isinstance(result, dict):
            raise BackendError(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                "Pilot result must be a JSON object",
            )
        return result

    def cancel(self, job_id: str) -> dict[str, Any]:
        process: subprocess.Popen[bytes] | None = None
        with self._lock:
            record = self._jobs.get(job_id)
            if record is None:
                raise BackendError(HTTPStatus.NOT_FOUND, "Job not found")
            if record["status"] in TERMINAL_STATUSES:
                return self._public_record(record)
            if record["status"] == "queued":
                record["status"] = "canceled"
                record["finished_at"] = _utc_now()
                self._write_record(record)
                return self._public_record(record)
            record["status"] = "canceling"
            self._write_record(record)
            process = self._processes.get(job_id)

        if process is not None and process.poll() is None:
            _terminate_process(process)
        return self.get(job_id)

    def resume(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            record = self._jobs.get(job_id)
            if record is None:
                raise BackendError(HTTPStatus.NOT_FOUND, "Job not found")
            if record.get("kind") != "optimization":
                raise BackendError(HTTPStatus.CONFLICT, "Only optimization jobs can resume")
            if record.get("status") not in {"interrupted", "failed", "canceled"}:
                raise BackendError(
                    HTTPStatus.CONFLICT,
                    "Only interrupted, failed, or canceled jobs can resume",
                )
            output_dir = Path(record["output_dir"])
            checkpoint = output_dir / "optimizer-checkpoint.pkl"
            if not checkpoint.is_file():
                raise BackendError(
                    HTTPStatus.CONFLICT,
                    "The job has no completed-generation optimizer checkpoint",
                )
            record["status"] = "queued"
            record["resume_from"] = str(output_dir)
            record["started_at"] = None
            record["finished_at"] = None
            record["exit_code"] = None
            record["error"] = None
            archive_path = Path(record["job_dir"]) / f"{job_id}-artifacts.zip"
            archive_path.unlink(missing_ok=True)
            self._write_record(record)
        self._queue.put(job_id)
        return self.get(job_id)

    def artifact_path(self, job_id: str) -> Path:
        with self._lock:
            record = self._jobs.get(job_id)
            if record is None:
                raise BackendError(HTTPStatus.NOT_FOUND, "Job not found")
            if record["status"] not in TERMINAL_STATUSES:
                raise BackendError(
                    HTTPStatus.CONFLICT,
                    "Artifacts are available after the job reaches a terminal status",
                )
            job_dir = Path(record["job_dir"])
            archive_path = job_dir / f"{job_id}-artifacts.zip"
            if archive_path.exists():
                return archive_path

        with self._archive_lock:
            if archive_path.exists():
                return archive_path
            temporary_path = archive_path.with_suffix(".zip.tmp")
            with zipfile.ZipFile(
                temporary_path, "w", compression=zipfile.ZIP_DEFLATED
            ) as archive:
                for path in sorted(job_dir.rglob("*")):
                    if not path.is_file() or path in {
                        archive_path,
                        temporary_path,
                    }:
                        continue
                    archive.write(path, path.relative_to(job_dir))
            temporary_path.replace(archive_path)
            return archive_path

    def sequence_scope_path(self, job_id: str) -> Path:
        with self._lock:
            record = self._jobs.get(job_id)
            if record is None:
                raise BackendError(HTTPStatus.NOT_FOUND, "Job not found")
            status = str(record["status"])
            scope_path = Path(record["output_dir"]) / "sequence_scope.svg"
        if scope_path.is_file():
            return scope_path
        if status in TERMINAL_STATUSES:
            raise BackendError(
                HTTPStatus.NOT_FOUND,
                "The job finished without generating sequence_scope.svg",
            )
        raise BackendError(
            HTTPStatus.CONFLICT,
            "Sequence scope is still being generated",
        )

    def summary(self) -> dict[str, int]:
        with self._lock:
            statuses = [record["status"] for record in self._jobs.values()]
        return {
            status: statuses.count(status)
            for status in sorted(set(statuses) | {"queued", "running"})
        }

    def shutdown(self) -> None:
        if self._stopping.is_set():
            return
        self._stopping.set()
        with self._lock:
            processes = list(self._processes.items())
            for job_id, _process in processes:
                record = self._jobs[job_id]
                record["status"] = "interrupting"
                self._write_record(record)
        for _job_id, process in processes:
            if process.poll() is None:
                _terminate_process(process)
        self._queue.put(None)
        self._worker.join(timeout=10)

    def _resolve_config(self, config_name: str) -> Path:
        if not config_name:
            raise BackendError(HTTPStatus.BAD_REQUEST, "config is required")
        relative = Path(config_name)
        if relative.is_absolute():
            raise BackendError(
                HTTPStatus.BAD_REQUEST, "config must be relative to the config root"
            )
        candidate = (self.settings.config_root / relative).resolve()
        if not candidate.is_relative_to(self.settings.config_root):
            raise BackendError(
                HTTPStatus.BAD_REQUEST, "config escapes the config root"
            )
        if candidate.suffix not in {".yaml", ".yml"} or not candidate.is_file():
            raise BackendError(HTTPStatus.NOT_FOUND, "Config not found")
        return candidate

    def _default_command(self, record: dict[str, Any]) -> list[str]:
        if record.get("kind") == "pilot":
            return [
                sys.executable,
                "-u",
                "-m",
                "morse",
                "pilot",
                "--config",
                record["config_path"],
                "--samples",
                str(record.get("pilot_samples") or 24),
                "--output-dir",
                record["output_dir"],
            ]
        command = [
            sys.executable,
            "-u",
            "-m",
            "morse",
            "run",
            "--config",
            record["config_path"],
            "--output-dir",
            record["output_dir"],
        ]
        if record.get("resume_from"):
            command.extend(["--resume-from", str(record["resume_from"])])
        return command

    def _worker_loop(self) -> None:
        while not self._stopping.is_set():
            job_id = self._queue.get()
            if job_id is None:
                return
            with self._lock:
                record = self._jobs.get(job_id)
                if record is None or record["status"] != "queued":
                    continue
            try:
                self._run_job(job_id)
            except Exception as exc:
                with self._lock:
                    record = self._jobs[job_id]
                    record["status"] = "failed"
                    record["finished_at"] = _utc_now()
                    record["error"] = (
                        f"Backend worker failed: {exc.__class__.__name__}: {exc}"
                    )
                    self._write_record(record)

    def _run_job(self, job_id: str) -> None:
        total_started = time.perf_counter()
        with self._lock:
            record = self._jobs[job_id]
            record["status"] = "running"
            record["started_at"] = _utc_now()
            self._write_record(record)
            log_path = Path(record["log_path"])

        exit_code: int | None = None
        launch_error: str | None = None
        with log_path.open("ab") as log:
            log.write(
                (
                    f"[morse-backend] started {record['started_at']} "
                    f"config={record['config']} run_name={record['run_name']}\n"
                ).encode()
            )
            log.flush()
            try:
                command = self._command_builder(record)
                process_started = time.perf_counter()
                process = subprocess.Popen(
                    command,
                    cwd=self.settings.workdir,
                    env={**os.environ, "PYTHONUNBUFFERED": "1"},
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    start_new_session=os.name != "nt",
                )
                with self._lock:
                    self._processes[job_id] = process
                    should_cancel = self._jobs[job_id]["status"] == "canceling"
                if should_cancel:
                    _terminate_process(process)
                exit_code = process.wait()
                record.setdefault("timings", {})[
                    "optimizer_process_seconds"
                    if record.get("kind") == "optimization"
                    else "pilot_process_seconds"
                ] = time.perf_counter() - process_started
            except Exception as exc:
                launch_error = f"{exc.__class__.__name__}: {exc}"
                log.write(f"[morse-backend] launch failed: {launch_error}\n".encode())
            finally:
                with self._lock:
                    self._processes.pop(job_id, None)
            try:
                self._write_provenance(record)
            except Exception as exc:
                if launch_error is None:
                    launch_error = (
                        "Failed to write provenance: "
                        f"{exc.__class__.__name__}: {exc}"
                    )
            record.setdefault("timings", {})["total_seconds"] = (
                time.perf_counter() - total_started
            )
            summary_path = Path(record["output_dir"]) / "summary.json"
            if summary_path.is_file():
                try:
                    summary = json.loads(summary_path.read_text())
                    first_calls = summary.get("scorer_timing", {}).get(
                        "startup_inclusive_first_calls_seconds"
                    )
                    if first_calls is not None:
                        record["timings"][
                            "scorer_startup_inclusive_seconds"
                        ] = float(first_calls)
                except (OSError, ValueError, TypeError):
                    pass

        with self._lock:
            record = self._jobs[job_id]
            if record["status"] == "canceling":
                record["status"] = "canceled"
            elif record["status"] == "interrupting":
                record["status"] = "interrupted"
                record["error"] = "Backend stopped before the job completed"
            elif launch_error is not None:
                record["status"] = "failed"
                record["error"] = launch_error
            elif exit_code == 0:
                record["status"] = "succeeded"
            else:
                record["status"] = "failed"
                record["error"] = _optimizer_exit_error(log_path, exit_code)
            record["exit_code"] = exit_code
            record["finished_at"] = _utc_now()
            self._write_record(record)

    def _write_provenance(self, record: dict[str, Any]) -> None:
        output_dir = Path(record["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)
        config_path = Path(record["config_path"])
        try:
            spec = load_run_spec(config_path)
        except (OSError, ValueError, TypeError):
            configured_scorers = []
        else:
            configured_scorers = [
                {
                    "name": scorer.name,
                    "mode": scorer.mode,
                    "image": scorer.image,
                    "command": scorer.command,
                    "output_prefix": scorer.output_prefix,
                    "namespace": scorer.namespace,
                    "metrics": spec.run_config.to_dict()["scorers"][index]["metrics"],
                }
                for index, scorer in enumerate(spec.run_config.scorers)
            ]
        provenance = {
            "schema_version": 1,
            "job_id": record["id"],
            "job_kind": record.get("kind", "optimization"),
            "created_at": record.get("created_at"),
            "runtime": runtime_metadata(),
            "configuration_fingerprint": record.get("configuration_fingerprint"),
            "resolved_config_sha256": hashlib.sha256(
                config_path.read_bytes()
            ).hexdigest(),
            "initialization": record.get("initialization"),
            "scorers": configured_scorers,
        }
        (output_dir / "provenance.json").write_text(
            json.dumps(provenance, indent=2, sort_keys=True) + "\n"
        )

    def _restore_jobs(self) -> None:
        for path in sorted(self.settings.results_root.glob("*/job.json")):
            try:
                record = json.loads(path.read_text())
                job_id = str(record["id"])
            except (OSError, ValueError, KeyError, TypeError):
                continue
            if not re.fullmatch(r"[0-9a-f]{32}", job_id):
                continue
            job_dir = path.parent.resolve()
            record["job_dir"] = str(job_dir)
            record["output_dir"] = str(job_dir / "results")
            record["log_path"] = str(job_dir / "run.log")
            record.setdefault("kind", "optimization")
            record.setdefault(
                "custom_configuration",
                bool(record.get("custom_design", False)),
            )
            record.setdefault("initialization", {"mode": "random"})
            record.setdefault("overrides", {})
            record.setdefault("resume_from", None)
            record.setdefault("pilot_samples", None)
            record.setdefault("timings", {})
            if record.get("status") not in TERMINAL_STATUSES:
                record["status"] = "interrupted"
                record["finished_at"] = _utc_now()
                record["error"] = "Backend restarted before the job completed"
                self._write_record(record)
            self._jobs[job_id] = record

    def _write_record(self, record: dict[str, Any]) -> None:
        path = Path(record["job_dir"]) / "job.json"
        temporary_path = path.with_suffix(".json.tmp")
        temporary_path.write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n"
        )
        temporary_path.replace(path)

    @staticmethod
    def _public_record(record: dict[str, Any]) -> dict[str, Any]:
        return {
            key: record.get(key)
            for key in (
                "id",
                "kind",
                "status",
                "config",
                "custom_configuration",
                "run_name",
                "configuration_fingerprint",
                "initialization",
                "resume_from",
                "pilot_samples",
                "created_at",
                "started_at",
                "finished_at",
                "exit_code",
                "error",
                "timings",
            )
        }


def _optimizer_exit_error(log_path: Path, exit_code: int | None) -> str:
    fallback = f"Optimizer exited with code {exit_code}"
    for line in reversed(_tail_text(log_path, 200).splitlines()):
        stripped = line.strip()
        if stripped.lower().startswith("error:"):
            detail = stripped.split(":", 1)[1].strip()
            if detail:
                return f"{fallback}: {detail}"
    return fallback


class MorseHTTPServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(
        self,
        server_address: tuple[str, int],
        manager: JobManager,
        token: str,
    ) -> None:
        super().__init__(server_address, MorseRequestHandler)
        self.manager = manager
        self.token = token


class MorseRequestHandler(BaseHTTPRequestHandler):
    server: MorseHTTPServer

    def do_GET(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler API
        self._handle(self._get)

    def do_POST(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler API
        self._handle(self._post)

    def do_DELETE(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler API
        self._handle(self._delete)

    def _handle(self, action: Callable[[], None]) -> None:
        try:
            action()
        except BackendError as exc:
            self._send_json(exc.status, {"error": exc.message})
        except Exception as exc:
            self.log_error("Unhandled backend error: %s", exc)
            self._send_json(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                {"error": "Internal backend error"},
            )

    def _get(self) -> None:
        parsed = urlsplit(self.path)
        if parsed.path == "/health":
            self._send_json(
                HTTPStatus.OK,
                {
                    "status": "ok",
                    "runtime": runtime_metadata(),
                    "jobs": self.server.manager.summary(),
                },
            )
            return
        self._require_auth()
        if parsed.path == "/v1/configs":
            self._send_json(
                HTTPStatus.OK, {"configs": self.server.manager.list_configs()}
            )
            return
        if parsed.path == "/v1/jobs":
            self._send_json(
                HTTPStatus.OK, {"jobs": self.server.manager.list_jobs()}
            )
            return
        job_id, suffix = _job_route(parsed.path)
        if suffix == "":
            self._send_json(HTTPStatus.OK, self.server.manager.get(job_id))
            return
        if suffix == "/log":
            query = parse_qs(parsed.query)
            try:
                lines = int(query.get("tail", ["200"])[0])
            except ValueError as exc:
                raise BackendError(
                    HTTPStatus.BAD_REQUEST, "tail must be an integer"
                ) from exc
            self._send_text(
                HTTPStatus.OK, self.server.manager.log_tail(job_id, lines)
            )
            return
        if suffix == "/metrics":
            self._send_json(
                HTTPStatus.OK,
                self.server.manager.metrics(job_id),
            )
            return
        if suffix == "/pilot":
            self._send_json(
                HTTPStatus.OK,
                self.server.manager.pilot_result(job_id),
            )
            return
        if suffix == "/artifacts":
            self._send_file(self.server.manager.artifact_path(job_id))
            return
        if suffix == "/sequence-scope":
            self._send_file(
                self.server.manager.sequence_scope_path(job_id),
                content_type="image/svg+xml; charset=utf-8",
                disposition="inline",
            )
            return
        raise BackendError(HTTPStatus.NOT_FOUND, "Endpoint not found")

    def _post(self) -> None:
        self._require_auth()
        path = urlsplit(self.path).path
        payload = self._read_json()
        if path == "/v1/jobs/preflight":
            self._send_json(
                HTTPStatus.OK,
                self.server.manager.preflight(payload),
            )
            return
        if path == "/v1/jobs/pilot":
            self._send_json(
                HTTPStatus.ACCEPTED,
                self.server.manager.submit_pilot(payload),
            )
            return
        if path.startswith("/v1/jobs/"):
            job_id, suffix = _job_route(path)
            if suffix != "/resume":
                raise BackendError(HTTPStatus.NOT_FOUND, "Endpoint not found")
            if payload:
                raise BackendError(
                    HTTPStatus.BAD_REQUEST,
                    "The resume endpoint does not accept a request body",
                )
            self._send_json(
                HTTPStatus.ACCEPTED,
                self.server.manager.resume(job_id),
            )
            return
        if path != "/v1/jobs":
            raise BackendError(HTTPStatus.NOT_FOUND, "Endpoint not found")
        self._send_json(
            HTTPStatus.ACCEPTED, self.server.manager.submit(payload)
        )

    def _delete(self) -> None:
        self._require_auth()
        job_id, suffix = _job_route(urlsplit(self.path).path)
        if suffix:
            raise BackendError(HTTPStatus.NOT_FOUND, "Endpoint not found")
        self._send_json(
            HTTPStatus.ACCEPTED, self.server.manager.cancel(job_id)
        )

    def _require_auth(self) -> None:
        authorization = self.headers.get("Authorization", "")
        prefix = "Bearer "
        supplied = (
            authorization[len(prefix) :] if authorization.startswith(prefix) else ""
        )
        if not supplied or not hmac.compare_digest(supplied, self.server.token):
            raise BackendError(HTTPStatus.UNAUTHORIZED, "Invalid bearer token")

    def _read_json(self) -> dict[str, Any]:
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError as exc:
            raise BackendError(
                HTTPStatus.BAD_REQUEST, "Invalid Content-Length"
            ) from exc
        if length < 1 or length > MAX_REQUEST_BYTES:
            raise BackendError(
                HTTPStatus.BAD_REQUEST,
                f"JSON body must be between 1 and {MAX_REQUEST_BYTES} bytes",
            )
        try:
            payload = json.loads(self.rfile.read(length))
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            raise BackendError(HTTPStatus.BAD_REQUEST, "Invalid JSON body") from exc
        if not isinstance(payload, dict):
            raise BackendError(HTTPStatus.BAD_REQUEST, "JSON body must be an object")
        return payload

    def _send_json(self, status: int, payload: Any) -> None:
        body = json.dumps(payload, separators=(",", ":")).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(body)

    def _send_text(self, status: int, text: str) -> None:
        body = text.encode()
        self.send_response(status)
        self.send_header("Content-Type", "text/plain; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(body)

    def _send_file(
        self,
        path: Path,
        *,
        content_type: str = "application/zip",
        disposition: str = "attachment",
    ) -> None:
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(path.stat().st_size))
        self.send_header(
            "Content-Disposition", f'{disposition}; filename="{path.name}"'
        )
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        with path.open("rb") as handle:
            while chunk := handle.read(1024 * 1024):
                self.wfile.write(chunk)

    def log_message(self, format: str, *args: Any) -> None:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        sys.stderr.write(
            f"{timestamp} {self.client_address[0]} {format % args}\n"
        )


def _job_route(path: str) -> tuple[str, str]:
    prefix = "/v1/jobs/"
    if not path.startswith(prefix):
        raise BackendError(HTTPStatus.NOT_FOUND, "Endpoint not found")
    remainder = path[len(prefix) :]
    raw_job_id, separator, raw_suffix = remainder.partition("/")
    job_id = unquote(raw_job_id)
    if not re.fullmatch(r"[0-9a-f]{32}", job_id):
        raise BackendError(HTTPStatus.NOT_FOUND, "Job not found")
    suffix = f"/{raw_suffix}" if separator else ""
    return job_id, suffix


def _tail_text(path: Path, lines: int) -> str:
    if not path.exists():
        return ""
    block_size = 8192
    chunks: list[bytes] = []
    newline_count = 0
    with path.open("rb") as handle:
        handle.seek(0, os.SEEK_END)
        position = handle.tell()
        while position > 0 and newline_count <= lines:
            read_size = min(block_size, position)
            position -= read_size
            handle.seek(position)
            chunk = handle.read(read_size)
            chunks.append(chunk)
            newline_count += chunk.count(b"\n")
    selected = (
        b"".join(reversed(chunks)).decode(errors="replace").splitlines()[-lines:]
    )
    return "\n".join(selected) + ("\n" if selected else "")


def _terminate_process(process: subprocess.Popen[bytes]) -> None:
    if process.poll() is not None:
        return
    if os.name == "nt":
        process.terminate()
    else:
        try:
            # Core's Bash batches own separate process groups. SIGINT lets
            # Python unwind their finally blocks and close those clients;
            # already running external services are not managed here.
            os.killpg(process.pid, signal.SIGINT)
        except (OSError, ProcessLookupError):
            process.terminate()
    try:
        process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        if os.name == "nt":
            process.kill()
        else:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except (OSError, ProcessLookupError):
                process.kill()
        process.wait(timeout=10)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="morse-backend",
        description="Authenticated asynchronous job API for MORSE.",
    )
    parser.add_argument(
        "--host",
        default=os.environ.get("MORSE_BACKEND_HOST", "127.0.0.1"),
    )
    parser.add_argument(
        "--port",
        type=int,
        default=int(os.environ.get("MORSE_BACKEND_PORT", "8000")),
    )
    parser.add_argument(
        "--config-root",
        type=Path,
        default=Path(
            os.environ.get("MORSE_BACKEND_CONFIG_ROOT", Path.cwd() / "configs")
        ),
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=Path(
            os.environ.get(
                "MORSE_BACKEND_RESULTS_ROOT",
                Path.cwd() / "results" / "backend",
            )
        ),
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=(
            None
            if not os.environ.get("MORSE_BACKEND_WORKDIR")
            else Path(os.environ["MORSE_BACKEND_WORKDIR"])
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    settings = BackendSettings.from_values(
        config_root=args.config_root,
        results_root=args.results_root,
        workdir=args.workdir,
    )
    manager = JobManager(settings)
    server = MorseHTTPServer(
        (args.host, args.port),
        manager,
        settings.token,
    )

    def stop_server(_signum: int, _frame: Any) -> None:
        threading.Thread(target=server.shutdown, daemon=True).start()

    previous_sigterm = signal.signal(signal.SIGTERM, stop_server)
    print(
        f"MORSE backend listening on http://{args.host}:{args.port} "
        f"configs={settings.config_root} results={settings.results_root}",
        file=sys.stderr,
        flush=True,
    )
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        signal.signal(signal.SIGTERM, previous_sigterm)
        server.shutdown()
        server.server_close()
        manager.shutdown()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
