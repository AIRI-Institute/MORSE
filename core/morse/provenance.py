"""Content-based provenance for inference and safe optimizer continuation.

Hash before inference, never substitute a missing asset with a fake scorer.
Custom scorers must declare additional imported code, data, and remote model
snapshots in ``provenance_paths`` when these cannot be inferred from their command.
"""
from __future__ import annotations

import copy
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from typing import Any

from morse.config import RunConfig, ScorerConfig

_SKIP = {".git", ".venv", "venv", "__pycache__", ".pytest_cache", ".ipynb_checkpoints"}


def content_digest(path: Path) -> dict[str, Any]:
    """Stream all file content; directory identity is independent of location."""
    path = path.expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Missing provenance input: {path}")
    files = [path] if path.is_file() else sorted(
        p for p in path.rglob("*")
        if p.is_file() and not any(x in _SKIP for x in p.relative_to(path).parts)
        and p.suffix not in {".pyc", ".pyo"}
    )
    digest = hashlib.sha256()
    total = 0
    single_file_hash = None
    for item in files:
        relative = item.relative_to(path).as_posix() if path.is_dir() else "file"
        file_hash = hashlib.sha256()
        with item.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                file_hash.update(chunk)
                total += len(chunk)
        digest.update(relative.encode() + b"\0" + file_hash.digest())
        single_file_hash = file_hash.hexdigest()
    return {"sha256": single_file_hash if path.is_file() else digest.hexdigest(), "files": len(files), "bytes": total}


def resolve_image_id(image: str) -> str:
    try:
        result = subprocess.run(
            ["docker", "image", "inspect", image, "--format", "{{.Id}}"],
            check=True, text=True, capture_output=True, timeout=60,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise RuntimeError(
            f"Cannot resolve immutable identity for image {image!r}; build or pull "
            "the image first, then retry. Docker must be available."
        ) from exc
    value = result.stdout.strip()
    if not value.startswith("sha256:"):
        raise ValueError(f"Docker returned an invalid image identity for {image!r}")
    return value


def scorer_provenance(scorer: ScorerConfig) -> dict[str, Any]:
    base = Path(scorer.workdir or Path.cwd()).expanduser().resolve()
    paths: dict[str, Path] = {}
    for index, value in enumerate(scorer.provenance_paths):
        raw = Path(value).expanduser()
        paths[f"declared:{index}"] = raw if raw.is_absolute() else base / raw
    if isinstance(scorer.command, str):
        from morse.scorers import BASH_EXECUTABLE
        # Shell code is deliberately opaque. Inspecting the Docker client or
        # guessing URLs/paths in a pipeline cannot verify the service's model,
        # source or state. Only explicitly declared host inputs are hashed.
        shell = Path(BASH_EXECUTABLE)
        environment = {"MORSE_PYTHON": os.environ.get("MORSE_PYTHON", sys.executable), **(scorer.env or {})}
        result = {
            "name": scorer.name, "mode": scorer.mode,
            "command": scorer.command, "workdir": str(base), "env": environment,
            "shell": {"path": str(shell), **content_digest(shell)},
            "inputs": {name: {"path": str(path.resolve()), **content_digest(path)}
                       for name, path in sorted(paths.items())},
            "external_runtime_identity_verified": False,
            "verification_scope": "Exact Bash command, working directory, declared environment and declared host input files only. External service/container code, model assets and state are not verified.",
        }
        if scorer.parameters:
            result["parameters"] = copy.deepcopy(scorer.parameters)
        return result
    command = list(scorer.command or [])
    if command:
        executable = command[0]
        if executable in {"python", "python3", "@python"}:
            executable = sys.executable
        resolved_executable = shutil.which(executable) or executable
        candidate = Path(resolved_executable)
        if not candidate.is_absolute():
            candidate = base / candidate
        if candidate.is_file():
            paths["executable"] = candidate
        for index, token in enumerate(command[1:], start=1):
            if token.startswith("-"):
                continue
            candidate = Path(token).expanduser()
            candidate = candidate if candidate.is_absolute() else base / candidate
            if candidate.is_file():
                paths[f"command:{index}"] = candidate
        if "-m" in command:
            index = command.index("-m")
            if index + 1 < len(command):
                module = command[index + 1]
                try:
                    spec = importlib.util.find_spec(module)
                except (ImportError, ModuleNotFoundError, ValueError):
                    spec = None
                if spec is not None and spec.origin and spec.origin not in {"built-in", "frozen"}:
                    origin = Path(spec.origin)
                    # A module may import siblings; hash that package's sources.
                    paths["python_module"] = origin.parent if "." in module else origin
    # Explicit env inputs (not the whole process environment) commonly carry
    # model/profile paths. Container-internal paths are covered by mounts/image.
    for name, value in (scorer.env or {}).items():
        raw = Path(value).expanduser()
        candidate = raw if raw.is_absolute() else base / raw
        if candidate.is_file() or candidate.is_dir():
            paths[f"env:{name}"] = candidate
    args = list(scorer.docker_args or [])
    for index, token in enumerate(args):
        value = None
        if token in {"-v", "--volume", "--mount"} and index + 1 < len(args):
            value = args[index + 1]
        elif token.startswith(("--volume=", "--mount=")):
            value = token.split("=", 1)[1]
        if value is None:
            continue
        if "type=bind" in value:
            fields = dict(part.split("=", 1) for part in value.split(",") if "=" in part)
            source = fields.get("source") or fields.get("src")
        elif ":" in value:
            source = value.split(":", 1)[0]
        else:
            source = None
        if source:
            candidate = Path(source).expanduser()
            candidate = candidate if candidate.is_absolute() else base / candidate
            if not candidate.exists():
                raise FileNotFoundError(f"Missing bind-mounted asset: {candidate}")
            paths[f"mount:{index}"] = candidate
    inputs = {name: {"path": str(path.resolve()), **content_digest(path)} for name, path in sorted(paths.items())}
    result: dict[str, Any] = {"name": scorer.name, "mode": scorer.mode, "inputs": inputs}
    if scorer.parameters:
        result["parameters"] = copy.deepcopy(scorer.parameters)
    if scorer.image:
        result["image"] = scorer.image
        result["image_id"] = resolve_image_id(scorer.image)
    return result


def runtime_provenance(config: RunConfig) -> dict[str, Any]:
    from morse import __version__
    versions = {}
    for name in ("numpy", "pymoo", "PyYAML", "cloudpickle", "torch", "ViennaRNA"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    core_source = content_digest(Path(__file__).parent)
    return {
        "schema_version": 1,
        "core_version": __version__,
        "core_source": core_source,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": versions,
        "scorers": [scorer_provenance(s) for s in config.scorers],
    }


def scoring_fingerprint(config: RunConfig, provenance: dict[str, Any]) -> str:
    payload = {"config": {key: config.to_dict()[key] for key in ("scorers", "manufacturability")}, "runtime": provenance}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
