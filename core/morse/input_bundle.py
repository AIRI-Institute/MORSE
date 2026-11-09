"""Text inputs accompanying a notebook run configuration."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
from typing import Any


def input_references(configuration):
    """Yield editable paths to sequence, library and initialization inputs."""
    construct = configuration.get("construct", {})
    if isinstance(construct, dict):
        if isinstance(construct.get("fasta"), str):
            yield construct, "fasta"
        changes = construct.get("changes", {})
        if not isinstance(changes, dict):
            raise ValueError("construct.changes must be a mapping")
        for change in changes.values():
            if not isinstance(change, dict):
                continue
            library = change.get("library")
            if isinstance(library, str):
                yield change, "library"
            elif isinstance(library, dict) and isinstance(library.get("path"), str):
                yield library, "path"
        parts = construct.get("parts", {})
        if not isinstance(parts, dict):
            raise ValueError("construct.parts must be a mapping")
        for part in parts.values():
            source = part.get("source", {}) if isinstance(part, dict) else {}
            if source.get("type") in {"fasta", "tsv"} and isinstance(source.get("path"), str):
                yield source, "path"
    optimizer = configuration.get("optimizer", {})
    if isinstance(optimizer, dict):
        for key in ("initial_population", "initial_population_audit"):
            if isinstance(optimizer.get(key), str):
                yield optimizer, key
    initialization = configuration.get("initialization", {})
    if isinstance(initialization, dict):
        for key in ("path", "audit", "profile_path"):
            if isinstance(initialization.get(key), str):
                yield initialization, key


def relative_input_name(name: str) -> Path:
    path = Path(name)
    if path.is_absolute() or not path.parts or any(part in {"", ".", ".."} for part in path.parts) or "\\" in name:
        raise ValueError("Input filenames must be relative paths without traversal")
    return path


def stage_input_bundle(configuration: dict[str, Any], files: Any, cache: Path) -> None:
    """Stage bounded text inputs immutably and rewrite only declared input paths."""
    if not isinstance(files, dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in files.items()):
        raise ValueError("input_files must map relative filenames to UTF-8 text")
    if sum(len(content.encode("utf-8")) for content in files.values()) > 10 * 1024 * 1024:
        raise ValueError("Notebook sequence inputs exceed the 10 MiB upload limit")
    references = list(input_references(configuration))
    names = {mapping[key] for mapping, key in references}
    if set(files) - names:
        raise ValueError("input_files contains files not referenced by the configuration")
    relative = {name: relative_input_name(name) for name in files}
    identity = hashlib.sha256(json.dumps(files, sort_keys=True, ensure_ascii=False).encode()).hexdigest()
    root = cache / identity
    root.mkdir(parents=True, exist_ok=True)
    for name, content in files.items():
        target = root / relative[name]
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if target.read_text(encoding="utf-8") != content:
                raise ValueError("Cached notebook input differs from its recorded content")
        else:
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=target.parent, delete=False) as handle:
                handle.write(content)
                pending = Path(handle.name)
            pending.replace(target)
    for mapping, key in references:
        if mapping[key] in files:
            mapping[key] = str((root / relative[mapping[key]]).resolve())
