"""Small, dependency-free readers for the frozen scientific bundle."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
from pathlib import Path

HERE = Path(os.environ.get("MORSE_ROOT", Path(__file__).resolve().parent.parent)).expanduser().resolve() / "paper"
REPO = HERE.parent
DATA = HERE / "data"


def read_json(path):
    path = Path(path)
    with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as stream:
        return json.load(stream)


def rows(path):
    path = Path(path)
    with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open(newline="")) as stream:
        return list(csv.DictReader(stream))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def write_csv(path, records):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not records:
        raise ValueError(f"Refusing to emit an empty scientific table: {path}")
    fields = list(dict.fromkeys(k for row in records for k in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def registry():
    return read_json(HERE / "registry.json")


def runs():
    return read_json(DATA / "analysis_runs.json.gz")
