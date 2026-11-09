from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path

from morse.codon import normalize_dna

REQUIRED_PART_COLUMNS = {"id", "type", "sequence", "set"}


@dataclass(frozen=True)
class PartRecord:
    id: str
    type: str
    sequence: str
    set_name: str
    metadata: dict[str, str]


def load_part_library(
    path: Path,
    expected_type: str,
    set_name: str | None = None,
    max_candidates: int | None = None,
) -> tuple[PartRecord, ...]:
    if max_candidates is not None and max_candidates < 1:
        raise ValueError("library.max_candidates must be >= 1")
    if not path.exists():
        raise FileNotFoundError(f"Parts library does not exist: {path}")

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = set(reader.fieldnames or [])
        missing = sorted(REQUIRED_PART_COLUMNS - columns)
        if missing:
            raise ValueError(
                f"Parts library {path} is missing columns: {', '.join(missing)}"
            )
        records = [
            _record_from_row(row, path)
            for row in reader
            if row.get("type", "").strip() == expected_type
            and (set_name is None or row.get("set", "").strip() == set_name)
        ]

    if not records:
        selection = f" set={set_name!r}" if set_name is not None else ""
        raise ValueError(f"Parts library {path} has no {expected_type!r} rows{selection}")
    _reject_duplicates(records, path)
    if max_candidates is not None:
        records = records[:max_candidates]
    return tuple(records)


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def _record_from_row(row: dict[str, str], path: Path) -> PartRecord:
    part_id = row.get("id", "").strip()
    part_type = row.get("type", "").strip()
    set_name = row.get("set", "").strip()
    if not part_id or not part_type or not set_name:
        raise ValueError(f"Parts library {path} contains an empty id, type, or set")
    sequence = normalize_dna(row.get("sequence", ""))
    if not sequence:
        raise ValueError(f"Parts library {path} contains an empty sequence for {part_id}")
    expected_hash = row.get("sequence_hash", "").strip()
    actual_hash = sequence_hash(sequence)
    if expected_hash and expected_hash != actual_hash:
        raise ValueError(f"Parts library sequence hash mismatch for {part_id}")
    metadata = {
        key: value.strip()
        for key, value in row.items()
        if key not in REQUIRED_PART_COLUMNS and value is not None and value.strip()
    }
    return PartRecord(part_id, part_type, sequence, set_name, metadata)


def _reject_duplicates(records: list[PartRecord], path: Path) -> None:
    ids = [record.id for record in records]
    if len(ids) != len(set(ids)):
        raise ValueError(f"Parts library {path} selection contains duplicate ids")
    sequences = [record.sequence for record in records]
    if len(sequences) != len(set(sequences)):
        raise ValueError(f"Parts library {path} selection contains duplicate sequences")
