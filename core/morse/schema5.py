"""Schema 5: explicit initiation/termination choices around a CDS body."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from morse.codon import normalize_dna
from morse.construct import PART_ORDER, STOP_CODONS, ConstructDesignSpace, construct_design_from_config
from morse.schema4 import _apply_changes, _mapping, _reject
from morse.schema4 import engine_config, resolve_paths as _resolve_paths


def _codon_choices(raw: Any, part: str) -> tuple[str, ...]:
    context = f"construct.{part}"
    if not isinstance(raw, list) or not raw or any(not isinstance(item, str) for item in raw):
        raise ValueError(f"{context} must be a non-empty list of codon strings")
    values = tuple(normalize_dna(item) for item in raw)
    if part == "start_codon":
        if any(len(value) != 3 or value in STOP_CODONS for value in values):
            raise ValueError(f"{context} must contain three-nucleotide non-stop codons")
    elif any(value not in STOP_CODONS for value in values):
        raise ValueError(f"{context} must contain UAA, UAG or UGA")
    if len(values) != len(set(values)):
        raise ValueError(f"{context} contains duplicate codons")
    return values


def resolve_paths(raw: dict[str, Any], base_dir: Path, *, output_base_dir: Path | None = None) -> None:
    """Resolve paths and record normalized RNA codon defaults in the snapshot."""
    construct = _mapping(raw.get("construct"), "construct")
    for part, default in (("start_codon", ["AUG"]), ("stop_codon", ["UAA"])):
        values = _codon_choices(construct.get(part, default), part)
        construct[part] = [value.replace("T", "U") for value in values]
    _resolve_paths(raw, base_dir, output_base_dir=output_base_dir)


def build_design(
    name: str, raw: dict[str, Any], base_dir: Path, *, source_sequence: str | None = None
) -> ConstructDesignSpace:
    _reject(raw, {"fasta", "coordinates", "changes", "allow_upstream_start_codons",
                  "start_codon", "stop_codon"}, "construct")
    changes = _mapping(raw.get("changes", {}), "construct.changes")
    if {"start_codon", "stop_codon"} & set(changes):
        raise ValueError(
            "Schema 5 defines codon choices only in construct.start_codon and "
            "construct.stop_codon; remove them from construct.changes"
        )
    _reject(changes, set(PART_ORDER) - {"start_codon", "stop_codon"}, "construct.changes")
    starts = _codon_choices(raw.get("start_codon", ["AUG"]), "start_codon")
    stops = _codon_choices(raw.get("stop_codon", ["UAA"]), "stop_codon")
    mutable = set(changes)
    if len(starts) > 1:
        mutable.add("start_codon")
    if len(stops) > 1:
        mutable.add("stop_codon")
    source = {key: value for key, value in raw.items()
              if key not in {"changes", "start_codon", "stop_codon"}}
    design = construct_design_from_config(
        name, source, base_dir, mutable, source_sequence=source_sequence,
        codon_choices=(starts, stops),
    )
    return _apply_changes(design, changes, base_dir)
