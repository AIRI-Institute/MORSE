"""RNA validation, FASTA input and homopolymer spans for scorer tools."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

RNA_BASES = frozenset("ACGU")


@dataclass(frozen=True)
class SequenceRun:
    base: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


def normalize_rna(sequence: str) -> str:
    normalized = "".join(sequence.upper().split()).replace("T", "U")
    invalid = sorted(set(normalized) - RNA_BASES)
    if invalid:
        raise ValueError(f"Sequence contains unsupported bases: {', '.join(invalid)}")
    return normalized


def validate_forbidden_motifs(motifs: Sequence[str]) -> tuple[str, ...]:
    if isinstance(motifs, (str, bytes)) or not isinstance(motifs, Sequence):
        raise ValueError("forbidden_motifs must be a list of sequences")
    validated: list[str] = []
    seen_motifs: set[str] = set()
    for index, motif in enumerate(motifs):
        if not isinstance(motif, str):
            raise ValueError(f"forbidden_motifs[{index}] must be a string")
        if not motif:
            raise ValueError(f"forbidden_motifs[{index}] must not be empty")
        invalid = sorted(set(motif) - RNA_BASES)
        if invalid:
            raise ValueError(
                f"forbidden_motifs[{index}] contains unsupported RNA bases: "
                + ", ".join(invalid)
            )
        if motif in seen_motifs:
            raise ValueError(f"forbidden_motifs contains a duplicate: {motif!r}")
        seen_motifs.add(motif)
        validated.append(motif)
    return tuple(validated)


def homopolymer_runs(sequence: str) -> list[SequenceRun]:
    if not sequence:
        return []
    runs: list[SequenceRun] = []
    start = 0
    for index in range(1, len(sequence) + 1):
        if index == len(sequence) or sequence[index] != sequence[start]:
            runs.append(SequenceRun(sequence[start], start, index))
            start = index
    return runs


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts)))
            header = line[1:].strip()
            if not header:
                raise ValueError("FASTA header must contain an identifier")
            name = header.split()[0]
            parts = []
        else:
            if name is None:
                raise ValueError("FASTA sequence appears before the first header")
            parts.append(line)
    if name is not None:
        records.append((name, "".join(parts)))
    if not records:
        raise ValueError(f"FASTA contains no records: {path}")
    return records
