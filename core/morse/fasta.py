from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from textwrap import wrap


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str


def read_fasta(path: Path) -> FastaRecord:
    header: str | None = None
    parts: list[str] = []

    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                raise ValueError(f"{path} contains more than one FASTA record")
            header = line[1:].strip()
        else:
            parts.append(line)

    if header is None:
        raise ValueError(f"{path} does not contain a FASTA header")

    sequence = "".join(parts).upper().replace("U", "T")
    return FastaRecord(header=header, sequence=sequence)


def write_fasta(path: Path, records: list[FastaRecord], width: int = 80) -> None:
    lines: list[str] = []
    for record in records:
        lines.append(f">{record.header}")
        lines.extend(wrap(record.sequence, width))
    path.write_text("\n".join(lines) + "\n")
