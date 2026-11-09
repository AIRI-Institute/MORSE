from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
import subprocess
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Mapping, Sequence

from morse.codon import (
    DNA_CODON_TABLE,
    SYNONYMOUS_CODONS,
    normalize_dna,
    translate_dna,
)
from morse.construct import ConstructDesignSpace

from morse_warm_starter.generation import GeneratedCandidate

LINEARDESIGN_REPOSITORY = (
    "https://github.com/LinearDesignSoftware/LinearDesign.git"
)
LINEARDESIGN_COMMIT = "f0126ca89a8b853088b4bccfd2cc8c378d3678be"
LINEARDESIGN_LICENSE_URL = (
    "https://github.com/LinearDesignSoftware/LinearDesign/blob/main/license.txt"
)


@dataclass(frozen=True)
class CaiProfile:
    profile_id: str
    taxon: str
    codon_counts: Mapping[str, float]
    sha256: str


def generate_lineardesign_candidates(
    design_space: ConstructDesignSpace,
    *,
    profile_path: str | Path,
    lambdas: Sequence[float],
    lineardesign_dir: str | Path,
    work_dir: str | Path,
    accept_license: bool,
) -> list[GeneratedCandidate]:
    profile = load_profile(profile_path)
    binary = install_lineardesign(
        lineardesign_dir,
        accept_license=accept_license,
    )
    csv_path = write_lineardesign_profile(
        profile,
        Path(work_dir) / f"{profile.profile_id}.lineardesign.csv",
    )
    candidates: list[GeneratedCandidate] = []
    seen: set[tuple[int, ...]] = set()
    for index, raw_lambda in enumerate(lambdas, start=1):
        lambda_value = float(raw_lambda)
        if not math.isfinite(lambda_value) or lambda_value < 0:
            raise ValueError("LinearDesign lambda values must be finite and >= 0")
        completed = subprocess.run(
            [str(binary), str(lambda_value), "1", str(csv_path)],
            input=design_space.amino_acids + "\n",
            text=True,
            capture_output=True,
            check=True,
            cwd=Path(lineardesign_dir),
        )
        coding_sequence, reported_mfe, reported_cai = _parse_output(
            completed.stdout + "\n" + completed.stderr,
            design_space.amino_acids,
        )
        source = design_space.source_construct
        construct = replace(
            source,
            start_codon=coding_sequence[:3],
            cds=coding_sequence[3:],
        )
        vector = tuple(design_space.encode(construct))
        if vector in seen:
            continue
        seen.add(vector)
        candidates.append(
            GeneratedCandidate(
                id=f"lineardesign-{index:03d}",
                vector=vector,
                sequence=construct.sequence,
                metadata={
                    "anchor": "lineardesign",
                    "lambda": lambda_value,
                    "reported_mfe": reported_mfe,
                    "reported_cai": reported_cai,
                    "profile_id": profile.profile_id,
                    "profile_sha256": profile.sha256,
                    "lineardesign_commit": LINEARDESIGN_COMMIT,
                },
            )
        )
    return candidates


def install_lineardesign(
    directory: str | Path,
    *,
    accept_license: bool,
) -> Path:
    if not accept_license:
        raise ValueError(
            "Read and accept the upstream LinearDesign academic/non-profit "
            f"license before installation: {LINEARDESIGN_LICENSE_URL}"
        )
    root = Path(directory).expanduser().resolve()
    binary = root / "bin" / "LinearDesign_2D"
    marker = root / ".morse_upstream_commit"
    if (
        binary.is_file()
        and marker.is_file()
        and marker.read_text().strip() == LINEARDESIGN_COMMIT
    ):
        return binary
    if root.exists():
        raise ValueError(
            f"LinearDesign directory exists but is not the pinned build: {root}"
        )
    subprocess.run(
        [
            "git",
            "clone",
            "--filter=blob:none",
            "--",
            LINEARDESIGN_REPOSITORY,
            str(root),
        ],
        check=True,
    )
    subprocess.run(
        ["git", "checkout", "--detach", LINEARDESIGN_COMMIT],
        cwd=root,
        check=True,
    )
    subprocess.run(["make", "lineardesign_2D"], cwd=root, check=True)
    if not binary.is_file():
        raise RuntimeError(f"LinearDesign build did not create {binary}")
    marker.write_text(LINEARDESIGN_COMMIT + "\n")
    return binary


def load_profile(path: str | Path) -> CaiProfile:
    source = Path(path).expanduser().resolve()
    payload = source.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    if source.suffix.lower() == ".json":
        raw = json.loads(payload)
        if not isinstance(raw, Mapping):
            raise ValueError("CAI profile JSON must contain an object")
        if int(raw.get("schema_version", 0)) != 1:
            raise ValueError("CAI profile must define schema_version: 1")
        profile_id = str(raw.get("profile_id", "")).strip()
        taxon = str(raw.get("taxon", "")).strip()
        counts = _validate_counts(raw.get("codon_counts"))
    elif source.suffix.lower() in {".csv", ".tsv"}:
        profile_id = source.stem
        taxon = "custom"
        counts = _counts_from_lineardesign_csv(
            payload.decode("utf-8-sig")
        )
    else:
        raise ValueError("CAI profile must be JSON, CSV, or TSV")
    if not profile_id or not taxon:
        raise ValueError("CAI profile must define profile_id and taxon")
    return CaiProfile(profile_id, taxon, counts, digest)


def write_lineardesign_profile(
    profile: CaiProfile,
    path: str | Path,
) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    lines = ["#,,"]
    for amino_acid in sorted(SYNONYMOUS_CODONS):
        codons = SYNONYMOUS_CODONS[amino_acid]
        total = math.fsum(profile.codon_counts[codon] for codon in codons)
        if total <= 0:
            raise ValueError(
                f"CAI profile has no frequency mass for {amino_acid}"
            )
        for codon in codons:
            lines.append(
                f"{codon.replace('T', 'U')},{amino_acid},"
                f"{profile.codon_counts[codon] / total:.12g}"
            )
    destination.write_text("\n".join(lines) + "\n")
    return destination


def _counts_from_lineardesign_csv(text: str) -> dict[str, float]:
    counts: dict[str, float] = {}
    for row in list(csv.reader(io.StringIO(text)))[1:]:
        if not row or not any(cell.strip() for cell in row):
            continue
        if len(row) != 3:
            raise ValueError(
                "LinearDesign CSV rows must contain codon, amino acid, and frequency"
            )
        codon = normalize_dna(row[0])
        expected = DNA_CODON_TABLE.get(codon)
        amino_acid = row[1].strip().upper()
        if expected is None or amino_acid not in {
            expected,
            "STOP" if expected == "*" else expected,
        }:
            raise ValueError(
                f"Amino acid {amino_acid!r} does not match codon {codon}"
            )
        counts[codon] = float(row[2])
    return _validate_counts(counts)


def _validate_counts(raw: Any) -> dict[str, float]:
    if not isinstance(raw, Mapping):
        raise ValueError("CAI profile must define codon_counts")
    counts: dict[str, float] = {}
    for raw_codon, raw_value in raw.items():
        codon = normalize_dna(str(raw_codon))
        if len(codon) != 3 or codon not in DNA_CODON_TABLE:
            raise ValueError(f"Unknown codon: {raw_codon!r}")
        value = float(raw_value)
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"Invalid frequency for {codon}")
        counts[codon] = value
    missing = sorted(set(DNA_CODON_TABLE) - set(counts))
    if missing:
        raise ValueError("CAI profile is missing codons: " + ", ".join(missing))
    return counts


def _parse_output(
    output: str,
    protein: str,
) -> tuple[str, float | None, float | None]:
    match = re.search(r"mRNA sequence:\s*([ACGUT]+)", output)
    if not match:
        raise RuntimeError(
            "LinearDesign output does not contain an mRNA sequence: "
            + output[-2000:]
        )
    sequence = normalize_dna(match.group(1))
    if translate_dna(sequence) != protein:
        raise RuntimeError("LinearDesign output does not preserve the protein")
    metrics = re.search(
        r"folding free energy:\s*([-+0-9.eE]+).*?mRNA CAI:\s*([-+0-9.eE]+)",
        output,
        re.DOTALL,
    )
    return (
        sequence,
        float(metrics.group(1)) if metrics else None,
        float(metrics.group(2)) if metrics else None,
    )
