from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping


DNA_CODON_TABLE: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

STOP_CODONS = frozenset(
    codon for codon, amino_acid in DNA_CODON_TABLE.items() if amino_acid == "*"
)
UNINFORMATIVE_CODONS = frozenset({"ATG", "TGG"})
SENSE_CODONS = frozenset(DNA_CODON_TABLE) - STOP_CODONS

SYNONYMOUS_CODONS: dict[str, tuple[str, ...]] = {}
for _amino_acid in sorted(set(DNA_CODON_TABLE.values()) - {"*"}):
    SYNONYMOUS_CODONS[_amino_acid] = tuple(
        sorted(codon for codon, value in DNA_CODON_TABLE.items() if value == _amino_acid)
    )


@dataclass(frozen=True)
class CaiProfile:
    profile_id: str
    taxon: str
    reference_kind: str
    zero_count_replacement: float
    codon_counts: Mapping[str, float]
    codon_weights: Mapping[str, float]
    sha256: str
    source: Mapping[str, Any]


@dataclass(frozen=True)
class CaiScore:
    cai: float
    log_cai: float
    minimum_weight: float
    low_weight_fraction: float
    informative_codons: int


class CaiScorer:
    def __init__(self, profile: CaiProfile, low_weight_threshold: float = 0.1) -> None:
        if not math.isfinite(low_weight_threshold) or not 0.0 < low_weight_threshold <= 1.0:
            raise ValueError("low_weight_threshold must satisfy 0 < value <= 1")
        self.profile = profile
        self.low_weight_threshold = low_weight_threshold

    def score_sequence(self, sequence: str) -> CaiScore:
        coding_sequence = normalize_dna(sequence)
        codons = split_codons(coding_sequence)
        if not codons:
            raise ValueError("Coding sequence must not be empty")

        stop_positions = [index for index, codon in enumerate(codons) if codon in STOP_CODONS]
        if stop_positions:
            positions = ", ".join(str(position + 1) for position in stop_positions)
            raise ValueError(
                "Coding sequence passed to CAI must not contain stop codons; "
                f"found at codon position(s) {positions}"
            )

        weights = [
            self.profile.codon_weights[codon]
            for index, codon in enumerate(codons)
            if index > 0 and codon not in UNINFORMATIVE_CODONS
        ]
        if not weights:
            raise ValueError(
                "CAI is undefined when the coding sequence contains no "
                "informative codons after the start codon"
            )

        log_cai = math.fsum(math.log(weight) for weight in weights) / len(weights)
        return CaiScore(
            cai=math.exp(log_cai),
            log_cai=log_cai,
            minimum_weight=min(weights),
            low_weight_fraction=sum(
                weight < self.low_weight_threshold for weight in weights
            )
            / len(weights),
            informative_codons=len(weights),
        )

    def score_outputs(self, sequence: str) -> dict[str, Any]:
        result = self.score_sequence(sequence)
        return {
            "cai": result.cai,
            "cai_log": result.log_cai,
            "cai_min_weight": result.minimum_weight,
            "cai_low_weight_fraction": result.low_weight_fraction,
            "cai_informative_codons": result.informative_codons,
            "cai_profile_id": self.profile.profile_id,
            "cai_profile_sha256": self.profile.sha256,
        }


def load_profile(path: Path) -> CaiProfile:
    payload = path.read_bytes()
    try:
        raw = json.loads(payload)
    except json.JSONDecodeError as exc:
        raise ValueError(f"CAI profile is not valid JSON: {path}") from exc
    if not isinstance(raw, dict):
        raise ValueError(f"CAI profile must contain a JSON object: {path}")
    return profile_from_mapping(raw, hashlib.sha256(payload).hexdigest())


def profile_from_mapping(raw: Mapping[str, Any], sha256: str | None = None) -> CaiProfile:
    if int(raw.get("schema_version", 0)) != 1:
        raise ValueError("CAI profile must define schema_version: 1")
    if int(raw.get("genetic_code", 0)) != 1:
        raise ValueError("CAI scorer supports only NCBI genetic code 1")

    profile_id = str(raw.get("profile_id", "")).strip()
    taxon = str(raw.get("taxon", "")).strip()
    reference_kind = str(raw.get("reference_kind", "")).strip()
    if not profile_id or not taxon or not reference_kind:
        raise ValueError("CAI profile must define profile_id, taxon, and reference_kind")

    zero_count_replacement = _finite_nonnegative(
        raw.get("zero_count_replacement", 0.5), "zero_count_replacement"
    )
    if zero_count_replacement <= 0:
        raise ValueError("zero_count_replacement must be greater than zero")

    counts_raw = raw.get("codon_counts")
    if not isinstance(counts_raw, dict):
        raise ValueError("CAI profile must define codon_counts as an object")
    normalized_counts: dict[str, float] = {}
    for raw_codon, raw_count in counts_raw.items():
        codon = normalize_dna(str(raw_codon))
        if len(codon) != 3 or codon not in DNA_CODON_TABLE:
            raise ValueError(f"Unknown codon in CAI profile: {raw_codon!r}")
        if codon in normalized_counts:
            raise ValueError(f"Duplicate codon in CAI profile after normalization: {codon}")
        normalized_counts[codon] = _finite_nonnegative(raw_count, f"codon_counts.{codon}")

    missing = sorted(SENSE_CODONS - set(normalized_counts))
    if missing:
        raise ValueError("CAI profile is missing sense codons: " + ", ".join(missing))

    weights: dict[str, float] = {}
    for codons in SYNONYMOUS_CODONS.values():
        if math.fsum(normalized_counts[codon] for codon in codons) <= 0:
            amino_acid = DNA_CODON_TABLE[codons[0]]
            raise ValueError(
                f"CAI profile has no observations for amino acid {amino_acid}: {codons}"
            )
        adjusted = {
            codon: (
                normalized_counts[codon]
                if normalized_counts[codon] > 0
                else zero_count_replacement
            )
            for codon in codons
        }
        maximum = max(adjusted.values())
        weights.update({codon: count / maximum for codon, count in adjusted.items()})

    source_raw = raw.get("source", {})
    if not isinstance(source_raw, dict):
        raise ValueError("CAI profile source must be an object")
    if sha256 is None:
        canonical = json.dumps(raw, sort_keys=True, separators=(",", ":"), allow_nan=False)
        sha256 = hashlib.sha256(canonical.encode()).hexdigest()

    return CaiProfile(
        profile_id=profile_id,
        taxon=taxon,
        reference_kind=reference_kind,
        zero_count_replacement=zero_count_replacement,
        codon_counts=normalized_counts,
        codon_weights=weights,
        sha256=sha256,
        source=dict(source_raw),
    )


def normalize_dna(sequence: str) -> str:
    normalized = "".join(sequence.upper().split()).replace("U", "T")
    invalid = sorted(set(normalized) - {"A", "C", "G", "T"})
    if invalid:
        raise ValueError(f"CDS contains unsupported bases: {', '.join(invalid)}")
    return normalized


def split_codons(sequence: str) -> list[str]:
    if len(sequence) % 3 != 0:
        raise ValueError(f"CDS length must be divisible by 3; got {len(sequence)} nt")
    return [sequence[index : index + 3] for index in range(0, len(sequence), 3)]


def _finite_nonnegative(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be a finite non-negative number")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a finite non-negative number") from exc
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"{label} must be a finite non-negative number")
    return result
