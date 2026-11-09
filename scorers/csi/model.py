from __future__ import annotations

import bisect
import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import Any

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
SENSE_CODONS = frozenset(DNA_CODON_TABLE) - STOP_CODONS


@dataclass(frozen=True)
class CsiProfile:
    profile_id: str
    organism: str
    context: str
    coefficient_method: str
    codon_stability_coefficients: Mapping[str, float]
    reference_csi: tuple[float, ...]
    reference_percentiles: tuple[float, ...]
    sha256: str
    source: Mapping[str, Any]


@dataclass(frozen=True)
class CsiScore:
    csi_raw: float
    csi_percentile: float
    window_p05: float
    window_minimum: float
    destabilizing_fraction: float
    codons: int


class CsiScorer:
    def __init__(
        self,
        profile: CsiProfile,
        *,
        window_codons: int = 15,
        destabilizing_threshold: float = 0.0,
    ) -> None:
        if isinstance(window_codons, bool) or int(window_codons) < 1:
            raise ValueError("window_codons must be an integer >= 1")
        if not math.isfinite(destabilizing_threshold):
            raise ValueError("destabilizing_threshold must be finite")
        self.profile = profile
        self.window_codons = int(window_codons)
        self.destabilizing_threshold = float(destabilizing_threshold)

    def score_sequence(self, sequence: str) -> CsiScore:
        coding_sequence = normalize_dna(sequence)
        codons = split_codons(coding_sequence)
        if len(codons) < 2:
            raise ValueError(
                "CSI requires at least one CDS codon after the explicit start codon"
            )
        stop_positions = [
            index for index, codon in enumerate(codons) if codon in STOP_CODONS
        ]
        if stop_positions:
            positions = ", ".join(str(position + 1) for position in stop_positions)
            raise ValueError(
                "Coding sequence passed to CSI must not contain stop codons; "
                f"found at codon position(s) {positions}"
            )

        values = [
            self.profile.codon_stability_coefficients[codon] for codon in codons[1:]
        ]
        csi_raw = math.fsum(values) / len(values)
        window_size = min(self.window_codons, len(values))
        window_scores = [
            math.fsum(values[start : start + window_size]) / window_size
            for start in range(len(values) - window_size + 1)
        ]
        return CsiScore(
            csi_raw=csi_raw,
            csi_percentile=self.profile_percentile(csi_raw),
            window_p05=quantile(window_scores, 0.05),
            window_minimum=min(window_scores),
            destabilizing_fraction=sum(
                value < self.destabilizing_threshold for value in values
            )
            / len(values),
            codons=len(values),
        )

    def profile_percentile(self, csi: float) -> float:
        values = self.profile.reference_csi
        percentiles = self.profile.reference_percentiles
        if csi <= values[0]:
            return percentiles[0]
        if csi >= values[-1]:
            return percentiles[-1]
        upper = bisect.bisect_right(values, csi)
        lower = upper - 1
        fraction = (csi - values[lower]) / (values[upper] - values[lower])
        return percentiles[lower] + fraction * (percentiles[upper] - percentiles[lower])

    def score_outputs(self, sequence: str) -> dict[str, Any]:
        result = self.score_sequence(sequence)
        return {
            "csi_raw": result.csi_raw,
            "csi_percentile": result.csi_percentile,
            "csi_window_p05": result.window_p05,
            "csi_window_min": result.window_minimum,
            "csi_destabilizing_fraction": result.destabilizing_fraction,
            "csi_codons": result.codons,
            "csi_profile_id": self.profile.profile_id,
            "csi_profile_sha256": self.profile.sha256,
        }


def load_profile(path: Path) -> CsiProfile:
    payload = path.read_bytes()
    try:
        raw = json.loads(payload)
    except json.JSONDecodeError as exc:
        raise ValueError(f"CSI profile is not valid JSON: {path}") from exc
    if not isinstance(raw, dict):
        raise TypeError(f"CSI profile must contain a JSON object: {path}")
    return profile_from_mapping(raw, hashlib.sha256(payload).hexdigest())


def profile_from_mapping(
    raw: Mapping[str, Any],
    sha256: str | None = None,
) -> CsiProfile:
    if int(raw.get("schema_version", 0)) != 1:
        raise ValueError("CSI profile must define schema_version: 1")
    if int(raw.get("genetic_code", 0)) != 1:
        raise ValueError("CSI scorer supports only NCBI genetic code 1")

    profile_id = str(raw.get("profile_id", "")).strip()
    organism = str(raw.get("organism", "")).strip()
    context = str(raw.get("context", "")).strip()
    coefficient_method = str(raw.get("coefficient_method", "")).strip()
    if not all((profile_id, organism, context, coefficient_method)):
        raise ValueError(
            "CSI profile must define profile_id, organism, context, and "
            "coefficient_method"
        )

    coefficients_raw = raw.get("codon_stability_coefficients")
    if not isinstance(coefficients_raw, dict):
        raise TypeError(
            "CSI profile must define codon_stability_coefficients as an object"
        )
    coefficients: dict[str, float] = {}
    for raw_codon, raw_value in coefficients_raw.items():
        codon = normalize_dna(str(raw_codon))
        if len(codon) != 3 or codon not in SENSE_CODONS:
            raise ValueError(f"Unknown sense codon in CSI profile: {raw_codon!r}")
        if codon in coefficients:
            raise ValueError(
                f"Duplicate codon in CSI profile after normalization: {codon}"
            )
        coefficients[codon] = _finite_number(
            raw_value, f"codon_stability_coefficients.{codon}"
        )
    missing = sorted(SENSE_CODONS - set(coefficients))
    extra = sorted(set(coefficients) - SENSE_CODONS)
    if missing or extra:
        detail = []
        if missing:
            detail.append("missing: " + ", ".join(missing))
        if extra:
            detail.append("extra: " + ", ".join(extra))
        raise ValueError(
            "CSI profile must define exactly 61 sense codons ("
            + "; ".join(detail)
            + ")"
        )

    cdf_raw = raw.get("reference_cdf")
    if not isinstance(cdf_raw, list) or len(cdf_raw) < 2:
        raise ValueError("CSI profile reference_cdf must contain at least two points")
    reference_csi: list[float] = []
    reference_percentiles: list[float] = []
    for index, point in enumerate(cdf_raw):
        if not isinstance(point, dict) or set(point) != {"csi", "percentile"}:
            raise ValueError(
                f"reference_cdf[{index}] must define only csi and percentile"
            )
        reference_csi.append(
            _finite_number(point["csi"], f"reference_cdf[{index}].csi")
        )
        percentile = _finite_number(
            point["percentile"], f"reference_cdf[{index}].percentile"
        )
        if not 0.0 <= percentile <= 100.0:
            raise ValueError("reference_cdf percentiles must be between 0 and 100")
        reference_percentiles.append(percentile)
    if reference_percentiles[0] != 0.0 or reference_percentiles[-1] != 100.0:
        raise ValueError("reference_cdf must start at percentile 0 and end at 100")
    if any(current <= previous for previous, current in pairwise(reference_csi)):
        raise ValueError("reference_cdf csi values must be strictly increasing")
    if any(
        current <= previous for previous, current in pairwise(reference_percentiles)
    ):
        raise ValueError("reference_cdf percentiles must be strictly increasing")

    source_raw = raw.get("source", {})
    if not isinstance(source_raw, dict):
        raise TypeError("CSI profile source must be an object")
    if sha256 is None:
        canonical = json.dumps(
            raw, sort_keys=True, separators=(",", ":"), allow_nan=False
        )
        sha256 = hashlib.sha256(canonical.encode()).hexdigest()
    return CsiProfile(
        profile_id=profile_id,
        organism=organism,
        context=context,
        coefficient_method=coefficient_method,
        codon_stability_coefficients=coefficients,
        reference_csi=tuple(reference_csi),
        reference_percentiles=tuple(reference_percentiles),
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


def quantile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def _finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{label} must be a finite number")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a finite number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be a finite number")
    return result
