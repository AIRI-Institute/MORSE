from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from morse.codon import (
    CodonDesignSpace,
    build_design_space,
    normalize_dna,
    translate_dna,
)
from morse.fasta import read_fasta
from morse.parts import load_part_library

STOP_CODONS = {"TAA", "TAG", "TGA"}
PART_ORDER = (
    "five_prime_utr",
    "start_codon",
    "cds",
    "stop_codon",
    "three_prime_utr",
    "poly_a",
)


@dataclass(frozen=True)
class SequenceVariant:
    name: str
    sequence: str
    metadata: dict[str, str] | None = None


@dataclass(frozen=True)
class VaccineConstruct:
    name: str
    five_prime_utr: SequenceVariant
    start_codon: str
    cds: str
    stop_codon: str
    three_prime_utr: SequenceVariant
    poly_a_length: int

    @property
    def coding_sequence(self) -> str:
        return self.start_codon + self.cds

    @property
    def coding_sequence_with_stop(self) -> str:
        return self.coding_sequence + self.stop_codon

    @property
    def sequence(self) -> str:
        return (
            self.five_prime_utr.sequence
            + self.coding_sequence_with_stop
            + self.three_prime_utr.sequence
            + "A" * self.poly_a_length
        )

    @property
    def segments(self) -> dict[str, str]:
        return {
            "five_prime_utr": self.five_prime_utr.sequence,
            "start_codon": self.start_codon,
            "cds": self.cds,
            "stop_codon": self.stop_codon,
            "three_prime_utr": self.three_prime_utr.sequence,
            "poly_a": "A" * self.poly_a_length,
        }

    @property
    def coordinates(self) -> dict[str, tuple[int, int]]:
        cursor = 0
        coordinates: dict[str, tuple[int, int]] = {}
        for name, sequence in self.segments.items():
            coordinates[name] = (cursor, cursor + len(sequence))
            cursor += len(sequence)
        return coordinates

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "five_prime_utr": {
                "name": self.five_prime_utr.name,
                "sequence": self.five_prime_utr.sequence,
                "metadata": dict(self.five_prime_utr.metadata or {}),
            },
            "start_codon": self.start_codon,
            "cds": self.cds,
            "stop_codon": self.stop_codon,
            "three_prime_utr": {
                "name": self.three_prime_utr.name,
                "sequence": self.three_prime_utr.sequence,
                "metadata": dict(self.three_prime_utr.metadata or {}),
            },
            "poly_a_length": self.poly_a_length,
            "sequence": self.sequence,
            "coordinates": {
                name: {"start": start, "end": end}
                for name, (start, end) in self.coordinates.items()
            },
        }

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> VaccineConstruct:
        five = _variant_from_mapping(raw["five_prime_utr"], "five_prime_utr")
        three = _variant_from_mapping(raw["three_prime_utr"], "three_prime_utr")
        return cls(
            name=str(raw["name"]),
            five_prime_utr=five,
            start_codon=normalize_dna(str(raw["start_codon"])),
            cds=normalize_dna(str(raw["cds"])),
            stop_codon=normalize_dna(str(raw["stop_codon"])),
            three_prime_utr=three,
            poly_a_length=int(raw["poly_a_length"]),
        )

    def fingerprint(self) -> str:
        payload = {
            "name": self.name,
            "five_prime_utr_name": self.five_prime_utr.name,
            "five_prime_utr": self.five_prime_utr.sequence,
            "start_codon": self.start_codon,
            "cds": self.cds,
            "stop_codon": self.stop_codon,
            "three_prime_utr_name": self.three_prime_utr.name,
            "three_prime_utr": self.three_prime_utr.sequence,
            "poly_a_length": self.poly_a_length,
        }
        return hashlib.sha256(_canonical_json(payload).encode()).hexdigest()


@dataclass(frozen=True)
class ConstructDesignSpace:
    name: str
    five_prime_utrs: tuple[SequenceVariant, ...]
    start_codons: tuple[str, ...]
    codons: CodonDesignSpace
    stop_codons: tuple[str, ...]
    three_prime_utrs: tuple[SequenceVariant, ...]
    poly_a_lengths: tuple[int, ...]
    allow_upstream_start_codons: bool = False
    mutable_parts: frozenset[str] = frozenset(PART_ORDER)
    # None preserves the schema-3/4 policy: check every allowed start triplet.
    upstream_start_codons: tuple[str, ...] | None = None

    @property
    def upstream_codon_choices(self) -> tuple[str, ...]:
        return self.upstream_start_codons or self.start_codons

    @property
    def n_var(self) -> int:
        return self.codons.n_var + 5

    @property
    def lower_bounds(self) -> list[int]:
        lower = self.original_vector
        cds_vector_start = 2
        cds_vector_end = cds_vector_start + self.codons.n_var
        if "five_prime_utr" in self.mutable_parts:
            lower[0] = 0
        if "start_codon" in self.mutable_parts:
            lower[1] = 0
        if "cds" in self.mutable_parts:
            lower[cds_vector_start:cds_vector_end] = self.codons.lower_bounds
        if "stop_codon" in self.mutable_parts:
            lower[cds_vector_end] = 0
        if "three_prime_utr" in self.mutable_parts:
            lower[cds_vector_end + 1] = 0
        if "poly_a" in self.mutable_parts:
            lower[cds_vector_end + 2] = 0
        return lower

    @property
    def upper_bounds(self) -> list[int]:
        upper = self.original_vector
        cds_vector_start = 2
        cds_vector_end = cds_vector_start + self.codons.n_var
        if "five_prime_utr" in self.mutable_parts:
            upper[0] = len(self.five_prime_utrs) - 1
        if "start_codon" in self.mutable_parts:
            upper[1] = len(self.start_codons) - 1
        if "cds" in self.mutable_parts:
            upper[cds_vector_start:cds_vector_end] = self.codons.upper_bounds
        if "stop_codon" in self.mutable_parts:
            upper[cds_vector_end] = len(self.stop_codons) - 1
        if "three_prime_utr" in self.mutable_parts:
            upper[cds_vector_end + 1] = len(self.three_prime_utrs) - 1
        if "poly_a" in self.mutable_parts:
            upper[cds_vector_end + 2] = len(self.poly_a_lengths) - 1
        return upper

    @property
    def original_vector(self) -> list[int]:
        return [0, 0, *self.codons.original_vector, 0, 0, 0]

    @property
    def source_construct(self) -> VaccineConstruct:
        return self.decode(self.original_vector)

    @property
    def amino_acids(self) -> str:
        return "M" + self.codons.amino_acids

    @property
    def mutable_variable_count(self) -> int:
        return sum(
            upper > lower
            for lower, upper in zip(self.lower_bounds, self.upper_bounds, strict=True)
        )

    @property
    def candidate_length_range(self) -> tuple[int, int]:
        source = self.source_construct
        five_lengths = (
            [len(variant.sequence) for variant in self.five_prime_utrs]
            if "five_prime_utr" in self.mutable_parts
            else [len(source.five_prime_utr.sequence)]
        )
        three_lengths = (
            [len(variant.sequence) for variant in self.three_prime_utrs]
            if "three_prime_utr" in self.mutable_parts
            else [len(source.three_prime_utr.sequence)]
        )
        poly_a_lengths = (
            list(self.poly_a_lengths)
            if "poly_a" in self.mutable_parts
            else [source.poly_a_length]
        )
        stop_lengths = (
            [len(stop) for stop in self.stop_codons]
            if "stop_codon" in self.mutable_parts else [len(source.stop_codon)]
        )
        fixed_length = len(source.start_codon) + len(source.cds)
        return (
            fixed_length + min(stop_lengths) + min(five_lengths) + min(three_lengths) + min(poly_a_lengths),
            fixed_length + max(stop_lengths) + max(five_lengths) + max(three_lengths) + max(poly_a_lengths),
        )

    def decode(self, vector: Iterable[int | float]) -> VaccineConstruct:
        values = [int(round(float(value))) for value in vector]
        if len(values) != self.n_var:
            raise ValueError(f"Expected {self.n_var} design variables, got {len(values)}")
        cds_vector_start = 2
        cds_vector_end = cds_vector_start + self.codons.n_var
        construct = VaccineConstruct(
            name=self.name,
            five_prime_utr=_select(values[0], self.five_prime_utrs),
            start_codon=_select(values[1], self.start_codons),
            cds=self.codons.decode(values[cds_vector_start:cds_vector_end]),
            stop_codon=_select(values[cds_vector_end], self.stop_codons),
            three_prime_utr=_select(
                values[cds_vector_end + 1], self.three_prime_utrs
            ),
            poly_a_length=_select(
                values[cds_vector_end + 2], self.poly_a_lengths
            ),
        )
        self.validate(construct)
        return construct

    def encode(self, construct: VaccineConstruct) -> list[int]:
        self.validate(construct)
        return [
            _variant_index(self.five_prime_utrs, construct.five_prime_utr),
            self.start_codons.index(construct.start_codon),
            *self.codons.encode(construct.cds),
            self.stop_codons.index(construct.stop_codon),
            _variant_index(self.three_prime_utrs, construct.three_prime_utr),
            self.poly_a_lengths.index(construct.poly_a_length),
        ]

    def validate(self, construct: VaccineConstruct) -> None:
        from morse.validators import StructuralValidator

        StructuralValidator().validate(construct, self).raise_for_errors()

    def changed_codon_count(self, construct: VaccineConstruct) -> int:
        return int(construct.start_codon != self.source_construct.start_codon) + (
            self.codons.changed_codon_count(construct.cds)
        )

    def fingerprint(self) -> str:
        return hashlib.sha256(_canonical_json(self._fingerprint_payload()).encode()).hexdigest()

    def portable_fingerprint(self) -> str:
        """Hash the complete design while excluding UTR library locations only."""
        payload = json.loads(_canonical_json(self._fingerprint_payload()))
        for key in ("five_prime_utrs", "three_prime_utrs"):
            for variant in payload[key]:
                if isinstance(variant["metadata"], dict):
                    variant["metadata"].pop("library", None)
        return hashlib.sha256(_canonical_json(payload).encode()).hexdigest()

    def _fingerprint_payload(self) -> dict[str, Any]:
        payload = {
            "schema_version": 3,
            "name": self.name,
            "five_prime_utrs": [variant.__dict__ for variant in self.five_prime_utrs],
            "start_codons": self.start_codons,
            "source_cds_body": self.codons.source_sequence,
            "stop_codons": self.stop_codons,
            "three_prime_utrs": [variant.__dict__ for variant in self.three_prime_utrs],
            "poly_a_lengths": self.poly_a_lengths,
            "allow_upstream_start_codons": self.allow_upstream_start_codons,
            "mutable_parts": sorted(self.mutable_parts),
        }
        # Keep legacy fingerprints byte-for-byte compatible, including resumes.
        if self.upstream_start_codons is not None:
            payload["upstream_start_codons"] = self.upstream_start_codons
        return payload

def construct_design_from_config(
    name: str,
    raw: dict[str, Any],
    base_dir: Path,
    mutable_parts: frozenset[str] | set[str],
    *,
    source_sequence: str | None = None,
    codon_choices: tuple[tuple[str, ...], tuple[str, ...]] | None = None,
) -> ConstructDesignSpace:
    if source_sequence is not None and ("coordinates" not in raw or "parts" in raw):
        raise ValueError("An explicit sequence requires construct.coordinates without construct.parts")
    from_coordinates = "fasta" in raw or "coordinates" in raw
    if codon_choices is not None and not from_coordinates:
        raise ValueError("Explicit codon choices require construct.fasta/coordinates")
    if from_coordinates:
        raw = _parts_from_fasta_coordinates(
            raw, base_dir, source_sequence=source_sequence, codon_choices=codon_choices
        )
    _reject_unknown_fields(
        raw,
        {"order", "allow_upstream_start_codons", "parts"},
        "construct",
    )
    if raw.get("order") != list(PART_ORDER):
        raise ValueError("construct.order must be exactly: " + ", ".join(PART_ORDER))
    parts = _mapping(raw.get("parts"), "construct.parts")
    missing = set(PART_ORDER) - set(parts)
    extra = set(parts) - set(PART_ORDER)
    if missing or extra:
        details = []
        if missing:
            details.append("missing " + ", ".join(sorted(missing)))
        if extra:
            details.append("unknown " + ", ".join(sorted(extra)))
        raise ValueError("construct.parts is invalid: " + "; ".join(details))

    unknown_mutable = set(mutable_parts) - set(PART_ORDER)
    if unknown_mutable:
        raise ValueError(
            "Unknown mutable construct parts: " + ", ".join(sorted(unknown_mutable))
        )
    selected_mutable = frozenset(mutable_parts)
    five = _parse_utr_part(parts["five_prime_utr"], "five_prime_utr", base_dir)
    three = _parse_utr_part(parts["three_prime_utr"], "three_prime_utr", base_dir)
    starts = _parse_start_part(parts["start_codon"])
    if not starts or any(len(start) != 3 or start in STOP_CODONS for start in starts):
        raise ValueError(
            "construct.parts.start_codon.source.sequences must contain "
            "non-stop three-nucleotide codons"
        )
    if len(starts) != len(set(starts)):
        raise ValueError("construct.parts.start_codon.source contains duplicates")
    source_codons = _parse_cds_part(parts["cds"], base_dir)
    if (not source_codons and not from_coordinates) or len(source_codons) % 3:
        raise ValueError("CDS body must be non-empty and divisible by three")
    if "*" in translate_dna(source_codons):
        raise ValueError("CDS body contains a stop codon")
    stops = _parse_stop_part(parts["stop_codon"])
    allowed_stops = STOP_CODONS | ({""} if from_coordinates else set())
    if not stops or any(stop not in allowed_stops for stop in stops):
        raise ValueError(
            "construct.parts.stop_codon.source.sequences must contain TAA, TAG, or TGA"
        )
    if len(stops) != len(set(stops)):
        raise ValueError("construct.parts.stop_codon.source contains duplicates")
    poly_a_lengths = _parse_poly_a_part(parts["poly_a"])
    minimum_poly_a = 0 if from_coordinates else 1
    if not poly_a_lengths or any(length < minimum_poly_a for length in poly_a_lengths):
        raise ValueError("construct.parts.poly_a source must contain positive lengths")
    if len(poly_a_lengths) != len(set(poly_a_lengths)):
        raise ValueError("construct.parts.poly_a source contains duplicate lengths")

    allow_upstream_start_codons = raw.get(
        "allow_upstream_start_codons", False
    )
    if not isinstance(allow_upstream_start_codons, bool):
        raise ValueError(
            "construct.allow_upstream_start_codons must be a boolean"
        )
    if not allow_upstream_start_codons:
        upstream_starts = ("ATG",) if codon_choices is not None else starts
        invalid_utrs = [
            variant.name
            for variant in five
            if any(start in variant.sequence for start in upstream_starts)
        ]
        if invalid_utrs:
            raise ValueError(
                "5'UTR variants contain configured upstream start codons: "
                + ", ".join(invalid_utrs)
            )

    design = ConstructDesignSpace(
        name=name,
        five_prime_utrs=five,
        start_codons=starts,
        codons=build_design_space(source_codons),
        stop_codons=stops,
        three_prime_utrs=three,
        poly_a_lengths=poly_a_lengths,
        allow_upstream_start_codons=allow_upstream_start_codons,
        mutable_parts=selected_mutable,
        upstream_start_codons=("ATG",) if codon_choices is not None else None,
    )
    design.validate(design.source_construct)
    return design


def _parts_from_fasta_coordinates(
    raw: dict[str, Any], base_dir: Path, *, source_sequence: str | None = None,
    codon_choices: tuple[tuple[str, ...], tuple[str, ...]] | None = None,
) -> dict[str, Any]:
    """Read one reference RNA and split it using explicit [start, end) ranges."""
    if "parts" in raw:
        raise ValueError("construct.fasta/coordinates cannot be combined with construct.parts")
    _reject_unknown_fields(
        raw,
        {"fasta", "coordinates", "order", "allow_upstream_start_codons"},
        "construct",
    )
    fasta = raw.get("fasta")
    if source_sequence is None and (not isinstance(fasta, str) or not fasta.strip()):
        raise ValueError("construct.fasta must be a non-empty path")
    coordinates = _mapping(raw.get("coordinates"), "construct.coordinates")
    input_order = ("five_prime_utr", "cds", "three_prime_utr", "poly_a")
    extra = sorted(str(part) for part in coordinates if part not in input_order)
    if "cds" not in coordinates or extra:
        details = []
        if "cds" not in coordinates:
            details.append("missing cds")
        if extra:
            details.append("unknown " + ", ".join(extra))
        raise ValueError("construct.coordinates is invalid: " + "; ".join(details))

    if source_sequence is None:
        path = Path(fasta).expanduser()
        if not path.is_absolute():
            path = base_dir / path
        sequence = normalize_dna(read_fasta(path).sequence)
    else:
        sequence = normalize_dna(source_sequence)
    segments: dict[str, str] = {}
    cursor = 0
    for part in input_order:
        if part not in coordinates:
            segments[part] = ""
            continue
        bounds = coordinates[part]
        context = f"construct.coordinates.{part}"
        if (
            not isinstance(bounds, list)
            or len(bounds) != 2
            or any(type(value) is not int for value in bounds)
        ):
            raise ValueError(f"{context} must be [start, end] with integer coordinates")
        start, end = bounds
        if not 0 <= start < end <= len(sequence):
            raise ValueError(
                f"{context} must satisfy 0 <= start < end <= {len(sequence)}"
            )
        if start != cursor:
            raise ValueError(
                f"{context} must start at {cursor}; parts must cover the full "
                "FASTA in order without gaps or overlaps"
            )
        segments[part] = sequence[start:end]
        cursor = end
    if cursor != len(sequence):
        raise ValueError(
            "construct.coordinates must cover the full FASTA: "
            f"{len(sequence) - cursor} trailing nucleotides are not assigned"
        )
    if set(segments["poly_a"]) - {"A"}:
        raise ValueError("construct.coordinates.poly_a must select only A nucleotides")

    coding = segments.pop("cds")
    if codon_choices is None:
        if not coding.startswith("ATG") or len(coding) % 3:
            raise ValueError(
                "construct.coordinates.cds must select a full CDS starting with "
                "ATG/AUG and a length divisible by three"
            )
        stop = coding[-3:] if coding[-3:] in STOP_CODONS else ""
        body = coding[3:-3] if stop else coding[3:]
        segments.update(start_codon="ATG", cds=body, stop_codon=stop)
    else:
        if not coding or len(coding) % 3:
            raise ValueError(
                "construct.coordinates.cds must select a non-empty CDS body "
                "with a length divisible by three"
            )
        segments.update(start_codon=codon_choices[0][0], cds=coding,
                        stop_codon=codon_choices[1][0])

    parts = {
        part: {"source": {"type": "fixed", "sequence": segments[part]}}
        for part in PART_ORDER
        if part != "poly_a"
    }
    parts["poly_a"] = {
        "source": {
            "type": "lengths",
            "nucleotide": "A",
            "values": [len(segments["poly_a"])],
        }
    }
    if codon_choices is not None:
        for part, choices in zip(("start_codon", "stop_codon"), codon_choices, strict=True):
            parts[part] = {"source": {"type": "inline", "sequences": list(choices)}}
    return {
        "order": raw.get("order", list(PART_ORDER)),
        "parts": parts,
        "allow_upstream_start_codons": raw.get("allow_upstream_start_codons", False),
    }


def _parse_cds_part(raw: Any, base_dir: Path) -> str:
    part = _mapping(raw, "construct.parts.cds")
    _reject_unknown_fields(part, {"source"}, "construct.parts.cds")

    source = _mapping(part.get("source"), "construct.parts.cds.source")
    source_type = str(source.get("type", "")).lower()
    if source_type == "fasta":
        _reject_unknown_fields(source, {"type", "path"}, "construct.parts.cds.source")
        if not source.get("path"):
            raise ValueError("construct.parts.cds.source.path is required")
        path = Path(str(source["path"]))
        if not path.is_absolute():
            path = base_dir / path
        return normalize_dna(read_fasta(path).sequence)
    if source_type == "fixed":
        _reject_unknown_fields(
            source, {"type", "sequence"}, "construct.parts.cds.source"
        )
        if source.get("sequence") is None:
            raise ValueError("construct.parts.cds.source.sequence is required")
        return normalize_dna(str(source["sequence"]))
    raise ValueError("construct.parts.cds.source.type must be fasta or fixed")


def _parse_start_part(raw: Any) -> tuple[str, ...]:
    return _parse_codon_choice_part(raw, "start_codon")


def _parse_codon_choice_part(raw: Any, label: str) -> tuple[str, ...]:
    context = f"construct.parts.{label}"
    part = _mapping(raw, context)
    _reject_unknown_fields(part, {"source"}, context)
    source = _mapping(part.get("source"), f"{context}.source")
    source_type = str(source.get("type", "")).lower()
    if source_type == "fixed":
        _reject_unknown_fields(source, {"type", "sequence"}, f"{context}.source")
        if source.get("sequence") is None:
            raise ValueError(f"{context}.source.sequence is required")
        values = [source["sequence"]]
    elif source_type == "inline":
        _reject_unknown_fields(source, {"type", "sequences"}, f"{context}.source")
        values = source.get("sequences")
        if isinstance(values, dict):
            values = list(values.values())
        if not isinstance(values, list) or not values:
            raise ValueError(f"{context}.source.sequences must be a non-empty list")
    else:
        raise ValueError(f"{context}.source.type must be fixed or inline")
    return tuple(normalize_dna(str(value)) for value in values)


def _parse_utr_part(raw: Any, label: str, base_dir: Path) -> tuple[SequenceVariant, ...]:
    context = f"construct.parts.{label}"
    part = _mapping(raw, context)
    _reject_unknown_fields(part, {"source"}, context)
    source = _mapping(part.get("source"), f"{context}.source")
    source_type = str(source.get("type", "")).lower()

    variants: list[SequenceVariant]
    if source_type == "fixed":
        _reject_unknown_fields(source, {"type", "name", "sequence"}, f"{context}.source")
        if source.get("sequence") is None:
            raise ValueError(f"{context}.source.sequence is required")
        variants = [
            SequenceVariant(
                str(source.get("name") or f"{label}_reference"),
                normalize_dna(str(source["sequence"])),
            )
        ]
    elif source_type == "inline":
        _reject_unknown_fields(source, {"type", "sequences"}, f"{context}.source")
        variants_raw = source.get("sequences")
        if isinstance(variants_raw, dict):
            variants = [
                SequenceVariant(str(name), normalize_dna(str(sequence)))
                for name, sequence in variants_raw.items()
            ]
        elif isinstance(variants_raw, list):
            variants = [
                _inline_variant(item, label, index)
                for index, item in enumerate(variants_raw, start=1)
            ]
        else:
            raise ValueError(f"{context}.source.sequences must be a mapping or list")
    elif source_type == "tsv":
        _reject_unknown_fields(
            source,
            {"type", "path", "set", "max_candidates"},
            f"{context}.source",
        )
        variants = _load_library_variants(source, label, base_dir)
    else:
        raise ValueError(f"{context}.source.type must be fixed, inline, or tsv")
    if not variants:
        raise ValueError(f"{context}.source cannot be empty")
    names = [variant.name for variant in variants]
    if len(names) != len(set(names)):
        raise ValueError(f"{label}.variants contains duplicate names")
    sequences = [variant.sequence for variant in variants]
    if len(sequences) != len(set(sequences)):
        raise ValueError(f"{label}.variants contains duplicate sequences")
    return tuple(variants)


def _inline_variant(raw: Any, label: str, index: int) -> SequenceVariant:
    if isinstance(raw, str):
        return SequenceVariant(f"{label}_{index}", normalize_dna(raw))
    mapping = _mapping(
        raw, f"construct.parts.{label}.source.sequences[{index - 1}]"
    )
    _reject_unknown_fields(
        mapping,
        {"name", "id", "sequence", "metadata"},
        f"construct.parts.{label}.source.sequences[{index - 1}]",
    )
    if mapping.get("sequence") is None:
        raise ValueError(f"inline {label} variant must define sequence")
    return _variant_from_mapping(
        {
            "name": str(mapping.get("name") or mapping.get("id") or f"{label}_{index}"),
            "sequence": mapping["sequence"],
            **(
                {"metadata": mapping["metadata"]}
                if mapping.get("metadata") is not None
                else {}
            ),
        },
        label,
    )


def _parse_stop_part(raw: Any) -> tuple[str, ...]:
    return _parse_codon_choice_part(raw, "stop_codon")


def _parse_poly_a_part(raw: Any) -> tuple[int, ...]:
    context = "construct.parts.poly_a"
    part = _mapping(raw, context)
    _reject_unknown_fields(part, {"source"}, context)
    source = _mapping(part.get("source"), f"{context}.source")
    source_type = str(source.get("type", "")).lower()
    if source_type == "lengths":
        _reject_unknown_fields(
            source, {"type", "nucleotide", "values"}, f"{context}.source"
        )
        if str(source.get("nucleotide", "A")).upper() != "A":
            raise ValueError("poly_a length source currently supports only nucleotide A")
        values = source.get("values")
        if not isinstance(values, list) or not values:
            raise ValueError(f"{context}.source.values must be a non-empty list")
        return tuple(int(value) for value in values)
    if source_type == "fixed":
        _reject_unknown_fields(source, {"type", "sequence"}, f"{context}.source")
        sequence = str(source.get("sequence", "")).upper()
        if not sequence or set(sequence) != {"A"}:
            raise ValueError("fixed poly_a sequence must contain only A")
        return (len(sequence),)
    if source_type == "inline":
        _reject_unknown_fields(source, {"type", "sequences"}, f"{context}.source")
        sequences = source.get("sequences")
        if not isinstance(sequences, list) or not sequences:
            raise ValueError("inline poly_a sequences must be a non-empty list")
        normalized = [str(sequence).upper() for sequence in sequences]
        if any(not sequence or set(sequence) != {"A"} for sequence in normalized):
            raise ValueError("inline poly_a sequences must contain only A")
        return tuple(len(sequence) for sequence in normalized)
    raise ValueError(f"{context}.source.type must be lengths, fixed, or inline")


def _load_library_variants(
    raw: Any, label: str, base_dir: Path
) -> list[SequenceVariant]:
    if not isinstance(raw, dict) or not raw.get("path"):
        raise ValueError(f"{label}.library must define path")
    path_value = str(raw["path"])
    path = Path(path_value)
    if not path.is_absolute():
        path = base_dir / path
    set_name = str(raw["set"]) if raw.get("set") is not None else None
    max_candidates = raw.get("max_candidates")
    if max_candidates is not None:
        if isinstance(max_candidates, bool):
            raise ValueError(f"{label}.library.max_candidates must be an integer")
        try:
            max_candidates = int(max_candidates)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{label}.library.max_candidates must be an integer"
            ) from exc
    records = load_part_library(path, label, set_name, max_candidates)
    return [
        SequenceVariant(
            record.id,
            record.sequence,
            {
                "library": path_value,
                "part_type": record.type,
                "part_set": record.set_name,
                **record.metadata,
            },
        )
        for record in records
    ]


def _variant_from_mapping(raw: Any, label: str) -> SequenceVariant:
    if not isinstance(raw, dict) or "name" not in raw or "sequence" not in raw:
        raise ValueError(f"{label} variant must define name and sequence")
    metadata_raw = raw.get("metadata", {})
    if not isinstance(metadata_raw, dict):
        raise ValueError(f"{label} variant metadata must be a mapping")
    metadata = {str(key): str(value) for key, value in metadata_raw.items()}
    return SequenceVariant(
        str(raw["name"]), normalize_dna(str(raw["sequence"])), metadata or None
    )


def _select(index: int, choices: tuple[Any, ...]) -> Any:
    if index < 0 or index >= len(choices):
        raise ValueError(f"Choice index {index} outside 0..{len(choices) - 1}")
    return choices[index]


def _variant_index(choices: tuple[SequenceVariant, ...], value: SequenceVariant) -> int:
    try:
        return choices.index(value)
    except ValueError as exc:
        raise ValueError(f"Unknown sequence variant {value.name!r}") from exc


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _mapping(value: Any, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"`{name}` must be a mapping")
    return value


def _reject_unknown_fields(
    raw: dict[str, Any], allowed: set[str], context: str
) -> None:
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ValueError(f"{context} contains unknown fields: " + ", ".join(unknown))
