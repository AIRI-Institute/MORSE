from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Protocol

from morse.codon import translate_dna
from morse.construct import (
    STOP_CODONS,
    ConstructDesignSpace,
    VaccineConstruct,
)

STRUCTURAL_MODULES = (
    "five_prime_utr",
    "start_codon",
    "cds",
    "stop_codon",
    "three_prime_utr",
    "poly_a",
    "assembly",
    "manufacturability",
)


@dataclass(frozen=True)
class ValidationFailure:
    code: str
    module: str
    message: str

    def to_dict(self) -> dict[str, str]:
        return {"code": self.code, "module": self.module, "message": self.message}


@dataclass(frozen=True)
class StructuralValidationResult:
    failures: tuple[ValidationFailure, ...] = ()

    @property
    def valid(self) -> bool:
        return not self.failures

    @property
    def binary(self) -> int:
        return int(self.valid)

    @property
    def module_validity(self) -> dict[str, bool]:
        failed_modules = {failure.module for failure in self.failures}
        return {module: module not in failed_modules for module in STRUCTURAL_MODULES}

    def to_dict(self) -> dict[str, object]:
        return {
            "valid": self.valid,
            "binary": self.binary,
            "modules": self.module_validity,
            "failures": [failure.to_dict() for failure in self.failures],
        }

    def raise_for_errors(self) -> None:
        if self.valid:
            return
        details = "; ".join(
            f"{failure.module}.{failure.code}: {failure.message}"
            for failure in self.failures
        )
        raise ValueError(f"Structural validation failed: {details}")


class Validator(Protocol):
    name: str

    def validate(
        self, construct: VaccineConstruct, design_space: ConstructDesignSpace
    ) -> StructuralValidationResult:
        ...


class StructuralValidator:
    name = "structural"

    def validate(
        self, construct: VaccineConstruct, design_space: ConstructDesignSpace
    ) -> StructuralValidationResult:
        failures: list[ValidationFailure] = []

        def fail(code: str, module: str, message: str) -> None:
            failures.append(ValidationFailure(code, module, message))

        if construct.name != design_space.name:
            fail("name_mismatch", "assembly", "construct name differs from design space")

        if construct.five_prime_utr not in design_space.five_prime_utrs:
            fail("unknown_variant", "five_prime_utr", "5′UTR variant is not allowed")
        _check_bases(construct.five_prime_utr.sequence, "five_prime_utr", fail)
        if not design_space.allow_upstream_start_codons:
            upstream_starts = [
                start
                for start in design_space.start_codons
                if start in construct.five_prime_utr.sequence
            ]
            if upstream_starts:
                fail(
                    "upstream_start_codon",
                    "five_prime_utr",
                    "5′UTR contains configured upstream start codon(s): "
                    + ", ".join(upstream_starts),
                )

        _check_bases(construct.start_codon, "start_codon", fail)
        if len(construct.start_codon) != 3:
            fail("length", "start_codon", "start codon must contain exactly 3 nt")
        if construct.start_codon in STOP_CODONS:
            fail("is_stop", "start_codon", "start codon cannot be a stop codon")
        if construct.start_codon not in design_space.start_codons:
            fail("not_allowed", "start_codon", "start codon is not allowed")

        _check_bases(construct.cds, "cds", fail)
        if len(construct.cds) != len(design_space.codons.source_sequence):
            fail("length_mismatch", "cds", "CDS body length differs from source")
        if len(construct.cds) % 3:
            fail("frame", "cds", "CDS body length is not divisible by three")
        if _has_only_dna(construct.cds) and len(construct.cds) % 3 == 0:
            translated = translate_dna(construct.cds)
            if "*" in translated:
                fail("internal_stop", "cds", "CDS body contains a stop codon")
            if "M" + translated != design_space.amino_acids:
                fail(
                    "protein_changed",
                    "cds",
                    "start codon and CDS body change the encoded amino-acid sequence",
                )

        _check_bases(construct.stop_codon, "stop_codon", fail)
        if construct.stop_codon not in STOP_CODONS and not (
            construct.stop_codon == "" and "" in design_space.stop_codons
        ):
            fail("not_stop", "stop_codon", "terminal codon is not a stop codon")
        if construct.stop_codon not in design_space.stop_codons:
            fail("not_allowed", "stop_codon", "stop codon is not allowed")

        if construct.three_prime_utr not in design_space.three_prime_utrs:
            fail("unknown_variant", "three_prime_utr", "3′UTR variant is not allowed")
        _check_bases(construct.three_prime_utr.sequence, "three_prime_utr", fail)

        if construct.poly_a_length not in design_space.poly_a_lengths:
            fail("length_not_allowed", "poly_a", "poly(A) length is not allowed")
        if construct.segments.get("poly_a") != "A" * construct.poly_a_length:
            fail("sequence_mismatch", "poly_a", "poly(A) sequence does not match its length")

        expected_segments = {
            "five_prime_utr": construct.five_prime_utr.sequence,
            "start_codon": construct.start_codon,
            "cds": construct.cds,
            "stop_codon": construct.stop_codon,
            "three_prime_utr": construct.three_prime_utr.sequence,
            "poly_a": "A" * construct.poly_a_length,
        }
        if construct.segments != expected_segments:
            fail("segments_mismatch", "assembly", "construct segments are inconsistent")
        expected_sequence = "".join(expected_segments.values())
        if construct.sequence != expected_sequence:
            fail("sequence_mismatch", "assembly", "assembled sequence is inconsistent")
        if construct.coordinates != _expected_coordinates(expected_segments):
            fail("coordinates", "assembly", "module coordinates are not contiguous")
        _check_bases(construct.sequence, "assembly", fail)

        return StructuralValidationResult(tuple(failures))


@dataclass(frozen=True)
class ManufacturabilityPolicy:
    enabled: bool = False
    forbidden_sequences: tuple[str, ...] = ()
    include_reverse_complements: bool = True
    gc_min: float = 0.30
    gc_max: float = 0.70
    local_gc_window: int = 50
    local_gc_min: float = 0.20
    local_gc_max: float = 0.80
    max_homopolymer: int = 8
    max_tandem_repeat_copies: int = 4
    max_tandem_repeat_unit: int = 6
    max_direct_repeat: int = 20
    max_inverted_repeat: int = 12
    palindrome_min_loop: int = 3

    def __post_init__(self) -> None:
        normalized = tuple(
            dict.fromkeys(
                _normalize_motif(value)
                for value in self.forbidden_sequences
                if str(value).strip()
            )
        )
        object.__setattr__(self, "forbidden_sequences", normalized)
        for low, high, label in (
            (self.gc_min, self.gc_max, "global GC"),
            (self.local_gc_min, self.local_gc_max, "local GC"),
        ):
            if not 0.0 <= float(low) <= float(high) <= 1.0:
                raise ValueError(
                    f"manufacturability {label} bounds must satisfy 0 <= min <= max <= 1"
                )
        for name in (
            "local_gc_window",
            "max_homopolymer",
            "max_tandem_repeat_copies",
            "max_tandem_repeat_unit",
            "max_direct_repeat",
            "max_inverted_repeat",
        ):
            if int(getattr(self, name)) < 1:
                raise ValueError(f"manufacturability.{name} must be positive")
        if int(self.palindrome_min_loop) < 0:
            raise ValueError(
                "manufacturability.palindrome_min_loop must be non-negative"
            )

    def to_dict(self) -> dict[str, object]:
        return {
            "enabled": self.enabled,
            "forbidden_sequences": list(self.forbidden_sequences),
            "include_reverse_complements": self.include_reverse_complements,
            "gc_min": self.gc_min,
            "gc_max": self.gc_max,
            "local_gc_window": self.local_gc_window,
            "local_gc_min": self.local_gc_min,
            "local_gc_max": self.local_gc_max,
            "max_homopolymer": self.max_homopolymer,
            "max_tandem_repeat_copies": self.max_tandem_repeat_copies,
            "max_tandem_repeat_unit": self.max_tandem_repeat_unit,
            "max_direct_repeat": self.max_direct_repeat,
            "max_inverted_repeat": self.max_inverted_repeat,
            "palindrome_min_loop": self.palindrome_min_loop,
        }


class ManufacturabilityValidator:
    name = "manufacturability"

    def __init__(self, policy: ManufacturabilityPolicy) -> None:
        self.policy = policy

    def validate(
        self, construct: VaccineConstruct, design_space: ConstructDesignSpace
    ) -> StructuralValidationResult:
        del design_space
        if not self.policy.enabled:
            return StructuralValidationResult()
        sequence = "".join(
            construct.segments[name]
            for name in (
                "five_prime_utr",
                "start_codon",
                "cds",
                "stop_codon",
                "three_prime_utr",
            )
        )
        failures: list[ValidationFailure] = []

        def fail(code: str, message: str) -> None:
            failures.append(
                ValidationFailure(code, "manufacturability", message)
            )

        gc = _gc(sequence)
        if gc < self.policy.gc_min:
            fail("global_gc_low", f"GC {gc:.4f} is below {self.policy.gc_min:g}")
        if gc > self.policy.gc_max:
            fail("global_gc_high", f"GC {gc:.4f} is above {self.policy.gc_max:g}")

        motifs = set(self.policy.forbidden_sequences)
        if self.policy.include_reverse_complements:
            motifs.update(_reverse_complement(motif) for motif in tuple(motifs))
        for motif in sorted(motifs):
            count = _overlapping_count(sequence, motif)
            if count:
                fail(
                    "forbidden_sequence",
                    f"forbidden sequence {motif} occurs {count} time(s)",
                )

        window = self.policy.local_gc_window
        if len(sequence) >= window:
            for start in range(len(sequence) - window + 1):
                value = _gc(sequence[start : start + window])
                if value < self.policy.local_gc_min:
                    fail(
                        "local_gc_low",
                        f"GC window at {start} is {value:.4f}, below {self.policy.local_gc_min:g}",
                    )
                    break
                if value > self.policy.local_gc_max:
                    fail(
                        "local_gc_high",
                        f"GC window at {start} is {value:.4f}, above {self.policy.local_gc_max:g}",
                    )
                    break

        longest_run = _longest_homopolymer(sequence)
        if longest_run > self.policy.max_homopolymer:
            fail(
                "homopolymer",
                f"longest homopolymer is {longest_run}, maximum is {self.policy.max_homopolymer}",
            )
        tandem = _largest_tandem_copy_count(
            sequence, self.policy.max_tandem_repeat_unit
        )
        if tandem > self.policy.max_tandem_repeat_copies:
            fail(
                "tandem_repeat",
                f"tandem repeat has {tandem} copies, maximum is {self.policy.max_tandem_repeat_copies}",
            )
        if _has_repeated_substring(sequence, self.policy.max_direct_repeat + 1):
            fail(
                "direct_repeat",
                f"direct repeat exceeds {self.policy.max_direct_repeat} nt",
            )
        if _has_inverted_repeat(
            sequence,
            self.policy.max_inverted_repeat + 1,
            self.policy.palindrome_min_loop,
        ):
            fail(
                "inverted_repeat",
                f"palindromic stem exceeds {self.policy.max_inverted_repeat} nt",
            )
        return StructuralValidationResult(tuple(failures))


class CompositeValidator:
    name = "composite"

    def __init__(self, *validators: Validator) -> None:
        self.validators = validators

    def validate(
        self, construct: VaccineConstruct, design_space: ConstructDesignSpace
    ) -> StructuralValidationResult:
        failures = tuple(
            failure
            for validator in self.validators
            for failure in validator.validate(construct, design_space).failures
        )
        return StructuralValidationResult(failures)


def validator_for_policy(policy: ManufacturabilityPolicy) -> Validator:
    if not policy.enabled:
        return StructuralValidator()
    return CompositeValidator(
        StructuralValidator(), ManufacturabilityValidator(policy)
    )

def _check_bases(
    sequence: str, module: str, fail: Callable[[str, str, str], None]
) -> None:
    invalid = sorted(set(sequence) - {"A", "C", "G", "T"})
    if invalid:
        fail("invalid_bases", module, f"invalid bases: {', '.join(invalid)}")


def _has_only_dna(sequence: str) -> bool:
    return not (set(sequence) - {"A", "C", "G", "T"})


def _expected_coordinates(segments: dict[str, str]) -> dict[str, tuple[int, int]]:
    cursor = 0
    coordinates: dict[str, tuple[int, int]] = {}
    for name, sequence in segments.items():
        coordinates[name] = (cursor, cursor + len(sequence))
        cursor += len(sequence)
    return coordinates


def _normalize_motif(value: str) -> str:
    motif = "".join(str(value).upper().split()).replace("U", "T")
    invalid = sorted(set(motif) - {"A", "C", "G", "T"})
    if not motif or invalid:
        raise ValueError(
            "manufacturability forbidden sequences must contain only A/C/G/T/U"
        )
    return motif


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _gc(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / max(len(sequence), 1)


def _overlapping_count(sequence: str, motif: str) -> int:
    return sum(
        sequence.startswith(motif, index)
        for index in range(max(0, len(sequence) - len(motif) + 1))
    )


def _longest_homopolymer(sequence: str) -> int:
    longest = current = 0
    previous = ""
    for base in sequence:
        current = current + 1 if base == previous else 1
        longest = max(longest, current)
        previous = base
    return longest


def _largest_tandem_copy_count(sequence: str, max_unit: int) -> int:
    best = 1
    for unit_length in range(1, max_unit + 1):
        for start in range(max(0, len(sequence) - 2 * unit_length + 1)):
            unit = sequence[start : start + unit_length]
            copies = 1
            cursor = start + unit_length
            while sequence[cursor : cursor + unit_length] == unit:
                copies += 1
                cursor += unit_length
            best = max(best, copies)
    return best


def _has_repeated_substring(sequence: str, length: int) -> bool:
    seen: set[str] = set()
    for start in range(max(0, len(sequence) - length + 1)):
        value = sequence[start : start + length]
        if value in seen:
            return True
        seen.add(value)
    return False


def _has_inverted_repeat(sequence: str, stem: int, min_loop: int) -> bool:
    earliest: dict[str, int] = {}
    for start in range(max(0, len(sequence) - stem + 1)):
        value = sequence[start : start + stem]
        partner_start = earliest.get(_reverse_complement(value))
        if (
            partner_start is not None
            and partner_start + stem + min_loop <= start
        ):
            return True
        earliest.setdefault(value, start)
    return False
