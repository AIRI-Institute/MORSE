from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

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

SYNONYMOUS_CODONS: dict[str, list[str]] = {}
for _codon, _aa in DNA_CODON_TABLE.items():
    SYNONYMOUS_CODONS.setdefault(_aa, []).append(_codon)
for _aa in SYNONYMOUS_CODONS:
    SYNONYMOUS_CODONS[_aa] = sorted(SYNONYMOUS_CODONS[_aa])


@dataclass(frozen=True)
class CodonDesignSpace:
    source_sequence: str
    amino_acids: str
    choices_by_position: tuple[tuple[str, ...], ...]
    original_indices: tuple[int, ...]

    @property
    def n_var(self) -> int:
        return len(self.choices_by_position)

    @property
    def lower_bounds(self) -> list[int]:
        return [0] * self.n_var

    @property
    def upper_bounds(self) -> list[int]:
        return [len(choices) - 1 for choices in self.choices_by_position]

    @property
    def original_vector(self) -> list[int]:
        return list(self.original_indices)

    def decode(self, vector: Iterable[int | float]) -> str:
        codons: list[str] = []
        for raw_value, choices in zip(vector, self.choices_by_position, strict=True):
            index = int(round(float(raw_value)))
            if index < 0 or index >= len(choices):
                raise ValueError(f"Codon choice index {index} outside 0..{len(choices) - 1}")
            codons.append(choices[index])
        return "".join(codons)

    def encode(self, sequence: str) -> list[int]:
        sequence = normalize_dna(sequence)
        self.validate_synonymous(sequence)
        vector: list[int] = []
        for codon, choices in zip(split_codons(sequence), self.choices_by_position, strict=True):
            try:
                vector.append(choices.index(codon))
            except ValueError as exc:
                raise ValueError(f"Codon {codon} is not available in the design space") from exc
        return vector

    def changed_codon_count(self, sequence: str) -> int:
        return sum(
            1
            for i in range(0, len(self.source_sequence), 3)
            if self.source_sequence[i : i + 3] != sequence[i : i + 3]
        )

    def validate_synonymous(self, sequence: str) -> None:
        if len(sequence) != len(self.source_sequence):
            raise ValueError("Candidate length differs from source CDS length")
        if translate_dna(sequence) != self.amino_acids:
            raise ValueError("Candidate changes translated amino-acid sequence")


def normalize_dna(sequence: str) -> str:
    normalized = "".join(sequence.upper().split()).replace("U", "T")
    invalid = sorted(set(normalized) - {"A", "C", "G", "T"})
    if invalid:
        raise ValueError(f"Invalid DNA/RNA bases: {', '.join(invalid)}")
    return normalized


def split_codons(sequence: str) -> list[str]:
    sequence = normalize_dna(sequence)
    if len(sequence) % 3 != 0:
        raise ValueError(f"CDS length must be divisible by 3; got {len(sequence)}")
    return [sequence[i : i + 3] for i in range(0, len(sequence), 3)]


def translate_dna(sequence: str) -> str:
    amino_acids: list[str] = []
    for codon in split_codons(sequence):
        try:
            amino_acids.append(DNA_CODON_TABLE[codon])
        except KeyError as exc:
            raise ValueError(f"Unknown codon {codon}") from exc
    return "".join(amino_acids)


def build_design_space(sequence: str, preserve_stop: bool = True) -> CodonDesignSpace:
    sequence = normalize_dna(sequence)
    codons = split_codons(sequence)
    amino_acids = translate_dna(sequence)
    choices_by_position: list[tuple[str, ...]] = []
    original_indices: list[int] = []

    for position, (codon, amino_acid) in enumerate(zip(codons, amino_acids, strict=True)):
        choices = SYNONYMOUS_CODONS[amino_acid]
        if preserve_stop and amino_acid == "*":
            choices = [codon]
        if position == 0 and amino_acid == "M":
            choices = [codon]
        choices_tuple = tuple(choices)
        choices_by_position.append(choices_tuple)
        original_indices.append(choices_tuple.index(codon))

    return CodonDesignSpace(
        source_sequence=sequence,
        amino_acids=amino_acids,
        choices_by_position=tuple(choices_by_position),
        original_indices=tuple(original_indices),
    )


def gc_fraction(sequence: str) -> float:
    sequence = normalize_dna(sequence)
    if not sequence:
        return 0.0
    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence)
