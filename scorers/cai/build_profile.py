#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterator

if __package__:
    from .model import (
        DNA_CODON_TABLE,
        STOP_CODONS,
        normalize_dna,
        profile_from_mapping,
        split_codons,
    )
else:
    from model import (
        DNA_CODON_TABLE,
        STOP_CODONS,
        normalize_dna,
        profile_from_mapping,
        split_codons,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a versioned CAI reference profile from CDS FASTA files."
    )
    parser.add_argument("fasta", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--profile-id", required=True)
    parser.add_argument("--taxon", required=True)
    parser.add_argument("--taxon-id", type=int)
    parser.add_argument("--reference-kind", default="highly_expressed_gene_set")
    parser.add_argument("--source-name", required=True)
    parser.add_argument("--source-version", required=True)
    parser.add_argument("--source-url")
    parser.add_argument("--zero-count-replacement", type=float, default=0.5)
    args = parser.parse_args()

    profile = build_profile(
        paths=args.fasta,
        profile_id=args.profile_id,
        taxon=args.taxon,
        taxon_id=args.taxon_id,
        reference_kind=args.reference_kind,
        source_name=args.source_name,
        source_version=args.source_version,
        source_url=args.source_url,
        zero_count_replacement=args.zero_count_replacement,
    )
    profile_from_mapping(profile)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(profile, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {args.output}: {profile['source']['cds_count']} CDS, "
        f"{profile['source']['codon_count']} codons"
    )
    return 0


def build_profile(
    paths: list[Path],
    profile_id: str,
    taxon: str,
    taxon_id: int | None,
    reference_kind: str,
    source_name: str,
    source_version: str,
    source_url: str | None,
    zero_count_replacement: float = 0.5,
) -> dict[str, object]:
    if zero_count_replacement <= 0:
        raise ValueError("zero_count_replacement must be greater than zero")
    counts = {codon: 0 for codon in DNA_CODON_TABLE}
    cds_count = 0
    codon_count = 0
    for path in paths:
        for record_id, raw_sequence in read_fasta(path):
            codons = _reference_codons(record_id, raw_sequence)
            for codon in codons:
                counts[codon] += 1
            cds_count += 1
            codon_count += len(codons)
    if cds_count == 0:
        raise ValueError("Reference FASTA set does not contain any CDS records")

    source: dict[str, object] = {
        "name": source_name,
        "version": source_version,
        "files": [str(path) for path in paths],
        "cds_count": cds_count,
        "codon_count": codon_count,
    }
    if source_url:
        source["url"] = source_url
    result: dict[str, object] = {
        "schema_version": 1,
        "profile_id": profile_id,
        "taxon": taxon,
        "genetic_code": 1,
        "reference_kind": reference_kind,
        "zero_count_replacement": zero_count_replacement,
        "source": source,
        "codon_counts": counts,
    }
    if taxon_id is not None:
        result["taxon_id"] = taxon_id
    return result


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    record_id: str | None = None
    chunks: list[str] = []
    for line_number, raw_line in enumerate(path.read_text().splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if record_id is not None:
                yield record_id, "".join(chunks)
            record_id = line[1:].split(maxsplit=1)[0]
            if not record_id:
                raise ValueError(f"Empty FASTA identifier in {path}:{line_number}")
            chunks = []
        elif record_id is None:
            raise ValueError(f"Sequence appears before FASTA header in {path}:{line_number}")
        else:
            chunks.append(line)
    if record_id is not None:
        yield record_id, "".join(chunks)


def _reference_codons(record_id: str, sequence: str) -> list[str]:
    cds = normalize_dna(sequence)
    codons = split_codons(cds)
    if not codons:
        raise ValueError(f"Reference CDS {record_id!r} is empty")
    if codons[0] != "ATG":
        raise ValueError(f"Reference CDS {record_id!r} does not start with ATG")
    if codons[-1] in STOP_CODONS:
        coding_codons = codons[:-1]
    else:
        coding_codons = codons
    internal_stops = [codon for codon in coding_codons if codon in STOP_CODONS]
    if internal_stops:
        raise ValueError(f"Reference CDS {record_id!r} contains an internal stop codon")
    return codons


if __name__ == "__main__":
    raise SystemExit(main())
