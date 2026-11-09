#!/usr/bin/env python3
"""Build a compact UTR parts library from GENCODE annotation and transcripts."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO
from urllib.parse import unquote

import numpy as np

GENCODE_ROOT = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
GENCODE_PAGE = "https://www.gencodegenes.org/human/"
from asset_paths import asset_directory

ASSET_DIR = asset_directory()
OUTPUT_COLUMNS = (
    "id",
    "type",
    "sequence",
    "set",
    "source",
    "release",
    "species",
    "gene",
    "gene_id",
    "transcript",
    "length_nt",
    "gc_fraction",
    "evidence_level",
    "source_url",
    "sequence_hash",
    "selection_rank",
    "notes",
)


@dataclass
class TranscriptAnnotation:
    transcript: str
    gene: str
    gene_id: str
    utr_length: int = 0
    segment_count: int = 0


@dataclass(frozen=True)
class UtrCandidate:
    transcript: str
    gene: str
    gene_id: str
    sequence: str
    segment_count: int


def build_library(
    annotation_path: Path,
    transcripts_path: Path,
    output_path: Path,
    *,
    release: int,
    utr_type: str,
    set_name: str,
    min_length: int,
    max_length: int,
    max_candidates: int | None,
    mane_only: bool = True,
) -> dict[str, int]:
    if utr_type not in {"three_prime_utr", "five_prime_utr"}:
        raise ValueError("utr_type must be three_prime_utr or five_prime_utr")
    if min_length < 1 or max_length < min_length:
        raise ValueError("Expected 1 <= min_length <= max_length")
    if max_candidates is not None and max_candidates < 1:
        raise ValueError("max_candidates must be >= 1")

    annotations = read_annotations(
        annotation_path,
        feature=utr_type.replace("_utr", "_UTR"),
        min_length=min_length,
        max_length=max_length,
        mane_only=mane_only,
    )
    if not annotations:
        raise ValueError("No transcripts match the annotation filters")
    extracted = extract_candidates(transcripts_path, annotations, utr_type)
    if not extracted:
        raise ValueError("No annotated transcripts were found in the transcript FASTA")
    unique = deduplicate(extracted)
    selected = (
        unique
        if max_candidates is None
        else select_representatives(unique, max_candidates)
    )
    write_library(output_path, selected, release, utr_type, set_name, mane_only)
    counts = {
        "annotated_transcripts": len(annotations),
        "extracted_transcripts": len(extracted),
        "unique_sequences": len(unique),
        "selected_sequences": len(selected),
    }
    write_manifest(
        output_path.with_suffix(".manifest.json"),
        annotation_path,
        transcripts_path,
        release=release,
        utr_type=utr_type,
        set_name=set_name,
        min_length=min_length,
        max_length=max_length,
        max_candidates=max_candidates,
        mane_only=mane_only,
        counts=counts,
    )
    return counts


def read_annotations(
    path: Path,
    *,
    feature: str,
    min_length: int,
    max_length: int,
    mane_only: bool,
) -> dict[str, TranscriptAnnotation]:
    annotations: dict[str, TranscriptAnnotation] = {}
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != feature:
                continue
            attributes = _parse_attributes(fields[8])
            tags = set(attributes.get("tag", "").split(","))
            if attributes.get("transcript_type") != "protein_coding":
                continue
            if mane_only and "MANE_Select" not in tags:
                continue
            transcript = _accession(attributes.get("transcript_id", ""))
            if not transcript:
                continue
            record = annotations.setdefault(
                transcript,
                TranscriptAnnotation(
                    transcript=transcript,
                    gene=attributes.get("gene_name", ""),
                    gene_id=_accession(attributes.get("gene_id", "")),
                ),
            )
            record.utr_length += int(fields[4]) - int(fields[3]) + 1
            record.segment_count += 1
    return {
        transcript: record
        for transcript, record in annotations.items()
        if min_length <= record.utr_length <= max_length
    }


def extract_candidates(
    path: Path,
    annotations: dict[str, TranscriptAnnotation],
    utr_type: str,
) -> list[UtrCandidate]:
    candidates: list[UtrCandidate] = []
    for transcript, raw_sequence in _read_fasta(path):
        accession = _accession(transcript.split("|", 1)[0])
        annotation = annotations.get(accession)
        if annotation is None:
            continue
        sequence = normalize_dna(raw_sequence)
        length = annotation.utr_length
        if length > len(sequence):
            raise ValueError(
                f"Annotated UTR length exceeds transcript length for {accession}"
            )
        utr = sequence[-length:] if utr_type == "three_prime_utr" else sequence[:length]
        candidates.append(
            UtrCandidate(
                transcript=accession,
                gene=annotation.gene,
                gene_id=annotation.gene_id,
                sequence=utr,
                segment_count=annotation.segment_count,
            )
        )
    return sorted(candidates, key=lambda item: item.transcript)


def deduplicate(candidates: list[UtrCandidate]) -> list[UtrCandidate]:
    unique: dict[str, UtrCandidate] = {}
    for candidate in candidates:
        unique.setdefault(candidate.sequence, candidate)
    return list(unique.values())


def select_representatives(
    candidates: list[UtrCandidate], max_candidates: int
) -> list[UtrCandidate]:
    if len(candidates) <= max_candidates:
        return candidates
    features = np.vstack([_sequence_features(item.sequence) for item in candidates])
    scale = features.std(axis=0)
    scale[scale == 0] = 1.0
    features = (features - features.mean(axis=0)) / scale

    centroid = features.mean(axis=0)
    selected = [int(np.argmin(np.sum((features - centroid) ** 2, axis=1)))]
    minimum_distance = np.sum((features - features[selected[0]]) ** 2, axis=1)
    minimum_distance[selected[0]] = -1.0
    while len(selected) < max_candidates:
        index = int(np.argmax(minimum_distance))
        selected.append(index)
        distance = np.sum((features - features[index]) ** 2, axis=1)
        minimum_distance = np.minimum(minimum_distance, distance)
        minimum_distance[selected] = -1.0
    return [candidates[index] for index in selected]


def write_library(
    path: Path,
    candidates: list[UtrCandidate],
    release: int,
    utr_type: str,
    set_name: str,
    mane_only: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for rank, candidate in enumerate(candidates, start=1):
            writer.writerow(
                {
                    "id": f"GENCODE_{release}_{candidate.transcript}",
                    "type": utr_type,
                    "sequence": candidate.sequence,
                    "set": set_name,
                    "source": "GENCODE",
                    "release": release,
                    "species": "Homo sapiens",
                    "gene": candidate.gene,
                    "gene_id": candidate.gene_id,
                    "transcript": candidate.transcript,
                    "length_nt": len(candidate.sequence),
                    "gc_fraction": f"{_gc_fraction(candidate.sequence):.6f}",
                    "evidence_level": "MANE_Select" if mane_only else "GENCODE",
                    "source_url": GENCODE_PAGE,
                    "sequence_hash": sequence_hash(candidate.sequence),
                    "selection_rank": rank,
                    "notes": f"full_native; utr_segments={candidate.segment_count}",
                }
            )


def write_manifest(
    path: Path,
    annotation_path: Path,
    transcripts_path: Path,
    **metadata: object,
) -> None:
    payload = {
        "schema_version": 1,
        **metadata,
        "inputs": {
            "annotation": {
                "file": annotation_path.name,
                "sha256": _file_hash(annotation_path),
            },
            "transcripts": {
                "file": transcripts_path.name,
                "sha256": _file_hash(transcripts_path),
            },
        },
        "source_urls": _source_urls(int(metadata["release"])),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _sequence_features(sequence: str) -> np.ndarray:
    length = len(sequence)
    features = [np.log1p(length)]
    features.extend(sequence.count(base) / length for base in "ACGT")
    bins = np.zeros(32, dtype=float)
    for index in range(max(0, length - 3)):
        kmer = sequence[index : index + 4]
        value = 0
        for base in kmer:
            value = value * 4 + "ACGT".index(base)
        bins[value % len(bins)] += 1
    if bins.sum():
        bins /= bins.sum()
    return np.concatenate((np.asarray(features), bins))


def _gc_fraction(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def normalize_dna(sequence: str) -> str:
    normalized = "".join(sequence.upper().split()).replace("U", "T")
    invalid = sorted(set(normalized) - {"A", "C", "G", "T"})
    if invalid:
        raise ValueError(f"Invalid DNA/RNA bases: {', '.join(invalid)}")
    return normalized


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _parse_attributes(raw: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in raw.split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        value = unquote(value)
        if key in attributes:
            attributes[key] += "," + value
        else:
            attributes[key] = value
    return attributes


def _read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    with _open_text(path) as handle:
        header = ""
        sequence: list[str] = []
        for line in handle:
            if line.startswith(">"):
                if header:
                    yield header, "".join(sequence)
                header = line[1:].strip()
                sequence = []
            else:
                sequence.append(line.strip())
        if header:
            yield header, "".join(sequence)


def _open_text(path: Path) -> TextIO:
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def _accession(value: str) -> str:
    return value.strip().split(".", 1)[0]


def _file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_urls(release: int) -> dict[str, str]:
    root = f"{GENCODE_ROOT}/release_{release}"
    prefix = f"gencode.v{release}"
    return {
        "annotation": f"{root}/{prefix}.annotation.gff3.gz",
        "transcripts": f"{root}/{prefix}.transcripts.fa.gz",
    }


def _download(url: str, path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    with urllib.request.urlopen(url) as response, temporary.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", type=Path)
    parser.add_argument("--transcripts", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--release", type=int, default=50)
    parser.add_argument(
        "--utr-type",
        choices=("three_prime_utr", "five_prime_utr"),
        default="three_prime_utr",
    )
    parser.add_argument("--set-name")
    parser.add_argument("--min-length", type=int, default=30)
    parser.add_argument("--max-length", type=int, default=500)
    parser.add_argument(
        "--max-candidates",
        type=int,
        default=None,
        help="Select a diverse subset of this size; omit to retain the full library.",
    )
    parser.add_argument("--all-transcripts", action="store_true")
    parser.add_argument("--download", action="store_true")
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path.home() / ".cache" / "morse" / "gencode",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    default_raw = ASSET_DIR / "data" / "raw" / f"gencode_v{args.release}"
    annotation = args.annotation or default_raw / f"gencode.v{args.release}.annotation.gff3.gz"
    transcripts = args.transcripts or default_raw / f"gencode.v{args.release}.transcripts.fa.gz"
    urls = _source_urls(args.release)
    if args.download:
        annotation = annotation or args.cache_dir / Path(urls["annotation"]).name
        transcripts = transcripts or args.cache_dir / Path(urls["transcripts"]).name
        _download(urls["annotation"], annotation)
        _download(urls["transcripts"], transcripts)
    if not annotation.is_file() or not transcripts.is_file():
        raise SystemExit("Run prepare_assets.py, provide --annotation/--transcripts, or use --download")
    set_name = args.set_name or (
        f"gencode_v{args.release}_{'all' if args.all_transcripts else 'mane'}_"
        f"{args.utr_type}_{args.min_length}_{args.max_length}nt"
    )
    output = args.output or (
        ASSET_DIR
        / "data"
        / "libraries"
        / f"gencode_v{args.release}_{'mane' if not args.all_transcripts else 'all'}_"
        f"{args.utr_type}_{args.max_candidates or 'all'}.tsv"
    )
    counts = build_library(
        annotation,
        transcripts,
        output,
        release=args.release,
        utr_type=args.utr_type,
        set_name=set_name,
        min_length=args.min_length,
        max_length=args.max_length,
        max_candidates=args.max_candidates,
        mane_only=not args.all_transcripts,
    )
    print(json.dumps(counts, sort_keys=True))


if __name__ == "__main__":
    main()
