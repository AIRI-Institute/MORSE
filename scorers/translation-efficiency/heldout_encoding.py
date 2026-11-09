"""Frozen common-holdout TE input encoding (including the terminal stop)."""
from __future__ import annotations
import numpy as np


def prepare_tracks(sequence: str, context: dict) -> np.ndarray:
    sequence = sequence.upper().replace("U", "T")
    if not sequence or set(sequence) - set("ACGT"):
        raise ValueError("TE requires an unambiguous transcript body")
    order = ("five_prime_utr", "start_codon", "cds", "stop_codon", "three_prime_utr")
    found, cursor, previous = {}, 0, -1
    for span in context.get("part_ranges", []):
        part = span["part"]
        if part not in order or part in found or order.index(part) <= previous:
            raise ValueError("Expected transcript-body parts in biological order, without poly(A)")
        start, end = int(span["start"]), int(span["end"])
        if start != cursor or end <= start or end > len(sequence) or span.get("part_start", 0) != 0:
            raise ValueError("Part ranges must cover the complete routed sequence")
        found[part] = (start, end)
        cursor, previous = end, order.index(part)
    if cursor != len(sequence) or any(k not in found for k in ("start_codon", "cds", "stop_codon")):
        raise ValueError("Explicit complete coding context is required")
    first, last = found["start_codon"][0], found["stop_codon"][1]
    if (found["start_codon"][1] - first != 3 or last - found["stop_codon"][0] != 3
            or sequence[first:first + 3] != "ATG" or sequence[last - 3:last] not in ("TAA", "TAG", "TGA")
            or (last - first) % 3):
        raise ValueError("Expected AUG and a canonical terminal stop in frame")
    if first > 1381 or len(sequence) - first > 11937:
        raise ValueError("Transcript exceeds the frozen adapter's input limits")
    result = np.zeros((len(sequence), 6), dtype=np.float32)
    columns = {base: i for i, base in enumerate("ACGT")}
    result[np.arange(len(sequence)), [columns[base] for base in sequence]] = 1.0
    result[first:last, 4] = 1.0
    return result
