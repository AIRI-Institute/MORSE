from __future__ import annotations

import json
import os
import sys
from typing import Any

WINDOW = int(os.environ.get("SEQUENCE_QUALITY_WINDOW", "60"))
STRIDE = int(os.environ.get("SEQUENCE_QUALITY_STRIDE", "10"))
if WINDOW < 1 or STRIDE < 1:
    raise ValueError("SEQUENCE_QUALITY_WINDOW and SEQUENCE_QUALITY_STRIDE must be positive")


def main() -> int:
    for line in sys.stdin:
        if not line.strip():
            continue
        try:
            request = json.loads(line)
            unknown = sorted(
                set(request) - {"request_id", "sequences", "sequence_contexts"}
            )
            if unknown:
                raise ValueError("Unknown request fields: " + ", ".join(unknown))
            sequences = request.get("sequences")
            if not isinstance(sequences, list) or not all(
                isinstance(sequence, str) for sequence in sequences
            ):
                raise ValueError("Request must define a list of sequences")
            sequence_contexts = request.get("sequence_contexts")
            if (
                not isinstance(sequence_contexts, list)
                or len(sequence_contexts) != len(sequences)
                or not all(isinstance(item, dict) for item in sequence_contexts)
            ):
                raise ValueError(
                    "Request must define one sequence context per sequence"
                )
            response = {
                "request_id": request.get("request_id"),
                "outputs": [score(sequence) for sequence in sequences],
            }
        except Exception as exc:
            response = {"request_id": _request_id(line), "error": str(exc)}
        print(json.dumps(response, separators=(",", ":")), flush=True)
    return 0


def score(raw_sequence: str) -> dict[str, Any]:
    sequence = _normalize(raw_sequence)
    windows = _windows(sequence)
    gc_windows = [_gc(window) for window in windows if window]
    return {
        "non_polya_max_homopolymer": _longest_homopolymer(sequence),
        "body_gc_window_min": min(gc_windows, default=0.0),
        "body_gc_window_max": max(gc_windows, default=0.0),
        "cpg_per_kb": 1000.0
        * _count_overlapping(sequence, "CG")
        / max(len(sequence), 1),
        "upa_per_kb": 1000.0
        * _count_overlapping(sequence, "TA")
        / max(len(sequence), 1),
        "repeated_12mer_count": _repeated_kmer_count(sequence, 12),
    }


def _normalize(sequence: str) -> str:
    normalized = "".join(sequence.upper().split()).replace("U", "T")
    invalid = sorted(set(normalized) - {"A", "C", "G", "T"})
    if invalid:
        raise ValueError(f"Unsupported bases: {', '.join(invalid)}")
    return normalized


def _gc(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / max(len(sequence), 1)


def _windows(sequence: str) -> list[str]:
    if len(sequence) <= WINDOW:
        return [sequence]
    starts = list(range(0, len(sequence) - WINDOW + 1, STRIDE))
    starts.append(len(sequence) - WINDOW)
    return [sequence[start : start + WINDOW] for start in sorted(set(starts))]


def _longest_homopolymer(sequence: str) -> int:
    longest = current = 0
    previous = ""
    for base in sequence:
        current = current + 1 if base == previous else 1
        longest = max(longest, current)
        previous = base
    return longest


def _count_overlapping(sequence: str, motif: str) -> int:
    return sum(
        sequence[index : index + len(motif)] == motif
        for index in range(max(len(sequence) - len(motif) + 1, 0))
    )


def _repeated_kmer_count(sequence: str, length: int) -> int:
    counts: dict[str, int] = {}
    for index in range(max(len(sequence) - length + 1, 0)):
        kmer = sequence[index : index + length]
        counts[kmer] = counts.get(kmer, 0) + 1
    return sum(count - 1 for count in counts.values() if count > 1)


def _request_id(line: str) -> Any:
    try:
        return json.loads(line).get("request_id")
    except Exception:
        return None


if __name__ == "__main__":
    raise SystemExit(main())
