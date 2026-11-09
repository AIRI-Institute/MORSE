"""Minimal JSONL scorer: transparent sequence metrics for an installation demo.

Replace score() with your inference function. These metrics do not predict
expression or immunogenicity. Model training is outside this repository.
"""
from __future__ import annotations
import json
import sys


def score(sequence: str) -> dict[str, float]:
    dna = sequence.upper().replace("U", "T")
    if not dna or set(dna) - set("ACGT"):
        raise ValueError("Expected a nonempty A/C/G/T (or RNA) sequence")
    gc = (dna.count("G") + dna.count("C")) / len(dna)
    return {
        "cpg_per_kb": 1000 * sum(dna[i:i + 2] == "CG" for i in range(len(dna) - 1)) / len(dna),
        "gc_distance": abs(gc - 0.5),
        "gc_fraction": gc,
    }


def main() -> None:
    for line in sys.stdin:
        request = json.loads(line)
        try:
            response = {"request_id": request["request_id"], "outputs": [score(s) for s in request["sequences"]]}
        except Exception as exc:
            response = {"request_id": request.get("request_id"), "error": str(exc)}
        print(json.dumps(response), flush=True)


if __name__ == "__main__":
    main()
