"""Minimal JSONL scorer: transparent sequence metrics for an installation demo.

Replace score() with your inference function. These metrics do not predict
expression or immunogenicity. Model training is outside this repository.
"""
from __future__ import annotations
import json
import math
import sys


def score(sequence: str, *, target_gc: float = 0.5) -> dict[str, float]:
    dna = sequence.upper().replace("U", "T")
    if not dna or set(dna) - set("ACGT"):
        raise ValueError("Expected a nonempty A/C/G/T (or RNA) sequence")
    gc = (dna.count("G") + dna.count("C")) / len(dna)
    return {
        "cpg_per_kb": 1000 * sum(dna[i:i + 2] == "CG" for i in range(len(dna) - 1)) / len(dna),
        "gc_distance": abs(gc - target_gc),
        "gc_fraction": gc,
    }


def main() -> None:
    for line in sys.stdin:
        request = json.loads(line)
        try:
            parameters = request.get("parameters", {})
            if not isinstance(parameters, dict) or set(parameters) - {"target_gc"}:
                raise ValueError("Supported parameters: target_gc")
            target_gc = parameters.get("target_gc", 0.5)
            if isinstance(target_gc, bool) or not isinstance(target_gc, (int, float)) or not math.isfinite(target_gc) or not 0 <= target_gc <= 1:
                raise ValueError("target_gc must be a finite number between 0 and 1")
            response = {"request_id": request["request_id"], "outputs": [score(s, target_gc=target_gc) for s in request["sequences"]]}
        except Exception as exc:
            response = {"request_id": request.get("request_id"), "error": str(exc)}
        print(json.dumps(response), flush=True)


if __name__ == "__main__":
    main()
