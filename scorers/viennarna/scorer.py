from __future__ import annotations

import json
import os
import sys
from typing import Any

ENABLE_PARTITION_FUNCTION = os.environ.get("VIENNARNA_ENABLE_PF", "1").lower() not in {
    "0",
    "false",
    "no",
    "off",
}
MAX_BP_SPAN = int(os.environ.get("VIENNARNA_MAX_BP_SPAN", "0") or 0)


def main() -> int:
    for line in sys.stdin:
        if not line.strip():
            continue
        try:
            request = json.loads(line)
            response = score_batch(request)
        except Exception as exc:
            response = _error_response(line, str(exc))
        print(json.dumps(response, sort_keys=True), flush=True)
    return 0


def score_batch(request: dict[str, Any]) -> dict[str, Any]:
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
    return {
        "request_id": request.get("request_id"),
        "outputs": [score_sequence(sequence) for sequence in sequences],
    }


def score_sequence(sequence: str) -> dict[str, Any]:
    dna_sequence = "".join(sequence.upper().split()).replace("U", "T")
    rna_sequence = dna_sequence.replace("T", "U")
    if not rna_sequence:
        raise ValueError("RNA sequence must not be empty")
    invalid = sorted(set(rna_sequence) - {"A", "C", "G", "U"})
    if invalid:
        raise ValueError(f"Invalid RNA bases: {', '.join(invalid)}")

    import RNA  # type: ignore

    model_details = RNA.md()
    if MAX_BP_SPAN > 0:
        model_details.max_bp_span = MAX_BP_SPAN
    fold_compound = RNA.fold_compound(rna_sequence, model_details)
    mfe_structure, mfe = fold_compound.mfe()
    gc_fraction = (dna_sequence.count("G") + dna_sequence.count("C")) / len(dna_sequence)

    outputs = {
        "mfe": float(mfe),
        "mfe_per_nt": float(mfe) / max(len(rna_sequence), 1),
        "gc_fraction": gc_fraction,
        "mfe_structure": mfe_structure,
    }
    if ENABLE_PARTITION_FUNCTION:
        fold_compound.exp_params_rescale(mfe)
        _, ensemble_free_energy = fold_compound.pf()
        outputs.update(
            {
                "ensemble_diversity": float(fold_compound.mean_bp_distance()),
                "ensemble_free_energy": float(ensemble_free_energy),
            }
        )

    return outputs


def _error_response(line: str, message: str) -> dict[str, Any]:
    try:
        raw = json.loads(line)
    except Exception:
        return {"request_id": None, "error": message}
    return {"request_id": raw.get("request_id"), "error": message}


if __name__ == "__main__":
    raise SystemExit(main())
