from __future__ import annotations

import json
import os
import sys
from typing import Any

ENABLE_PARTITION_FUNCTION = os.environ.get("VIENNARNA_ENABLE_PF", "0").lower() not in {
    "0",
    "false",
    "no",
    "off",
}
MAX_BP_SPAN = int(os.environ.get("VIENNARNA_MAX_BP_SPAN", "300") or 300)


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
        set(request) - {"request_id", "sequences", "sequence_contexts", "parameters"}
    )
    if unknown:
        raise ValueError("Unknown request fields: " + ", ".join(unknown))
    sequences = request.get("sequences")
    if not isinstance(sequences, list) or not all(
        isinstance(sequence, str) for sequence in sequences
    ):
        raise ValueError("Request must define a list of sequences")
    if "parameters" in request:
        if not isinstance(request["parameters"], dict):
            raise ValueError("parameters must be a mapping")
        parameters = _resolve_parameters(request["parameters"])
        outputs = [score_sequence(sequence, parameters=parameters) for sequence in sequences]
    else:
        outputs = [score_sequence(sequence) for sequence in sequences]
    return {
        "request_id": request.get("request_id"),
        "outputs": outputs,
    }


def _resolve_parameters(parameters: dict[str, Any] | None = None) -> dict[str, Any]:
    if parameters is None:
        parameters = {}
    if not isinstance(parameters, dict):
        raise ValueError("parameters must be a mapping")
    unknown = set(parameters) - {"max_bp_span", "compute_ensemble"}
    if unknown:
        raise ValueError("Unknown MFE parameters: " + ", ".join(sorted(map(str, unknown))))
    span = parameters.get("max_bp_span", MAX_BP_SPAN)
    ensemble = parameters.get("compute_ensemble", ENABLE_PARTITION_FUNCTION)
    if type(span) is not int or not 0 <= span <= 100000:
        raise ValueError("max_bp_span must be an integer between 0 and 100000 (nt)")
    if type(ensemble) is not bool:
        raise ValueError("compute_ensemble must be true or false")
    return {"max_bp_span": span, "compute_ensemble": ensemble}


def score_sequence(sequence: str, parameters: dict[str, Any] | None = None) -> dict[str, Any]:
    dna_sequence = "".join(sequence.upper().split()).replace("U", "T")
    rna_sequence = dna_sequence.replace("T", "U")
    if not rna_sequence:
        raise ValueError("RNA sequence must not be empty")
    invalid = sorted(set(rna_sequence) - {"A", "C", "G", "U"})
    if invalid:
        raise ValueError(f"Invalid RNA bases: {', '.join(invalid)}")
    settings = _resolve_parameters(parameters)

    try:
        import RNA  # type: ignore
    except ModuleNotFoundError as exc:
        if exc.name != "RNA":
            raise
        raise RuntimeError(
            "ViennaRNA is not installed in the scorer environment. "
            "From the MORSE checkout run: uv sync --extra mfe. "
            "Alternatively install ViennaRNA==2.7.2 with that environment's pip."
        ) from exc

    model_details = RNA.md()
    if settings["max_bp_span"] > 0:
        model_details.max_bp_span = settings["max_bp_span"]
    fold_compound = RNA.fold_compound(rna_sequence, model_details)
    mfe_structure, mfe = fold_compound.mfe()
    gc_fraction = (dna_sequence.count("G") + dna_sequence.count("C")) / len(dna_sequence)

    outputs = {
        "mfe": float(mfe),
        "mfe_per_nt": float(mfe) / max(len(rna_sequence), 1),
        "gc_fraction": gc_fraction,
        "mfe_structure": mfe_structure,
    }
    if settings["compute_ensemble"]:
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
