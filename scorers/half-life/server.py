from __future__ import annotations

import copy
import json
import os
import sys
import traceback
from typing import Any


def handle_request(
    predictor: Any, request: dict[str, Any]
) -> dict[str, Any]:
    unknown = sorted(
        set(request) - {"request_id", "sequences", "sequence_contexts", "parameters", "params"}
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
        raise ValueError("Request must define one sequence context per sequence")
    parameters = resolve_parameters(predictor, request)
    configured = copy.copy(predictor)
    configured.max_batch_size = parameters["batch_size"]
    configured.max_tokens = parameters["max_tokens"]
    return {
        "request_id": request.get("request_id"),
        "parameters": parameters,
        "outputs": configured.predict(
            sequences,
            sequence_contexts=sequence_contexts,
        ),
    }


def resolve_parameters(predictor: Any, request: dict[str, Any]) -> dict[str, Any]:
    if "parameters" in request and "params" in request:
        raise ValueError("Use either parameters or params, not both")
    raw = request.get("parameters", request.get("params", {}))
    if not isinstance(raw, dict):
        raise ValueError("parameters must be a mapping")
    unknown = set(raw) - {"batch_size", "max_tokens"}
    if unknown:
        raise ValueError("Unknown Orthrus parameters: " + ", ".join(sorted(map(str, unknown))))
    size = raw.get("batch_size", predictor.max_batch_size)
    tokens = raw.get("max_tokens", getattr(predictor, "max_tokens", None))
    if type(size) is not int or size < 1:
        raise ValueError("batch_size must be a positive integer")
    if tokens is not None and (type(tokens) is not int or tokens < 1):
        raise ValueError("max_tokens must be a positive integer or null")
    return {"batch_size": size, "max_tokens": tokens}


def main() -> int:
    from model import EXPECTED_TASK

    if os.environ.get("ORTHRUS_PORTABLE", "0") == "1":
        from portable_orthrus import PortableOrthrusPredictor

        predictor = PortableOrthrusPredictor()
    else:
        from model import OrthrusPredictor

        predictor = OrthrusPredictor()
    print(
        f"Orthrus scorer ready task={EXPECTED_TASK} "
        f"model={predictor.model_id} device={predictor.device}",
        file=sys.stderr,
    )
    for line in sys.stdin:
        if not line.strip():
            continue
        try:
            request = json.loads(line)
            response = handle_request(predictor, request)
        except Exception as exc:  # scorer failures must be visible to MORSE
            traceback.print_exc(file=sys.stderr)
            response = {"request_id": _request_id(line), "error": str(exc)}
        print(json.dumps(response, separators=(",", ":")), flush=True)
    return 0


def _request_id(line: str) -> Any:
    try:
        return json.loads(line).get("request_id")
    except Exception:
        return None


if __name__ == "__main__":
    raise SystemExit(main())
