"""Persistent JSONL server for the versioned feasibility preset."""

from __future__ import annotations

import json
import os
import sys
import traceback
from typing import Any, Mapping

from .feasibility import (
    DEFAULT_PRESET,
    FLOAT_PARAMETERS,
    JSON_PARAMETERS,
    STRING_PARAMETERS,
    STRATEGY_DEFAULTS,
    SequenceFeasibilityScorer,
)


def handle_request(
    request: dict[str, Any],
    scorer: SequenceFeasibilityScorer,
) -> dict[str, Any]:
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
        "outputs": [scorer.score_sequence(sequence) for sequence in sequences],
    }


def error_response(line: str, message: str) -> dict[str, Any]:
    try:
        raw = json.loads(line)
    except Exception:
        return {"request_id": None, "error": message}
    request_id = raw.get("request_id") if isinstance(raw, dict) else None
    return {"request_id": request_id, "error": message}


def from_environment(
    environment: Mapping[str, str] | None = None,
) -> SequenceFeasibilityScorer:
    env = os.environ if environment is None else environment
    preset = env.get("SEQUENCE_FEASIBILITY_PRESET", DEFAULT_PRESET)
    defaults = SequenceFeasibilityScorer(preset).parameters
    overrides = {}
    for key in {**STRATEGY_DEFAULTS, **defaults}:
        name = "SEQUENCE_FEASIBILITY_" + key.upper()
        if name not in env:
            continue
        raw = env[name]
        if raw.strip() == "preset":
            continue
        try:
            overrides[key] = (
                json.loads(raw)
                if key in JSON_PARAMETERS
                else raw if key in STRING_PARAMETERS
                else float(raw) if key in FLOAT_PARAMETERS else int(raw)
            )
        except (ValueError, TypeError) as exc:
            raise ValueError(f"Invalid startup parameter {name}") from exc
    return SequenceFeasibilityScorer(preset, **overrides)


def main() -> int:
    try:
        scorer = from_environment()
    except Exception:
        traceback.print_exc(file=sys.stderr)
        return 1
    print(
        f"Sequence feasibility scorer ready preset={scorer.profile_id} config_sha256={scorer.configuration_sha256}",
        file=sys.stderr,
        flush=True,
    )
    for line in sys.stdin:
        if not line.strip():
            continue
        try:
            request = json.loads(line)
            if not isinstance(request, dict):
                raise ValueError("Request must be a JSON object")
            response = handle_request(request, scorer)
        except Exception as exc:
            response = error_response(line, str(exc))
        print(
            json.dumps(
                response, sort_keys=True, separators=(",", ":"), allow_nan=False
            ),
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
