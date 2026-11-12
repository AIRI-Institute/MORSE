"""Persistent JSONL server for the versioned feasibility preset."""

from __future__ import annotations

import copy
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
    if not isinstance(request, dict):
        raise ValueError("Request must be a JSON object")
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
    effective = _for_request(request, scorer)
    return {
        "request_id": request.get("request_id"),
        "outputs": [effective.score_sequence(sequence) for sequence in sequences],
        "parameters": {"preset": effective.profile_id, **effective.parameters},
    }


def _for_request(
    request: dict[str, Any], baseline: SequenceFeasibilityScorer,
) -> SequenceFeasibilityScorer:
    if "parameters" in request and "params" in request:
        raise ValueError("Use either parameters or params, not both")
    parameters = request.get("parameters", request.get("params", {}))
    if not isinstance(parameters, dict):
        raise ValueError("parameters must be a mapping")
    if not parameters:
        return baseline
    parameters = copy.deepcopy(parameters)
    preset = parameters.pop("preset", baseline.profile_id)
    if not isinstance(preset, str):
        raise ValueError("preset must be a string")
    if preset == baseline.profile_id:
        defaults = baseline.parameters
    else:
        defaults = getattr(baseline, "_startup_overrides", None)
        if defaults is None:
            original = SequenceFeasibilityScorer(baseline.profile_id).parameters
            defaults = {
                key: value for key, value in baseline.parameters.items()
                if key not in original or value != original[key]
            }
    return SequenceFeasibilityScorer(preset, **{**defaults, **parameters})


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
    scorer = SequenceFeasibilityScorer(preset, **overrides)
    scorer._startup_overrides = copy.deepcopy(overrides)
    return scorer


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
