from __future__ import annotations

import json
import math
import os
import sys
import traceback
from pathlib import Path
from typing import Any

if __package__:
    from .model import CaiScorer, load_profile
else:
    from model import CaiScorer, load_profile


DEFAULT_PROFILE = Path(__file__).with_name("profiles") / "homo_sapiens_kazusa_93487_cds.json"


def main() -> int:
    try:
        scorer = from_environment()
    except Exception:
        traceback.print_exc(file=sys.stderr)
        return 1

    print(
        f"CAI scorer ready profile={scorer.profile.profile_id} "
        f"sha256={scorer.profile.sha256} low_weight_threshold={scorer.low_weight_threshold:g}",
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
            traceback.print_exc(file=sys.stderr)
            response = _error_response(line, str(exc))
        print(
            json.dumps(response, separators=(",", ":"), sort_keys=True, allow_nan=False),
            flush=True,
        )
    return 0


def handle_request(request: dict[str, Any], scorer: CaiScorer) -> dict[str, Any]:
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
        "outputs": [effective.score_outputs(sequence) for sequence in sequences],
        "parameters": {
            "profile_path": getattr(effective, "profile_path", None),
            "low_weight_threshold": effective.low_weight_threshold,
        },
    }


def from_environment() -> CaiScorer:
    path = Path(os.environ.get("CAI_PROFILE", str(DEFAULT_PROFILE)))
    scorer = CaiScorer(
        load_profile(path),
        low_weight_threshold=float(os.environ.get("CAI_LOW_WEIGHT_THRESHOLD", "0.1")),
    )
    scorer.profile_path = str(path.resolve())
    return scorer


def _for_request(request: dict[str, Any], baseline: CaiScorer) -> CaiScorer:
    if "parameters" in request and "params" in request:
        raise ValueError("Use either parameters or params, not both")
    parameters = request.get("parameters", request.get("params", {}))
    if not isinstance(parameters, dict):
        raise ValueError("parameters must be a mapping")
    unknown = set(parameters) - {"profile_path", "low_weight_threshold"}
    if unknown:
        raise ValueError("Unknown CAI parameters: " + ", ".join(sorted(map(str, unknown))))
    if not parameters:
        return baseline
    threshold = parameters.get("low_weight_threshold", baseline.low_weight_threshold)
    if type(threshold) not in (int, float) or not math.isfinite(threshold):
        raise ValueError("low_weight_threshold must be a finite number")
    profile = baseline.profile
    path = getattr(baseline, "profile_path", None)
    if "profile_path" in parameters:
        raw_path = parameters["profile_path"]
        if not isinstance(raw_path, str) or not raw_path.strip():
            raise ValueError("profile_path must be a nonempty string")
        selected = _resolve_profile_path(raw_path)
        profile = load_profile(selected)
        path = str(selected.resolve())
    scorer = CaiScorer(profile, low_weight_threshold=threshold)
    scorer.profile_path = path
    return scorer


def _resolve_profile_path(value: str) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute() or path.is_file():
        return path
    if path.parts[:2] == ("scorers", "cai"):
        path = Path(*path.parts[2:])
    local = Path(__file__).parent / path
    if local.is_file():
        return local
    bundled = DEFAULT_PROFILE.parent / path
    return bundled if bundled.is_file() else path


def _error_response(line: str, message: str) -> dict[str, Any]:
    try:
        raw = json.loads(line)
    except Exception:
        return {"request_id": None, "error": message}
    request_id = raw.get("request_id") if isinstance(raw, dict) else None
    return {"request_id": request_id, "error": message}


if __name__ == "__main__":
    raise SystemExit(main())
