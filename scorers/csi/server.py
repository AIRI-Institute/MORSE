from __future__ import annotations

import json
import math
import os
import sys
import traceback
from pathlib import Path
from typing import Any

if __package__:
    from .model import CsiScorer, load_profile
else:
    from model import CsiScorer, load_profile

DEFAULT_PROFILE = (
    Path(__file__).with_name("profiles") / "homo_sapiens_consensus_wu2019_v1.json"
)


def main() -> int:
    try:
        scorer = from_environment()
    except Exception:  # noqa: BLE001 - startup failures must be reported on stderr
        traceback.print_exc(file=sys.stderr)
        return 1

    print(
        f"CSI scorer ready profile={scorer.profile.profile_id} "
        f"sha256={scorer.profile.sha256} window_codons={scorer.window_codons}",
        file=sys.stderr,
        flush=True,
    )
    for line in sys.stdin:
        if not line.strip():
            continue
        try:
            request = json.loads(line)
            if not isinstance(request, dict):
                raise TypeError("Request must be a JSON object")
            response = handle_request(request, scorer)
        except Exception as exc:  # noqa: BLE001 - keep the JSONL process alive
            traceback.print_exc(file=sys.stderr)
            response = _error_response(line, str(exc))
        print(
            json.dumps(
                response, separators=(",", ":"), sort_keys=True, allow_nan=False
            ),
            flush=True,
        )
    return 0


def handle_request(request: dict[str, Any], scorer: CsiScorer) -> dict[str, Any]:
    if not isinstance(request, dict):
        raise ValueError("Request must be a JSON object")
    unknown = sorted(set(request) - {"request_id", "sequences", "sequence_contexts", "parameters", "params"})
    if unknown:
        raise ValueError("Unknown request fields: " + ", ".join(unknown))
    sequences = request.get("sequences")
    if not isinstance(sequences, list) or not all(
        isinstance(sequence, str) for sequence in sequences
    ):
        raise ValueError("Request must define a list of sequences")
    contexts = request.get("sequence_contexts")
    if contexts is not None and (
        not isinstance(contexts, list) or len(contexts) != len(sequences)
    ):
        raise ValueError("sequence_contexts must contain one context per sequence")
    effective = _for_request(request, scorer)
    return {
        "request_id": request.get("request_id"),
        "outputs": [effective.score_outputs(sequence) for sequence in sequences],
        "parameters": {
            "profile_path": getattr(effective, "profile_path", None),
            "window_codons": effective.window_codons,
            "destabilizing_threshold": effective.destabilizing_threshold,
        },
    }


def from_environment() -> CsiScorer:
    path = Path(os.environ.get("CSI_PROFILE", str(DEFAULT_PROFILE)))
    scorer = CsiScorer(
        load_profile(path),
        window_codons=int(os.environ.get("CSI_WINDOW_CODONS", "15")),
        destabilizing_threshold=float(os.environ.get("CSI_DESTABILIZING_THRESHOLD", "0.0")),
    )
    scorer.profile_path = str(path.resolve())
    return scorer


def _for_request(request: dict[str, Any], baseline: CsiScorer) -> CsiScorer:
    if "parameters" in request and "params" in request:
        raise ValueError("Use either parameters or params, not both")
    parameters = request.get("parameters", request.get("params", {}))
    if not isinstance(parameters, dict):
        raise ValueError("parameters must be a mapping")
    unknown = set(parameters) - {"profile_path", "window_codons", "destabilizing_threshold"}
    if unknown:
        raise ValueError("Unknown CSI parameters: " + ", ".join(sorted(map(str, unknown))))
    if not parameters:
        return baseline
    window = parameters.get("window_codons", baseline.window_codons)
    threshold = parameters.get("destabilizing_threshold", baseline.destabilizing_threshold)
    if type(window) is not int or window < 1:
        raise ValueError("window_codons must be a positive integer")
    if type(threshold) not in (int, float) or not math.isfinite(threshold):
        raise ValueError("destabilizing_threshold must be a finite number")
    profile = baseline.profile
    path = getattr(baseline, "profile_path", None)
    if "profile_path" in parameters:
        raw_path = parameters["profile_path"]
        if not isinstance(raw_path, str) or not raw_path.strip():
            raise ValueError("profile_path must be a nonempty string")
        selected = _resolve_profile_path(raw_path)
        profile = load_profile(selected)
        path = str(selected.resolve())
    scorer = CsiScorer(profile, window_codons=window, destabilizing_threshold=threshold)
    scorer.profile_path = path
    return scorer


def _resolve_profile_path(value: str) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute() or path.is_file():
        return path
    if path.parts[:2] == ("scorers", "csi"):
        path = Path(*path.parts[2:])
    local = Path(__file__).parent / path
    if local.is_file():
        return local
    bundled = DEFAULT_PROFILE.parent / path
    return bundled if bundled.is_file() else path


def _error_response(line: str, message: str) -> dict[str, Any]:
    try:
        raw = json.loads(line)
    except json.JSONDecodeError:
        return {"request_id": None, "error": message}
    request_id = raw.get("request_id") if isinstance(raw, dict) else None
    return {"request_id": request_id, "error": message}


if __name__ == "__main__":
    raise SystemExit(main())
