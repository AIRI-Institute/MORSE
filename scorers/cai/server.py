from __future__ import annotations

import json
import os
import sys
import traceback
from pathlib import Path
from typing import Any

from model import CaiScorer, load_profile


DEFAULT_PROFILE = Path(__file__).with_name("profiles") / "homo_sapiens_kazusa_93487_cds.json"


def main() -> int:
    try:
        profile_path = Path(os.environ.get("CAI_PROFILE", str(DEFAULT_PROFILE)))
        low_weight_threshold = float(os.environ.get("CAI_LOW_WEIGHT_THRESHOLD", "0.1"))
        scorer = CaiScorer(
            load_profile(profile_path),
            low_weight_threshold=low_weight_threshold,
        )
    except Exception:
        traceback.print_exc(file=sys.stderr)
        return 1

    print(
        f"CAI scorer ready profile={scorer.profile.profile_id} "
        f"sha256={scorer.profile.sha256} low_weight_threshold={low_weight_threshold:g}",
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
        "outputs": [scorer.score_outputs(sequence) for sequence in sequences],
    }


def _error_response(line: str, message: str) -> dict[str, Any]:
    try:
        raw = json.loads(line)
    except Exception:
        return {"request_id": None, "error": message}
    request_id = raw.get("request_id") if isinstance(raw, dict) else None
    return {"request_id": request_id, "error": message}


if __name__ == "__main__":
    raise SystemExit(main())
