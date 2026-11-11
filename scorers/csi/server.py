from __future__ import annotations

import json
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
        profile_path = Path(os.environ.get("CSI_PROFILE", str(DEFAULT_PROFILE)))
        scorer = CsiScorer(
            load_profile(profile_path),
            window_codons=int(os.environ.get("CSI_WINDOW_CODONS", "15")),
            destabilizing_threshold=float(
                os.environ.get("CSI_DESTABILIZING_THRESHOLD", "0.0")
            ),
        )
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
    unknown = sorted(set(request) - {"request_id", "sequences", "sequence_contexts"})
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
    return {
        "request_id": request.get("request_id"),
        "outputs": [scorer.score_outputs(sequence) for sequence in sequences],
    }


def _error_response(line: str, message: str) -> dict[str, Any]:
    try:
        raw = json.loads(line)
    except json.JSONDecodeError:
        return {"request_id": None, "error": message}
    request_id = raw.get("request_id") if isinstance(raw, dict) else None
    return {"request_id": request_id, "error": message}


if __name__ == "__main__":
    raise SystemExit(main())
