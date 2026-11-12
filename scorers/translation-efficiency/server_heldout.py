"""Inference-only JSONL adapter for the paper's common-holdout TE checkpoint."""
from __future__ import annotations
import copy
import hashlib
import json
import os
from pathlib import Path
import sys
import traceback
import numpy as np

CHECKPOINT_SHA = "70913fc29ba9c073100913b56f59542520ec18a473fae358ff19bd57ecddd604"


def _asset_root(environ: dict[str, str] | None = None, source_file: str | Path | None = None) -> Path:
    """Resolve assets in a checkout, installed package, or flat /app container."""
    env = os.environ if environ is None else environ
    if env.get("MORSE_ASSETS_DIR"):
        return Path(env["MORSE_ASSETS_DIR"]).expanduser()
    if env.get("MORSE_ROOT"):
        return Path(env["MORSE_ROOT"]).expanduser() / ".cache/morse/assets"
    source = Path(source_file or __file__).resolve()
    for ancestor in source.parents:
        if (ancestor / "assets/registry.json").is_file() and (ancestor / "scorers").is_dir():
            return ancestor / ".cache/morse/assets"
    return Path.home() / ".cache/morse/assets"


def inference_paths(environ: dict[str, str] | None = None, source_file: str | Path | None = None) -> tuple[Path, Path]:
    env = os.environ if environ is None else environ
    checkpoint = env.get("ORTHRUS_CHECKPOINT")
    backbone = env.get("ORTHRUS_MODEL_SOURCE")
    # Explicit mounts need no source-layout or home-directory lookup at all.
    if checkpoint and backbone:
        return Path(checkpoint).expanduser(), Path(backbone).expanduser()
    root = _asset_root(env, source_file)
    return (Path(checkpoint).expanduser() if checkpoint else root / "orthrus-te-held-out/best.pt",
            Path(backbone).expanduser() if backbone else root / "orthrus-backbone")


class Predictor:
    def __init__(self) -> None:
        import torch
        from transformers import AutoModel
        from model import OrthrusRegressor, _resolve_device
        path, backbone_path = inference_paths()
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest() != CHECKPOINT_SHA:
            raise ValueError("The held-out TE adapter requires the frozen common-holdout checkpoint")
        self.checkpoint = torch.load(path, map_location="cpu", weights_only=False)
        if self.checkpoint.get("task") != "translation_efficiency_common_holdout":
            raise ValueError("Wrong checkpoint task")
        backbone = AutoModel.from_pretrained(
            str(backbone_path),
            trust_remote_code=True, local_files_only=True)
        self.model = OrthrusRegressor(backbone)
        self.model.load_state_dict(self.checkpoint["model_state"], strict=True)
        self.device = _resolve_device()
        torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "2")))
        self.model.to(self.device).eval()
        self.batch_size = int(os.environ.get("ORTHRUS_BATCH_SIZE", "4"))
        self.max_tokens = None
        if self.batch_size < 1:
            raise ValueError("ORTHRUS_BATCH_SIZE must be positive")

    def predict(self, sequences: list[str], contexts: list[dict]) -> list[dict]:
        import torch
        from heldout_encoding import prepare_tracks
        from model import PreparedCandidate, _token_batches
        if len(sequences) != len(contexts):
            raise ValueError("One context is required per sequence")
        result = []
        prepared = [PreparedCandidate(prepare_tracks(seq, ctx))
                    for seq, ctx in zip(sequences, contexts)]
        with torch.inference_mode():
            for batch in _token_batches(list(enumerate(prepared)),
                                        max_batch_size=self.batch_size,
                                        max_tokens=self.max_tokens):
                tracks = [candidate.tracks for _, candidate in batch]
                lengths = torch.tensor([len(row) for row in tracks], dtype=torch.long)
                inputs = torch.zeros(len(tracks), int(lengths.max()), 6)
                for i, row in enumerate(tracks):
                    inputs[i, :len(row)] = torch.from_numpy(row)
                with torch.autocast(self.device.type, dtype=torch.bfloat16, enabled=self.device.type == "cuda"):
                    z = self.model(inputs.to(self.device), lengths.to(self.device))
                values = z.float().cpu().numpy() * self.checkpoint["target_std"] + self.checkpoint["target_mean"]
                if not np.isfinite(values).all():
                    raise ValueError("Nonfinite TE score")
                result.extend({"translation_efficiency_score": float(value)} for value in values)
        return result


def handle_request(predictor, request: dict) -> dict:
    if not isinstance(request, dict):
        raise ValueError("Request must be a JSON object")
    unknown = set(request) - {"request_id", "sequences", "sequence_contexts", "parameters", "params"}
    if unknown:
        raise ValueError("Unknown request fields: " + ", ".join(sorted(map(str, unknown))))
    sequences = request.get("sequences")
    contexts = request.get("sequence_contexts")
    if not isinstance(sequences, list) or not all(isinstance(s, str) for s in sequences):
        raise ValueError("sequences must be a list of strings")
    if not isinstance(contexts, list) or len(contexts) != len(sequences) or not all(isinstance(c, dict) for c in contexts):
        raise ValueError("One explicit sequence context is required per sequence")
    if "parameters" in request and "params" in request:
        raise ValueError("Use either parameters or params, not both")
    raw = request.get("parameters", request.get("params", {}))
    if not isinstance(raw, dict):
        raise ValueError("parameters must be a mapping")
    unknown = set(raw) - {"batch_size", "max_tokens"}
    if unknown:
        raise ValueError("Unknown Orthrus parameters: " + ", ".join(sorted(map(str, unknown))))
    parameters = {"batch_size": raw.get("batch_size", predictor.batch_size),
                  "max_tokens": raw.get("max_tokens", getattr(predictor, "max_tokens", None))}
    if type(parameters["batch_size"]) is not int or parameters["batch_size"] < 1:
        raise ValueError("batch_size must be a positive integer")
    tokens = parameters["max_tokens"]
    if tokens is not None and (type(tokens) is not int or tokens < 1):
        raise ValueError("max_tokens must be a positive integer or null")
    configured = copy.copy(predictor)
    configured.batch_size = parameters["batch_size"]
    configured.max_tokens = tokens
    return {"request_id": request.get("request_id"), "parameters": parameters,
            "outputs": configured.predict(sequences, contexts)}


def main() -> int:
    predictor = Predictor()
    print("Held-out TE ready; continuous CDS mask includes stop", file=sys.stderr, flush=True)
    for line in sys.stdin:
        if not line.strip():
            continue
        request = {}
        try:
            request = json.loads(line)
            response = handle_request(predictor, request)
        except Exception as error:
            traceback.print_exc(file=sys.stderr)
            response = {"request_id": request.get("request_id") if isinstance(request, dict) else None, "error": str(error)}
        print(json.dumps(response, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
