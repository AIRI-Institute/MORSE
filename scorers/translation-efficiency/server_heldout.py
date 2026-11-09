"""Inference-only JSONL adapter for the paper's common-holdout TE checkpoint."""
from __future__ import annotations
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
        if self.batch_size < 1:
            raise ValueError("ORTHRUS_BATCH_SIZE must be positive")

    def predict(self, sequences: list[str], contexts: list[dict]) -> list[dict]:
        import torch
        from heldout_encoding import prepare_tracks
        if len(sequences) != len(contexts):
            raise ValueError("One context is required per sequence")
        result = []
        with torch.inference_mode():
            for offset in range(0, len(sequences), self.batch_size):
                tracks = [prepare_tracks(seq, ctx) for seq, ctx in zip(
                    sequences[offset:offset + self.batch_size], contexts[offset:offset + self.batch_size])]
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


def main() -> int:
    predictor = Predictor()
    print("Held-out TE ready; continuous CDS mask includes stop", file=sys.stderr, flush=True)
    for line in sys.stdin:
        if not line.strip():
            continue
        request = {}
        try:
            request = json.loads(line)
            unknown = set(request) - {"request_id", "sequences", "sequence_contexts"}
            if unknown:
                raise ValueError(f"Unknown request fields: {sorted(unknown)}")
            response = {"request_id": request.get("request_id"), "outputs": predictor.predict(
                request["sequences"], request["sequence_contexts"])}
        except Exception as error:
            traceback.print_exc(file=sys.stderr)
            response = {"request_id": request.get("request_id"), "error": str(error)}
        print(json.dumps(response, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
