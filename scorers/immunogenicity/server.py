from __future__ import annotations

import json
import hashlib
import os
import sys
import traceback
from pathlib import Path
from typing import Any

import numpy as np

DEFAULT_WINDOW_SIZE = 27
DEFAULT_WINDOW_STRIDE = 3


def main() -> int:
    scorer = RinalmoWindowScorer.from_env()
    for line in sys.stdin:
        if not line.strip():
            continue
        try:
            request = json.loads(line)
            response = scorer.handle_request(request)
        except Exception as exc:  # noqa: BLE001 - scorer must serialize failures to the caller.
            traceback.print_exc(file=sys.stderr)
            response = {
                "request_id": _safe_request_id(line),
                "error": f"{exc.__class__.__name__}: {exc}",
            }
        print(json.dumps(response), flush=True)
    return 0


class RinalmoWindowScorer:
    def __init__(
        self,
        checkpoint_path: Path,
        window_size: int,
        stride: int,
        offsets: list[int],
        batch_size: int,
    ) -> None:
        import torch
        from multimolecule import RnaTokenizer

        from runtime_model import RinalmoSequenceClassifier, normalize_rna

        self.normalize_rna = normalize_rna
        self.window_size = window_size
        self.stride = stride
        self.offsets = offsets
        self.batch_size = batch_size
        self.torch = torch
        requested_device = os.environ.get("RINALMO_DEVICE", "auto").strip().lower()
        if requested_device == "auto":
            requested_device = "cuda" if torch.cuda.is_available() else "cpu"
        if requested_device.startswith("cuda") and not torch.cuda.is_available():
            raise RuntimeError(
                f"RINALMO_DEVICE={requested_device} requested but CUDA is unavailable"
            )
        self.device = torch.device(requested_device)

        expected = os.environ.get("MORSE_EXPECTED_CHECKPOINT_SHA256")
        if expected:
            digest = hashlib.sha256()
            with checkpoint_path.open("rb") as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(block)
            if digest.hexdigest() != expected:
                raise ValueError("RiNALMo checkpoint SHA-256 mismatch")
        checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
        model_config = checkpoint.get("model_config")
        if model_config is None:
            config_path = Path(
                os.environ.get(
                    "RINALMO_CONFIG",
                    str(Path(__file__).resolve().parent / "model_assets/rinalmo_mega_config.json"),
                )
            )
            model_config = json.loads(config_path.read_text())
        model = RinalmoSequenceClassifier(
            dropout=0.0,
            backbone_config=model_config,
        )
        model.load_state_dict(checkpoint["model_state_dict"])
        model.to(self.device)
        model.eval()
        self.model = model
        tokenizer_path = os.environ.get(
            "RINALMO_TOKENIZER",
            str(Path(__file__).resolve().parent / "model_assets/tokenizer"),
        )
        self.tokenizer = RnaTokenizer.from_pretrained(
            tokenizer_path,
            local_files_only=True,
        )
        print(
            "Loaded Rinalmo window scorer "
            f"device={self.device} window_size={window_size} stride={stride} "
            f"offsets={offsets or 'default'} batch_size={batch_size}",
            file=sys.stderr,
            flush=True,
        )

    @classmethod
    def from_env(cls) -> "RinalmoWindowScorer":
        checkpoint_path = Path(
            os.environ.get("RINALMO_CHECKPOINT", "/models/best_model.pt")
        )
        if not checkpoint_path.is_file():
            raise FileNotFoundError(
                f"RiNALMo checkpoint not found: {checkpoint_path}. "
                "Set RINALMO_CHECKPOINT or mount the model at /models/best_model.pt."
            )
        return cls(
            checkpoint_path=checkpoint_path,
            window_size=_parse_positive_int("RINALMO_WINDOW_SIZE", DEFAULT_WINDOW_SIZE),
            stride=_parse_positive_int("RINALMO_WINDOW_STRIDE", DEFAULT_WINDOW_STRIDE),
            offsets=_parse_offsets(os.environ.get("RINALMO_WINDOW_OFFSETS", "")),
            batch_size=_parse_positive_int("RINALMO_BATCH_SIZE", 256),
        )

    def handle_request(self, request: dict[str, Any]) -> dict[str, Any]:
        unknown = sorted(
            set(request) - {"request_id", "sequences", "sequence_contexts"}
        )
        if unknown:
            raise ValueError("Unknown request fields: " + ", ".join(unknown))
        sequences = request.get("sequences")
        if not isinstance(sequences, list) or not all(
            isinstance(sequence, str) for sequence in sequences
        ):
            raise ValueError("request.sequences must be a list of strings")

        unique_windows: list[str] = []
        window_indices: dict[str, int] = {}
        candidate_indices: list[np.ndarray] = []
        for raw_sequence in sequences:
            sequence = self.normalize_rna(raw_sequence)
            windows = self._windows(sequence)
            indices = []
            for window in windows:
                index = window_indices.get(window)
                if index is None:
                    index = len(unique_windows)
                    window_indices[window] = index
                    unique_windows.append(window)
                indices.append(index)
            candidate_indices.append(np.asarray(indices, dtype=np.int64))

        unique_logits = self._score_sequences(unique_windows)
        unique_probs = 1.0 / (1.0 + np.exp(-np.clip(unique_logits, -60, 60)))

        outputs = []
        for indices in candidate_indices:
            candidate_logits = unique_logits[indices]
            candidate_probs = unique_probs[indices]
            outputs.append(
                _aggregate_outputs(candidate_logits, candidate_probs)
            )

        return {"request_id": request.get("request_id"), "outputs": outputs}

    def _windows(self, sequence: str) -> list[str]:
        return _sliding_windows(sequence, self.window_size, self.stride, self.offsets)

    def _score_sequences(self, sequences: list[str]) -> np.ndarray:
        logits_all: list[np.ndarray] = []
        torch = self.torch
        with torch.no_grad():
            for index in range(0, len(sequences), self.batch_size):
                batch = sequences[index : index + self.batch_size]
                tok = self.tokenizer(
                    batch,
                    return_tensors="pt",
                    padding=True,
                    return_special_tokens_mask=True,
                )
                tok = {key: value.to(self.device) for key, value in tok.items()}
                with torch.autocast(
                    device_type="cuda",
                    dtype=torch.bfloat16,
                    enabled=self.device.type == "cuda",
                ):
                    out = self.model(**tok)
                logits_all.append(out["logits"].float().cpu().numpy().reshape(-1))
        if not logits_all:
            return np.asarray([], dtype=float)
        return np.concatenate(logits_all).astype(float)


def _aggregate_outputs(logits: np.ndarray, probs: np.ndarray) -> dict[str, float | int]:
    if len(logits) == 0:
        raise ValueError("Cannot aggregate an empty set of RiNALMo windows")
    return {
        "rinalmo_window_count": int(len(logits)),
        "rinalmo_window_logit_mean": float(np.mean(logits)),
        "rinalmo_window_logit_median": float(np.median(logits)),
        "rinalmo_window_logit_max": float(np.max(logits)),
        "rinalmo_window_logit_p95": float(np.percentile(logits, 95)),
        "rinalmo_window_probability_mean": float(np.mean(probs)),
        "rinalmo_window_probability_max": float(np.max(probs)),
        "rinalmo_window_probability_p95": float(np.percentile(probs, 95)),
    }


def _parse_offsets(raw: str) -> list[int]:
    if not raw.strip():
        return []
    offsets = [int(part.strip()) for part in raw.split(",") if part.strip()]
    if any(offset < 0 for offset in offsets):
        raise ValueError("RINALMO_WINDOW_OFFSETS values must be >= 0")
    return offsets


def _parse_positive_int(name: str, default: int) -> int:
    value = int(os.environ.get(name, str(default)))
    if value < 1:
        raise ValueError(f"{name} must be >= 1")
    return value


def _sliding_windows(
    sequence: str,
    window_size: int,
    stride: int,
    offsets: list[int] | None = None,
) -> list[str]:
    if window_size < 1:
        raise ValueError("window_size must be >= 1")
    if stride < 1:
        raise ValueError("stride must be >= 1")
    if offsets and any(offset < 0 for offset in offsets):
        raise ValueError("offsets must be >= 0")
    if len(sequence) <= window_size:
        return [sequence]

    max_start = len(sequence) - window_size
    starts: set[int] = set()
    for offset in offsets or [0]:
        if offset <= max_start:
            starts.update(range(offset, max_start + 1, stride))
    starts.add(max_start)
    return [sequence[start : start + window_size] for start in sorted(starts)]


def _safe_request_id(line: str) -> Any:
    try:
        return json.loads(line).get("request_id")
    except Exception:  # noqa: BLE001
        return None


if __name__ == "__main__":
    raise SystemExit(main())
