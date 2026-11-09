from __future__ import annotations

import os
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import torch
import torch.nn as nn


EXPECTED_TASK = "translation_efficiency"
OUTPUT_KEY = "translation_efficiency_score"
DEFAULT_CHECKPOINT = Path("/models/model_final.pt")
DEFAULT_BACKBONE = Path("/models/orthrus-backbone")

_BASE_TO_COLUMN = {"A": 0, "C": 1, "G": 2, "T": 3}


class MORSEHead(nn.Module):
    def __init__(self, hidden: int) -> None:
        super().__init__()
        self.layers = nn.Sequential(
            nn.LayerNorm(hidden),
            nn.Linear(hidden, hidden // 2),
            nn.GELU(),
            nn.Dropout(0.15),
            nn.Linear(hidden // 2, 1),
        )

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        return self.layers(inputs)


class OrthrusRegressor(nn.Module):
    def __init__(self, backbone: nn.Module) -> None:
        super().__init__()
        self.backbone = backbone
        hidden = int(backbone.config.ssm_model_dim)
        self.head = MORSEHead(hidden)

    def forward(self, inputs: torch.Tensor, lengths: torch.Tensor) -> torch.Tensor:
        representation = self.backbone.representation(
            inputs, lengths, channel_last=True
        )
        return self.head(representation).squeeze(-1)


@dataclass(frozen=True)
class PreparedCandidate:
    tracks: np.ndarray

    @property
    def length(self) -> int:
        return int(self.tracks.shape[0])


class OrthrusPredictor:
    def __init__(
        self,
        checkpoint_path: Path | None = None,
        backbone_path: Path | None = None,
    ) -> None:
        self.checkpoint_path = checkpoint_path or Path(
            os.environ.get("ORTHRUS_CHECKPOINT", str(DEFAULT_CHECKPOINT))
        )
        if not self.checkpoint_path.exists():
            raise FileNotFoundError(
                f"Checkpoint not found: {self.checkpoint_path}. "
                "Mount prepared assets at /models or set ORTHRUS_CHECKPOINT."
            )
        verify_checkpoint(self.checkpoint_path)
        checkpoint = torch.load(
            self.checkpoint_path, map_location="cpu", weights_only=False
        )
        if checkpoint.get("task") != EXPECTED_TASK:
            raise ValueError(
                f"Expected task {EXPECTED_TASK!r}, got {checkpoint.get('task')!r}"
            )
        model_config = dict(checkpoint["model"])
        recipe = dict(checkpoint["recipe"])
        if recipe.get("head") != "layernorm":
            raise ValueError("This scorer supports the selected LayerNorm head only")

        source = backbone_path or Path(
            os.environ.get("ORTHRUS_MODEL_SOURCE", str(DEFAULT_BACKBONE))
        )
        if not source.exists():
            raise FileNotFoundError(
                f"Orthrus backbone assets not found: {source}. "
                "Install orthrus-backbone with assets/manage.py, then set "
                "ORTHRUS_MODEL_SOURCE to that directory (mount it at /backbone in Docker)."
            )
        from transformers import AutoModel

        backbone = AutoModel.from_pretrained(
            str(source),
            trust_remote_code=True,
            local_files_only=True,
        )
        self.model = OrthrusRegressor(backbone)
        self.model.load_state_dict(checkpoint["model_state"], strict=True)
        self.device = _resolve_device()
        self.model.to(self.device).eval()
        self.target_mean = float(checkpoint["target_mean"])
        self.target_std = float(checkpoint["target_std"])
        self.cds_track_mode = str(
            checkpoint.get("input_tracks", {}).get("cds")
            or "continuous_cds_mask"
        )
        self.max_batch_size = _positive_int("ORTHRUS_BATCH_SIZE", 8)
        self.max_tokens = _positive_int("ORTHRUS_MAX_TOKENS", 32_000)
        self.model_id = str(model_config["model_id"])

    def predict(
        self,
        sequences: list[str],
        *,
        sequence_contexts: list[dict[str, Any]] | None = None,
    ) -> list[dict[str, Any]]:
        if sequence_contexts is None or len(sequence_contexts) != len(sequences):
            raise ValueError("sequence_contexts must match sequences in length")
        prepared = [
            prepare_candidate(
                sequence,
                cds_track_mode=self.cds_track_mode,
                sequence_context=context,
            )
            for sequence, context in zip(
                sequences, sequence_contexts, strict=True
            )
        ]
        predictions = np.empty(len(prepared), dtype=np.float32)
        indexed = sorted(enumerate(prepared), key=lambda item: item[1].length)
        with torch.inference_mode():
            for batch in _token_batches(
                indexed,
                max_batch_size=self.max_batch_size,
                max_tokens=self.max_tokens,
            ):
                row_indices, batch_candidates = zip(*batch)
                inputs, lengths = _collate(list(batch_candidates))
                with torch.autocast(
                    device_type=self.device.type,
                    dtype=torch.bfloat16,
                    enabled=self.device.type == "cuda",
                ):
                    predicted_z = self.model(
                        inputs.to(self.device, non_blocking=True),
                        lengths.to(self.device, non_blocking=True),
                    )
                values = (
                    predicted_z.float().cpu().numpy() * self.target_std
                    + self.target_mean
                )
                predictions[np.asarray(row_indices)] = values

        return [
            {OUTPUT_KEY: float(predictions[index])}
            for index in range(len(prepared))
        ]


def prepare_candidate(
    raw_sequence: str,
    *,
    cds_track_mode: str,
    sequence_context: dict[str, Any],
) -> PreparedCandidate:
    sequence = _normalize_sequence(raw_sequence)
    coding_start, coding_end = _resolve_context_coding_interval(
        sequence_context, sequence
    )
    if not (0 <= coding_start < coding_end <= len(sequence)):
        raise ValueError(
            f"Invalid coding interval [{coding_start}, {coding_end}) for sequence "
            f"length {len(sequence)}"
        )

    tracks = np.zeros((len(sequence), 6), dtype=np.float32)
    for position, base in enumerate(sequence):
        tracks[position, _BASE_TO_COLUMN[base]] = 1.0
    if cds_track_mode == "codon_start":
        tracks[coding_start:coding_end:3, 4] = 1.0
    elif cds_track_mode == "continuous_cds_mask":
        tracks[coding_start:coding_end, 4] = 1.0
    else:
        raise ValueError(f"Unsupported CDS track mode: {cds_track_mode}")
    # Synthetic IVT constructs have no natural exon-junction annotation.
    tracks[:, 5] = 0.0
    return PreparedCandidate(tracks=tracks)


def _resolve_context_coding_interval(
    context: dict[str, Any], sequence: str
) -> tuple[int, int]:
    ranges = context.get("part_ranges")
    if not isinstance(ranges, list):
        raise ValueError("sequence context must define part_ranges")
    by_part = {
        str(item.get("part")): item
        for item in ranges
        if isinstance(item, dict)
    }
    start = by_part.get("start_codon")
    cds = by_part.get("cds")
    stop = by_part.get("stop_codon")
    if start is None or cds is None or stop is None:
        raise ValueError(
            "Orthrus sequence context must include start_codon, cds, and stop_codon"
        )
    coding_start = int(start.get("start", -1))
    start_end = int(start.get("end", -1))
    body_start = int(cds.get("start", -1))
    body_end = int(cds.get("end", -1))
    stop_start = int(stop.get("start", -1))
    stop_end = int(stop.get("end", -1))
    if (
        start_end - coding_start != 3
        or start_end != body_start
        or body_end != stop_start
        or stop_end - stop_start != 3
        or (body_end - body_start) % 3 != 0
        or not 0 <= coding_start < stop_start < stop_end <= len(sequence)
    ):
        raise ValueError("Orthrus sequence context has non-contiguous coding parts")
    if sequence[coding_start:start_end] in {"TAA", "TAG", "TGA"}:
        raise ValueError("Orthrus start codon cannot be a stop codon")
    if sequence[stop_start:stop_end] not in {"TAA", "TAG", "TGA"}:
        raise ValueError("Orthrus stop_codon range must contain a stop codon")
    return coding_start, stop_start


def _normalize_sequence(sequence: str) -> str:
    normalized = "".join(sequence.upper().split()).replace("U", "T")
    if not normalized:
        raise ValueError("Sequence must not be empty")
    unknown = sorted(set(normalized) - set(_BASE_TO_COLUMN))
    if unknown:
        raise ValueError(f"Unsupported nucleotide symbols: {''.join(unknown)}")
    return normalized


def _token_batches(
    indexed: list[tuple[int, PreparedCandidate]],
    *,
    max_batch_size: int,
    max_tokens: int,
) -> list[list[tuple[int, PreparedCandidate]]]:
    batches: list[list[tuple[int, PreparedCandidate]]] = []
    current: list[tuple[int, PreparedCandidate]] = []
    current_max = 0
    for item in indexed:
        candidate_max = max(current_max, item[1].length)
        candidate_size = len(current) + 1
        if current and (
            candidate_size > max_batch_size
            or candidate_max * candidate_size > max_tokens
        ):
            batches.append(current)
            current = []
            current_max = 0
        current.append(item)
        current_max = max(current_max, item[1].length)
    if current:
        batches.append(current)
    return batches


def _collate(
    candidates: list[PreparedCandidate],
) -> tuple[torch.Tensor, torch.Tensor]:
    lengths = torch.tensor([candidate.length for candidate in candidates])
    inputs = torch.zeros(
        (len(candidates), int(lengths.max()), 6), dtype=torch.float32
    )
    for row, candidate in enumerate(candidates):
        inputs[row, : candidate.length] = torch.from_numpy(candidate.tracks)
    return inputs, lengths


def _resolve_device() -> torch.device:
    requested = os.environ.get("ORTHRUS_DEVICE", "auto").strip().lower()
    if requested == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    device = torch.device(requested)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is not available")
    return device


def _positive_int(name: str, default: int) -> int:
    value = int(os.environ.get(name, str(default)))
    if value < 1:
        raise ValueError(f"{name} must be positive")
    return value


def verify_checkpoint(path: Path) -> None:
    expected = os.environ.get("MORSE_EXPECTED_CHECKPOINT_SHA256")
    if not expected:
        return
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    if result.hexdigest() != expected:
        raise ValueError("Checkpoint SHA-256 does not match the declared inference asset")
