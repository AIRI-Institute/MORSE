from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Any

import torch
import torch.nn as nn
import torch.nn.functional as F

from model import EXPECTED_TASK, OUTPUT_KEY, MORSEHead, prepare_candidate, verify_checkpoint


class PortableMamba(nn.Module):
    """State-dict-compatible Mamba-1 inference using only PyTorch operations."""

    def __init__(
        self,
        d_model: int,
        *,
        d_state: int = 16,
        d_conv: int = 4,
        expand: int = 2,
    ) -> None:
        super().__init__()
        self.d_model = d_model
        self.d_state = d_state
        self.d_conv = d_conv
        self.expand = expand
        self.d_inner = expand * d_model
        self.dt_rank = math.ceil(d_model / 16)
        self.in_proj = nn.Linear(d_model, self.d_inner * 2, bias=False)
        self.conv1d = nn.Conv1d(
            self.d_inner,
            self.d_inner,
            d_conv,
            groups=self.d_inner,
            padding=d_conv - 1,
            bias=True,
        )
        self.x_proj = nn.Linear(
            self.d_inner, self.dt_rank + 2 * d_state, bias=False
        )
        self.dt_proj = nn.Linear(self.dt_rank, self.d_inner, bias=True)
        self.A_log = nn.Parameter(torch.empty(self.d_inner, d_state))
        self.D = nn.Parameter(torch.empty(self.d_inner))
        self.out_proj = nn.Linear(self.d_inner, d_model, bias=False)
        with torch.no_grad():
            transition = torch.arange(1, d_state + 1, dtype=torch.float32)
            self.A_log.copy_(transition.log().repeat(self.d_inner, 1))
            self.D.fill_(1.0)

    def forward(self, hidden_states: torch.Tensor) -> torch.Tensor:
        batch, length, _ = hidden_states.shape
        projected_inputs = self.in_proj(hidden_states).transpose(1, 2)
        x_value, z_value = projected_inputs.chunk(2, dim=1)
        x_value = F.silu(self.conv1d(x_value)[..., :length])
        flat = x_value.transpose(1, 2).reshape(batch * length, self.d_inner)
        projected = self.x_proj(flat)
        delta, b_value, c_value = torch.split(
            projected,
            [self.dt_rank, self.d_state, self.d_state],
            dim=-1,
        )
        delta = (
            (self.dt_proj.weight @ delta.t())
            .t()
            .reshape(batch, length, self.d_inner)
            .permute(0, 2, 1)
        )
        b_value = (
            b_value.reshape(batch, length, self.d_state)
            .permute(0, 2, 1)
            .contiguous()
        )
        c_value = (
            c_value.reshape(batch, length, self.d_state)
            .permute(0, 2, 1)
            .contiguous()
        )
        output = _selective_scan_reference(
            x_value,
            delta,
            -torch.exp(self.A_log.float()),
            b_value,
            c_value,
            self.D.float(),
            z_value,
            self.dt_proj.bias.float(),
        )
        return self.out_proj(output.transpose(1, 2))


class PortableMambaBlock(nn.Module):
    def __init__(
        self,
        dimension: int,
        *,
        d_state: int,
        d_conv: int,
        expand: int,
    ) -> None:
        super().__init__()
        self.mixer = PortableMamba(
            dimension,
            d_state=d_state,
            d_conv=d_conv,
            expand=expand,
        )
        self.norm = nn.LayerNorm(dimension)

    def forward(
        self,
        hidden_states: torch.Tensor,
        residual: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        residual = (
            hidden_states + residual if residual is not None else hidden_states
        )
        return self.mixer(self.norm(residual.to(self.norm.weight.dtype))), residual


class PortableOrthrusBackbone(nn.Module):
    def __init__(
        self,
        *,
        n_tracks: int,
        dimension: int,
        layers: int,
        d_state: int,
        d_conv: int,
        expand: int,
    ) -> None:
        super().__init__()
        self.embedding = nn.Linear(n_tracks, dimension)
        self.layers = nn.ModuleList(
            PortableMambaBlock(
                dimension,
                d_state=d_state,
                d_conv=d_conv,
                expand=expand,
            )
            for _ in range(layers)
        )
        self.norm_f = nn.LayerNorm(dimension)

    def representation(
        self, inputs: torch.Tensor, lengths: torch.Tensor
    ) -> torch.Tensor:
        hidden = self.embedding(inputs)
        residual: torch.Tensor | None = None
        for layer in self.layers:
            hidden, residual = layer(hidden, residual)
        merged = hidden + residual if residual is not None else hidden
        hidden = self.norm_f(merged.to(self.norm_f.weight.dtype))
        positions = torch.arange(hidden.shape[1], device=hidden.device)[None, :]
        mask = positions < lengths[:, None]
        return (hidden * mask.unsqueeze(-1)).sum(dim=1) / lengths.unsqueeze(
            -1
        ).float()


class PortableOrthrusRegressor(nn.Module):
    def __init__(
        self, backbone: PortableOrthrusBackbone, dimension: int
    ) -> None:
        super().__init__()
        self.backbone = backbone
        self.head = MORSEHead(dimension)

    def forward(
        self, inputs: torch.Tensor, lengths: torch.Tensor
    ) -> torch.Tensor:
        return self.head(self.backbone.representation(inputs, lengths)).squeeze(
            -1
        )


class PortableOrthrusPredictor:
    def __init__(self, checkpoint_path: str | Path | None = None) -> None:
        self.checkpoint_path = Path(
            checkpoint_path
            or os.environ.get("ORTHRUS_CHECKPOINT", "/models/model_final.pt")
        ).expanduser()
        if not self.checkpoint_path.is_file():
            raise FileNotFoundError(
                f"Orthrus checkpoint not found: {self.checkpoint_path}"
            )
        verify_checkpoint(self.checkpoint_path)
        checkpoint = torch.load(
            self.checkpoint_path, map_location="cpu", weights_only=False
        )
        if checkpoint.get("task") != EXPECTED_TASK:
            raise ValueError(
                f"Expected task {EXPECTED_TASK!r}, got {checkpoint.get('task')!r}"
            )
        state = checkpoint.get("model_state")
        if not isinstance(state, dict):
            raise ValueError("Orthrus checkpoint does not contain model_state")
        architecture = _infer_architecture(state)
        self.model = PortableOrthrusRegressor(
            PortableOrthrusBackbone(**architecture),
            architecture["dimension"],
        )
        self.model.load_state_dict(state, strict=True)
        self.device = _resolve_device()
        self.model.to(self.device).eval()
        self.target_mean = float(checkpoint["target_mean"])
        self.target_std = float(checkpoint["target_std"])
        self.cds_track_mode = str(
            checkpoint.get("input_tracks", {}).get("cds")
            or "continuous_cds_mask"
        )
        self.max_batch_size = _positive_int("ORTHRUS_BATCH_SIZE", 2)
        self.model_id = str(checkpoint.get("model", {}).get("model_id", "orthrus"))

    def predict(
        self,
        sequences: list[str],
        *,
        sequence_contexts: list[dict[str, Any]] | None = None,
    ) -> list[dict[str, float]]:
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
        predictions = [0.0] * len(prepared)
        indexed = sorted(enumerate(prepared), key=lambda item: item[1].length)
        with torch.inference_mode():
            for start in range(0, len(indexed), self.max_batch_size):
                batch = indexed[start : start + self.max_batch_size]
                row_indices, candidates = zip(*batch, strict=True)
                lengths = torch.tensor(
                    [candidate.length for candidate in candidates],
                    dtype=torch.long,
                    device=self.device,
                )
                inputs = torch.zeros(
                    (
                        len(candidates),
                        int(lengths.max().item()),
                        6,
                    ),
                    dtype=torch.float32,
                    device=self.device,
                )
                for row, candidate in enumerate(candidates):
                    inputs[row, : candidate.length] = torch.from_numpy(
                        candidate.tracks
                    ).to(self.device)
                values = (
                    self.model(inputs, lengths).float().cpu() * self.target_std
                    + self.target_mean
                )
                for row_index, value in zip(
                    row_indices, values.tolist(), strict=True
                ):
                    predictions[row_index] = float(value)
        return [{OUTPUT_KEY: value} for value in predictions]


def _infer_architecture(state: dict[str, torch.Tensor]) -> dict[str, int]:
    embedding = state["backbone.embedding.weight"]
    dimension, n_tracks = map(int, embedding.shape)
    layer_indices = {
        int(key.split(".")[2])
        for key in state
        if key.startswith("backbone.layers.")
    }
    prefix = "backbone.layers.0.mixer."
    d_state = int(state[prefix + "A_log"].shape[1])
    d_conv = int(state[prefix + "conv1d.weight"].shape[2])
    expand = int(
        state[prefix + "in_proj.weight"].shape[0] // (2 * dimension)
    )
    return {
        "n_tracks": n_tracks,
        "dimension": dimension,
        "layers": len(layer_indices),
        "d_state": d_state,
        "d_conv": d_conv,
        "expand": expand,
    }


def _selective_scan_reference(
    inputs: torch.Tensor,
    delta: torch.Tensor,
    a_value: torch.Tensor,
    b_value: torch.Tensor,
    c_value: torch.Tensor,
    d_value: torch.Tensor,
    z_value: torch.Tensor,
    delta_bias: torch.Tensor,
) -> torch.Tensor:
    input_dtype = inputs.dtype
    inputs_float = inputs.float()
    delta_float = F.softplus(delta.float() + delta_bias[..., None].float())
    state = a_value.new_zeros(
        (inputs.shape[0], inputs.shape[1], a_value.shape[1])
    )
    outputs: list[torch.Tensor] = []
    for index in range(inputs.shape[2]):
        step = delta_float[:, :, index]
        delta_a = torch.exp(torch.einsum("bd,dn->bdn", step, a_value))
        delta_b_u = torch.einsum(
            "bd,bn,bd->bdn",
            step,
            b_value[:, :, index].float(),
            inputs_float[:, :, index],
        )
        state = delta_a * state + delta_b_u
        outputs.append(
            torch.einsum("bdn,bn->bd", state, c_value[:, :, index].float())
        )
    output = torch.stack(outputs, dim=2)
    output = output + inputs_float * d_value.reshape(1, -1, 1)
    return (output * F.silu(z_value.float())).to(input_dtype)


def _resolve_device() -> torch.device:
    requested = os.environ.get("ORTHRUS_DEVICE", "auto").strip().lower()
    if requested == "auto":
        requested = "cuda" if torch.cuda.is_available() else "cpu"
    device = torch.device(requested)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is not available")
    return device


def _positive_int(name: str, default: int) -> int:
    value = int(os.environ.get(name, str(default)))
    if value < 1:
        raise ValueError(f"{name} must be positive")
    return value
