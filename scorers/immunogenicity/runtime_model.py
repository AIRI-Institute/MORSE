"""Inference-only RiNALMo classifier compatible with the frozen IMrnaDS weights.

This module contains no model fitting, dataset loading, loss, or evaluation code.
The MultiMolecule dependency retains its own AGPL license; see THIRD_PARTY_NOTICES.
"""
from __future__ import annotations

from typing import Any
import torch
from torch import nn


class RinalmoSequenceClassifier(nn.Module):
    def __init__(self, *, backbone_config: dict[str, Any], dropout: float = 0.0) -> None:
        super().__init__()
        from multimolecule import RiNALMoConfig, RiNALMoModel
        self.backbone = RiNALMoModel(RiNALMoConfig.from_dict(backbone_config))
        hidden = next((getattr(self.backbone.config, key, None)
                       for key in ("hidden_size", "d_model", "embedding_size")
                       if getattr(self.backbone.config, key, None) is not None), None)
        if hidden is None:
            raise ValueError("RiNALMo configuration does not declare its hidden size")
        self.dropout = nn.Dropout(dropout)
        self.classifier = nn.Linear(int(hidden), 1)

    def forward(self, input_ids: torch.Tensor, attention_mask: torch.Tensor,
                special_tokens_mask: torch.Tensor | None = None,
                **kwargs: Any) -> dict[str, torch.Tensor]:
        output = self.backbone(input_ids=input_ids, attention_mask=attention_mask, **kwargs)
        mask = attention_mask.bool()
        if special_tokens_mask is not None:
            mask = mask & ~special_tokens_mask.bool()
        weights = mask.to(output.last_hidden_state.dtype).unsqueeze(-1)
        pooled = (output.last_hidden_state * weights).sum(1) / weights.sum(1).clamp_min(1.0)
        return {"logits": self.classifier(self.dropout(pooled)).squeeze(-1)}


def normalize_rna(sequence: str) -> str:
    normalized = sequence.strip().upper().replace("T", "U")
    bad = sorted(set(normalized) - set("ACGU"))
    if bad:
        raise ValueError(f"Unsupported RNA characters: {bad!r}")
    if not normalized:
        raise ValueError("RNA sequence must not be empty")
    return normalized
