from __future__ import annotations

import contextlib
import json
import os
import sys
import tempfile
import traceback
from pathlib import Path
from typing import Any

import numpy as np


ENV_TO_CONDITION_FILE = {
    "HEK293T": "HEK293T_normal_10552_RPKM_HEK293T.npz",
    "A549": "A549_normal_10552_RPKM_A549.npz",
    "HeLa": "HeLa_normal_10552_RPKM_HeLa.npz",
    # The author's Hela/HeLa files have identical registered bytes. Use the
    # canonical spelling so external preparation also works on macOS/APFS.
    "Hela": "HeLa_normal_10552_RPKM_HeLa.npz",
}


class RiboDecodeScorer:
    def __init__(self) -> None:
        import torch
        from TranslationModel.inference import InferenceModel, base_dir

        self.torch = torch
        self.base_dir = Path(base_dir)
        self.output_prefix = os.environ.get("RIBODECODE_OUTPUT_PREFIX", "ribodecode")
        self.batch_size = int(os.environ.get("RIBODECODE_BATCH_SIZE", "256"))
        model_config = os.environ.get(
            "RIBODECODE_MODEL_CONFIG", str(self.base_dir / "model_config.json")
        )
        best_model = os.environ.get(
            "RIBODECODE_BEST_MODEL", str(self.base_dir / "best_model.p")
        )
        condition_path = self._condition_path()

        with contextlib.redirect_stdout(sys.stderr):
            self.model = InferenceModel(model_config, best_model)

        self.model.predictorModel.eval()
        self.condition_array = self._load_condition_array(condition_path)
        self.mrna_count = np.array([[4.5]], dtype=np.float32)
        print(
            "RiboDecode scorer ready "
            f"device={self.model.device} condition={condition_path} batch_size={self.batch_size}",
            file=sys.stderr,
            flush=True,
        )

    def score(self, raw_sequences: list[str]) -> list[dict[str, Any]]:
        sequences = [
            _normalize_coding_sequence(sequence) for sequence in raw_sequences
        ]
        log_scores = self._predict_log_scores(sequences)
        outputs: list[dict[str, Any]] = []
        if len(log_scores) != len(sequences):
            raise RuntimeError(f"Expected {len(sequences)} scores, got {len(log_scores)}")
        for log_score in log_scores:
            log_score = float(log_score)
            translation = float(np.expm1(np.clip(log_score, -60.0, 60.0)) / 5.0)
            outputs.append(
                {
                    f"{self.output_prefix}_translation_log": log_score,
                    f"{self.output_prefix}_translation": translation,
                }
            )
        return outputs

    def _predict_log_scores(self, sequences: list[str]) -> list[float]:
        scores: list[float] = []
        torch = self.torch
        condition_tensor = torch.tensor(
            self.condition_array, dtype=torch.float32
        ).to(self.model.device)
        count_tensor = torch.tensor(self.mrna_count, dtype=torch.float32).to(self.model.device)

        with torch.no_grad():
            for start in range(0, len(sequences), self.batch_size):
                batch = sequences[start : start + self.batch_size]
                cds_oh = [self.model.my_process_line_for_cds(sequence, 4500) for sequence in batch]
                input_cds = torch.tensor(np.vstack(cds_oh), dtype=torch.float32)
                cds_sequence = torch.transpose(input_cds, 2, 1)[:, :4, :].to(self.model.device)
                condition_batch = torch.cat([condition_tensor] * len(batch), dim=0)
                count_batch = torch.cat([count_tensor] * len(batch), dim=0)
                prediction = self.model.predictorModel(cds_sequence, condition_batch, count_batch)
                scores.extend(prediction.float().cpu().numpy()[:, 0].tolist())
        return scores

    def _condition_path(self) -> Path:
        csv_path = os.environ.get("RIBODECODE_ENV_CSV")
        if csv_path:
            return _condition_npz_from_csv(Path(csv_path))

        env = os.environ.get("RIBODECODE_ENV", "HEK293T")
        if env not in ENV_TO_CONDITION_FILE:
            choices = ", ".join(sorted(ENV_TO_CONDITION_FILE))
            raise ValueError(f"Unsupported RIBODECODE_ENV={env!r}; expected one of {choices}")
        return self.base_dir / "conditions" / ENV_TO_CONDITION_FILE[env]

    def _load_condition_array(self, condition_path: Path) -> np.ndarray:
        npzfile = np.load(condition_path, allow_pickle=True)
        condition = np.expand_dims(npzfile["mRNA"].astype(np.float32), 0)
        condition = np.log1p(condition * 5.0)
        return condition[:, :10552]


def main() -> int:
    try:
        scorer = RiboDecodeScorer()
    except Exception:
        traceback.print_exc(file=sys.stderr)
        return 1

    for line in sys.stdin:
        if not line.strip():
            continue
        try:
            request = json.loads(line)
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
            response = {
                "request_id": request.get("request_id"),
                "outputs": scorer.score(sequences),
            }
        except Exception as exc:
            traceback.print_exc(file=sys.stderr)
            response = {"request_id": _request_id_from_line(line), "error": str(exc)}
        print(json.dumps(response, separators=(",", ":")), flush=True)
    return 0


def _normalize_coding_sequence(sequence: str) -> str:
    coding_sequence = "".join(sequence.upper().split()).replace("U", "T")
    if len(coding_sequence) > 4500:
        raise ValueError(
            "RiboDecode supports coding sequences up to 4500 nt; "
            f"got {len(coding_sequence)} nt"
        )
    if len(coding_sequence) % 3 != 0:
        raise ValueError(
            "Coding-sequence length must be divisible by 3; "
            f"got {len(coding_sequence)} nt"
        )
    invalid = sorted(set(coding_sequence) - {"A", "C", "G", "T", "N"})
    if invalid:
        raise ValueError(
            "Coding sequence contains unsupported bases: " + "".join(invalid)
        )
    return coding_sequence


def _condition_npz_from_csv(csv_path: Path) -> Path:
    import pandas as pd

    df = pd.read_csv(csv_path, index_col=0)
    if df.shape != (10552, 1):
        raise ValueError(
            "Custom RiboDecode CSV must contain exactly 10,552 genes and one RPKM column; "
            f"got {df.shape}: {csv_path}"
        )
    if not df.index.is_unique:
        raise ValueError(f"Custom RiboDecode CSV contains duplicate gene IDs: {csv_path}")
    template_path = Path(os.environ.get("RIBODECODE_ENV_TEMPLATE", "/assets/env_file.csv"))
    if template_path.is_file():
        template = pd.read_csv(template_path, index_col=0)
        if list(df.index.astype(str)) != list(template.index.astype(str)):
            raise ValueError(
                "Custom RiboDecode CSV gene IDs or their order differ from the official template"
            )
    mrna_counts = np.asarray(df.iloc[:, 0], dtype=np.float32)
    if not np.isfinite(mrna_counts).all() or (mrna_counts < 0).any():
        raise ValueError("Custom RiboDecode RPKM values must be finite and non-negative")
    temp_dir = Path(tempfile.mkdtemp(prefix="ribodecode_condition_"))
    npz_path = temp_dir / "condition.npz"
    np.savez(npz_path, mRNA=mrna_counts)
    return npz_path


def _request_id_from_line(line: str) -> Any:
    try:
        return json.loads(line).get("request_id")
    except Exception:
        return None


if __name__ == "__main__":
    raise SystemExit(main())
