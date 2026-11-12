from __future__ import annotations

import contextlib
import copy
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
        self.cell_context = os.environ.get("RIBODECODE_ENV", "HEK293T")
        csv_path = os.environ.get("RIBODECODE_ENV_CSV")
        self.context_file = str(Path(csv_path).resolve()) if csv_path else None
        if self.context_file is None:
            self.cell_context = _validate_cell_context(self.cell_context)
        self._condition_cache: dict[tuple[Any, ...], np.ndarray] = {}
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
        self._condition_cache[self._condition_key(self.cell_context, self.context_file)] = self.condition_array
        self.mrna_count = np.array([[4.5]], dtype=np.float32)
        print(
            "RiboDecode scorer ready "
            f"device={self.model.device} condition={condition_path} batch_size={self.batch_size}",
            file=sys.stderr,
            flush=True,
        )

    def handle_request(self, request: dict[str, Any]) -> dict[str, Any]:
        if not isinstance(request, dict):
            raise ValueError("Request must be a JSON object")
        unknown = set(request) - {"request_id", "sequences", "sequence_contexts", "parameters", "params"}
        if unknown:
            raise ValueError("Unknown request fields: " + ", ".join(sorted(unknown)))
        sequences = request.get("sequences")
        if not isinstance(sequences, list) or not all(isinstance(sequence, str) for sequence in sequences):
            raise ValueError("Request must define a list of sequences")
        overrides = _request_parameters(request)
        parameters = self.resolve_parameters(overrides)
        scorer = copy.copy(self)
        scorer.batch_size = parameters["batch_size"]
        scorer.cell_context = parameters.get("cell_context")
        scorer.context_file = parameters["context_file"]
        if "cell_context" in overrides or "context_file" in overrides:
            key = self._condition_key(scorer.cell_context, scorer.context_file)
            if key not in self._condition_cache:
                self._condition_cache[key] = self._load_condition_array(scorer._condition_path())
            scorer.condition_array = self._condition_cache[key]
        return {"request_id": request.get("request_id"), "outputs": scorer.score(sequences), "parameters": parameters}

    def resolve_parameters(self, parameters: dict[str, Any]) -> dict[str, Any]:
        unknown = set(parameters) - {"cell_context", "context_file", "batch_size"}
        if unknown:
            raise ValueError("Unknown RiboDecode parameters: " + ", ".join(sorted(unknown)))
        # An explicit built-in context takes precedence over a startup CSV.
        context_file = parameters.get("context_file", None if "cell_context" in parameters else self.context_file)
        if context_file is not None:
            if not isinstance(context_file, str) or not context_file.strip():
                raise ValueError("context_file must be a non-empty path or null")
            context_file = str(Path(context_file).resolve())
        batch_size = parameters.get("batch_size", self.batch_size)
        if type(batch_size) is not int or batch_size < 1:
            raise ValueError("batch_size must be a positive integer")
        settings = {"context_file": context_file, "batch_size": batch_size}
        if context_file is None:
            settings["cell_context"] = _validate_cell_context(parameters.get("cell_context", self.cell_context))
        return settings

    @staticmethod
    def _condition_key(cell_context: str | None, context_file: str | None) -> tuple[Any, ...]:
        if context_file is None:
            return ("builtin", cell_context)
        stat = Path(context_file).stat()
        return ("csv", context_file, stat.st_mtime_ns, stat.st_size)

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
        if self.context_file:
            return _condition_npz_from_csv(Path(self.context_file))
        return self.base_dir / "conditions" / ENV_TO_CONDITION_FILE[self.cell_context]

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
            response = scorer.handle_request(request)
        except Exception as exc:
            traceback.print_exc(file=sys.stderr)
            response = {"request_id": _request_id_from_line(line), "error": str(exc)}
        print(json.dumps(response, separators=(",", ":")), flush=True)
    return 0


def _validate_cell_context(value: Any) -> str:
    if not isinstance(value, str) or value not in ENV_TO_CONDITION_FILE:
        raise ValueError("cell_context must be one of: HEK293T, A549, HeLa")
    return "HeLa" if value == "Hela" else value


def _request_parameters(request: dict[str, Any]) -> dict[str, Any]:
    if "parameters" in request and "params" in request:
        raise ValueError("Use only one of parameters and params")
    parameters = request.get("parameters", request.get("params", {}))
    if not isinstance(parameters, dict) or any(not isinstance(key, str) for key in parameters):
        raise ValueError("parameters must be a mapping with string keys")
    return parameters


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
