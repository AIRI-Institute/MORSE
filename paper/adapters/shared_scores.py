"""Experiment-owned official VaxPress addon; no third-party source patches."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
import time
from functools import lru_cache
from pathlib import Path

from Bio.Seq import Seq
from vaxpress.scoring import ScoringFunction


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_verified_module(name, path, expected):
    if sha256(path) != expected:
        raise ValueError(f"Source checksum mismatch: {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=4)
def get_evaluator(config_path):
    return SharedEvaluator(config_path)


class SharedEvaluator:
    def __init__(self, config_path):
        self.cfg = json.loads(Path(config_path).read_text())
        cfg = self.cfg
        self.mfe = load_verified_module(
            "shared_pinned_viennarna", cfg["viennarna_model_path"],
            cfg["viennarna_model_sha256"])
        self.mfe.ENABLE_PARTITION_FUNCTION = False
        self.mfe.MAX_BP_SPAN = cfg["max_bp_span"]
        csi_module = load_verified_module(
            "shared_pinned_csi", cfg["csi_model_path"], cfg["csi_model_sha256"])
        if sha256(cfg["csi_profile_path"]) != cfg["csi_profile_sha256"]:
            raise ValueError("CSI profile checksum mismatch")
        self.csi = csi_module.CsiScorer(
            csi_module.load_profile(Path(cfg["csi_profile_path"])))
        self.protein = str(Seq(cfg["initial_cds"]).translate())

    def score(self, sequence):
        cds = sequence.upper().replace("U", "T")
        cfg = self.cfg
        if (len(cds) != cfg["cds_length"] or not cds.startswith("ATG")
                or str(Seq(cds).translate()) != self.protein):
            raise ValueError("Candidate does not preserve the fixed MLANA protein")
        construct = cfg["five_prime_utr"] + cds + cfg["stop"] + cfg["three_prime_utr"]
        assert len(construct) == cfg["transcript_body_length"]
        fold = self.mfe.score_sequence(construct)
        csi = self.csi.score_sequence(cds)
        result = {"cds": cds, "mfe_per_nt": fold["mfe_per_nt"],
                  "mfe": fold["mfe"], "csi_percentile": csi.csi_percentile,
                  "csi_raw": csi.csi_raw, "csi_codons": csi.codons,
                  "transcript_body_length": len(construct)}
        assert all(math.isfinite(result[key]) for key in (
            "mfe_per_nt", "mfe", "csi_percentile", "csi_raw"))
        return result


class SharedFitness(ScoringFunction):
    name = "shared"
    description = "Pinned transcript-context MFE/nt and human CSI"
    uses_folding = False
    arguments = [
        ("weight", dict(type=float, default=1.0)),
        ("alpha", dict(type=float, default=0.5)),
        ("config", dict(type=str, required=True)),
        ("log-path", dict(type=str, required=True)),
    ]

    def __init__(self, weight, alpha, config, log_path, _length_cds):
        if weight != 1.0 or not 0 <= alpha <= 1:
            raise ValueError("Expected unit addon weight and alpha in [0, 1]")
        self.alpha, self.config, self.log_path = alpha, config, log_path
        cfg = json.loads(Path(config).read_text())
        if _length_cds != cfg["cds_length"]:
            raise ValueError("VaxPress input must be CDS including ATG, without stop")

    def score(self, seqs):
        evaluator = get_evaluator(self.config)
        begin = time.monotonic()
        records = [evaluator.score(seq) for seq in seqs]
        cfg = evaluator.cfg
        mfe = [self.alpha * (-row["mfe_per_nt"] / cfg["mfe_fitness_scale"])
               for row in records]
        csi = [(1 - self.alpha) * row["csi_percentile"] / cfg["csi_fitness_scale"]
               for row in records]
        elapsed = time.monotonic() - begin
        with Path(self.log_path).open("a") as handle:
            for index, row in enumerate(records):
                handle.write(json.dumps(dict(row, alpha=self.alpha,
                    fitness=mfe[index] + csi[index], index_in_batch=index,
                    batch_size=len(records), batch_scoring_seconds=elapsed)) + "\n")
        return {"shared_mfe": mfe, "shared_csi": csi}, {
            "shared_mfe_per_nt": [row["mfe_per_nt"] for row in records],
            "shared_csi_percentile": [row["csi_percentile"] for row in records]}
