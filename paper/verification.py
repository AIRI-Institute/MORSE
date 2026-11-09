"""Content and numerical verification; no scorers or network are invoked."""
from __future__ import annotations

from collections import Counter
import math
from pathlib import Path

from .io import DATA, HERE, read_json, registry, runs, sha256


def verify(*, raw=False, raw_root=None):
    manifest = read_json(HERE / "provenance/bundle_manifest.json")
    missing, mismatches = [], []
    for name, expected in manifest["files"].items():
        path = HERE / name
        if not path.is_file():
            missing.append(name)
        elif sha256(path) != expected["sha256"] or path.stat().st_size != expected["bytes"]:
            mismatches.append(name)
    if missing or mismatches:
        raise ValueError(f"Frozen bundle changed: missing={missing}; checksum mismatch={mismatches}")
    records = runs()
    groups = Counter((r["figure"], r["condition"]) for r in records)
    if len(records) != 40 or len(groups) != 8 or set(groups.values()) != {5}:
        raise ValueError("Expected exactly 40 independent runs/sweeps in 8 groups of 5")
    expected = {r["run_id"]: r for r in registry()["runs"]}
    for run in records:
        frozen = expected[run["run_id"]]
        for key in ["generations", "recorded_evaluations"]:
            if run.get(key) != frozen.get(key):
                raise ValueError(f"Run record mismatch {run['run_id']}:{key}")
        if run["figure"] in {"F2", "F3B", "F3C"}:
            for selected in run["selected"]:
                candidate = selected["candidate"]
                if not candidate["raw"].get("feasibility_pass", False):
                    raise ValueError(f"Infeasible selected candidate in {run['run_id']}")
    # Native individual values are summarized only after panel-specific scaling.
    for filename, section in [("figure_display_coordinates_0105.json", "figure3"),
                              ("figure_display_coordinates_0105.json", "figure5")]:
        d = read_json(DATA / filename)[section]
        for profile in d["radar_profiles"]:
            values = profile["normalized_values"]
            for j, observed_mean in enumerate(profile.get("mean", profile.get("display_mean"))):
                recomputed = sum(row[j] for row in values) / len(values)
                if not math.isclose(recomputed, observed_mean, abs_tol=1e-12):
                    raise ValueError(f"Radar aggregation mismatch: {section}/{profile['role']}")
    raw_count = 0
    if raw:
        root = Path(raw_root or HERE / "artifacts/raw")
        for name, digest in read_json(HERE / "provenance/raw_sources.json")["sources"].items():
            path = root / name
            if not path.is_file():
                raise FileNotFoundError(f"Missing archived run source {path}; extract the paper raw-source release bundle")
            if sha256(path) != digest:
                raise ValueError(f"Raw source checksum mismatch: {name}")
            raw_count += 1
    return {"status": "pass", "frozen_files": len(manifest["files"]), "independent_runs_or_sweeps": 40,
            "conditions": 8, "repeats_per_condition": 5, "raw_files_verified": raw_count,
            "model_inference": False, "optimization_rerun": False}
