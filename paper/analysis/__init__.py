"""Recompute paper tables and verify selections against the complete archives."""
from __future__ import annotations

from collections import defaultdict
import gzip
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys

from ..io import DATA, HERE, read_json, registry, rows, runs, write_csv, write_json


def summarize(values):
    values = [float(v) for v in values]
    return {"n": len(values), "mean": statistics.mean(values),
            "sample_sd": statistics.stdev(values) if len(values) > 1 else None}


def tables(output):
    """Generate numeric CSV tables from immutable individual observations."""
    import yaml
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    records = runs()
    result = []
    configurations = []
    environments = read_json(HERE / "provenance/scorer-environments.json")
    for r in registry()["runs"]:
        if r["config"]:
            cfg = yaml.safe_load((HERE / r["config"]).read_text())
            for scorer in cfg["scorers"]:
                profile = environments["configs"][Path(r["config"]).stem][scorer["name"]]
                logical_scorer = environments["profiles"][profile]["logical_scorer"]
                for metric, policy in scorer["metrics"].items():
                    policy = dict(policy)
                    value_type = policy.pop("value_type", "float")
                    # Keep the manuscript table's established terminology and
                    # numeric types while run YAML uses the simpler interface.
                    if policy["role"] == "optimize":
                        policy["role"] = "objective"
                    configurations.append({"study": r["study"], "condition": r["condition"], "seed": r["seed"],
                        "scorer": logical_scorer, "output": metric,
                        "scope": json.dumps(scorer.get("scope"), sort_keys=True), "type": value_type, **policy})
    write_csv(output / "SupplementaryTable1_scorers.csv", configurations)
    protocols = [{k: r.get(k) for k in ["study", "run_id", "condition", "seed", "cap", "generations", "attempted_evaluations", "recorded_evaluations"]} for r in registry()["runs"]]
    write_csv(output / "SupplementaryTable2_protocols.csv", protocols)
    bootstrap = rows(DATA / "SupplementaryFigure3_Mauger_U/paired_bootstrap_20261016/comparison_results.csv")
    write_csv(output / "SupplementaryTable3_paired_bootstrap.csv", bootstrap)
    policy = read_json(DATA / "Sequence_feasibility_manufacturing_strategy_v1.json")
    write_json(output / "SupplementaryTable4_feasibility.json", policy)
    f2 = [r for r in records if r["figure"] == "F2"]
    endpoint = []
    for r in f2:
        c = next(c for c in r["curves"] if c["generation"] == 220)
        endpoint.append({"run_id": r["run_id"], "condition": r["condition"], "seed": r["seed"], **c})
    write_csv(output / "SupplementaryTable5_generation220.csv", endpoint)
    grouped = defaultdict(list)
    for row in endpoint:
        for metric, value in row.items():
            if metric not in {"seed", "generation"} and isinstance(value, (int, float)):
                grouped[(row["condition"], metric)].append(value)
    write_csv(output / "SupplementaryTable5_summary.csv", [dict(condition=c, metric=k, **summarize(v)) for (c, k), v in grouped.items()])
    utr = []
    for r in f2:
        for distribution in r["u_distributions"]:
            utr.append({"run_id": r["run_id"], "condition": r["condition"], "seed": r["seed"],
                        "stage": distribution["stage"], "generation": distribution["generation"], "region": distribution["region"],
                        "mean_max_u_run": statistics.mean(distribution["max_u_runs"]),
                        "fraction_with_u_ge10": sum(x > 0 for x in distribution["long_run_counts"]) / 64,
                        "mean_u_ge10_per_kb": statistics.mean(distribution["long_runs_per_kb"])})
    write_csv(output / "SupplementaryTable6_uridine.csv", utr)
    numeric, ham, utrs = [], [], []
    for r in records:
        base = r["baseline"]
        for selection in r["selected"]:
            c = selection["candidate"]
            common = {"run_id": r["run_id"], "figure": r["figure"], "condition": r["condition"], "seed": r["seed"], "role": selection["role"], "candidate_id": c["candidate_id"]}
            for metric, value in c["raw"].items():
                if isinstance(value, (int, float)) and metric in base["raw"]:
                    numeric.append({**common, "metric": metric, "initial": base["raw"][metric], "selected": value, "change": value - base["raw"][metric]})
            bcds, cds = base["cds"], c["cds"]
            changed = sum(bcds[i:i+3] != cds[i:i+3] for i in range(0, len(cds), 3))
            ham.append({**common, "body_codons": len(cds) // 3, "changed_codons": changed, "codon_hamming_fraction": changed / (len(cds) // 3)})
            if r["figure"] == "F3B":
                utrs.append({**common, "five_prime_utr": c["five_prime_utr"], "three_prime_utr": c["three_prime_utr"]})
    write_csv(output / "SupplementaryTable7_native_comparisons.csv", numeric)
    write_csv(output / "SupplementaryTable8_additional_scores.csv", rows(DATA / "Figure4_additional_scores.csv"))
    write_csv(output / "SupplementaryTable9_UTRs.csv", utrs)
    write_csv(output / "SupplementaryTable10_codon_hamming.csv", ham)
    write_csv(output / "Run_termination_records.csv", protocols)
    # All summaries preserve the independent-run unit; roles within a run overlap.
    grouped = defaultdict(list)
    for row in numeric:
        grouped[(row["figure"], row["condition"], row["role"], row["metric"])].append(row["selected"])
    write_csv(output / "Native_score_summary.csv", [dict(figure=f, condition=c, role=r, metric=m, **summarize(v)) for (f, c, r, m), v in grouped.items()])
    return [p.name for p in sorted(output.glob("*")) if p.is_file()]


def _events(record, raw_root):
    import numpy as np
    sources = [p for p in record["sources"] if p.endswith((".jsonl", ".jsonl.gz"))]
    candidates, generations = [], set()
    for source in sources:
        path = Path(raw_root) / source
        with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as stream:
            event = 0
            for line in stream:
                r = json.loads(line)
                if record["figure"] == "F3A" and "output_mfe_per_nt" not in r:
                    event += 1
                    candidates.append({"id": r.get("candidate_id",f'alpha{float(r.get("alpha",0))}_event{event}'), "sequence": r["cds"], "index": r.get("event_index",event),
                                       "valid": True, "raw": {k: v for k, v in r.items() if isinstance(v, (int, float))}})
                elif r.get("type", "evaluation") == "evaluation":
                    event += 1
                    candidates.append({"id": r["candidate_id"], "sequence": r["sequence"], "index": r.get("eval_index", event),
                                       "valid": r.get("valid", False), "raw": {k[7:]: v for k, v in r.items() if k.startswith("output_") and isinstance(v, (int, float, bool))}})
                elif r.get("type") == "population":
                    generations.add(r["generation"])
    if record["figure"] == "F3A":
        candidates = candidates[:12800]
    return candidates, generations


def audit_raw(raw_root):
    """Recompute all selected IDs from full in-budget archived evaluations."""
    import numpy as np
    from ..verification import verify
    verify(raw=True, raw_root=raw_root)
    registered = registry()["runs"]
    pools = {r["run_id"]: _events(r, raw_root)[0] for r in registered}
    metrics5 = ["ribonn_mean_te", "translation_efficiency_score", "ribodecode_translation_log"]
    unique = {}
    for r in registered:
        if r["figure"] == "F4":
            for c in pools[r["run_id"]]:
                if c["valid"]:
                    unique[c["sequence"]] = [c["raw"][k] for k in metrics5]
    if len(unique) != 104239:
        raise ValueError(f"Fig5 CDF population mismatch: {len(unique)} vs 104239")
    reference = np.sort(np.asarray(list(unique.values())), axis=0)
    normalizations = read_json(DATA / "selection_normalization.json")
    from .compare import select
    checks = []
    for r in registered:
        pool = pools[r["run_id"]]
        if len(pool) != r["recorded_evaluations"]:
            raise ValueError(f"Recorded event mismatch: {r['run_id']} {len(pool)} vs {r['recorded_evaluations']}")
        selected = {role:candidate["id"] for role,candidate in select(r,pool,normalizations,reference).items()}
        expected = {c["role"]: c["candidate_id"] for c in r["expected_selected"]}
        if selected != expected:
            raise ValueError(f"Recomputed selections differ for {r['run_id']}: actual={selected}; frozen={expected}")
        checks.append({"run_id": r["run_id"], "recorded_evaluations": r["recorded_evaluations"], "selected": selected, "status": "pass"})
    return {"unique_fig5_transcripts": len(unique), "checks": checks, "all_selected_candidate_ids_recomputed": True}


def analyze(output, *, raw=False, raw_root=None, bootstrap=False,runs_root=None):
    from ..verification import verify
    verify()
    output = Path(output).resolve()
    result = {"status": "pass", "tables": tables(output), "source": "frozen individual numeric observations", "new_model_inference": False}
    from .mauger import calculate
    result["mauger_reconstruction"]=calculate(output,permutations=True)
    if raw:
        result["archive_audit"] = audit_raw(raw_root or HERE / "artifacts/raw")
    if runs_root:
        from .compare import collect
        result["fresh_run_comparison"]=collect(runs_root,output/"fresh")
    if bootstrap:
        subprocess.run([sys.executable, str(HERE / "analysis/bootstrap.py"), "--input-dir", str(DATA / "SupplementaryFigure3_Mauger_U/paired_bootstrap_20261016/inputs"), "--output-dir", str(output / "bootstrap"), "--resamples", "200000", "--seed", "20261016"], check=True)
        result["bootstrap"] = "200000 paired reconstructing resamples, seed20261016"
        expected_bootstrap=rows(DATA/"SupplementaryFigure3_Mauger_U/paired_bootstrap_20261016/comparison_results.csv")
        actual_bootstrap=rows(output/"bootstrap/comparison_results.csv")
        if len(actual_bootstrap)!=len(expected_bootstrap):raise ValueError("Bootstrap comparison count mismatch")
        for observed,expected in zip(actual_bootstrap,expected_bootstrap):
            for key in ["difference","bca_95_lower","bca_95_upper","bca_bonferroni_lower","bca_bonferroni_upper"]:
                if not math.isclose(float(observed[key]),float(expected[key]),abs_tol=1e-12):raise ValueError(f"Recomputed bootstrap mismatch: {observed['component']}:{key}")
        result["bootstrap_intervals_match_archived"]=True
    write_json(output / "analysis_receipt.json", result)
    return result
