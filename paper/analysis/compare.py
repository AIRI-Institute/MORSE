"""Select and compare completed fresh runs without declaring trajectory equality."""
from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path
import tempfile

import numpy as np

from ..io import DATA, HERE, read_json, registry, runs, sha256, write_csv, write_json

FIG5_METRICS = ["ribonn_mean_te", "translation_efficiency_score", "ribodecode_translation_log"]
ALPHAS = [0, .25, .5, .75, 1]


def _finite_metrics(pool, metrics):
    for candidate in pool:
        for metric in metrics:
            value = candidate["raw"].get(metric)
            if not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"Missing/nonfinite selection output {metric} for {candidate['id']}")


def select(record, pool, normalizations, reference):
    """The archived per-study specialist, hypervolume, and tie-breaking rules."""
    if not pool:
        return {}
    if record["figure"] == "F3A":
        _finite_metrics([pool[0]], ["mfe_per_nt", "csi_percentile"])
        bm, bc = pool[0]["raw"]["mfe_per_nt"], pool[0]["raw"]["csi_percentile"]
        eligible = [c for c in pool if c["valid"]]
        _finite_metrics(eligible, ["mfe_per_nt", "csi_percentile"])
        selected = {}
        for role in ["mfe_per_nt", "csi_percentile", "balanced"]:
            def allowed(c):
                m, s = c["raw"]["mfe_per_nt"], c["raw"]["csi_percentile"]
                return s >= bc if role == "mfe_per_nt" else m <= bm if role == "csi_percentile" else s >= bc and m <= bm
            def key(c):
                m = (bm - c["raw"]["mfe_per_nt"]) / .5
                s = (c["raw"]["csi_percentile"] - bc) / 100
                return (m, s, -c["index"]) if role == "mfe_per_nt" else (s, m, -c["index"]) if role == "csi_percentile" else (min(m, s), m + s, -c["index"])
            choices = [c for c in eligible if allowed(c)]
            if choices:
                selected[role] = max(choices, key=key)
        return selected
    pool = [c for c in pool if c["valid"] and c["raw"].get("feasibility_pass", True)]
    if not pool:
        return {}
    if record["figure"] == "F4":
        _finite_metrics(pool, FIG5_METRICS)
        reference = np.asarray(reference)
        if reference.ndim != 2 or reference.shape[1] != 3 or not len(reference) or not np.all(np.isfinite(reference)):
            raise ValueError("Fig5 selection requires a finite nonempty three-column empirical reference")
        x = np.asarray([[c["raw"][k] for k in FIG5_METRICS] for c in pool])
        ranks = np.column_stack([np.searchsorted(reference[:, j], x[:, j], side="right") / len(reference) for j in range(3)])
        selected = {key: pool[max(range(len(pool)), key=lambda i: (x[i, j], min(np.delete(ranks[i], j)), np.mean(np.delete(ranks[i], j)), -pool[i]["index"]))] for j, key in enumerate(FIG5_METRICS)}
        selected["maximin"] = pool[max(range(len(pool)), key=lambda i: (min(ranks[i]), np.mean(ranks[i]), -pool[i]["index"]))]
        return selected
    norm = normalizations[record["figure"]]
    metrics = norm["metrics"]
    _finite_metrics(pool, metrics)
    x = np.asarray([[c["raw"][k] for k in metrics] for c in pool])
    z = (x * norm["cost_signs"] - norm["cost_ideal"]) / norm["cost_span"]
    hv = np.prod(np.maximum(0, norm["reference"] - z), axis=1)
    selected = {key: pool[min(range(len(pool)), key=lambda i: (z[i, j], pool[i]["id"]))] for j, key in enumerate(metrics)}
    selected["best_hv"] = pool[min(range(len(pool)), key=lambda i: (-hv[i], pool[i]["id"]))]
    return selected


def _read_rows(path):
    with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as stream:
        for line in stream:
            if line.strip():
                yield json.loads(line)


def _paper_metric_aliases(record):
    """Relate new scorer.metric names to the frozen observation columns."""
    config = record.get("config")
    if not config:
        return {}
    migrated = read_json(HERE / "provenance/schema4_migration.json")["configs"]
    return migrated.get(Path(config).stem, {}).get("metric_names", {})


def _historical_output_values(values, aliases):
    result = {}
    for name, value in values.items():
        target = aliases.get(name, name)
        if target in result:
            raise ValueError(f"Ambiguous fresh paper output mapping for {target}")
        result[target] = value
    return result


def _protocol_mode(directory, record):
    plan_path = directory / "plan.json"
    if not plan_path.is_file():
        return "direct-morse-run" if (directory / "resolved_run.yaml").is_file() else "unverified-protocol"
    plan = read_json(plan_path)
    if plan.get("run_id") != record["run_id"]:
        raise ValueError(f"Plan identity differs from its run directory: {directory}")
    if plan.get("config_sha256"):
        config = directory / "resolved.yaml"
        if not config.is_file() or sha256(config) != plan["config_sha256"]:
            raise ValueError(f"Resolved configuration changed after planning: {directory}")
    if plan.get("protocol_sha256"):
        protocol = directory / "vaxpress_protocol.json"
        if not protocol.is_file() or sha256(protocol) != plan["protocol_sha256"]:
            raise ValueError(f"VaxPress protocol changed after planning: {directory}")
    if record.get("initialization_policy"):
        # Older stored launcher evidence is readable, but no historical runtime
        # configuration is used to start a new MORSE run.
        import yaml
        cfg = yaml.safe_load((directory / "resolved.yaml").read_text())
        supplied = bool(cfg["optimizer"].get("initial_population"))
        if supplied != bool(record["retained_initial_population"]):
            return "superseded-initialization-protocol"
        presets = read_json(HERE / "provenance/scorer-environments.json")
        mapping = presets["configs"][Path(record["config"]).stem]
        for scorer in cfg["scorers"]:
            policy = presets["profiles"][mapping[scorer["name"]]]["recipe"]["paper_runtime"]
            if policy["scorer"] in {"ribodecode", "ribonn", "orthrus-te-full-fit", "orthrus-te-held-out"} and scorer.get("mode") == "persistent-container":
                docker = scorer.get("docker_args", [])
                gpu = any(a == "--gpus" or a.startswith("--gpus=") for a in docker)
                if gpu != (policy["reference_device"] == "cuda"):
                    return "superseded-device-protocol"
    return plan.get("mode", "unverified-protocol")


def _read_completed(directory, record):
    """Read all event counts, retaining only the declared scientific budget."""
    vaxpress = record["condition"] == "vaxpress"
    if vaxpress:
        paths = [directory / f"alpha_{a}/evaluations.jsonl" for a in ALPHAS]
        summaries = [directory / f"alpha_{a}/summary.json" for a in ALPHAS]
    else:
        results = directory if (directory / "run_log.jsonl.gz").is_file() else directory / "results"
        paths = [results / "run_log.jsonl.gz"]
        summaries = [results / "summary.json"]
    missing = [str(p.relative_to(directory)) for p in [*paths, *summaries] if not p.is_file()]
    if missing:
        return None, {"run_id": record["run_id"], "status": "incomplete-artifacts", "missing": missing}
    pool, events, last_generation = [], 0, 0
    metric_aliases = {} if vaxpress else _paper_metric_aliases(record)
    alpha_counts = []
    for path in paths:
        local = 0
        for row in _read_rows(path):
            if vaxpress:
                local += 1
                events += 1
                # Every weight has its own 2,560-event quota. A longer first
                # sweep must never displace observations from a later weight.
                if local > 2560:
                    continue
                pool.append({"id": f'alpha{float(row["alpha"])}_event{local}', "sequence": row["cds"],
                             "sequence_scope": "coding_sequence", "index": local, "valid": True,
                             "raw": {k: v for k, v in row.items() if isinstance(v, (int, float))}})
            elif row.get("type", "evaluation") == "evaluation":
                events += 1
                generation = int(row["generation"])
                last_generation = max(last_generation, generation)
                if record["figure"] == "F3A" and events > 12800:
                    continue
                if record.get("cap") is not None and generation > record["cap"]:
                    continue
                pool.append({"id": row["candidate_id"], "sequence": row["sequence"], "sequence_scope": "full_transcript",
                             "index": row["eval_index"], "valid": bool(row.get("valid", False)),
                             "raw": _historical_output_values(
                                 {k[7:]: v for k, v in row.items() if k.startswith("output_") and isinstance(v, (int, float, bool))},
                                 metric_aliases)})
            elif row.get("type") == "population":
                last_generation = max(last_generation, int(row["generation"]))
        if vaxpress:
            alpha_counts.append(local)
    summary_values = [read_json(p) for p in summaries]
    if vaxpress:
        for alpha, actual, summary in zip(ALPHAS, alpha_counts, summary_values):
            if summary.get("seed") != record["seed"] or float(summary.get("alpha", -1)) != alpha:
                raise ValueError(f"VaxPress summary seed/weight differs from the frozen recipe in {directory}")
            if int(summary.get("events", -1)) != actual:
                raise ValueError(f"VaxPress summary/event count differs in {directory}")
        budget_complete = all(n >= 2560 for n in alpha_counts)
    else:
        optimizer = summary_values[0].get("optimizer", {})
        if optimizer.get("seed") != record["seed"]:
            raise ValueError(f"Completed run seed differs from the frozen recipe in {directory}")
        if int(optimizer.get("generations", -1)) != last_generation:
            raise ValueError(f"Completed generation count differs between log and summary in {directory}")
        budget_complete = len(pool) >= 12800 if record["figure"] == "F3A" else True
    mode = _protocol_mode(directory, record)
    provenance_files = [directory / "plan.json", directory / "resolved.yaml", directory / "vaxpress_protocol.json", directory / "results/run_provenance.json", directory / "resolved_run.yaml", directory / "run_provenance.json"]
    sources = [*paths, *summaries, *(p for p in provenance_files if p.is_file())]
    completion_status = "completed" if vaxpress else summary_values[0].get("status", "completed")
    details = {"run_id": record["run_id"], "status": completion_status, "protocol_mode": mode,
               "recorded_events_total": events, "in_budget_events": len(pool), "ignored_out_of_budget_events": events - len(pool),
               "expected_in_budget_events": record["recorded_evaluations"], "last_generation": last_generation or None,
               "expected_last_generation": record["generations"], "event_budget_complete": budget_complete,
               "alpha_event_counts": alpha_counts if vaxpress else None,
               "source_files": [{"path": str(p), "sha256": sha256(p)} for p in sources]}
    if not budget_complete:
        details["comparison_mode"] = "incomplete-evaluation-budget"
    elif mode == "diagnostic-generation-override":
        details["comparison_mode"] = "diagnostic-generation-override"
    elif record.get("cap") is not None and last_generation > record["cap"]:
        details["comparison_mode"] = "generation-budget-exceeded"
    else:
        details["comparison_mode"] = mode
    return pool, details


def collect(root, output):
    root, output = Path(root).expanduser().resolve(), Path(output).expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Refusing to replace existing comparison evidence: {output}; choose another output directory")
    registered = registry()["runs"]
    found, summaries, skipped, preselected = {}, {}, [], {}
    norms = read_json(DATA / "selection_normalization.json")
    for record in registered:
        directory = root / record["run_id"]
        if not directory.is_dir():
            continue
        pool, details = _read_completed(directory, record)
        if pool is None:
            skipped.append(details)
        else:
            rid = record["run_id"]
            # Only Fig5 needs cross-run pools for its empirical reference. For
            # other studies retain the selected candidates, not every archive.
            if record["figure"] == "F4":
                found[rid] = pool
            else:
                preselected[rid] = select(record, pool, norms, None)
                found[rid] = None
            summaries[rid] = details
    if not found:
        raise ValueError(f"No completed paper runs found under {root}; completed summaries and event logs are required")
    fig5 = [r for r in registered if r["figure"] == "F4"]
    fresh_fig5 = [r for r in fig5 if r["run_id"] in found]
    reference, reference_policy = None, "not used: no fresh Fig5 runs"
    full_fig5 = len(fig5) == 5 and len(fresh_fig5) == 5 and all(summaries[r["run_id"]]["comparison_mode"] == "paper-protocol" for r in fig5)
    if full_fig5:
        unique = {c["sequence"]: [c["raw"][k] for k in FIG5_METRICS] for r in fig5 for c in found[r["run_id"]] if c["valid"]}
        if not unique:
            raise ValueError("Fresh Fig5 runs contain no valid scored transcripts")
        reference = np.sort(np.asarray(list(unique.values())), axis=0)
        reference_policy = "pooled unique transcripts from all five completed fresh Fig5 paper-protocol archives"
    elif fresh_fig5:
        with np.load(DATA / "Figure5_empirical_reference.npz") as source:
            reference = source["values"]
        reference_policy = "frozen historical CDF for partial, diagnostic, or unverified Fig5 comparison"
    archived = {r["run_id"]: {s["role"]: s["candidate"] for s in r["selected"]} for r in runs()}
    comparisons, selected_payload, fasta = [], {}, []
    for record in registered:
        rid = record["run_id"]
        if rid not in found:
            continue
        selected = preselected[rid] if rid in preselected else select(record, found[rid], norms, reference)
        selected_payload[rid] = selected
        expected = {s["role"]: s for s in record["expected_selected"]}
        ids = {role: candidate["id"] for role, candidate in selected.items()}
        sequence_matches, metric_matches = [], []
        for role, candidate in selected.items():
            previous = archived[rid][role]
            expected_sequence = previous["coding_sequence"] if candidate.get("sequence_scope") == "coding_sequence" else previous["sequence"]
            same_sequence = candidate["sequence"] == expected_sequence
            sequence_matches.append(same_sequence)
            candidate["sequence_sha256"] = hashlib.sha256(candidate["sequence"].encode()).hexdigest()
            old = expected[role]
            for metric, expected_value in old["raw"].items():
                actual = candidate["raw"].get(metric)
                same_value = isinstance(actual, (int, float)) and math.isfinite(actual) and math.isclose(actual, expected_value, rel_tol=1e-7, abs_tol=1e-9)
                metric_matches.append(same_value)
                comparisons.append({"run_id": rid, "mode": summaries[rid]["comparison_mode"], "role": role,
                    "candidate_id": candidate["id"], "expected_candidate_id": old["candidate_id"], "sequence_matches": same_sequence,
                    "metric": metric, "actual": actual, "expected": expected_value,
                    "difference": actual - expected_value if isinstance(actual, (int, float)) else None, "within_tolerance": same_value})
            fasta.append(f'>{rid} role={role} candidate={candidate["id"]} scope={candidate.get("sequence_scope", "full_transcript")}')
            fasta.extend(candidate["sequence"][i:i + 80] for i in range(0, len(candidate["sequence"]), 80))
        summaries[rid].update(selected=ids, missing_roles=sorted(set(expected) - set(selected)),
            all_selected_ids_match=ids == {k: v["candidate_id"] for k, v in expected.items()},
            all_selected_sequences_match=len(selected) == len(expected) and all(sequence_matches),
            all_expected_metrics_match=len(selected) == len(expected) and bool(metric_matches) and all(metric_matches))
    result = {"status": "collected", "runs_found": len(found), "runs_expected": len(registered), "skipped_incomplete_runs": skipped,
        "fig5_reference_policy": reference_policy, "runs": list(summaries.values()),
        "numeric_tolerance": {"relative": 1e-7, "absolute": 1e-9},
        "interpretation": "Candidate IDs alone do not establish equality. Sequence identity, numeric differences, budgets and protocol modes are reported separately; diagnostic runs are not full paper reproduction."}
    # Publish a complete comparison directory only after collection succeeds.
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".morse-comparison-", dir=output.parent) as temporary:
        stage = Path(temporary) / "comparison"
        stage.mkdir()
        write_json(stage / "fresh_selected_candidates.json", selected_payload)
        if comparisons:
            write_csv(stage / "fresh_score_comparison.csv", comparisons)
        (stage / "fresh_selected_sequences.fasta").write_text("\n".join(fasta) + ("\n" if fasta else ""))
        write_json(stage / "fresh_comparison_receipt.json", result)
        if output.exists():
            output.rmdir()
        stage.rename(output)
    return result
