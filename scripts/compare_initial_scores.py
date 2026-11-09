#!/usr/bin/env python3
"""Audit the same 64 initial transcripts against archived score observations.

No optimization/inference is run. Matching uses full transcript sequences, never
candidate numbering. Only generation one of JSONL archives is read. The optional
checkpoint fallback deserializes a trusted local pickle and may read a large file.
Exit codes: 0 pass, 1 mismatch, 2 invalid input, 3 partial/unavailable evidence.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import pickle
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from paper.analysis.compare import _historical_output_values, _paper_metric_aliases, _read_rows
from paper.io import HERE, read_json, registry

# Fixed declared defaults. Overrides require a named metric in --tolerances JSON;
# the report retains both the selected tolerance and its source. Nothing adapts
# a tolerance to observed differences.
TOLERANCES = {
    "mfe": (1e-5, 0.0), "mfe_per_nt": (1e-7, 0.0),
    "cds_mfe": (1e-5, 0.0), "cds_mfe_per_nt": (1e-7, 0.0),
    "translation_efficiency_score": (1e-5, 1e-5),
    "half_life_score": (1e-5, 1e-5), "orthrus_te_fullfit": (1e-5, 1e-5),
    "ribodecode_translation_log": (1e-5, 1e-5), "ribodecode_translation": (1e-4, 1e-5),
    "ribonn_mean_te": (1e-5, 1e-5), "ribonn_model_sd": (1e-5, 1e-5),
    "rinalmo_window_logit_p95": (1e-5, 1e-5),
    **{name: (1e-12, 1e-12) for name in (
        "csi_percentile", "csi_raw", "csi_destabilizing_fraction", "csi_window_min", "csi_window_p05",
        "cai", "cai_log", "cai_min_weight", "cai_low_weight_fraction", "native_frequency_cai",
        "gc_fraction", "global_gc_fraction", "local_gc_min", "local_gc_max", "cds_gc_fraction",
        "cpg_per_kb", "upa_per_kb", "body_gc_window_min", "body_gc_window_max",
    )},
}
_EXACT_COUNTS = re.compile(r"(?:.*(?:_count|_codons|_length)|analyzed_length|scored_length|changed_codons|poly_a_length|structural_binary|max_homopolymer_[acgu]|threshold_.*_(?:observed|min|max|equals))$")
_TOP_FIELDS = {"valid", "structural_valid", "structural_binary", "structural_failures", "changed_codons", "poly_a_length", "start_codon", "stop_codon", "gc_fraction"}


def tolerances_from_file(path):
    if path is None:
        return {}
    raw = read_json(path)
    if not isinstance(raw, dict):
        raise ValueError("Tolerance file must map exact metric names to {atol, rtol}")
    result = {}
    for name, value in raw.items():
        if not isinstance(value, dict) or set(value) != {"atol", "rtol"}:
            raise ValueError(f"Metric {name} must specify exactly atol and rtol")
        if isinstance(value["atol"], bool) or isinstance(value["rtol"], bool):
            raise ValueError(f"Metric {name} tolerances must be numeric, not boolean")
        atol, rtol = float(value["atol"]), float(value["rtol"])
        if not math.isfinite(atol) or not math.isfinite(rtol) or min(atol, rtol) < 0:
            raise ValueError(f"Metric {name} has an invalid tolerance")
        field = str(name)
        if field.removeprefix("output_") not in _TOP_FIELDS:
            field = field.removeprefix("output_")
        result[field] = (atol, rtol)
    return result


def first_generation(path):
    """Stop at the first later-generation evaluation; do not scan the archive."""
    rows, consumed = [], 0
    for row in _read_rows(Path(path)):
        consumed += 1
        if row.get("type", "evaluation") != "evaluation":
            continue
        generation = int(row.get("generation", -1))
        if generation > 1 and rows:
            break
        if generation == 1:
            rows.append(row)
        if len(rows) > 10000:
            raise ValueError("Unexpectedly large generation one; expected a 64-member paper population")
    if not rows:
        raise ValueError(f"No generation-one evaluation observations in {path}")
    return rows, {"path": str(Path(path).resolve()), "records_read": consumed,
                  "generation_one_records": len(rows), "complete_file_hash_verified": False}


def _reconstruct_archive_sequences(rows, record):
    if all("sequence" in row for row in rows):
        return rows, "full transcripts recorded in the canonical archive"
    if record["figure"] != "F3A" or record["condition"] != "morse":
        raise ValueError("Canonical archive lacks full transcripts and has no supported fixed-part reconstruction")
    from morse.run_spec import load_run_spec
    from morse.validators import StructuralValidator
    path = HERE / record["config"]
    space = load_run_spec(path).design_space
    source = space.source_construct
    result = []
    for row in rows:
        coding = row["cds"].upper().replace("U", "T")
        if len(coding) != len(source.coding_sequence):
            raise ValueError("Archived benchmark CDS length differs from the frozen construct")
        construct = replace(source, start_codon=coding[:3], cds=coding[3:])
        if not StructuralValidator().validate(construct, space).valid:
            raise ValueError("Archived benchmark CDS fails structural validation in its frozen fixed parts")
        result.append({**row, "sequence": construct.sequence})
    return result, "reconstructed from archived coding sequence and frozen fixed UTR/start/stop/poly(A) construct; unavailable archived fields remain unavailable"


def checkpoint_generation(path):
    from morse.config import OutputSpec
    from morse.optimizer import EvaluationRecord, _evaluation_rows
    payload = pickle.loads(Path(path).read_bytes())
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise ValueError("Unsupported trusted checkpoint format")
    records = payload.get("records", [])
    if not records or not all(isinstance(r, EvaluationRecord) for r in records):
        raise ValueError("Trusted checkpoint has no compatible evaluation records")
    selected = [r for r in records if r.generation == 1]
    names = sorted({name for r in selected for name in r.objective_values})
    specs = [OutputSpec(name, "float", "objective", direction="min") for name in names]
    rows = _evaluation_rows(selected, specs)
    if not rows:
        raise ValueError("Trusted checkpoint contains no generation-one observations")
    return rows, {"path": str(Path(path).resolve()), "trusted_checkpoint_deserialized": True,
                  "generation_one_records": len(rows), "total_checkpoint_records": len(records)}


def _fields(row, aliases):
    if any(name.startswith("output_") for name in row):
        outputs = _historical_output_values(
            {name[7:]: value for name, value in row.items() if name.startswith("output_")}, aliases)
    else:
        # Legacy F3A prefix tables log only these native score observations.
        outputs = {name: row[name] for name in ("mfe_per_nt", "csi_percentile", "csi_raw") if name in row}
    result = {}
    for name, value in outputs.items():
        # Keep scoped outputs distinct from transcript metadata and reject any
        # collision introduced by escaping a conflicting native output name.
        target = f"output_{name}" if name in _TOP_FIELDS or name.startswith(("structural_module_", "threshold_")) else name
        if target in result:
            raise ValueError(f"Ambiguous output observation name for {target}")
        result[target] = value
    for name, value in row.items():
        if name in _TOP_FIELDS or name.startswith("structural_module_"):
            target = name
        elif name.startswith("threshold_"):
            metric, separator, suffix = name[10:].rpartition("_")
            target = f"threshold_{aliases.get(metric, metric)}_{suffix}" if separator else name
        else:
            continue
        if target in result:
            raise ValueError(f"Ambiguous metadata observation mapping for {target}")
        result[target] = value
    return result


def _semantic(value, name):
    if isinstance(value, str) and (name.endswith("_json") or name == "structural_failures"):
        try:
            return json.loads(value)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid structured observation {name}") from exc
    return value


def _same(left, right):
    if isinstance(left, bool) or isinstance(right, bool):
        return type(left) is type(right) and left == right
    if isinstance(left, dict) and isinstance(right, dict):
        return left.keys() == right.keys() and all(_same(left[k], right[k]) for k in left)
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(_same(a, b) for a, b in zip(left, right))
    return left == right


def _tolerance(name, override):
    if name in override:
        a, r = override[name]
        return {"atol": a, "rtol": r, "source": "explicit override"}
    native_name = name.removeprefix("output_")
    if native_name in TOLERANCES:
        a, r = TOLERANCES[native_name]
        return {"atol": a, "rtol": r, "source": "fixed declared metric policy"}
    if _EXACT_COUNTS.fullmatch(name):
        return {"atol": 0.0, "rtol": 0.0, "source": "exact native-count/length/threshold policy"}
    return None


def _safe_value(value):
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {k: _safe_value(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_safe_value(v) for v in value]
    return value


def _number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def compare_rows(archived, fresh, *, expected_size=64, overrides=None, metric_aliases=None):
    aliases = metric_aliases or {}
    if len(set(aliases.values())) != len(aliases):
        raise ValueError("Paper metric aliases must map one-to-one; scorer outputs cannot be merged")
    normalized_overrides = {}
    for name, value in (overrides or {}).items():
        target = aliases.get(name, name)
        if name in aliases and target in _TOP_FIELDS:
            target = f"output_{target}"
        if target in normalized_overrides:
            raise ValueError(f"Ambiguous tolerance override mapping for {target}")
        normalized_overrides[target] = value
    overrides = normalized_overrides
    expected_counts = Counter(row["sequence"] for row in archived)
    observed_counts = Counter(row["sequence"] for row in fresh)
    expected = {row["sequence"]: _fields(row, aliases) for row in archived}
    observed = {row["sequence"]: _fields(row, aliases) for row in fresh}
    paired = sorted(expected.keys() & observed.keys())
    fields = sorted({name for item in [*expected.values(), *observed.values()] for name in item})
    metrics = {}
    for name in fields:
        present = [(sequence, expected[sequence][name], observed[sequence][name]) for sequence in paired
                   if name in expected[sequence] and name in observed[sequence]]
        missing_archive = sum(name not in expected[sequence] for sequence in paired)
        missing_fresh = sum(name not in observed[sequence] for sequence in paired)
        nullable_bound = name.startswith("threshold_") and name.endswith(("_min", "_max", "_equals"))
        boolean = name in {"valid", "structural_valid", "feasibility_pass"} or (name.startswith(("structural_module_", "threshold_")) and name.endswith("_valid"))
        policy = _tolerance(name, overrides)
        numeric = (policy is not None and not nullable_bound and not boolean) or (bool(present) and all(_number(a) and _number(b) for _, a, b in present))
        policy = policy if numeric else None
        exact, accepted, absolute, relative, zero_mismatches, examples = 0, 0, [], [], 0, []
        invalid = 0
        for sequence, old, new in present:
            old, new = _semantic(old, name), _semantic(new, name)
            same = _same(old, new)
            exact += int(same)
            within = same
            if boolean and not (isinstance(old, bool) and isinstance(new, bool)):
                within = False
                invalid += 1
            elif nullable_bound and not name.endswith("_equals") and any(value is not None and (not _number(value) or not math.isfinite(value)) for value in (old, new)):
                within = False
                invalid += 1
            elif numeric:
                if not (_number(old) and _number(new)) or not math.isfinite(old) or not math.isfinite(new):
                    within = False
                    invalid += 1
                else:
                    difference = abs(new - old)
                    absolute.append(difference)
                    if old:
                        relative.append(difference / abs(old))
                    elif difference:
                        zero_mismatches += 1
                    within = policy is not None and difference <= policy["atol"] + policy["rtol"] * abs(old)
            accepted += int(within)
            if not within and len(examples) < 5:
                examples.append({"sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(), "archived": _safe_value(old), "fresh": _safe_value(new)})
        if missing_archive or missing_fresh:
            status = "missing_observations"
        elif invalid:
            status = "invalid_observations"
        elif numeric and policy is None:
            status = "unknown_numeric_tolerance"
        elif accepted != len(present) or not present:
            status = "mismatch"
        else:
            status = "pass"
        known_mismatch = bool(invalid) or (not numeric or policy is not None) and accepted != len(present)
        metrics[name] = {"status": status, "has_known_mismatch": known_mismatch, "paired_observations": len(present), "exact_matches": exact,
            "within_tolerance": accepted, "missing_archived": missing_archive, "missing_fresh": missing_fresh, "invalid_observations": invalid,
            "max_absolute_difference": max(absolute, default=None), "max_relative_difference": max(relative, default=None),
            "zero_reference_mismatches": zero_mismatches, "tolerance": policy if numeric else "exact semantic equality",
            "mismatch_examples": examples}
    population_ok = len(archived) == len(fresh) == expected_size and expected_counts == observed_counts and all(n == 1 for n in expected_counts.values())
    definite_mismatch = not population_ok or any(row["has_known_mismatch"] or row["status"] == "mismatch" for row in metrics.values())
    partial = not metrics or any(row["status"] != "pass" for row in metrics.values())
    status = "fail" if definite_mismatch else "partial" if partial else "pass"
    digest = lambda sequence: hashlib.sha256(sequence.encode()).hexdigest()
    return {"status": status, "population": {"expected_size": expected_size, "archived_records": len(archived), "fresh_records": len(fresh),
        "full_sequence_multiset_matches": expected_counts == observed_counts, "paired_unique_transcripts": len(paired),
        "duplicate_archived_sequences": {digest(s): n for s, n in expected_counts.items() if n > 1},
        "duplicate_fresh_sequences": {digest(s): n for s, n in observed_counts.items() if n > 1},
        "missing_sequence_sha256": [digest(s) for s in expected.keys() - observed.keys()],
        "extra_sequence_sha256": [digest(s) for s in observed.keys() - expected.keys()]},
        "metrics": metrics, "unresolved_fields": [name for name, value in metrics.items() if value["status"] in {"missing_observations", "invalid_observations", "unknown_numeric_tolerance"}],
        "note": "Unknown or missing observations prevent a pass; tolerances are declared before examining differences and are never relaxed automatically."}


def audit(run_root, raw_root, *, run_id=None, overrides=None, allow_checkpoint=False):
    root, raw_root = Path(run_root).expanduser().resolve(), Path(raw_root).expanduser().resolve()
    known = {r["run_id"]: r for r in registry()["runs"]}
    if run_id is not None and (root / run_id).is_dir():
        root = root / run_id
    if root.name == "results":
        root = root.parent
    if run_id is None:
        run_id = read_json(root / "plan.json").get("run_id") if (root / "plan.json").is_file() else root.name
    if run_id not in known or known[run_id]["condition"] == "vaxpress":
        raise ValueError("Select a canonical Core paper run ID (VaxPress has a different initialization protocol)")
    record = known[run_id]
    sources = [name for name in record["sources"] if name.endswith((".jsonl", ".jsonl.gz"))]
    if len(sources) != 1:
        raise ValueError(f"Expected one canonical Core evaluation stream for {run_id}")
    archived, archive_read = first_generation(raw_root / sources[0])
    archived, reconstruction = _reconstruct_archive_sequences(archived, record)
    log = root / "results/run_log.jsonl.gz"
    if not log.is_file():
        log = root / "run_log.jsonl.gz"
    if log.is_file():
        fresh, fresh_read = first_generation(log)
    else:
        checkpoint = root / "results/optimizer-checkpoint.pkl"
        if not checkpoint.is_file():
            checkpoint = root / "optimizer-checkpoint.pkl"
        if not allow_checkpoint or not checkpoint.is_file():
            raise ValueError("No fresh evaluation log. For an older failed run, explicitly allow a trusted checkpoint with --allow-checkpoint")
        fresh, fresh_read = checkpoint_generation(checkpoint)
    aliases = _paper_metric_aliases(record)
    result = compare_rows(archived, fresh, overrides=overrides, metric_aliases=aliases)
    result.update(run_id=run_id, archive_read=archive_read, fresh_read=fresh_read,
        metric_name_mapping=aliases,
        archive_expected_file_sha256=record["sources"][sources[0]], archive_sequence_source=reconstruction,
        claim="Numerical inference replay on the same first-generation inputs; this does not establish a full optimization trajectory replay.")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--raw-root", type=Path, default=HERE / "artifacts/raw")
    parser.add_argument("--run-id")
    parser.add_argument("--output", type=Path, required=True, help="New JSON report path; existing files are never replaced")
    parser.add_argument("--tolerances", type=Path, help="JSON metric -> {atol, rtol}; explicit overrides are recorded")
    parser.add_argument("--allow-checkpoint", action="store_true", help="Permit deserialization of a trusted local checkpoint when no evaluation log exists")
    args = parser.parse_args(argv)
    try:
        if args.output.exists():
            raise ValueError(f"Refusing to replace existing numerical replay report: {args.output}")
        result = audit(args.run_root, args.raw_root, run_id=args.run_id,
                       overrides=tolerances_from_file(args.tolerances), allow_checkpoint=args.allow_checkpoint)
        result["tolerance_override_file"] = str(args.tolerances.resolve()) if args.tolerances else None
        result["audit_script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result["tolerance_override_sha256"] = hashlib.sha256(args.tolerances.read_bytes()).hexdigest() if args.tolerances else None
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        print(json.dumps({"status": result["status"], "run_id": result["run_id"], "paired_transcripts": result["population"]["paired_unique_transcripts"], "unresolved_fields": result["unresolved_fields"], "report": str(args.output.resolve())}, indent=2))
        return {"pass": 0, "fail": 1, "partial": 3}[result["status"]]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"initial score audit: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
