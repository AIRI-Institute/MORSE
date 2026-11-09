from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np

from morse.optimizer import VaccineOptimizationProblem
from morse.run_spec import RunSpec
from morse.validators import validator_for_policy


def run_pilot(
    spec: RunSpec,
    scorer: Any,
    *,
    sample_count: int,
    seed: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Score the reference and a reproducible random design-space sample."""

    requested = int(sample_count)
    if not 2 <= requested <= 512:
        raise ValueError("pilot sample_count must be between 2 and 512")
    output_dir.mkdir(parents=True, exist_ok=True)
    design = spec.design_space
    design_space_size = math.prod(
        high - low + 1
        for low, high in zip(
            design.lower_bounds,
            design.upper_bounds,
            strict=True,
        )
    )
    actual_count = min(requested, design_space_size)
    if actual_count < 2:
        raise ValueError("The design space needs at least two candidates for a pilot")
    progress_path = output_dir / "progress.json"
    progress_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "kind": "pilot",
                "stage": "sampling",
                "evaluated_candidates": 0,
                "total_candidates": actual_count,
                "progress_fraction": 0.0,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )

    vectors = _sample_vectors(
        design.lower_bounds,
        design.upper_bounds,
        design.original_vector,
        actual_count,
        seed,
    )
    constructs = [design.decode(vector) for vector in vectors]
    problem = VaccineOptimizationProblem(
        design,
        scorer,
        list(spec.run_config.outputs),
        validator=validator_for_policy(spec.run_config.manufacturability),
        constraint_handling=spec.run_config.optimizer.constraint_handling,
    )
    started = time.perf_counter()
    problem.score_constructs(constructs)
    elapsed = time.perf_counter() - started
    records = [problem.cache[construct.fingerprint()] for construct in constructs]
    valid = [record for record in records if record.valid]
    structurally_valid = [
        record for record in records if record.structural_validation.valid
    ]
    objectives = [
        output for output in spec.run_config.outputs if output.role == "objective"
    ]
    thresholds = [
        output for output in spec.run_config.outputs if output.role == "threshold"
    ]
    metric_rows = []
    for output in spec.run_config.outputs:
        numeric = [
            float(record.outputs[output.name])
            for record in structurally_valid
            if isinstance(record.outputs.get(output.name), (int, float))
            and not isinstance(record.outputs.get(output.name), bool)
            and math.isfinite(float(record.outputs[output.name]))
        ]
        metric_rows.append(
            {
                "name": output.name,
                "role": output.role,
                "direction": output.direction,
                "count": len(numeric),
                "min": min(numeric) if numeric else None,
                "median": float(np.median(numeric)) if numeric else None,
                "max": max(numeric) if numeric else None,
                "reference": records[0].outputs.get(output.name),
            }
        )

    threshold_rows = []
    for output in thresholds:
        checked_records = [
            (
                record,
                next(
                    check
                    for check in record.threshold_results
                    if check.output == output.name
                ),
            )
            for record in records
            if record.structural_validation.valid
        ]
        checks = [check for _record, check in checked_records]
        reference_check = next(
            (
                check
                for check in records[0].threshold_results
                if check.output == output.name
            ),
            None,
        )
        threshold_rows.append(
            {
                "name": output.name,
                "min": output.min_value,
                "max": output.max_value,
                "equals": output.equals_value,
                "pass_fraction": (
                    sum(check.valid for check in checks) / len(checks) if checks else 0.0
                ),
                "reference_pass": (
                    reference_check.valid if reference_check is not None else None
                ),
                "failing_examples": [
                    {
                        "candidate_id": record.candidate_id,
                        "start_codon": record.construct.start_codon,
                        "changed_codons": record.changed_codons,
                        "observed": check.observed,
                        "reason": check.reason,
                    }
                    for record, check in checked_records
                    if not check.valid
                ][:3],
            }
        )

    correlations = []
    for left_index, left in enumerate(objectives):
        for right in objectives[left_index + 1 :]:
            pairs = [
                (record.outputs.get(left.name), record.outputs.get(right.name))
                for record in structurally_valid
            ]
            numeric_pairs = [
                (float(a), float(b))
                for a, b in pairs
                if isinstance(a, (int, float))
                and not isinstance(a, bool)
                and isinstance(b, (int, float))
                and not isinstance(b, bool)
                and math.isfinite(float(a))
                and math.isfinite(float(b))
            ]
            correlation = None
            if len(numeric_pairs) >= 3:
                left_values, right_values = zip(*numeric_pairs, strict=True)
                if np.std(left_values) > 0 and np.std(right_values) > 0:
                    correlation = float(np.corrcoef(left_values, right_values)[0, 1])
            correlations.append(
                {
                    "left": left.name,
                    "right": right.name,
                    "pearson": correlation,
                }
            )

    optimizer = spec.run_config.optimizer
    expected_evaluations = optimizer.population_size + (
        optimizer.n_offsprings or optimizer.population_size
    ) * max(optimizer.generations - 1, 0)
    mean_seconds = elapsed / len(records)
    result = {
        "schema_version": 1,
        "kind": "pilot",
        "sample_count": len(records),
        "seed": int(seed),
        "elapsed_seconds": elapsed,
        "mean_wall_seconds_per_candidate": mean_seconds,
        "structurally_valid_fraction": len(structurally_valid) / len(records),
        "hard_feasible_fraction": len(valid) / len(records),
        "expected_optimizer_evaluations": expected_evaluations,
        "estimated_optimizer_seconds": mean_seconds * expected_evaluations,
        "reference": {
            "candidate_id": records[0].candidate_id,
            "structural_valid": records[0].structural_validation.valid,
            "hard_feasible": records[0].valid,
            "outputs": records[0].outputs,
        },
        "structural_failure_examples": [
            {
                "candidate_id": record.candidate_id,
                "start_codon": record.construct.start_codon,
                "changed_codons": record.changed_codons,
                "failures": [
                    failure.to_dict()
                    for failure in record.structural_validation.failures
                ],
            }
            for record in records
            if not record.structural_validation.valid
        ][:3],
        "metrics": metric_rows,
        "thresholds": threshold_rows,
        "objective_correlations": correlations,
        "warnings": _pilot_warnings(metric_rows, threshold_rows, correlations),
        "blocking_issues": _pilot_blocking_issues(metric_rows, threshold_rows),
    }
    (output_dir / "pilot.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    progress_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "kind": "pilot",
                "stage": "complete",
                "evaluated_candidates": len(records),
                "total_candidates": len(records),
                "progress_fraction": 1.0,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return result


def _sample_vectors(
    lower: list[int],
    upper: list[int],
    original: list[int],
    count: int,
    seed: int,
) -> list[list[int]]:
    rng = np.random.default_rng(seed)
    vectors = [list(original)]
    seen = {tuple(original)}
    attempts = 0
    while len(vectors) < count:
        attempts += 1
        vector = [
            int(rng.integers(low, high + 1)) if high > low else int(low)
            for low, high in zip(lower, upper, strict=True)
        ]
        key = tuple(vector)
        if key not in seen:
            seen.add(key)
            vectors.append(vector)
        if attempts > count * 10_000:
            raise ValueError(
                "The design space is too small for the requested unique pilot sample"
            )
    return vectors


def _pilot_warnings(
    metrics: list[dict[str, Any]],
    thresholds: list[dict[str, Any]],
    correlations: list[dict[str, Any]],
) -> list[str]:
    warnings: list[str] = []
    for metric in metrics:
        if (
            metric["role"] == "objective"
            and metric["count"] > 0
            and metric["min"] == metric["max"]
        ):
            warnings.append(f"Objective {metric['name']} is constant in the pilot sample")
    for threshold in thresholds:
        if threshold["pass_fraction"] == 0:
            warnings.append(f"No pilot candidate passes threshold {threshold['name']}")
    for row in correlations:
        correlation = row["pearson"]
        if correlation is not None and abs(correlation) >= 0.95:
            warnings.append(
                f"Objectives {row['left']} and {row['right']} are strongly correlated "
                f"in the pilot sample (r={correlation:.3f})"
            )
    return warnings


def _pilot_blocking_issues(
    metrics: list[dict[str, Any]], thresholds: list[dict[str, Any]]
) -> list[dict[str, str]]:
    issues: list[dict[str, str]] = []
    for metric in metrics:
        if (
            metric["role"] == "objective"
            and metric["count"] > 0
            and metric["min"] == metric["max"]
        ):
            issues.append(
                {
                    "type": "constant objective",
                    "metric": str(metric["name"]),
                    "action": "Change its scope or role before the full run",
                }
            )
    for threshold in thresholds:
        if threshold["pass_fraction"] == 0:
            issues.append(
                {
                    "type": "infeasible threshold",
                    "metric": str(threshold["name"]),
                    "action": "Adjust the threshold before the full run",
                }
            )
    return issues
