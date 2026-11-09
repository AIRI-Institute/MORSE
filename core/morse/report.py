from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

from morse.config import RunConfig
from morse.run_log import RUN_LOG_FILENAME, read_run_log_tables
from morse.visualization import scope_label


def write_report(
    output_dir: Path, plot_paths: list[Path], run_config: RunConfig, config_path: Path
) -> Path:
    summary = json.loads((output_dir / "summary.json").read_text())
    report_path = output_dir / "report.md"
    part_lengths = summary["input"].get("part_max_lengths")
    present_parts = (
        {part for part, length in part_lengths.items() if length > 0}
        if isinstance(part_lengths, dict) else None
    )
    complete = present_parts is None or len(present_parts) == 6
    module_labels = [
        ("five_prime_utr", "5′UTR", "five_prime_utr_variants"),
        ("start_codon", "start codon", "start_codon_variants"),
        ("three_prime_utr", "3′UTR", "three_prime_utr_variants"),
        ("stop_codon", "stop", "stop_codon_variants"),
        ("poly_a", "poly(A)", "poly_a_variants"),
    ]
    variants = ", ".join(
        f"{label} {summary['input'][key]}" for part, label, key in module_labels
        if present_parts is None or part in present_parts
    )

    lines = [
        "# Full-length mRNA vaccine optimization report" if complete else "# RNA sequence optimization report",
        "",
        f"Run status: **{summary.get('status', 'completed')}**. {summary.get('completion_message', '')}",
        "",
        "## Input sequence",
        "",
        f"- Construct: `{summary['input']['name']}`",
        f"- Run config: `{summary['input']['manifest']}`",
        f"- Reference length: {summary['input']['length_nt']} nt",
        f"- CDS body: {summary['input']['cds_body_codons']} codons; "
        f"encoded protein: {summary['input']['protein_amino_acids']} amino acids",
        f"- Module variants: {variants}",
        "",
        "![Sequence and scope plan](sequence_scope.svg)",
        "",
        "## Configuration",
        "",
        f"- YAML config: `{config_path}`",
        f"- Run name: `{summary['run']['run_name']}`",
        f"- Run index: {summary['run']['run_index']}",
        f"- Resumed from: `{summary['run']['resumed_from'] or 'none'}`",
        "- Structural validator: "
        f"`{summary.get('structural_validation', {}).get('validator', 'not recorded')}`",
        "- Scorers:",
        *_format_scorers(run_config),
        "- Scorer contract: JSON requests in, ordered `outputs` in the JSON response.",
        "  YAML `metrics` assigns each metric a role: optimize, threshold, or log.",
        f"- Scorer invocations: {summary['scorer_timing']['invocations']}",
        f"- Total scorer time: {summary['scorer_timing']['total_seconds']:.3f} s",
        "- Mean scorer time/candidate: "
        f"{summary['scorer_timing']['mean_seconds_per_candidate']:.3f} s",
        f"- Run-chain scorer time: {summary['scorer_timing']['chain_total_seconds']:.3f} s",
        "- Run-chain mean scorer time/candidate: "
        f"{summary['scorer_timing']['chain_mean_seconds_per_candidate']:.3f} s",
        *_format_timing_breakdown(output_dir),
        "",
        "Run chain:",
        "",
    ]

    for run in summary.get("run_chain", {}).get("runs", []):
        lines.append(
            f"- {run['run_index']}. `{run['run_name']}` generations "
            f"{run['generation_start']}..{run['generation_end']} dir `{run['run_dir']}`"
        )

    lines.extend(
        [
        "",
        "## Optimization setup",
        "",
        f"- Algorithm: {summary['optimizer']['algorithm']}",
        f"- Population size: {summary['optimizer']['population_size']}",
        f"- Generations completed: {summary['optimizer']['generations']}",
        "- Generations requested: "
        f"{summary['optimizer'].get('requested_generations', summary['optimizer']['generations'])}",
        f"- Seed: {summary['optimizer']['seed']}",
        *_format_optimizer_parameters(summary["optimizer"]),
        *_format_initial_population(
            summary["optimizer"].get("initial_population", {})
        ),
        *_format_early_stopping(summary["optimizer"].get("early_stopping", {})),
        *_format_hypervolume(summary["optimizer"].get("hypervolume", {})),
        "",
        "Objectives:",
        "",
        ]
    )

    for objective in summary["optimizer"]["objectives"]:
        lines.append(f"- `{objective['output']}` direction `{objective['direction']}`")

    lines.extend(["", "Thresholds:", ""])
    thresholds = summary["output_policy"]["thresholds"]
    if thresholds:
        for threshold in thresholds:
            lines.append(
                f"- `{threshold['output']}` min `{threshold['min']}` max "
                f"`{threshold['max']}` equals `{threshold.get('equals')}`"
            )
    else:
        lines.append("- none")

    structural = summary.get("structural_validation")
    if structural:
        lines.extend(
            [
                "",
                "Structural validation:",
                "",
                f"- Validator: `{structural['validator']}`",
                f"- Passed before scoring: {structural['accepted']}",
                f"- Rejected before scoring: {structural['rejected']}",
                f"- Run-chain rejected before scoring: {structural['chain_rejected']}",
            ]
        )

    lines.extend(
        [
            "",
            "Logged outputs:",
            "",
        ]
    )
    logs = summary["output_policy"]["logs"]
    if logs:
        for item in logs:
            lines.append(f"- `{item['output']}`")
    else:
        lines.append("- none")

    lines.extend(
        [
            "",
            f"- Accepted candidates: {summary['output_policy']['accepted']}",
            f"- Rejected candidates: {summary['output_policy']['rejected']}",
            f"- Run-chain accepted candidates: {summary['output_policy']['chain_accepted']}",
            f"- Run-chain rejected candidates: {summary['output_policy']['chain_rejected']}",
            "",
            "## Result summary",
            "",
            f"- Final-population Pareto candidates: {summary['pareto_count']}",
            f"- Global feasible Pareto archive: {summary.get('pareto_archive_count', summary['pareto_count'])}",
            "",
            "Original sequence:",
            "",
            _format_record(summary["original"], present_parts=present_parts),
            "",
            "Best candidates by objective:",
            "",
        ]
    )

    for output, record in summary["best_by_objective"].items():
        lines.extend([f"### Best by `{output}`", "", _format_record(record, present_parts=present_parts), ""])

    reference_comparison = summary.get("reference_comparison", {})
    if reference_comparison:
        lines.extend(["## Global archive versus reference", ""])
        for output, comparison in reference_comparison.items():
            lines.append(
                f"- `{output}`: reference `{comparison['reference']:.6g}`, "
                f"best evaluated `{comparison['best_evaluated']:.6g}`, "
                f"beats reference `{comparison['beats_reference']}` "
                f"(`{comparison['candidate_id']}`)"
            )
        lines.append("")

    applicability = summary.get("model_applicability", {})
    if applicability.get("warning"):
        lines.extend(
            [
                "## Model applicability warning",
                "",
                f"- {applicability['warning']}",
                "- Non-AUG candidates require explicit experimental justification.",
                "",
            ]
        )

    lines.extend(
        [
            "## Artifacts",
            "",
            "- `run_log.jsonl.gz`: compressed event log with evaluations, "
            "population members, initial-population decisions, and scorer timings",
            "- `run_state.json`: final population and run-chain metadata; NSGA-II "
            "and hard-threshold GA runs can resume",
            "- `resolved_run.yaml`: complete resolved config used for this run",
            "- `sequence_scope.svg`: construct, optimizer, and scorer input scopes",
            "- `pareto.csv`: final non-dominated valid candidates",
            "- `pareto.fasta`: final non-dominated "
            + ("full-length RNA sequences" if complete else "RNA sequences"),
            "- `pareto_constructs.json`: module-resolved Pareto constructs with coordinates",
            "- `pareto_archive.csv`: global hard-feasible non-dominated candidates from all evaluations",
            "- `pareto_archive.fasta`: global hard-feasible Pareto sequences",
            "- `pareto_archive_constructs.json`: module-resolved global archive constructs",
            "- `reference.fasta`: unmodified source construct used as the experimental control",
        ]
    )
    if (output_dir / "provenance.json").is_file():
        lines.append(
            "- `provenance.json`: runtime, configuration, initialization, "
            "and scorer provenance"
        )
    for plot_path in plot_paths:
        lines.append(f"- `{plot_path.name}`")

    lines.extend(["", "## Plots", ""])
    for plot_path in plot_paths:
        lines.extend([f"![{plot_path.stem}]({plot_path.name})", ""])
        vector_path = plot_path.with_suffix(".pdf")
        if vector_path.is_file():
            lines.extend([f"[Vector PDF]({vector_path.name})", ""])
        if plot_path.name == "metric_winners_spider_plots.png":
            lines.extend([
                "Half-spider plot: initial sequence and the existing objective winners "
                "and best-average candidate share the same semicircular axes. Repeated "
                "candidates are drawn once. Each axis retains this run's direction-aware "
                "min–max desirability across the final Pareto set and initial sequence "
                "(0–100; 100 is better), not a CDF or an external "
                "reference percentile. Colours follow scorer order in the resolved run "
                "configuration; all outputs of one scorer share its colour.",
                "",
            ])

    report_path.write_text("\n".join(lines).rstrip() + "\n")
    return report_path


def _format_optimizer_parameters(optimization: dict[str, Any]) -> list[str]:
    parameters = optimization.get("parameters", {})
    crossover = parameters.get("crossover", {})
    mutation = parameters.get("mutation", {})
    lines = [
        "- Offspring size: "
        f"{parameters.get('n_offsprings') or 'population size (pymoo default)'}",
        f"- Eliminate duplicates: {parameters.get('eliminate_duplicates', True)}",
        "- Crossover: "
        f"`{crossover.get('type', 'not recorded')}` probability "
        f"`{crossover.get('probability', 'not recorded')}`",
        "- Mutation: "
        f"`{mutation.get('type', 'not recorded')}` per-variable probability "
        f"`{mutation.get('probability', 'not recorded')}`",
    ]
    reference = parameters.get("reference_directions")
    if reference:
        lines.append(
            "- Reference directions: "
            f"`{reference['method']}` count `{reference['count']}` dimensions "
            f"`{reference['dimensions']}`"
        )
        for key in ("n_partitions", "n_points", "scaling", "seed"):
            if key in reference:
                lines.append(f"  - `{key}`: `{reference[key]}`")
    constraint_handling = parameters.get("constraint_handling")
    if constraint_handling:
        lines.append(
            "- Constraint handling: "
            f"`{constraint_handling.get('method')}`, normalization "
            f"`{constraint_handling.get('normalization')}`"
        )
        if constraint_handling.get("method") == "adaptive_epsilon":
            lines.append(
                "  - threshold epsilon: "
                f"`{constraint_handling.get('initial_epsilon')}` → `0` by generation "
                f"`{constraint_handling.get('relaxation_generations')}` using "
                f"`{constraint_handling.get('curve')}` decay"
            )
        if constraint_handling.get("hard_feasible_archive"):
            lines.append(
                "  - hard-feasible candidates are archived during threshold relaxation"
            )
    return lines


def _format_initial_population(initial_population: dict[str, Any]) -> list[str]:
    if not initial_population:
        return []
    source = initial_population.get("source", "reference")
    if source == "reference":
        return [
            "- Initial population: reference construct plus "
            f"`{initial_population.get('random_fill', 0)}` random candidates"
        ]
    if source == "resume":
        return [
            "- Initial population: resumed population "
            f"(`{initial_population.get('accepted', 0)}` candidates)"
        ]
    lines = [
        "- Initial population file: "
        f"`{initial_population.get('path')}`; provided "
        f"`{initial_population.get('provided', 0)}`, accepted "
        f"`{initial_population.get('accepted', 0)}`, rejected "
        f"`{initial_population.get('rejected', 0)}`, random fill "
        f"`{initial_population.get('random_fill', 0)}`"
    ]
    rejected_ids = initial_population.get("rejected_ids", [])
    relaxed_accepted = int(initial_population.get("relaxed_accepted", 0) or 0)
    if relaxed_accepted:
        lines.append(
            "- Initial candidates admitted by generation-1 threshold epsilon: "
            f"`{relaxed_accepted}`"
        )
    if rejected_ids:
        lines.append(
            "- Rejected initial-population IDs: "
            + ", ".join(f"`{candidate_id}`" for candidate_id in rejected_ids)
        )
    return lines


def _format_scorers(run_config: RunConfig) -> list[str]:
    lines: list[str] = []
    for scorer in run_config.scorers:
        parts = [f"name `{scorer.name}`", f"mode `{scorer.mode}`"]
        if scorer.image:
            parts.append(f"image `{scorer.image}`")
        if scorer.command:
            command = scorer.command if isinstance(scorer.command, str) else ' '.join(scorer.command)
            parts.append(f"command `{command}`")
        parts.append(f"scope `{scope_label(scorer.scope)}`")
        lines.append("  - " + ", ".join(parts))
    return lines


def _format_timing_breakdown(output_dir: Path) -> list[str]:
    path = output_dir / RUN_LOG_FILENAME
    if not path.exists():
        return []
    by_scorer: dict[str, dict[str, float]] = defaultdict(
        lambda: {"seconds": 0.0, "candidates": 0.0, "invocations": 0.0}
    )
    for row in read_run_log_tables(path)["scorer_timing"]:
        scorer = row.get("scorer", "scorer")
        by_scorer[scorer]["seconds"] += float(row.get("elapsed_seconds", 0) or 0)
        by_scorer[scorer]["candidates"] += float(row.get("batch_size", 0) or 0)
        by_scorer[scorer]["invocations"] += 1
    if not by_scorer:
        return []
    lines = ["- Scorer timing by component:"]
    for scorer, stats in sorted(by_scorer.items()):
        mean = stats["seconds"] / stats["candidates"] if stats["candidates"] else 0.0
        lines.append(
            f"  - `{scorer}`: {stats['seconds']:.3f} s total, "
            f"{mean:.3f} s/candidate, {int(stats['invocations'])} invocations"
        )
    return lines


def _format_early_stopping(early_stopping: dict[str, Any]) -> list[str]:
    if not early_stopping or not early_stopping.get("enabled"):
        return ["- Early stopping: disabled"]
    reason = early_stopping.get("stop_reason") or "not triggered"
    return [
        "- Early stopping: enabled "
        f"metric `{early_stopping.get('metric')}` "
        f"min_generations `{early_stopping.get('min_generations')}` "
        f"patience `{early_stopping.get('patience')}` "
        f"rel_tol `{early_stopping.get('rel_tol')}` "
        f"abs_tol `{early_stopping.get('abs_tol', 0.0)}`",
        f"- Stop reason: {reason}",
    ]


def _format_hypervolume(hypervolume: dict[str, Any]) -> list[str]:
    if not hypervolume or not hypervolume.get("available"):
        return [
            "- Hypervolume: unavailable (requires at least two objectives and "
            "a selection-feasible population)"
        ]
    return [
        "- Hypervolume: normalized against the first selection-feasible "
        "population; reference point "
        f"`{hypervolume.get('reference_point', 'not recorded')}`",
        f"- Final hypervolume: `{float(hypervolume['final']):.6g}`",
        f"- Best hypervolume: `{float(hypervolume['best']):.6g}`",
    ]


def _format_record(record: dict[str, Any], *, present_parts: set[str] | None = None) -> str:
    lines = [
        f"- Candidate: `{record['candidate_id']}`",
        f"- Valid: {record['valid']}",
        f"- Length: {record['length_nt']} nt",
    ]
    for part, detail in (
        ("five_prime_utr", f"- 5′UTR: `{record['modules']['five_prime_utr']}`"),
        ("three_prime_utr", f"- 3′UTR: `{record['modules']['three_prime_utr']}`"),
        ("stop_codon", f"- Stop codon: `{record['modules']['stop_codon']}`"),
        ("poly_a", f"- poly(A): {record['modules']['poly_a_length']} nt"),
    ):
        if present_parts is None or part in present_parts:
            lines.append(detail)
    lines.extend([
        f"- Changed codons: {record['changed_codons']}",
        f"- GC fraction: {record['gc_fraction']:.4f}",
    ])
    for key, value in record["outputs"].items():
        if key.rsplit(".", 1)[-1] == "mfe_structure":
            continue
        if isinstance(value, float):
            lines.append(f"- Output `{key}`: {value:.6g}")
        else:
            lines.append(f"- Output `{key}`: {value}")
    failed = [threshold for threshold in record["thresholds"] if not threshold["valid"]]
    for threshold in failed:
        lines.append(f"- Threshold failure `{threshold['output']}`: {threshold['reason']}")
    structural = record.get("structural_validation", {})
    for failure in structural.get("failures", []):
        lines.append(
            f"- Structural failure `{failure['module']}.{failure['code']}`: "
            f"{failure['message']}"
        )
    return "\n".join(lines)
