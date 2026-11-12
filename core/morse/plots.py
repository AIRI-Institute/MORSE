from __future__ import annotations

import csv
import json
import os
import sys
import tempfile
import textwrap
from collections import Counter, defaultdict
from contextvars import ContextVar
from functools import wraps
from pathlib import Path
from statistics import median
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "morse-matplotlib"))

import matplotlib

_ipython = sys.modules.get("IPython")
_shell = getattr(_ipython, "get_ipython", lambda: None)()
if "matplotlib.pyplot" not in sys.modules and getattr(_shell, "kernel", None) is None:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from morse.publication_style import (
    COLORS, FIXED_FILL, GRID, INITIAL, INK, MUTED, RC_STYLE,
    palette_color, run_color_mapping,
)
from morse.run_log import RUN_LOG_FILENAME, read_run_log_tables


_RUN_COLORS: ContextVar[dict[str, dict[str, str]] | None] = ContextVar("morse_plot_colors", default=None)


def _read_run_colors(output_dir: Path, summary: dict[str, Any] | None = None):
    import yaml
    scorers = None
    for filename in ("resolved_run.yaml", "experiment.yaml"):
        path = output_dir / filename
        if path.is_file():
            config = yaml.safe_load(path.read_text())
            if isinstance(config, dict) and isinstance(config.get("scorers"), list):
                scorers = config["scorers"]
                if config.get("schema_version") == 4:
                    scorers = [{**item, "namespace": item["name"]} for item in scorers]
                break
    if scorers is None:
        scorers = (summary or {}).get("scorers")
    if not isinstance(scorers, list):
        # Standalone diagnostic plots may have only a summary, with no saved
        # scorer configuration. Preserve its declared objective order then.
        scorers = [{"name": item["output"], "outputs": {item["output"]: {}}}
                   for item in (summary or {}).get("optimizer", {}).get("objectives", [])]
    return run_color_mapping(scorers)


def _publication_plot(function):
    @wraps(function)
    def styled(*args, **kwargs):
        token = None
        if _RUN_COLORS.get() is None:
            summary = kwargs.get("summary")
            if summary is None:
                summary = next((item for item in args[1:] if isinstance(item, dict)), None)
            token = _RUN_COLORS.set(_read_run_colors(Path(args[0]), summary))
        try:
            with plt.rc_context(RC_STYLE):
                return function(*args, **kwargs)
        finally:
            if token is not None:
                _RUN_COLORS.reset(token)
    return styled


def _output_color(name: str, index: int = 0) -> str:
    key = name.removeprefix("output_")
    return (_RUN_COLORS.get() or {}).get("outputs", {}).get(key, palette_color(index))


def _scorer_color(name: str, index: int = 0) -> str:
    return (_RUN_COLORS.get() or {}).get("scorers", {}).get(name, palette_color(index))


def _save_plot(fig, path: Path) -> None:
    """Keep PNG report contracts and export a vector sibling without style leaks."""
    for axis in fig.axes:
        if not axis.axison:
            continue
        axis.set_axisbelow(True)
        axis.grid(False)
        axis.grid(axis="y", color=GRID, linewidth=0.5, alpha=0.5)
        for name, spine in axis.spines.items():
            if name in {"top", "right"}:
                spine.set_visible(False)
            spine.set_linewidth(0.6)
            spine.set_color(MUTED)
    if fig._suptitle is not None:
        fig._suptitle.set(x=0.065, ha="left", fontsize=11, fontweight="bold")
    scorer_count = len((_RUN_COLORS.get() or {}).get("scorers", {}))
    if scorer_count > len(COLORS):
        fig.text(0.99, 0.005, f"{scorer_count} scorer bindings · extended palette",
                 ha="right", color=MUTED, fontsize=6.5)
    fig.savefig(path, dpi=220, bbox_inches="tight", pad_inches=0.12, facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.12,
                facecolor="white", metadata={"CreationDate": None, "ModDate": None})


@_publication_plot
def plot_run(output_dir: Path) -> list[Path]:
    summary = json.loads((output_dir / "summary.json").read_text())
    run_log = read_run_log_tables(output_dir / RUN_LOG_FILENAME)
    history_rows = run_log["population_member"]
    evaluations = run_log["evaluation"]
    pareto_rows = _read_csv(output_dir / "pareto.csv")
    scorer_timings = run_log["scorer_timing"]

    plot_paths = [
        _plot_hypervolume_history(output_dir, summary),
        _plot_constraint_flow(output_dir, evaluations, history_rows, summary),
        _plot_threshold_metric_history(output_dir, evaluations, summary),
        _plot_population_history(output_dir, history_rows, summary, metric_group="objective"),
        _plot_population_history(output_dir, history_rows, summary, metric_group="measured"),
        _plot_single_objective_histogram(output_dir, evaluations, summary),
        _plot_two_objective_pareto(output_dir, summary, pareto_rows),
        _plot_metric_winners_radar(output_dir, summary, pareto_rows),
        _plot_module_frequencies(output_dir, history_rows, summary),
        _plot_codon_changes(output_dir, summary, evaluations),
        _plot_scorer_timings(output_dir, scorer_timings),
        _plot_scorer_timing_breakdown(output_dir, scorer_timings),
    ]
    return [path for path in plot_paths if path is not None]


@_publication_plot
def _plot_hypervolume_history(
    output_dir: Path,
    summary: dict[str, Any],
) -> Path | None:
    hypervolume = summary.get("optimizer", {}).get("hypervolume", {})
    history = hypervolume.get("history", []) if isinstance(hypervolume, dict) else []
    rows = [
        item
        for item in history
        if isinstance(item, dict)
        and _is_float(item.get("generation"))
        and _is_float(item.get("hypervolume"))
    ]
    if not rows:
        return None

    generations = [int(float(item["generation"])) for item in rows]
    current = [_float(item["hypervolume"]) for item in rows]
    best = [
        _float(item.get("best_hypervolume", item["hypervolume"]))
        for item in rows
    ]
    reference_point = hypervolume.get("reference_point", "not recorded")

    fig, axis = plt.subplots(figsize=(9, 4.2))
    axis.plot(
        generations,
        current,
        marker="o",
        markersize=3,
        markevery=max(1, len(generations) // 10),
        linewidth=1.4,
        color=INK,
        label="current",
    )
    axis.plot(
        generations,
        best,
        linestyle="--",
        linewidth=1.1,
        color=MUTED,
        label="best so far",
    )
    axis.set_title("Normalized hypervolume across generations")
    axis.set_xlabel("Generation")
    axis.set_ylabel("Hypervolume")
    axis.grid(alpha=0.25)
    axis.legend(frameon=False)
    axis.text(
        0.01,
        0.01,
        "Selection-feasible Pareto set · "
        f"reference point {reference_point}",
        transform=axis.transAxes,
        fontsize=8,
        color=MUTED,
    )
    _add_run_boundaries(axis, summary)
    fig.tight_layout()

    path = output_dir / "hypervolume_history.png"
    _save_plot(fig, path)
    plt.close(fig)
    return path


@_publication_plot
def _plot_population_history(
    output_dir: Path,
    rows: list[dict[str, Any]],
    summary: dict[str, Any] | None = None,
    metric_group: str = "all",
) -> Path | None:
    if not rows:
        return None

    metrics = _population_metrics(rows, summary, metric_group=metric_group)
    available = [(key, label) for key, label in metrics if _has_numeric(rows, key)]
    if not available:
        return None

    fig, axes = plt.subplots(len(available), 1, figsize=(9, 2.6 * len(available)), sharex=True)
    if len(available) == 1:
        axes = [axes]

    generation_key = "chain_generation" if "chain_generation" in rows[0] else "generation"
    generations = sorted({int(float(row[generation_key])) for row in rows})
    grouped: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[int(float(row[generation_key]))].append(row)

    for index, (axis, (metric, label)) in enumerate(zip(axes, available, strict=True)):
        mins: list[float] = []
        meds: list[float] = []
        maxs: list[float] = []
        for generation in generations:
            values = [
                _float(row[metric])
                for row in grouped[generation]
                if _is_float(row.get(metric))
            ]
            mins.append(min(values))
            meds.append(float(median(values)))
            maxs.append(max(values))
        color = _output_color(metric, index)
        axis.fill_between(generations, mins, maxs, color=color, alpha=0.12, label="min-max")
        axis.plot(generations, meds, color=color, linewidth=1.5, label="median")
        axis.plot(generations, mins, color=color, linewidth=0.8, linestyle="--", alpha=0.8, label="min")
        axis.set_ylabel(label)
        axis.grid(True, alpha=0.25)
        if summary is not None:
            _add_run_boundaries(axis, summary)
            if metric_group == "objective":
                _annotate_best_objective_delta(axis, metric, summary, generations[-1])
    axes[-1].set_xlabel("Generation")
    axes[0].legend(loc="best", fontsize=8)
    title, filename = _population_plot_metadata(metric_group)
    fig.suptitle(title)
    fig.tight_layout(rect=(0, 0, 1, 0.98))

    path = output_dir / filename
    _save_plot(fig, path)
    plt.close(fig)
    return path


@_publication_plot
def _plot_constraint_flow(
    output_dir: Path,
    evaluations: list[dict[str, Any]],
    population_rows: list[dict[str, Any]],
    summary: dict[str, Any],
) -> Path | None:
    grouped_evaluations = _group_by_generation(evaluations)
    if not grouped_evaluations:
        return None

    grouped_population = _group_by_generation(population_rows)
    generations = sorted(set(grouped_evaluations) | set(grouped_population))
    feasible: list[int] = []
    threshold_rejected: list[int] = []
    structural_rejected: list[int] = []
    population_invalid: list[int] = []
    threshold_epsilons: list[float] = []

    for generation in generations:
        rows = grouped_evaluations.get(generation, [])
        structural_failures = sum(
            not _parse_bool(row.get("structural_valid", True)) for row in rows
        )
        valid = sum(_parse_bool(row.get("valid", False)) for row in rows)
        feasible.append(valid)
        structural_rejected.append(structural_failures)
        threshold_rejected.append(len(rows) - valid - structural_failures)
        population_invalid.append(
            sum(
                not _parse_bool(row.get("valid", False))
                for row in grouped_population.get(generation, [])
            )
        )
        epsilon_values = [
            _float(row["threshold_epsilon"])
            for row in grouped_population.get(generation, [])
            if _is_float(row.get("threshold_epsilon"))
        ]
        threshold_epsilons.append(
            float(np.median(epsilon_values)) if epsilon_values else 0.0
        )

    x = np.asarray(generations, dtype=float)
    feasible_array = np.asarray(feasible, dtype=float)
    threshold_array = np.asarray(threshold_rejected, dtype=float)
    structural_array = np.asarray(structural_rejected, dtype=float)
    population_invalid_array = np.asarray(population_invalid, dtype=float)
    width = min(0.8, max(0.18, 8.0 / max(len(generations), 1)))

    fig, axis = plt.subplots(figsize=(10, 4.8))
    axis.bar(
        x,
        feasible_array,
        width=width,
        color=COLORS[1],
        alpha=0.80,
        label="passed all constraints",
    )
    axis.bar(
        x,
        threshold_array,
        width=width,
        bottom=feasible_array,
        color=COLORS[0],
        alpha=0.65,
        label="threshold rejected",
    )
    axis.bar(
        x,
        structural_array,
        width=width,
        bottom=feasible_array + threshold_array,
        color=INITIAL,
        label="structurally rejected",
    )
    axis.plot(
        x,
        population_invalid_array,
        color=INK,
        marker="D",
        markersize=3,
        markevery=max(1, len(generations) // 10),
        linewidth=1.0,
        markerfacecolor="white",
        label=(
            "invalid after "
            f"{summary.get('optimizer', {}).get('algorithm', 'optimizer')} selection"
        ),
        zorder=5,
    )

    rejected = threshold_array + structural_array
    totals = feasible_array + rejected
    annotation_generations = set(_readable_generation_ticks(generations))
    for generation, total, count in zip(x, totals, rejected, strict=True):
        if count <= 0 or generation not in annotation_generations:
            continue
        axis.annotate(
            f"−{int(count)}",
            (generation, total),
            xytext=(0, -5),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=8,
            color=MUTED,
        )

    axis.set_xlabel("Generation")
    axis.set_ylabel("Candidates")
    axis.set_title("Candidate constraint flow by generation")
    axis.set_xticks(_readable_generation_ticks(generations))
    axis.set_ylim(bottom=0)
    axis.grid(axis="y", alpha=0.25)
    secondary_axis = None
    if any(epsilon > 0.0 for epsilon in threshold_epsilons):
        secondary_axis = axis.twinx()
        secondary_axis.plot(
            x,
            np.asarray(threshold_epsilons, dtype=float),
            color=COLORS[1],
            linewidth=1.8,
            linestyle="--",
            label="selection threshold ε",
        )
        secondary_axis.set_ylabel("Normalized threshold ε", color=COLORS[1])
        secondary_axis.tick_params(axis="y", labelcolor=COLORS[1])
        secondary_axis.set_ylim(bottom=0)
    handles, labels = axis.get_legend_handles_labels()
    if secondary_axis is not None:
        extra_handles, extra_labels = secondary_axis.get_legend_handles_labels()
        handles.extend(extra_handles)
        labels.extend(extra_labels)
    axis.legend(handles, labels, loc="best", fontsize=8)
    _add_run_boundaries(axis, summary)
    fig.tight_layout()

    path = output_dir / "constraint_flow_history.png"
    _save_plot(fig, path)
    plt.close(fig)
    return path


@_publication_plot
def _plot_threshold_metric_history(
    output_dir: Path,
    evaluations: list[dict[str, Any]],
    summary: dict[str, Any],
) -> Path | None:
    grouped = _group_by_generation(evaluations)
    if not grouped:
        return None

    thresholds = [
        item
        for item in summary.get("output_policy", {}).get("thresholds", [])
        if _has_numeric(evaluations, f"threshold_{item['output']}_observed")
    ]
    if not thresholds:
        return None

    ncols = 2 if len(thresholds) > 1 else 1
    nrows = (len(thresholds) + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(11, 3.8 * nrows),
        squeeze=False,
        sharex=False,
    )

    for axis, threshold in zip(axes.flat, thresholds, strict=False):
        output = str(threshold["output"])
        value_key = f"threshold_{output}_observed"
        valid_key = f"threshold_{output}_valid"
        generations: list[int] = []
        mins: list[float] = []
        p10s: list[float] = []
        medians: list[float] = []
        p90s: list[float] = []
        maxs: list[float] = []
        failed: list[int] = []
        failure_points: list[float | None] = []
        totals: list[int] = []

        for generation in sorted(grouped):
            rows = grouped[generation]
            values = [
                _float(row[value_key])
                for row in rows
                if _is_float(row.get(value_key))
            ]
            if not values:
                continue
            generations.append(generation)
            mins.append(float(np.min(values)))
            p10s.append(float(np.percentile(values, 10)))
            medians.append(float(np.median(values)))
            p90s.append(float(np.percentile(values, 90)))
            maxs.append(float(np.max(values)))
            failed_values = [
                _float(row[value_key])
                for row in rows
                if _is_float(row.get(value_key))
                and not _parse_bool(row.get(valid_key, True))
            ]
            failed.append(len(failed_values))
            failure_points.append(
                _violation_marker_value(failed_values, threshold)
                if failed_values
                else None
            )
            totals.append(len(values))

        if not generations:
            axis.set_visible(False)
            continue

        bounds = _threshold_numeric_bounds(threshold)
        scale_values = [*mins, *maxs, *[value for value in bounds if value is not None]]
        y_min, y_max = _padded_limits(scale_values)
        axis.set_ylim(y_min, y_max)
        _shade_allowed_threshold_region(axis, threshold, y_min, y_max)

        axis.vlines(
            generations,
            mins,
            maxs,
            color=MUTED,
            alpha=0.55,
            linewidth=0.6,
            label="min–max",
        )
        axis.fill_between(
            generations,
            p10s,
            p90s,
            color=_output_color(output),
            alpha=0.13,
            label="10–90%",
        )
        axis.plot(
            generations,
            medians,
            color=_output_color(output),
            marker="o",
            markersize=2.5,
            markevery=max(1, len(generations) // 10),
            linewidth=1.4,
            label="median",
        )
        _draw_threshold_lines(axis, threshold)

        annotation_generations = set(
            generations
            if len(generations) <= 12
            else [
                generations[int(index)]
                for index in np.linspace(
                    0, len(generations) - 1, num=6, dtype=int
                )
            ]
        )
        for generation, marker_value, failed_count, total in zip(
            generations, failure_points, failed, totals, strict=True
        ):
            if failed_count <= 0 or marker_value is None:
                continue
            axis.scatter(
                [generation],
                [marker_value],
                color=INK,
                edgecolor="white",
                linewidth=0.6,
                s=30,
                zorder=6,
            )
            if generation in annotation_generations:
                axis.annotate(
                    f"{failed_count}/{total}",
                    (generation, marker_value),
                    xytext=(0, 6),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    color=MUTED,
                )

        axis.set_title(_metric_label(output))
        axis.set_xlabel("Generation")
        axis.set_ylabel(_threshold_axis_label(output))
        axis.set_xticks(_readable_generation_ticks(generations))
        axis.grid(alpha=0.25)
        _add_run_boundaries(axis, summary)

    for axis in axes.flat[len(thresholds) :]:
        axis.set_visible(False)

    handles, labels = axes.flat[0].get_legend_handles_labels()
    if handles:
        unique_legend = dict(zip(labels, handles, strict=True))
        fig.legend(
            unique_legend.values(),
            unique_legend.keys(),
            loc="upper center",
            bbox_to_anchor=(0.5, 0.965),
            ncol=min(4, len(handles)),
            fontsize=8,
        )
    fig.suptitle("Threshold metric history across evaluated candidates", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.90))

    path = output_dir / "threshold_metric_history.png"
    _save_plot(fig, path)
    plt.close(fig)
    return path


@_publication_plot
def _plot_single_objective_histogram(
    output_dir: Path,
    rows: list[dict[str, Any]],
    summary: dict[str, Any],
) -> Path | None:
    objectives = summary.get("optimizer", {}).get("objectives", [])
    if len(objectives) != 1:
        return None

    objective = objectives[0]
    output = str(objective["output"])
    output_key = f"output_{output}"
    values = [
        _float(row[output_key])
        for row in rows
        if _is_finite_float(row.get(output_key))
    ]
    if not values:
        return None

    direction = str(objective["direction"])
    bin_count = min(30, max(5, round(np.sqrt(len(values)))))
    fig, axis = plt.subplots(figsize=(8.5, 5.2))
    axis.hist(
        values,
        bins=bin_count,
        color=_output_color(output),
        edgecolor="white",
        linewidth=0.6,
        alpha=0.72,
    )

    original_value = summary.get("original", {}).get("outputs", {}).get(output)
    if _is_finite_float(original_value):
        axis.axvline(
            _float(original_value),
            color=INITIAL,
            linewidth=1.2,
            linestyle="--",
            label="original",
        )

    best_value = (
        summary.get("best_by_objective", {})
        .get(output, {})
        .get("outputs", {})
        .get(output)
    )
    if _is_finite_float(best_value):
        axis.axvline(
            _float(best_value),
            color=_output_color(output),
            linewidth=1.5,
            label="best feasible",
        )

    axis.set_xlabel(_radar_metric_label(output, direction))
    axis.set_ylabel("Evaluated candidates")
    axis.set_title("Single-objective value distribution")
    axis.grid(axis="y", alpha=0.25)
    handles, _ = axis.get_legend_handles_labels()
    if handles:
        axis.legend(loc="best")
    fig.tight_layout()

    path = output_dir / "objective_histogram.png"
    _save_plot(fig, path)
    plt.close(fig)
    return path


@_publication_plot
def _plot_two_objective_pareto(
    output_dir: Path,
    summary: dict[str, Any],
    rows: list[dict[str, Any]],
) -> Path | None:
    objectives = summary.get("optimizer", {}).get("objectives", [])
    if len(objectives) != 2 or not rows:
        return None

    x_objective, y_objective = objectives
    x_metric = str(x_objective["output"])
    y_metric = str(y_objective["output"])
    x_key = f"output_{x_metric}"
    y_key = f"output_{y_metric}"
    if not _has_numeric(rows, x_key) or not _has_numeric(rows, y_key):
        return None

    x = np.asarray([_float(row[x_key]) for row in rows])
    y = np.asarray([_float(row[y_key]) for row in rows])
    fig, axis = plt.subplots(figsize=(8, 6))
    axis.scatter(
        x,
        y,
        s=22,
        color=INK,
        edgecolor="white",
        linewidth=0.5,
        alpha=0.72,
        label="Pareto candidates",
    )

    original_outputs = summary.get("original", {}).get("outputs", {})
    if _is_finite_float(original_outputs.get(x_metric)) and _is_finite_float(
        original_outputs.get(y_metric)
    ):
        axis.scatter(
            [_float(original_outputs[x_metric])],
            [_float(original_outputs[y_metric])],
            marker="x",
            s=42,
            color=INITIAL,
            linewidth=1.4,
            label="original",
            zorder=5,
        )

    annotations: dict[tuple[float, float], list[str]] = defaultdict(list)
    for metric, record in summary.get("best_by_objective", {}).items():
        outputs = record.get("outputs", {})
        if _is_finite_float(outputs.get(x_metric)) and _is_finite_float(
            outputs.get(y_metric)
        ):
            point = (_float(outputs[x_metric]), _float(outputs[y_metric]))
            annotations[point].append(_radar_metric_name(metric))

    x_midpoint = (float(np.min(x)) + float(np.max(x))) / 2
    for (point_x, point_y), metrics in annotations.items():
        annotate_left = point_x > x_midpoint
        axis.annotate(
            "best " + "\n+ ".join(metrics),
            (point_x, point_y),
            textcoords="offset points",
            xytext=(-8 if annotate_left else 8, 10),
            ha="right" if annotate_left else "left",
            va="bottom",
            fontsize=8,
            annotation_clip=False,
            arrowprops={"arrowstyle": "-", "color": MUTED, "linewidth": 0.6},
        )
    for index, (metric, record) in enumerate(summary.get("best_by_objective", {}).items()):
        outputs = record.get("outputs", {})
        if _is_finite_float(outputs.get(x_metric)) and _is_finite_float(outputs.get(y_metric)):
            axis.scatter([_float(outputs[x_metric])], [_float(outputs[y_metric])],
                         color=_output_color(metric, index), s=36, edgecolor="white",
                         linewidth=0.7, zorder=6, label="best " + _radar_metric_name(metric))

    axis.set_xlabel(
        _radar_metric_label(x_metric, str(x_objective["direction"]))
    )
    axis.set_ylabel(
        _radar_metric_label(y_metric, str(y_objective["direction"]))
    )
    axis.set_title("Final two-objective Pareto front")
    axis.grid(alpha=0.25)
    axis.margins(x=0.12, y=0.1)
    axis.legend(loc="best")
    fig.tight_layout()

    path = output_dir / "pareto_front.png"
    _save_plot(fig, path)
    plt.close(fig)
    return path


@_publication_plot
def _plot_metric_winners_radar(
    output_dir: Path, summary: dict[str, Any], rows: list[dict[str, Any]]
) -> Path | None:
    chart_data = _metric_winner_radar_data(summary, rows)
    if chart_data is None:
        # One or two objectives are already represented by the histogram/Pareto
        # charts; a semicircle with only its two end spokes has no visible area.
        return None
    objectives, panels, original_scores = chart_data
    angles = np.linspace(np.pi, 0.0, len(objectives))
    initial_id = str(summary.get("original", {}).get("candidate_id", "original"))
    curves = []
    by_id = {}
    initial_roles = []
    for index, panel in enumerate(panels):
        if panel["is_original"]:
            continue
        if panel["candidate_id"] == initial_id:
            initial_roles.append(panel["title"])
            continue
        identifier = panel["candidate_id"]
        if identifier in by_id:
            by_id[identifier]["titles"].append(panel["title"])
            continue
        color = _output_color(objectives[index]["output"], index) if index < len(objectives) else INK
        curve = {"panel": panel, "titles": [panel["title"]], "color": color}
        curves.append(curve)
        by_id[identifier] = curve
    color_uses = Counter()
    for curve in curves:
        curve["linestyle"] = ("-", "--", "-.", ":")[color_uses[curve["color"]] % 4]
        color_uses[curve["color"]] += 1

    # Facets only reduce visual overlap. Each facet uses exactly the same
    # objective axes, normalization and initial reference from the old helper.
    groups = [curves] if len(curves) <= 7 else [curves[i:i+5] for i in range(0, len(curves), 5)]
    ncols = min(2, len(groups))
    nrows = (len(groups) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(10.5 * ncols, 5.4 * nrows), squeeze=False)
    arc = np.linspace(0, np.pi, 181)
    for panel_index, (axis, group) in enumerate(zip(axes.flat, groups, strict=False)):
        for radius in (25, 50, 75, 100):
            axis.plot(radius*np.cos(arc), radius*np.sin(arc), color=GRID, linewidth=0.55, zorder=0)
            axis.text(-radius, -4.5, str(radius), color=MUTED, ha="center", va="top", fontsize=7)
        axis.plot([-100, 100], [0, 0], color=GRID, linewidth=0.55, zorder=0)
        for index, (angle, objective) in enumerate(zip(angles, objectives, strict=True)):
            axis.plot([0, 100*np.cos(angle)], [0, 100*np.sin(angle)], color=GRID, linewidth=0.55, zorder=0)
            horizontal = "right" if np.cos(angle) < -0.7 else "left" if np.cos(angle) > 0.7 else "center"
            label = _radar_metric_label(objective["output"], objective["direction"])
            axis.text(111*np.cos(angle), 111*np.sin(angle), label,
                      color=_output_color(objective["output"], index), fontsize=8,
                      ha=horizontal, va="center")
        values = np.asarray(original_scores, dtype=float)
        initial_label = f"Initial · {initial_id}"
        if initial_roles:
            roles = [title.removeprefix("Best by ").replace("Best average", "average")
                     for title in dict.fromkeys(initial_roles)]
            initial_label += "\nAlso best: " + ", ".join(roles)
        axis.plot(values*np.cos(angles), values*np.sin(angles), color=INITIAL,
                  linestyle="--", linewidth=1.25, label=initial_label, zorder=2)
        for curve in group:
            panel = curve["panel"]
            values = np.asarray(panel["scores"], dtype=float)
            roles = [title.removeprefix("Best by ").replace("Best average", "average")
                     for title in dict.fromkeys(curve["titles"])]
            label = textwrap.fill(f"Best: {', '.join(roles)} · {panel['candidate_id']}",
                                  width=58, break_long_words=False, break_on_hyphens=False)
            axis.plot(values*np.cos(angles), values*np.sin(angles), color=curve["color"],
                      linestyle=curve["linestyle"], linewidth=1.5, label=label, zorder=3)
        top = max(100.0, float(np.max(111*np.sin(angles)))) + 12
        axis.set(xlim=(-145, 145), ylim=(-17, top), aspect="equal")
        axis.axis("off")
        if len(groups) > 1:
            axis.set_title(f"{chr(97 + panel_index)}  Selected candidates", loc="left", fontweight="bold")
        axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02),
                    ncol=2 if len(group) > 1 else 1, fontsize=8, handlelength=2.6,
                    columnspacing=2.0, labelspacing=0.8)
    for axis in axes.flat[len(groups):]:
        axis.set_visible(False)
    fig.suptitle("Selected candidates across optimized metrics", x=0.065, ha="left", fontweight="bold")
    fig.text(0.065, 0.025,
             "Final Pareto + initial min–max desirability (0–100; 100 is better). Arrows indicate raw objective direction.\n"
             "Shared candidates are drawn once; colours follow scorer order in the run configuration.",
             ha="left", fontsize=8, color=MUTED)
    fig.subplots_adjust(top=0.93, bottom=0.26, left=0.08, right=0.92, hspace=0.6)
    path = output_dir / "metric_winners_spider_plots.png"
    _save_plot(fig, path)
    plt.close(fig)
    return path


def _metric_winner_radar_data(
    summary: dict[str, Any], rows: list[dict[str, Any]]
) -> tuple[list[dict[str, str]], list[dict[str, Any]], list[float]] | None:
    objectives = [
        {"output": str(item["output"]), "direction": str(item["direction"])}
        for item in summary.get("optimizer", {}).get("objectives", [])
        if item.get("direction") in {"min", "max"}
    ]
    if len(objectives) < 3 or not rows:
        return None

    candidates: list[dict[str, Any]] = []
    for row in rows:
        outputs = {
            objective["output"]: _float(row[f"output_{objective['output']}"])
            for objective in objectives
            if _is_finite_float(row.get(f"output_{objective['output']}"))
        }
        if len(outputs) == len(objectives):
            candidates.append(
                {
                    "candidate_id": str(row.get("candidate_id", "candidate")),
                    "outputs": outputs,
                }
            )
    if not candidates:
        return None

    original = summary.get("original", {})
    original_outputs = original.get("outputs", {})
    if not all(_is_finite_float(original_outputs.get(item["output"])) for item in objectives):
        return None
    original_record = {
        "candidate_id": str(original.get("candidate_id", "original")),
        "outputs": {
            item["output"]: _float(original_outputs[item["output"]])
            for item in objectives
        },
    }

    bounds: dict[str, tuple[float, float]] = {}
    for objective in objectives:
        output = objective["output"]
        values = [record["outputs"][output] for record in candidates]
        values.append(original_record["outputs"][output])
        bounds[output] = (min(values), max(values))

    def with_scores(record: dict[str, Any]) -> dict[str, Any]:
        scores = [
            _normalize_radar_score(
                record["outputs"][objective["output"]],
                *bounds[objective["output"]],
                objective["direction"],
            )
            for objective in objectives
        ]
        return {
            **record,
            "scores": scores,
            "mean_score": float(np.mean(scores)),
        }

    scored_candidates = [with_scores(record) for record in candidates]
    scored_original = with_scores(original_record)
    panels: list[dict[str, Any]] = []
    for index, objective in enumerate(objectives):
        winner = max(scored_candidates, key=lambda record: record["scores"][index])
        panels.append(
            {
                **winner,
                "title": f"Best by {_radar_metric_name(objective['output'])}",
                "is_original": False,
            }
        )

    balanced = max(scored_candidates, key=lambda record: record["mean_score"])
    panels.append({**balanced, "title": "Best average", "is_original": False})
    panels.append(
        {**scored_original, "title": "Original sequence", "is_original": True}
    )
    return objectives, panels, scored_original["scores"]


def _normalize_radar_score(
    value: float, lower: float, upper: float, direction: str
) -> float:
    if lower == upper:
        # A constant objective is tied at the run-wide optimum for every candidate.
        return 100.0
    if direction == "min":
        score = (upper - value) / (upper - lower)
    else:
        score = (value - lower) / (upper - lower)
    return float(np.clip(score * 100.0, 0.0, 100.0))


def _radar_metric_name(output: str) -> str:
    labels = {
        "mfe_per_nt": "MFE/nt",
        "rinalmo_window_logit_p95": "RiNALMo p95",
        "ribodecode_translation_log": "RiboDecode",
    }
    return labels.get(output, _metric_label(output))


def _radar_metric_label(output: str, direction: str) -> str:
    arrow = "↓" if direction == "min" else "↑"
    return f"{_radar_metric_name(output)} {arrow}"


@_publication_plot
def _plot_module_frequencies(
    output_dir: Path, rows: list[dict[str, Any]], summary: dict[str, Any] | None = None
) -> Path | None:
    if not rows:
        return None
    generation_key = "chain_generation" if "chain_generation" in rows[0] else "generation"
    final_generation = max(int(float(row[generation_key])) for row in rows)
    final_rows = [
        row for row in rows if int(float(row[generation_key])) == final_generation
    ]
    modules = [
        ("five_prime_utr", "5′UTR"),
        ("three_prime_utr", "3′UTR"),
        ("stop_codon", "Stop codon"),
        ("poly_a_length", "poly(A) length"),
    ]
    part_lengths = (summary or {}).get("input", {}).get("part_max_lengths")
    if isinstance(part_lengths, dict):
        modules = [(key, title) for key, title in modules
                   if part_lengths.get("poly_a" if key == "poly_a_length" else key, 0) > 0]
    if not modules:
        for suffix in (".png", ".pdf"):
            (output_dir / ("final_population_module_choices" + suffix)).unlink(missing_ok=True)
        return None

    fig, axes = (plt.subplots(2, 2, figsize=(11, 7)) if len(modules) == 4 else
                 plt.subplots(1, len(modules), figsize=(4 * len(modules), 3.8), squeeze=False))
    for axis, (key, title) in zip(axes.flat, modules, strict=True):
        counts = Counter(
            str(row[key]) for row in final_rows if row.get(key) not in {None, ""}
        )
        labels = sorted(counts)
        axis.barh(labels, [counts[label] for label in labels], color=COLORS[0], alpha=0.60, height=0.55)
        axis.set_title(title)
        axis.set_xlabel("Candidates in final population")
        axis.grid(axis="x", alpha=0.25)
    fig.suptitle(f"Module choices at generation {final_generation}")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = output_dir / "final_population_module_choices.png"
    _save_plot(fig, path)
    plt.close(fig)
    return path


def _annotate_best_objective_delta(
    axis: Any,
    metric: str,
    summary: dict[str, Any],
    x_position: int,
) -> None:
    if not metric.startswith("output_"):
        return
    output_name = metric.removeprefix("output_")
    original_value = summary.get("original", {}).get("outputs", {}).get(output_name)
    best_record = summary.get("best_by_objective", {}).get(output_name)
    if not isinstance(original_value, int | float) or not best_record:
        return
    best_value = best_record.get("outputs", {}).get(output_name)
    if not isinstance(best_value, int | float):
        return

    delta = float(best_value) - float(original_value)
    axis.scatter(
        [x_position],
        [float(best_value)],
        marker="o",
        s=24,
        color=_output_color(metric),
        edgecolor="white",
        linewidth=0.4,
        zorder=6,
    )
    axis.annotate(
        f"delta {_format_delta(delta)}",
        xy=(x_position, float(best_value)),
        xytext=(-8, 0),
        textcoords="offset points",
        ha="right",
        va="center",
        fontsize=9,
        color=INK,
        bbox={
            "boxstyle": "round,pad=0.2",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.75,
        },
        zorder=7,
    )


@_publication_plot
def _plot_codon_changes(
    output_dir: Path, summary: dict[str, Any], evaluations: list[dict[str, Any]]
) -> Path | None:
    by_id = {row["candidate_id"]: row for row in evaluations}
    original_id = summary["original"]["candidate_id"]
    if original_id not in by_id:
        return None
    original_sequence = by_id[original_id].get(
        "coding_sequence", by_id[original_id]["sequence"]
    )

    labels: list[str] = []
    winner_metrics: list[str] = []
    matrix: list[list[int]] = []
    output_order = {name: index for index, name in
                    enumerate((_RUN_COLORS.get() or {}).get("outputs", {}))}
    ordered_winners = sorted(summary["best_by_objective"].items(),
                             key=lambda item: output_order.get(item[0], len(output_order)))
    for metric, record in ordered_winners:
        candidate_id = record["candidate_id"]
        if candidate_id not in by_id:
            continue
        sequence = by_id[candidate_id].get(
            "coding_sequence", by_id[candidate_id]["sequence"]
        )
        labels.append(f"best {metric}")
        winner_metrics.append(metric)
        matrix.append(
            [
                int(original_sequence[i : i + 3] != sequence[i : i + 3])
                for i in range(0, len(original_sequence), 3)
            ]
        )

    if not matrix:
        return None

    data = np.asarray(matrix)
    fig, axis = plt.subplots(figsize=(10.5, max(2.1, 1.3 + 0.42 * len(labels))))
    codons = data.shape[1]
    axis.text(codons+1.5, -0.43, "Changed / total", color=MUTED, fontsize=8,
              ha="left", va="center")
    for index, (row, metric) in enumerate(zip(data, winner_metrics, strict=True)):
        axis.add_patch(Rectangle((0.5, index-0.10), codons, 0.20, facecolor=FIXED_FILL, linewidth=0))
        for position in np.flatnonzero(row):
            axis.add_patch(Rectangle((position+0.5, index-0.10), 1, 0.20,
                                     facecolor=_output_color(metric, index), linewidth=0))
        axis.text(codons+1.5, index, f"{int(row.sum())}/{codons}", va="center", color=MUTED, fontsize=8)
    axis.set_yticks(range(len(labels)))
    axis.set_yticklabels(labels)
    for index, label in enumerate(axis.get_yticklabels()):
        label.set_color(_output_color(winner_metrics[index], index))
    axis.set_xlabel("Codon position (start + CDS, 1-based)")
    axis.set_title("Coding sequence (start + CDS): codon changes")
    axis.set_xticks(np.linspace(1, codons, min(8, codons), dtype=int))
    axis.set(xlim=(0.5, codons+max(8, codons*0.09)), ylim=(len(labels)-0.55, -0.65))
    axis.tick_params(axis="y", length=0)
    axis.spines["left"].set_visible(False)
    fig.tight_layout()

    path = output_dir / "codon_changes_best_candidates.png"
    _save_plot(fig, path)
    plt.close(fig)
    return path


@_publication_plot
def _plot_scorer_timings(output_dir: Path, rows: list[dict[str, Any]]) -> Path | None:
    if not rows or not _has_numeric(rows, "elapsed_seconds"):
        return None

    invocation = np.arange(1, len(rows) + 1)
    elapsed = np.asarray([_float(row["elapsed_seconds"]) for row in rows])
    batch_size = np.asarray([_float(row["batch_size"]) for row in rows])
    seconds_per_candidate = np.asarray(
        [_float(row["seconds_per_candidate"]) for row in rows]
    )
    scorer_names = [row.get("scorer", "scorer") for row in rows]
    scorer_to_color = {
        scorer: _scorer_color(scorer, index)
        for index, scorer in enumerate(dict.fromkeys(scorer_names))
    }
    colors = [scorer_to_color[scorer] for scorer in scorer_names]

    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    axes[0].bar(invocation, elapsed, color=colors, alpha=0.72)
    axes[0].set_ylabel("Elapsed seconds")
    axes[0].set_title("Scorer invocation timing")
    axes[0].grid(True, axis="y", alpha=0.25)
    handles = [
        plt.Line2D([0], [0], color=color, linewidth=2.2, label=scorer)
        for scorer, color in scorer_to_color.items()
    ]
    axes[0].legend(handles=handles, loc="best", fontsize=7)

    twin = axes[0].twinx()
    twin.plot(invocation, batch_size, color=INK, marker="o", markersize=2.5,
              markevery=max(1, len(invocation)//12), linewidth=0.9)
    twin.set_ylabel("Batch size")

    axes[1].scatter(
        invocation,
        seconds_per_candidate,
        color=colors,
        s=12,
    )
    axes[1].plot(invocation, seconds_per_candidate, color=MUTED, linewidth=0.8, alpha=0.5)
    axes[1].set_xlabel("Scorer invocation")
    axes[1].set_ylabel("Seconds per candidate")
    axes[1].grid(True, alpha=0.25)
    _add_timing_run_boundaries(axes, rows)

    fig.tight_layout()
    path = output_dir / "scorer_timings.png"
    _save_plot(fig, path)
    plt.close(fig)
    return path


@_publication_plot
def _plot_scorer_timing_breakdown(output_dir: Path, rows: list[dict[str, Any]]) -> Path | None:
    if not rows or not _has_numeric(rows, "elapsed_seconds"):
        return None

    by_scorer: dict[str, dict[str, float]] = defaultdict(
        lambda: {"seconds": 0.0, "candidates": 0.0, "invocations": 0.0}
    )
    for row in rows:
        scorer = row.get("scorer", "scorer")
        by_scorer[scorer]["seconds"] += _float(row["elapsed_seconds"])
        by_scorer[scorer]["candidates"] += _float(row["batch_size"])
        by_scorer[scorer]["invocations"] += 1

    scorers = list(by_scorer)
    total_seconds = np.asarray([by_scorer[name]["seconds"] for name in scorers])
    seconds_per_candidate = np.asarray(
        [
            by_scorer[name]["seconds"] / by_scorer[name]["candidates"]
            if by_scorer[name]["candidates"]
            else 0.0
            for name in scorers
        ]
    )

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    colors = [_scorer_color(name, index) for index, name in enumerate(scorers)]
    axes[0].bar(scorers, total_seconds, color=colors, alpha=0.76)
    axes[0].set_title("Total scorer time")
    axes[0].set_ylabel("Seconds")
    axes[0].tick_params(axis="x", rotation=20)
    axes[0].grid(True, axis="y", alpha=0.25)

    axes[1].bar(scorers, seconds_per_candidate, color=colors, alpha=0.76)
    axes[1].set_title("Mean time per candidate")
    axes[1].set_ylabel("Seconds / candidate")
    axes[1].tick_params(axis="x", rotation=20)
    axes[1].grid(True, axis="y", alpha=0.25)

    fig.suptitle("Scorer timing breakdown")
    fig.tight_layout()
    path = output_dir / "scorer_timing_breakdown.png"
    _save_plot(fig, path)
    plt.close(fig)
    return path


def _population_metrics(
    rows: list[dict[str, Any]], summary: dict[str, Any] | None, metric_group: str = "all"
) -> list[tuple[str, str]]:
    metrics: list[tuple[str, str]] = []
    if summary is not None:
        include_objectives = metric_group in {"all", "objective"}
        include_measured = metric_group in {"all", "measured"}
        if include_objectives:
            for objective in summary.get("optimizer", {}).get("objectives", []):
                output = objective["output"]
                metrics.append((f"output_{output}", f"{output} ({objective['direction']})"))
        if include_measured:
            for item in summary.get("output_policy", {}).get("thresholds", []):
                output = item["output"]
                metrics.append((f"output_{output}", output))
            for item in summary.get("output_policy", {}).get("logs", []):
                output = item["output"]
                if output.endswith("_structure") or output.endswith("_count"):
                    continue
                metrics.append((f"output_{output}", output))
    else:
        for key in rows[0]:
            if key.startswith("output_"):
                metrics.append((key, key.removeprefix("output_")))

    if metric_group in {"all", "measured"}:
        metrics.append(("changed_codons", "Changed codons"))
    deduped: list[tuple[str, str]] = []
    seen: set[str] = set()
    for key, label in metrics:
        if key in seen:
            continue
        seen.add(key)
        deduped.append((key, _metric_label(label)))
    return deduped[:7]


def _population_plot_metadata(metric_group: str) -> tuple[str, str]:
    if metric_group == "objective":
        return "Optimized parameter dynamics", "population_objective_history.png"
    if metric_group == "measured":
        return "Measured parameter dynamics", "population_measured_history.png"
    return "Population metric dynamics", "population_metric_history.png"


def _metric_label(label: str) -> str:
    labels = {
        "mfe": "MFE, kcal/mol",
        "mfe_per_nt": "MFE per nt, kcal/mol",
        "ensemble_diversity": "Ensemble diversity",
        "gc_fraction": "GC fraction",
        "ensemble_free_energy": "Ensemble free energy, kcal/mol",
        "changed_codons": "Changed codons",
    }
    return labels.get(label, label.replace("_", " "))


def _format_delta(value: float) -> str:
    abs_value = abs(value)
    if abs_value >= 100:
        text = f"{value:+.0f}"
    elif abs_value >= 10:
        text = f"{value:+.1f}"
    elif abs_value >= 1:
        text = f"{value:+.2f}"
    else:
        text = f"{value:+.3f}"
    return text


def _add_run_boundaries(axis: Any, summary: dict[str, Any]) -> None:
    runs = summary.get("run_chain", {}).get("runs", [])
    if len(runs) < 2:
        return
    for index, run in enumerate(runs[1:], start=2):
        start = run.get("generation_start")
        if start is None:
            continue
        x = float(start) - 0.5
        axis.axvline(x, color=MUTED, linewidth=1, linestyle=":", alpha=0.7)
        axis.text(
            x + 0.1,
            0.98,
            f"resume {index}",
            transform=axis.get_xaxis_transform(),
            rotation=90,
            va="top",
            ha="left",
            fontsize=8,
            color=MUTED,
        )


def _add_timing_run_boundaries(axes: Any, rows: list[dict[str, Any]]) -> None:
    if not rows or "run_name" not in rows[0]:
        return
    run_names = {row.get("run_name", "") for row in rows} - {""}
    if len(run_names) < 2:
        return
    first_run = rows[0].get("run_name", "")
    seen: set[str] = {first_run}
    resume_index = 1
    for index, row in enumerate(rows, start=1):
        run_name = row.get("run_name", "")
        if not run_name or run_name in seen:
            continue
        seen.add(run_name)
        resume_index += 1
        for axis in axes:
            axis.axvline(index - 0.5, color=MUTED, linewidth=1, linestyle=":", alpha=0.7)
            axis.text(
                index - 0.4,
                0.98,
                f"resume {resume_index}",
                transform=axis.get_xaxis_transform(),
                rotation=90,
                va="top",
                ha="left",
                fontsize=8,
                color=MUTED,
            )


def _group_by_generation(
    rows: list[dict[str, Any]],
) -> dict[int, list[dict[str, Any]]]:
    grouped: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        generation = _row_generation(row)
        if generation is not None:
            grouped[generation].append(row)
    return grouped


def _row_generation(row: dict[str, Any]) -> int | None:
    for key in ("chain_generation", "generation"):
        if _is_float(row.get(key)):
            return int(float(row[key]))
    return None


def _readable_generation_ticks(generations: list[int]) -> list[int]:
    if len(generations) <= 12:
        return generations
    indices = np.linspace(0, len(generations) - 1, num=10, dtype=int)
    return sorted({generations[int(index)] for index in indices})


def _threshold_numeric_bounds(
    threshold: dict[str, Any],
) -> tuple[float | None, float | None, float | None]:
    return (
        _optional_float(threshold.get("min")),
        _optional_float(threshold.get("max")),
        _optional_float(threshold.get("equals")),
    )


def _optional_float(value: Any) -> float | None:
    if value in {None, ""}:
        return None
    if isinstance(value, bool):
        return float(value)
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _padded_limits(values: list[float]) -> tuple[float, float]:
    lower = float(np.min(values))
    upper = float(np.max(values))
    span = upper - lower
    if np.isclose(span, 0.0):
        padding = max(abs(lower) * 0.1, 0.1)
    else:
        padding = span * 0.08
    y_min = lower - padding
    y_max = upper + padding
    if lower >= 0 and upper <= 1:
        y_min = max(0.0, y_min)
        y_max = min(1.0, y_max)
    if np.isclose(y_min, y_max):
        y_max = y_min + 1.0
    return y_min, y_max


def _shade_allowed_threshold_region(
    axis: Any,
    threshold: dict[str, Any],
    y_min: float,
    y_max: float,
) -> None:
    minimum, maximum, equals = _threshold_numeric_bounds(threshold)
    if equals is not None:
        return
    allowed_min = minimum if minimum is not None else y_min
    allowed_max = maximum if maximum is not None else y_max
    axis.axhspan(
        max(y_min, allowed_min),
        min(y_max, allowed_max),
        color=FIXED_FILL,
        alpha=0.7,
        zorder=0,
    )


def _draw_threshold_lines(axis: Any, threshold: dict[str, Any]) -> None:
    minimum, maximum, equals = _threshold_numeric_bounds(threshold)
    if minimum is not None:
        axis.axhline(
            minimum,
            color=MUTED,
            linewidth=0.9,
            linestyle="--",
            label="threshold",
        )
    if maximum is not None:
        axis.axhline(
            maximum,
            color=MUTED,
            linewidth=0.9,
            linestyle="--",
            label="threshold",
        )
    if equals is not None:
        axis.axhline(
            equals,
            color=MUTED,
            linewidth=0.9,
            linestyle="--",
            label="threshold",
        )


def _violation_marker_value(
    failed_values: list[float], threshold: dict[str, Any]
) -> float:
    minimum, maximum, equals = _threshold_numeric_bounds(threshold)
    below = [value for value in failed_values if minimum is not None and value < minimum]
    if below:
        return min(below)
    above = [value for value in failed_values if maximum is not None and value > maximum]
    if above:
        return max(above)
    if equals is not None:
        return max(failed_values, key=lambda value: abs(value - equals))
    return failed_values[0]


def _threshold_axis_label(output: str) -> str:
    if output in {"body_gc_window_min", "body_gc_window_max"}:
        return "GC fraction"
    if output == "non_polya_max_homopolymer":
        return "Nucleotides"
    return _metric_label(output)


def _read_csv(path: Path) -> list[dict[str, Any]]:
    if not path.exists() or not path.read_text().strip():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _has_numeric(rows: list[dict[str, Any]], key: str) -> bool:
    return any(_is_float(row.get(key)) for row in rows)


def _is_float(value: Any) -> bool:
    if value is None:
        return False
    try:
        float(value)
    except ValueError:
        return False
    return True


def _is_finite_float(value: Any) -> bool:
    return _is_float(value) and bool(np.isfinite(_float(value)))


def _float(value: Any) -> float:
    return float(value)


def _parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}
