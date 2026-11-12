"""Read saved MORSE results without loading optimizer checkpoints."""
from __future__ import annotations

import copy
import csv
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
import statistics
import sys
from typing import Any


@dataclass(frozen=True)
class Candidate:
    """One selected construct, with RNA sequences and its recorded scores."""

    id: str
    sequence: str
    parts: dict[str, Any]
    metrics: dict[str, Any]

    @property
    def candidate_id(self) -> str:
        return self.id


def _read_json(path: Path, default=None):
    return json.loads(path.read_text()) if path.is_file() else copy.deepcopy(default)


def _value(value: str) -> Any:
    if value == "":
        return None
    if value.lower() in {"true", "false"}:
        return value.lower() == "true"
    try:
        return json.loads(value)
    except (ValueError, TypeError):
        return value


def _finite(value) -> bool:
    return not isinstance(value, bool) and isinstance(value, (float, int)) and math.isfinite(value)


def _rna(sequence: str) -> str:
    return sequence.upper().replace("T", "U")


def _ensure_notebook_figure_formatter() -> None:
    """Enable returned figures in fresh kernels without changing user formats."""
    ipython = sys.modules.get("IPython")
    shell = getattr(ipython, "get_ipython", lambda: None)()
    if getattr(shell, "kernel", None) is None:
        return
    from matplotlib.figure import Figure
    formatters = shell.display_formatter.formatters
    for mime, formatter in formatters.items():
        if mime.startswith("image/") or mime == "application/pdf":
            try:
                formatter.lookup_by_type(Figure)
            except KeyError:
                continue
            return
    from IPython.core.pylabtools import print_figure
    formatters["image/png"].for_type(Figure, lambda figure: print_figure(figure, fmt="png"))


class OptimizationResult:
    """A lightweight view of a completed or interrupted run directory.

    ``pareto`` and ``history`` are lists of dictionaries and can be passed to
    ``pandas.DataFrame``. Properties reread saved artifacts, so ``load_run`` can
    also be used to inspect a run that is still writing progress.
    """

    def __init__(self, output_dir: str | Path):
        self.output_dir = Path(output_dir).expanduser().resolve()
        if not self.output_dir.is_dir():
            raise FileNotFoundError(f"Run directory does not exist: {self.output_dir}")
        if not any((self.output_dir / name).is_file() for name in (
            "summary.json", "progress.json", "optimizer-checkpoint.pkl", "run_state.json"
        )):
            raise ValueError(f"No MORSE run artifacts found in {self.output_dir}")

    @property
    def summary(self) -> dict[str, Any]:
        summary = _read_json(self.output_dir / "summary.json")
        if summary is not None:
            return summary
        progress = _read_json(self.output_dir / "progress.json", {})
        return {"status": "incomplete", "output_dir": str(self.output_dir), "progress": progress}

    @property
    def pareto(self) -> list[dict[str, Any]]:
        path = self.output_dir / "pareto.csv"
        if not path.is_file():
            return []
        with path.open(newline="") as handle:
            rows = [{key: _value(value) for key, value in row.items()}
                    for row in csv.DictReader(handle)]
        for row in rows:
            if row.get("sequence"):
                row["sequence"] = _rna(row["sequence"])
        return rows

    def _population_rows(self):
        path = self.output_dir / "run_log.jsonl.gz"
        if path.is_file():
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                for line in handle:
                    event = json.loads(line)
                    if event.get("type") == "population_member":
                        yield event

    def _objective_specs(self) -> list[dict[str, Any]]:
        return self.summary.get("optimizer", {}).get("objectives", [])

    def _has_utr(self, field: str) -> bool:
        lengths = self.summary.get("input", {}).get("part_max_lengths", {})
        if field in lengths:
            return lengths[field] > 0
        for name in ("resolved_run.yaml", "experiment.yaml"):
            path = self.output_dir / name
            if path.is_file():
                import yaml
                config = yaml.safe_load(path.read_text())
                construct = config.get("construct", {}) if isinstance(config, dict) else {}
                coordinates = construct.get("coordinates")
                if isinstance(coordinates, dict):
                    return field in coordinates or field in construct.get("changes", {})
                break
        # Older artifacts may not retain enough information to identify absent UTRs.
        return True

    def _legacy_history(self) -> list[dict[str, Any]]:
        grouped: dict[int, list[dict[str, Any]]] = {}
        for row in self._population_rows():
            grouped.setdefault(int(row["generation"]), []).append(row)
        specs = {item["output"]: item for item in self._objective_specs()}
        history = []
        for generation, rows in sorted(grouped.items()):
            outputs = sorted({key[7:] for row in rows for key, value in row.items()
                              if key.startswith("output_") and _finite(value)})
            objectives, logged = [], []
            for name in outputs:
                values = [row["output_" + name] for row in rows if _finite(row.get("output_" + name))]
                if not values:
                    continue
                direction = specs.get(name, {}).get("direction")
                metric = {"name": name, "count": len(values), "mean": statistics.fmean(values),
                          "std": statistics.pstdev(values), "median": statistics.median(values),
                          "observed_min": min(values), "observed_max": max(values),
                          "best": min(values) if direction == "min" else max(values)}
                if direction:
                    metric["direction"] = direction
                    objectives.append(metric)
                else:
                    logged.append(metric)
            frequencies = {}
            for field in ("five_prime_utr", "three_prime_utr"):
                counts: dict[str, int] = {}
                for row in rows:
                    if row.get(field) is not None:
                        name = str(row[field])
                        counts[name] = counts.get(name, 0) + 1
                frequencies[field] = counts
            feasible = sum(row.get("valid") is True or str(row.get("valid")).lower() == "true"
                           for row in rows)
            history.append({"generation": generation, "population_size": len(rows),
                            "hard_feasible_count": feasible, "hard_feasible_fraction": feasible / len(rows),
                            "objectives": objectives, "logged_metrics": logged,
                            "utr_frequencies": frequencies})
        return history

    @property
    def history(self) -> list[dict[str, Any]]:
        progress = _read_json(self.output_dir / "progress.json", {})
        history = copy.deepcopy(progress.get("history", []))
        if not history and progress.get("generation") is not None:
            history = [copy.deepcopy(progress)]
        needs_legacy = not history or any(
            "utr_frequencies" not in row or any("std" not in metric for group in
                ("objectives", "thresholds", "logged_metrics") for metric in row.get(group, []))
            for row in history)
        if needs_legacy:
            legacy = self._legacy_history()
            if not history:
                history = legacy
            elif legacy:
                by_generation = {row["generation"]: row for row in legacy}
                # Completed continuation logs use global generation indices.
                runs = self.summary.get("run_chain", {}).get("runs", [])
                offset = int(runs[-1].get("generation_start", 1)) - 1 if runs else 0
                for row in history:
                    source = by_generation.get(int(row["generation"]) + offset)
                    if source is None:
                        continue
                    for key in ("population_size", "hard_feasible_count", "utr_frequencies"):
                        row.setdefault(key, copy.deepcopy(source[key]))
                    metrics = {item["name"]: item for group in ("objectives", "logged_metrics")
                               for item in source.get(group, [])}
                    for group in ("objectives", "thresholds", "logged_metrics"):
                        for metric in row.get(group, []):
                            if metric["name"] in metrics:
                                metric.setdefault("std", metrics[metric["name"]]["std"])
        for row in history:
            if "population_size" not in row and progress.get("population_size") is not None:
                row["population_size"] = progress["population_size"]
        saved_hv = self.summary.get("optimizer", {}).get("hypervolume", {}).get("history", [])
        if not history:
            history = copy.deepcopy(saved_hv)
        elif saved_hv:
            by_generation = {int(row["generation"]): row for row in saved_hv}
            for row in history:
                for key, value in by_generation.get(int(row["generation"]), {}).items():
                    if key in {"hypervolume", "best_hypervolume"}:
                        row.setdefault(key, value)
        return history

    def _candidate(self, row: dict[str, Any]) -> Candidate:
        identifier = str(row["candidate_id"])
        entries = _read_json(self.output_dir / "pareto_constructs.json", [])
        matching = [item for item in entries if item["candidate_id"] == identifier]
        if len(matching) > 1:
            raise ValueError(f"Duplicate construct record: {identifier}")
        parts = copy.deepcopy(matching[0]["construct"]) if matching else {}
        sequence = str(row.get("sequence") or parts.get("sequence") or "")
        if not sequence:
            raise ValueError(f"Candidate has no saved sequence: {identifier}")
        if parts.get("sequence") and _rna(parts["sequence"]) != _rna(sequence):
            raise ValueError(f"Saved candidate sequence disagrees with its parts: {identifier}")
        for key in ("sequence", "start_codon", "cds", "stop_codon"):
            if isinstance(parts.get(key), str):
                parts[key] = _rna(parts[key])
        for key in ("five_prime_utr", "three_prime_utr"):
            if isinstance(parts.get(key), dict) and "sequence" in parts[key]:
                parts[key]["sequence"] = _rna(parts[key]["sequence"])
        metrics = {key[7:]: value for key, value in row.items() if key.startswith("output_")}
        return Candidate(identifier, _rna(sequence), parts, metrics)

    def select(self, method: str = "hvc", *, metric: str | None = None) -> Candidate:
        """Select from the final feasible Pareto set by HVC or one objective.

        HVC uses min-max scaling over that set plus the starting construct and
        reference 1.1 on every cost axis. Equal contributions use candidate ID.
        """
        rows = self.pareto
        if not rows:
            raise ValueError("No feasible Pareto candidates are available")
        ids = [str(row["candidate_id"]) for row in rows]
        if len(ids) != len(set(ids)):
            raise ValueError("Pareto candidates must have unique candidate IDs")
        for row in rows:
            flags = [value for key, value in row.items() if key in {"valid", "structural_valid"}
                     or key.startswith("threshold_") and key.endswith("_valid")]
            if any(value is not True for value in flags):
                raise ValueError("Pareto set contains a non-feasible candidate")
        specs = self._objective_specs()
        if metric is not None:
            spec = next((item for item in specs if item["output"] == metric), None)
            if spec is None or spec.get("direction") not in {"min", "max"}:
                raise ValueError(f"Unknown optimization metric: {metric}")
            values = [row.get("output_" + metric) for row in rows]
            if not all(_finite(value) for value in values):
                raise ValueError(f"Missing or non-finite objective values: {metric}")
            sign = 1 if spec["direction"] == "min" else -1
            index = min(range(len(rows)), key=lambda i: (sign * values[i], ids[i]))
            return self._candidate(rows[index])
        if method != "hvc":
            raise ValueError("Selection method must be 'hvc', or provide metric='scorer.metric'")
        if not specs or any(item.get("direction") not in {"min", "max"} for item in specs):
            raise ValueError("Saved result does not define optimization objectives")
        import numpy as np
        from pymoo.indicators.hv import HV
        names = [item["output"] for item in specs]
        original = self.summary.get("original", {}).get("outputs", {})
        native = [[row.get("output_" + name) for name in names] for row in rows]
        reference_values = [original.get(name) for name in names]
        if not all(_finite(value) for row in native for value in row) or not all(_finite(value) for value in reference_values):
            raise ValueError("HVC selection requires finite Pareto and starting-construct objective values")
        signs = np.asarray([1 if item["direction"] == "min" else -1 for item in specs])
        costs = np.asarray(native, dtype=float) * signs
        pool = np.vstack([costs, np.asarray(reference_values) * signs])
        ideal = pool.min(axis=0)
        span = pool.max(axis=0) - ideal
        span = np.where(span == 0, 1.0, span)
        if not np.isfinite(span).all():
            raise ValueError("HVC normalization has non-finite objective spans")
        points = (costs - ideal) / span
        hv = HV(ref_point=np.full(len(names), 1.1))
        total = float(hv(points))
        contributions = [total - float(hv(np.delete(points, i, axis=0))) if len(rows) > 1 else total
                         for i in range(len(rows))]
        if not all(math.isfinite(value) and value >= -1e-12 for value in contributions):
            raise ValueError("Invalid hypervolume contributions")
        index = min(range(len(rows)), key=lambda i: (-max(0, contributions[i]), ids[i]))
        return self._candidate(rows[index])

    def plot(self, kind: str = "metrics", *, metrics: list[str] | None = None, top: int = 5):
        """Return a Matplotlib figure from recorded generations, without scoring."""
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
        from morse.publication_style import GRID, INK, MUTED, RC_STYLE, palette_color
        from morse.plots import _read_run_colors
        colors = _read_run_colors(self.output_dir, self.summary).get("outputs", {})
        history = self.history
        if not history:
            raise ValueError("No population history is available")
        x = [row["generation"] for row in history]
        with plt.rc_context(RC_STYLE):
            if kind == "metrics":
                names = list(dict.fromkeys(item["name"] for row in history
                             for group in ("objectives", "thresholds", "logged_metrics")
                             for item in row.get(group, []) if item.get("mean") is not None))
                if metrics is not None:
                    unknown = set(metrics) - set(names)
                    if unknown:
                        raise ValueError("Unknown recorded metrics: " + ", ".join(sorted(unknown)))
                    names = list(metrics)
                if not names:
                    raise ValueError("No numeric metric history is available")
                columns = min(2, len(names))
                fig, axes = plt.subplots(math.ceil(len(names) / columns), columns,
                                         figsize=(6 * columns, 3.2 * math.ceil(len(names) / columns)), squeeze=False)
                for index, name in enumerate(names):
                    axis = axes.flat[index]
                    observations = [next((item for group in ("objectives", "thresholds", "logged_metrics")
                                          for item in row.get(group, []) if item["name"] == name), {})
                                    for row in history]
                    means = [item.get("mean") if _finite(item.get("mean")) else math.nan for item in observations]
                    color = colors.get(name, palette_color(index))
                    axis.plot(x, means, color=color, label="Population mean")
                    if any(_finite(item.get("std")) for item in observations):
                        std = [item.get("std") if _finite(item.get("std")) else math.nan for item in observations]
                        axis.fill_between(x, [m-s for m, s in zip(means, std)], [m+s for m, s in zip(means, std)],
                                          color=color, alpha=.18, linewidth=0, label="±1 SD")
                    best = [item.get("best") if _finite(item.get("best")) else math.nan for item in observations]
                    directed = any(item.get("direction") in {"min", "max"} for item in observations)
                    if directed and any(_finite(value) for value in best):
                        axis.plot(x, best, color=color, linestyle="--", linewidth=1, label="Best")
                    axis.set(title=name, xlabel="Generation")
                    axis.legend(fontsize=8)
                for axis in list(axes.flat)[len(names):]:
                    axis.set_visible(False)
            elif kind == "hypervolume":
                if not any(_finite(row.get("hypervolume")) for row in history):
                    raise ValueError("Hypervolume history is not available for this run")
                fig, axis = plt.subplots(figsize=(7, 3.3))
                axis.plot(x, [row.get("hypervolume", math.nan) for row in history], color=INK, label="Hypervolume")
                if any(_finite(row.get("best_hypervolume")) for row in history):
                    axis.plot(x, [row.get("best_hypervolume", math.nan) for row in history], color=MUTED,
                              linestyle="--", label="Best so far")
                axis.set(xlabel="Generation", ylabel="Normalized hypervolume", title="Hypervolume")
                axis.legend()
            elif kind == "utrs":
                if isinstance(top, bool) or not isinstance(top, int) or top < 1:
                    raise ValueError("top must be a positive integer")
                fields = [(field, title) for field, title in (("five_prime_utr", "5′UTR"), ("three_prime_utr", "3′UTR"))
                          if self._has_utr(field) and any(row.get("utr_frequencies", {}).get(field) for row in history)]
                if not fields:
                    raise ValueError("UTR frequency history is not available for this run")
                fig, axes = plt.subplots(1, len(fields), figsize=(6 * len(fields), 3.5), squeeze=False)
                for axis, (field, title) in zip(axes.flat, fields):
                    distributions = [row.get("utr_frequencies", {}).get(field, {}) for row in history]
                    variants = {name for counts in distributions for name in counts}
                    ordered = sorted(variants, key=lambda name: (-max(counts.get(name, 0) for counts in distributions), name))
                    shown = ordered[:top]
                    denominators = [row.get("population_size") or sum(counts.values()) for row, counts in zip(history, distributions)]
                    for index, variant in enumerate(shown):
                        axis.plot(x, [100 * counts.get(variant, 0) / n if n else math.nan
                                      for counts, n in zip(distributions, denominators)],
                                  color=palette_color(index), label=variant)
                    if len(ordered) > top:
                        axis.plot(x, [100 * sum(v for k, v in counts.items() if k not in shown) / n if n else math.nan
                                      for counts, n in zip(distributions, denominators)], color=MUTED, label="Other")
                    axis.set(title=title, xlabel="Generation", ylabel="Share of population (%)", ylim=(-2, 102))
                    axis.legend(fontsize=7)
            else:
                raise ValueError("Plot kind must be 'metrics', 'hypervolume', or 'utrs'")
            for axis in fig.axes:
                if axis.get_visible():
                    axis.grid(axis="y", color=GRID, linewidth=.5)
                    axis.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
            fig.tight_layout()
            _ensure_notebook_figure_formatter()
            return fig

    def report(self) -> Path:
        """Display the saved report in a notebook and return its file path."""
        path = self.output_dir / "report.md"
        if not path.is_file():
            raise FileNotFoundError(f"Report is not available yet: {path}")
        try:
            from IPython import get_ipython
        except ImportError:
            return path
        if get_ipython() is not None:
            from morse.notebook_tools import discover_result_artifacts, display_result_artifacts
            display_result_artifacts(discover_result_artifacts(self.output_dir))
        return path
