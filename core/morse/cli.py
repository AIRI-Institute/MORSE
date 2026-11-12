from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from morse.config import RunConfig, ScorerConfig
from morse.run_spec import RunSpec, load_run_spec, write_run_spec
from morse.scorers import (
    CompositeScorer,
    ConstructCandidate,
    ScorerBinding,
    make_raw_scorer,
)
from morse.validators import StructuralValidator
from morse.visualization import write_sequence_scope_svg

LOGGER = logging.getLogger("morse")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    _configure_logging(args)
    try:
        return int(args.func(args))
    except Exception as exc:
        if LOGGER.isEnabledFor(logging.DEBUG):
            LOGGER.exception("Command failed")
        else:
            LOGGER.error("%s", exc)
        return 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="morse-core",
        description="Configurable RNA sequence optimization.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    init = subparsers.add_parser("init", help="Create a portable CPU experiment from a full coding FASTA.")
    init.add_argument("--cds", type=Path, required=True)
    init.add_argument("--output-dir", type=Path, required=True)
    init.add_argument("--five-prime-utr", type=Path, action="append", default=[])
    init.add_argument("--three-prime-utr", type=Path, action="append", default=[])
    init.add_argument("--seed", type=int, default=43)
    init.add_argument("--generations", type=int, default=10)
    init.add_argument("--population-size", type=int, default=16)
    add_logging_args(init)
    init.set_defaults(func=cmd_init)

    score = subparsers.add_parser(
        "score", help="Score the reference construct from one run config."
    )
    add_config_arg(score)
    add_input_root_arg(score)
    add_logging_args(score)
    score.set_defaults(func=cmd_score)

    pilot = subparsers.add_parser(
        "pilot",
        help="Score the reference plus a reproducible random sample before optimization.",
    )
    add_config_arg(pilot)
    add_input_root_arg(pilot)
    pilot.add_argument("--samples", type=int, default=24)
    pilot.add_argument("--output-dir", type=Path, required=True)
    add_logging_args(pilot)
    pilot.set_defaults(func=cmd_pilot)

    run = subparsers.add_parser(
        "run", help="Execute one experiment; output_dir may be set in schema-4 YAML."
    )
    add_config_arg(run)
    add_input_root_arg(run)
    run.add_argument("--output-dir", type=Path)
    resume = run.add_mutually_exclusive_group()
    resume.add_argument(
        "--resume",
        action="store_true",
        help="Resume from the checkpoint already present in --output-dir.",
    )
    resume.add_argument(
        "--resume-from",
        type=Path,
        help=(
            "Interrupted run directory with an exact checkpoint, or a completed "
            "run whose final population should seed a continuation."
        ),
    )
    add_logging_args(run)
    run.set_defaults(func=cmd_run)

    validate = subparsers.add_parser(
        "validate", help="Validate a config and all referenced inputs without running."
    )
    add_config_arg(validate)
    add_input_root_arg(validate)
    add_logging_args(validate)
    validate.set_defaults(func=cmd_validate)

    preview = subparsers.add_parser(
        "preview", help="Render construct, optimizer, and scorer scopes without running scorers."
    )
    add_config_arg(preview)
    add_input_root_arg(preview)
    preview.add_argument("--output", type=Path, default=Path("sequence_scope.svg"))
    add_logging_args(preview)
    preview.set_defaults(func=cmd_preview)

    report = subparsers.add_parser(
        "report", help="Regenerate plots and report for an existing run directory."
    )
    report.add_argument("output_dir", type=Path)
    add_logging_args(report)
    report.set_defaults(func=cmd_report)

    return parser


def add_config_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", type=Path, required=True)


def add_input_root_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--input-root",
        type=Path,
        default=None,
        help=(
            "Base directory for relative FASTA, profile, and initial-population "
            "paths; defaults to the config directory."
        ),
    )


def add_logging_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--log-level",
        choices=("error", "warning", "info", "debug"),
        default="info",
    )
    parser.add_argument(
        "--log-format",
        choices=("text", "json"),
        default="text",
    )
    parser.add_argument("--log-file", type=Path, default=None)


def cmd_init(args: argparse.Namespace) -> int:
    from morse.experiment import create_experiment
    path = create_experiment(
        args.cds, args.output_dir, five_prime_utrs=args.five_prime_utr,
        three_prime_utrs=args.three_prime_utr, seed=args.seed,
        generations=args.generations, population_size=args.population_size,
    )
    LOGGER.info("Created experiment: %s", path)
    return 0


def cmd_score(args: argparse.Namespace) -> int:
    from morse.execution import score_reference

    spec = load_run_spec(args.config, input_root=args.input_root)
    print(json.dumps(score_reference(spec, scorer_factory=make_scorer), indent=2, sort_keys=True))
    return 0


def cmd_pilot(args: argparse.Namespace) -> int:
    from morse.pilot import run_pilot

    spec = load_run_spec(args.config, input_root=args.input_root)
    scorer = make_scorer(spec.run_config)
    try:
        result = run_pilot(
            spec,
            scorer,
            sample_count=args.samples,
            seed=spec.run_config.optimizer.seed,
            output_dir=args.output_dir,
        )
    finally:
        close = getattr(scorer, "close", None)
        if callable(close):
            close()
    print(json.dumps(result, sort_keys=True))
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    config_path = args.config.expanduser().resolve()
    spec = load_run_spec(config_path, input_root=args.input_root)
    configured_output = spec.resolved.get("output_dir")
    selected_output = args.output_dir or (Path(configured_output) if configured_output else None)
    if selected_output is None:
        raise ValueError("Provide output_dir in the YAML or --output-dir")
    output_dir = selected_output.expanduser().resolve()
    resume_from = (
        output_dir
        if args.resume
        else (
            args.resume_from.expanduser().resolve()
            if args.resume_from is not None
            else None
        )
    )
    result = _run_optimize_and_report(
        spec,
        output_dir,
        run_name=spec.name,
        resume_from=resume_from,
        source_config_path=config_path,
    )
    LOGGER.info("Run directory: %s", result.output_dir)
    LOGGER.info("Report: %s", result.output_dir / "report.md")
    if result.summary.get("status") == "completed_no_feasible":
        LOGGER.error("Search completed without hard-feasible candidates; complete observations are retained in %s", result.output_dir)
        return 1
    return 0


def cmd_validate(args: argparse.Namespace) -> int:
    spec = load_run_spec(args.config, input_root=args.input_root)
    LOGGER.info(
        "Valid config %s: %d objectives, %d thresholds",
        spec.name,
        sum(output.role == "objective" for output in spec.run_config.outputs),
        sum(output.role == "threshold" for output in spec.run_config.outputs),
    )
    if spec.resolved.get("schema_version") == 4:
        from morse.schema4 import scope_expression
        construct = spec.resolved["construct"]
        coordinates = construct["coordinates"]
        LOGGER.info("RNA: %d nt; parts: %s", len(spec.design_space.source_construct.sequence),
                    ", ".join(f"{part}[{start}:{end}]" for part, (start, end) in coordinates.items()))
        LOGGER.info("Changes: %s", ", ".join(construct.get("changes", {})))
        for scorer in spec.resolved["scorers"]:
            policies = []
            for name, policy in scorer["metrics"].items():
                details = [str(policy["role"])]
                if policy.get("direction"):
                    details.append(str(policy["direction"]))
                details.extend(f"{key}={policy[key]}" for key in ("min", "max", "equals") if key in policy)
                policies.append(f"{name}: {' '.join(details)}")
            LOGGER.info("Scorer %s: %s; metrics: %s", scorer["name"],
                        scope_expression(scorer["scope"]),
                        "; ".join(policies))
        LOGGER.info("Search: %s; population=%d, generations=%d, seed=%d",
                    spec.run_config.optimizer.type, spec.run_config.optimizer.population_size,
                    spec.run_config.optimizer.generations, spec.run_config.optimizer.seed)
        if spec.resolved.get("output_dir"):
            LOGGER.info("Results: %s", spec.resolved["output_dir"])
    return 0


def cmd_preview(args: argparse.Namespace) -> int:
    spec = load_run_spec(args.config, input_root=args.input_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    path = write_sequence_scope_svg(args.output, spec)
    LOGGER.info("Sequence plan: %s", path.resolve())
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    from morse.plots import plot_run
    from morse.report import write_report

    config_path = args.output_dir / "resolved_run.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Cannot regenerate report: missing {config_path}")
    spec = load_run_spec(config_path)
    run_config = spec.run_config
    diagram_path = write_sequence_scope_svg(
        args.output_dir / "sequence_scope.svg", spec
    )
    LOGGER.info("Sequence plan: %s", diagram_path)
    plot_paths = plot_run(args.output_dir)
    report_path = write_report(args.output_dir, plot_paths, run_config, config_path)
    for path in plot_paths:
        LOGGER.info("Plot: %s", path)
    LOGGER.info("Report: %s", report_path)
    return 0


def _run_optimize_and_report(
    spec: RunSpec,
    output_dir: Path,
    run_name: str,
    resume_from: Path | None,
    source_config_path: Path | None = None,
) -> Any:
    """Compatibility wrapper around the shared execution path."""
    from morse.execution import execute_run

    return execute_run(
        spec, output_dir, run_name, resume_from, source_config_path,
        scorer_factory=make_scorer,
    )


def make_scorer(run_config: RunConfig) -> Any:
    bindings = [
        ScorerBinding(
            _make_one_scorer(scorer_config),
            scorer_config.scope,
            output_prefix=scorer_config.output_prefix,
            namespace=scorer_config.namespace,
        )
        for scorer_config in run_config.scorers
    ]
    return CompositeScorer(bindings)


def _make_one_scorer(scorer_config: ScorerConfig) -> Any:
    """Compatibility entry point used by existing report and paper helpers."""
    return make_raw_scorer(scorer_config, legacy_python=sys.executable)


class _JsonLogFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        payload: dict[str, Any] = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "level": record.levelname.lower(),
            "message": record.getMessage(),
        }
        if record.exc_info:
            payload["exception"] = self.formatException(record.exc_info)
        return json.dumps(payload, sort_keys=True)


def _configure_logging(args: argparse.Namespace) -> None:
    level = getattr(logging, str(args.log_level).upper())
    formatter: logging.Formatter = (
        _JsonLogFormatter()
        if args.log_format == "json"
        else logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    )
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if args.log_file is not None:
        log_path = args.log_file.expanduser().resolve()
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_path, encoding="utf-8"))
    LOGGER.handlers.clear()
    LOGGER.setLevel(level)
    LOGGER.propagate = False
    for handler in handlers:
        handler.setLevel(level)
        handler.setFormatter(formatter)
        LOGGER.addHandler(handler)
