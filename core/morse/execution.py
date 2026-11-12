"""Shared execution path for the command line and Python API."""
from __future__ import annotations

import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Any

from morse.config import RunConfig, ScorerConfig
from morse.run_spec import RunSpec, write_run_spec
from morse.scorers import CompositeScorer, ConstructCandidate, ScorerBinding, make_raw_scorer
from morse.validators import StructuralValidator
from morse.visualization import write_sequence_scope_svg

LOGGER = logging.getLogger("morse")


def score_reference(spec: RunSpec, *, scorer_factory: Any = None) -> dict[str, Any]:
    """Score one reference with the common Python/CLI contract."""
    run_config = spec.run_config
    design_space = spec.design_space
    structural_validation = StructuralValidator().validate(
        design_space.source_construct, design_space
    )
    scorer = (scorer_factory or make_scorer)(run_config)
    try:
        result = scorer.score([ConstructCandidate("input", design_space.source_construct)])[0]
    finally:
        close = getattr(scorer, "close", None)
        if callable(close):
            close()
    outputs = dict(result.outputs)
    for metric in run_config.outputs:
        if metric.name in outputs:
            outputs[metric.name] = metric.coerce(outputs[metric.name])
        elif metric.role != "log":
            raise ValueError(f"Missing required metric {metric.name!r}")
    thresholds = [metric.check_threshold(outputs) for metric in run_config.outputs if metric.role == "threshold"]
    from dataclasses import asdict
    return {
        "id": result.candidate_id,
        "structural_validation": structural_validation.to_dict(),
        "outputs": outputs,
        "thresholds": [asdict(item) for item in thresholds],
        "feasible": structural_validation.valid and all(item.valid for item in thresholds),
        "config_outputs": [
            {
                "name": spec.name,
                "type": spec.type,
                "role": spec.role,
                "direction": spec.direction,
                "min": spec.min_value,
                "max": spec.max_value,
                "equals": spec.equals_value,
            }
            for spec in run_config.outputs
        ],
    }


def execute_run(
    spec: RunSpec,
    output_dir: Path,
    run_name: str,
    resume_from: Path | None,
    source_config_path: Path | None = None,
    on_progress: Any = None,
    scorer_factory: Any = None,
) -> Any:
    from morse.optimizer import run_optimization
    from morse.validators import validator_for_policy

    if spec.resolved.get("schema_version") == 4:
        from dataclasses import replace
        spec = replace(spec, resolved={**spec.resolved, "output_dir": str(output_dir.resolve())})

    run_config = spec.run_config
    if (
        resume_from is not None
        and not (resume_from / "optimizer-checkpoint.pkl").is_file()
        and run_config.optimizer.constraint_handling.method == "adaptive_epsilon"
    ):
        raise ValueError(
            "Resume is not supported with adaptive_epsilon constraint handling"
        )

    from morse.provenance import runtime_provenance, scoring_fingerprint

    provenance = runtime_provenance(run_config)
    fingerprint = scoring_fingerprint(run_config, provenance)
    provenance["scoring_fingerprint"] = fingerprint
    # Preserve the original evidence if an incompatible resume is rejected.
    if resume_from is not None and (resume_from / "run_provenance.json").is_file():
        previous = json.loads((resume_from / "run_provenance.json").read_text())
        if previous.get("scoring_fingerprint") != fingerprint:
            raise ValueError("Cannot resume: scorer source, model assets, image, or runtime changed")

    _validate_resume_before_write(spec, output_dir, resume_from, fingerprint)

    output_dir.mkdir(parents=True, exist_ok=True)
    if source_config_path is not None:
        (output_dir / "experiment.yaml").write_bytes(
            source_config_path.read_bytes()
        )
    config_snapshot_path = output_dir / "resolved_run.yaml"
    write_run_spec(config_snapshot_path, spec)
    if resume_from is None and run_config.initialization.mode != "random":
        spec = _prepare_configured_initial_population(
            spec,
            config_snapshot_path,
            output_dir,
        )
        run_config = spec.run_config
        write_run_spec(config_snapshot_path, spec)
    diagram_path = write_sequence_scope_svg(output_dir / "sequence_scope.svg", spec)
    LOGGER.info("Sequence plan: %s", diagram_path)
    provenance_path = output_dir / "run_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    # Execute the exact image inspected for provenance, even if its tag moves.
    from dataclasses import replace
    runtime_config = replace(run_config, scorers=[
        replace(config, image=record.get("image_id", config.image))
        for config, record in zip(run_config.scorers, provenance.get("scorers", [{}] * len(run_config.scorers)))
    ])
    scorer = (scorer_factory or make_scorer)(runtime_config)
    try:
        result = run_optimization(
            input_manifest=spec.path,
            output_dir=output_dir,
            scorer=scorer,
            output_specs=run_config.outputs,
            optimizer=run_config.optimizer,
            run_name=run_name,
            resume_from=resume_from,
            scoring_fingerprint=fingerprint,
            design_space=spec.design_space,
            validator=validator_for_policy(run_config.manufacturability),
            on_progress=on_progress,
        )
    finally:
        effective_parameters = getattr(scorer, "resolved_parameters", {})
        if effective_parameters:
            (output_dir / "scorer_parameters.json").write_text(
                json.dumps(effective_parameters, indent=2, sort_keys=True, allow_nan=False) + "\n"
            )
        close = getattr(scorer, "close", None)
        if callable(close):
            close()
    write_run_report(output_dir, run_config, config_snapshot_path)
    return result


def write_run_report(output_dir: Path, run_config: RunConfig, config_path: Path) -> Path:
    """Render recorded observations without launching scorers again."""
    from morse.plots import plot_run
    from morse.report import write_report

    plot_paths = plot_run(output_dir)
    write_report(output_dir, plot_paths, run_config, config_path)
    return output_dir / "report.md"


def _validate_resume_before_write(
    spec: RunSpec, output_dir: Path, resume_from: Path | None, fingerprint: str
) -> None:
    """Reject incompatible continuation before replacing configuration evidence."""
    if resume_from is None:
        return
    from morse.optimizer import (
        CHECKPOINT_FILENAME,
        _optimizer_checkpoint_fingerprint,
        _read_optimizer_checkpoint,
        _validate_design_space_fingerprint,
    )

    checkpoint = resume_from / CHECKPOINT_FILENAME
    state_path = resume_from / "run_state.json"
    if checkpoint.is_file() and not state_path.is_file():
        if output_dir.resolve() != resume_from.resolve():
            raise ValueError("An interrupted optimizer checkpoint must resume in its original output directory")
        expected = _optimizer_checkpoint_fingerprint(
            fingerprint, spec.design_space, spec.run_config.optimizer, spec.run_config.outputs
        )
        _read_optimizer_checkpoint(checkpoint, expected)
        return
    if not state_path.is_file():
        raise FileNotFoundError(f"Cannot resume: missing {state_path}")
    state = json.loads(state_path.read_text())
    _validate_design_space_fingerprint(state, spec.design_space, state_path)
    if state.get("scoring_fingerprint") != fingerprint:
        raise ValueError("Cannot resume: scorer configuration differs from the source run")
    if spec.run_config.optimizer.initial_population is not None:
        raise ValueError("resume_from and optimizer.initial_population cannot be used together")


def _prepare_configured_initial_population(
    spec: RunSpec,
    config_snapshot_path: Path,
    output_dir: Path,
) -> RunSpec:
    initialization = spec.run_config.initialization
    optimizer = spec.run_config.optimizer
    if optimizer.initial_population is not None:
        raise ValueError(
            "Configure either initialization or optimizer.initial_population, not both"
        )
    manifest_path = output_dir / "initial-population.json"
    command = [
        sys.executable,
        "-u",
        "-m",
        "morse_warm_starter.cli",
        "generate",
        "--config",
        str(config_snapshot_path),
        "--output",
        str(manifest_path),
        "--population-size",
        str(optimizer.population_size),
        "--seed",
        str(optimizer.seed),
        "--anchors-only",
    ]
    if initialization.mode == "linear_design":
        command.extend(
            [
                "--lineardesign-profile",
                str(initialization.profile_path),
                "--lineardesign-dir",
                initialization.lineardesign_dir,
                "--accept-lineardesign-license",
            ]
        )
        for lambda_value in initialization.lambdas:
            command.extend(["--lineardesign-lambda", str(lambda_value)])
    elif initialization.mode == "balanced_warm_start":
        objective_pairs = [
            (scorer.name, output.name)
            for scorer in spec.run_config.scorers
            for output in scorer.outputs
            if output.role == "objective"
        ]
        desired = max(
            2,
            int(round(optimizer.population_size * initialization.warm_fraction)),
        )
        generated_slots = max(1, desired - 1)
        base, remainder = divmod(generated_slots, len(objective_pairs))
        for index, (scorer_name, output_name) in enumerate(objective_pairs):
            count = base + (1 if index < remainder else 0)
            if count:
                command.extend(
                    ["--greedy", f"{scorer_name}:{output_name}:{count}"]
                )
        command.extend(["--passes", str(initialization.passes)])
    else:  # pragma: no cover - guarded by InitializationConfig
        return spec
    LOGGER.info("Preparing %s initial population", initialization.mode)
    try:
        completed = subprocess.run(
            command,
            check=True,
            text=True,
            capture_output=True,
        )
        if completed.stdout.strip():
            LOGGER.debug("Warm starter stdout:\n%s", completed.stdout.rstrip())
        if completed.stderr.strip():
            LOGGER.debug("Warm starter stderr:\n%s", completed.stderr.rstrip())
    except FileNotFoundError as exc:
        raise RuntimeError(
            "MORSE Warm Starter is not installed in this environment"
        ) from exc
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or "").strip()
        raise RuntimeError(
            "Warm Starter failed" + (f": {detail}" if detail else "")
        ) from exc
    return spec.with_initial_population(manifest_path)


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
