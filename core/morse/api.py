"""Python entry points for the same configurations and engine as the CLI."""
from __future__ import annotations

import copy
from collections.abc import Callable, Mapping, Sequence
from dataclasses import replace
from datetime import datetime
from pathlib import Path
import re
from typing import Any

from morse.results import OptimizationResult


def load_config(path: str | Path) -> dict[str, Any]:
    """Read YAML into an ordinary dict, resolving paths against its directory."""
    import yaml
    from morse.run_spec import resolve_run_mapping

    path = Path(path).expanduser().resolve()
    return resolve_run_mapping(yaml.safe_load(path.read_text()), path=path)


def _prepare(sequence, params, base_dir, *, require_mutable=False):
    from morse.run_spec import parse_run_spec

    if not isinstance(params, Mapping):
        raise TypeError("params must be a configuration dict; use morse.load_config(path)")
    if sequence is not None and not isinstance(sequence, str):
        raise TypeError("rna_seq must be an RNA/DNA string or None")
    base = Path(base_dir).expanduser().resolve() if base_dir is not None else Path.cwd()
    return parse_run_spec(
        copy.deepcopy(dict(params)), path=base / "experiment.yaml", input_root=base,
        sequence=sequence, require_mutable=require_mutable,
    )


def validate(rna_seq: str | None = None, params: Mapping[str, Any] | None = None,
             *, base_dir: str | Path | None = None) -> dict[str, Any]:
    """Check configuration, coordinates and structure without calling scorers."""
    from morse.validators import StructuralValidator

    spec = _prepare(rna_seq, params, base_dir)
    construct = spec.design_space.source_construct
    validation = StructuralValidator().validate(construct, spec.design_space)
    return {
        "valid": validation.valid, "name": spec.name,
        "sequence_length": len(construct.sequence),
        "coordinates": construct.to_dict()["coordinates"],
        "mutable_parts": sorted(spec.design_space.mutable_parts),
        "structural_validation": validation.to_dict(),
    }


def score(rna_seq: str | None = None, params: Mapping[str, Any] | None = None,
          *, base_dir: str | Path | None = None) -> dict[str, Any]:
    """Score the supplied sequence, or the config's FASTA, without optimizing."""
    from morse.execution import score_reference

    return score_reference(_prepare(rna_seq, params, base_dir))


def _observer(output_dir, verbose, plots, on_progress):
    from morse.monitor import ProgressDisplay

    if on_progress is not None and not callable(on_progress):
        raise TypeError("on_progress must be callable or None")
    display = ProgressDisplay(output_dir, verbose=verbose, plots=plots)

    def update(payload):
        display.update(payload)
        if on_progress is not None:
            on_progress(copy.deepcopy(payload))

    return display, update


def optimize(
    rna_seq: str | None = None,
    params: Mapping[str, Any] | None = None,
    *,
    output_dir: str | Path | None = None,
    verbose: int = 1,
    plots: Sequence[str] | None = None,
    base_dir: str | Path | None = None,
    on_progress: Callable[[dict[str, Any]], None] | None = None,
) -> OptimizationResult:
    """Run optimization synchronously, persisting inputs, progress and results.

    ``rna_seq`` overrides ``construct.fasta`` and uses the configured coordinates.
    ``verbose`` is 0 (quiet), 1 (status), or 2 (status and live plots).
    Plain dict paths are relative to ``base_dir`` or the working directory.
    """
    from morse.execution import execute_run
    from morse.fasta import FastaRecord, write_fasta
    from morse.run_spec import write_run_spec

    spec = _prepare(rna_seq, params, base_dir, require_mutable=True)
    selected = output_dir or spec.resolved.get("output_dir")
    if selected is None:
        name = re.sub(r"[^A-Za-z0-9_-]", "_", spec.name) or "run"
        selected = Path("outputs") / f"{name}_{datetime.now():%Y%m%d_%H%M%S_%f}"
    target = Path(selected).expanduser().resolve()
    display, observer = _observer(target, verbose, plots, on_progress)
    if target.exists() and (not target.is_dir() or any(target.iterdir())):
        raise FileExistsError(f"Result directory is not empty: {target}; use a new directory or morse.resume()")
    resolved = spec.to_dict()
    if resolved["schema_version"] == 4:
        resolved["output_dir"] = str(target)
    if rna_seq is not None:
        previous = Path(resolved["construct"]["fasta"])
        fasta = target / "input.fasta"
        resolved["construct"]["fasta"] = str(fasta)
        sources = tuple(replace(source, path=fasta) if source.path == previous else source
                        for source in spec.part_sources)
        spec = replace(spec, part_sources=sources)
    spec = replace(spec, path=target / "experiment.yaml", resolved=resolved)
    target.mkdir(parents=True, exist_ok=True)
    if rna_seq is not None:
        write_fasta(fasta, [FastaRecord(spec.name, spec.design_space.source_construct.sequence.replace("T", "U"))])
    write_run_spec(spec.path, spec)
    display.start()
    try:
        execute_run(spec, target, spec.name, None, on_progress=observer)
    except BaseException:
        _stopped(display, target)
        raise
    display.finish("Completed")
    return OptimizationResult(target)


def _stopped(display, target):
    if (target / "optimizer-checkpoint.pkl").is_file():
        message = "Stopped; continue from the last saved generation with morse.resume()"
    else:
        message = "Stopped before the first generation was saved"
    display.finish(message, success=False)


def load_run(path: str | Path) -> OptimizationResult:
    """Open a completed or partial run without rerunning scorers."""
    return OptimizationResult(path)


def resume(
    path: str | Path, *, verbose: int = 1, plots: Sequence[str] | None = None,
    on_progress: Callable[[dict[str, Any]], None] | None = None,
) -> OptimizationResult:
    """Continue an interrupted run in place from its last complete generation."""
    from morse.execution import execute_run
    from morse.run_spec import load_run_spec

    target = Path(path).expanduser().resolve()
    result = OptimizationResult(target)
    display, observer = _observer(target, verbose, plots, on_progress)
    if (target / "run_state.json").is_file():
        # A finished run is already a result, not a request for a new experiment.
        if not (target / "report.md").is_file():
            from morse.execution import write_run_report
            config = target / "resolved_run.yaml"
            spec = load_run_spec(config)
            write_run_report(target, spec.run_config, config)
        display.finish("Already completed")
        return result
    if not (target / "optimizer-checkpoint.pkl").is_file():
        raise FileNotFoundError(f"No saved generation to resume in {target}")
    config = target / "resolved_run.yaml"
    if not config.is_file():
        raise FileNotFoundError(f"Cannot resume without {config}")
    spec = load_run_spec(config)
    display.start()
    try:
        execute_run(spec, target, spec.name, target, on_progress=observer)
    except BaseException:
        _stopped(display, target)
        raise
    display.finish("Completed")
    return OptimizationResult(target)
