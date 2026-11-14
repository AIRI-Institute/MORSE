from __future__ import annotations

import copy
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import yaml

from morse.config import RunConfig, expand_env_vars, parse_run_config
from morse.construct import (
    PART_ORDER,
    ConstructDesignSpace,
    construct_design_from_config,
)
from morse.router import resolve_scope_ranges


_OPERATORS = {
    "five_prime_utr": "choose_one",
    "start_codon": "choose_one",
    "cds": "synonymous_codons",
    "stop_codon": "choose_one",
    "three_prime_utr": "choose_one",
    "poly_a": "choose_one",
}


@dataclass(frozen=True)
class PartSourceSpec:
    part: str
    kind: str
    path: Path
    set_name: str | None = None
    max_candidates: int | None = None


@dataclass(frozen=True)
class RunSpec:
    path: Path
    name: str
    design_space: ConstructDesignSpace
    run_config: RunConfig
    part_sources: tuple[PartSourceSpec, ...]
    resolved: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return copy.deepcopy(self.resolved)

    def with_initial_population(self, path: Path) -> RunSpec:
        resolved = copy.deepcopy(self.resolved)
        resolved["optimizer"]["initial_population"] = str(path.resolve())
        return replace(
            self,
            run_config=replace(
                self.run_config,
                optimizer=replace(
                    self.run_config.optimizer,
                    initial_population=str(path.resolve()),
                ),
            ),
            resolved=resolved,
        )


def load_run_spec(path: Path, *, input_root: Path | None = None) -> RunSpec:
    path = Path(path).resolve()
    return parse_run_spec(yaml.safe_load(path.read_text()), path=path, input_root=input_root)


def resolve_run_mapping(
    raw: dict[str, Any], *, path: Path, input_root: Path | None = None
) -> dict[str, Any]:
    """Copy and resolve config paths without opening its sequence or libraries."""
    path = Path(path).resolve()
    base_dir = (
        input_root.expanduser().resolve()
        if input_root is not None
        else path.parent
    )
    raw = expand_env_vars(copy.deepcopy(raw))
    if not isinstance(raw, dict):
        raise ValueError(f"Run config must be a YAML mapping: {path}")
    schema_version = int(raw.get("schema_version", 0))
    if schema_version not in {3, 4, 5}:
        if schema_version == 2:
            raise ValueError(
                "Run config schema_version 2 is no longer accepted: migrate "
                "start_codon from a CDS annotation to a physical construct part"
            )
        raise ValueError("Run config must define schema_version: 3, 4 or 5")
    unknown_top_level = sorted(
        set(raw)
        - {
            "schema_version",
            "name",
            "construct",
            "optimizer",
            "scorers",
            "manufacturability",
            "initialization",
            *( {"output_dir"} if schema_version in {4, 5} else set()),
        }
    )
    if unknown_top_level:
        raise ValueError(
            "Run config contains unknown fields: " + ", ".join(unknown_top_level)
        )

    resolved = copy.deepcopy(raw)
    resolved["schema_version"] = schema_version
    name = str(resolved.get("name") or path.stem)
    resolved["name"] = name
    if schema_version in {4, 5}:
        if schema_version == 5:
            from morse.schema5 import resolve_paths
        else:
            from morse.schema4 import resolve_paths
        resolve_paths(resolved, base_dir, output_base_dir=path.parent)
    else:
        _resolve_part_paths(resolved, base_dir)
    _resolve_initial_population_path(resolved, base_dir)
    _resolve_initialization_paths(resolved, base_dir)
    _resolve_scorer_paths(resolved, base_dir)
    return resolved


def parse_run_spec(
    raw: dict[str, Any],
    *,
    path: Path | None = None,
    input_root: Path | None = None,
    sequence: str | None = None,
    require_mutable: bool = True,
) -> RunSpec:
    """Build the CLI's run specification from a mapping and optional input RNA."""
    path = Path(path).resolve() if path is not None else Path.cwd() / "experiment.yaml"
    if sequence is not None:
        if not isinstance(sequence, str):
            raise TypeError("sequence must be a string")
        raw = copy.deepcopy(raw)
        if not isinstance(raw, dict):
            raise ValueError(f"Run config must be a YAML mapping: {path}")
        construct = _mapping(raw.get("construct"), "construct")
        if "coordinates" not in construct or "parts" in construct:
            raise ValueError(
                "An explicit sequence requires construct.coordinates; "
                "construct.parts cannot be used with a sequence override"
            )
        # This path describes the input for snapshots. The parser neither reads
        # nor writes it; the sequence is passed directly to the design builder.
        construct["fasta"] = str(path.parent / "input.fasta")
    resolved = resolve_run_mapping(raw, path=path, input_root=input_root)
    schema_version = resolved["schema_version"]
    name = resolved["name"]
    base_dir = input_root.expanduser().resolve() if input_root is not None else path.parent
    part_sources = _parse_part_sources(resolved, path)

    if schema_version in {4, 5}:
        if schema_version == 5:
            from morse.schema5 import build_design, engine_config
        else:
            from morse.schema4 import build_design, engine_config
        design_space = build_design(name, _mapping(resolved.get("construct"), "construct"),
                                    base_dir, source_sequence=sequence)
        runtime_raw = engine_config(resolved)
    else:
        optimizer = _mapping(resolved.get("optimizer"), "optimizer")
        mutable_parts = (
            frozenset()
            if not require_mutable and optimizer.get("variables") in (None, {})
            else _parse_variables(optimizer.get("variables"))
        )
        design_space = construct_design_from_config(
            name, _mapping(resolved.get("construct"), "construct"),
            base_dir, mutable_parts=mutable_parts, source_sequence=sequence,
        )
        runtime_raw = resolved
    if require_mutable and design_space.mutable_variable_count == 0:
        key = "construct.changes / codon choices" if schema_version == 5 else (
            "construct.changes" if schema_version == 4 else "optimizer.variables"
        )
        raise ValueError(f"{key} does not expose any mutable choices")

    run_config = parse_run_config(runtime_raw)
    _validate_scopes(design_space, run_config)
    return RunSpec(
        path=path,
        name=name,
        design_space=design_space,
        run_config=run_config,
        part_sources=part_sources,
        resolved=resolved,
    )


def write_run_spec(path: Path, spec: RunSpec) -> None:
    path.write_text(yaml.safe_dump(spec.to_dict(), sort_keys=False))


def _parse_variables(raw: Any) -> frozenset[str]:
    variables = _mapping(raw, "optimizer.variables")
    if not variables:
        raise ValueError("optimizer.variables cannot be empty")
    unknown = set(variables) - set(PART_ORDER)
    if unknown:
        raise ValueError(
            "optimizer.variables contains unknown parts: "
            + ", ".join(sorted(unknown))
        )
    for part, raw_spec in variables.items():
        spec = raw_spec if isinstance(raw_spec, dict) else {"operator": raw_spec}
        unknown_fields = sorted(set(spec) - {"operator"})
        if unknown_fields:
            raise ValueError(
                f"optimizer.variables.{part} contains unknown fields: "
                + ", ".join(unknown_fields)
            )
        operator = str(spec.get("operator", ""))
        if operator != _OPERATORS[part]:
            raise ValueError(
                f"optimizer.variables.{part}.operator must be {_OPERATORS[part]}"
            )
    return frozenset(variables)


def _validate_scopes(design: ConstructDesignSpace, run_config: RunConfig) -> None:
    source = design.source_construct
    shortest = replace(
        source,
        stop_codon=(min(design.stop_codons, key=len)
                    if "stop_codon" in design.mutable_parts else source.stop_codon),
        five_prime_utr=(
            min(design.five_prime_utrs, key=lambda item: len(item.sequence))
            if "five_prime_utr" in design.mutable_parts
            else source.five_prime_utr
        ),
        three_prime_utr=(
            min(design.three_prime_utrs, key=lambda item: len(item.sequence))
            if "three_prime_utr" in design.mutable_parts
            else source.three_prime_utr
        ),
        poly_a_length=(
            min(design.poly_a_lengths)
            if "poly_a" in design.mutable_parts
            else source.poly_a_length
        ),
    )
    for scorer in run_config.scorers:
        config = scorer.scope
        try:
            resolve_scope_ranges(source, config)
            resolve_scope_ranges(shortest, config)
        except ValueError as exc:
            raise ValueError(
                f"scorer {scorer.name!r} has an invalid scope: {exc}"
            ) from exc


def _resolve_part_paths(raw: dict[str, Any], base_dir: Path) -> None:
    construct = _mapping(raw.get("construct"), "construct")
    if "fasta" in construct or "coordinates" in construct:
        value = construct.get("fasta")
        if not isinstance(value, str) or not value.strip():
            raise ValueError("construct.fasta must be a non-empty path")
        path = Path(value).expanduser()
        construct["fasta"] = str(
            path if path.is_absolute() else (base_dir / path).resolve()
        )
        return
    parts = _mapping(construct.get("parts"), "construct.parts")
    for part, raw_part in parts.items():
        source = _part_source(raw_part, str(part))
        if str(source.get("type", "")).lower() not in {"fasta", "tsv"}:
            continue
        if not source.get("path"):
            continue
        path = Path(str(source["path"]))
        source["path"] = str(path if path.is_absolute() else (base_dir / path).resolve())


def _resolve_initial_population_path(raw: dict[str, Any], base_dir: Path) -> None:
    optimizer = _mapping(raw.get("optimizer"), "optimizer")
    for key in ("initial_population", "initial_population_audit"):
        value = optimizer.get(key)
        if value is None:
            continue
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"optimizer.{key} must be a non-empty path")
        path = Path(value)
        optimizer[key] = str(path if path.is_absolute() else (base_dir / path).resolve())


def _resolve_initialization_paths(raw: dict[str, Any], base_dir: Path) -> None:
    initialization = raw.get("initialization")
    if initialization is None:
        return
    value = _mapping(initialization, "initialization")
    for field in ("profile_path", "lineardesign_dir"):
        raw_path = value.get(field)
        if raw_path is None:
            continue
        if not isinstance(raw_path, str) or not raw_path.strip():
            raise ValueError(f"initialization.{field} must be a non-empty path")
        path = Path(raw_path).expanduser()
        value[field] = str(path if path.is_absolute() else (base_dir / path).resolve())


def _parse_part_sources(
    raw: dict[str, Any], config_path: Path
) -> tuple[PartSourceSpec, ...]:
    construct = _mapping(raw.get("construct"), "construct")
    if raw.get("schema_version") in {4, 5}:
        changes = construct.get("changes", {})
        sources = []
        for part in PART_ORDER:
            rule = changes.get(part, {})
            rule = rule if isinstance(rule, dict) else {}
            library = rule.get("library")
            sources.append(PartSourceSpec(
                part=part,
                kind="tsv" if str(library).endswith(".tsv") else "fasta",
                path=Path(library or construct["fasta"]),
                set_name=rule.get("set"),
                max_candidates=rule.get("max_candidates"),
            ))
        return tuple(sources)
    if "fasta" in construct or "coordinates" in construct:
        return tuple(
            PartSourceSpec(part=part, kind="fasta", path=Path(construct["fasta"]))
            for part in PART_ORDER
        )
    parts = _mapping(construct.get("parts"), "construct.parts")
    sources: list[PartSourceSpec] = []
    for part in PART_ORDER:
        source = _part_source(parts[part], part)
        kind = str(source.get("type", "")).lower()
        source_path = (
            Path(str(source["path"]))
            if kind in {"fasta", "tsv"}
            else config_path
        )
        max_candidates = source.get("max_candidates")
        sources.append(
            PartSourceSpec(
                part=part,
                kind=kind,
                path=source_path,
                set_name=(
                    str(source["set"]) if source.get("set") is not None else None
                ),
                max_candidates=(
                    int(max_candidates) if max_candidates is not None else None
                ),
            )
        )
    return tuple(sources)


def _part_source(raw_part: Any, part: str) -> dict[str, Any]:
    part_mapping = _mapping(raw_part, f"construct.parts.{part}")
    return _mapping(part_mapping.get("source"), f"construct.parts.{part}.source")


def _mapping(value: Any, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a mapping")
    return value


def _resolve_scorer_paths(raw: dict[str, Any], base_dir: Path) -> None:
    """Resolve explicit scorer paths against the experiment, never shell cwd."""
    for scorer in raw.get("scorers", []):
        if not isinstance(scorer, dict):
            continue
        if isinstance(scorer.get("command"), str) and scorer.get("workdir") is None:
            scorer["workdir"] = str(base_dir.resolve())
        if scorer.get("workdir") is not None:
            value = Path(str(scorer["workdir"])).expanduser()
            scorer["workdir"] = str(value if value.is_absolute() else (base_dir / value).resolve())
        if isinstance(scorer.get("provenance_paths"), list):
            scorer["provenance_paths"] = [
                str(p if p.is_absolute() else (base_dir / p).resolve())
                for p in (Path(str(value)).expanduser() for value in scorer["provenance_paths"])
            ]
