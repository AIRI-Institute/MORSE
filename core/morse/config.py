from __future__ import annotations

import copy
from dataclasses import dataclass, field, replace
from math import cos, isfinite, pi
import os
import re
from typing import Any

from morse.scorer_scope import (
    ScorerScopeConfig,
    parse_scorer_scope,
    scorer_scope_to_dict,
)
from morse.validators import ManufacturabilityPolicy

_ENV_VAR_PATTERN = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)(?::-(.*?))?\}")
_CONTAINER_SCORER_MODES = {"persistent-container"}
_PROCESS_SCORER_MODES = {"persistent-process"}
_BASH_SCORER_MODE = "bash-command"
_SCORER_MODES = _CONTAINER_SCORER_MODES | _PROCESS_SCORER_MODES | {_BASH_SCORER_MODE}
_OUTPUT_PREFIX_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,63}$")


@dataclass(frozen=True)
class ScorerConfig:
    name: str
    mode: str
    image: str | None
    timeout_seconds: int
    scope: ScorerScopeConfig
    outputs: tuple[OutputSpec, ...] = ()
    command: str | list[str] | None = None
    env: dict[str, str] | None = None
    workdir: str | None = None
    docker_args: list[str] | None = None
    output_prefix: str | None = None
    provenance_paths: tuple[str, ...] = ()
    namespace: str | None = None
    parameters: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "parameters", parse_scorer_parameters(self.parameters))


@dataclass(frozen=True)
class EarlyStoppingConfig:
    enabled: bool = False
    metric: str = "hypervolume"
    min_generations: int = 0
    patience: int = 0
    rel_tol: float = 0.0
    abs_tol: float = 0.0


@dataclass(frozen=True)
class ConstraintHandlingConfig:
    method: str = "feasibility_first"
    normalization: str = "none"
    initial_epsilon: float = 0.0
    relaxation_generations: int = 0
    curve: str = "linear"

    def __post_init__(self) -> None:
        if self.method not in {"feasibility_first", "adaptive_epsilon"}:
            raise ValueError(
                "optimizer.constraint_handling.method must be "
                "feasibility_first or adaptive_epsilon"
            )
        if self.normalization not in {"none", "relative"}:
            raise ValueError(
                "optimizer.constraint_handling.normalization must be none or relative"
            )
        if self.curve not in {"linear", "cosine"}:
            raise ValueError(
                "optimizer.constraint_handling.curve must be linear or cosine"
            )
        if self.initial_epsilon < 0:
            raise ValueError(
                "optimizer.constraint_handling.initial_epsilon must be >= 0"
            )
        if self.relaxation_generations < 0:
            raise ValueError(
                "optimizer.constraint_handling.relaxation_generations must be >= 0"
            )
        if self.method == "feasibility_first" and (
            self.initial_epsilon != 0.0 or self.relaxation_generations != 0
        ):
            raise ValueError(
                "feasibility_first constraint handling does not accept "
                "initial_epsilon or relaxation_generations"
            )
        if self.method == "adaptive_epsilon":
            if self.normalization != "relative":
                raise ValueError(
                    "adaptive_epsilon constraint handling requires "
                    "normalization: relative"
                )
            if self.initial_epsilon <= 0:
                raise ValueError(
                    "adaptive_epsilon constraint handling requires initial_epsilon > 0"
                )
            if self.relaxation_generations < 2:
                raise ValueError(
                    "adaptive_epsilon constraint handling requires "
                    "relaxation_generations >= 2"
                )

    def epsilon(self, generation: int) -> float:
        if self.method != "adaptive_epsilon":
            return 0.0
        if generation <= 1:
            return self.initial_epsilon
        if generation >= self.relaxation_generations:
            return 0.0
        progress = (generation - 1) / (self.relaxation_generations - 1)
        if self.curve == "cosine":
            fraction = 0.5 * (1.0 + cos(pi * progress))
        else:
            fraction = 1.0 - progress
        return self.initial_epsilon * fraction


@dataclass(frozen=True)
class ReferenceDirectionsConfig:
    method: str
    n_partitions: int | None = None
    n_points: int | None = None
    scaling: float | None = None
    seed: int | None = None

    def __post_init__(self) -> None:
        if self.method not in {"das-dennis", "energy"}:
            raise ValueError(
                "optimizer.reference_directions.method must be das-dennis or energy"
            )
        if self.method == "das-dennis":
            if self.n_partitions is None or self.n_partitions < 1:
                raise ValueError(
                    "das-dennis reference directions require n_partitions >= 1"
                )
            if self.n_points is not None or self.seed is not None:
                raise ValueError(
                    "das-dennis reference directions do not accept n_points or seed"
                )
            if self.scaling is not None and not 0.0 < self.scaling <= 1.0:
                raise ValueError(
                    "reference_directions.scaling must satisfy 0 < scaling <= 1"
                )
        else:
            if self.n_points is None or self.n_points < 2:
                raise ValueError("energy reference directions require n_points >= 2")
            if self.n_partitions is not None or self.scaling is not None:
                raise ValueError(
                    "energy reference directions do not accept n_partitions or scaling"
                )


@dataclass(frozen=True)
class OptimizerConfig:
    population_size: int
    generations: int
    seed: int
    type: str = "nsga2"
    n_offsprings: int | None = None
    eliminate_duplicates: bool = True
    crossover_type: str = "uniform"
    crossover_probability: float = 0.9
    mutation_type: str = "integer_reset"
    mutation_probability: float | None = None
    reference_directions: ReferenceDirectionsConfig | None = None
    initial_population: str | None = None
    initial_population_policy: str = "filter"
    initial_population_audit: str | None = None
    constraint_handling: ConstraintHandlingConfig = ConstraintHandlingConfig()
    early_stopping: EarlyStoppingConfig = EarlyStoppingConfig()

    def __post_init__(self) -> None:
        if self.type not in {"ga", "nsga2", "nsga3"}:
            raise ValueError("optimizer.type must be ga, nsga2, or nsga3")
        if self.population_size < 2:
            raise ValueError("optimizer.population_size must be >= 2")
        if self.generations < 1:
            raise ValueError("optimizer.generations must be >= 1")
        if self.n_offsprings is not None and self.n_offsprings < 1:
            raise ValueError("optimizer.n_offsprings must be >= 1")
        if self.crossover_type != "uniform":
            raise ValueError("optimizer.crossover.type currently supports only uniform")
        if not 0.0 <= self.crossover_probability <= 1.0:
            raise ValueError(
                "optimizer.crossover.probability must be between 0 and 1"
            )
        if self.mutation_type != "integer_reset":
            raise ValueError(
                "optimizer.mutation.type currently supports only integer_reset"
            )
        if self.mutation_probability is not None and not (
            0.0 <= self.mutation_probability <= 1.0
        ):
            raise ValueError(
                "optimizer.mutation.probability must be auto or between 0 and 1"
            )
        if self.type == "nsga3" and self.reference_directions is None:
            raise ValueError(
                "optimizer.reference_directions is required for optimizer.type nsga3"
            )
        if self.type != "nsga3" and self.reference_directions is not None:
            raise ValueError(
                "optimizer.reference_directions is supported only for optimizer.type nsga3"
            )
        if self.initial_population is not None and not self.initial_population.strip():
            raise ValueError("optimizer.initial_population must be a non-empty path")
        if self.initial_population_policy not in {"filter", "retain_all"}:
            raise ValueError("optimizer.initial_population_policy must be filter or retain_all")
        if self.initial_population_policy == "retain_all":
            if self.initial_population is None:
                raise ValueError("retain_all requires optimizer.initial_population")
            if self.constraint_handling.method != "feasibility_first":
                raise ValueError("retain_all requires strict feasibility_first constraints")
        if self.initial_population_audit is not None and not self.initial_population_audit.strip():
            raise ValueError("optimizer.initial_population_audit must be a non-empty path")
        if self.initial_population_audit and self.constraint_handling.method != "feasibility_first":
            raise ValueError("initial_population_audit requires strict feasibility_first constraints")
        if self.constraint_handling.relaxation_generations > self.generations:
            raise ValueError(
                "optimizer.constraint_handling.relaxation_generations must be "
                "<= optimizer.generations"
            )
        if (
            self.early_stopping.enabled
            and self.early_stopping.metric == "best_feasible_objective"
            and self.type != "ga"
        ):
            raise ValueError(
                "best_feasible_objective early stopping is supported only for optimizer.type ga"
            )
        if (
            self.early_stopping.enabled
            and self.constraint_handling.method == "adaptive_epsilon"
            and self.early_stopping.min_generations
            < self.constraint_handling.relaxation_generations
        ):
            raise ValueError(
                "optimizer.early_stopping.min_generations must be >= "
                "constraint_handling.relaxation_generations"
            )


@dataclass(frozen=True)
class InitializationConfig:
    mode: str = "random"
    warm_fraction: float = 0.25
    passes: int = 1
    profile_path: str | None = None
    lambdas: tuple[float, ...] = ()
    lineardesign_dir: str = "~/.cache/morse/LinearDesign"
    accept_lineardesign_license: bool = False

    def __post_init__(self) -> None:
        if self.mode not in {"random", "balanced_warm_start", "linear_design"}:
            raise ValueError(
                "initialization.mode must be random, balanced_warm_start, or linear_design"
            )
        if not 0.1 <= self.warm_fraction <= 0.75:
            raise ValueError("initialization.warm_fraction must be between 0.1 and 0.75")
        if not 1 <= self.passes <= 5:
            raise ValueError("initialization.passes must be between 1 and 5")
        if self.mode == "linear_design":
            if not self.profile_path:
                raise ValueError("linear_design initialization requires profile_path")
            if not 1 <= len(self.lambdas) <= 32:
                raise ValueError("linear_design initialization requires 1 to 32 lambdas")
            if any(not 0.0 <= value <= 100.0 for value in self.lambdas):
                raise ValueError("linear_design lambdas must be between 0 and 100")
            if not self.lineardesign_dir:
                raise ValueError("linear_design initialization requires lineardesign_dir")
            if not self.accept_lineardesign_license:
                raise ValueError(
                    "linear_design initialization requires explicit upstream license acceptance"
                )

    def to_dict(self) -> dict[str, Any]:
        if self.mode == "random":
            return {"mode": "random"}
        if self.mode == "balanced_warm_start":
            return {
                "mode": self.mode,
                "warm_fraction": self.warm_fraction,
                "passes": self.passes,
            }
        return {
            "mode": self.mode,
            "profile_path": self.profile_path,
            "lambdas": list(self.lambdas),
            "lineardesign_dir": self.lineardesign_dir,
            "accept_lineardesign_license": self.accept_lineardesign_license,
        }


@dataclass(frozen=True)
class ThresholdResult:
    output: str
    valid: bool
    observed: Any | None
    min_value: float | None
    max_value: float | None
    equals_value: Any | None
    reason: str


@dataclass(frozen=True)
class OutputSpec:
    name: str
    type: str
    role: str
    direction: str | None = None
    min_value: float | None = None
    max_value: float | None = None
    equals_value: Any | None = None
    description: str = ""

    def coerce(self, value: Any) -> Any:
        if self.type in {"float", "number"}:
            try:
                number = float(value)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ValueError(f"Metric {self.name!r} must be a finite number") from exc
            if not isfinite(number):
                raise ValueError(f"Metric {self.name!r} must be finite; got {number}")
            return number
        if self.type in {"int", "integer"}:
            try:
                return int(value)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ValueError(f"Metric {self.name!r} must be a finite integer") from exc
        if self.type in {"str", "string", "structure"}:
            return str(value)
        if self.type in {"bool", "boolean", "binary"}:
            return _parse_bool(value)
        raise ValueError(f"Unsupported output type for {self.name}: {self.type}")

    def objective_value(self, outputs: dict[str, Any]) -> float:
        value = float(self.coerce(_required_output(outputs, self.name)))
        if self.direction == "min":
            return value
        if self.direction == "max":
            return -value
        raise ValueError(f"Objective output {self.name} must define direction min or max")

    def check_threshold(self, outputs: dict[str, Any]) -> ThresholdResult:
        observed = self.coerce(_required_output(outputs, self.name))
        reasons: list[str] = []
        if self.min_value is not None and float(observed) < self.min_value:
            reasons.append(
                f"{self.name}={_format_value(observed)} below min {self.min_value:.6g}"
            )
        if self.max_value is not None and float(observed) > self.max_value:
            reasons.append(
                f"{self.name}={_format_value(observed)} above max {self.max_value:.6g}"
            )
        expected = (
            None if self.equals_value is None else self.coerce(self.equals_value)
        )
        if expected is not None and observed != expected:
            reasons.append(
                f"{self.name}={_format_value(observed)} does not equal "
                f"{_format_value(expected)}"
            )
        return ThresholdResult(
            output=self.name,
            valid=not reasons,
            observed=observed,
            min_value=self.min_value,
            max_value=self.max_value,
            equals_value=expected,
            reason="; ".join(reasons),
        )


@dataclass(frozen=True)
class RunConfig:
    scorers: list[ScorerConfig]
    optimizer: OptimizerConfig
    outputs: list[OutputSpec]
    manufacturability: ManufacturabilityPolicy = ManufacturabilityPolicy()
    initialization: InitializationConfig = InitializationConfig()

    def to_dict(self) -> dict[str, Any]:
        scorer_dicts = [
            {
                key: value
                for key, value in {
                    "name": scorer.name,
                    "mode": None if scorer.mode == _BASH_SCORER_MODE else scorer.mode,
                    "image": scorer.image,
                    "timeout_seconds": scorer.timeout_seconds,
                    "command": scorer.command,
                    "env": scorer.env,
                    "workdir": scorer.workdir,
                    "docker_args": scorer.docker_args,
                    "output_prefix": scorer.output_prefix,
                    "namespace": scorer.namespace,
                    "provenance_paths": list(scorer.provenance_paths) or None,
                    "parameters": copy.deepcopy(scorer.parameters) or None,
                    "scope": scorer_scope_to_dict(scorer.scope),
                    "metrics": _metrics_to_dict(scorer.outputs) if scorer.outputs else None,
                }.items()
                if value is not None
            }
            for scorer in self.scorers
        ]
        return {
            "scorers": scorer_dicts,
            "manufacturability": self.manufacturability.to_dict(),
            "initialization": self.initialization.to_dict(),
            "optimizer": {
                "type": self.optimizer.type,
                "population_size": self.optimizer.population_size,
                "generations": self.optimizer.generations,
                "seed": self.optimizer.seed,
                "n_offsprings": self.optimizer.n_offsprings,
                "eliminate_duplicates": self.optimizer.eliminate_duplicates,
                "crossover": {
                    "type": self.optimizer.crossover_type,
                    "probability": self.optimizer.crossover_probability,
                },
                "mutation": {
                    "type": self.optimizer.mutation_type,
                    "probability": (
                        "auto"
                        if self.optimizer.mutation_probability is None
                        else self.optimizer.mutation_probability
                    ),
                },
                "reference_directions": _reference_directions_to_dict(
                    self.optimizer.reference_directions
                ),
                "initial_population": self.optimizer.initial_population,
                "initial_population_policy": self.optimizer.initial_population_policy,
                "initial_population_audit": self.optimizer.initial_population_audit,
                "constraint_handling": _constraint_handling_to_dict(
                    self.optimizer.constraint_handling
                ),
                "early_stopping": {
                    "enabled": self.optimizer.early_stopping.enabled,
                    "metric": self.optimizer.early_stopping.metric,
                    "min_generations": self.optimizer.early_stopping.min_generations,
                    "patience": self.optimizer.early_stopping.patience,
                    "rel_tol": self.optimizer.early_stopping.rel_tol,
                    "abs_tol": self.optimizer.early_stopping.abs_tol,
                },
            },
        }


def parse_run_config(raw: dict[str, Any]) -> RunConfig:
    scorers = _parse_scorers(raw)
    optimizer_raw = _mapping(raw.get("optimizer"), "optimizer")
    _reject_unknown_fields(
        optimizer_raw,
        {
            "type",
            "population_size",
            "generations",
            "seed",
            "n_offsprings",
            "eliminate_duplicates",
            "crossover",
            "mutation",
            "reference_directions",
            "initial_population",
            "initial_population_policy",
            "initial_population_audit",
            "constraint_handling",
            "early_stopping",
            "variables",
        },
        "optimizer",
    )
    crossover_type, crossover_probability = _parse_crossover(optimizer_raw.get("crossover"))
    mutation_type, mutation_probability = _parse_mutation(optimizer_raw.get("mutation"))

    outputs = _flatten_scorer_outputs(scorers)
    if not outputs:
        raise ValueError("Config must define at least one metric under `scorers[].metrics`")
    if not any(spec.role == "objective" for spec in outputs):
        raise ValueError("Config must define at least one metric with role: optimize")

    optimizer = OptimizerConfig(
        population_size=int(optimizer_raw.get("population_size", 32)),
        generations=int(optimizer_raw.get("generations", 12)),
        seed=int(optimizer_raw.get("seed", 7)),
        type=str(optimizer_raw.get("type", "nsga2")).lower(),
        n_offsprings=(
            None
            if optimizer_raw.get("n_offsprings") is None
            else int(optimizer_raw["n_offsprings"])
        ),
        eliminate_duplicates=_parse_bool(
            optimizer_raw.get("eliminate_duplicates", True)
        ),
        crossover_type=crossover_type,
        crossover_probability=crossover_probability,
        mutation_type=mutation_type,
        mutation_probability=mutation_probability,
        reference_directions=_parse_reference_directions(
            optimizer_raw.get("reference_directions")
        ),
        initial_population=(
            None
            if optimizer_raw.get("initial_population") is None
            else str(optimizer_raw["initial_population"])
        ),
        initial_population_policy=str(optimizer_raw.get("initial_population_policy", "filter")),
        initial_population_audit=(None if optimizer_raw.get("initial_population_audit") is None else str(optimizer_raw["initial_population_audit"])),
        constraint_handling=_parse_constraint_handling(
            optimizer_raw.get("constraint_handling")
        ),
        early_stopping=_parse_early_stopping(optimizer_raw.get("early_stopping")),
    )
    objective_count = sum(spec.role == "objective" for spec in outputs)
    if optimizer.type == "ga" and objective_count != 1:
        raise ValueError(
            "optimizer.type ga requires exactly one output with role: objective"
        )
    if (
        optimizer.constraint_handling.method == "adaptive_epsilon"
        and not any(spec.role == "threshold" for spec in outputs)
    ):
        raise ValueError(
            "adaptive_epsilon constraint handling requires at least one "
            "output with role: threshold"
        )
    return RunConfig(
        scorers=scorers,
        optimizer=optimizer,
        outputs=outputs,
        manufacturability=_parse_manufacturability(
            raw.get("manufacturability")
        ),
        initialization=_parse_initialization(raw.get("initialization")),
    )


def _parse_initialization(raw: Any) -> InitializationConfig:
    if raw is None:
        return InitializationConfig()
    value = _mapping(raw, "initialization")
    allowed = {
        "mode",
        "warm_fraction",
        "passes",
        "profile_path",
        "lambdas",
        "lineardesign_dir",
        "accept_lineardesign_license",
    }
    _reject_unknown_fields(value, allowed, "initialization")
    lambdas_raw = value.get("lambdas", [])
    if not isinstance(lambdas_raw, list):
        raise ValueError("initialization.lambdas must be a list")
    return InitializationConfig(
        mode=str(value.get("mode", "random")).strip().lower(),
        warm_fraction=float(value.get("warm_fraction", 0.25)),
        passes=int(value.get("passes", 1)),
        profile_path=(
            None
            if value.get("profile_path") is None
            else str(value["profile_path"])
        ),
        lambdas=tuple(float(item) for item in lambdas_raw),
        lineardesign_dir=str(value.get("lineardesign_dir", "~/.cache/morse/LinearDesign")),
        accept_lineardesign_license=_parse_bool(
            value.get("accept_lineardesign_license", False)
        ),
    )


def _parse_manufacturability(raw: Any) -> ManufacturabilityPolicy:
    if raw is None:
        return ManufacturabilityPolicy()
    value = _mapping(raw, "manufacturability")
    allowed = {
        "enabled",
        "forbidden_sequences",
        "include_reverse_complements",
        "gc_min",
        "gc_max",
        "local_gc_window",
        "local_gc_min",
        "local_gc_max",
        "max_homopolymer",
        "max_tandem_repeat_copies",
        "max_tandem_repeat_unit",
        "max_direct_repeat",
        "max_inverted_repeat",
        "palindrome_min_loop",
    }
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise ValueError(
            "Unknown manufacturability fields: " + ", ".join(unknown)
        )
    forbidden = value.get("forbidden_sequences", [])
    if not isinstance(forbidden, list) or not all(
        isinstance(item, str) for item in forbidden
    ):
        raise ValueError("manufacturability.forbidden_sequences must be a list of strings")
    return ManufacturabilityPolicy(
        enabled=_parse_bool(value.get("enabled", True)),
        forbidden_sequences=tuple(forbidden),
        include_reverse_complements=_parse_bool(
            value.get("include_reverse_complements", True)
        ),
        gc_min=float(value.get("gc_min", 0.30)),
        gc_max=float(value.get("gc_max", 0.70)),
        local_gc_window=int(value.get("local_gc_window", 50)),
        local_gc_min=float(value.get("local_gc_min", 0.20)),
        local_gc_max=float(value.get("local_gc_max", 0.80)),
        max_homopolymer=int(value.get("max_homopolymer", 8)),
        max_tandem_repeat_copies=int(
            value.get("max_tandem_repeat_copies", 4)
        ),
        max_tandem_repeat_unit=int(value.get("max_tandem_repeat_unit", 6)),
        max_direct_repeat=int(value.get("max_direct_repeat", 20)),
        max_inverted_repeat=int(value.get("max_inverted_repeat", 12)),
        palindrome_min_loop=int(value.get("palindrome_min_loop", 3)),
    )


def _parse_scorers(raw: dict[str, Any]) -> list[ScorerConfig]:
    scorers_raw = raw.get("scorers")
    if not isinstance(scorers_raw, list) or not scorers_raw:
        raise ValueError("`scorers` must be a non-empty list")
    scorers = [
        _parse_scorer_config(_mapping(item, f"scorers[{index}]"), default_name=f"scorer_{index}")
        for index, item in enumerate(scorers_raw, start=1)
    ]
    names = [scorer.name for scorer in scorers]
    if len(names) != len(set(names)):
        raise ValueError("scorers contains duplicate names")
    return scorers


def _parse_scorer_config(raw: dict[str, Any], default_name: str) -> ScorerConfig:
    allowed_fields = {
        "name",
        "mode",
        "image",
        "timeout_seconds",
        "outputs",
        "metrics",
        "command",
        "env",
        "workdir",
        "docker_args",
        "output_prefix",
        "provenance_paths",
        "scope",
        "namespace",
        "parameters",
        "params",
    }
    unknown = sorted(set(raw) - allowed_fields)
    if unknown:
        raise ValueError("Unknown scorer fields: " + ", ".join(unknown))
    mode = str(raw.get("mode") or ("bash-command" if isinstance(raw.get("command"), str) else ""))
    if mode not in _SCORER_MODES:
        allowed = ", ".join(sorted(_SCORER_MODES))
        raise ValueError(f"scorer mode must be one of: {allowed}")
    image = raw.get("image")
    if mode in _CONTAINER_SCORER_MODES and (image is None or str(image).strip() == ""):
        raise ValueError("container scorer must define image")
    command = _parse_command(raw.get("command"), mode)
    scorer_name = str(raw.get("name", default_name))
    output_prefix_raw = raw.get("output_prefix")
    output_prefix = None
    if output_prefix_raw is not None:
        output_prefix = str(output_prefix_raw).strip()
        if not _OUTPUT_PREFIX_PATTERN.fullmatch(output_prefix):
            raise ValueError(
                "scorer.output_prefix must start with a letter and contain only "
                "letters, numbers, underscore, hyphen, or dot"
            )

    env_raw = raw.get("env")
    env = None
    if env_raw is not None:
        env_map = _mapping(env_raw, "scorer.env")
        env = {str(key): str(value) for key, value in env_map.items()}

    docker_args_raw = raw.get("docker_args")
    docker_args = None
    if docker_args_raw is not None:
        if mode not in _CONTAINER_SCORER_MODES:
            raise ValueError("scorer.docker_args is supported only for container modes")
        if not isinstance(docker_args_raw, list) or not all(
            isinstance(item, str) for item in docker_args_raw
        ):
            raise ValueError("scorer.docker_args must be a list of strings")
        docker_args = list(docker_args_raw)

    provenance_paths = raw.get("provenance_paths", [])
    if not isinstance(provenance_paths, list) or not all(isinstance(p, str) and p for p in provenance_paths):
        raise ValueError("scorer.provenance_paths must be a list of non-empty paths")

    if raw.get("scope") is None:
        raise ValueError(f"scorer {scorer_name!r} must define scope")
    outputs = tuple(_parse_scorer_outputs(raw, scorer_name))
    timeout_seconds = float(raw.get("timeout_seconds", 180)) if mode == _BASH_SCORER_MODE else int(raw.get("timeout_seconds", 180))
    if mode == _BASH_SCORER_MODE and (not isfinite(timeout_seconds) or timeout_seconds <= 0):
        raise ValueError("Bash scorer timeout_seconds must be positive and finite")
    return ScorerConfig(
        name=scorer_name,
        mode=mode,
        image=None if image is None else str(image),
        timeout_seconds=timeout_seconds,
        scope=parse_scorer_scope(
            raw["scope"],
            f"scorers.{scorer_name}.scope",
        ),
        outputs=outputs,
        output_prefix=output_prefix,
        command=command,
        env=env,
        workdir=None if raw.get("workdir") is None else str(raw["workdir"]),
        docker_args=docker_args,
        provenance_paths=tuple(provenance_paths),
        namespace=None if raw.get("namespace") is None else str(raw["namespace"]),
        parameters=scorer_parameters_from_mapping(raw),
    )


def scorer_parameters_from_mapping(raw: Any) -> dict[str, Any]:
    """Read calculation settings, accepting the shorter spelling as an alias."""
    if "params" in raw and "parameters" in raw:
        raise ValueError("Use scorer.parameters or scorer.params, not both")
    return parse_scorer_parameters(raw.get("parameters", raw.get("params", {})))


def parse_scorer_parameters(raw: Any) -> dict[str, Any]:
    """Copy a scorer's finite JSON object without coercing keys or values."""
    if not isinstance(raw, dict):
        raise ValueError("scorer.parameters must be a mapping")
    active: set[int] = set()

    def copy_value(value: Any, location: str) -> Any:
        if value is None or type(value) in {str, bool, int}:
            return value
        if type(value) is float:
            if not isfinite(value):
                raise ValueError(f"{location} must be finite")
            return value
        if isinstance(value, (dict, list)):
            if id(value) in active:
                raise ValueError(f"{location} must not contain circular references")
            active.add(id(value))
            try:
                if isinstance(value, dict):
                    result = {}
                    for key, item in value.items():
                        if not isinstance(key, str) or not key.strip():
                            raise ValueError(f"{location} keys must be non-empty strings")
                        result[key] = copy_value(item, f"{location}.{key}")
                    return result
                return [copy_value(item, f"{location}[{index}]")
                        for index, item in enumerate(value)]
            finally:
                active.remove(id(value))
        raise ValueError(
            f"{location} must contain only JSON values "
            "(objects, arrays, strings, finite numbers, booleans or null)"
        )

    return copy_value(raw, "scorer.parameters")


def _parse_command(raw_command: Any, mode: str) -> str | list[str] | None:
    if mode == _BASH_SCORER_MODE:
        if not isinstance(raw_command, str) or not raw_command.strip():
            raise ValueError("scorer.command must be a non-empty Bash string")
        return raw_command
    if mode not in _PROCESS_SCORER_MODES:
        if raw_command is not None:
            raise ValueError("scorer.command is supported only for process modes")
        return None
    if raw_command is None:
        raise ValueError("process scorer must define command")
    if not isinstance(raw_command, list) or not raw_command:
        raise ValueError("process scorer command must be a non-empty list of strings")
    if not all(isinstance(item, str) and item for item in raw_command):
        raise ValueError("process scorer command must be a non-empty list of strings")
    return list(raw_command)


def _parse_scorer_outputs(raw: dict[str, Any], scorer_name: str) -> list[OutputSpec]:
    if "metrics" in raw and "outputs" in raw:
        raise ValueError(f"scorer {scorer_name!r} cannot define both metrics and legacy outputs")
    if "metrics" in raw:
        metrics = raw["metrics"]
        if not isinstance(metrics, dict) or not metrics:
            raise ValueError(f"scorer {scorer_name!r} must define a non-empty metrics mapping")
        return [_parse_metric_spec(name, spec, context=f"scorers.{scorer_name}.metrics.{name}")
                for name, spec in metrics.items()]
    outputs_raw = raw.get("outputs")
    if outputs_raw is None:
        return []
    if not isinstance(outputs_raw, dict) or not outputs_raw:
        raise ValueError(f"scorer {scorer_name!r} must define a non-empty `outputs` mapping")
    return [
        _parse_output_spec(name, spec, context=f"scorers.{scorer_name}.outputs.{name}")
        for name, spec in outputs_raw.items()
    ]


def _parse_metric_spec(name: str, raw_spec: Any, context: str) -> OutputSpec:
    if not isinstance(name, str) or not name.strip():
        raise ValueError("Metric names must be non-empty strings")
    spec = _mapping(raw_spec, context)
    _reject_unknown_fields(spec, {"role", "direction", "min", "max", "equals", "value_type", "description"}, context)
    role = str(spec.get("role", "")).lower()
    if role not in {"optimize", "threshold", "log"}:
        raise ValueError(f"{context}.role must be optimize, threshold, or log")
    if role != "optimize" and "direction" in spec:
        raise ValueError(f"{context}: direction is allowed only for role optimize")
    if role != "threshold" and any(key in spec for key in ("min", "max", "equals")):
        raise ValueError(f"{context}: min/max/equals bounds are allowed only for role threshold")
    value_type = str(spec.get("value_type", "float")).lower()
    if value_type not in {"float", "number", "int", "integer", "bool", "boolean", "binary", "str", "string", "structure"}:
        raise ValueError(f"{context}.value_type must be float, int, bool, string, or structure")
    translated = {key: value for key, value in spec.items() if key != "value_type"}
    translated.update(type=value_type, role="objective" if role == "optimize" else role)
    result = _parse_output_spec(name, translated, context)
    if role == "threshold":
        for bound in (result.min_value, result.max_value):
            if bound is not None and not isfinite(bound):
                raise ValueError(f"{context}: threshold bounds must be finite")
        if result.min_value is not None and result.max_value is not None and result.min_value > result.max_value:
            raise ValueError(f"{context}: threshold min must not exceed max")
    return result


def _parse_output_spec(name: str, raw_spec: Any, context: str) -> OutputSpec:
    spec = _mapping(raw_spec, context)
    _reject_unknown_fields(
        spec,
        {"type", "role", "direction", "min", "max", "equals", "description"},
        context,
    )
    role = str(spec.get("role", "log")).lower()
    if role not in {"objective", "threshold", "log"}:
        raise ValueError(f"{context}.role must be objective, threshold, or log")

    direction = spec.get("direction")
    if direction is not None:
        direction = str(direction).lower()
    if role == "objective" and direction not in {"min", "max"}:
        raise ValueError(f"{context} objective must define direction min or max")

    min_value = spec.get("min")
    max_value = spec.get("max")
    equals_value = spec.get("equals")
    if role == "threshold" and all(
        value is None for value in (min_value, max_value, equals_value)
    ):
        raise ValueError(f"{context} threshold must define min, max, and/or equals")

    result = OutputSpec(
        name=name,
        type=str(spec.get("type", "float")).lower(),
        role=role,
        direction=direction,
        min_value=None if min_value is None else float(min_value),
        max_value=None if max_value is None else float(max_value),
        equals_value=equals_value,
        description=str(spec.get("description", "")),
    )
    for bound in (result.min_value, result.max_value):
        if bound is not None and not isfinite(bound):
            raise ValueError(f"{context}: threshold bounds must be finite")
    if result.min_value is not None and result.max_value is not None and result.min_value > result.max_value:
        raise ValueError(f"{context}: threshold min must not exceed max")
    if result.equals_value is not None and result.type in {"float", "number", "int", "integer"}:
        result.coerce(result.equals_value)
    return result


def _flatten_scorer_outputs(scorers: list[ScorerConfig]) -> list[OutputSpec]:
    outputs: list[OutputSpec] = []
    seen: dict[str, str] = {}
    for scorer in scorers:
        for spec in scorer.outputs:
            output_name = qualified_output_name(scorer, spec.name)
            if output_name in seen:
                raise ValueError(
                    f"Output {output_name!r} is declared by both scorer {seen[output_name]!r} "
                    f"and scorer {scorer.name!r}"
                )
            seen[output_name] = scorer.name
            outputs.append(replace(spec, name=output_name))
    return outputs


def qualified_output_name(scorer: ScorerConfig, output_name: str) -> str:
    """Return the run-global name for one scorer-local output."""
    if scorer.namespace is not None:
        return f"{scorer.namespace}.{output_name}"
    if scorer.output_prefix is None:
        return output_name
    return f"{scorer.output_prefix}_{output_name}"


def _outputs_to_dict(outputs: list[OutputSpec] | tuple[OutputSpec, ...]) -> dict[str, Any]:
    return {
        spec.name: {
            key: value
            for key, value in {
                "type": spec.type,
                "role": spec.role,
                "direction": spec.direction,
                "min": spec.min_value,
                "max": spec.max_value,
                "equals": spec.equals_value,
                "description": spec.description or None,
            }.items()
            if value is not None
        }
        for spec in outputs
    }


def _metrics_to_dict(outputs: list[OutputSpec] | tuple[OutputSpec, ...]) -> dict[str, Any]:
    """Serialize public metric roles while keeping internal objective semantics."""
    return {
        spec.name: {
            key: value for key, value in {
                "role": "optimize" if spec.role == "objective" else spec.role,
                "direction": spec.direction if spec.role == "objective" else None,
                "min": spec.min_value if spec.role == "threshold" else None,
                "max": spec.max_value if spec.role == "threshold" else None,
                "equals": spec.equals_value if spec.role == "threshold" else None,
                "value_type": spec.type if spec.type != "float" else None,
                "description": spec.description or None,
            }.items() if value is not None
        } for spec in outputs
    }


def _parse_crossover(raw: Any) -> tuple[str, float]:
    if raw is None:
        return "uniform", 0.9
    spec = _mapping(raw, "optimizer.crossover")
    _reject_unknown_fields(spec, {"type", "probability"}, "optimizer.crossover")
    return (
        str(spec.get("type", "uniform")).lower(),
        float(spec.get("probability", 0.9)),
    )


def _parse_mutation(raw: Any) -> tuple[str, float | None]:
    if raw is None:
        return "integer_reset", None
    spec = _mapping(raw, "optimizer.mutation")
    _reject_unknown_fields(spec, {"type", "probability"}, "optimizer.mutation")
    probability_raw = spec.get("probability", "auto")
    probability = (
        None
        if probability_raw is None
        or (
            isinstance(probability_raw, str)
            and probability_raw.strip().lower() == "auto"
        )
        else float(probability_raw)
    )
    return str(spec.get("type", "integer_reset")).lower(), probability


def _parse_reference_directions(raw: Any) -> ReferenceDirectionsConfig | None:
    if raw is None:
        return None
    spec = _mapping(raw, "optimizer.reference_directions")
    _reject_unknown_fields(
        spec,
        {"method", "n_partitions", "n_points", "scaling", "seed"},
        "optimizer.reference_directions",
    )
    method = str(spec.get("method", "")).strip().lower()
    return ReferenceDirectionsConfig(
        method=method,
        n_partitions=(
            None if spec.get("n_partitions") is None else int(spec["n_partitions"])
        ),
        n_points=None if spec.get("n_points") is None else int(spec["n_points"]),
        scaling=None if spec.get("scaling") is None else float(spec["scaling"]),
        seed=None if spec.get("seed") is None else int(spec["seed"]),
    )


def _reference_directions_to_dict(
    config: ReferenceDirectionsConfig | None,
) -> dict[str, Any] | None:
    if config is None:
        return None
    return {
        key: value
        for key, value in {
            "method": config.method,
            "n_partitions": config.n_partitions,
            "n_points": config.n_points,
            "scaling": config.scaling,
            "seed": config.seed,
        }.items()
        if value is not None
    }


def _parse_constraint_handling(raw: Any) -> ConstraintHandlingConfig:
    if raw is None:
        return ConstraintHandlingConfig()
    spec = _mapping(raw, "optimizer.constraint_handling")
    _reject_unknown_fields(
        spec,
        {"method", "normalization", "initial_epsilon", "relaxation_generations", "curve"},
        "optimizer.constraint_handling",
    )
    return ConstraintHandlingConfig(
        method=str(spec.get("method", "feasibility_first")).strip().lower(),
        normalization=str(spec.get("normalization", "none")).strip().lower(),
        initial_epsilon=float(spec.get("initial_epsilon", 0.0)),
        relaxation_generations=int(spec.get("relaxation_generations", 0)),
        curve=str(spec.get("curve", "linear")).strip().lower(),
    )


def _constraint_handling_to_dict(
    config: ConstraintHandlingConfig,
) -> dict[str, Any]:
    return {
        "method": config.method,
        "normalization": config.normalization,
        "initial_epsilon": config.initial_epsilon,
        "relaxation_generations": config.relaxation_generations,
        "curve": config.curve,
    }


def _parse_early_stopping(raw: Any) -> EarlyStoppingConfig:
    if raw is None:
        return EarlyStoppingConfig(enabled=False)
    spec = _mapping(raw, "optimizer.early_stopping")
    _reject_unknown_fields(
        spec,
        {"enabled", "metric", "min_generations", "patience", "rel_tol", "abs_tol"},
        "optimizer.early_stopping",
    )
    enabled = _parse_bool(spec.get("enabled", True))
    metric = str(spec.get("metric", "hypervolume")).lower()
    if metric not in {"hypervolume", "best_feasible_objective"}:
        raise ValueError(
            "optimizer.early_stopping.metric must be hypervolume or "
            "best_feasible_objective"
        )

    min_generations = int(spec.get("min_generations", 0))
    patience = int(spec.get("patience", 0))
    rel_tol = float(spec.get("rel_tol", 0.0))
    abs_tol = float(spec.get("abs_tol", 0.0))
    if enabled:
        if min_generations < 1:
            raise ValueError("optimizer.early_stopping.min_generations must be >= 1")
        if patience < 1:
            raise ValueError("optimizer.early_stopping.patience must be >= 1")
        if rel_tol < 0:
            raise ValueError("optimizer.early_stopping.rel_tol must be >= 0")
        if abs_tol < 0:
            raise ValueError("optimizer.early_stopping.abs_tol must be >= 0")

    return EarlyStoppingConfig(
        enabled=enabled,
        metric=metric,
        min_generations=min_generations,
        patience=patience,
        rel_tol=rel_tol,
        abs_tol=abs_tol,
    )


def _parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, int) and value in {0, 1}:
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"1", "true", "yes", "y", "on"}:
            return True
        if normalized in {"0", "false", "no", "n", "off"}:
            return False
    raise ValueError(f"Expected a boolean value, got {value!r}")


def _mapping(value: Any, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"`{name}` must be a mapping")
    return value


def _reject_unknown_fields(
    raw: dict[str, Any], allowed: set[str], context: str
) -> None:
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ValueError(f"{context} contains unknown fields: " + ", ".join(unknown))


def _required_output(outputs: dict[str, Any], name: str) -> Any:
    if name not in outputs:
        raise RuntimeError(f"Scorer did not return required output {name!r}")
    return outputs[name]


def expand_env_vars(value: Any) -> Any:
    if isinstance(value, str):
        return _ENV_VAR_PATTERN.sub(_env_replacement, value)
    if isinstance(value, list):
        return [expand_env_vars(item) for item in value]
    if isinstance(value, dict):
        # A shell command belongs to Bash, including ${local}, $X, pipelines
        # and command substitutions. Expanding it here would silently erase
        # shell-local variables before the process ever starts.
        return {key: item if key == "command" and isinstance(item, str) else expand_env_vars(item)
                for key, item in value.items()}
    return value


def _env_replacement(match: re.Match[str]) -> str:
    name = match.group(1)
    default = match.group(2)
    if name in os.environ:
        return os.environ[name]
    if default is not None:
        return default
    return ""


def _format_value(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)
