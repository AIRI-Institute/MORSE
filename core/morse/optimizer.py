from __future__ import annotations

import copy
import csv
import hashlib
import json
import os
import pickle
import time
import warnings
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable, Protocol

import numpy as np
import cloudpickle
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.core.callback import Callback
from pymoo.core.population import Population
from pymoo.core.problem import Problem
from pymoo.core.sampling import Sampling
from pymoo.core.mutation import Mutation
from pymoo.core.termination import Termination
from pymoo.indicators.hv import HV
from pymoo.operators.crossover.ux import UniformCrossover
from pymoo.optimize import minimize
from pymoo.termination import get_termination
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from pymoo.util.ref_dirs import get_reference_directions

from morse.codon import gc_fraction
from morse.config import (
    ConstraintHandlingConfig,
    EarlyStoppingConfig,
    OptimizerConfig,
    OutputSpec,
    ReferenceDirectionsConfig,
    ThresholdResult,
)
from morse.construct import (
    ConstructDesignSpace,
    VaccineConstruct,
)
from morse.fasta import FastaRecord, write_fasta
from morse.initial_population import (
    InitialPopulationCandidate,
    load_initial_population,
    verify_initial_population,
)
from morse.run_log import RUN_LOG_FILENAME, read_run_log_tables, write_run_log
from morse.scorers import ConstructCandidate, ScoreResult
from morse.validators import (
    StructuralValidationResult,
    StructuralValidator,
    Validator,
)

INVALID_OBJECTIVE_PENALTY = 1e12
PROGRESS_FILENAME = "progress.json"
CHECKPOINT_FILENAME = "optimizer-checkpoint.pkl"


class Scorer(Protocol):
    def score(self, candidates: list[ConstructCandidate]) -> list[ScoreResult]:
        ...


@dataclass
class EvaluationRecord:
    eval_index: int
    candidate_id: str
    construct: VaccineConstruct
    outputs: dict[str, Any]
    objective_values: dict[str, float]
    changed_codons: int
    gc_fraction: float
    valid: bool
    structural_validation: StructuralValidationResult
    generation: int | None = None
    threshold_results: list[ThresholdResult] = field(default_factory=list)

    @property
    def sequence(self) -> str:
        return self.construct.sequence

    @property
    def key(self) -> str:
        return self.construct.fingerprint()


@dataclass
class ScorerTimingRecord:
    invocation_index: int
    scorer: str
    batch_size: int
    elapsed_seconds: float
    seconds_per_candidate: float
    candidate_ids: list[str]


@dataclass
class ResumeData:
    final_population: list[list[int]]
    cached_records: list[EvaluationRecord]
    chain_evaluation_rows: list[dict[str, Any]]
    chain_population_rows: list[dict[str, Any]]
    chain_scorer_timing_rows: list[dict[str, Any]]
    chain_initial_population_rows: list[dict[str, Any]]
    previous_runs: list[dict[str, Any]]
    completed_generations: int


@dataclass
class RunResult:
    output_dir: Path
    run_state_path: Path
    summary: dict[str, Any]


class ConstructSampling(Sampling):
    def __init__(
        self,
        original_vector: list[int],
        seed: int | None = None,
        seed_vectors: list[list[int]] | None = None,
        include_original: bool = True,
        excluded_vectors: list[list[int]] | None = None,
    ) -> None:
        super().__init__()
        self.original_vector = np.asarray(original_vector, dtype=int)
        self.seed_vectors = [
            np.asarray(vector, dtype=int) for vector in (seed_vectors or [])
        ]
        self.include_original = include_original
        self.excluded_vectors = {
            tuple(int(value) for value in vector)
            for vector in (excluded_vectors or [])
        }
        self.rng = np.random.default_rng(seed)

    def _do(self, problem: Problem, n_samples: int, **kwargs: Any) -> np.ndarray:
        xl = np.asarray(problem.xl, dtype=int)
        xu = np.asarray(problem.xu, dtype=int)
        samples = np.empty((n_samples, problem.n_var), dtype=int)
        row = 0
        seen: set[tuple[int, ...]] = set()
        for vector in self.seed_vectors:
            if row >= n_samples:
                break
            clipped = np.clip(vector, xl, xu).astype(int)
            key = tuple(int(value) for value in clipped)
            if key in seen:
                continue
            samples[row] = clipped
            seen.add(key)
            row += 1
        if row < n_samples and self.include_original:
            key = tuple(int(value) for value in self.original_vector)
            if key not in seen:
                samples[row] = self.original_vector
                seen.add(key)
                row += 1
        while row < n_samples:
            for _ in range(10_000):
                candidate = self.rng.integers(xl, xu + 1)
                key = tuple(int(value) for value in candidate)
                if key in seen or key in self.excluded_vectors:
                    continue
                samples[row] = candidate
                seen.add(key)
                row += 1
                break
            else:
                raise RuntimeError(
                    "Could not randomly fill the initial population without reusing "
                    "a rejected supplied candidate"
                )
        return samples


class IntegerResetMutation(Mutation):
    def __init__(self, probability: float | None = None) -> None:
        super().__init__()
        self.probability = probability

    def _do(
        self, problem: Problem, X: np.ndarray, random_state: Any = None, **kwargs: Any
    ) -> np.ndarray:
        result = np.asarray(X, dtype=int).copy()
        probability = (
            1.0 / max(problem.n_var, 1)
            if self.probability is None
            else self.probability
        )
        if probability <= 0.0:
            return result
        for row in range(len(result)):
            changed = False
            for column in range(problem.n_var):
                if (
                    random_state.random() >= probability
                    or problem.xl[column] == problem.xu[column]
                ):
                    continue
                result[row, column] = _different_random_value(
                    int(result[row, column]),
                    int(problem.xl[column]),
                    int(problem.xu[column]),
                    random_state,
                )
                changed = True
            if not changed:
                mutable = np.flatnonzero(np.asarray(problem.xu) > np.asarray(problem.xl))
                if len(mutable):
                    column = int(random_state.choice(mutable))
                    result[row, column] = _different_random_value(
                        int(result[row, column]),
                        int(problem.xl[column]),
                        int(problem.xu[column]),
                        random_state,
                    )
        return result


class VaccineOptimizationProblem(Problem):
    def __init__(
        self,
        design_space: ConstructDesignSpace,
        scorer: Scorer,
        output_specs: list[OutputSpec],
        cached_records: list[EvaluationRecord] | None = None,
        validator: Validator | None = None,
        constraint_handling: ConstraintHandlingConfig = ConstraintHandlingConfig(),
    ) -> None:
        self.output_specs = output_specs
        self.objectives = [spec for spec in output_specs if spec.role == "objective"]
        self.thresholds = [spec for spec in output_specs if spec.role == "threshold"]
        if not self.objectives:
            raise ValueError("At least one output must have role: objective")

        super().__init__(
            n_var=design_space.n_var,
            n_obj=len(self.objectives),
            n_ieq_constr=1 + len(self.thresholds),
            xl=np.asarray(design_space.lower_bounds, dtype=int),
            xu=np.asarray(design_space.upper_bounds, dtype=int),
            vtype=int,
            requires_kwargs=True,
        )
        self.design_space = design_space
        self.scorer = scorer
        self.validator = validator or StructuralValidator()
        self.constraint_handling = constraint_handling
        self.cache: dict[str, EvaluationRecord] = {}
        self.records: list[EvaluationRecord] = []
        self.scorer_timings: list[ScorerTimingRecord] = []
        self._next_eval_index = 1
        if cached_records:
            self.preload_cache(cached_records)

    def preload_cache(self, records: list[EvaluationRecord]) -> None:
        max_candidate_index = 0
        for record in records:
            self.cache[record.key] = record
            if record.candidate_id.startswith("cand_"):
                try:
                    max_candidate_index = max(
                        max_candidate_index, int(record.candidate_id.removeprefix("cand_"))
                    )
                except ValueError:
                    pass
        self._next_eval_index = max(self._next_eval_index, max_candidate_index + 1)

    def _evaluate(self, X: np.ndarray, out: dict[str, Any], *args: Any, **kwargs: Any) -> None:
        X = np.atleast_2d(X)
        constructs = [self.design_space.decode(np.rint(row).astype(int)) for row in X]
        self.score_constructs(constructs)
        generation = _algorithm_generation(kwargs.get("algorithm"))
        out["F"] = np.asarray(
            [
                [
                    self.cache[item.fingerprint()].objective_values[spec.name]
                    for spec in self.objectives
                ]
                for item in constructs
            ],
            dtype=float,
        )
        out["G"] = np.asarray(
            [
                self.constraint_values(
                    self.cache[item.fingerprint()], generation=generation
                )
                for item in constructs
            ],
            dtype=float,
        )

    def constraint_values(
        self,
        record: EvaluationRecord,
        *,
        generation: int,
    ) -> list[float]:
        return _constraint_values(
            record,
            len(self.thresholds),
            epsilon=self.threshold_epsilon(generation),
            normalization=self.constraint_handling.normalization,
        )

    def threshold_epsilon(self, generation: int) -> float:
        return self.constraint_handling.epsilon(generation)

    def refresh_population_constraints(self, population: Any, generation: int) -> None:
        if population is None or len(population) == 0:
            return
        constraints: list[list[float]] = []
        for row in np.atleast_2d(population.get("X")):
            construct = self.design_space.decode(np.rint(row).astype(int))
            record = self.cache[construct.fingerprint()]
            constraints.append(self.constraint_values(record, generation=generation))
        population.set("G", np.asarray(constraints, dtype=float))
        population.set("CV", None)

    def score_constructs(self, constructs: list[VaccineConstruct]) -> None:
        unseen: list[ConstructCandidate] = []
        validation_by_candidate: dict[str, StructuralValidationResult] = {}
        for construct in constructs:
            key = construct.fingerprint()
            if key in self.cache:
                continue
            candidate_id = f"cand_{self._next_eval_index:06d}"
            self._next_eval_index += 1
            structural_validation = self.validator.validate(construct, self.design_space)
            if not structural_validation.valid:
                record = self._structurally_invalid_record(
                    candidate_id, construct, structural_validation
                )
                self.cache[key] = record
                self.records.append(record)
                continue
            validation_by_candidate[candidate_id] = structural_validation
            unseen.append(ConstructCandidate(candidate_id, construct))

        if not unseen:
            return

        scorer_start = time.perf_counter()
        score_results = self.scorer.score(unseen)
        elapsed_seconds = time.perf_counter() - scorer_start
        component_timings = getattr(self.scorer, "last_timings", None)
        if component_timings:
            for timing in component_timings:
                self.scorer_timings.append(
                    ScorerTimingRecord(
                        invocation_index=len(self.scorer_timings) + 1,
                        scorer=str(timing.scorer),
                        batch_size=int(timing.batch_size),
                        elapsed_seconds=float(timing.elapsed_seconds),
                        seconds_per_candidate=float(timing.seconds_per_candidate),
                        candidate_ids=list(timing.candidate_ids),
                    )
                )
        else:
            self.scorer_timings.append(
                ScorerTimingRecord(
                    invocation_index=len(self.scorer_timings) + 1,
                    scorer=_scorer_name(self.scorer),
                    batch_size=len(unseen),
                    elapsed_seconds=elapsed_seconds,
                    seconds_per_candidate=elapsed_seconds / len(unseen),
                    candidate_ids=[candidate.candidate_id for candidate in unseen],
                )
            )

        for result, candidate in zip(score_results, unseen, strict=True):
            construct = candidate.construct
            outputs = self._coerce_outputs(result.outputs)
            threshold_results = [spec.check_threshold(outputs) for spec in self.thresholds]
            valid = all(check.valid for check in threshold_results)
            objective_values = self._objective_values(outputs)

            record = EvaluationRecord(
                eval_index=int(candidate.candidate_id.removeprefix("cand_")),
                candidate_id=result.candidate_id,
                construct=construct,
                outputs=outputs,
                objective_values=objective_values,
                changed_codons=self.design_space.changed_codon_count(construct),
                gc_fraction=float(outputs.get("gc_fraction", gc_fraction(construct.sequence))),
                valid=valid,
                structural_validation=validation_by_candidate[candidate.candidate_id],
                threshold_results=threshold_results,
            )
            self.cache[record.key] = record
            self.records.append(record)

    def _coerce_outputs(self, outputs: dict[str, Any]) -> dict[str, Any]:
        coerced = dict(outputs)
        for spec in self.output_specs:
            if spec.name in coerced:
                coerced[spec.name] = spec.coerce(coerced[spec.name])
        return coerced

    def _objective_values(self, outputs: dict[str, Any]) -> dict[str, float]:
        return {spec.name: spec.objective_value(outputs) for spec in self.objectives}

    def _structurally_invalid_record(
        self,
        candidate_id: str,
        construct: VaccineConstruct,
        validation: StructuralValidationResult,
    ) -> EvaluationRecord:
        threshold_results = [
            ThresholdResult(
                output=spec.name,
                valid=False,
                observed=None,
                min_value=spec.min_value,
                max_value=spec.max_value,
                equals_value=spec.equals_value,
                reason="not evaluated because structural validation failed",
            )
            for spec in self.thresholds
        ]
        return EvaluationRecord(
            eval_index=int(candidate_id.removeprefix("cand_")),
            candidate_id=candidate_id,
            construct=construct,
            outputs={},
            objective_values={
                spec.name: INVALID_OBJECTIVE_PENALTY for spec in self.objectives
            },
            changed_codons=_safe_changed_codon_count(self.design_space, construct),
            gc_fraction=_safe_gc_fraction(construct.sequence),
            valid=False,
            structural_validation=validation,
            threshold_results=threshold_results,
        )


class NormalizedHypervolumeTracker:
    """Track comparable hypervolume values using the first feasible scale."""

    reference_value = 1.1

    def __init__(self) -> None:
        self._ideal: np.ndarray | None = None
        self._span: np.ndarray | None = None

    def value(self, algorithm: Any) -> float | None:
        population = getattr(algorithm, "pop", None)
        if population is None:
            return None
        try:
            values = np.asarray(population.get("F"), dtype=float)
        except Exception:
            return None
        try:
            raw_constraints = population.get("G")
        except Exception:
            raw_constraints = None
        if values.ndim != 2 or values.size == 0 or values.shape[1] < 2:
            return None
        valid_mask = np.all(np.isfinite(values), axis=1)
        if raw_constraints is not None:
            constraints = np.asarray(raw_constraints, dtype=float)
            if constraints.size:
                if constraints.ndim == 1:
                    constraints = constraints.reshape(-1, 1)
                if constraints.shape[0] != values.shape[0]:
                    return None
                valid_mask &= np.all(constraints <= 0, axis=1)
        values = values[valid_mask]
        if values.size == 0:
            return None
        return self.normalized_value(values)

    def normalized_value(self, values: np.ndarray) -> float:
        if self._ideal is None or self._span is None:
            ideal = np.min(values, axis=0)
            nadir = np.max(values, axis=0)
            span = nadir - ideal
            fallback_span = np.maximum(np.abs(ideal) * 0.01, 1.0)
            span = np.where(span > 1e-12, span, fallback_span)
            self._ideal = ideal
            self._span = span

        normalized = (values - self._ideal) / self._span
        normalized = np.minimum(normalized, self.reference_value)
        reference = np.full(normalized.shape[1], self.reference_value)
        return float(HV(ref_point=reference)(normalized))


class PopulationHistory(Callback):
    def __init__(
        self,
        problem: VaccineOptimizationProblem,
        *,
        progress_path: Path | None = None,
        total_generations: int | None = None,
        hypervolume_tracker: NormalizedHypervolumeTracker | None = None,
    ) -> None:
        super().__init__()
        self.problem = problem
        self.rows: list[dict[str, Any]] = []
        self.progress_path = progress_path
        self.total_generations = total_generations
        self.started_at = time.monotonic()
        self.progress_history: list[dict[str, Any]] = []
        self.latest_snapshot: dict[str, Any] | None = None
        self.hypervolume_tracker = (
            hypervolume_tracker or NormalizedHypervolumeTracker()
        )
        self.best_hypervolume: float | None = None

    def notify(self, algorithm: Any) -> None:
        generation = int(algorithm.n_gen)
        threshold_epsilon = self.problem.threshold_epsilon(generation)
        population_records: list[EvaluationRecord] = []
        selection_feasible_count = 0
        for record in self.problem.records:
            if record.generation is None:
                record.generation = generation
        for rank, individual in enumerate(algorithm.pop, start=1):
            construct = self.problem.design_space.decode(np.rint(individual.X).astype(int))
            record = self.problem.cache[construct.fingerprint()]
            population_records.append(record)
            selection_feasible_count += int(
                np.all(np.asarray(individual.G, dtype=float) <= 0.0)
            )
            row: dict[str, Any] = {
                "generation": generation,
                "rank": rank,
                "candidate_id": record.candidate_id,
                "changed_codons": record.changed_codons,
                "gc_fraction": record.gc_fraction,
                "five_prime_utr": record.construct.five_prime_utr.name,
                "three_prime_utr": record.construct.three_prime_utr.name,
                "stop_codon": record.construct.stop_codon,
                "poly_a_length": record.construct.poly_a_length,
                "threshold_epsilon": threshold_epsilon,
                "selection_feasible": bool(
                    np.all(np.asarray(individual.G, dtype=float) <= 0.0)
                ),
                "selection_constraint_violation": float(
                    np.asarray(individual.CV, dtype=float).reshape(-1)[0]
                ),
            }
            _append_record_outputs(row, record, self.problem.objectives)
            self.rows.append(row)
        if self.progress_path is not None:
            hypervolume = self.hypervolume_tracker.value(algorithm)
            if hypervolume is not None:
                self.best_hypervolume = (
                    hypervolume
                    if self.best_hypervolume is None
                    else max(self.best_hypervolume, hypervolume)
                )
            snapshot = _generation_progress_snapshot(
                self.problem,
                population_records,
                generation=generation,
                total_generations=self.total_generations,
                selection_feasible_count=selection_feasible_count,
                threshold_epsilon=threshold_epsilon,
                elapsed_seconds=time.monotonic() - self.started_at,
                hypervolume=hypervolume,
                best_hypervolume=self.best_hypervolume,
            )
            self.latest_snapshot = snapshot
            self.progress_history.append(_compact_progress_snapshot(snapshot))
            _write_progress_json(
                self.progress_path,
                {
                    "schema_version": 1,
                    **snapshot,
                    "history": self.progress_history,
                },
            )

    def hypervolume_summary(self) -> dict[str, Any]:
        history = [
            {
                "generation": int(item["generation"]),
                "hypervolume": item.get("hypervolume"),
                "best_hypervolume": item.get("best_hypervolume"),
            }
            for item in self.progress_history
            if item.get("hypervolume") is not None
        ]
        final_hypervolume = (
            None if not history else history[-1]["hypervolume"]
        )
        return {
            "available": bool(history),
            "normalization": "first selection-feasible population",
            "reference_point": self.hypervolume_tracker.reference_value,
            "final": final_hypervolume,
            "best": self.best_hypervolume,
            "history": history,
        }


def _generation_progress_snapshot(
    problem: VaccineOptimizationProblem,
    population_records: list[EvaluationRecord],
    *,
    generation: int,
    total_generations: int | None,
    selection_feasible_count: int,
    threshold_epsilon: float,
    elapsed_seconds: float,
    hypervolume: float | None,
    best_hypervolume: float | None,
) -> dict[str, Any]:
    population_size = len(population_records)
    structural_valid_count = sum(
        record.structural_validation.valid for record in population_records
    )
    threshold_valid_count = sum(
        all(check.valid for check in record.threshold_results)
        for record in population_records
    )
    hard_feasible_count = sum(
        record.structural_validation.valid
        and all(check.valid for check in record.threshold_results)
        for record in population_records
    )

    objectives = []
    for spec in problem.objectives:
        values = [
            record.outputs[spec.name]
            for record in population_records
            if spec.name in record.outputs
        ]
        objectives.append(
            {
                "name": spec.name,
                "direction": spec.direction,
                **_numeric_progress_summary(values, direction=spec.direction),
            }
        )

    thresholds = []
    for spec in problem.thresholds:
        checks = [
            check
            for record in population_records
            for check in record.threshold_results
            if check.output == spec.name
        ]
        thresholds.append(
            {
                "name": spec.name,
                "min": spec.min_value,
                "max": spec.max_value,
                "equals": spec.equals_value,
                "passed": sum(check.valid for check in checks),
                "population_size": population_size,
                "pass_fraction": _fraction(
                    sum(check.valid for check in checks),
                    population_size,
                ),
                **_numeric_progress_summary(
                    [
                        check.observed
                        for check in checks
                        if check.observed is not None
                    ]
                ),
            }
        )

    logged_metrics = [
        {
            "name": spec.name,
            **_numeric_progress_summary([
                record.outputs[spec.name]
                for record in population_records
                if spec.name in record.outputs
            ]),
        }
        for spec in problem.output_specs
        if spec.role == "log" and spec.type in {"float", "int"}
    ]

    timings: dict[str, dict[str, Any]] = {}
    for timing in problem.scorer_timings:
        aggregate = timings.setdefault(
            timing.scorer,
            {
                "scorer": timing.scorer,
                "calls": 0,
                "candidates": 0,
                "elapsed_seconds": 0.0,
            },
        )
        aggregate["calls"] += 1
        aggregate["candidates"] += timing.batch_size
        aggregate["elapsed_seconds"] += timing.elapsed_seconds

    return {
        "generation": generation,
        "total_generations": total_generations,
        "progress_fraction": (
            None
            if total_generations is None
            else min(generation / max(total_generations, 1), 1.0)
        ),
        "evaluated_candidates": len(problem.records),
        "population_size": population_size,
        "structural_valid_count": structural_valid_count,
        "structural_valid_fraction": _fraction(
            structural_valid_count,
            population_size,
        ),
        "threshold_valid_count": threshold_valid_count,
        "threshold_valid_fraction": _fraction(
            threshold_valid_count,
            population_size,
        ),
        "hard_feasible_count": hard_feasible_count,
        "hard_feasible_fraction": _fraction(
            hard_feasible_count,
            population_size,
        ),
        "selection_feasible_count": selection_feasible_count,
        "selection_feasible_fraction": _fraction(
            selection_feasible_count,
            population_size,
        ),
        "threshold_epsilon": threshold_epsilon,
        "hypervolume": hypervolume,
        "best_hypervolume": best_hypervolume,
        "elapsed_seconds": elapsed_seconds,
        "objectives": objectives,
        "thresholds": thresholds,
        "logged_metrics": logged_metrics,
        "scorer_timings": list(timings.values()),
        "utr_frequencies": {
            field: dict(sorted(Counter(
                getattr(record.construct, field).name
                for record in population_records
            ).items()))
            for field in ("five_prime_utr", "three_prime_utr")
        },
        "updated_at_unix": time.time(),
    }


def _numeric_progress_summary(
    values: list[Any],
    *,
    direction: str | None = None,
) -> dict[str, float | None]:
    numeric: list[float] = []
    for value in values:
        try:
            parsed = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(parsed):
            numeric.append(parsed)
    if not numeric:
        return {
            "count": 0,
            "mean": None,
            "std": None,
            "best": None,
            "median": None,
            "observed_min": None,
            "observed_max": None,
        }
    best = min(numeric) if direction == "min" else max(numeric)
    return {
        "count": len(numeric),
        "mean": float(np.mean(np.asarray(numeric, dtype=float))),
        "std": float(np.std(np.asarray(numeric, dtype=float), ddof=0)),
        "best": best,
        "median": float(np.median(np.asarray(numeric, dtype=float))),
        "observed_min": min(numeric),
        "observed_max": max(numeric),
    }


def _compact_progress_snapshot(snapshot: dict[str, Any]) -> dict[str, Any]:
    return {
        key: snapshot[key]
        for key in (
            "generation",
            "total_generations",
            "population_size",
            "evaluated_candidates",
            "structural_valid_count",
            "threshold_valid_count",
            "hard_feasible_count",
            "hard_feasible_fraction",
            "selection_feasible_count",
            "selection_feasible_fraction",
            "threshold_epsilon",
            "hypervolume",
            "best_hypervolume",
            "elapsed_seconds",
            "objectives",
            "thresholds",
            "logged_metrics",
            "utr_frequencies",
        )
    }


def _write_progress_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(path.suffix + ".tmp")
    temporary_path.write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n"
    )
    temporary_path.replace(path)


def _fraction(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


class _ThresholdAwareEvolutionMixin:
    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self._hard_feasible_archive: Population | None = None
        self._hard_feasible_archive_reinserted = False

    def _infill(self) -> Any:
        generation = _algorithm_generation(self)
        self._update_hard_feasible_archive(self.pop)
        self.problem.refresh_population_constraints(self.pop, generation)
        return super()._infill()

    def _advance(self, infills: Any = None, **kwargs: Any) -> None:
        population = self.pop
        if infills is not None:
            population = Population.merge(population, infills)
        self._update_hard_feasible_archive(population)
        generation = _algorithm_generation(self)
        if (
            self.problem.threshold_epsilon(generation) <= 0.0
            and not self._hard_feasible_archive_reinserted
            and self._hard_feasible_archive is not None
        ):
            archived = [
                copy.deepcopy(individual)
                for individual in self._hard_feasible_archive
                if not any(
                    np.array_equal(individual.X, current.X)
                    for current in population
                )
            ]
            if archived:
                population = Population.merge(
                    population,
                    Population.create(*archived),
                )
            self._hard_feasible_archive_reinserted = True
        self.problem.refresh_population_constraints(population, generation)
        self.pop = self.survival.do(
            self.problem,
            population,
            n_survive=self.pop_size,
            algorithm=self,
            random_state=self.random_state,
            **kwargs,
        )

    def _update_hard_feasible_archive(self, population: Any) -> None:
        if (
            population is None
            or self._hard_feasible_archive_reinserted
            or self.problem.constraint_handling.method != "adaptive_epsilon"
        ):
            return
        candidates = (
            []
            if self._hard_feasible_archive is None
            else [copy.deepcopy(individual) for individual in self._hard_feasible_archive]
        )
        candidates.extend(
            copy.deepcopy(individual)
            for individual in population
            if self._is_hard_feasible(individual)
        )
        if not candidates:
            return

        unique: dict[tuple[int, ...], Any] = {}
        for individual in candidates:
            key = tuple(int(value) for value in np.rint(individual.X))
            unique[key] = individual
        candidates = list(unique.values())
        objective_values = np.asarray(
            [np.asarray(individual.F, dtype=float) for individual in candidates],
            dtype=float,
        )
        if self.problem.n_obj == 1:
            keep = [int(np.argmin(objective_values[:, 0]))]
        else:
            keep = [
                int(index)
                for index in NonDominatedSorting().do(
                    objective_values,
                    only_non_dominated_front=True,
                )
            ]
        self._hard_feasible_archive = Population.create(
            *(candidates[index] for index in keep)
        )

    def _is_hard_feasible(self, individual: Any) -> bool:
        construct = self.problem.design_space.decode(
            np.rint(individual.X).astype(int)
        )
        record = self.problem.cache.get(construct.fingerprint())
        return bool(record is not None and record.valid)


class ThresholdAwareGA(_ThresholdAwareEvolutionMixin, GA):
    pass


class ThresholdAwareNSGA2(_ThresholdAwareEvolutionMixin, NSGA2):
    pass


class ThresholdAwareNSGA3(_ThresholdAwareEvolutionMixin, NSGA3):
    pass


class HypervolumeEarlyStoppingTermination(Termination):
    def __init__(
        self,
        max_generations: int,
        early_stopping: EarlyStoppingConfig,
        *,
        hypervolume_tracker: NormalizedHypervolumeTracker | None = None,
    ) -> None:
        super().__init__()
        self.max_generations = max_generations
        self.early_stopping = early_stopping
        self.completed_generations = 0
        self.stop_reason: str | None = None
        self.history: list[dict[str, Any]] = []
        self.hypervolume_tracker = (
            hypervolume_tracker or NormalizedHypervolumeTracker()
        )
        self._best_hv: float | None = None
        self._last_significant_improvement_generation: int | None = None

    def _update(self, algorithm: Any) -> float:
        generation = int(getattr(algorithm, "n_gen", 0) or 0)
        self.completed_generations = max(self.completed_generations, generation)
        if generation >= self.max_generations:
            self.stop_reason = f"reached max generations ({self.max_generations})"
            return 1.0

        if generation < self.early_stopping.min_generations:
            return min(generation / max(self.max_generations, 1), 0.999)

        hv = self._hypervolume(algorithm)
        if hv is None:
            return min(generation / max(self.max_generations, 1), 0.999)

        significant = False
        if self._best_hv is None:
            self._best_hv = hv
            self._last_significant_improvement_generation = generation
            significant = True
        else:
            improvement = hv - self._best_hv
            relative = improvement / max(abs(self._best_hv), 1e-12)
            if (
                improvement > 0
                and improvement >= self.early_stopping.abs_tol
                and relative >= self.early_stopping.rel_tol
            ):
                self._best_hv = hv
                self._last_significant_improvement_generation = generation
                significant = True

        last_improvement = self._last_significant_improvement_generation or generation
        stale_generations = generation - last_improvement
        self.history.append(
            {
                "generation": generation,
                "hypervolume": hv,
                "best_hypervolume": self._best_hv,
                "significant_improvement": significant,
                "stale_generations": stale_generations,
            }
        )
        if stale_generations >= self.early_stopping.patience:
            self.stop_reason = (
                "early stopping: hypervolume plateau "
                f"for {stale_generations} generations "
                f"(rel_tol={self.early_stopping.rel_tol:g}, "
                f"abs_tol={self.early_stopping.abs_tol:g})"
            )
            return 1.0

        return min(generation / max(self.max_generations, 1), 0.999)

    def _hypervolume(self, algorithm: Any) -> float | None:
        return self.hypervolume_tracker.value(algorithm)

    def _normalized_hypervolume(self, values: np.ndarray) -> float:
        return self.hypervolume_tracker.normalized_value(values)

    def status(self) -> dict[str, Any]:
        return {
            "enabled": True,
            "metric": self.early_stopping.metric,
            "min_generations": self.early_stopping.min_generations,
            "patience": self.early_stopping.patience,
            "rel_tol": self.early_stopping.rel_tol,
            "abs_tol": self.early_stopping.abs_tol,
            "completed_generations": self.completed_generations,
            "stopped_early": bool(
                self.stop_reason and self.stop_reason.startswith("early stopping")
            ),
            "stop_reason": self.stop_reason,
            "hypervolume_history": self.history,
        }


class BestFeasibleObjectiveEarlyStoppingTermination(Termination):
    def __init__(
        self,
        max_generations: int,
        early_stopping: EarlyStoppingConfig,
    ) -> None:
        super().__init__()
        self.max_generations = max_generations
        self.early_stopping = early_stopping
        self.completed_generations = 0
        self.stop_reason: str | None = None
        self.history: list[dict[str, Any]] = []
        self._best_value: float | None = None
        self._last_significant_improvement_generation: int | None = None

    def _update(self, algorithm: Any) -> float:
        generation = _algorithm_generation(algorithm)
        self.completed_generations = max(self.completed_generations, generation)
        current_best = self._best_hard_feasible_value(algorithm)
        problem = algorithm.problem
        significant = False

        if current_best is not None:
            if self._best_value is None:
                self._best_value = current_best
                self._last_significant_improvement_generation = generation
                significant = True
            else:
                improvement = self._best_value - current_best
                required = max(
                    self.early_stopping.abs_tol,
                    self.early_stopping.rel_tol
                    * max(abs(self._best_value), 1e-12),
                )
                if improvement > 0.0 and improvement >= required:
                    self._best_value = current_best
                    self._last_significant_improvement_generation = generation
                    significant = True

        last_improvement = self._last_significant_improvement_generation
        stale_generations = (
            None if last_improvement is None else generation - last_improvement
        )
        self.history.append(
            {
                "generation": generation,
                "best_feasible_objective": current_best,
                "best_seen_feasible_objective": self._best_value,
                "significant_improvement": significant,
                "stale_generations": stale_generations,
                "threshold_epsilon": problem.threshold_epsilon(generation),
            }
        )

        if generation >= self.max_generations:
            self.stop_reason = f"reached max generations ({self.max_generations})"
            return 1.0
        if (
            generation >= self.early_stopping.min_generations
            and stale_generations is not None
            and stale_generations >= self.early_stopping.patience
        ):
            self.stop_reason = (
                "early stopping: best hard-feasible objective plateau "
                f"for {stale_generations} generations "
                f"(rel_tol={self.early_stopping.rel_tol:g}, "
                f"abs_tol={self.early_stopping.abs_tol:g})"
            )
            return 1.0
        return min(generation / max(self.max_generations, 1), 0.999)

    def _best_hard_feasible_value(self, algorithm: Any) -> float | None:
        problem = algorithm.problem
        population = getattr(algorithm, "pop", None)
        if population is None:
            return None
        values: list[float] = []
        objective = problem.objectives[0]
        for row in np.atleast_2d(population.get("X")):
            construct = problem.design_space.decode(np.rint(row).astype(int))
            record = problem.cache[construct.fingerprint()]
            if record.valid:
                values.append(record.objective_values[objective.name])
        return min(values) if values else None

    def status(self) -> dict[str, Any]:
        return {
            "enabled": True,
            "metric": self.early_stopping.metric,
            "min_generations": self.early_stopping.min_generations,
            "patience": self.early_stopping.patience,
            "rel_tol": self.early_stopping.rel_tol,
            "abs_tol": self.early_stopping.abs_tol,
            "completed_generations": self.completed_generations,
            "stopped_early": bool(
                self.stop_reason and self.stop_reason.startswith("early stopping")
            ),
            "stop_reason": self.stop_reason,
            "objective_history": self.history,
        }


def _build_algorithm(
    problem: VaccineOptimizationProblem,
    optimizer: OptimizerConfig,
    sampling: ConstructSampling,
) -> tuple[Any, dict[str, Any]]:
    common_parameters = {
        "pop_size": optimizer.population_size,
        "n_offsprings": optimizer.n_offsprings,
        "sampling": sampling,
        "crossover": UniformCrossover(prob=optimizer.crossover_probability),
        "mutation": IntegerResetMutation(optimizer.mutation_probability),
        "eliminate_duplicates": optimizer.eliminate_duplicates,
    }
    algorithm_names = {
        "ga": "pymoo GA",
        "nsga2": "pymoo NSGA-II",
        "nsga3": "pymoo NSGA-III",
    }
    metadata: dict[str, Any] = {
        "type": optimizer.type,
        "name": algorithm_names[optimizer.type],
        "parameters": {
            "n_offsprings": optimizer.n_offsprings,
            "eliminate_duplicates": optimizer.eliminate_duplicates,
            "crossover": {
                "type": optimizer.crossover_type,
                "probability": optimizer.crossover_probability,
            },
            "mutation": {
                "type": optimizer.mutation_type,
                "probability": (
                    "auto"
                    if optimizer.mutation_probability is None
                    else optimizer.mutation_probability
                ),
            },
            "constraint_handling": {
                "method": optimizer.constraint_handling.method,
                "normalization": optimizer.constraint_handling.normalization,
                "initial_epsilon": optimizer.constraint_handling.initial_epsilon,
                "relaxation_generations": (
                    optimizer.constraint_handling.relaxation_generations
                ),
                "curve": optimizer.constraint_handling.curve,
                "hard_feasible_archive": (
                    optimizer.constraint_handling.method == "adaptive_epsilon"
                ),
            },
        },
    }
    if optimizer.type == "ga":
        if problem.n_obj != 1:
            raise ValueError(
                "optimizer.type ga requires exactly one output with role: objective"
            )
        return ThresholdAwareGA(**common_parameters), metadata
    if optimizer.type == "nsga2":
        return ThresholdAwareNSGA2(**common_parameters), metadata

    if problem.n_obj < 2:
        raise ValueError("optimizer.type nsga3 requires at least 2 objectives")
    reference_config = optimizer.reference_directions
    if reference_config is None:  # guarded by OptimizerConfig, kept for type safety
        raise ValueError("optimizer.type nsga3 requires reference_directions")
    reference_directions = _make_reference_directions(
        reference_config,
        n_objectives=problem.n_obj,
        default_seed=optimizer.seed,
    )
    if len(reference_directions) < problem.n_obj:
        raise ValueError(
            "NSGA-III requires at least one reference direction per objective "
            f"({len(reference_directions)} < {problem.n_obj})"
        )
    if optimizer.population_size < len(reference_directions):
        raise ValueError(
            "optimizer.population_size must be >= the number of NSGA-III reference "
            f"directions ({optimizer.population_size} < {len(reference_directions)})"
        )
    metadata["parameters"]["reference_directions"] = {
        **_reference_directions_metadata(reference_config, optimizer.seed),
        "count": len(reference_directions),
        "dimensions": int(reference_directions.shape[1]),
    }
    return (
        ThresholdAwareNSGA3(
            ref_dirs=reference_directions,
            **common_parameters,
        ),
        metadata,
    )


def _make_reference_directions(
    config: ReferenceDirectionsConfig,
    *,
    n_objectives: int,
    default_seed: int,
) -> np.ndarray:
    if config.method == "das-dennis":
        kwargs: dict[str, Any] = {"n_partitions": config.n_partitions}
        if config.scaling is not None:
            kwargs["scaling"] = config.scaling
        directions = get_reference_directions(
            "das-dennis", n_objectives, **kwargs
        )
    else:
        directions = get_reference_directions(
            "energy",
            n_objectives,
            config.n_points,
            seed=default_seed if config.seed is None else config.seed,
        )
    directions = np.asarray(directions, dtype=float)
    if directions.ndim != 2 or directions.shape[1] != n_objectives:
        raise ValueError(
            "Generated NSGA-III reference directions have unexpected shape "
            f"{directions.shape}; expected (*, {n_objectives})"
        )
    if len(directions) == 0 or not np.all(np.isfinite(directions)):
        raise ValueError("Generated NSGA-III reference directions are empty or non-finite")
    return directions


def _reference_directions_metadata(
    config: ReferenceDirectionsConfig, default_seed: int
) -> dict[str, Any]:
    return {
        key: value
        for key, value in {
            "method": config.method,
            "n_partitions": config.n_partitions,
            "n_points": config.n_points,
            "scaling": config.scaling,
            "seed": (
                (default_seed if config.seed is None else config.seed)
                if config.method == "energy"
                else None
            ),
        }.items()
        if value is not None
    }


def _optimizer_checkpoint_fingerprint(
    scoring_fingerprint: str,
    design_space: ConstructDesignSpace,
    optimizer: OptimizerConfig,
    output_specs: list[OutputSpec],
) -> str:
    optimizer_payload = asdict(optimizer)
    # The generated warm-start manifest is an input to generation zero only.
    # Its output-directory path is intentionally excluded so the canonical
    # config can resume the exact checkpoint after a process/runtime restart.
    optimizer_payload.pop("initial_population", None)
    audit_path = optimizer_payload.pop("initial_population_audit", None)
    if audit_path:
        optimizer_payload["initial_population_audit_sha256"] = hashlib.sha256(Path(audit_path).read_bytes()).hexdigest()
    if optimizer.initial_population_policy == "retain_all" and optimizer.initial_population:
        optimizer_payload["initial_population_sha256"] = hashlib.sha256(Path(optimizer.initial_population).read_bytes()).hexdigest()
    payload = {
        "schema_version": 1,
        "scoring_fingerprint": scoring_fingerprint,
        "design_space_fingerprint": design_space.fingerprint(),
        "optimizer": optimizer_payload,
        "outputs": [asdict(output) for output in output_specs],
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _write_optimizer_checkpoint(
    path: Path,
    fingerprint: str,
    algorithm: Any,
    problem: VaccineOptimizationProblem,
    *,
    elapsed_seconds: float,
) -> None:
    """Atomically persist the exact pymoo state after a completed generation."""

    callback = getattr(algorithm, "callback", None)
    saved_problem = algorithm.problem
    saved_callback_problem = getattr(callback, "problem", None)
    algorithm.problem = None
    if isinstance(callback, PopulationHistory):
        callback.problem = None  # type: ignore[assignment]
    try:
        algorithm_bytes = cloudpickle.dumps(
            algorithm,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    finally:
        algorithm.problem = saved_problem
        if isinstance(callback, PopulationHistory):
            callback.problem = saved_callback_problem
    payload = {
        "schema_version": 1,
        "fingerprint": fingerprint,
        "algorithm": algorithm_bytes,
        "records": list(problem.records),
        "scorer_timings": list(problem.scorer_timings),
        "elapsed_seconds": float(elapsed_seconds),
        "generation": int(getattr(algorithm, "n_gen", 0) or 0),
    }
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(pickle.dumps(payload, protocol=pickle.HIGHEST_PROTOCOL))
    os.replace(temporary, path)


def _read_optimizer_checkpoint(path: Path, fingerprint: str) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Cannot resume: missing optimizer checkpoint {path}")
    try:
        payload = pickle.loads(path.read_bytes())
    except Exception as exc:
        raise ValueError(f"Cannot read optimizer checkpoint {path}: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise ValueError("Optimizer checkpoint has an unsupported schema")
    if payload.get("fingerprint") != fingerprint:
        raise ValueError(
            "Optimizer checkpoint is incompatible with the current construct, "
            "optimizer, scorer configuration, or model assets"
        )
    if not isinstance(payload.get("algorithm"), bytes):
        raise ValueError("Optimizer checkpoint does not contain an algorithm state")
    if not isinstance(payload.get("records"), list):
        raise ValueError("Optimizer checkpoint does not contain evaluation records")
    if not isinstance(payload.get("scorer_timings"), list):
        raise ValueError("Optimizer checkpoint does not contain scorer timings")
    return payload


def run_optimization(
    input_manifest: Path,
    output_dir: Path,
    scorer: Scorer,
    output_specs: list[OutputSpec],
    optimizer: OptimizerConfig,
    run_name: str,
    design_space: ConstructDesignSpace,
    resume_from: Path | None = None,
    scoring_fingerprint: str | None = None,
    validator: Validator | None = None,
    on_progress: Callable[[dict[str, Any]], None] | None = None,
    _interrupt_after_generations: int | None = None,
) -> RunResult:
    if on_progress is not None and not callable(on_progress):
        raise TypeError("on_progress must be callable or None")
    checkpoint_path = output_dir / CHECKPOINT_FILENAME
    exact_resume = bool(
        resume_from is not None
        and (resume_from / CHECKPOINT_FILENAME).is_file()
        and not (resume_from / "run_state.json").is_file()
    )
    if exact_resume and output_dir.resolve() != resume_from.resolve():
        raise ValueError(
            "An interrupted optimizer checkpoint must resume in its original output directory"
        )
    if (
        resume_from is not None
        and not exact_resume
        and optimizer.constraint_handling.method == "adaptive_epsilon"
    ):
        raise ValueError(
            "Resume is not supported with adaptive_epsilon constraint handling"
        )
    if resume_from is not None and not exact_resume and optimizer.initial_population is not None:
        raise ValueError(
            "resume_from and optimizer.initial_population cannot be used together"
        )
    pop_size = optimizer.population_size
    generations = optimizer.generations
    seed = optimizer.seed
    validator = validator or StructuralValidator()
    scoring_fingerprint = scoring_fingerprint or _output_specs_fingerprint(output_specs)
    audit_manifest = (load_initial_population(Path(optimizer.initial_population_audit), design_space, pop_size)
                      if optimizer.initial_population_audit else None)
    if audit_manifest is not None and (
        len(audit_manifest.candidates) != pop_size or any(candidate.sequence is None for candidate in audit_manifest.candidates)
    ):
        raise ValueError("Initial population audit requires the complete population and stored full DNA sequences")
    output_dir.mkdir(parents=True, exist_ok=True)

    checkpoint_fingerprint = _optimizer_checkpoint_fingerprint(
        scoring_fingerprint,
        design_space,
        optimizer,
        output_specs,
    )
    checkpoint_payload = (
        _read_optimizer_checkpoint(
            resume_from / CHECKPOINT_FILENAME,
            checkpoint_fingerprint,
        )
        if exact_resume and resume_from is not None
        else None
    )

    resume_data = (
        _load_resume_data(
            resume_from, output_specs, design_space, scoring_fingerprint, validator
        )
        if resume_from and not exact_resume
        else None
    )
    initial_seed_vectors: list[list[int]] = []
    excluded_initial_vectors: list[list[int]] = []
    initial_population_rows: list[dict[str, Any]] = []
    initial_population_summary: dict[str, Any]
    if resume_data:
        initial_seed_vectors.extend(resume_data.final_population)
    problem = VaccineOptimizationProblem(
        design_space,
        scorer,
        output_specs,
        cached_records=(
            list(checkpoint_payload["records"])
            if checkpoint_payload is not None
            else (resume_data.cached_records if resume_data else None)
        ),
        validator=validator,
        constraint_handling=optimizer.constraint_handling,
    )
    if checkpoint_payload is not None:
        problem.records = list(checkpoint_payload["records"])
        problem.scorer_timings = list(checkpoint_payload["scorer_timings"])
    if optimizer.type == "ga" and problem.n_obj != 1:
        raise ValueError(
            "optimizer.type ga requires exactly one output with role: objective"
        )
    if optimizer.type == "nsga3" and problem.n_obj < 2:
        raise ValueError("optimizer.type nsga3 requires at least 2 objectives")
    if resume_data:
        initial_population_summary = {
            "source": "resume",
            "path": None,
            "provided": len(initial_seed_vectors),
            "accepted": len(initial_seed_vectors),
            "rejected": 0,
            "rejected_ids": [],
            "random_fill": max(0, pop_size - len(initial_seed_vectors)),
        }
    elif optimizer.initial_population is not None:
        (
            initial_seed_vectors,
            excluded_initial_vectors,
            initial_population_rows,
            initial_population_summary,
        ) = _prepare_initial_population(
            Path(optimizer.initial_population),
            design_space,
            problem,
            pop_size,
            policy=optimizer.initial_population_policy,
        )
    else:
        initial_population_summary = {
            "source": "reference",
            "path": None,
            "provided": 1,
            "accepted": 1,
            "rejected": 0,
            "rejected_ids": [],
            "random_fill": max(0, pop_size - 1),
        }
    new_algorithm, algorithm_metadata = _build_algorithm(
        problem,
        optimizer,
        ConstructSampling(
            design_space.original_vector,
            seed=seed,
            seed_vectors=initial_seed_vectors,
            include_original=(
                resume_data is None and optimizer.initial_population is None
            ),
            excluded_vectors=excluded_initial_vectors,
        ),
    )
    if checkpoint_payload is not None:
        algorithm = cloudpickle.loads(checkpoint_payload["algorithm"])
        algorithm.problem = problem
        history = algorithm.callback
        if not isinstance(history, PopulationHistory):
            raise ValueError("Optimizer checkpoint contains an incompatible callback")
        history.problem = problem
        history.progress_path = output_dir / PROGRESS_FILENAME
        history.total_generations = generations
        elapsed = float(checkpoint_payload.get("elapsed_seconds", 0.0))
        history.started_at = time.monotonic() - max(0.0, elapsed)
        hypervolume_tracker = history.hypervolume_tracker
    else:
        algorithm = new_algorithm
        hypervolume_tracker = NormalizedHypervolumeTracker()
        history = PopulationHistory(
            problem,
            progress_path=output_dir / PROGRESS_FILENAME,
            total_generations=generations,
            hypervolume_tracker=hypervolume_tracker,
        )
    original_construct = design_space.source_construct
    problem.score_constructs([original_construct])
    termination: Any
    early_stopping = optimizer.early_stopping
    if early_stopping.enabled:
        if early_stopping.metric == "best_feasible_objective":
            termination = BestFeasibleObjectiveEarlyStoppingTermination(
                generations, early_stopping
            )
        else:
            termination = HypervolumeEarlyStoppingTermination(
                generations,
                early_stopping,
                hypervolume_tracker=hypervolume_tracker,
            )
    else:
        termination = get_termination("n_gen", generations)
    if checkpoint_payload is None:
        algorithm.setup(
            problem,
            termination=termination,
            seed=seed,
            callback=history,
            save_history=False,
            verbose=False,
        )
    generations_this_process = 0
    while algorithm.has_next():
        algorithm.next()
        if audit_manifest is not None and int(algorithm.n_gen) == 2:
            # pymoo increments n_gen after the generation-one callback. The
            # audit reads actual scoring order and does not feed the sampler.
            audit = verify_initial_population(algorithm, problem, audit_manifest)
            audit["initialization_policy"] = optimizer.initial_population_policy if optimizer.initial_population else "native"
            audit["audit_only"] = optimizer.initial_population is None
            (output_dir / "initial_population_check.json").write_text(json.dumps(audit, indent=2) + "\n")
        generations_this_process += 1
        _write_optimizer_checkpoint(
            checkpoint_path,
            checkpoint_fingerprint,
            algorithm,
            problem,
            elapsed_seconds=time.monotonic() - history.started_at,
        )
        if on_progress is not None and history.latest_snapshot is not None:
            # Presentation observers may contain IPython handles or user closures.
            # Keep them outside the algorithm, and notify only after a safe resume
            # point exists. Mutating the payload must not affect the optimizer.
            payload = {
                "schema_version": 1,
                **history.latest_snapshot,
                "history": history.progress_history,
            }
            status = getattr(algorithm.termination, "status", None)
            if callable(status):
                early_status = status()
                for key in ("hypervolume_history", "objective_history"):
                    if early_status.get(key):
                        early_status["stale_generations"] = early_status[key][-1].get("stale_generations")
                payload["early_stopping"] = early_status
            else:
                payload["early_stopping"] = {"enabled": False}
            on_progress(copy.deepcopy(payload))
        if (
            _interrupt_after_generations is not None
            and generations_this_process >= _interrupt_after_generations
        ):
            raise RuntimeError("simulated interruption after checkpoint")
    result = algorithm.result()
    result.algorithm = algorithm
    completed_generations = _completed_generation_count(result, history, generations)
    used_termination = getattr(getattr(result, "algorithm", None), "termination", termination)
    if isinstance(
        used_termination,
        (
            HypervolumeEarlyStoppingTermination,
            BestFeasibleObjectiveEarlyStoppingTermination,
        ),
    ):
        used_termination.completed_generations = completed_generations
        early_stopping_status = used_termination.status()
    else:
        early_stopping_status = {
            "enabled": False,
            "completed_generations": completed_generations,
            "stopped_early": False,
            "stop_reason": None,
        }

    final_population_vectors = _final_population_vectors(result)
    pareto_vectors = np.asarray(final_population_vectors, dtype=int)
    pareto_records = _nondominated_valid_records_for_vectors(problem, pareto_vectors)
    best_by_objective = _best_records_by_objective(pareto_records, problem.objectives)
    pareto_archive_records = _nondominated_valid_records(
        problem.records,
        problem.objectives,
    )
    archive_best_by_objective = _best_records_by_objective(
        pareto_archive_records,
        problem.objectives,
    )

    run_log_path = output_dir / RUN_LOG_FILENAME
    pareto_path = output_dir / "pareto.csv"
    pareto_fasta_path = output_dir / "pareto.fasta"
    pareto_constructs_path = output_dir / "pareto_constructs.json"
    pareto_archive_path = output_dir / "pareto_archive.csv"
    pareto_archive_fasta_path = output_dir / "pareto_archive.fasta"
    pareto_archive_constructs_path = output_dir / "pareto_archive_constructs.json"
    reference_fasta_path = output_dir / "reference.fasta"
    run_state_path = output_dir / "run_state.json"
    summary_path = output_dir / "summary.json"

    evaluation_rows = _evaluation_rows(problem.records, problem.objectives)
    population_rows = [dict(row) for row in history.rows]
    scorer_timing_rows = _scorer_timing_rows(problem.scorer_timings)
    run_index = len(resume_data.previous_runs) + 1 if resume_data else 1
    prior_generations = resume_data.completed_generations if resume_data else 0
    current_run_info = {
        "run_index": run_index,
        "run_name": run_name,
        "run_dir": str(output_dir),
        "optimizer": algorithm_metadata,
        "generations": completed_generations,
        "requested_generations": generations,
        "population_size": pop_size,
        "seed": seed,
        "generation_start": prior_generations + 1,
        "generation_end": prior_generations + completed_generations,
        "resumed_from": str(resume_from) if resume_from else None,
        "initial_population": initial_population_summary,
        "early_stopping": early_stopping_status,
        "hypervolume": history.hypervolume_summary(),
    }
    chain_runs = [*(resume_data.previous_runs if resume_data else []), current_run_info]
    chain_evaluation_rows = _combine_rows(
        resume_data.chain_evaluation_rows if resume_data else [],
        _tag_rows(
            _offset_generation_rows(evaluation_rows, prior_generations),
            run_name,
            run_index,
        ),
    )
    chain_population_rows = _combine_rows(
        resume_data.chain_population_rows if resume_data else [],
        _tag_rows(
            _offset_generation_rows(population_rows, prior_generations),
            run_name,
            run_index,
        ),
    )
    chain_scorer_timing_rows = _combine_rows(
        resume_data.chain_scorer_timing_rows if resume_data else [],
        _tag_rows(scorer_timing_rows, run_name, run_index),
    )
    chain_initial_population_rows = _combine_rows(
        resume_data.chain_initial_population_rows if resume_data else [],
        _tag_rows(initial_population_rows, run_name, run_index),
    )
    write_run_log(
        run_log_path,
        chain_evaluation_rows,
        chain_population_rows,
        chain_scorer_timing_rows,
        initial_population_rows=chain_initial_population_rows,
    )
    _write_record_csv(pareto_path, pareto_records, problem.objectives)
    pareto_constructs_path.write_text(
        json.dumps(
            [
                {
                    "candidate_id": record.candidate_id,
                    "construct": record.construct.to_dict(),
                    "outputs": record.outputs,
                }
                for record in pareto_records
            ],
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    write_fasta(
        pareto_fasta_path,
        [
            FastaRecord(
                header=f"{rec.candidate_id} changed_codons={rec.changed_codons} "
                + " ".join(
                    f"{spec.name}={float(rec.outputs[spec.name]):.6g}"
                    for spec in problem.objectives
                ),
                sequence=rec.sequence,
            )
            for rec in pareto_records
        ],
    )
    _write_record_csv(
        pareto_archive_path,
        pareto_archive_records,
        problem.objectives,
    )
    pareto_archive_constructs_path.write_text(
        json.dumps(
            [
                {
                    "candidate_id": record.candidate_id,
                    "construct": record.construct.to_dict(),
                    "outputs": record.outputs,
                }
                for record in pareto_archive_records
            ],
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    write_fasta(
        pareto_archive_fasta_path,
        [
            FastaRecord(
                header=f"{rec.candidate_id} changed_codons={rec.changed_codons} "
                + " ".join(
                    f"{spec.name}={float(rec.outputs[spec.name]):.6g}"
                    for spec in problem.objectives
                ),
                sequence=rec.sequence,
            )
            for rec in pareto_archive_records
        ],
    )

    for metric, best_record in best_by_objective.items():
        path = output_dir / f"best_by_{metric}.fasta"
        write_fasta(
            path,
            [
                FastaRecord(
                    header=(
                        f"{best_record.candidate_id} best_by={metric} "
                        f"changed_codons={best_record.changed_codons}"
                    ),
                    sequence=best_record.sequence,
                )
            ],
        )

    original_record = problem.cache[original_construct.fingerprint()]
    write_fasta(
        reference_fasta_path,
        [
            FastaRecord(
                header=(
                    f"reference source_construct changed_codons="
                    f"{original_record.changed_codons}"
                ),
                sequence=original_record.sequence,
            )
        ],
    )
    summary = _build_summary(
        input_manifest=input_manifest,
        design_space=design_space,
        original_record=original_record,
        pareto_records=pareto_records,
        best_by_objective=best_by_objective,
        pareto_archive_records=pareto_archive_records,
        archive_best_by_objective=archive_best_by_objective,
        output_specs=output_specs,
        pop_size=pop_size,
        generations=completed_generations,
        requested_generations=generations,
        seed=seed,
        early_stopping=early_stopping_status,
        hypervolume=history.hypervolume_summary(),
        records=problem.records,
        scorer_timings=problem.scorer_timings,
        run_name=run_name,
        run_index=run_index,
        resume_from=resume_from,
        chain_runs=chain_runs,
        chain_evaluation_rows=chain_evaluation_rows,
        chain_population_rows=chain_population_rows,
        chain_scorer_timing_rows=chain_scorer_timing_rows,
        chain_initial_population_rows=chain_initial_population_rows,
        validator_name=validator.name,
        algorithm_metadata=algorithm_metadata,
        initial_population=initial_population_summary,
    )
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    _write_run_state(
        run_state_path,
        result,
        run_name=run_name,
        output_dir=output_dir,
        resume_from=resume_from,
        chain_runs=chain_runs,
        completed_generations=prior_generations + completed_generations,
        early_stopping=early_stopping_status,
        design_space_fingerprint=design_space.fingerprint(),
        scoring_fingerprint=scoring_fingerprint,
        algorithm_metadata=algorithm_metadata,
    )

    return RunResult(
        output_dir=output_dir,
        run_state_path=run_state_path,
        summary=summary,
    )


def _prepare_initial_population(
    path: Path,
    design_space: ConstructDesignSpace,
    problem: VaccineOptimizationProblem,
    population_size: int,
    *, policy: str = "filter",
) -> tuple[
    list[list[int]],
    list[list[int]],
    list[dict[str, Any]],
    dict[str, Any],
]:
    manifest = load_initial_population(path, design_space, population_size)
    if policy == "retain_all" and (
        problem.constraint_handling.method != "feasibility_first"
        or problem.threshold_epsilon(1) != 0
        or len(manifest.candidates) != population_size
    ):
        raise ValueError("retain_all requires a complete supplied population and strict feasibility_first constraints")
    rows: list[dict[str, Any]] = []
    decoded: list[
        tuple[InitialPopulationCandidate, VaccineConstruct, dict[str, Any]]
    ] = []
    seen_vectors: set[tuple[int, ...]] = set()
    seen_constructs: set[str] = set()
    excluded_vectors: list[list[int]] = []

    for candidate in manifest.candidates:
        row: dict[str, Any] = {
            "initial_population_id": candidate.id,
            "vector": list(candidate.vector),
            "metadata": candidate.metadata,
            "source_path": str(path),
        }
        rows.append(row)
        if candidate.vector in seen_vectors:
            row.update(
                {
                    "status": "rejected_duplicate",
                    "accepted": False,
                    "reasons": ["duplicate design vector"],
                }
            )
            continue
        seen_vectors.add(candidate.vector)
        try:
            construct = design_space.decode(candidate.vector)
        except ValueError as exc:
            row.update(
                {
                    "status": "rejected_structural",
                    "accepted": False,
                    "reasons": [str(exc)],
                }
            )
            excluded_vectors.append(list(candidate.vector))
            continue
        if candidate.sequence is not None and candidate.sequence != construct.sequence:
            raise ValueError(
                f"Initial population candidate {candidate.id!r} sequence does not "
                "match its design vector"
            )
        construct_key = construct.fingerprint()
        if construct_key in seen_constructs:
            row.update(
                {
                    "status": "rejected_duplicate",
                    "accepted": False,
                    "reasons": ["duplicate assembled construct"],
                }
            )
            continue
        seen_constructs.add(construct_key)
        decoded.append((candidate, construct, row))

    if policy == "retain_all" and len(decoded) != population_size:
        raise ValueError("retain_all cannot retain duplicate or structurally invalid initial candidates")
    if decoded:
        problem.score_constructs([construct for _, construct, _ in decoded])

    accepted_vectors: list[list[int]] = []
    relaxed_accepted = 0
    retained_infeasible = 0
    for candidate, construct, row in decoded:
        record = problem.cache[construct.fingerprint()]
        row["candidate_id"] = record.candidate_id
        row["structural_valid"] = record.structural_validation.valid
        row["threshold_valid"] = all(
            result.valid for result in record.threshold_results
        )
        row["threshold_epsilon"] = problem.threshold_epsilon(1)
        selection_feasible = all(
            value <= 0.0
            for value in problem.constraint_values(record, generation=1)
        )
        row["selection_feasible"] = selection_feasible
        if policy == "retain_all":
            row["thresholds_relaxed"] = False
        if record.structural_validation.valid and selection_feasible:
            is_relaxed = not record.valid
            row.update(
                {
                    "status": (
                        "accepted_relaxed_threshold" if is_relaxed else "accepted"
                    ),
                    "accepted": True,
                    "reasons": [],
                }
            )
            accepted_vectors.append(list(candidate.vector))
            relaxed_accepted += int(is_relaxed)
        else:
            structural_reasons = [
                f"{failure.module}.{failure.code}: {failure.message}"
                for failure in record.structural_validation.failures
            ]
            threshold_reasons = [
                result.reason
                for result in record.threshold_results
                if not result.valid and result.reason
            ]
            if policy == "retain_all":
                if not record.structural_validation.valid:
                    raise ValueError("retain_all cannot retain a structurally invalid initial candidate")
                row.update(status="accepted_infeasible_initial", accepted=True, reasons=[],
                           initial_threshold_reasons=threshold_reasons, thresholds_relaxed=False)
                accepted_vectors.append(list(candidate.vector))
                retained_infeasible += 1
                continue
            row.update(
                {
                    "status": (
                        "rejected_structural"
                        if structural_reasons
                        else "rejected_threshold"
                    ),
                    "accepted": False,
                    "reasons": structural_reasons or threshold_reasons,
                }
            )
            excluded_vectors.append(list(candidate.vector))

    rejected_ids = [
        str(row["initial_population_id"])
        for row in rows
        if not row.get("accepted", False)
    ]
    if rejected_ids:
        warnings.warn(
            "Initial population candidates rejected and replaced with random "
            f"candidates: {', '.join(rejected_ids)}. See initial_population "
            "events in run_log.jsonl.gz for details.",
            UserWarning,
            stacklevel=2,
        )

    summary = {
        "source": "file",
        "path": str(path),
        "provided": len(manifest.candidates),
        "accepted": len(accepted_vectors),
        "rejected": len(rejected_ids),
        "rejected_ids": rejected_ids,
        "random_fill": max(0, population_size - len(accepted_vectors)),
        "metadata": manifest.metadata,
    }
    if relaxed_accepted:
        summary["relaxed_accepted"] = relaxed_accepted
    if policy == "retain_all":
        summary.update(retained_threshold_infeasible=retained_infeasible,
                       retention_policy="retain_all", thresholds_relaxed=False)
    if manifest.relocated:
        summary["fingerprint_relocation"] = {"historical": manifest.design_space_fingerprint,
                                            "resolved": manifest.resolved_design_space_fingerprint,
                                            "portable_design_verified": True}

    return (
        accepted_vectors,
        excluded_vectors,
        rows,
        summary,
    )


def _nondominated_valid_records_for_vectors(
    problem: VaccineOptimizationProblem, vectors: np.ndarray
) -> list[EvaluationRecord]:
    records: list[EvaluationRecord] = []
    seen: set[str] = set()
    for vector in vectors:
        construct = problem.design_space.decode(np.rint(vector).astype(int))
        key = construct.fingerprint()
        if key in seen:
            continue
        seen.add(key)
        record = problem.cache[key]
        if record.valid:
            records.append(record)
    if not records:
        return []
    objective_matrix = np.asarray(
        [
            [record.objective_values[spec.name] for spec in problem.objectives]
            for record in records
        ],
        dtype=float,
    )
    nondominated = NonDominatedSorting().do(
        objective_matrix, only_non_dominated_front=True
    )
    return [records[int(index)] for index in nondominated]


def _nondominated_valid_records(
    records: list[EvaluationRecord], objectives: list[OutputSpec]
) -> list[EvaluationRecord]:
    unique = {
        record.construct.fingerprint(): record
        for record in records
        if record.valid
    }
    candidates = list(unique.values())
    if not candidates:
        return []
    objective_matrix = np.asarray(
        [
            [record.objective_values[spec.name] for spec in objectives]
            for record in candidates
        ],
        dtype=float,
    )
    if len(objectives) == 1:
        nondominated = [int(np.argmin(objective_matrix[:, 0]))]
    else:
        nondominated = NonDominatedSorting().do(
            objective_matrix,
            only_non_dominated_front=True,
        )
    return [candidates[int(index)] for index in nondominated]


def _completed_generation_count(
    result: Any, history: PopulationHistory, requested_generations: int
) -> int:
    if history.rows:
        return max(int(row["generation"]) for row in history.rows)
    algorithm = getattr(result, "algorithm", None)
    generation = getattr(algorithm, "n_gen", None)
    if generation is not None:
        return int(generation)
    return requested_generations


def _best_records_by_objective(
    records: list[EvaluationRecord], objectives: list[OutputSpec]
) -> dict[str, EvaluationRecord]:
    best: dict[str, EvaluationRecord] = {}
    if not records:
        return best
    for spec in objectives:
        if spec.direction == "min":
            best[spec.name] = min(records, key=lambda rec: float(rec.outputs[spec.name]))
        else:
            best[spec.name] = max(records, key=lambda rec: float(rec.outputs[spec.name]))
    return best


def _append_record_outputs(
    row: dict[str, Any], record: EvaluationRecord, objectives: list[OutputSpec]
) -> None:
    row["valid"] = record.valid
    row["structural_valid"] = record.structural_validation.valid
    row["structural_binary"] = record.structural_validation.binary
    for module, valid in record.structural_validation.module_validity.items():
        row[f"structural_module_{module}_valid"] = valid
    row["structural_failures"] = json.dumps(
        [failure.to_dict() for failure in record.structural_validation.failures],
        sort_keys=True,
    )
    for check in record.threshold_results:
        prefix = f"threshold_{check.output}"
        row[f"{prefix}_valid"] = check.valid
        row[f"{prefix}_observed"] = check.observed
        row[f"{prefix}_min"] = check.min_value
        row[f"{prefix}_max"] = check.max_value
        row[f"{prefix}_equals"] = check.equals_value
        row[f"{prefix}_reason"] = check.reason
    for spec in objectives:
        row[f"objective_{spec.name}"] = record.objective_values[spec.name]
    for key, value in sorted(record.outputs.items()):
        row[f"output_{key}"] = value


def _evaluation_rows(
    records: list[EvaluationRecord], objectives: list[OutputSpec]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        row: dict[str, Any] = {
            "eval_index": record.eval_index,
            "generation": record.generation,
            "candidate_id": record.candidate_id,
            "sequence": record.sequence,
            "start_codon": record.construct.start_codon,
            "cds": record.construct.cds,
            "coding_sequence": record.construct.coding_sequence,
            "coding_sequence_with_stop": record.construct.coding_sequence_with_stop,
            "five_prime_utr": record.construct.five_prime_utr.name,
            "three_prime_utr": record.construct.three_prime_utr.name,
            "stop_codon": record.construct.stop_codon,
            "poly_a_length": record.construct.poly_a_length,
            "construct_json": json.dumps(record.construct.to_dict(), sort_keys=True),
            "changed_codons": record.changed_codons,
            "gc_fraction": record.gc_fraction,
        }
        _append_record_outputs(row, record, objectives)
        rows.append(row)
    return rows


def _scorer_timing_rows(
    records: list[ScorerTimingRecord],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        rows.append(
            {
                "invocation_index": record.invocation_index,
                "scorer": record.scorer,
                "batch_size": record.batch_size,
                "elapsed_seconds": record.elapsed_seconds,
                "seconds_per_candidate": record.seconds_per_candidate,
                "candidate_ids": ";".join(record.candidate_ids),
            }
        )
    return rows


def _write_record_csv(
    path: Path, records: list[EvaluationRecord], objectives: list[OutputSpec]
) -> None:
    _write_csv(path, _evaluation_rows(records, objectives))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _build_summary(
    input_manifest: Path,
    design_space: ConstructDesignSpace,
    original_record: EvaluationRecord,
    pareto_records: list[EvaluationRecord],
    best_by_objective: dict[str, EvaluationRecord],
    pareto_archive_records: list[EvaluationRecord],
    archive_best_by_objective: dict[str, EvaluationRecord],
    output_specs: list[OutputSpec],
    pop_size: int,
    generations: int,
    requested_generations: int,
    seed: int,
    early_stopping: dict[str, Any],
    hypervolume: dict[str, Any],
    records: list[EvaluationRecord],
    scorer_timings: list[ScorerTimingRecord],
    run_name: str,
    run_index: int,
    resume_from: Path | None,
    chain_runs: list[dict[str, Any]],
    chain_evaluation_rows: list[dict[str, Any]],
    chain_population_rows: list[dict[str, Any]],
    chain_scorer_timing_rows: list[dict[str, Any]],
    chain_initial_population_rows: list[dict[str, Any]],
    validator_name: str,
    algorithm_metadata: dict[str, Any],
    initial_population: dict[str, Any],
) -> dict[str, Any]:
    rejected_records = [record for record in records if not record.valid]
    structural_rejected_records = [
        record for record in records if not record.structural_validation.valid
    ]
    objectives = [spec for spec in output_specs if spec.role == "objective"]
    thresholds = [spec for spec in output_specs if spec.role == "threshold"]
    logs = [spec for spec in output_specs if spec.role == "log"]
    total_scorer_seconds = sum(record.elapsed_seconds for record in scorer_timings)
    total_scored_candidates = sum(record.batch_size for record in scorer_timings)
    first_scorer_invocations: dict[str, ScorerTimingRecord] = {}
    for timing in scorer_timings:
        current = first_scorer_invocations.get(timing.scorer)
        if current is None or timing.invocation_index < current.invocation_index:
            first_scorer_invocations[timing.scorer] = timing
    startup_inclusive_seconds = sum(
        timing.elapsed_seconds for timing in first_scorer_invocations.values()
    )
    chain_scorer_seconds = sum(
        float(row.get("elapsed_seconds", 0) or 0) for row in chain_scorer_timing_rows
    )
    chain_scored_candidates = sum(
        int(float(row.get("batch_size", 0) or 0)) for row in chain_scorer_timing_rows
    )
    chain_rejected = sum(
        1 for row in chain_evaluation_rows if not _parse_bool(row.get("valid", True))
    )
    chain_structural_rejected = sum(
        1
        for row in chain_evaluation_rows
        if not _parse_bool(row.get("structural_valid", True))
    )
    return {
        "status": "completed" if pareto_archive_records else "completed_no_feasible",
        "completion_message": (
            "Configured search completed."
            if pareto_archive_records
            else "Configured search completed without a hard-feasible evaluated candidate; observations and constraint failures are retained."
        ),
        "run": {
            "run_name": run_name,
            "run_index": run_index,
            "resumed_from": str(resume_from) if resume_from else None,
        },
        "run_chain": {
            "runs": chain_runs,
            "total_runs": len(chain_runs),
            "total_generations": sum(int(run.get("generations", 0)) for run in chain_runs),
        },
        "input": {
            "manifest": str(input_manifest),
            "name": design_space.name,
            "length_nt": len(design_space.source_construct.sequence),
            "start_codon_variants": len(design_space.start_codons),
            "cds_body_codons": design_space.codons.n_var,
            "protein_amino_acids": len(design_space.amino_acids),
            "five_prime_utr_variants": len(design_space.five_prime_utrs),
            "three_prime_utr_variants": len(design_space.three_prime_utrs),
            "stop_codon_variants": len(design_space.stop_codons),
            "poly_a_variants": len(design_space.poly_a_lengths),
            "part_max_lengths": {
                "five_prime_utr": max(len(item.sequence) for item in design_space.five_prime_utrs),
                "start_codon": max(map(len, design_space.start_codons)),
                "cds": len(design_space.codons.source_sequence),
                "stop_codon": max(map(len, design_space.stop_codons)),
                "three_prime_utr": max(len(item.sequence) for item in design_space.three_prime_utrs),
                "poly_a": max(design_space.poly_a_lengths),
            },
        },
        "optimizer": {
            "algorithm": algorithm_metadata["name"],
            "type": algorithm_metadata["type"],
            "parameters": algorithm_metadata["parameters"],
            "population_size": pop_size,
            "generations": generations,
            "requested_generations": requested_generations,
            "seed": seed,
            "initial_population": initial_population,
            "early_stopping": early_stopping,
            "hypervolume": hypervolume,
            "objectives": [
                {"output": spec.name, "type": spec.type, "direction": spec.direction}
                for spec in objectives
            ],
        },
        "output_policy": {
            "thresholds": [
                {
                    "output": spec.name,
                    "type": spec.type,
                    "min": spec.min_value,
                    "max": spec.max_value,
                    "equals": spec.equals_value,
                }
                for spec in thresholds
            ],
            "logs": [{"output": spec.name, "type": spec.type} for spec in logs],
            "accepted": len(records) - len(rejected_records),
            "rejected": len(rejected_records),
            "chain_accepted": len(chain_evaluation_rows) - chain_rejected,
            "chain_rejected": chain_rejected,
        },
        "structural_validation": {
            "validator": validator_name,
            "accepted": len(records) - len(structural_rejected_records),
            "rejected": len(structural_rejected_records),
            "chain_accepted": len(chain_evaluation_rows) - chain_structural_rejected,
            "chain_rejected": chain_structural_rejected,
        },
        "scorer_timing": {
            "invocations": len(scorer_timings),
            "total_seconds": total_scorer_seconds,
            "total_candidates": total_scored_candidates,
            "mean_seconds_per_candidate": (
                total_scorer_seconds / total_scored_candidates
                if total_scored_candidates
                else 0.0
            ),
            "max_invocation_seconds": max(
                (record.elapsed_seconds for record in scorer_timings), default=0.0
            ),
            "startup_inclusive_first_calls_seconds": startup_inclusive_seconds,
            "subsequent_calls_seconds": max(
                0.0,
                total_scorer_seconds - startup_inclusive_seconds,
            ),
            "chain_total_seconds": chain_scorer_seconds,
            "chain_total_candidates": chain_scored_candidates,
            "chain_mean_seconds_per_candidate": (
                chain_scorer_seconds / chain_scored_candidates
                if chain_scored_candidates
                else 0.0
            ),
        },
        "chain_outputs": {
            "evaluated_rows": len(chain_evaluation_rows),
            "population_rows": len(chain_population_rows),
            "scorer_timing_rows": len(chain_scorer_timing_rows),
            "initial_population_rows": len(chain_initial_population_rows),
        },
        "original": _record_summary(original_record),
        "pareto_count": len(pareto_records),
        "best_by_objective": {
            output: _record_summary(record) for output, record in best_by_objective.items()
        },
        "pareto_archive_count": len(pareto_archive_records),
        "archive_best_by_objective": {
            output: _record_summary(record)
            for output, record in archive_best_by_objective.items()
        },
        "reference_comparison": {
            spec.name: _reference_comparison(
                original_record,
                archive_best_by_objective[spec.name],
                spec,
            )
            for spec in objectives
            if spec.name in archive_best_by_objective
        },
        "model_applicability": _model_applicability_summary(
            records,
            output_specs,
        ),
    }


def _reference_comparison(
    reference: EvaluationRecord,
    best: EvaluationRecord,
    objective: OutputSpec,
) -> dict[str, Any]:
    reference_value = float(reference.outputs[objective.name])
    best_value = float(best.outputs[objective.name])
    improvement = (
        reference_value - best_value
        if objective.direction == "min"
        else best_value - reference_value
    )
    return {
        "direction": objective.direction,
        "reference": reference_value,
        "best_evaluated": best_value,
        "difference_in_preferred_direction": improvement,
        "beats_reference": improvement > 0.0,
        "candidate_id": best.candidate_id,
    }


def _model_applicability_summary(
    records: list[EvaluationRecord],
    output_specs: list[OutputSpec],
) -> dict[str, Any]:
    non_aug = sorted(
        {record.construct.start_codon for record in records} - {"ATG"}
    )
    affected_output_names = {
        "ribodecode_translation_log",
        "translation_efficiency_score",
        "half_life_score",
    }
    selected_outputs = sorted(
        {spec.name for spec in output_specs if spec.name in affected_output_names}
    )
    applies = bool(non_aug and selected_outputs)
    return {
        "status": "out of distribution" if applies else "in distribution",
        "non_aug_start_codons": non_aug,
        "affected_outputs": selected_outputs if applies else [],
        "warning": (
            "Non-AUG candidates are outside the documented training distribution "
            "of one or more selected model outputs"
            if applies
            else None
        ),
    }


def _scorer_name(scorer: Scorer) -> str:
    name = getattr(scorer, "name", None)
    if name:
        return str(name)
    image = getattr(scorer, "image", None)
    if image:
        return f"{scorer.__class__.__name__}:{image}"
    return scorer.__class__.__name__


def _load_resume_data(
    resume_from: Path,
    output_specs: list[OutputSpec],
    design_space: ConstructDesignSpace,
    scoring_fingerprint: str,
    validator: Validator,
) -> ResumeData:
    state_path = resume_from / "run_state.json"
    if not state_path.exists():
        raise FileNotFoundError(f"Cannot resume: missing {state_path}")
    state = json.loads(state_path.read_text())
    _validate_design_space_fingerprint(state, design_space, state_path)
    if state.get("scoring_fingerprint") != scoring_fingerprint:
        raise ValueError("Cannot resume: scorer configuration differs from the source run")
    tables = read_run_log_tables(resume_from / RUN_LOG_FILENAME)
    chain_evaluation_rows = tables["evaluation"]
    chain_population_rows = tables["population_member"]
    chain_scorer_timing_rows = tables["scorer_timing"]
    chain_initial_population_rows = tables["initial_population"]
    cached_records = [
        _record_from_log_row(row, output_specs, design_space, validator)
        for row in chain_evaluation_rows
        if row.get("construct_json")
    ]
    chain_runs = state.get("run_chain", {}).get("runs", [])
    if not chain_runs:
        chain_runs = [
            {
                "run_index": 1,
                "run_name": state.get("run_name", resume_from.name),
                "run_dir": str(resume_from),
                "generations": state.get("completed_generations", 0),
                "generation_start": 1,
                "generation_end": state.get("completed_generations", 0),
            }
        ]
    return ResumeData(
        final_population=state.get("final_population", []),
        cached_records=cached_records,
        chain_evaluation_rows=chain_evaluation_rows,
        chain_population_rows=chain_population_rows,
        chain_scorer_timing_rows=chain_scorer_timing_rows,
        chain_initial_population_rows=chain_initial_population_rows,
        previous_runs=chain_runs,
        completed_generations=int(
            state.get("completed_generations")
            or sum(int(run.get("generations", 0)) for run in chain_runs)
        ),
    )


def _record_from_log_row(
    row: dict[str, Any],
    output_specs: list[OutputSpec],
    design_space: ConstructDesignSpace,
    validator: Validator,
) -> EvaluationRecord:
    outputs = {
        key.removeprefix("output_"): _parse_log_value(value)
        for key, value in row.items()
        if key.startswith("output_") and value not in {"", None}
    }
    threshold_results: list[ThresholdResult] = []
    for spec in output_specs:
        if spec.role != "threshold":
            continue
        if spec.name in outputs:
            threshold_results.append(spec.check_threshold(outputs))
            continue
        prefix = f"threshold_{spec.name}"
        valid_key = f"{prefix}_valid"
        if valid_key not in row:
            continue
        threshold_results.append(
            ThresholdResult(
                output=spec.name,
                valid=_parse_bool(row.get(valid_key, True)),
                observed=_parse_optional_output(
                    spec, row.get(f"{prefix}_observed")
                ),
                min_value=_parse_optional_float(row.get(f"{prefix}_min")),
                max_value=_parse_optional_float(row.get(f"{prefix}_max")),
                equals_value=spec.equals_value,
                reason=str(row.get(f"{prefix}_reason", "")),
            )
        )
    threshold_valid = (
        all(check.valid for check in threshold_results)
        if threshold_results
        else _parse_bool(row.get("valid", True))
    )
    objective_values: dict[str, float] = {}
    for spec in output_specs:
        if spec.role != "objective":
            continue
        objective_key = f"objective_{spec.name}"
        if spec.name in outputs:
            objective_values[spec.name] = spec.objective_value(outputs)
        elif objective_key in row and row[objective_key] not in {"", None}:
            objective_values[spec.name] = float(row[objective_key])
        elif not threshold_valid:
            objective_values[spec.name] = INVALID_OBJECTIVE_PENALTY
        else:
            objective_values[spec.name] = spec.objective_value(outputs)
    construct_raw = json.loads(str(row["construct_json"]))
    construct = VaccineConstruct.from_dict(construct_raw)
    structural_validation = validator.validate(construct, design_space)
    valid = structural_validation.valid and threshold_valid
    return EvaluationRecord(
        eval_index=int(float(row.get("eval_index", 0) or 0)),
        candidate_id=str(row["candidate_id"]),
        construct=construct,
        outputs=outputs,
        objective_values=objective_values,
        changed_codons=int(
            float(row.get("changed_codons") or design_space.changed_codon_count(construct))
        ),
        gc_fraction=float(row.get("gc_fraction") or gc_fraction(construct.sequence)),
        valid=valid,
        structural_validation=structural_validation,
        generation=_parse_optional_int(row.get("generation")),
        threshold_results=threshold_results,
    )


def _write_run_state(
    path: Path,
    result: Any,
    run_name: str,
    output_dir: Path,
    resume_from: Path | None,
    chain_runs: list[dict[str, Any]],
    completed_generations: int,
    early_stopping: dict[str, Any],
    design_space_fingerprint: str,
    scoring_fingerprint: str,
    algorithm_metadata: dict[str, Any],
) -> None:
    path.write_text(
        json.dumps(
            {
                "run_name": run_name,
                "run_dir": str(output_dir),
                "resumed_from": str(resume_from) if resume_from else None,
                "completed_generations": completed_generations,
                "run_chain": {"runs": chain_runs},
                "early_stopping": early_stopping,
                "design_space_fingerprint": design_space_fingerprint,
                "scoring_fingerprint": scoring_fingerprint,
                "optimizer": algorithm_metadata,
                "final_population": _final_population_vectors(result),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def _final_population_vectors(result: Any) -> list[list[int]]:
    population = getattr(getattr(result, "algorithm", None), "pop", None)
    if population is None:
        return []
    vectors: list[list[int]] = []
    for individual in population:
        vectors.append([int(round(float(value))) for value in individual.X])
    return vectors


def _validate_design_space_fingerprint(
    state: dict[str, Any], design_space: ConstructDesignSpace, state_path: Path
) -> None:
    stored = state.get("design_space_fingerprint")
    if stored is None:
        raise ValueError(
            f"Cannot resume from {state_path}: run_state.json does not contain "
            "`design_space_fingerprint`; recreate the source run with the current optimizer."
        )
    current = design_space.fingerprint()
    if stored != current:
        raise ValueError(
            "Cannot resume: construct manifest or module design space differs from the source run "
            f"(run_state fingerprint {stored}, current fingerprint {current})"
        )


def _threshold_violation(
    check: ThresholdResult,
    *,
    normalization: str = "none",
) -> float:
    if check.observed is None:
        return INVALID_OBJECTIVE_PENALTY
    if check.equals_value is not None and check.observed != check.equals_value:
        return 1.0
    violation = 0.0
    observed = float(check.observed)
    if check.min_value is not None and observed < check.min_value:
        violation = max(violation, check.min_value - observed)
    if check.max_value is not None and observed > check.max_value:
        violation = max(violation, observed - check.max_value)
    if violation > 0.0 and normalization == "relative":
        violation /= _threshold_scale(check)
    return violation


def _threshold_scale(check: ThresholdResult) -> float:
    if (
        check.min_value is not None
        and check.max_value is not None
        and check.max_value > check.min_value
    ):
        return max(check.max_value - check.min_value, 1e-12)
    bounds = [
        abs(value)
        for value in (check.min_value, check.max_value)
        if value is not None
    ]
    scale = max(bounds, default=1.0)
    return scale if scale > 1e-12 else 1.0


def _constraint_values(
    record: EvaluationRecord,
    threshold_count: int,
    *,
    epsilon: float = 0.0,
    normalization: str = "none",
) -> list[float]:
    values = [0.0 if record.structural_validation.valid else 1.0]
    values.extend(
        _threshold_violation(check, normalization=normalization) - epsilon
        for check in record.threshold_results
    )
    values.extend(
        [INVALID_OBJECTIVE_PENALTY]
        * max(0, threshold_count - len(record.threshold_results))
    )
    return values[: 1 + threshold_count]


def _algorithm_generation(algorithm: Any) -> int:
    if algorithm is None:
        return 1
    raw = getattr(algorithm, "n_gen", None)
    if raw is None:
        raw = getattr(algorithm, "n_iter", None)
    return max(int(raw or 1), 1)


def _safe_changed_codon_count(
    design_space: ConstructDesignSpace, construct: VaccineConstruct
) -> int:
    try:
        return design_space.changed_codon_count(construct)
    except (ValueError, IndexError):
        return 0


def _safe_gc_fraction(sequence: str) -> float:
    try:
        return gc_fraction(sequence)
    except ValueError:
        return 0.0


def _different_random_value(
    current: int, lower: int, upper: int, random_state: Any
) -> int:
    draw = int(random_state.integers(lower, upper))
    return draw + 1 if draw >= current else draw


def _output_specs_fingerprint(output_specs: list[OutputSpec]) -> str:
    payload = [
        {
            "name": spec.name,
            "type": spec.type,
            "role": spec.role,
            "direction": spec.direction,
            "min": spec.min_value,
            "max": spec.max_value,
            "equals": spec.equals_value,
        }
        for spec in output_specs
    ]
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _tag_rows(rows: list[dict[str, Any]], run_name: str, run_index: int) -> list[dict[str, Any]]:
    tagged: list[dict[str, Any]] = []
    for row in rows:
        new_row = dict(row)
        new_row["run_name"] = run_name
        new_row["run_index"] = run_index
        tagged.append(new_row)
    return tagged


def _offset_generation_rows(
    rows: list[dict[str, Any]], generation_offset: int
) -> list[dict[str, Any]]:
    shifted: list[dict[str, Any]] = []
    for row in rows:
        new_row = dict(row)
        generation = int(float(new_row.get("generation", 0) or 0))
        new_row["chain_generation"] = generation + generation_offset
        shifted.append(new_row)
    return shifted


def _combine_rows(
    previous_rows: list[dict[str, Any]], current_rows: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    return [dict(row) for row in previous_rows] + [dict(row) for row in current_rows]


def _parse_log_value(value: Any) -> Any:
    if isinstance(value, bool | int | float):
        return value
    if value is None:
        return value
    text = str(value)
    if text == "":
        return text
    for parser in (int, float):
        try:
            return parser(text)
        except ValueError:
            continue
    if text in {"True", "False"}:
        return text == "True"
    return text


def _parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).lower() in {"1", "true", "yes"}


def _parse_optional_float(value: Any) -> float | None:
    if value in {None, ""}:
        return None
    return float(value)


def _parse_optional_int(value: Any) -> int | None:
    if value in {None, ""}:
        return None
    return int(float(value))


def _parse_optional_output(spec: OutputSpec, value: Any) -> Any | None:
    if value in {None, ""}:
        return None
    return spec.coerce(value)


def _record_summary(record: EvaluationRecord) -> dict[str, Any]:
    return {
        "candidate_id": record.candidate_id,
        "length_nt": len(record.sequence),
        "modules": {
            "five_prime_utr": record.construct.five_prime_utr.name,
            "three_prime_utr": record.construct.three_prime_utr.name,
            "stop_codon": record.construct.stop_codon,
            "poly_a_length": record.construct.poly_a_length,
        },
        "changed_codons": record.changed_codons,
        "gc_fraction": record.gc_fraction,
        "valid": record.valid,
        "structural_validation": record.structural_validation.to_dict(),
        "outputs": {
            key: value
            for key, value in record.outputs.items()
            if isinstance(value, int | float | str | bool)
        },
        "thresholds": [
            {
                "output": check.output,
                "valid": check.valid,
                "observed": check.observed,
                "min": check.min_value,
                "max": check.max_value,
                "equals": check.equals_value,
                "reason": check.reason,
            }
            for check in record.threshold_results
        ],
    }
