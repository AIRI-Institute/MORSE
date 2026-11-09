from __future__ import annotations

import json
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from morse.construct import ConstructDesignSpace

INITIAL_POPULATION_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class InitialPopulationCandidate:
    id: str
    vector: tuple[int, ...]
    sequence: str | None
    metadata: dict[str, Any]


@dataclass(frozen=True)
class InitialPopulationManifest:
    path: Path
    design_space_fingerprint: str
    candidates: tuple[InitialPopulationCandidate, ...]
    metadata: dict[str, Any]
    resolved_design_space_fingerprint: str | None = None
    relocated: bool = False
    source_sha256: str = ""


def load_initial_population(
    path: Path,
    design_space: ConstructDesignSpace,
    population_size: int,
) -> InitialPopulationManifest:
    if not path.exists():
        raise FileNotFoundError(f"Initial population file does not exist: {path}")
    try:
        source_bytes = path.read_bytes()
        raw = json.loads(source_bytes)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Initial population file is not valid JSON: {path}") from exc
    if not isinstance(raw, dict):
        raise ValueError("Initial population file must contain a JSON object")
    if raw.get("schema_version") != INITIAL_POPULATION_SCHEMA_VERSION:
        raise ValueError(
            "Initial population file must define "
            f"schema_version: {INITIAL_POPULATION_SCHEMA_VERSION}"
        )

    fingerprint = raw.get("design_space_fingerprint")
    if not isinstance(fingerprint, str) or not fingerprint:
        raise ValueError(
            "Initial population file must define design_space_fingerprint"
        )
    current_fingerprint = design_space.fingerprint()
    portable = "portable_design_fingerprint" in raw
    if fingerprint != current_fingerprint and not portable:
        raise ValueError(
            "Initial population design space differs from the run config "
            f"(file fingerprint {fingerprint}, current fingerprint "
            f"{current_fingerprint})"
        )
    if portable and raw.get("portable_design_fingerprint") != design_space.portable_fingerprint():
        raise ValueError("Portable initial population design differs from the run config; all design choices and variant ordering must remain identical")

    candidates_raw = raw.get("candidates")
    if not isinstance(candidates_raw, list):
        raise ValueError("Initial population file candidates must be a list")
    if len(candidates_raw) > population_size:
        raise ValueError(
            "Initial population contains more candidates than population_size "
            f"({len(candidates_raw)} > {population_size})"
        )

    lower = design_space.lower_bounds
    upper = design_space.upper_bounds
    candidates: list[InitialPopulationCandidate] = []
    seen_ids: set[str] = set()
    for index, item in enumerate(candidates_raw):
        if not isinstance(item, dict):
            raise ValueError(
                f"Initial population candidates[{index}] must be an object"
            )
        candidate_id = item.get("id")
        if not isinstance(candidate_id, str) or not candidate_id.strip():
            raise ValueError(
                f"Initial population candidates[{index}].id must be a non-empty string"
            )
        if candidate_id in seen_ids:
            raise ValueError(
                f"Initial population contains duplicate candidate id {candidate_id!r}"
            )
        seen_ids.add(candidate_id)

        vector_raw = item.get("vector")
        if not isinstance(vector_raw, list) or len(vector_raw) != design_space.n_var:
            raise ValueError(
                f"Initial population candidate {candidate_id!r} vector must contain "
                f"{design_space.n_var} integers"
            )
        if any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in vector_raw
        ):
            raise ValueError(
                f"Initial population candidate {candidate_id!r} vector must contain "
                "integers"
            )
        vector = tuple(int(value) for value in vector_raw)
        for position, (value, minimum, maximum) in enumerate(
            zip(vector, lower, upper, strict=True)
        ):
            if value < minimum or value > maximum:
                raise ValueError(
                    f"Initial population candidate {candidate_id!r} variable "
                    f"{position}={value} is outside {minimum}..{maximum}"
                )

        sequence = item.get("sequence")
        if sequence is not None and (
            not isinstance(sequence, str) or not sequence.strip()
        ):
            raise ValueError(
                f"Initial population candidate {candidate_id!r} sequence must be "
                "a non-empty string"
            )
        if portable:
            if sequence is None or sequence != design_space.decode(vector).sequence:
                raise ValueError(f"Portable initial population candidate {candidate_id!r} must exactly match its stored sequence")
        metadata = item.get("metadata", {})
        if not isinstance(metadata, dict):
            raise ValueError(
                f"Initial population candidate {candidate_id!r} metadata must be an object"
            )
        candidates.append(
            InitialPopulationCandidate(
                id=candidate_id,
                vector=vector,
                sequence=None if sequence is None else sequence.upper().replace("U", "T"),
                metadata=dict(metadata),
            )
        )

    metadata = raw.get("metadata", {})
    if not isinstance(metadata, dict):
        raise ValueError("Initial population metadata must be an object")
    return InitialPopulationManifest(
        path=path,
        design_space_fingerprint=fingerprint,
        candidates=tuple(candidates),
        metadata=dict(metadata),
        resolved_design_space_fingerprint=current_fingerprint,
        relocated=fingerprint != current_fingerprint,
        source_sha256=hashlib.sha256(source_bytes).hexdigest(),
    )


def verify_initial_population(algorithm: Any, problem: Any, manifest: InitialPopulationManifest) -> dict[str, Any]:
    """Audit generation one without supplying vectors or drawing random numbers."""
    expected = [candidate.sequence for candidate in manifest.candidates]
    if hashlib.sha256(manifest.path.read_bytes()).hexdigest() != manifest.source_sha256:
        raise ValueError("Initial population audit file changed after loading")
    if any(sequence is None for sequence in expected):
        raise ValueError("Initial population audit requires stored full DNA sequences")
    observed = [problem.design_space.decode([int(round(x)) for x in individual.X]).sequence
                for individual in algorithm.pop]
    evaluated = [record.sequence for record in problem.records]
    if len(observed) != len(expected) or sorted(observed) != sorted(expected):
        raise ValueError("Generation 1 differs from the frozen initial population audit")
    if evaluated != expected:
        raise ValueError("Generation 1 evaluation order differs from the frozen initial population audit")
    if problem.constraint_handling.method != "feasibility_first" or problem.threshold_epsilon(1) != 0:
        raise ValueError("Initial population audit requires unchanged strict feasibility_first constraints")
    passing = 0
    for individual in algorithm.pop:
        construct = problem.design_space.decode([int(round(x)) for x in individual.X])
        record = problem.cache[construct.fingerprint()]
        if not record.structural_validation.valid:
            raise ValueError("Structural validation failed in audited generation 1")
        constraints = list(problem.constraint_values(record, generation=1))
        if list(individual.G) != constraints:
            raise ValueError("Generation 1 constraint values differ from the hard thresholds")
        passing += int(all(value <= 0 for value in constraints))
    digest = lambda sequences: hashlib.sha256("\n".join(sequences).encode()).hexdigest()
    return {"verified": True, "generation": 1, "population_size": len(expected),
            "passing": passing, "threshold_epsilon": 0, "evaluation_order_verified": True,
            "expected_sequence_set_sha256": digest(sorted(expected)),
            "observed_sequence_set_sha256": digest(sorted(observed)),
            "expected_evaluation_order_sha256": digest(expected),
            "observed_evaluation_order_sha256": digest(evaluated),
            "source_sha256": manifest.source_sha256}
