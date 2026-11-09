from __future__ import annotations

import heapq
from dataclasses import dataclass
from itertools import islice, product
from math import isfinite, prod
import random
from typing import Any, Iterable, Protocol

from morse.config import OutputSpec, ScorerConfig
from morse.construct import (
    PART_ORDER,
    ConstructDesignSpace,
    VaccineConstruct,
)
from morse.router import resolve_scope_ranges
from morse.scorer_scope import ScorerScopeConfig

from morse_warm_starter.scoring import ScorerSessions


@dataclass(frozen=True)
class AnchorRequest:
    mode: str
    scorer: str
    objective: str
    count: int = 1

    def __post_init__(self) -> None:
        if self.mode not in {"exhaustive", "greedy"}:
            raise ValueError("Anchor mode must be exhaustive or greedy")
        if not self.scorer or not self.objective:
            raise ValueError("Anchor scorer and objective must be non-empty")
        if self.count < 1:
            raise ValueError("Anchor count must be >= 1")


@dataclass(frozen=True)
class GenerationOptions:
    population_size: int
    seed: int
    include_original: bool = True
    anchors_only: bool = False
    max_combinations: int = 100_000
    batch_size: int = 512
    passes: int = 1
    candidate_pool: int = 100
    min_hamming_distance: int = 1

    def __post_init__(self) -> None:
        if self.population_size < 1:
            raise ValueError("population_size must be >= 1")
        if self.max_combinations < 1:
            raise ValueError("max_combinations must be >= 1")
        if self.batch_size < 1:
            raise ValueError("batch_size must be >= 1")
        if self.passes < 1:
            raise ValueError("passes must be >= 1")
        if self.candidate_pool < 1:
            raise ValueError("candidate_pool must be >= 1")
        if self.min_hamming_distance < 0:
            raise ValueError("min_hamming_distance must be >= 0")


@dataclass(frozen=True)
class GeneratedCandidate:
    id: str
    vector: tuple[int, ...]
    sequence: str
    metadata: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "vector": list(self.vector),
            "sequence": self.sequence,
            "metadata": self.metadata,
        }


class ScoringSession(Protocol):
    def score(
        self,
        constructs: list[VaccineConstruct],
        *,
        objective: str,
        batch_size: int,
    ) -> list[dict[str, Any]]:
        ...


class SessionProvider(Protocol):
    def get(self, config: ScorerConfig) -> ScoringSession:
        ...

    def close(self) -> None:
        ...


def generate_population(
    *,
    design_space: ConstructDesignSpace,
    scorer_configs: tuple[ScorerConfig, ...],
    requests: list[AnchorRequest],
    options: GenerationOptions,
    sessions: SessionProvider | None = None,
    prebuilt_candidates: Iterable[GeneratedCandidate] = (),
) -> dict[str, Any]:
    prebuilt = list(prebuilt_candidates)
    requested_slots = sum(request.count for request in requests)
    requested_slots += len(prebuilt)
    if options.include_original:
        requested_slots += 1
    if requested_slots > options.population_size:
        raise ValueError(
            "Requested anchor count exceeds population_size "
            f"({requested_slots} > {options.population_size})"
        )
    if options.candidate_pool < max(
        (request.count for request in requests),
        default=1,
    ):
        raise ValueError("candidate_pool must be >= every requested anchor count")

    configs = {config.name: config for config in scorer_configs}
    owned_sessions = sessions is None
    session_provider = sessions or ScorerSessions()
    generated: list[GeneratedCandidate] = []
    seen_vectors: set[tuple[int, ...]] = set()
    seen_sequences: set[str] = set()
    anchor_results: list[dict[str, Any]] = []

    def add(candidate: GeneratedCandidate) -> bool:
        if candidate.vector in seen_vectors or candidate.sequence in seen_sequences:
            return False
        seen_vectors.add(candidate.vector)
        seen_sequences.add(candidate.sequence)
        generated.append(candidate)
        return True

    try:
        if options.include_original:
            source = design_space.source_construct
            add(
                GeneratedCandidate(
                    id="original-001",
                    vector=tuple(design_space.original_vector),
                    sequence=source.sequence,
                    metadata={"anchor": "original"},
                )
            )

        accepted_prebuilt = 0
        for candidate in prebuilt:
            construct = design_space.decode(candidate.vector)
            if construct.sequence != candidate.sequence:
                raise ValueError(
                    f"Prebuilt candidate {candidate.id!r} sequence does not match its vector"
                )
            accepted_prebuilt += int(add(candidate))
        if prebuilt:
            anchor_results.append(
                {
                    "mode": "lineardesign",
                    "requested": len(prebuilt),
                    "proposed": len(prebuilt),
                    "accepted": accepted_prebuilt,
                    "duplicates": len(prebuilt) - accepted_prebuilt,
                }
            )

        for request in requests:
            config, objective = _resolve_request(request, configs)
            session = session_provider.get(config)
            if request.mode == "exhaustive":
                proposals = _exhaustive(
                    design_space,
                    config,
                    objective,
                    request,
                    session,
                    options,
                )
            else:
                proposals = _greedy(
                    design_space,
                    config,
                    objective,
                    request,
                    session,
                    options,
                )
            accepted = 0
            for rank, proposal in enumerate(proposals, start=1):
                accepted += int(
                    add(
                        GeneratedCandidate(
                            id=(
                                f"{request.mode}-{_safe_name(request.scorer)}-"
                                f"{_safe_name(request.objective)}-{rank:03d}"
                            ),
                            vector=proposal.vector,
                            sequence=proposal.construct.sequence,
                            metadata={
                                "anchor": request.mode,
                                "scorer": request.scorer,
                                "objective": request.objective,
                                "direction": objective.direction,
                                "objective_value": proposal.value,
                                "variable_positions": list(proposal.positions),
                                **proposal.metadata,
                            },
                        )
                    )
                )
            anchor_results.append(
                {
                    "mode": request.mode,
                    "scorer": request.scorer,
                    "objective": request.objective,
                    "requested": request.count,
                    "proposed": len(proposals),
                    "accepted": accepted,
                    "duplicates": len(proposals) - accepted,
                }
            )

        if not options.anchors_only:
            _fill_random(
                design_space,
                options,
                generated,
                seen_vectors,
                seen_sequences,
            )
    finally:
        if owned_sessions:
            session_provider.close()

    return {
        "schema_version": 1,
        "design_space_fingerprint": design_space.fingerprint(),
        "metadata": {
            "generator": "morse-warm-starter",
            "population_size": options.population_size,
            "seed": options.seed,
            "include_original": options.include_original,
            "anchors_only": options.anchors_only,
            "anchor_requests": [
                {
                    "mode": request.mode,
                    "scorer": request.scorer,
                    "objective": request.objective,
                    "count": request.count,
                }
                for request in requests
            ],
            "anchor_results": anchor_results,
            "generated_candidates": len(generated),
            "random_candidates": sum(
                candidate.metadata.get("anchor") == "random"
                for candidate in generated
            ),
        },
        "candidates": [candidate.to_dict() for candidate in generated],
    }


@dataclass(frozen=True)
class _Proposal:
    vector: tuple[int, ...]
    construct: VaccineConstruct
    value: float
    positions: tuple[int, ...]
    metadata: dict[str, Any]


def _resolve_request(
    request: AnchorRequest,
    configs: dict[str, ScorerConfig],
) -> tuple[ScorerConfig, OutputSpec]:
    try:
        config = configs[request.scorer]
    except KeyError as exc:
        raise ValueError(f"Unknown scorer {request.scorer!r}") from exc
    outputs = {output.name: output for output in config.outputs}
    try:
        objective = outputs[request.objective]
    except KeyError as exc:
        raise ValueError(
            f"Scorer {request.scorer!r} does not declare output "
            f"{request.objective!r}"
        ) from exc
    if objective.role != "objective" or objective.direction not in {"min", "max"}:
        raise ValueError(
            f"Output {request.objective!r} must be a min/max objective"
        )
    return config, objective


def _exhaustive(
    design: ConstructDesignSpace,
    config: ScorerConfig,
    objective: OutputSpec,
    request: AnchorRequest,
    session: ScoringSession,
    options: GenerationOptions,
) -> list[_Proposal]:
    positions = tuple(_variable_positions_for_scope(design, config.scope))
    lower = design.lower_bounds
    upper = design.upper_bounds
    widths = [upper[position] - lower[position] + 1 for position in positions]
    combinations_count = prod(widths) if widths else 1
    if combinations_count > options.max_combinations:
        raise ValueError(
            f"Exhaustive anchor {request.scorer}:{request.objective} expands to "
            f"{combinations_count:,} combinations, exceeding "
            f"max_combinations={options.max_combinations:,}"
        )

    values: Iterable[tuple[int, ...]]
    if positions:
        values = product(
            *(range(lower[position], upper[position] + 1) for position in positions)
        )
    else:
        values = [()]
    base = design.original_vector
    heap: list[tuple[float, int, _Proposal]] = []
    serial = 0
    for batch_values in _batched(values, options.batch_size):
        vectors: list[list[int]] = []
        constructs: list[VaccineConstruct] = []
        for combination in batch_values:
            vector = list(base)
            for position, value in zip(positions, combination, strict=True):
                vector[position] = value
            vectors.append(vector)
            constructs.append(design.decode(vector))
        scored = session.score(
            constructs,
            objective=objective.name,
            batch_size=options.batch_size,
        )
        for vector, construct, outputs in zip(
            vectors,
            constructs,
            scored,
            strict=True,
        ):
            value = _finite_objective(outputs, objective, request)
            proposal = _Proposal(
                vector=tuple(vector),
                construct=construct,
                value=value,
                positions=positions,
                metadata={"combination_count": combinations_count},
            )
            serial += 1
            item = (_desirability(value, objective), serial, proposal)
            if len(heap) < options.candidate_pool:
                heapq.heappush(heap, item)
            elif item[:2] > heap[0][:2]:
                heapq.heapreplace(heap, item)

    ranked = [
        item[2]
        for item in sorted(heap, key=lambda item: item[:2], reverse=True)
    ]
    return _select_diverse(
        ranked,
        request.count,
        options.min_hamming_distance,
    )


def _greedy(
    design: ConstructDesignSpace,
    config: ScorerConfig,
    objective: OutputSpec,
    request: AnchorRequest,
    session: ScoringSession,
    options: GenerationOptions,
) -> list[_Proposal]:
    positions = tuple(_variable_positions_for_scope(design, config.scope))
    lower = design.lower_bounds
    upper = design.upper_bounds
    current = list(design.original_vector)
    completed_passes = 0

    for _ in range(options.passes):
        changed = False
        for position in positions:
            vectors = []
            for value in range(lower[position], upper[position] + 1):
                vector = list(current)
                vector[position] = value
                vectors.append(vector)
            ranked = _score_vectors(
                design,
                vectors,
                positions,
                objective,
                request,
                session,
                options.batch_size,
                preferred=tuple(current),
            )
            if ranked and ranked[0].vector != tuple(current):
                current = list(ranked[0].vector)
                changed = True
        completed_passes += 1
        if not changed:
            break

    neighbors = [list(current)]
    for position in positions:
        for value in range(lower[position], upper[position] + 1):
            if value == current[position]:
                continue
            vector = list(current)
            vector[position] = value
            neighbors.append(vector)
    ranked = _score_vectors(
        design,
        neighbors,
        positions,
        objective,
        request,
        session,
        options.batch_size,
    )
    ranked = [
        _Proposal(
            vector=proposal.vector,
            construct=proposal.construct,
            value=proposal.value,
            positions=proposal.positions,
            metadata={"passes": completed_passes},
        )
        for proposal in ranked[: options.candidate_pool]
    ]
    return _select_diverse(
        ranked,
        request.count,
        options.min_hamming_distance,
    )


def _score_vectors(
    design: ConstructDesignSpace,
    vectors: list[list[int]],
    positions: tuple[int, ...],
    objective: OutputSpec,
    request: AnchorRequest,
    session: ScoringSession,
    batch_size: int,
    preferred: tuple[int, ...] | None = None,
) -> list[_Proposal]:
    constructs = [design.decode(vector) for vector in vectors]
    outputs = session.score(
        constructs,
        objective=objective.name,
        batch_size=batch_size,
    )
    proposals = [
        _Proposal(
            vector=tuple(vector),
            construct=construct,
            value=_finite_objective(result, objective, request),
            positions=positions,
            metadata={},
        )
        for vector, construct, result in zip(
            vectors,
            constructs,
            outputs,
            strict=True,
        )
    ]
    proposals.sort(
        key=lambda proposal: (
            _desirability(proposal.value, objective),
            preferred is not None and proposal.vector == preferred,
        ),
        reverse=True,
    )
    return proposals


def _fill_random(
    design: ConstructDesignSpace,
    options: GenerationOptions,
    generated: list[GeneratedCandidate],
    seen_vectors: set[tuple[int, ...]],
    seen_sequences: set[str],
) -> None:
    lower = design.lower_bounds
    upper = design.upper_bounds
    distinct_vectors = prod(
        maximum - minimum + 1
        for minimum, maximum in zip(lower, upper, strict=True)
    )
    if options.population_size > distinct_vectors:
        raise ValueError(
            "population_size exceeds the number of design vectors "
            f"({options.population_size} > {distinct_vectors})"
        )
    rng = random.Random(options.seed)
    attempts = 0
    random_index = 1
    while len(generated) < options.population_size:
        vector = tuple(
            rng.randint(minimum, maximum)
            for minimum, maximum in zip(lower, upper, strict=True)
        )
        attempts += 1
        if vector in seen_vectors:
            if attempts > 100_000:
                raise RuntimeError("Could not generate a unique random population")
            continue
        construct = design.decode(vector)
        if construct.sequence in seen_sequences:
            if attempts > 100_000:
                raise RuntimeError("Could not generate a unique random population")
            continue
        seen_vectors.add(vector)
        seen_sequences.add(construct.sequence)
        generated.append(
            GeneratedCandidate(
                id=f"random-{random_index:03d}",
                vector=vector,
                sequence=construct.sequence,
                metadata={"anchor": "random", "seed": options.seed},
            )
        )
        random_index += 1


def _variable_positions_for_scope(
    design: ConstructDesignSpace,
    scope: ScorerScopeConfig,
) -> list[int]:
    positions: set[int] = set()
    cds_vector_start = 2
    cds_vector_end = cds_vector_start + design.codons.n_var

    def add_part(
        part: str,
        start: int | None = None,
        end: int | None = None,
    ) -> None:
        if part == "five_prime_utr":
            positions.add(0)
        elif part == "start_codon":
            positions.add(1)
        elif part == "cds":
            local_start = 0 if start is None else start
            local_end = len(design.codons.source_sequence) if end is None else end
            if local_end <= local_start:
                return
            first_codon = local_start // 3
            last_codon = (local_end - 1) // 3
            positions.update(
                range(
                    cds_vector_start + first_codon,
                    cds_vector_start + last_codon + 1,
                )
            )
        elif part == "stop_codon":
            positions.add(cds_vector_end)
        elif part == "three_prime_utr":
            positions.add(cds_vector_end + 1)
        elif part == "poly_a":
            positions.add(cds_vector_end + 2)

    source = design.source_construct
    for resolved in resolve_scope_ranges(source, scope):
        for part in PART_ORDER:
            part_start, part_end = source.coordinates[part]
            overlap_start = max(resolved.sequence_start, part_start)
            overlap_end = min(resolved.sequence_end, part_end)
            if overlap_end > overlap_start:
                add_part(
                    part,
                    overlap_start - part_start,
                    overlap_end - part_start,
                )

    lower = design.lower_bounds
    upper = design.upper_bounds
    return sorted(
        position
        for position in positions
        if 0 <= position < design.n_var and lower[position] < upper[position]
    )


def _select_diverse(
    ranked: list[_Proposal],
    count: int,
    min_hamming_distance: int,
) -> list[_Proposal]:
    selected: list[_Proposal] = []
    deferred: list[_Proposal] = []
    for proposal in ranked:
        if all(
            _hamming_distance(proposal.vector, other.vector)
            >= min_hamming_distance
            for other in selected
        ):
            selected.append(proposal)
        else:
            deferred.append(proposal)
        if len(selected) >= count:
            return selected
    for proposal in deferred:
        selected.append(proposal)
        if len(selected) >= count:
            break
    return selected


def _finite_objective(
    outputs: dict[str, Any],
    objective: OutputSpec,
    request: AnchorRequest,
) -> float:
    try:
        value = float(outputs[objective.name])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            f"Scorer {request.scorer!r} returned no numeric "
            f"{objective.name!r} value"
        ) from exc
    if not isfinite(value):
        raise ValueError(
            f"Scorer {request.scorer!r} returned non-finite "
            f"{objective.name!r}: {value}"
        )
    return value


def _desirability(value: float, objective: OutputSpec) -> float:
    return value if objective.direction == "max" else -value


def _hamming_distance(
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> int:
    return sum(a != b for a, b in zip(left, right, strict=True))


def _safe_name(value: str) -> str:
    return "".join(
        character if character.isalnum() or character in {"-", "_"} else "_"
        for character in value
    )


def _batched(
    values: Iterable[tuple[int, ...]],
    size: int,
) -> Iterable[list[tuple[int, ...]]]:
    iterator = iter(values)
    while batch := list(islice(iterator, size)):
        yield batch
