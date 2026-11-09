from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from morse.construct import PART_ORDER, VaccineConstruct
from morse.scorer_scope import ScopeAnchor, ScorerScopeConfig


@dataclass(frozen=True)
class ResolvedScopeSpan:
    part: str | None
    sequence_start: int
    sequence_end: int


@dataclass(frozen=True)
class RoutedPartRange:
    part: str
    routed_start: int
    routed_end: int
    part_start: int
    part_end: int

    def to_dict(self) -> dict[str, int | str]:
        return {
            "part": self.part,
            "start": self.routed_start,
            "end": self.routed_end,
            "part_start": self.part_start,
            "part_end": self.part_end,
        }


@dataclass(frozen=True)
class RoutedSequence:
    sequence: str
    part_ranges: tuple[RoutedPartRange, ...]

    def context(self) -> dict[str, object]:
        return {"part_ranges": [item.to_dict() for item in self.part_ranges]}


class ScorerInputRouter:
    """Resolve a scorer scope inside Core and return only its nucleotide sequence."""

    def route(
        self,
        construct: VaccineConstruct,
        config: ScorerScopeConfig,
    ) -> str:
        return self.route_with_context(construct, config).sequence

    def route_with_context(
        self,
        construct: VaccineConstruct,
        config: ScorerScopeConfig,
    ) -> RoutedSequence:
        resolved = self.resolve_scope_ranges(construct, config)
        chunks: list[str] = []
        part_ranges: list[RoutedPartRange] = []
        routed_cursor = 0
        for span in resolved:
            chunk = construct.sequence[span.sequence_start : span.sequence_end]
            chunks.append(chunk)
            for part in PART_ORDER:
                part_start, part_end = construct.coordinates[part]
                overlap_start = max(span.sequence_start, part_start)
                overlap_end = min(span.sequence_end, part_end)
                if overlap_end <= overlap_start:
                    continue
                offset = overlap_start - span.sequence_start
                part_ranges.append(
                    RoutedPartRange(
                        part=part,
                        routed_start=routed_cursor + offset,
                        routed_end=routed_cursor + offset + overlap_end - overlap_start,
                        part_start=overlap_start - part_start,
                        part_end=overlap_end - part_start,
                    )
                )
            routed_cursor += len(chunk)
        return RoutedSequence("".join(chunks), tuple(part_ranges))

    def resolve_scope_ranges(
        self,
        construct: VaccineConstruct,
        config: ScorerScopeConfig,
    ) -> tuple[ResolvedScopeSpan, ...]:
        scope = config.scope
        coordinates = _scope_coordinates(construct)
        if scope == "full_sequence":
            return (
                ResolvedScopeSpan(
                    None,
                    0,
                    len(construct.sequence),
                ),
            )
        if scope == "parts":
            resolved: list[ResolvedScopeSpan] = []
            for span in config.spans:
                segment_start, segment_end = coordinates[span.part]
                segment_length = segment_end - segment_start
                local_start = _resolve_part_index(
                    span.start,
                    segment_length,
                    default=0,
                )
                local_end = _resolve_part_index(
                    span.end,
                    segment_length,
                    default=segment_length,
                )
                if (
                    local_start < 0
                    or local_end <= local_start
                    or local_end > segment_length
                ):
                    raise ValueError(
                        f"scorer scope {span.part}[{local_start}, {local_end}) exceeds "
                        f"part length {segment_length}"
                    )
                resolved.append(
                    ResolvedScopeSpan(
                        part=span.part,
                        sequence_start=segment_start + local_start,
                        sequence_end=segment_start + local_end,
                    )
                )
            _validate_resolved_order(resolved)
            return tuple(resolved)
        if scope == "window":
            if config.start is None or config.end is None:
                raise ValueError("Window scorer scope requires start and end anchors")
            sequence_start = _resolve_anchor(coordinates, config.start)
            sequence_end = _resolve_anchor(coordinates, config.end)
            if (
                sequence_start < 0
                or sequence_end <= sequence_start
                or sequence_end > len(construct.sequence)
            ):
                raise ValueError(
                    "scorer window "
                    f"[{sequence_start}, {sequence_end}) exceeds assembled "
                    f"sequence length {len(construct.sequence)}"
                )
            return (
                ResolvedScopeSpan(
                    part=None,
                    sequence_start=sequence_start,
                    sequence_end=sequence_end,
                ),
            )
        raise ValueError(f"Unsupported scorer input scope: {scope}")


DEFAULT_SCORER_INPUT_ROUTER = ScorerInputRouter()


def _scope_coordinates(construct: VaccineConstruct) -> dict[str, tuple[int, int]]:
    coordinates = dict(construct.coordinates)
    coordinates["full_cds"] = (
        coordinates["start_codon"][0], coordinates["stop_codon"][1]
    )
    return coordinates


def resolve_scope_ranges(
    construct: VaccineConstruct,
    config: ScorerScopeConfig,
) -> tuple[ResolvedScopeSpan, ...]:
    return DEFAULT_SCORER_INPUT_ROUTER.resolve_scope_ranges(construct, config)


def _resolve_part_index(
    value: int | None,
    length: int,
    *,
    default: int,
) -> int:
    if value is None:
        return default
    return length + value if value < 0 else value


def _resolve_anchor(
    coordinates: Mapping[str, tuple[int, int]],
    anchor: ScopeAnchor,
) -> int:
    segment_start, segment_end = coordinates[anchor.part]
    base = segment_start if anchor.edge == "start" else segment_end
    return base + anchor.offset


def _validate_resolved_order(spans: list[ResolvedScopeSpan]) -> None:
    for previous, current in zip(spans, spans[1:]):
        if current.sequence_start < previous.sequence_end:
            raise ValueError(
                "scorer scope spans must follow construct order without overlap: "
                f"[{previous.sequence_start}, {previous.sequence_end}) then "
                f"[{current.sequence_start}, {current.sequence_end})"
            )
