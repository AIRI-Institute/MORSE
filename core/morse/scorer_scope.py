from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any


PART_NAMES = {
    "five_prime_utr",
    "start_codon",
    "cds",
    "stop_codon",
    "three_prime_utr",
    "poly_a",
}

PART_ALIASES = {
    "full_cds": "full_cds",
    "five_prime_utr": "five_prime_utr",
    "5utr": "five_prime_utr",
    "5'utr": "five_prime_utr",
    "5′utr": "five_prime_utr",
    "utr5": "five_prime_utr",
    "start": "start_codon",
    "start_codon": "start_codon",
    "cds": "cds",
    "stop": "stop_codon",
    "stop_codon": "stop_codon",
    "three_prime_utr": "three_prime_utr",
    "3utr": "three_prime_utr",
    "3'utr": "three_prime_utr",
    "3′utr": "three_prime_utr",
    "utr3": "three_prime_utr",
    "poly_a": "poly_a",
    "polya": "poly_a",
}

_SCOPE_TERM = re.compile(
    r"^(?P<part>[A-Za-z0-9_\u2032']+)"
    r"(?:\[\s*(?P<start>-?\d*)\s*:\s*(?P<end>-?\d*)\s*\])?$"
)
_SCOPE_ANCHOR_EXPRESSION = re.compile(
    r"^(?P<part>[A-Za-z0-9_\u2032']+)\."
    r"(?P<edge>start|end)(?P<offset>[+-]\d+)?$"
)


@dataclass(frozen=True)
class ScopeSpan:
    part: str
    start: int | None = None
    end: int | None = None


@dataclass(frozen=True)
class ScopeAnchor:
    part: str
    edge: str
    offset: int = 0


@dataclass(frozen=True)
class ScorerScopeConfig:
    scope: str
    spans: tuple[ScopeSpan, ...] = ()
    start: ScopeAnchor | None = None
    end: ScopeAnchor | None = None


def parse_scorer_scope(raw: Any, context: str) -> ScorerScopeConfig:
    mapping = _mapping(raw, context)
    scope_type = str(mapping.get("type", "")).lower()
    if scope_type == "full_sequence":
        _reject_unknown_fields(mapping, {"type"}, context)
        return ScorerScopeConfig(scope="full_sequence")
    if scope_type == "window":
        _reject_unknown_fields(mapping, {"type", "start", "end"}, context)
        return ScorerScopeConfig(
            scope="window",
            start=_parse_anchor(mapping.get("start"), f"{context}.start"),
            end=_parse_anchor(mapping.get("end"), f"{context}.end"),
        )
    if scope_type != "parts":
        raise ValueError(
            f"{context}.type must be full_sequence, parts, or window"
        )

    _reject_unknown_fields(mapping, {"type", "spans"}, context)
    spans_raw = mapping.get("spans")
    if not isinstance(spans_raw, list) or not spans_raw:
        raise ValueError(f"{context}.spans must be a non-empty list")
    spans: list[ScopeSpan] = []
    for index, raw_span in enumerate(spans_raw):
        span_context = f"{context}.spans[{index}]"
        span = _mapping(raw_span, span_context)
        _reject_unknown_fields(span, {"part", "range"}, span_context)
        part = str(span.get("part", ""))
        if part not in PART_NAMES | {"full_cds"}:
            raise ValueError(
                f"{span_context}.part must be one of: "
                + ", ".join(sorted(PART_NAMES))
            )
        local_range = (
            None
            if span.get("range") is None
            else _parse_range(span["range"], f"{span_context}.range")
        )
        spans.append(
            ScopeSpan(
                part=part,
                start=None if local_range is None else local_range[0],
                end=None if local_range is None else local_range[1],
            )
        )
    return ScorerScopeConfig(scope="parts", spans=tuple(spans))


def scorer_scope_to_dict(config: ScorerScopeConfig) -> dict[str, Any]:
    if config.scope == "full_sequence":
        return {"type": "full_sequence"}
    if config.scope == "parts":
        return {
            "type": "parts",
            "spans": [
                {
                    "part": span.part,
                    **(
                        {"range": [span.start, span.end]}
                        if span.start is not None or span.end is not None
                        else {}
                    ),
                }
                for span in config.spans
            ],
        }
    if config.scope == "window":
        if config.start is None or config.end is None:
            raise ValueError("Window scope requires start and end anchors")
        return {
            "type": "window",
            "start": _anchor_to_dict(config.start),
            "end": _anchor_to_dict(config.end),
        }
    raise ValueError(f"Unsupported scorer scope: {config.scope}")


def normalize_scorer_scope(raw: Any, context: str) -> dict[str, Any]:
    return scorer_scope_to_dict(parse_scorer_scope(raw, context))


def parse_scope_expression(expression: str) -> dict[str, Any]:
    """Compile the compact notebook scope expression to the YAML contract."""

    text = str(expression).strip().lower()
    if not text:
        raise ValueError("Scorer scope expression must not be empty")
    if text == "full_sequence":
        return {"type": "full_sequence"}
    if text.startswith("window("):
        return scorer_scope_to_dict(_parse_window_expression(text))

    spans: list[ScopeSpan] = []
    for index, raw_term in enumerate(text.split("+")):
        term = re.sub(r"\s+", "", raw_term)
        if not term:
            raise ValueError(
                "Invalid scorer scope expression: `+` must separate two parts"
            )
        match = _SCOPE_TERM.fullmatch(term)
        if match is None:
            raise ValueError(
                f"Invalid scorer scope term #{index + 1}: {raw_term.strip()!r}. "
                "Use forms such as `cds`, `cds[:30]`, or `cds[-30:]`"
            )
        raw_part = match.group("part")
        if raw_part.isdigit():
            raise ValueError(
                "A number cannot be concatenated to a scope part. "
                "Use a slice such as `cds[:30]` for the first 30 nt or "
                "`window(cds, downstream=30)` for downstream context"
            )
        part = _normalize_part_name(raw_part, f"term #{index + 1}")
        start_text = match.group("start")
        end_text = match.group("end")
        if start_text is None and end_text is None:
            start = end = None
        else:
            start = int(start_text) if start_text else 0
            end = int(end_text) if end_text else None
            if start == 0 and end is None:
                start = None
            _validate_range_shape(start, end, f"term #{index + 1}")
        spans.append(ScopeSpan(part=part, start=start, end=end))
    return scorer_scope_to_dict(
        ScorerScopeConfig(scope="parts", spans=tuple(spans))
    )


def scorer_scope_expression(raw: Any, context: str = "scope") -> str:
    """Render a normalized scope using the compact notebook syntax."""

    config = parse_scorer_scope(raw, context)
    if config.scope == "full_sequence":
        return "full_sequence"
    if config.scope == "parts":
        terms = []
        for span in config.spans:
            term = span.part
            if span.start is not None or span.end is not None:
                start = "" if span.start in {None, 0} else str(span.start)
                end = "" if span.end is None else str(span.end)
                term += f"[{start}:{end}]"
            terms.append(term)
        return " + ".join(terms)
    if config.scope == "window":
        if config.start is None or config.end is None:
            raise ValueError("Window scope requires start and end anchors")
        if (
            config.start.part == config.end.part
            and config.start.edge == "start"
            and config.end.edge == "end"
            and config.start.offset <= 0
            and config.end.offset >= 0
        ):
            arguments = [config.start.part]
            if config.start.offset:
                arguments.append(f"upstream={-config.start.offset}")
            if config.end.offset:
                arguments.append(f"downstream={config.end.offset}")
            return f"window({', '.join(arguments)})"
        return (
            "window("
            + _anchor_expression(config.start)
            + ", "
            + _anchor_expression(config.end)
            + ")"
        )
    raise ValueError(f"Unsupported scorer scope: {config.scope}")


def scope_parts(raw: Any, context: str = "scope") -> tuple[str, ...]:
    """Return the parts explicitly selected or used as window anchors."""

    config = parse_scorer_scope(raw, context)
    if config.scope == "full_sequence":
        return tuple(sorted(PART_NAMES))
    if config.scope == "parts":
        return tuple(dict.fromkeys(span.part for span in config.spans))
    if config.scope == "window":
        if config.start is None or config.end is None:
            raise ValueError("Window scope requires start and end anchors")
        return tuple(dict.fromkeys((config.start.part, config.end.part)))
    raise ValueError(f"Unsupported scorer scope: {config.scope}")


def _parse_range(raw: Any, context: str) -> tuple[int | None, int | None]:
    if not isinstance(raw, list) or len(raw) != 2:
        raise ValueError(f"{context} must be a two-item [start, end] list")
    parsed: list[int | None] = []
    for value in raw:
        if value is None:
            parsed.append(None)
            continue
        if isinstance(value, bool):
            raise ValueError(f"{context} values must be integers or null")
        try:
            parsed.append(int(value))
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{context} values must be integers or null"
            ) from exc
    start, end = parsed
    if start is None and end is not None:
        start = 0
    if start == 0 and end is None:
        start = None
    _validate_range_shape(start, end, context)
    return start, end


def _validate_range_shape(
    start: int | None,
    end: int | None,
    context: str,
) -> None:
    if start is None and end is None:
        return
    if start is not None and end is not None:
        if start >= 0 and end >= 0 and end <= start:
            raise ValueError(f"{context} must satisfy start < end")
        if start < 0 and end < 0 and end <= start:
            raise ValueError(f"{context} must satisfy start < end")


def _parse_anchor(raw: Any, context: str) -> ScopeAnchor:
    mapping = _mapping(raw, context)
    _reject_unknown_fields(mapping, {"part", "edge", "offset"}, context)
    part = _normalize_part_name(str(mapping.get("part", "")), context)
    edge = str(mapping.get("edge", "")).strip().lower()
    if edge not in {"start", "end"}:
        raise ValueError(f"{context}.edge must be start or end")
    offset_raw = mapping.get("offset", 0)
    if isinstance(offset_raw, bool):
        raise ValueError(f"{context}.offset must be an integer")
    try:
        offset = int(offset_raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context}.offset must be an integer") from exc
    return ScopeAnchor(part=part, edge=edge, offset=offset)


def _anchor_to_dict(anchor: ScopeAnchor) -> dict[str, Any]:
    return {
        "part": anchor.part,
        "edge": anchor.edge,
        **({"offset": anchor.offset} if anchor.offset else {}),
    }


def _parse_window_expression(text: str) -> ScorerScopeConfig:
    if not text.endswith(")"):
        raise ValueError(
            "Invalid window scope. Use `window(cds, upstream=30, downstream=30)`"
        )
    arguments = [item.strip() for item in text[7:-1].split(",")]
    if not arguments or not arguments[0]:
        raise ValueError("Window scope requires a part name")
    if len(arguments) == 2 and all("=" not in item for item in arguments):
        return ScorerScopeConfig(
            scope="window",
            start=_parse_anchor_expression(arguments[0], "window start"),
            end=_parse_anchor_expression(arguments[1], "window end"),
        )
    part = _normalize_part_name(arguments[0], "window")
    values = {"upstream": 0, "downstream": 0}
    seen: set[str] = set()
    for argument in arguments[1:]:
        name, separator, raw_value = argument.partition("=")
        name = name.strip()
        if not separator or name not in values or not raw_value.strip().isdigit():
            raise ValueError(
                "Window options must be non-negative `upstream=N` or "
                "`downstream=N` values"
            )
        if name in seen:
            raise ValueError(f"Window option `{name}` is repeated")
        seen.add(name)
        values[name] = int(raw_value)
    return ScorerScopeConfig(
        scope="window",
        start=ScopeAnchor(part, "start", -values["upstream"]),
        end=ScopeAnchor(part, "end", values["downstream"]),
    )


def _normalize_part_name(value: str, context: str) -> str:
    part = PART_ALIASES.get(value.strip().lower())
    if part is None:
        raise ValueError(
            f"Unknown sequence part in {context}: {value!r}. Expected one of: "
            + ", ".join(sorted(PART_NAMES))
        )
    return part


def _parse_anchor_expression(value: str, context: str) -> ScopeAnchor:
    match = _SCOPE_ANCHOR_EXPRESSION.fullmatch(re.sub(r"\s+", "", value))
    if match is None:
        raise ValueError(
            f"Invalid {context}: {value!r}. "
            "Use an anchor such as `cds.start-30` or `cds.end+30`"
        )
    return ScopeAnchor(
        part=_normalize_part_name(match.group("part"), context),
        edge=match.group("edge"),
        offset=int(match.group("offset") or 0),
    )


def _anchor_expression(anchor: ScopeAnchor) -> str:
    offset = "" if not anchor.offset else f"{anchor.offset:+d}"
    return f"{anchor.part}.{anchor.edge}{offset}"


def _mapping(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"`{context}` must be a mapping")
    return dict(value)


def _reject_unknown_fields(
    mapping: Mapping[str, Any],
    allowed: set[str],
    context: str,
) -> None:
    unknown = sorted(set(mapping) - allowed)
    if unknown:
        raise ValueError(f"Unknown {context} fields: " + ", ".join(unknown))
