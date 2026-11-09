from __future__ import annotations

from dataclasses import dataclass
from html import escape
from math import prod
import json
from pathlib import Path

from morse.construct import PART_ORDER
from morse.run_spec import PartSourceSpec, RunSpec
from morse.router import resolve_scope_ranges
from morse.scorer_scope import ScorerScopeConfig
from morse.publication_style import (
    INK, MUTED, GRID, FIXED_FILL, MUTABLE_FILL, MUTABLE_EDGE, palette_color,
)


_PART_LABELS = {
    "full_cds": "CDS (full)",
    "five_prime_utr": "5′UTR",
    "start_codon": "START",
    "cds": "CDS BODY",
    "stop_codon": "STOP",
    "three_prime_utr": "3′UTR",
    "poly_a": "poly(A)",
}


@dataclass(frozen=True)
class OptimizerPartSummary:
    part: str
    variant_count: int
    min_length: int
    max_length: int
    source: PartSourceSpec
    mutable_positions: int | None = None


def write_sequence_scope_svg(path: Path, spec: RunSpec) -> Path:
    path.write_text(render_sequence_scope_svg(spec))
    return path


def render_sequence_scope_svg(spec: RunSpec) -> str:
    """Publication-style reference coordinates and actual routed scorer inputs.

    The nucleotide axis is proportional, including three-base codons. Narrow
    labels use separate callout lanes rather than enlarging a sequence segment.
    """
    design, run_config = spec.design_space, spec.run_config
    construct = design.source_construct
    length = len(construct.sequence)
    visible_parts = tuple(
        part for part in PART_ORDER
        if construct.coordinates[part][1] > construct.coordinates[part][0]
    )
    optimizer_parts = optimizer_part_summaries(spec)
    summaries = {item.part: item for item in optimizer_parts}
    sources = {source.part: source for source in spec.part_sources}
    total_variants = prod(item.variant_count for item in optimizer_parts)
    optimizer_name = {"ga": "GA", "nsga2": "NSGA-II", "nsga3": "NSGA-III"}[run_config.optimizer.type]
    minimum, maximum = design.candidate_length_range
    labels = {
        "five_prime_utr": "5′UTR", "start_codon": construct.start_codon.replace("T", "U"),
        "cds": "CDS", "stop_codon": "stop " + construct.stop_codon.replace("T", "U"),
        "three_prime_utr": "3′UTR", "poly_a": f"A({construct.poly_a_length})",
    }
    width, left, right = 1280, 364, 1244
    track_width, track_height = right - left, 46

    def x_position(position):
        return left + track_width * position / length

    def wrapped(text, columns=43, font_size=11):
        # Approximate Arial advance widths, including wide uppercase names and
        # unbroken identifiers. All original text remains in the wrapped rows.
        def advance(char):
            if char in "W@%":
                return font_size
            if char in "Mwm":
                return font_size * .86
            if char in "ilI.,:;'′|! ":
                return font_size * .30
            return font_size * (.70 if char.isupper() else .59)
        limit, rows, current, size = columns * 6.2, [], "", 0.0
        for word in str(text).split():
            word_size = sum(map(advance, word))
            if current and size + advance(" ") + word_size > limit:
                rows.append(current)
                current, size = "", 0.0
            for char in ((" " if current else "") + word):
                char_size = advance(char)
                if current and size + char_size > limit:
                    rows.append(current)
                    current, size = "", 0.0
                current += char
                size += char_size
        return rows + ([current] if current else []) or [""]

    title_rows = wrapped(design.name, 188, font_size=22)
    header_bottom = 33 + (len(title_rows) - 1) * 26
    # Label the narrow parts above the true-scale track, with collision-free
    # lanes. Wide parts retain their names and lengths inside the block.
    callouts, occupied = [], []
    for part in visible_parts:
        start, end = construct.coordinates[part]
        part_width = x_position(end) - x_position(start)
        if part_width >= 104:
            continue
        text = f"{labels[part]} · {end - start} nt"
        label_width = max(72, len(text) * 6.7 + 12)
        center = min(right - label_width / 2, max(left + label_width / 2, x_position((start + end) / 2)))
        interval = (center - label_width / 2 - 8, center + label_width / 2 + 8)
        lane = next((i for i, spans in enumerate(occupied) if all(interval[1] <= a or interval[0] >= b for a, b in spans)), len(occupied))
        if lane == len(occupied):
            occupied.append([])
        occupied[lane].append(interval)
        callouts.append((part, center, lane, text))
    sequence_y = header_bottom + 96 + max(0, len(occupied) - 1) * 25
    scorer_header_y = sequence_y + track_height + 77
    y = scorer_header_y + 20
    scorer_rows = []
    for index, scorer in enumerate(run_config.scorers):
        spans = resolve_scope_ranges(construct, scorer.scope)
        ranges = _merge_ranges([(span.sequence_start, span.sequence_end) for span in spans])
        names = wrapped(scorer.name, 47, font_size=13)
        description = scope_label(scorer.scope).replace("CDS BODY", "CDS").replace("START", "start").replace("STOP", "stop").replace("poly(A)", "A(n)")
        detail = wrapped(description, 45)
        routed = sum(end - start for start, end in ranges)
        bounds = wrapped(" + ".join(f"[{start}, {end})" for start, end in ranges) + f" · {routed} nt", 45)
        row_height = max(58, 17 * len(names) + 14 * (len(detail) + len(bounds)) + 15)
        scorer_rows.append((index, scorer, ranges, names, detail, bounds, y, row_height))
        y += row_height
    table_header_y = y + 34
    table_y = table_header_y + 48
    part_rows = []
    for part in visible_parts:
        start, end = construct.coordinates[part]
        item = summaries.get(part)
        if item:
            details = list(_source_details(item))
            choices = _variant_label(item.variant_count)
        else:
            source = sources[part]
            details = [f"{source.kind.upper() if source.kind in {'tsv', 'fasta'} else source.kind} · {_display_path(source.path)}"]
            if source.set_name:
                details.append(f"set: {source.set_name}")
            choices = "1 reference variant"
        source_rows = [line for detail in details for line in wrapped(detail, 60)]
        row_height = max(48, 14 * len(source_rows) + 18)
        part_rows.append((part, item, start, end, choices, source_rows, table_y, row_height))
        table_y += row_height
    height = table_y + 52
    metadata = {
        "reference_sequence_dna": construct.sequence,
        "coordinate_convention": "0-based half-open reference coordinates; track proportional to nucleotide length",
        "reference_length_nt": length, "candidate_length_range_nt": [minimum, maximum],
        "mutable_variable_count": design.mutable_variable_count, "search_space_candidates": str(total_variants),
        "parts": [{"part": part, "start": start, "end": end, "mutable": part in design.mutable_parts,
                   "variant_count": str(item.variant_count if item else 1), "source": _display_path(sources[part].path)}
                  for part, item, start, end, *_ in part_rows],
        "scorers": [{"name": scorer.name, "color": palette_color(index),
                     "ranges": [{"start": start, "end": end} for start, end in ranges]}
                    for index, scorer, ranges, *_ in scorer_rows],
    }
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title description">',
        '<title id="title">mRNA design space and scorer input scopes</title>',
        '<desc id="description">Reference construct with fixed and variable parts, exact nucleotide boundaries, search-space counts and routed scorer inputs. RNA codons are displayed; the reference sequence is preserved as DNA in metadata.</desc>',
        '<metadata id="morse-scope-metadata">' + escape(json.dumps(metadata, ensure_ascii=False)) + '</metadata>',
        '<style>',
        f'text{{font-family:Arial,Helvetica,sans-serif;fill:{INK}}}',
        f'.title{{font-size:22px;font-weight:600}}.subtitle{{font-size:13px;fill:{MUTED}}}',
        f'.section{{font-size:12px;font-weight:600;letter-spacing:1px;fill:{INK}}}',
        f'.label{{font-size:13px;font-weight:500}}.detail{{font-size:11px;fill:{MUTED}}}',
        f'.column{{font-size:10px;font-weight:600;fill:{MUTED}}}',
        f'.inside{{font-size:13px;text-anchor:middle}}.axis{{font-size:10px;fill:{MUTED}}}',
        f'.guide{{stroke:{GRID};stroke-width:0.8}}.callout{{fill:none;stroke:{MUTED};stroke-width:0.8}}',
        '.scope{fill:none;stroke-width:1.8;stroke-linecap:butt}',
        '</style>',
        f'<rect class="publication-background" width="{width}" height="{height}" fill="#FFFFFF"/>',
    ]
    for i, title in enumerate(title_rows):
        lines.append(f'<text class="title" x="28" y="{33 + i * 26}">{escape(title)}</text>')
    lines.append(f'<text class="subtitle" x="28" y="{header_bottom + 24}">Reference {length:,} nt · candidate range {minimum:,}–{maximum:,} nt · {design.mutable_variable_count} mutable variables</text>')
    legend_y = header_bottom + 40
    for x, fill, edge, label in [(28, FIXED_FILL, GRID, "Fixed part"), (158, MUTABLE_FILL, MUTABLE_EDGE, "Variable part (optimizer)")]:
        lines.extend([f'<rect x="{x}" y="{legend_y}" width="18" height="12" rx="2" fill="{fill}" stroke="{edge}"/>', f'<text class="detail" x="{x + 26}" y="{legend_y + 10}">{label}</text>'])
    lines.append(f'<text class="section" x="28" y="{sequence_y + 18}">REFERENCE mRNA</text>')
    lines.append(f'<text class="detail" x="28" y="{sequence_y + 36}">RNA codon labels · proportional nt axis</text>')
    for part in visible_parts:
        start, end = construct.coordinates[part]
        x = x_position(start)
        part_width = x_position(end) - x
        mutable = part in design.mutable_parts
        fill, edge = (MUTABLE_FILL, MUTABLE_EDGE) if mutable else (FIXED_FILL, GRID)
        if mutable and part in {"five_prime_utr", "three_prime_utr", "poly_a"} and summaries[part].variant_count > 1 and part_width >= 24:
            for offset in (6, 3):
                lines.append(f'<rect class="library-stack" x="{x + offset:.2f}" y="{sequence_y - offset:.2f}" width="{part_width - offset:.2f}" height="{track_height}" rx="3" fill="{fill}" stroke="{edge}" stroke-width="0.7" aria-hidden="true"/>')
        lines.append(f'<rect class="sequence-part" data-part="{part}" data-start="{start}" data-end="{end}" data-mutable="{str(mutable).lower()}" x="{x:.4f}" y="{sequence_y}" width="{part_width:.4f}" height="{track_height}" fill="{fill}" stroke="{edge}" stroke-width="{1.2 if mutable else 0.8}"><title>{escape(labels[part])}: [{start}, {end}), {end - start} nt; {"variable" if mutable else "fixed"}</title></rect>')
        if part_width >= 104:
            center = x + part_width / 2
            lines.extend([f'<text class="inside" x="{center:.2f}" y="{sequence_y + 20}">{escape(labels[part])}</text>', f'<text class="detail" text-anchor="middle" x="{center:.2f}" y="{sequence_y + 36}">{end - start} nt</text>'])
    for part, center, lane, text in callouts:
        start, end = construct.coordinates[part]
        actual = x_position((start + end) / 2)
        label_y = sequence_y - 18 - lane * 25
        lines.extend([f'<path class="callout" d="M {actual:.2f},{sequence_y - 1} L {actual:.2f},{label_y + 7} L {center:.2f},{label_y + 7}"/>', f'<text class="detail" text-anchor="middle" x="{center:.2f}" y="{label_y}">{escape(text)}</text>'])
    axis_y = sequence_y + track_height + 16
    for position in sorted({0, length, *(value for values in construct.coordinates.values() for value in values)}):
        x = x_position(position)
        lines.append(f'<line class="guide" x1="{x:.4f}" x2="{x:.4f}" y1="{sequence_y + track_height + 2}" y2="{axis_y - 6}"/>')
    lines.extend([f'<text class="axis" x="{left}" y="{axis_y + 7}">0</text>', f'<text class="axis" text-anchor="end" x="{right}" y="{axis_y + 7}">{length:,} nt</text>', f'<text class="detail" text-anchor="middle" x="{(left + right) / 2}" y="{axis_y + 7}">0-based half-open coordinates; exact boundaries below</text>'])
    lines.append(f'<text class="section" x="28" y="{scorer_header_y}">SCORER INPUT SCOPES</text>')
    for index, scorer, ranges, names, details, bounds, row_y, row_height in scorer_rows:
        color = palette_color(index)
        lane_y = row_y + row_height / 2
        lines.append(f'<g class="scorer-row" data-scorer-name="{escape(scorer.name, quote=True)}" data-scorer-index="{index}" data-color="{color}">')
        lines.append(f'<line x1="28" x2="38" y1="{row_y + 13}" y2="{row_y + 13}" stroke="{color}" stroke-width="2"/>')
        text_y = row_y + 17
        for name in names:
            lines.append(f'<text class="label" x="46" y="{text_y}">{escape(name)}</text>')
            text_y += 17
        for detail in details + bounds:
            lines.append(f'<text class="detail" x="46" y="{text_y}">{escape(detail)}</text>')
            text_y += 14
        lines.append(f'<line class="guide" x1="{left}" x2="{right}" y1="{lane_y:.2f}" y2="{lane_y:.2f}"/>')
        for start, end in ranges:
            x1, x2 = x_position(start), x_position(end)
            lines.append(f'<g class="routed-scope" data-start="{start}" data-end="{end}"><line class="scope" x1="{x1:.4f}" x2="{x2:.4f}" y1="{lane_y:.2f}" y2="{lane_y:.2f}" stroke="{color}"/><line class="scope" x1="{x1:.4f}" x2="{x1:.4f}" y1="{lane_y - 5:.2f}" y2="{lane_y + 5:.2f}" stroke="{color}"/><line class="scope" x1="{x2:.4f}" x2="{x2:.4f}" y1="{lane_y - 5:.2f}" y2="{lane_y + 5:.2f}" stroke="{color}"/></g>')
        lines.append('</g>')
    lines.extend([f'<text class="section" x="28" y="{table_header_y}">SEARCH SPACE AND PART SOURCES</text>', f'<text class="detail" x="28" y="{table_header_y + 18}">{optimizer_name} · ≈{_format_count(total_variants)} candidate combinations · reference coordinates below</text>'])
    for x, label in [(28, "Part"), (158, "Editable"), (245, "Variants & length"), (476, "Source")]:
        lines.append(f'<text class="column" x="{x}" y="{table_header_y + 39}">{escape(label)}</text>')
    lines.append(f'<text class="column" text-anchor="end" x="{right}" y="{table_header_y + 39}">Reference interval</text>')
    for part, item, start, end, choices, details, row_y, row_height in part_rows:
        mutable = part in design.mutable_parts
        color = MUTABLE_EDGE if mutable else MUTED
        lines.append(f'<line class="guide" x1="28" x2="{right}" y1="{row_y}" y2="{row_y}"/>')
        lines.append(f'<text class="label" x="28" y="{row_y + 20}">{escape(labels[part])}</text>')
        lines.append(f'<text class="detail" x="158" y="{row_y + 20}" style="fill:{color}">{"Variable" if mutable else "Fixed"}</text>')
        lines.append(f'<text class="detail" x="245" y="{row_y + 20}">{escape(choices)}</text>')
        if item:
            lines.append(f'<text class="detail" x="245" y="{row_y + 35}">{_length_label(item.min_length, item.max_length)}</text>')
        for i, detail in enumerate(details):
            lines.append(f'<text class="detail" x="476" y="{row_y + 20 + i * 14}">{escape(detail)}</text>')
        lines.append(f'<text class="detail" text-anchor="end" x="{right}" y="{row_y + 20}">[{start}, {end}) · {end - start} nt</text>')
    lines.append(f'<text class="detail" x="28" y="{height - 18}">Scorer lines show the configured model inputs. Intervals use [start, end) reference nucleotide coordinates.</text>')
    lines.append('</svg>\n')
    return "\n".join(lines)


def optimizer_part_summaries(spec: RunSpec) -> tuple[OptimizerPartSummary, ...]:
    design = spec.design_space
    source_by_part = {source.part: source for source in spec.part_sources}
    lengths_by_part = {
        "five_prime_utr": [len(item.sequence) for item in design.five_prime_utrs],
        "start_codon": [len(item) for item in design.start_codons],
        "cds": [len(design.codons.source_sequence)],
        "stop_codon": [len(item) for item in design.stop_codons],
        "three_prime_utr": [len(item.sequence) for item in design.three_prime_utrs],
        "poly_a": list(design.poly_a_lengths),
    }
    counts = {
        "five_prime_utr": len(design.five_prime_utrs),
        "start_codon": len(design.start_codons),
        "cds": prod(len(choices) for choices in design.codons.choices_by_position),
        "stop_codon": len(design.stop_codons),
        "three_prime_utr": len(design.three_prime_utrs),
        "poly_a": len(design.poly_a_lengths),
    }
    return tuple(
        OptimizerPartSummary(
            part=part,
            variant_count=counts[part],
            min_length=min(lengths_by_part[part]),
            max_length=max(lengths_by_part[part]),
            source=source_by_part[part],
            mutable_positions=(
                sum(len(choices) > 1 for choices in design.codons.choices_by_position)
                if part == "cds"
                else None
            ),
        )
        for part in PART_ORDER
        if part in design.mutable_parts
    )


def scope_label(config: ScorerScopeConfig) -> str:
    if config.scope == "full_sequence":
        return "full sequence"
    if config.scope == "parts":
        labels = []
        for span in config.spans:
            label = _PART_LABELS[span.part]
            if span.start is not None or span.end is not None:
                start = "" if span.start in {None, 0} else str(span.start)
                end = "" if span.end is None else str(span.end)
                label += f"[{start}:{end}]"
            labels.append(label)
        return " + ".join(labels)
    if config.scope == "window":
        if config.start is None or config.end is None:
            return "invalid window"
        start_offset = (
            "" if not config.start.offset else f"{config.start.offset:+d}"
        )
        end_offset = "" if not config.end.offset else f"{config.end.offset:+d}"
        return (
            f"{_PART_LABELS[config.start.part]}.{config.start.edge}{start_offset}"
            " → "
            f"{_PART_LABELS[config.end.part]}.{config.end.edge}{end_offset}"
        )
    return config.scope


def _source_details(summary: OptimizerPartSummary) -> tuple[str, ...]:
    source = summary.source
    if summary.part == "cds":
        return (
            f"synonymous codons · {summary.mutable_positions} mutable positions",
            f"FASTA · {_display_path(source.path)}"
            if source.kind == "fasta"
            else f"reference CDS body · {source.path.name}",
        )
    if source.kind == "tsv":
        selection = f"set: {source.set_name}" if source.set_name else "all matching rows"
        if source.max_candidates is not None:
            selection += f" · limit: {source.max_candidates}"
        return (f"TSV · {_display_path(source.path)}", selection)
    source_label = {
        "inline": "inline variants",
        "fixed": "fixed sequence",
        "lengths": "length choices",
    }.get(source.kind, source.kind)
    return (f"{source_label} · {source.path.name}",)


def _display_path(path: Path) -> str:
    parts = path.parts
    for marker in ("examples", "configs", "data", "results", "scorers"):
        if marker in parts:
            return "/".join(parts[parts.index(marker) :])
    return path.name


def _variant_label(count: int) -> str:
    if count >= 1_000_000:
        return f"≈{_format_count(count)} variants"
    return f"{count} {'variant' if count == 1 else 'variants'}"


def _format_count(count: int) -> str:
    if count < 1_000_000:
        return f"{count:,}"
    digits = str(count)
    leading = int(digits[:3]) + (int(digits[3]) >= 5)
    exponent = len(digits) - 1
    if leading == 1000:
        leading = 100
        exponent += 1
    significant = str(leading)
    decimals = significant[1:].rstrip("0")
    coefficient = significant[0] + (f".{decimals}" if decimals else "")
    return f"{coefficient}×10^{exponent}"


def _length_label(min_length: int, max_length: int) -> str:
    return f"{min_length} nt" if min_length == max_length else f"{min_length}–{max_length} nt"


def _merge_ranges(ranges: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not ranges:
        return []
    merged: list[list[int]] = []
    for start, end in sorted(ranges):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]
