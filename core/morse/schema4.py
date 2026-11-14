"""Public schema-4 conventions translated to the established engine parts."""
from __future__ import annotations

import copy
import re
from dataclasses import replace
from pathlib import Path
from typing import Any

from morse.construct import (
    PART_ORDER, ConstructDesignSpace, SequenceVariant,
    construct_design_from_config,
)
from morse.codon import normalize_dna
from morse.parts import load_part_library
from morse.scorer_scope import parse_scope_expression, scorer_scope_expression


def _public_expression_to_internal(text: str) -> str:
    text = text.lower()
    text = re.sub(r"\bcds_body\b", "__body__", text)
    text = re.sub(r"\bcds\b", "full_cds", text)
    return text.replace("__body__", "cds")


def scope_to_internal(raw: Any, *, stop_present: bool = True) -> dict[str, Any]:
    """Translate public CDS/body names, without changing schema-3 semantics."""
    del stop_present  # full_cds is a logical range and handles an absent stop itself.
    if isinstance(raw, str):
        return parse_scope_expression(_public_expression_to_internal(raw))
    if isinstance(raw, list):
        if not raw:
            raise ValueError("scorer.scope cannot be empty")
        spans = []
        for item in raw:
            if isinstance(item, str):
                parsed = scope_to_internal(item)
                if parsed["type"] != "parts":
                    raise ValueError("A scope list must contain parts or part slices")
                spans.extend(parsed["spans"])
            elif isinstance(item, dict):
                spans.append(_map_scope_names(item, to_internal=True))
            else:
                raise ValueError("A scope list must contain part names or ranges")
        return {"type": "parts", "spans": spans}
    if isinstance(raw, dict):
        return _map_scope_names(copy.deepcopy(raw), to_internal=True)
    raise ValueError("scorer.scope must be a part, expression, list or scope mapping")


def _map_scope_names(value: Any, *, to_internal: bool) -> Any:
    if isinstance(value, list):
        return [_map_scope_names(item, to_internal=to_internal) for item in value]
    if not isinstance(value, dict):
        return value
    result = {}
    mapping = {"cds": "full_cds", "cds_body": "cds"} if to_internal else {"cds": "cds_body", "full_cds": "cds"}
    for key, item in value.items():
        result[key] = mapping.get(item, item) if key == "part" and isinstance(item, str) else _map_scope_names(item, to_internal=to_internal)
    return result


def scope_from_internal(raw: dict[str, Any]) -> str | list[Any] | dict[str, Any]:
    """Convert an existing engine/manifest scope to public schema-4 names."""
    result = _map_scope_names(copy.deepcopy(raw), to_internal=False)
    if result.get("type") == "full_sequence":
        return "full_sequence"
    if result.get("type") == "parts":
        terms = []
        for span in result.get("spans", []):
            if set(span) != {"part"}:
                return result
            terms.append(span["part"])
        compact = []
        index = 0
        while index < len(terms):
            if terms[index:index + 3] == ["start_codon", "cds_body", "stop_codon"]:
                compact.append("cds")
                index += 3
            else:
                compact.append(terms[index])
                index += 1
        terms = compact
        return terms[0] if len(terms) == 1 else terms
    return result


legacy_scope_to_public_v4 = scope_from_internal


def scope_expression(raw: Any) -> str:
    """Editable public expression for strings, lists and detailed mappings."""
    if isinstance(raw, str):
        scope_to_internal(raw)
        return raw
    internal = scope_to_internal(raw)
    text = scorer_scope_expression(internal)
    text = re.sub(r"\bcds\b", "cds_body", text)
    return text.replace("full_cds", "cds")


def resolve_paths(raw: dict[str, Any], base_dir: Path, *, output_base_dir: Path | None = None) -> None:
    construct = _mapping(raw.get("construct"), "construct")
    construct["fasta"] = _path(construct.get("fasta"), base_dir, "construct.fasta")
    changes = _mapping(construct.get("changes", {}), "construct.changes")
    for part, choice in changes.items():
        if isinstance(choice, dict) and "library" in choice:
            choice["library"] = _path(choice["library"], base_dir, f"construct.changes.{part}.library")
    if raw.get("output_dir") is not None:
        raw["output_dir"] = _path(raw["output_dir"], output_base_dir or base_dir, "output_dir")


def build_design(
    name: str, raw: dict[str, Any], base_dir: Path, *, source_sequence: str | None = None
) -> ConstructDesignSpace:
    _reject(raw, {"fasta", "coordinates", "changes", "allow_upstream_start_codons"}, "construct")
    changes = _mapping(raw.get("changes", {}), "construct.changes")
    _reject(changes, set(PART_ORDER), "construct.changes")
    source_raw = {key: value for key, value in raw.items() if key != "changes"}
    design = construct_design_from_config(name, source_raw, base_dir, set(changes),
                                         source_sequence=source_sequence)
    return _apply_changes(design, changes, base_dir)


def _apply_changes(
    design: ConstructDesignSpace, changes: dict[str, Any], base_dir: Path
) -> ConstructDesignSpace:
    """Apply the shared library/synonymous rules after parsing reference parts."""
    updates: dict[str, Any] = {}
    source = design.source_construct
    for part, choice in changes.items():
        if part == "cds":
            if choice != "synonymous":
                raise ValueError("construct.changes.cds must be synonymous")
        elif part in {"five_prime_utr", "three_prime_utr"}:
            rule = _mapping(choice, f"construct.changes.{part}")
            _reject(rule, {"library", "set", "max_candidates"}, f"construct.changes.{part}")
            library = Path(_path(rule.get("library"), base_dir, f"construct.changes.{part}.library"))
            reference = getattr(source, part)
            variants = _library(library, part, rule)
            matching = [item for item in variants if item.sequence == reference.sequence]
            if matching:
                selected = matching[0]
                variants = [selected] + [item for item in variants if item is not selected]
            else:
                if any(item.name == reference.name for item in variants):
                    reference = replace(reference, name=f"{part}_input")
                variants = [reference] + variants
            updates["five_prime_utrs" if part == "five_prime_utr" else "three_prime_utrs"] = tuple(variants)
        elif part == "poly_a":
            rule = _mapping(choice, "construct.changes.poly_a")
            _reject(rule, {"lengths"}, "construct.changes.poly_a")
            lengths = rule.get("lengths")
            if not isinstance(lengths, list) or not lengths or any(type(value) is not int or value < 0 for value in lengths):
                raise ValueError("construct.changes.poly_a.lengths must contain non-negative integers")
            if len(lengths) != len(set(lengths)):
                raise ValueError("construct.changes.poly_a.lengths contains duplicates")
            updates["poly_a_lengths"] = tuple([source.poly_a_length] + [value for value in lengths if value != source.poly_a_length])
        elif part in {"start_codon", "stop_codon"}:
            rule = _mapping(choice, f"construct.changes.{part}")
            _reject(rule, {"sequences"}, f"construct.changes.{part}")
            values = rule.get("sequences")
            if not isinstance(values, list) or not values:
                raise ValueError(f"construct.changes.{part}.sequences must be a non-empty list")
            variants = [normalize_dna(str(item)) for item in values]
            reference = getattr(source, part)
            updates["start_codons" if part == "start_codon" else "stop_codons"] = tuple([reference] + [item for item in variants if item != reference])
    design = replace(design, **updates)
    # Reuse physical-part content checks for choices that do not pass the old parser.
    if any(len(value) != 3 or value in {"TAA", "TAG", "TGA"} for value in design.start_codons):
        raise ValueError("Start choices must be three-nucleotide non-stop codons")
    if any(value not in {"", "TAA", "TAG", "TGA"} for value in design.stop_codons):
        raise ValueError("Stop choices must be terminal stop codons")
    if any(len(values) != len(set(values)) for values in (design.start_codons, design.stop_codons)):
        raise ValueError("Codon choices contain duplicates")
    for variant in design.five_prime_utrs:
        design.validate(replace(design.source_construct, five_prime_utr=variant))
    for variant in design.three_prime_utrs:
        design.validate(replace(design.source_construct, three_prime_utr=variant))
    design.validate(design.source_construct)
    return design


def _library(path: Path, part: str, rule: dict[str, Any]) -> list[SequenceVariant]:
    if path.suffix.lower() in {".tsv", ".tab"}:
        records = load_part_library(path, part, rule.get("set"), rule.get("max_candidates"))
        return [SequenceVariant(item.id, item.sequence, {"library": str(path), "part_type": item.type, "part_set": item.set_name, **item.metadata}) for item in records]
    if rule.get("set") is not None:
        raise ValueError("A library set is supported only for TSV libraries")
    records = []
    header = None
    chunks = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append(SequenceVariant(header, normalize_dna("".join(chunks))))
            header = line[1:].split()[0] if line[1:].strip() else ""
            chunks = []
        else:
            if header is None:
                raise ValueError(f"FASTA library {path} has sequence before a header")
            chunks.append(line)
    if header is not None:
        records.append(SequenceVariant(header, normalize_dna("".join(chunks))))
    if not records or any(not item.name or not item.sequence for item in records):
        raise ValueError(f"FASTA library {path} must contain named, non-empty records")
    if len({item.name for item in records}) != len(records) or len({item.sequence for item in records}) != len(records):
        raise ValueError(f"FASTA library {path} contains duplicate names or sequences")
    limit = rule.get("max_candidates")
    if limit is not None:
        if type(limit) is not int or limit < 1:
            raise ValueError("library.max_candidates must be a positive integer")
        records = records[:limit]
    return records


def engine_config(raw: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(raw)
    schema = result.get("schema_version", 4)
    optimizer = _mapping(result.get("optimizer"), "optimizer")
    if "type" in optimizer or "variables" in optimizer:
        raise ValueError(f"Schema {schema} uses optimizer.algorithm and construct.changes")
    optimizer["type"] = optimizer.pop("algorithm", "nsga2")
    names = []
    for scorer in result.get("scorers", []):
        if not isinstance(scorer, dict):
            raise ValueError("Each scorer must be a mapping")
        if not isinstance(scorer.get("command"), str) or not scorer["command"].strip():
            raise ValueError(f"Schema-{schema} scorers require a prepared Bash command string")
        if any(key in scorer for key in ("outputs", "mode", "image", "docker_args")):
            raise ValueError(f"Schema {schema} uses prepared command strings and metrics, not legacy scorer runtime fields")
        name = scorer.get("name")
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]{0,63}", name):
            raise ValueError("Scorer names must start with a letter and use letters, numbers, underscore or hyphen")
        if name in names:
            raise ValueError(f"Duplicate scorer name: {name}")
        names.append(name)
        if "output_prefix" in scorer:
            raise ValueError(f"Schema {schema} namespaces metrics automatically as scorer.metric")
        if "namespace" in scorer:
            raise ValueError(f"Schema {schema} derives metric namespaces from scorer names")
        scorer["namespace"] = name
        scorer["scope"] = scope_to_internal(scorer.get("scope"))
    return result


def _mapping(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be a mapping")
    return value


def _reject(raw: dict[str, Any], allowed: set[str], context: str) -> None:
    unknown = set(raw) - allowed
    if unknown:
        raise ValueError(f"{context} contains unknown fields: " + ", ".join(sorted(map(str, unknown))))


def _path(value: Any, base_dir: Path, context: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{context} must be a non-empty path")
    path = Path(value).expanduser()
    return str(path.resolve() if path.is_absolute() else (base_dir / path).resolve())
