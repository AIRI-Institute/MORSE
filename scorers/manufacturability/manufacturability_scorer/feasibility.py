"""Versioned, deterministic transcript-body feasibility checks (no dependencies)."""

from __future__ import annotations

import copy
import hashlib
import json
import math
from bisect import bisect_left, bisect_right
from collections import Counter
from pathlib import Path
from typing import Any

from .sequence import homopolymer_runs, normalize_rna, validate_forbidden_motifs

CHECKS = (
    "global_gc",
    "local_gc",
    "homopolymers",
    "direct_repeats",
    "tandem_repeats",
    "forbidden_motifs",
    "restriction_sites",
    "capping_start",
    "inverted_repeats",
)
DEFAULT_PRESET = "sequence_feasibility_v1"
PRESETS = {
    name: Path(__file__).parent / "profiles" / f"{name}.json"
    for name in (DEFAULT_PRESET, "manufacturing_strategy_v1")
}
JSON_PARAMETERS = frozenset(
    {"enabled_checks", "forbidden_motifs", "exempt_regions", "restriction_enzymes"}
)
STRING_PARAMETERS = frozenset({"capping_method", "custom_5prime_prefix"})
RESTRICTION_SITES = {"BsaI": "GGUCUC", "BsmBI": "CGUCUC"}
STRATEGY_DEFAULTS = {
    "restriction_enzymes": [],
    "capping_method": "none",
    "custom_5prime_prefix": "",
    "inverted_repeat_stem": 20,
    "inverted_repeat_window": 200,
    "inverted_repeat_min_gap": 0,
}
FLOAT_PARAMETERS = frozenset(
    {"global_gc_min", "global_gc_max", "local_gc_min", "local_gc_max"}
)
RNA_COMPLEMENT = str.maketrans("ACGU", "UGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RNA_COMPLEMENT)[::-1]


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


class SequenceFeasibilityScorer:
    def __init__(self, preset: str = DEFAULT_PRESET, **overrides: Any) -> None:
        if preset not in PRESETS:
            raise ValueError(f"Unknown feasibility preset: {preset!r}")
        payload = PRESETS[preset].read_bytes()
        profile = json.loads(payload)
        if profile["profile_id"] != preset or profile["schema_version"] != 1:
            raise ValueError("Invalid feasibility preset identity or schema")
        parameters = profile["parameters"]
        requested_checks = overrides.get("enabled_checks", [])
        if not isinstance(requested_checks, list):
            raise ValueError("enabled_checks must be a JSON list")
        # Preserve the released v1 configuration/hash unless strategy options
        # are explicitly requested. New profiles include these defaults in full.
        if set(overrides) & set(STRATEGY_DEFAULTS) or any(
            rule in requested_checks for rule in CHECKS[6:]
        ):
            parameters = {**copy.deepcopy(STRATEGY_DEFAULTS), **parameters}
        unknown = set(overrides) - set(parameters)
        if unknown:
            raise ValueError(
                "Unknown feasibility parameters: " + ", ".join(sorted(unknown))
            )
        parameters.update(copy.deepcopy(overrides))
        self._validate(parameters)
        self._parameters = parameters
        self.profile_id = preset
        self.profile_sha256 = hashlib.sha256(payload).hexdigest()
        self.configuration_json = _json(
            {
                "preset": preset,
                "preset_sha256": self.profile_sha256,
                "parameters": parameters,
            }
        )
        self.configuration_sha256 = hashlib.sha256(
            self.configuration_json.encode()
        ).hexdigest()

    @property
    def parameters(self) -> dict[str, Any]:
        return copy.deepcopy(self._parameters)

    @staticmethod
    def _validate(p: dict[str, Any]) -> None:
        for key, value in p.items():
            if key in JSON_PARAMETERS:
                if not isinstance(value, list):
                    raise ValueError(f"{key} must be a JSON list")
            elif key in STRING_PARAMETERS:
                if not isinstance(value, str):
                    raise ValueError(f"{key} must be a string")
            elif key == "inverted_repeat_min_gap":
                if type(value) is not int or value < 0:
                    raise ValueError(f"{key} must be a nonnegative integer")
            elif key in FLOAT_PARAMETERS:
                if (
                    isinstance(value, bool)
                    or not isinstance(value, (int, float))
                    or not math.isfinite(value)
                    or not 0 <= value <= 1
                ):
                    raise ValueError(f"{key} must be a finite fraction in [0, 1]")
                p[key] = float(value)
            elif type(value) is not int or value < 1:
                raise ValueError(f"{key} must be a positive integer")
        for prefix in ("global_gc", "local_gc"):
            if p[f"{prefix}_min"] > p[f"{prefix}_max"]:
                raise ValueError(f"{prefix}_min must not exceed {prefix}_max")
        enabled = p["enabled_checks"]
        if not all(isinstance(v, str) and v in CHECKS for v in enabled) or len(
            set(enabled)
        ) != len(enabled):
            raise ValueError("enabled_checks must contain unique known check names")
        p["enabled_checks"] = [rule for rule in CHECKS if rule in enabled]
        p["forbidden_motifs"] = sorted(validate_forbidden_motifs(p["forbidden_motifs"]))
        if "restriction_enzymes" in p:
            enzymes = p["restriction_enzymes"]
            if not all(
                isinstance(e, str) and e in RESTRICTION_SITES for e in enzymes
            ) or len(set(enzymes)) != len(enzymes):
                raise ValueError(
                    "restriction_enzymes must contain unique supported names: BsaI, BsmBI"
                )
            p["restriction_enzymes"] = sorted(enzymes)
            if p["capping_method"] not in ("none", "cleancap_ag", "custom"):
                raise ValueError("capping_method must be none, cleancap_ag, or custom")
            prefix = p["custom_5prime_prefix"]
            if prefix and any(base not in "ACGU" for base in prefix):
                raise ValueError("custom_5prime_prefix must be exact uppercase RNA")
            if p["capping_method"] == "custom" and not prefix:
                raise ValueError("Custom capping requires a nonempty custom_5prime_prefix")
            if p["capping_method"] != "custom" and prefix:
                raise ValueError("custom_5prime_prefix is only used with custom capping")
            if p["inverted_repeat_window"] < (
                2 * p["inverted_repeat_stem"] + p["inverted_repeat_min_gap"]
            ):
                raise ValueError(
                    "inverted_repeat_window must fit both arms and minimum gap"
                )
        if p["direct_repeat_min_copies"] < 2:
            raise ValueError("direct_repeat_min_copies must be at least 2")
        if not 2 <= p["tandem_min_period"] <= p["tandem_max_period"]:
            raise ValueError("Tandem periods must satisfy 2 <= min <= max")
        if p["tandem_min_copies"] < 2:
            raise ValueError("tandem_min_copies must be at least 2")
        for region in p["exempt_regions"]:
            if not isinstance(region, dict) or set(region) != {
                "start",
                "end",
                "checks",
                "reason",
            }:
                raise ValueError(
                    "Each exempt region requires start, end, checks, and reason"
                )
            start, end = region["start"], region["end"]
            if type(start) is not int or type(end) is not int or not 0 <= start < end:
                raise ValueError("Exempt regions require integer 0 <= start < end")
            if not isinstance(region["reason"], str) or not region["reason"].strip():
                raise ValueError("Exempt regions require a nonempty reason")
            rules = region["checks"]
            if (
                not isinstance(rules, list)
                or not rules
                or not all(isinstance(v, str) and v in CHECKS[1:] for v in rules)
                or len(set(rules)) != len(rules)
            ):
                raise ValueError(
                    "Exempt checks must be unique known local checks; global_gc cannot be exempted"
                )
            region["checks"] = [rule for rule in CHECKS if rule in rules]
        p["exempt_regions"].sort(
            key=lambda r: (r["start"], r["end"], r["checks"], r["reason"])
        )

    def score_sequence(self, raw_sequence: str) -> dict[str, Any]:
        if not isinstance(raw_sequence, str):
            raise ValueError("Sequence must be a string")
        sequence = normalize_rna(raw_sequence)
        n = len(sequence)
        if not n:
            raise ValueError("Routed sequence must not be empty")
        p = self._parameters
        if any(region["end"] > n for region in p["exempt_regions"]):
            raise ValueError("Exempt region extends beyond the routed sequence")
        events: list[dict[str, Any]] = []
        prefix_gc = [0]
        for base in sequence:
            prefix_gc.append(prefix_gc[-1] + (base in "GC"))
        global_gc = prefix_gc[-1] / n
        if not p["global_gc_min"] <= global_gc <= p["global_gc_max"]:
            events.append(
                {
                    "check": "global_gc",
                    "spans": [[0, n]],
                    "value": global_gc,
                    "bounds": [p["global_gc_min"], p["global_gc_max"]],
                }
            )
        width = min(n, p["local_gc_window"])
        starts = list(range(0, n - width + 1, p["local_gc_stride"]))
        if starts[-1] != n - width:
            starts.append(n - width)
        gc_values = []
        for start in starts:
            value = (prefix_gc[start + width] - prefix_gc[start]) / width
            gc_values.append(value)
            if not p["local_gc_min"] <= value <= p["local_gc_max"]:
                events.append(
                    {
                        "check": "local_gc",
                        "spans": [[start, start + width]],
                        "value": value,
                        "bounds": [p["local_gc_min"], p["local_gc_max"]],
                    }
                )
        max_runs = dict.fromkeys("ACGU", 0)
        for run in homopolymer_runs(sequence):
            max_runs[run.base] = max(max_runs[run.base], run.length)
            limit = p[f"homopolymer_max_{run.base.lower()}"]
            if run.length > limit:
                events.append(
                    {
                        "check": "homopolymers",
                        "spans": [[run.start, run.end]],
                        "base": run.base,
                        "length": run.length,
                        "max_length": limit,
                    }
                )
        k = p["direct_repeat_length"]
        occurrences: dict[str, list[int]] = {}
        for start in range(n - k + 1):
            motif = sequence[start : start + k]
            occurrences.setdefault(motif, []).append(start)
        for motif, positions in occurrences.items():
            # Equal-length greedy interval selection gives the maximum number
            # of pairwise nonoverlapping copies. Retain ALL occurrences for
            # exemption checking so an unselected occurrence cannot be hidden.
            selected = []
            for start in positions:
                if not selected or start >= selected[-1] + k:
                    selected.append(start)
            if len(selected) >= p["direct_repeat_min_copies"]:
                events.append(
                    {
                        "check": "direct_repeats",
                        "spans": [[s, s + k] for s in positions],
                        "motif": motif,
                        "copies": len(selected),
                    }
                )
        for start, end, period in tandem_runs(
            sequence,
            p["tandem_min_period"],
            p["tandem_max_period"],
            p["tandem_min_length"],
            p["tandem_min_copies"],
        ):
            events.append(
                {
                    "check": "tandem_repeats",
                    "spans": [[start, end]],
                    "period": period,
                    "motif": sequence[start : start + period],
                    "length": end - start,
                    "complete_copies": (end - start) // period,
                }
            )
        for motif in p["forbidden_motifs"]:
            start = sequence.find(motif)
            while start >= 0:
                events.append(
                    {
                        "check": "forbidden_motifs",
                        "spans": [[start, start + len(motif)]],
                        "motif": motif,
                    }
                )
                start = sequence.find(motif, start + 1)
        for enzyme in p.get("restriction_enzymes", []):
            motif = RESTRICTION_SITES[enzyme]
            for strand, word in (("+", motif), ("-", reverse_complement(motif))):
                start = sequence.find(word)
                while start >= 0:
                    events.append({
                        "check": "restriction_sites",
                        "spans": [[start, start + len(word)]],
                        "enzyme": enzyme, "strand": strand, "motif": word,
                    })
                    start = sequence.find(word, start + 1)
        method = p.get("capping_method", "none")
        required = ""
        if method == "cleancap_ag":
            required = "AG"
        elif method == "custom":
            required = p["custom_5prime_prefix"]
        if required and not sequence.startswith(required):
            events.append({
                "check": "capping_start", "spans": [[0, min(n, len(required))]],
                "method": method, "required_prefix": required,
                "observed_prefix": sequence[:len(required)],
            })
        if "inverted_repeat_stem" in p:
            stem = p["inverted_repeat_stem"]
            for left, right in inverted_repeat_pairs(
                sequence, stem, p["inverted_repeat_window"], p["inverted_repeat_min_gap"]
            ):
                events.append({
                    "check": "inverted_repeats",
                    "spans": [[left, left + stem], [right, right + stem]],
                    "stem_length": stem, "gap": right - left - stem,
                    "window": p["inverted_repeat_window"],
                })
        raw_counts = Counter(event["check"] for event in events)
        active, exempt = [], []
        for event in events:
            if event["check"] not in p["enabled_checks"]:
                continue
            reasons = []
            for start, end in event["spans"]:
                matches = [
                    region["reason"]
                    for region in p["exempt_regions"]
                    if event["check"] in region["checks"]
                    and region["start"] <= start
                    and end <= region["end"]
                ]
                if not matches:
                    break
                reasons.extend(matches)
            else:
                exempt.append({**event, "reasons": sorted(set(reasons))})
                continue
            active.append(event)
        active_counts = Counter(event["check"] for event in active)
        return {
            "feasibility_pass": not active,
            "feasibility_violation_count": len(active),
            "feasibility_failed_check_count": len(active_counts),
            "feasibility_failed_checks_json": _json(
                [rule for rule in CHECKS if active_counts[rule]]
            ),
            "feasibility_violations_json": _json(active),
            "feasibility_exemptions_json": _json(exempt),
            "feasibility_raw_counts_json": _json(
                {rule: raw_counts[rule] for rule in CHECKS}
            ),
            "analyzed_length": n,
            "global_gc_fraction": global_gc,
            "local_gc_min": min(gc_values),
            "local_gc_max": max(gc_values),
            "local_gc_window_count": len(starts),
            "max_homopolymer_length": max(max_runs.values()),
            **{
                f"max_homopolymer_{base.lower()}": value
                for base, value in max_runs.items()
            },
            **{f"{rule}_violation_count": active_counts[rule] for rule in CHECKS},
            "feasibility_profile_id": self.profile_id,
            "feasibility_profile_sha256": self.profile_sha256,
            "feasibility_config_sha256": self.configuration_sha256,
            "feasibility_config_json": self.configuration_json,
        }


def inverted_repeat_pairs(sequence: str, stem: int, window: int, min_gap: int = 0):
    """Exact nonoverlapping reverse-complement arms within an inclusive span.

    Zero gap includes reverse-complement palindromes of length 2*stem.
    Enumerates qualifying arm pairs, not predicted thermodynamic structures.
    """
    positions: dict[str, list[int]] = {}
    for i in range(len(sequence) - stem + 1):
        positions.setdefault(sequence[i:i + stem], []).append(i)
    for left in range(len(sequence) - 2 * stem - min_gap + 1):
        matches = positions.get(reverse_complement(sequence[left:left + stem]), [])
        lo = bisect_left(matches, left + stem + min_gap)
        hi = bisect_right(matches, left + window - stem)
        for right in matches[lo:hi]:
            yield left, right


def tandem_runs(
    sequence: str, min_period: int, max_period: int, min_length: int, min_copies: int
) -> list[tuple[int, int, int]]:
    """Maximal exact periodic intervals with a primitive (nonreducible) motif.

    A partial trailing copy contributes to length, but not complete-copy count.
    Comparisons are linear in sequence length for each tested period.
    """
    result = []
    for period in range(min_period, min(max_period, len(sequence) // min_copies) + 1):
        match_start = None
        for index in range(period, len(sequence) + 1):
            if index < len(sequence) and sequence[index] == sequence[index - period]:
                if match_start is None:
                    match_start = index
            elif match_start is not None:
                start, end = match_start - period, index
                match_start = None
                if end - start < max(min_length, min_copies * period):
                    continue
                motif = sequence[start : start + period]
                if any(
                    period % d == 0 and motif == motif[:d] * (period // d)
                    for d in range(1, period)
                ):
                    continue
                result.append((start, end, period))
    return sorted(result)
