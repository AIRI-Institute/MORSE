"""FASTA interface; preset overrides use the same names as the manifest."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .feasibility import DEFAULT_PRESET, PRESETS, SequenceFeasibilityScorer
from .sequence import normalize_rna, read_fasta


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("--preset", choices=tuple(PRESETS), default=DEFAULT_PRESET)
    parser.add_argument(
        "--parameters", type=Path, help="JSON object with named preset overrides"
    )
    parser.add_argument(
        "--poly-a-length",
        type=int,
        default=0,
        help="Explicit annotated terminal poly(A) length to exclude",
    )
    args = parser.parse_args(argv)
    if args.poly_a_length < 0:
        parser.error("--poly-a-length must be nonnegative")
    parameters = (
        {} if args.parameters is None else json.loads(args.parameters.read_text())
    )
    if not isinstance(parameters, dict):
        parser.error("--parameters must contain a JSON object")
    scorer = SequenceFeasibilityScorer(args.preset, **parameters)
    for name, raw in read_fasta(args.fasta):
        sequence = normalize_rna(raw)
        if args.poly_a_length:
            if (
                args.poly_a_length >= len(sequence)
                or sequence[-args.poly_a_length :] != "A" * args.poly_a_length
            ):
                parser.error(
                    "Annotated poly(A) must be a terminal A-only segment shorter than the sequence"
                )
            sequence = sequence[: -args.poly_a_length]
        print(
            json.dumps(
                {"id": name, "outputs": scorer.score_sequence(sequence)},
                sort_keys=True,
                allow_nan=False,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
