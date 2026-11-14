from __future__ import annotations

import argparse
import json
from pathlib import Path

from morse.run_spec import load_run_spec

from morse_warm_starter.generation import (
    AnchorRequest,
    GenerationOptions,
    generate_population,
)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except Exception as exc:
        parser.exit(1, f"error: {exc}\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="morse-warm-start",
        description=(
            "Generate a scorer-guided initial-population manifest for MORSE Core."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser(
        "generate",
        help="Generate anchors and optionally fill the remaining population randomly.",
    )
    generate.add_argument("--config", type=Path, required=True)
    generate.add_argument("--output", type=Path, required=True)
    generate.add_argument("--population-size", type=int, default=None)
    generate.add_argument("--seed", type=int, default=None)
    generate.add_argument(
        "--exhaustive",
        action="append",
        default=[],
        metavar="SCORER:OBJECTIVE[:COUNT]",
        help=(
            "Exhaustively optimize mutable variables in a scorer's input scope; "
            "repeat for multiple anchor groups."
        ),
    )
    generate.add_argument(
        "--greedy",
        action="append",
        default=[],
        metavar="SCORER:OBJECTIVE[:COUNT]",
        help=(
            "Coordinate-descent optimization in a scorer's input scope; "
            "repeat for multiple anchor groups."
        ),
    )
    generate.add_argument(
        "--no-original",
        action="store_true",
        help="Do not include the reference construct as an anchor.",
    )
    generate.add_argument(
        "--anchors-only",
        action="store_true",
        help=(
            "Write only generated anchors; Core will randomly fill missing slots. "
            "By default this CLI writes a full population."
        ),
    )
    generate.add_argument("--max-combinations", type=int, default=100_000)
    generate.add_argument("--batch-size", type=int, default=512)
    generate.add_argument("--passes", type=int, default=1)
    generate.add_argument("--candidate-pool", type=int, default=100)
    generate.add_argument("--min-hamming-distance", type=int, default=1)
    generate.add_argument(
        "--lineardesign-profile",
        type=Path,
        default=None,
        help="Human, mouse, or custom CAI profile in JSON/LinearDesign CSV format.",
    )
    generate.add_argument(
        "--lineardesign-lambda",
        type=float,
        action="append",
        default=[],
        help="LinearDesign lambda; repeat the flag to create a lambda grid.",
    )
    generate.add_argument(
        "--lineardesign-dir",
        type=Path,
        default=Path("~/.cache/morse/LinearDesign").expanduser(),
    )
    generate.add_argument(
        "--accept-lineardesign-license",
        action="store_true",
        help="Confirm that the upstream LinearDesign license was reviewed and accepted.",
    )
    generate.set_defaults(func=cmd_generate)
    return parser


def cmd_generate(args: argparse.Namespace) -> int:
    spec = load_run_spec(args.config)
    population_size = (
        spec.run_config.optimizer.population_size
        if args.population_size is None
        else args.population_size
    )
    seed = (
        spec.run_config.optimizer.seed
        if args.seed is None
        else args.seed
    )
    requests = [
        *(_parse_request("exhaustive", value) for value in args.exhaustive),
        *(_parse_request("greedy", value) for value in args.greedy),
    ]
    lineardesign_candidates = []
    if args.lineardesign_lambda:
        if args.lineardesign_profile is None:
            raise ValueError(
                "--lineardesign-profile is required with --lineardesign-lambda"
            )
        from morse_warm_starter.lineardesign import (
            generate_lineardesign_candidates,
        )

        lineardesign_candidates = generate_lineardesign_candidates(
            spec.design_space,
            profile_path=args.lineardesign_profile,
            lambdas=args.lineardesign_lambda,
            lineardesign_dir=args.lineardesign_dir,
            work_dir=args.output.parent / "lineardesign-profile",
            accept_license=args.accept_lineardesign_license,
            preserve_start_codon=spec.resolved["schema_version"] == 5,
        )
    manifest = generate_population(
        design_space=spec.design_space,
        scorer_configs=spec.run_config.scorers,
        requests=requests,
        options=GenerationOptions(
            population_size=population_size,
            seed=seed,
            include_original=not args.no_original,
            anchors_only=args.anchors_only,
            max_combinations=args.max_combinations,
            batch_size=args.batch_size,
            passes=args.passes,
            candidate_pool=args.candidate_pool,
            min_hamming_distance=args.min_hamming_distance,
        ),
        prebuilt_candidates=lineardesign_candidates,
    )
    manifest["metadata"]["config"] = str(args.config.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "candidates": len(manifest["candidates"]),
                "population_size": population_size,
                "design_space_fingerprint": manifest[
                    "design_space_fingerprint"
                ],
            },
            sort_keys=True,
        )
    )
    return 0


def _parse_request(mode: str, value: str) -> AnchorRequest:
    parts = value.split(":")
    if len(parts) not in {2, 3} or not parts[0] or not parts[1]:
        raise ValueError(
            f"{mode} request must be SCORER:OBJECTIVE[:COUNT], got {value!r}"
        )
    try:
        count = 1 if len(parts) == 2 else int(parts[2])
    except ValueError as exc:
        raise ValueError(
            f"{mode} request count must be an integer, got {value!r}"
        ) from exc
    return AnchorRequest(
        mode=mode,
        scorer=parts[0],
        objective=parts[1],
        count=count,
    )


if __name__ == "__main__":
    raise SystemExit(main())
