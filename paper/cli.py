"""Public command-line interface for paper evidence and experiments."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .io import HERE, registry


def main(argv=None):
    parser = argparse.ArgumentParser(prog="morse paper", description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    p = subs.add_parser("verify", help="Verify frozen bundle and scientific invariants")
    p.add_argument("--raw", action="store_true", help="Also verify all 155 archived run-source hashes")
    p.add_argument("--raw-root", type=Path, default=HERE / "artifacts/raw")
    p = subs.add_parser("figures", help="Regenerate all eight figures from frozen numerical observations")
    p.add_argument("--output", type=Path, default=Path("results/paper"))
    p.add_argument("--format", choices=["pdf", "png", "both"], default="both")
    p = subs.add_parser("analyze", help="Recompute scientific tables and optionally audit all raw evaluations")
    p.add_argument("--output", type=Path, default=Path("results/paper/analysis"))
    p.add_argument("--raw", action="store_true")
    p.add_argument("--raw-root", type=Path, default=HERE / "artifacts/raw")
    p.add_argument("--bootstrap", action="store_true", help="Recompute the full 200,000 paired bootstrap")
    p.add_argument("--runs-root",type=Path,help="Also compare freshly executed run artifacts with frozen selections/scores/budgets")
    p = subs.add_parser("baseline", help="Run the separate native VaxPress benchmark")
    p.add_argument("tool", choices=["vaxpress"])
    p.add_argument("--seed", type=int, choices=[7,43,101,137,211])
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--execute", action="store_true")
    p.add_argument("--runtime", choices=["container", "process"], default="container")
    p.add_argument("--image-map", type=Path)
    p.add_argument("--scoring-events", type=int, help="Diagnostic sweep total divisible by five; default12800")
    p = subs.add_parser("upstream", help="Regenerate an upstream design with the separately installed authors' tool")
    p.add_argument("tool", choices=["lineardesign", "ribodecode"])
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--source", type=Path, help="External LinearDesign checkout at the frozen revision")
    p.add_argument("--seed", type=int, choices=[7, 43, 101, 137, 211], default=7)
    p.add_argument("--image", help="Installed native RiboDecode generation image")
    p.add_argument("--device", default="cuda:0")
    p.add_argument("--execute", action="store_true")
    p.add_argument("--pilot", action="store_true", help="RiboDecode diagnostic: one epoch and 400 final samples")
    p = subs.add_parser("archive", help="Pack or verify/extract the optional raw-source release artifact")
    p.add_argument("action",choices=["pack","extract"])
    p.add_argument("--file",type=Path,required=True)
    p.add_argument("--raw-root",type=Path,default=HERE/"artifacts/raw")
    p = subs.add_parser("rescore",help="Replay auxiliary measurements using the real ready scorers")
    p.add_argument("assay",choices=["mauger","fig3-logged"])
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--execute",action="store_true")
    p.add_argument("--limit",type=int,help="Diagnostic target prefix; original model-batch context is retained")
    p.add_argument("--batch-size",type=int,help="Diagnostic override of frozen per-scorer request grouping")
    p.add_argument("--scorer",action="append",help="Diagnostic selection of scorer IDs (repeatable)")
    p.add_argument("--device",default="cuda:0")
    p.add_argument("--runtime",choices=["container","process"],default="container")
    p.add_argument("--asset-root",type=Path,help="Explicit host asset override; self-contained published containers otherwise use bundled assets")
    p.add_argument("--image-map",type=Path,help="Image overrides JSON, or assets/runtime-lock.json for published immutable images")
    p.add_argument("--source-overlay",action="store_true",help="Overlay current adapter source; self-contained images keep their bundled models and packages")
    subs.add_parser("list", help="List the 40 independent runs/sweeps and frozen costs")
    args = parser.parse_args(argv)
    try:
        if args.command == "list":
            for r in registry()["runs"]:
                print(f'{r["study"]:6s} {r["run_id"]:25s} seed={r["seed"]:3d} events={r["recorded_evaluations"]} generations={r["generations"]}')
        elif args.command == "verify":
            from .verification import verify
            print(json.dumps(verify(raw=args.raw, raw_root=args.raw_root), indent=2))
        elif args.command == "figures":
            from .figures import build
            print(json.dumps(build(args.output, args.format), indent=2))
        elif args.command == "analyze":
            from .analysis import analyze
            print(json.dumps(analyze(args.output, raw=args.raw, raw_root=args.raw_root, bootstrap=args.bootstrap,runs_root=args.runs_root), indent=2))
        elif args.command == "baseline":
            from .runtime import run_baseline
            print(json.dumps(run_baseline(args), indent=2))
        elif args.command == "upstream":
            from .upstream import run_upstream
            print(json.dumps(run_upstream(args), indent=2))
        elif args.command == "archive":
            from .archive import manage
            print(json.dumps(manage(args.action,args.file,args.raw_root),indent=2))
        elif args.command == "rescore":
            from .rescore import run
            result=run(args)
            print(json.dumps(result,indent=2))
            if result.get("status")=="completed_with_differences":return 1
        return 0
    except (ValueError, RuntimeError, FileNotFoundError) as exc:
        print(f"paper: {exc}", file=sys.stderr)
        return 2
