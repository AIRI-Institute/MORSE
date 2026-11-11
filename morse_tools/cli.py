"""Thin entry point; scientific execution remains in the existing Core."""
from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def repository_root() -> Path:
    return Path(os.environ.get("MORSE_ROOT", Path(__file__).resolve().parents[1])).expanduser().resolve()


def configure_environment() -> Path:
    root = repository_root()
    default_assets = root / ".cache/morse/assets" if (root / "pyproject.toml").is_file() else Path.home() / ".cache/morse/assets"
    os.environ.setdefault("MORSE_ASSETS_DIR", str(default_assets))
    return root


def _paper_bundle_available(root: Path) -> bool:
    return (root / "paper/registry.json").is_file() and all(
        (root / "paper" / name).is_dir() for name in ("data", "configs", "inputs", "reference", "adapters")
    )


def _require_checkout_resources(root: Path, command: str) -> None:
    available = _paper_bundle_available(root) if command == "paper" else any((root / "scorers").glob("*/morse-scorer.yaml"))
    if not available:
        raise SystemExit(
            f"morse {command} requires the complete monorepo checkout and its data files. "
            f"Resources were not found under {root}. Install with uv from that checkout, "
            "or set MORSE_ROOT to its absolute path. A standalone wheel supports "
            "independent experiments and native CPU scorers; the full paper bundle requires a checkout."
        )


def doctor(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="morse doctor", description="Check the active research environment.")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--gpu", action="store_true", help="also inspect NVIDIA GPU availability")
    args = parser.parse_args(argv)
    root = configure_environment()
    modules = {name: importlib.util.find_spec(name) is not None for name in ["morse", "morse_warm_starter", "manufacturability_scorer", "numpy", "pymoo", "yaml", "matplotlib"]}
    versions = {}
    for name in ["morse-research", "morse-core", "morse-warm-starter", "numpy", "pymoo", "matplotlib", "scipy"]:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    result = {"python": sys.version.split()[0], "repository": str(root), "assets_directory": os.environ["MORSE_ASSETS_DIR"], "modules": modules, "versions": versions, "docker_available": shutil.which("docker") is not None, "paper_manifest_present": (root / "paper" / "registry.json").is_file(), "paper_bundle_available": _paper_bundle_available(root)}
    if args.gpu:
        if shutil.which("nvidia-smi"):
            try:
                completed = subprocess.run(["nvidia-smi", "--query-gpu=index,name,memory.total,memory.free", "--format=csv,noheader"], capture_output=True, text=True, timeout=15)
                result["gpus"] = completed.stdout.strip().splitlines() if completed.returncode == 0 else []
            except (OSError, subprocess.TimeoutExpired) as exc:
                result["gpus"] = []
                result["gpu_error"] = str(exc)
        else:
            result["gpus"] = []
    result["ready"] = all(modules.values()) and (bool(result.get("gpus")) if args.gpu else True)
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"MORSE {versions['morse-research'] or 'source'} | Python {result['python']}")
        print(f"Repository: {root}")
        print(f"Ready scorers cache: {result['assets_directory']}")
        for name, available in modules.items():
            print(f"  {'OK' if available else 'MISSING'} {name}")
        if args.gpu:
            for gpu in result["gpus"]:
                print(f"  GPU {gpu}")
        print("Ready" if result["ready"] else "Install missing dependencies before running experiments.")
    return 0 if result["ready"] else 1


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    root = configure_environment()
    if not argv or argv[0] in {"-h", "--help"}:
        print("MORSE research tools\n\nCommands:\n  doctor          Check installation and optional GPU availability\n  init            Create a portable experiment from a CDS FASTA\n  run             Run a custom YAML experiment\n  score           Score a reference construct\n  validate        Validate a YAML experiment\n  preview         Draw construct and scorer scopes\n  report          Regenerate an individual run report\n  pilot           Estimate experiment runtime\n  paper           Verify, analyze, plot, or rerun the paper\n  assets          List, install, or verify ready model assets\n  scorers         Inspect scorer interfaces and generate YAML templates\n\nUse morse COMMAND --help for command options.")
        return 0
    command, rest = argv[0], argv[1:]
    if command in {"init", "run", "score", "validate", "preview", "report", "pilot"}:
        from morse.cli import main as core_main
        return core_main([command, *rest]) or 0
    if command == "paper":
        _require_checkout_resources(root, command)
        from paper.cli import main as paper_main
        return paper_main(rest) or 0
    if command == "assets":
        from assets.manage import main as assets_main
        return assets_main(rest) or 0
    if command == "doctor":
        return doctor(rest)
    if command == "scorers":
        from morse_tools.scorers import main as scorers_main
        return scorers_main(rest, root=root)
    raise SystemExit(f"Unknown command {command!r}. Run morse --help.")


if __name__ == "__main__":
    raise SystemExit(main())
