"""Inspect scorer manifests without importing models or executing scorer code."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

import yaml

from morse.scorer_manifest import load_scorer_manifest, request_parameter_specs
from morse.schema4 import scope_expression, scope_from_internal


def _catalog(root: Path) -> dict[str, tuple[Path, dict[str, Any]]]:
    catalog = {}
    for path in sorted((root / "scorers").glob("*/morse-scorer*.yaml")):
        manifest = load_scorer_manifest(path)
        scorer_id = manifest["id"]
        if scorer_id in catalog:
            raise ValueError(f"Duplicate scorer id in catalog: {scorer_id}")
        catalog[scorer_id] = path, manifest
    if not catalog:
        raise ValueError(
            "The scorer catalog requires the complete monorepo checkout. "
            "Set MORSE_ROOT to its path, or use --manifest /path/to/morse-scorer.yaml."
        )
    return catalog


def _metrics(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    result = {}
    for name, output in manifest["outputs"].items():
        metric = {
            key: copy.deepcopy(value)
            for key, value in output["default_policy"].items()
            if key in {"role", "direction", "min", "max", "equals"}
        }
        if metric["role"] == "objective":
            metric["role"] = "optimize"
        if output["type"] != "float":
            metric["value_type"] = output["type"]
        result[name] = metric
    return result


def scorer_template(manifest: dict[str, Any], command: str | None) -> dict[str, Any]:
    """Return a pasteable section using the scorer's declared default policies."""
    if command is not None and not command.strip():
        raise ValueError("--command must be a non-empty Bash command")
    section = {
        "name": manifest["id"],
        "command": command if command is not None else "REPLACE_WITH_YOUR_SCORER_COMMAND",
    }
    defaults = {key: copy.deepcopy(spec["default"])
                for key, spec in request_parameter_specs(manifest).items() if "default" in spec}
    if defaults:
        section["parameters"] = defaults
    section.update({
        "timeout_seconds": manifest["runtime"]["timeout_seconds"],
        "scope": scope_from_internal(manifest["scope"]["default"]),
        "metrics": _metrics(manifest),
    })
    return {"scorers": [section]}


def scorer_interface(manifest: dict[str, Any]) -> dict[str, Any]:
    scope = copy.deepcopy(manifest["scope"])
    scope["default"] = scope_from_internal(scope["default"])
    scope["allowed_parts"] = ["cds_body" if part == "cds" else part for part in scope.get("allowed_parts", [])]
    metrics = {}
    for name, output in manifest["outputs"].items():
        metrics[name] = {
            "value_type": output["type"],
            "description": output["description"],
            "default": _metrics(manifest)[name],
        }
        for field in ("unit", "value_range"):
            if field in output:
                metrics[name][field] = output[field]
    parameter_specs = request_parameter_specs(manifest)
    defaults = {key: copy.deepcopy(spec["default"])
                for key, spec in parameter_specs.items() if "default" in spec}
    request_example = {
        "request_id": 1,
        "sequences": [manifest["self_test"]["sequence"]],
        "sequence_contexts": [manifest["self_test"]["sequence_context"]],
    }
    if defaults:
        request_example["parameters"] = defaults
    return {
        "id": manifest["id"], "name": manifest["name"], "version": manifest["version"],
        "description": manifest["description"], "scope": scope,
        "metrics": metrics,
        "parameters": parameter_specs,
        "protocol": {
            "stdin": "One JSON request followed by EOF per batch",
            "stdout": "One JSON object with the same request_id and an outputs list in input order",
            "stderr": "Diagnostics",
            "parameters": "Optional JSON object containing scorer-specific run settings",
        },
        "request_example": request_example,
    }


def main(argv: list[str], *, root: Path) -> int:
    parser = argparse.ArgumentParser(prog="morse scorers", description="Inspect scorer interfaces and print YAML sections for prepared scorers.")
    sub = parser.add_subparsers(dest="action")
    sub.add_parser("list", help="list the local scorer catalog")
    server = sub.add_parser("serve", help="optionally expose one persistent JSONL scorer over local HTTP before a run")
    server.add_argument("--command", required=True, help="Bash command for one existing JSONL scorer")
    server.add_argument("--port", type=int, required=True)
    server.add_argument("--timeout", type=float, default=1800, help="scoring timeout in seconds")
    for action in ("describe", "template"):
        item = sub.add_parser(action, help="show metrics and input scope" if action == "describe" else "print a scorer YAML section")
        item.add_argument("scorer", nargs="?", help="scorer id from the local catalog")
        item.add_argument("--manifest", type=Path, help="use a custom morse-scorer.yaml instead of a catalog id")
        if action == "describe":
            item.add_argument("--json", action="store_true", help="print the interface as JSON")
        else:
            item.add_argument("--command", help="Bash command to an already prepared script, service or container")
    args = parser.parse_args(argv)
    try:
        if args.action == "serve":
            from morse_tools.scorer_server import serve
            return serve(args.command, args.port, timeout_seconds=args.timeout)
        if args.action in {None, "list"}:
            for scorer_id, (path, manifest) in _catalog(root).items():
                print(f"{scorer_id}\t{manifest['version']}\t{path.relative_to(root)}")
            return 0
        if bool(args.scorer) == bool(args.manifest):
            parser.error("provide either a scorer id or --manifest PATH")
        if args.manifest:
            manifest = load_scorer_manifest(args.manifest.expanduser())
        else:
            catalog = _catalog(root)
            if args.scorer not in catalog:
                parser.error(f"unknown scorer {args.scorer!r}; use 'morse scorers' to list ids")
            _, manifest = catalog[args.scorer]
        if args.action == "template":
            if args.command is None:
                print("# Set command to your prepared scorer, or pass --command to this CLI.")
            print(yaml.safe_dump(scorer_template(manifest, args.command), sort_keys=False), end="")
        else:
            interface = scorer_interface(manifest)
            if args.json:
                print(json.dumps(interface, indent=2, ensure_ascii=False))
            else:
                print(f"{manifest['name']} ({manifest['id']} {manifest['version']})")
                print(manifest["description"])
                print(f"\nInput scope: {scope_expression(interface['scope']['default'])}")
                print("\nMetrics:")
                for name, metric in interface["metrics"].items():
                    policy = metric["default"]
                    label = policy["role"]
                    if policy.get("direction"):
                        label += " " + policy["direction"]
                    print(f"  {name} [{metric['value_type']}; {label}]")
                    print(f"    {metric['description']}")
                    if metric.get("value_range"):
                        print(f"    Range: {metric['value_range']}")
                print("\nProtocol: JSON request on stdin, then EOF; JSON response on stdout.")
                print("Response: same request_id; outputs has one metric object per sequence, in input order.")
                if interface["parameters"]:
                    print("\nParameters:")
                    for name, spec in interface["parameters"].items():
                        default = f"; default: {json.dumps(spec['default'])}" if "default" in spec else ""
                        print(f"  {name} [{spec['type']}{default}]")
                        print(f"    {spec['description']}")
                    print("Set these in the scorer's YAML parameters section. Omitted values retain the scorer's startup settings.")
                print("Scorer dependencies must already be installed.")
        return 0
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2
