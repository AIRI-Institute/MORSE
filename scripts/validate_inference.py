#!/usr/bin/env python3
"""Compare an image's inference adapter against the current source overlay.

Models and Docker images are supplied locally; no model training or downloads are
performed. Only this command's uniquely named containers are removed on timeout.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import uuid


MODELS = {
    "ribonn": {"folder": "ribonn", "asset": "ribonn-human-ensemble", "env": {"RIBONN_MODEL_DIR": "/models", "RIBONN_DEVICE": "cuda"}, "python": "python", "server": "server.py"},
    "orthrus-te-full-fit": {"folder": "translation-efficiency", "asset": "orthrus-te-full-fit", "env": {"ORTHRUS_CHECKPOINT": "/models/model_final.pt", "ORTHRUS_MODEL_SOURCE": "/backbone", "ORTHRUS_DEVICE": "cuda"}, "python": "python", "server": "server.py"},
    "orthrus-half-life": {"folder": "half-life", "asset": "orthrus-half-life", "env": {"ORTHRUS_CHECKPOINT": "/models/model_final.pt", "ORTHRUS_MODEL_SOURCE": "/backbone", "ORTHRUS_DEVICE": "cuda"}, "python": "python", "server": "server.py"},
    "rinalmo-imrnads": {"folder": "immunogenicity", "asset": "rinalmo-imrnads", "env": {"RINALMO_CHECKPOINT": "/models/best_model.pt", "RINALMO_DEVICE": "cuda", "RINALMO_BATCH_SIZE": "64"}, "python": "python", "server": "server.py"},
}


def evaluate(root, assets, model, image, request, gpu, overlay, output, timeout, ribonn_source=None):
    item = MODELS[model]
    name = "morse-inference-check-" + uuid.uuid4().hex[:12]
    command = ["docker", "run", "--rm", "--interactive", "--name", name, "--label", "morse.validation=inference", "--network=none", "--gpus", f"device={gpu}", "--cpus=2", "--memory=10g", "--pids-limit=256", "--volume", f"{assets / item['asset']}:/models:ro"]
    if model.startswith("orthrus"):
        command += ["--volume", f"{assets / 'orthrus-backbone'}:/backbone:ro"]
    for key, value in {"OMP_NUM_THREADS": "2", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "2", "HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1", **item["env"]}.items():
        command += ["--env", f"{key}={value}"]
    if overlay:
        command += ["--volume", f"{root / 'scorers' / item['folder']}:/app:ro"]
    if model == "ribonn":
        source = ribonn_source or root / "scorers/ribonn/third_party/ribonn"
        command += ["--mount", f"type=bind,src={source},dst=/app/third_party/ribonn,readonly"]
    command += ["--entrypoint", item["python"], image, "-u", f"/app/{item['server']}"]
    started = time.monotonic()
    try:
        completed = subprocess.run(command, input=json.dumps(request) + "\n", text=True, capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        subprocess.run(["docker", "rm", "-f", name], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        raise
    label = "current" if overlay else "archived"
    (output / f"{model}-{label}.stderr.log").write_text(completed.stderr)
    (output / f"{model}-{label}.stdout.log").write_text(completed.stdout)
    if completed.returncode:
        raise RuntimeError(f"{model}/{label}: container exited {completed.returncode}; inspect logs")
    lines = [line for line in completed.stdout.splitlines() if line.strip().startswith("{")]
    if len(lines) != 1:
        raise RuntimeError(f"{model}/{label}: expected exactly one JSONL response")
    response = json.loads(lines[0])
    if "error" in response:
        raise RuntimeError(f"{model}/{label}: {response['error']}")
    if response.get("request_id") != request["request_id"] or len(response.get("outputs", [])) != len(request["sequences"]):
        raise RuntimeError(f"{model}/{label}: malformed output contract")
    return {"elapsed_seconds": time.monotonic() - started, "response": response}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=MODELS, required=True)
    parser.add_argument("--image", required=True)
    parser.add_argument("--assets", type=Path, required=True)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gpu", type=int, default=0)
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--atol", type=float, default=1e-6)
    parser.add_argument("--ribonn-source", type=Path,
                        help="Exact external inference source (default prepared source asset, then vendored extraction)")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    args.output.mkdir(parents=True, exist_ok=True)
    fixture = json.loads(args.request.read_text())
    request = fixture.get("request", fixture)
    ribonn_source, verified_source = None, None
    if args.model == "ribonn":
        sys.path.insert(0, str(root))
        from assets.manage import REGISTRY, verify_asset
        source_asset = next(item for item in json.loads(REGISTRY.read_text())["assets"]
                            if item["id"] == "ribonn-inference-source")
        prepared = args.assets.resolve() / "ribonn-inference-source"
        ribonn_source = (args.ribonn_source or (prepared if prepared.is_dir()
                          else root / "scorers/ribonn/third_party/ribonn")).resolve()
        verified_source = verify_asset(source_asset, ribonn_source)
    archived = evaluate(root, args.assets.resolve(), args.model, args.image, request, args.gpu, False, args.output, args.timeout, ribonn_source)
    current = evaluate(root, args.assets.resolve(), args.model, args.image, request, args.gpu, True, args.output, args.timeout, ribonn_source)
    differences = []
    for before, after in zip(archived["response"]["outputs"], current["response"]["outputs"], strict=True):
        if set(before) != set(after):
            raise RuntimeError("Output fields changed")
        for key in before:
            if isinstance(before[key], (float, int)) and not isinstance(before[key], bool):
                difference = abs(float(before[key]) - float(after[key]))
                if not __import__("math").isfinite(difference):
                    raise RuntimeError("Nonfinite output")
                differences.append({"field": key, "absolute_difference": difference})
            elif before[key] != after[key]:
                raise RuntimeError(f"Output metadata changed: {key}")
    maximum = max((x["absolute_difference"] for x in differences), default=0.0)
    revision = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    receipt = {"model": args.model, "source_revision": revision, "image": args.image, "request_sha256": hashlib.sha256(args.request.read_bytes()).hexdigest(), "sequences": len(request["sequences"]), "archived": archived, "current": current, "maximum_absolute_difference": maximum, "tolerance": args.atol, "passed": maximum <= args.atol}
    receipt["comparison"] = "Image adapter versus current source overlay in the same supplied image; not an assertion of historical-runtime identity"
    if verified_source is not None:
        receipt["verified_external_inference_source"] = verified_source
    (args.output / f"{args.model}-parity.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({k: receipt[k] for k in ["model", "sequences", "maximum_absolute_difference", "tolerance", "passed"]}))
    return 0 if receipt["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
