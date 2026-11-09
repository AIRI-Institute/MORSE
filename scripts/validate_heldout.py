#!/usr/bin/env python3
"""Compare the current held-out TE adapter with frozen paper observations.

Uses a local compatible Docker image and ready checkpoints. It performs no
downloads, fitting, commits or remote operations. Numerical discrepancies are
reported against the supplied frozen observations without adjusting tolerance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import uuid

CHECKPOINT_SHA = "70913fc29ba9c073100913b56f59542520ec18a473fae358ff19bd57ecddd604"
FIELD = "translation_efficiency_score"
MODEL = "orthrus-te-held-out"


def sha256(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def fixture_request(fixture: dict) -> dict:
    request = fixture.get("request")
    if not isinstance(request, dict):
        raise ValueError("Fixture must contain a request object")
    sequences = request.get("sequences")
    contexts = request.get("sequence_contexts")
    observations = fixture.get("archived_observations")
    if not isinstance(sequences, list) or not sequences or not all(isinstance(s, str) and s for s in sequences):
        raise ValueError("Fixture must contain nonempty sequence strings")
    if not isinstance(contexts, list) or len(contexts) != len(sequences) or not all(isinstance(c, dict) for c in contexts):
        raise ValueError("Fixture requires one explicit context per sequence")
    if not isinstance(observations, list) or len(observations) != len(sequences):
        raise ValueError("Fixture requires one archived observation per sequence")
    if "request_id" not in request:
        raise ValueError("Fixture request must carry request_id")
    for observation in observations:
        if not isinstance(observation, dict) or not isinstance(observation.get("raw"), dict):
            raise ValueError("Malformed archived observation")
        value = observation["raw"].get(FIELD)
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
            raise ValueError("Archived TE observations must be finite numbers")
    return request


def compare_outputs(fixture: dict, response: dict, atol: float) -> list[dict]:
    request = fixture_request(fixture)
    outputs = response.get("outputs")
    if response.get("request_id") != request["request_id"]:
        raise ValueError("Response request_id does not match the fixture")
    if "error" in response:
        raise ValueError(f"Scorer returned an error: {response['error']}")
    if not isinstance(outputs, list) or len(outputs) != len(request["sequences"]):
        raise ValueError("Response must contain exactly one output per input sequence")
    differences = []
    for index, (observation, output) in enumerate(zip(fixture["archived_observations"], outputs)):
        if not isinstance(output, dict) or set(output) != {FIELD}:
            raise ValueError("Held-out scorer returned unexpected output fields")
        actual = output[FIELD]
        if isinstance(actual, bool) or not isinstance(actual, (int, float)) or not math.isfinite(actual):
            raise ValueError("Scorer returned a nonfinite or nonnumeric TE score")
        expected = float(observation["raw"][FIELD])
        signed = float(actual) - expected
        absolute = abs(signed)
        differences.append({
            "index": index, "candidate_id": observation.get("candidate_id"),
            "sequence_sha256": hashlib.sha256(request["sequences"][index].encode()).hexdigest(),
            "field": FIELD, "archived_value": expected, "current_value": float(actual),
            "signed_difference": signed, "absolute_difference": absolute,
            "relative_difference": absolute / abs(expected) if expected else None,
            "within_absolute_tolerance": absolute <= atol,
        })
    return differences


def _text(value: str | bytes | None) -> str:
    return value.decode(errors="replace") if isinstance(value, bytes) else value or ""


def request_plan(request: dict, mode: str, batch_size: int) -> list[dict]:
    count = len(request["sequences"])
    if mode == "combined":
        groups = [list(range(count))]
    elif mode == "baseline-and-candidates":
        groups = [[0]] + ([list(range(1, count))] if count > 1 else [])
    elif mode == "individual":
        groups = [[index] for index in range(count)]
    else:
        raise ValueError(f"Unknown request mode: {mode}")
    plan = []
    for query_index, indices in enumerate(groups):
        request_id = request["request_id"] if mode == "combined" else f"{request['request_id']}::{mode}::{query_index}"
        query = {"request_id": request_id,
                 "sequences": [request["sequences"][index] for index in indices],
                 "sequence_contexts": [request["sequence_contexts"][index] for index in indices]}
        batches = []
        for offset in range(0, len(indices), batch_size):
            batch_indices = indices[offset:offset + batch_size]
            lengths = [len(request["sequences"][index]) for index in batch_indices]
            maximum = max(lengths)
            batches.append({"fixture_indices": batch_indices, "sequence_lengths": lengths,
                            "collated_tensor_shape": [len(lengths), maximum, 6],
                            "right_padding_nt": [maximum - length for length in lengths]})
        plan.append({"request": query, "fixture_indices": indices, "adapter_batches": batches})
    return plan


def merge_responses(request: dict, plan: list[dict], responses: list[dict]) -> dict:
    if len(plan) != len(responses):
        raise ValueError("Expected one JSONL response per physical request")
    outputs = [None] * len(request["sequences"])
    for item, response in zip(plan, responses):
        query = item["request"]
        if not isinstance(response, dict) or response.get("request_id") != query["request_id"]:
            raise ValueError("Physical response request_id does not match its request")
        if "error" in response:
            raise ValueError(f"Scorer returned an error: {response['error']}")
        values = response.get("outputs")
        if not isinstance(values, list) or len(values) != len(item["fixture_indices"]):
            raise ValueError("Physical request output count changed")
        for index, value in zip(item["fixture_indices"], values):
            outputs[index] = value
    if any(value is None for value in outputs):
        raise ValueError("A fixture input has no physical response")
    return {"request_id": request["request_id"], "outputs": outputs}


def evaluate(root: Path, assets: Path, image: str, request: dict, gpu: int,
             batch_size: int, output: Path, timeout: int, request_mode: str = "combined") -> dict:
    plan = request_plan(request, request_mode, batch_size)
    name = "morse-heldout-check-" + uuid.uuid4().hex[:12]
    command = [
        "docker", "run", "--rm", "--interactive", "--name", name,
        "--label", "morse.validation=heldout-paper-observations", "--network=none",
        "--gpus", f"device={gpu}", "--cpus=2", "--memory=10g", "--pids-limit=256",
        "--volume", f"{assets / MODEL}:/models:ro",
        "--volume", f"{assets / 'orthrus-backbone'}:/backbone:ro",
        "--volume", f"{root / 'scorers/translation-efficiency'}:/app:ro",
        "--workdir", "/app",
    ]
    env = {
        "OMP_NUM_THREADS": "2", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "2",
        "HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1", "PYTHONUNBUFFERED": "1",
        "ORTHRUS_CHECKPOINT": "/models/best.pt", "ORTHRUS_MODEL_SOURCE": "/backbone",
        "ORTHRUS_DEVICE": "cuda", "ORTHRUS_BATCH_SIZE": str(batch_size),
        "ORTHRUS_PORTABLE": "0", "MORSE_EXPECTED_CHECKPOINT_SHA256": CHECKPOINT_SHA,
    }
    for key, value in env.items():
        command += ["--env", f"{key}={value}"]
    command += ["--entrypoint", "python", image, "-u", "/app/server_heldout.py"]
    started = time.monotonic()
    stdout, stderr = "", ""
    try:
        payload = "".join(json.dumps(item["request"], allow_nan=False) + "\n" for item in plan)
        completed = subprocess.run(command, input=payload,
                                   text=True, capture_output=True, timeout=timeout)
        stdout, stderr = completed.stdout, completed.stderr
        if completed.returncode:
            raise RuntimeError(f"Held-out container exited {completed.returncode}; inspect saved logs")
    except subprocess.TimeoutExpired as error:
        stdout, stderr = _text(error.stdout), _text(error.stderr)
        raise RuntimeError(f"Held-out inference exceeded the declared {timeout}-second timeout") from error
    finally:
        # The name belongs only to this invocation; do not inspect/remove other jobs.
        try:
            subprocess.run(["docker", "rm", "-f", name], capture_output=True, timeout=15)
        except (OSError, subprocess.TimeoutExpired):
            pass
        (output / f"{MODEL}-current.stdout.log").write_text(stdout)
        (output / f"{MODEL}-current.stderr.log").write_text(stderr)
    lines = [line for line in stdout.splitlines() if line.strip().startswith("{")]
    if len(lines) != len(plan):
        raise ValueError("Expected one JSONL response per physical request; inspect saved stdout")
    responses = [json.loads(line) for line in lines]
    response = merge_responses(request, plan, responses)
    return {"elapsed_seconds": time.monotonic() - started, "response": response,
            "raw_responses": responses, "request_mode": request_mode,
            "physical_requests": [{"request_id": item["request"]["request_id"],
                                   "fixture_indices": item["fixture_indices"],
                                   "adapter_batches": item["adapter_batches"]} for item in plan],
            "shape_provenance": "Inferred from the frozen adapter's contiguous batching and zero-padding collate contract; not a GPU profiler trace",
            "container_name": name, "adapter": "/app/server_heldout.py", "runtime_env": env}


def main(argv: list[str] | None = None) -> int:
    root = Path(os.environ.get("MORSE_ROOT", str(Path(__file__).resolve().parents[1]))).expanduser().resolve()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True, help="Existing local image with compatible Orthrus/Mamba dependencies")
    parser.add_argument("--assets", required=True, type=Path, help="Ready asset cache root")
    parser.add_argument("--request", type=Path, default=root / "validation/fixtures/fig5_inference.json")
    parser.add_argument("--output", required=True, type=Path, help="Directory for an immutable validation receipt and logs")
    parser.add_argument("--gpu", type=int, default=0)
    parser.add_argument("--batch-size", type=int, default=4, help="Frozen adapter default; always recorded in the receipt")
    parser.add_argument("--request-mode", choices=["combined", "baseline-and-candidates", "individual"], default="combined",
                        help="combined preserves the diagnostic default; baseline-and-candidates sends the Core reference alone then candidates; individual sends each input separately in one persistent process")
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--atol", type=float, default=1e-4, help="Absolute tolerance fixed before execution; default 1e-4")
    args = parser.parse_args(argv)
    if args.gpu < 0 or args.batch_size < 1 or args.timeout < 1 or not math.isfinite(args.atol) or args.atol < 0:
        parser.error("GPU index/tolerance must be nonnegative; batch size/timeout must be positive")
    output = args.output.expanduser().resolve()
    receipt_path = output / f"{MODEL}-parity.json"
    if any((output / name).exists() for name in [receipt_path.name, f"{MODEL}-current.stdout.log", f"{MODEL}-current.stderr.log"]):
        parser.error("Refusing to replace existing validation evidence; choose a new --output directory")
    output.mkdir(parents=True, exist_ok=True)
    receipt = {"model": MODEL, "comparison": "current-source overlay versus archived paper raw observations",
               "image": args.image, "tolerance": {"kind": "absolute", "value": args.atol},
               "batch_size": args.batch_size, "gpu_index": args.gpu, "timeout_seconds": args.timeout,
               "request_mode": args.request_mode,
               "request_boundary_note": "Core evaluates its reference construct alone before search (optimizer.py score_constructs([original_construct])). Combined and individual modes are explicit batch-shape diagnostics; no model or encoding is modified.",
               "source_overlay": True, "portable_fullfit_encoder_used": False, "passed": False,
               "numerical_note": "Different batch shapes, GPU kernels and BF16 arithmetic can change scores; discrepancies remain recorded and never cause automatic tolerance relaxation."}
    try:
        fixture = json.loads(args.request.read_text())
        request = fixture_request(fixture)
        receipt.update(request_sha256=sha256(args.request), sequences=len(request["sequences"]),
                       archived_run_id=fixture.get("run_id"))
        sys.path.insert(0, str(root))
        from assets.manage import REGISTRY, verify_asset
        registry = json.loads(REGISTRY.read_text())
        asset = next(item for item in registry["assets"] if item["id"] == MODEL)
        verified = verify_asset(asset, args.assets.expanduser().resolve() / MODEL)
        if verified["files"][0]["sha256"] != CHECKPOINT_SHA:
            raise ValueError("Registry does not identify the frozen held-out checkpoint")
        receipt["checkpoint_sha256"] = CHECKPOINT_SHA
        receipt["verified_backbone_files"] = verified.get("dependencies", [{}])[0].get("files", [])
        receipt["source_sha256"] = {name: sha256(root / "scorers/translation-efficiency" / name)
                                    for name in ["server_heldout.py", "heldout_encoding.py", "model.py"]}
        receipt["source_revision"] = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True, timeout=15).strip()
        receipt["resolved_image_id"] = subprocess.check_output(["docker", "image", "inspect", "--format", "{{.Id}}", args.image], text=True, timeout=15).strip()
        current = evaluate(root, args.assets.expanduser().resolve(), args.image, request,
                           args.gpu, args.batch_size, output, args.timeout, args.request_mode)
        differences = compare_outputs(fixture, current["response"], args.atol)
        maximum = max(item["absolute_difference"] for item in differences)
        receipt.update(status="completed", current=current, differences=differences,
                       maximum_absolute_difference=maximum,
                       passed=all(item["within_absolute_tolerance"] for item in differences))
    except (OSError, ValueError, KeyError, StopIteration, RuntimeError, subprocess.SubprocessError) as error:
        receipt.update(status="execution_error", error=f"{type(error).__name__}: {error}")
    receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: receipt.get(key) for key in ["model", "status", "request_mode", "sequences", "maximum_absolute_difference", "tolerance", "passed", "error"]}))
    if receipt["status"] == "execution_error":
        return 2
    return 0 if receipt["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
