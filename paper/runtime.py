"""Historical auxiliary scorer settings and the separate native VaxPress baseline.

MORSE experiments run directly from YAML; this module does not launch or wrap them.
"""
from __future__ import annotations

import copy
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import uuid

from .io import DATA, HERE, REPO, read_json, registry, sha256, write_json


def asset_root(value=None):
    return Path(value or os.environ.get("MORSE_ASSETS_DIR", REPO / ".cache/morse/assets")).expanduser().resolve()


def _device(device):
    if device == "cpu":
        return None
    if device.startswith("cuda:") and device[5:].isdigit():
        return int(device[5:])
    raise ValueError("Device must be cpu or cuda:N")


def _expand(value):
    if isinstance(value, str):
        return value.replace("@paper/", str(HERE) + "/").replace("@repo/", str(REPO) + "/").replace("@python", sys.executable)
    return value


def _image_overrides(path):
    """Accept an explicit legacy image map or the published runtime catalog."""
    if path is None:
        return {}, False
    data = read_json(path)
    if not isinstance(data, dict):
        raise ValueError("Image map must be a JSON object")
    if "published_images" not in data:
        return data, False
    if not isinstance(data["published_images"], list):
        raise ValueError("published_images must be a list")
    aliases = {"translation-efficiency": "orthrus-te-full-fit",
               "translation-efficiency-held-out": "orthrus-te-held-out"}
    images = {}
    for item in data["published_images"]:
        reference = item.get("immutable_reference")
        if not isinstance(reference, str) or not re.fullmatch(r"ghcr\.io/[A-Za-z0-9_./-]+@sha256:[0-9a-f]{64}", reference):
            raise ValueError("Published image catalog requires immutable GHCR references")
        identifier = item["scorer"]
        images[aliases.get(identifier, identifier)] = reference
    return images, True


def _uses_baked_assets(args, logical):
    """Only an explicitly capable published image may bypass host asset setup."""
    if args.runtime != "container" or args.asset_root is not None or args.image_map is None:
        return False
    data = read_json(args.image_map)
    if not isinstance(data, dict):
        return False
    aliases = {"translation-efficiency": "orthrus-te-full-fit",
               "translation-efficiency-held-out": "orthrus-te-held-out"}
    entries = data.get("published_images", [])
    if not isinstance(entries, list):
        return False
    selected = None
    for item in entries:
        if isinstance(item, dict) and aliases.get(item.get("scorer"), item.get("scorer")) == logical:
            selected = item
    # Match the same final override selected by _image_overrides.
    return selected is not None and selected.get("self_contained") is True


def _overlay_source_files(logical):
    """Overlay adapter files while preserving assets already present in /app."""
    folder = "translation-efficiency" if logical.startswith("orthrus-te") else logical
    names = {
        "orthrus-te-full-fit": ("model.py", "server.py", "portable_orthrus.py"),
        "orthrus-te-held-out": ("model.py", "server.py", "portable_orthrus.py", "server_heldout.py", "heldout_encoding.py"),
        "half-life": ("model.py", "server.py", "portable_orthrus.py"),
        "immunogenicity": ("server.py", "runtime_model.py"),
        "ribonn": ("encoding.py", "model.py", "server.py"),
        "ribodecode": ("server.py",),
    }[logical]
    return [REPO / "scorers" / folder / name for name in names]


def _source_file_mounts(logical):
    return [argument for path in _overlay_source_files(logical)
            for argument in ("--volume", f"{path}:/app/{path.name}:ro")]


def _overlay_asset_mounts(logical, assets):
    """Plan current-source mounts without installing or loading any assets."""
    if logical == "ribonn":
        source = REPO / "scorers/ribonn/third_party/ribonn"
        return ["--mount", f"type=bind,src={source},dst=/app/third_party/ribonn,readonly"], {}
    if logical == "ribodecode":
        package = assets / "ribodecode-translationmodel"
        return ["--mount", f"type=bind,src={package / 'package'},dst=/opt/ribodecode,readonly",
                "--mount", f"type=bind,src={package / 'env_file.csv'},dst=/assets/env_file.csv,readonly"], {
                    "PYTHONPATH": "/opt/ribodecode", "RIBODECODE_ENV_TEMPLATE": "/assets/env_file.csv"}
    return [], {}


def _resolve_scorer(scorer, args, tag):
    scorer = copy.deepcopy(scorer)
    runtime = scorer.pop("paper_runtime")
    logical = runtime["scorer"]
    assets = asset_root(args.asset_root)
    baked = _uses_baked_assets(args, logical)
    env = {k: _expand(v) for k, v in scorer.get("env", {}).items()}
    reference_device = runtime.get("reference_device")
    if reference_device not in {None, "cpu", "cuda"}:
        raise ValueError(f"Unknown reference device for {logical}: {reference_device}")
    # --device selects the available GPU, never promotes a historical CPU
    # scorer onto CUDA. RiboDecode otherwise auto-selects a visible GPU.
    gpu = None if reference_device == "cpu" else _device(args.device)
    local_device = "cpu" if gpu is None else (args.device if args.runtime == "process" else "cuda:0")
    env.pop("CUDA_VISIBLE_DEVICES", None)
    if gpu is None:
        env["CUDA_VISIBLE_DEVICES"] = "-1"
    if logical.startswith("orthrus-te"):
        checkpoint = "best.pt" if logical.endswith("held-out") else "model_final.pt"
        env.update(ORTHRUS_DEVICE=local_device, ORTHRUS_MODEL_SOURCE=str(assets / "orthrus-backbone"),
                   ORTHRUS_CHECKPOINT=str(assets / logical / checkpoint),
                   ORTHRUS_PORTABLE="1" if args.runtime=="process" or gpu is None else "0")
    elif logical == "ribonn":
        env.update(RIBONN_DEVICE=local_device, RIBONN_MODEL_DIR=str(assets / "ribonn-human-ensemble"))
    elif logical == "ribodecode":
        env["RIBODECODE_DEVICE"] = local_device
    if args.runtime == "process" or logical in {"csi", "cai"}:
        folder = {"orthrus-te-full-fit": "translation-efficiency", "orthrus-te-held-out": "translation-efficiency"}.get(logical, logical)
        server = "server_heldout.py" if logical == "orthrus-te-held-out" else "scorer.py" if logical == "viennarna" else "server.py"
        scorer.update(mode="persistent-process", command=[sys.executable, server], workdir=str(REPO / "scorers" / folder))
        if logical == "manufacturability":
            scorer["command"] = [sys.executable, "-m", "manufacturability_scorer.feasibility_server"]
        scorer.pop("image", None)
        scorer.pop("docker_args", None)
    else:
        overrides, published = _image_overrides(args.image_map)
        if published and logical not in overrides:
            raise ValueError(f"No published image for {logical} in the supplied runtime lock")
        image = overrides.get(logical, runtime["historical_image"])
        if not image:
            raise ValueError(f"No Docker image configured for {logical}")
        scorer.update(mode="persistent-container", image=image)
        docker = ["--name", f"morse-paper-{tag}-{logical}", "--network=none", "--pids-limit=256", "--log-driver=none"]
        if gpu is not None and logical in {"ribonn", "ribodecode", "orthrus-te-full-fit", "orthrus-te-held-out"}:
            docker += ["--gpus", f"device={gpu}"]
        if baked and logical.startswith("orthrus-te"):
            checkpoint = "best.pt" if logical.endswith("held-out") else "model_final.pt"
            env.update(ORTHRUS_CHECKPOINT="/models/" + checkpoint, ORTHRUS_MODEL_SOURCE="/backbone")
        elif logical.startswith("orthrus-te"):
            if logical.endswith("held-out"):
                docker += ["--volume", f"{assets / logical / 'best.pt'}:/te/best.pt:ro"]
                env["ORTHRUS_CHECKPOINT"] = "/te/best.pt"
            else:
                docker += ["--volume", f"{assets / logical}:/models:ro"]
                env["ORTHRUS_CHECKPOINT"] = "/models/model_final.pt"
            docker += ["--volume", f"{assets / 'orthrus-backbone'}:/backbone:ro"]
            env["ORTHRUS_MODEL_SOURCE"] = "/backbone"
        elif logical == "ribonn":
            if not baked:
                docker += ["--volume", f"{assets / 'ribonn-human-ensemble'}:/models:ro"]
            env["RIBONN_MODEL_DIR"] = "/models"
        if args.source_overlay:
            folder={"orthrus-te-full-fit":"translation-efficiency","orthrus-te-held-out":"translation-efficiency"}.get(logical,logical)
            if baked:
                docker.extend(_source_file_mounts(logical))
            else:
                docker += ["--volume", f"{REPO/'scorers'/folder}:/app:ro"]
            docker += ["--volume",f"{HERE/'adapters/container_entrypoints'}:/morse-entrypoints:ro","--entrypoint",f"/morse-entrypoints/{logical}"]
            if not baked:
                mounts, current_env = _overlay_asset_mounts(logical, assets)
                docker.extend(mounts)
                env.update(current_env)
        scorer["docker_args"] = docker
        scorer.pop("command", None)
        scorer.pop("workdir", None)
    for key, value in env.items():
        if isinstance(value, str) and value.startswith("@input/"):
            name = value.split("/", 1)[1]
            replacement = HERE / "inputs" / name
            if not replacement.is_file():
                raise ValueError(f"Unresolved scorer input {logical}:{key}={name}")
            env[key] = str(replacement)
    scorer["env"] = env
    return scorer


def vaxpress(record, args):
    cfg = read_json(HERE / "inputs/vaxpress_protocol.json")
    for key, relative in [("csi_model_path", "csi/model.py"), ("csi_profile_path", "csi/profiles/homo_sapiens_consensus_wu2019_v1.json"), ("viennarna_model_path", "viennarna/scorer.py")]:
        cfg[key] = str(REPO / "scorers" / relative)
    output = args.output.resolve() / record["run_id"]
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Refusing to replace existing experiment evidence: {output}; choose a new --output directory")
    output.mkdir(parents=True, exist_ok=True)
    path = output / "vaxpress_protocol.json"
    write_json(path, cfg)
    commands = []
    overrides, published = _image_overrides(args.image_map)
    if published and "vaxpress" not in overrides:
        raise ValueError("The published scorer catalog excludes the external VaxPress baseline; supply its prepared image map")
    image = overrides.get("vaxpress", cfg["vaxpress_image"])
    total_events=12800 if args.scoring_events is None else args.scoring_events
    if total_events<=0 or total_events%5:raise ValueError("--scoring-events must be positive and divisible by five weights")
    per_weight=total_events//5
    for alpha in cfg["alphas"]:
        target = output / f"alpha_{alpha}"
        command = [sys.executable, str(HERE / "adapters/run_vaxpress_shared.py"), "--config", str(path), "--output", str(target),
                   "--seed", str(record["seed"]), "--alpha", str(alpha), "--budget", str(per_weight), "--processes", "4"]
        if args.runtime == "container":
            command = ["docker", "run", "--rm", "--network=none", "--volume", f"{REPO}:{REPO}:ro", "--volume", f"{output}:{output}:rw", "--entrypoint", "python3", image, *command[1:]]
        commands.append(command)
    receipt = {"run_id": record["run_id"], "commands": commands,"protocol_sha256":sha256(path),"mode":"diagnostic-scoring-event-override" if args.scoring_events else "paper-protocol",
               "total_scoring_events":total_events,"events_per_weight":per_weight,"budget": f"5 weights * {per_weight} scoring events = {total_events} per independent sweep", "upstream_commit": cfg["vaxpress_commit"]}
    write_json(output / "plan.json", receipt)
    if args.execute:
        for command in commands:
            subprocess.run(command, check=True)
        receipt["status"] = "completed"
        write_json(output / "execution.json", receipt)
    return receipt


def run_baseline(args):
    """Run the separate native VaxPress benchmark, never a MORSE optimizer wrapper."""
    selected = [r for r in registry()["runs"] if r["condition"] == "vaxpress" and (args.seed is None or r["seed"] == args.seed)]
    receipts = [vaxpress(record, args) for record in selected]
    result = {"executed": bool(args.execute), "baseline": "vaxpress", "runs": receipts}
    write_json(args.output / "campaign.json", result)
    return result
