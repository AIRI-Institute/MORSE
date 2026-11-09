"""Ready-scorer replay of the two retrospective paper panels.

All sequence projection, JSONL communication and output association use Core's
ScorerInputRouter/ScorerBinding/CompositeScorer. Predictors are never fitted here.
"""
from __future__ import annotations

from collections import defaultdict
from copy import copy, deepcopy
import hashlib
import json
import math
from pathlib import Path
import shlex
import shutil
import statistics
import subprocess
import time
import uuid

from .io import DATA, HERE, REPO, read_json, rows, sha256, write_csv, write_json
from .runtime import (_device, _overlay_asset_mounts, _overlay_source_files,
                      _resolve_scorer, _source_file_mounts, _uses_baked_assets, asset_root)

IMAGES = {
    "viennarna": "sha256:04fc0b293a581cc2a6905a1e21ffe843bdcf1d7cc7929d82842db8f69c93f0c6",
    "orthrus-te-held-out": "sha256:a98eb79bb3fa10f4fbe8faab73e53959e6d4e2c3f2fa84b854498bcb050880b8",
    "ribonn": "sha256:bc3d16cd3fc58680c4a5335dcb9310723441a7ba398a66e851d63deb5a6b2f6b",
    "ribodecode": "sha256:70083cd8ae45f6f622bf9d70c0ac4c131cd6c0635d63b04d2ff83af23425cc91",
    "manufacturability": "sha256:ef2ee6c012f2df193ec1d0326f9bd0fe00a6c18f36a22a584190aad4fd4f522a",
    "half-life": "sha256:f38800a12100ff848ece332d28e8937fd4dcb7daf3bf825d91c64175250e1df4",
    "immunogenicity": "sha256:1b664479eb83967cd090e2aace92196f90ff542cc4682b75d4c5222dc1bc4f06",
    "irna": "sha256:43609c166502635b5ad2c9d6d49cafabe785721559d9df5441672314d6707696",
}
ASSETS = {"orthrus-te-full-fit": ["orthrus-te-full-fit"],
          "orthrus-te-held-out": ["orthrus-te-held-out"], "half-life": ["orthrus-half-life"],
          "immunogenicity": ["rinalmo-imrnads"], "ribonn": ["ribonn-human-ensemble"],
          "ribodecode": ["ribodecode-translationmodel"]}
BODY = ["five_prime_utr", "start_codon", "cds", "stop_codon", "three_prime_utr"]
CDS = ["start_codon", "cds"]
CUDA_SCORERS = {"orthrus-te-held-out", "ribonn", "half-life", "immunogenicity"}
RD_BOOTSTRAP = "import torch,runpy; torch.set_num_threads(1); torch.set_num_interop_threads(1); runpy.run_path('/app/server.py',run_name='__main__')"


def panel(name):
    """Load exact constructs and independent expected numeric observations."""
    from morse.construct import SequenceVariant, VaccineConstruct
    if name == "mauger":
        source = DATA / "SupplementaryFigure3_Mauger_U/constructs_U.json"
        records = read_json(source)
        candidates = []
        for r in records:
            cds = r["cds_with_stop"].upper().replace("U", "T")
            c = VaccineConstruct(r["id"], SequenceVariant("mauger_5utr", r["utr5"]), cds[:3], cds[3:-3], cds[-3:], SequenceVariant("mauger_3utr", r["utr3"]), 0)
            if c.sequence != r["transcript_body"]:
                raise ValueError(f"Mauger construct boundary mismatch: {r['id']}")
            candidates.append((r["id"], c))
        if len(candidates) != 39 or len(set(i for i, _ in candidates)) != 39:
            raise ValueError("Expected exactly 39 distinct Mauger IDs")
        expected = {r["id"]: r for r in rows(DATA / "SupplementaryFigure3_Mauger_U/measurements_and_scores.csv")}
        return candidates, expected, None, {"input_sha256": sha256(source), "n": 39}
    source = DATA / "rescore_fig3_panel.json.gz"
    frozen = read_json(source)
    candidates = [(r["id"], VaccineConstruct.from_dict(r["construct"])) for r in frozen["constructs"]]
    for identifier, c in candidates:
        body = c.sequence[:-c.poly_a_length]
        if identifier != "f3p_" + hashlib.sha256(body.replace("T", "U").encode()).hexdigest()[:20]:
            raise ValueError(f"Fig3 construct identifier mismatch: {identifier}")
    expected = {r["construct_id"]: r for r in rows(DATA / "Figure3_logged_scores_0100/analysis/extra_scores.csv")}
    membership = frozen["population_membership"]
    counts = defaultdict(int)
    for row in membership:
        counts[row["generation"]] += 1
    if len(candidates) != 5960 or set(expected) != {i for i, _ in candidates} or set(counts) != set(range(1, 357)) or set(counts.values()) != {64}:
        raise ValueError("Fig3 requires 5960 unique constructs and 64 members in all 356 generations")
    return candidates, expected, membership, {"input_sha256": sha256(source), "source_raw_sha256": frozen["source_raw_sha256"], "n": 5960, "population_memberships": len(membership)}


def recipes(assay):
    """Frozen assay settings; request grouping is part of model provenance."""
    import yaml
    base = yaml.safe_load((HERE / "configs/f2_filtered_s7.yaml").read_text())
    source_presets = read_json(HERE / "provenance/scorer-environments.json")
    profile = source_presets["configs"]["f2_filtered_s7"]["sequence_feasibility"]
    public_scorer = next(s for s in base["scorers"] if s["name"] == "sequence_feasibility")
    feasibility = {k: deepcopy(v) for k, v in public_scorer.items() if k != "command"}
    # Auxiliary replay keeps its historical projection and local metric names.
    feasibility["scope"] = {"type": "parts", "spans": [{"part": p} for p in BODY]}
    feasibility.update(deepcopy(source_presets["profiles"][profile]["recipe"]))
    specs = []
    def add(identifier, logical, parts, outputs, env=None, batch=4):
        reference_device = "cuda" if logical in CUDA_SCORERS else "cpu"
        environment = dict(env or {})
        # Mauger's auxiliary scorers and its copied CPU RiboDecode run used one
        # thread; TE/RiboNN and the later Fig3 panel used two torch/OMP threads.
        threads = "1" if assay == "mauger" and logical not in {"orthrus-te-held-out", "ribonn"} else "2"
        environment.update(OMP_NUM_THREADS=threads, MKL_NUM_THREADS=threads,
                           OPENBLAS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
        if logical == "ribodecode":
            environment.update(BLIS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
        scorer = {"name": identifier, "mode": "persistent-container", "image": IMAGES.get(logical), "timeout_seconds": 1800,
                  "scope": {"type": "parts", "spans": [{"part": p} for p in parts]},
                  "outputs": {raw: {"type": "float", "role": "log"} for raw in outputs.values()}, "env": environment,
                  "paper_runtime": {"scorer": logical, "historical_image": IMAGES.get(logical), "reference_device": reference_device}}
        specs.append({"id": identifier, "logical": logical, "scorer": scorer, "outputs": outputs, "request_batch": batch, "reference_device": reference_device})
    cai_profile = HERE / "inputs/lineardesign_cai_profile.json" if assay == "mauger" else REPO / "scorers/cai/profiles/homo_sapiens_kazusa_93487_cds.json"
    add("cai", "cai", CDS, {"cai_ld" if assay == "mauger" else "cai": "cai"}, {"CAI_PROFILE": str(cai_profile)}, 16)
    if assay == "mauger":
        add("csi", "csi", CDS, {"csi": "csi_percentile"}, {"CSI_PROFILE": "profiles/homo_sapiens_consensus_wu2019_v1.json"}, 16)
        add("viennarna-body", "viennarna", BODY, {"mfe_body": "mfe_per_nt", "gc": "gc_fraction"}, {"VIENNARNA_ENABLE_PF": "0", "VIENNARNA_MAX_BP_SPAN": "200"}, 16)
        feasibility["name"] = "manufacturability"
        feasibility["paper_runtime"]["reference_device"] = "cpu"
        feasibility.setdefault("env", {}).update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
        specs.append({"id": "manufacturability", "logical": "manufacturability", "scorer": feasibility,
                      "outputs": {"feasibility": "feasibility_violation_count"}, "request_batch": 16, "reference_device": "cpu"})
    add("viennarna-cds", "viennarna", CDS, {"mfe_cds" if assay == "mauger" else "cds_mfe": "mfe"}, {"VIENNARNA_ENABLE_PF": "0", "VIENNARNA_MAX_BP_SPAN": "0"}, 16 if assay == "mauger" else 4)
    add("orthrus-te-held-out", "orthrus-te-held-out", BODY, {"te" if assay == "mauger" else "orthrus_te_common_holdout": "translation_efficiency_score"}, {"ORTHRUS_BATCH_SIZE": "4", "ORTHRUS_MAX_TOKENS": "32000"}, 4)
    add("ribonn", "ribonn", BODY, {"ribonn" if assay == "mauger" else "ribonn_mean_te": "ribonn_mean_te"}, {"RIBONN_BATCH_SIZE": "16" if assay == "mauger" else "4"}, 16 if assay == "mauger" else 4)
    add("half-life", "half-life", BODY, {"half_life" if assay == "mauger" else "half_life_score": "half_life_score"}, {"ORTHRUS_BATCH_SIZE": "2", "ORTHRUS_MAX_TOKENS": "4096"}, 16 if assay == "mauger" else 4)
    if assay == "mauger":
        add("ribodecode", "ribodecode", CDS, {"ribodecode": "ribodecode_translation_log"}, {"RIBODECODE_ENV": "HEK293T", "RIBODECODE_BATCH_SIZE": "32"}, 32)
    add("irna", "irna", BODY, {"mirna" if assay == "mauger" else "irna_max_profile_risk": "irna_max_profile_risk"},
        {"IRNA_SPECIES": "human", "IRNA_WORKERS": "2", "IRNA_MAX_LOG_KD": "-1.0", "IRNA_ORF_MAX_LOG_KD": "-1.5"}, 16 if assay == "mauger" else 64)
    add("immunogenicity", "immunogenicity", BODY, {"immunogenicity" if assay == "mauger" else "rinalmo_window_logit_p95": "rinalmo_window_logit_p95"},
        {"RINALMO_WINDOW_SIZE": "27", "RINALMO_WINDOW_STRIDE": "3", "RINALMO_WINDOW_OFFSETS": "", "RINALMO_BATCH_SIZE": "256"}, 4)
    if assay == "mauger":
        for spec in specs:
            spec["historical_panel_context"] = (
                "Separate 39-row Mauger job; no preceding Leppek rows"
                if spec["logical"] in {"ribonn", "ribodecode"} else
                "Mauger suffix of Leppek24+Mauger39. The 24 excluded Leppek rows do not change TE/RiNALMo groups4. "
                "Half-life sorts each outer request and uses batches2: even prefix24 and equal 1815-nt Mauger lengths preserve every model pair and tensor shape. "
                "Remaining scorers operate independently per sequence.")
    return specs


def effective_device(spec, requested):
    """Select a GPU index without moving historical CPU inference onto CUDA."""
    _device(requested)  # Validate the user selection even for a CPU-only scorer.
    reference = spec["scorer"]["paper_runtime"]["reference_device"]
    if reference not in {"cpu", "cuda"}:
        raise ValueError(f"Unknown reference device for {spec['id']}: {reference}")
    return "cpu" if reference == "cpu" else requested


def resolve(spec, args, directory):
    """Reuse optimization runtime resolution, then add the three auxiliary adapters."""
    opts = copy(args)
    opts.device = effective_device(spec, args.device)
    opts.source_overlay = False  # Correct source entrypoint is installed below for every logical scorer.
    raw = _resolve_scorer(spec["scorer"], opts, uuid.uuid4().hex[:8])
    logical = spec["logical"]
    baked = _uses_baked_assets(args, logical)
    assets = asset_root(args.asset_root)
    gpu = _device(opts.device)
    device = opts.device if args.runtime == "process" else "cpu" if gpu is None else "cuda:0"
    env = raw.setdefault("env", {})
    for key,value in {"OMP_NUM_THREADS":"2","OPENBLAS_NUM_THREADS":"1","MKL_NUM_THREADS":"2","NUMEXPR_NUM_THREADS":"1","HF_HUB_OFFLINE":"1","TRANSFORMERS_OFFLINE":"1","HF_HOME":"/tmp/hf"}.items():env.setdefault(key,value)
    if logical == "half-life":
        env.update(ORTHRUS_CHECKPOINT=str(assets / "orthrus-half-life/model_final.pt"), ORTHRUS_MODEL_SOURCE=str(assets / "orthrus-backbone"), ORTHRUS_DEVICE=device,
                   ORTHRUS_PORTABLE="1" if args.runtime == "process" or gpu is None else "0")
    elif logical == "immunogenicity":
        env.update(RINALMO_CHECKPOINT=str(assets / "rinalmo-imrnads/best_model.pt"), RINALMO_DEVICE=device,
                   MORSE_EXPECTED_CHECKPOINT_SHA256="923d1e9f843d12fe29a9515a95b6a32eae286017790dfeef1c82215446958f0a")
    elif logical == "irna":
        env["IRNA_PANEL_ROOT"] = str(REPO / "scorers/irna/panels/2026-08-21") if args.runtime == "process" else "/app/panels/2026-08-21"
        if args.runtime == "process":
            raw["command"] = ["Rscript", str(REPO / "scorers/irna/server.R")]
    if logical == "ribodecode" and raw["mode"] == "persistent-process":
        bootstrap = RD_BOOTSTRAP.replace("'/app/server.py'", repr(str(REPO / "scorers/ribodecode/server.py")))
        raw["command"] = [raw["command"][0], "-u", "-c", bootstrap]
    if raw["mode"] == "persistent-container":
        docker = raw.setdefault("docker_args", [])
        if logical in {"half-life", "immunogenicity"}:
            if gpu is not None:
                docker.extend(["--gpus", f"device={gpu}"])
            if logical=="half-life":
                if not baked:
                    docker.extend(["--volume",f"{assets/'orthrus-half-life'}:/models:ro","--volume",f"{assets/'orthrus-backbone'}:/backbone:ro"])
                env.update(ORTHRUS_CHECKPOINT="/models/model_final.pt",ORTHRUS_MODEL_SOURCE="/backbone")
            else:
                if not baked:
                    docker.extend(["--volume",f"{assets/'rinalmo-imrnads'}:/models:ro"])
                env["RINALMO_CHECKPOINT"]="/models/best_model.pt"
        # Auxiliary utility images originally used mounted adapters; provide the
        # checked-out inference-only source and an explicit executable launcher.
        overlay = args.source_overlay or (not baked and logical in {"half-life", "immunogenicity", "ribodecode"})
        # Keep the frozen RiboDecode CPU thread bootstrap without mounting its
        # host package or replacing the baked adapter when no overlay is asked.
        if overlay or logical == "ribodecode":
            folder = "translation-efficiency" if logical.startswith("orthrus-te") else logical
            source = REPO / "scorers" / folder
            launcher = directory / f"entrypoint-{spec['id']}"
            command = "Rscript /app/server.R" if logical == "irna" else "python3 /app/scorer.py" if logical == "viennarna" else "python3 -m manufacturability_scorer.feasibility_server" if logical == "manufacturability" else "python3 /app/server_heldout.py" if logical == "orthrus-te-held-out" else "python3 /app/server.py"
            if logical == "ribodecode":
                command = "python3 -u -c " + shlex.quote(RD_BOOTSTRAP)
            launcher.write_text("#!/bin/sh\ncd /app\nexec " + command + "\n")
            launcher.chmod(0o755)
            if overlay:
                if baked:
                    docker.extend(_source_file_mounts(logical))
                    raw["provenance_paths"] = [str(path) for path in _overlay_source_files(logical)]
                else:
                    docker.extend(["--volume", f"{source}:/app:ro"])
                    raw["provenance_paths"] = [str(p) for p in source.rglob("*") if p.is_file() and p.suffix in {".py",".R",".json",".txt"} and not any(part in {"__pycache__",".venv","tests"} for part in p.relative_to(source).parts)]
            docker.extend(["--volume", f"{launcher}:/morse-rescore-entrypoint:ro", "--entrypoint", "/morse-rescore-entrypoint"])
        if args.source_overlay and not baked:
            mounts, current_env = _overlay_asset_mounts(logical, assets)
            docker.extend(mounts)
            env.update(current_env)
    return raw


def preflight(specs, args):
    from assets.manage import verify_asset
    manifest = {r["id"]: r for r in read_json(REPO / "assets/registry.json")["assets"]}
    cache = asset_root(args.asset_root)
    if _device(args.device) is None and any(s["logical"]=="orthrus-te-held-out" for s in specs):
        raise RuntimeError("The frozen held-out TE adapter requires a CUDA/Mamba runtime; use --device cuda:N or select supported light scorers for a CPU diagnostic")
    verified = []
    external_specs = [spec for spec in specs if not _uses_baked_assets(args, spec["logical"])]
    current_containers = {spec["logical"] for spec in external_specs if args.source_overlay
                          and spec["resolved"]["mode"] == "persistent-container"}
    for spec in specs:
        if _uses_baked_assets(args, spec["logical"]):
            for identifier in ASSETS.get(spec["logical"], [spec["logical"]]):
                verified.append({"id": identifier, "scorer": spec["logical"], "status": "bundled",
                                 "source": "immutable-published-image", "image": spec["resolved"]["image"]})
    for identifier in sorted({a for spec in external_specs for a in ASSETS.get(spec["logical"], [])}):
        runtime = identifier == "ribodecode-translationmodel" and "ribodecode" in current_containers
        try:
            verified.append(verify_asset(manifest[identifier], cache / identifier,
                                         dependency_root=cache, require_runtime=runtime))
        except (OSError, ValueError) as error:
            if not runtime:
                raise
            raise RuntimeError("Current RiboDecode source overlay requires the verified external package. "
                               "Prepare it with 'morse assets install --asset ribodecode-translationmodel' "
                               f"using the same asset root. {error}") from error
    if "ribonn" in current_containers:
        verified.append(verify_asset(manifest["ribonn-inference-source"],
                                     REPO / "scorers/ribonn/third_party/ribonn"))
    # Profiles and the human scanMiR panel are versioned source assets. Verify
    # those actual source files, because the scorers read the bundled copies.
    for identifier in {"cai-profiles", "csi-profiles", "scanmir-panels", "sequence-feasibility-profiles"}:
        if identifier not in manifest:
            continue
        for item in manifest[identifier].get("files", []):
            source = REPO / item["repository_path"]
            if sha256(source) != item["sha256"]:
                raise ValueError(f"Bundled scorer profile changed: {source}")
    for spec in specs:
        raw = spec["resolved"]
        if raw["mode"] == "persistent-container":
            try:
                result = subprocess.run(["docker", "image", "inspect", raw["image"], "--format", "{{.Id}}"], capture_output=True, text=True)
            except FileNotFoundError as exc:
                raise RuntimeError("Docker is required for container replay; install it or use --runtime process with all real dependencies") from exc
            if result.returncode:
                raise RuntimeError(f"Missing {spec['logical']} image {raw['image']}; build/install the real scorer image and supply --image-map")
            spec["actual_image_id"] = result.stdout.strip()
        elif shutil.which(raw["command"][0]) is None:
            raise RuntimeError(f"Missing scorer interpreter {raw['command'][0]} for {spec['id']}")
        if spec["logical"] in {"viennarna","irna"}:
            if spec["logical"]=="viennarna":
                program=["python3","-c","import RNA; print(RNA.__version__)"]
            else:
                program=["Rscript","-e",'cat(as.character(packageVersion("scanMiR")),as.character(packageVersion("scanMiRData")),sep=" ")']
            if raw["mode"]=="persistent-container":
                command=["docker","run","--rm","--pull=never","--network=none","--entrypoint",program[0],raw["image"],*program[1:]]
            else:command=[raw["command"][0],*program[1:]]
            process=subprocess.run(command,capture_output=True,text=True,timeout=60)
            if process.returncode:raise RuntimeError(f"Missing real {spec['logical']} runtime dependency: {process.stderr.strip()}")
            observed=process.stdout.strip()
            expected="2.7.2" if spec["logical"]=="viennarna" else "1.16.0 1.16.0"
            if observed!=expected:raise ValueError(f"Wrong {spec['logical']} runtime version: expected {expected}, observed {observed}")
            spec["runtime_version"]=observed
    return verified


def aggregate_fig3(values, membership, output):
    """Rebuild the exact 64-member populations, preserving repeated membership."""
    if len(values) != 5960:
        return {"status": "not-run", "reason": "Partial diagnostics cannot represent full historical populations"}
    protocol = read_json(DATA / "Figure3_logged_scores_0100/analysis/normalization.json")
    available = [m for m in protocol["metrics"] if all(m["key"] in v for v in values.values())]
    normalized = {}
    for metric in available:
        key = metric["key"]
        scores = [v[key] for v in values.values()]
        low, high = min(scores), max(scores)
        if high == low:
            raise ValueError(f"Constant retrospective output requires an explicit normalization policy: {key}")
        normalized[key] = {"minimum": low, "maximum": high, "direction": metric["direction"]}
    by_generation = defaultdict(list)
    for row in membership:
        by_generation[row["generation"]].append(row["construct_id"])
    source = {int(r["generation"]): r for r in rows(DATA / "Figure3_logged_scores_0100/analysis/extra_population_trajectories.csv")}
    trajectories, comparisons = [], []
    for generation in range(1, 357):
        members = by_generation[generation]
        if len(members) != 64:
            raise ValueError(f"Incomplete original population at generation {generation}")
        row = {"generation": generation, "population_n": 64}
        for metric in available:
            key = metric["key"]
            native = [values[identifier][key] for identifier in members]
            scale = normalized[key]
            q = [(v-scale["minimum"])/(scale["maximum"]-scale["minimum"]) if metric["direction"] == "max" else (scale["maximum"]-v)/(scale["maximum"]-scale["minimum"]) for v in native]
            row[key + "_raw_mean"] = statistics.fmean(native)
            row[key + "_relative_mean"] = statistics.fmean(q)
            row[key + "_sd_population"] = statistics.stdev(native)
            comparisons.append({"generation": generation, "metric": key, "native_mean": row[key+"_raw_mean"], "expected_native_mean": float(source[generation][key]), "difference": row[key+"_raw_mean"]-float(source[generation][key]),
                                "relative_mean": row[key+"_relative_mean"], "expected_relative_mean": float(source[generation][key+"_normalized"])})
        trajectories.append(row)
    write_csv(output / "Figure3_modular_design.d_series.csv", trajectories)
    write_csv(output / "population_comparison.csv", comparisons)
    write_json(output / "recomputed_display_normalization.json", normalized)
    return {"status": "completed", "generations": 356, "members_per_generation": 64, "metrics": [m["key"] for m in available], "all_seven_metrics": len(available) == 7,
            "selection_unchanged": True, "normalization": "direction-aware retrospective min/max across all 5960 freshly rescored unique constructs"}


def _mauger_rank_outputs(values, expected, output):
    if len(values) != 39 or any(len(v) != 12 for v in values.values()):
        return {"status": "not-run", "reason": "All 39 constructs and 12 output measurements are required for the fixed composites"}
    import numpy as np
    from scipy.stats import rankdata
    protocol = read_json(DATA / "SupplementaryFigure3_Mauger_U/analysis_protocol.json")
    keys = protocol["individual_scorers_retained"]
    identifiers = list(values)
    ranks = {key: rankdata([values[i][key] for i in identifiers]) for key in keys}
    products = {}
    for label, members in protocol["combinations"].items():
        products[label] = [math.prod(int(round(2*ranks[key][i])) for key in members) for i in range(39)]
        ranks[label] = rankdata(np.asarray(products[label], dtype=np.int64))
    y = np.asarray([float(expected[i]["U_bar_height_pt"]) for i in identifiers])
    yr = rankdata(y)
    def unit(a):
        z = a-a.mean()
        norm=np.linalg.norm(z)
        if not np.isfinite(z).all() or not math.isfinite(float(norm)) or norm<=0:raise ValueError("Undefined correlation: constant or nonfinite rank input")
        return z/norm
    rng = np.random.default_rng(20261129)
    permuted = np.empty((100000, 39));permuted[0] = unit(yr)
    for start in range(1, 100000, 4096):
        count = min(4096, 100000-start)
        permuted[start:start+count] = rng.permuted(np.broadcast_to(permuted[0], (count, 39)), axis=1)
    if not np.isfinite(permuted).all() or np.max(np.abs(permuted))>1:raise ValueError("Invalid bounded permutation vectors")
    stats = []
    for key in keys + ["all translation", "All Relevant"]:
        oriented=unit(ranks[key]);rho = float(np.sum(unit(yr)*oriented));null = np.einsum("ij,j->i",permuted,oriented,optimize=False)
        if not np.isfinite(null).all() or np.max(np.abs(null))>1+1e-12:raise ValueError(f"Invalid permutation correlation for {key}")
        p = float(np.mean(np.round(np.abs(null), 12) >= round(abs(rho), 12)))
        stats.append({"scorer": key, "n": 39, "spearman_rho": rho, "p_value": p, "random_permutations": 99999, "permutation_seed": 20261129})
    order = np.argsort(y, kind="stable")
    matrix = np.asarray([100*(yr[order]-1)/38] + [100*(ranks[key][order]-1)/38 for key in keys+["all translation", "All Relevant"]])
    write_csv(output / "mauger_recomputed_statistics.csv", stats)
    write_json(output / "mauger_rank_heatmap.json", {"ids": [identifiers[i] for i in order], "rows": ["experiment", *keys, "all translation", "All Relevant"], "percentiles": matrix.tolist(), "exact_twice_rank_products": products})
    return {"status": "completed", "constructs": 39, "components": 12, "composites": 2, "permutation_tests": 14}


def run(args):
    from morse.cli import _make_one_scorer
    from morse.config import _parse_scorer_config
    from morse.scorers import CompositeScorer, ConstructCandidate, ScorerBinding
    from .verification import verify
    verify()
    if args.limit is not None and args.limit <= 0:
        raise ValueError("--limit must be positive")
    if args.batch_size is not None and args.batch_size <= 0:
        raise ValueError("--batch-size must be positive")
    directory = args.output.resolve()
    if directory.exists() and any(directory.iterdir()):
        raise ValueError(f"Refusing to overwrite rescore evidence: {directory}")
    all_candidates, expected, membership, evidence = panel(args.assay)
    target = all_candidates[:args.limit] if args.limit else all_candidates
    target_ids = {identifier for identifier, _ in target}
    specs = recipes(args.assay)
    if args.scorer:
        unknown = set(args.scorer)-{s["id"] for s in specs}
        if unknown:
            raise ValueError(f"Unknown scorer IDs {sorted(unknown)}; available: {[s['id'] for s in specs]}")
        specs = [s for s in specs if s["id"] in args.scorer]
    directory.mkdir(parents=True, exist_ok=True)
    for spec in specs:
        spec["resolved"] = resolve(spec, args, directory)
    diagnostic = args.limit is not None or args.batch_size is not None or bool(args.scorer) or (args.device == "cpu" and any(s["reference_device"] == "cuda" for s in specs))
    plan = {"assay": args.assay, "mode": "diagnostic" if diagnostic else "paper-replay", "target_constructs": len(target), "canonical_constructs": len(all_candidates), "input_evidence": evidence,
            "runtime": args.runtime, "device": args.device, "scorers": specs, "batch_protocol": "Preserve canonical row order and frozen neural model batch context, including canonical siblings beyond --limit. Mauger excludes the preceding Leppek24 only where model batches remain identical; see each scorer's historical_panel_context. An explicit --batch-size changes request grouping and is diagnostic.",
            "predictor_training": False, "optimization": False}
    write_json(directory / "plan.json", plan)
    if not args.execute:
        return plan
    values = {identifier: {} for identifier, _ in target}
    native_values = {identifier: {} for identifier, _ in target}
    comparisons, timings, checks = [], [], []
    try:
        assets = preflight(specs, args)
        write_json(directory/"asset_verification.json",assets)
        write_json(directory/"plan.json",plan)
        for spec in specs:
            raw = spec["resolved"]
            config = _parse_scorer_config(raw, spec["id"])
            component = _make_one_scorer(config)
            scorer = CompositeScorer([ScorerBinding(component, config.scope)])
            request_size = args.batch_size or spec["request_batch"]
            phase_start = time.monotonic()
            scored = 0
            try:
                with (directory / (spec["id"] + ".jsonl")).open("x") as stream:
                    for offset in range(0, len(target), request_size):
                        batch = all_candidates[offset:offset+request_size]
                        records = [ConstructCandidate(identifier, construct) for identifier, construct in batch]
                        results = scorer.score(records)
                        routed=scorer.bindings[0].routed_sequences(records)
                        input_digests={r.candidate_id:hashlib.sha256(json.dumps({"sequence":q.sequence,"context":q.context()},sort_keys=True).encode()).hexdigest() for r,q in zip(records,routed)}
                        for result in results:
                            if spec["logical"] == "irna" and (result.outputs.get("irna_panel_id") != "human-core_975-2026-08-21" or result.outputs.get("irna_species") != "human"):
                                raise ValueError("scanMiR returned a different species/panel than the paper")
                            identifier = result.candidate_id
                            stream.write(json.dumps({"id": identifier, "context_only": identifier not in target_ids,"routed_input_sha256":input_digests[identifier], "outputs": result.outputs}, allow_nan=False)+"\n")
                            if identifier not in target_ids:
                                continue
                            for key, native_key in spec["outputs"].items():
                                if native_key not in result.outputs:
                                    raise ValueError(f"Missing {spec['id']} required output {native_key} for {identifier}")
                                native = float(result.outputs[native_key])
                                if not math.isfinite(native):
                                    raise ValueError(f"Non-finite {spec['id']} output for {identifier}")
                                sign = -1 if args.assay == "mauger" and key in {"mfe_body", "mfe_cds", "mirna", "immunogenicity", "feasibility"} else 1
                                values[identifier][key] = sign*native
                                native_values[identifier][key] = native
                                reference = sign*float(expected[identifier][key])
                                neural = spec["logical"] in {"orthrus-te-held-out", "ribonn", "half-life", "ribodecode", "immunogenicity"}
                                rtol, atol = (1e-5, 1e-5) if neural else (1e-10, 1e-12)
                                comparisons.append({"id": identifier, "scorer": spec["id"], "metric": key, "native_output": native_key, "actual_native": native, "expected_native": reference, "difference": native-reference,
                                    "orientation": sign, "rtol": rtol, "atol": atol, "matches": math.isclose(native, reference, rel_tol=rtol, abs_tol=atol)})
                        scored += len(batch)
                        timings.append({"scorer": spec["id"], "offset": offset, "request_size": len(batch), "ids": [identifier for identifier, _ in batch], "seconds": scorer.last_timings[0].elapsed_seconds})
                        stream.flush()
                        write_json(directory/"progress.json",{"scorer":spec["id"],"targets_completed":min(offset+len(batch),len(target)),"target_total":len(target),"scored_including_context":scored,"elapsed_seconds":time.monotonic()-phase_start})
            finally:
                scorer.close()
            source_files=[Path(p) for p in raw.get("provenance_paths",[])]
            if raw.get("workdir"):
                source_directory=Path(raw["workdir"])
                source_files += [p for p in source_directory.rglob("*") if p.is_file() and p.suffix in {".py",".R",".json",".txt"} and not any(part in {"__pycache__",".venv","tests"} for part in p.relative_to(source_directory).parts)]
            checks.append({"scorer": spec["id"], "target_constructs": len(target), "scored_with_context": scored, "seconds": time.monotonic()-phase_start, "request_batch": request_size,
                           "actual_image_id":spec.get("actual_image_id"),"runtime_version":spec.get("runtime_version"),"source_sha256":{str(p.relative_to(REPO)):sha256(p) for p in source_files}})
        write_csv(directory / "native_scores.csv", [dict(id=identifier, **scores) for identifier, scores in native_values.items()])
        if args.assay=="mauger":write_csv(directory/"oriented_scores.csv",[dict(id=identifier,**scores) for identifier,scores in values.items()])
        write_csv(directory / "comparison.csv", comparisons)
        write_json(directory / "request_timings.json", timings)
        aggregate = aggregate_fig3(values, membership, directory) if args.assay == "fig3-logged" else _mauger_rank_outputs(values, expected, directory)
        result = {"status": "completed" if all(c["matches"] for c in comparisons) else "completed_with_differences", "mode": plan["mode"], "assay": args.assay, "target_constructs": len(target), "metric_comparisons": len(comparisons),
                  "matched_comparisons": sum(c["matches"] for c in comparisons), "maximum_absolute_difference": max(abs(c["difference"]) for c in comparisons), "all_native_scores_match": all(c["matches"] for c in comparisons),
                  "scorers": checks, "assets_verified": [a["id"] for a in assets], "aggregation": aggregate, "predictor_training": False,
                  "scope": "Ready-scorer inference on frozen constructs; no optimization, reselection or model fitting"}
        write_json(directory / "rescore_receipt.json", result)
        return result
    except Exception as exc:
        write_json(directory / "failure.json", {"status": "failed", "assay": args.assay, "error": str(exc), "completed_scorers": checks})
        raise
