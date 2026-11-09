"""Launch separately obtained upstream design tools without bundling their code."""
from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess

from .io import HERE, REPO, read_json, registry, sha256, write_json
from .runtime import _device

LD_COMMIT = "f0126ca89a8b853088b4bccfd2cc8c378d3678be"


def run_upstream(args):
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.tool == "lineardesign":
        if args.source is None:
            raise ValueError("--source must name your separately obtained LinearDesign checkout; see paper/UPSTREAM.md")
        source = args.source.resolve()
        revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
        if revision != LD_COMMIT:
            raise ValueError(f"LinearDesign source revision must be {LD_COMMIT}; observed {revision}")
        command = [str(source / "bin/LinearDesign_2D"), "3", "0", str(source / "codon_usage_freq_table_human.csv")]
        protein = "".join(line.strip() for line in (HERE / "inputs/mlana_protein.fasta").read_text().splitlines() if not line.startswith(">"))
        result = {"tool": "LinearDesign", "revision": revision, "lambda": 3, "command": command, "protein": protein, "model_training": False}
        if args.execute:
            if not Path(command[0]).is_file():
                subprocess.run(["make"], cwd=source, check=True)
            process = subprocess.run(command, input=protein+"\n", capture_output=True, text=True, check=True, cwd=source)
            (args.output / "stdout.txt").write_text(process.stdout)
            (args.output / "stderr.txt").write_text(process.stderr)
            sequences = re.findall(r"[ACGU]{100,}", process.stdout)
            if not sequences:
                raise ValueError("LinearDesign produced no identifiable RNA sequence; inspect stdout.txt")
            sequence = max(sequences, key=len).replace("U", "T")
            from morse.codon import translate_dna
            if translate_dna(sequence).rstrip("*") != protein:
                raise ValueError("LinearDesign output does not preserve MLANA protein")
            (args.output / "design.fasta").write_text(">MLANA_LinearDesign_lambda3\n"+sequence+"\n")
            result.update(status="completed", sequence_sha256=sha256(args.output / "design.fasta"))
    else:
        config = read_json(HERE / "inputs/native_rd.json")
        image = args.image or config["image"]
        gpu = _device(args.device)
        if gpu is None:
            raise ValueError("Frozen native RiboDecode generation requires CUDA; its authors' script calls CUDA directly")
        command = ["docker", "run", "--rm", "--network=none", "--gpus", f"device={gpu}",
            "--env", f"PYTHONHASHSEED={args.seed}", "--env", "CUBLAS_WORKSPACE_CONFIG=:4096:8",
            "--env", "CUDA_VISIBLE_DEVICES=0", "--volume", f"{HERE}:/paper:ro", "--volume", f"{args.output}:/work:rw",
            "--workdir", "/work", "--entrypoint", "python3.8", image, "/paper/adapters/seeded_rd.py",
            "--seed", str(args.seed), "--config", "/paper/inputs/native_rd.json"]
        if args.pilot:
            command.append("--pilot")
        result = {"tool": "RiboDecode native generation", "seed": args.seed, "command": command,
                  "epochs": 1 if args.pilot else 10, "final_samples": 400 if args.pilot else 40000,
                  "mode": "diagnostic" if args.pilot else "paper-protocol", "predictor_training": False}
        if args.execute:
            subprocess.run(command, check=True)
            from morse.codon import translate_dna
            candidates = {}
            path = args.output / f"results_natural/HEK293T_optim_mfe_dist-optim/{0 if args.pilot else 9}/samples/optim_results.txt"
            protein = translate_dna(config["initial_sequence"][:-3])
            for line in path.read_text().splitlines():
                sequence, score, _ = line.split("\t")
                sequence = sequence.upper().replace("U", "T")
                if len(sequence) != 357 or not sequence.startswith("ATG") or sequence[-3:] not in {"TAA", "TAG", "TGA"}:
                    continue
                if translate_dna(sequence[:-3]) != protein:
                    continue
                coding = sequence[:-3]
                if coding not in candidates or float(score) > candidates[coding]["native_score"]:
                    candidates[coding] = {"coding_sequence": coding, "native_score": float(score), "original_stop": sequence[-3:]}
            if not candidates:
                raise ValueError("No protein-preserving native RiboDecode candidate")
            selected = min(candidates.values(), key=lambda row: (-row["native_score"], row["coding_sequence"]))
            selected.update(seed=args.seed, source_sha256=sha256(path), selection="highest native score before feasibility screening; lexicographic CDS breaks ties", standardized_stop="TAA")
            write_json(args.output / "selected_native.json", selected)
            (args.output / "design.fasta").write_text(">MLANA_RiboDecode_native\n"+selected["coding_sequence"]+"TAA\n")
            result.update(status="completed", selection=selected)
    write_json(args.output / "upstream_receipt.json", result)
    return result
