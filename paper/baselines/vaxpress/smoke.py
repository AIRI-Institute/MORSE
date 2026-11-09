"""Prepare or execute the original native VaxPress driver for exactly 64 events."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

REPO = Path(__file__).resolve().parents[3]


def prepare(output: Path, image: str, repo: Path = REPO):
    output, repo = output.expanduser().resolve(), repo.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to replace baseline smoke evidence: {output}")
    config = json.loads((repo / "paper/inputs/vaxpress_protocol.json").read_text())
    for key, relative in [("csi_model", "scorers/csi/model.py"),
                          ("csi_profile", "scorers/csi/profiles/homo_sapiens_consensus_wu2019_v1.json"),
                          ("viennarna_model", "scorers/viennarna/scorer.py")]:
        actual = hashlib.sha256((repo / relative).read_bytes()).hexdigest()
        if actual != config[key + "_sha256"]:
            raise ValueError(f"Frozen shared scorer source/profile changed: {relative}")
        config[key + "_path"] = "/morse/" + relative
    command = ["docker", "run", "--rm", "--network=none", "--cpus=4", "--memory=2g", "--volume", f"{repo}:/morse:ro",
        "--volume", f"{output}:/output:rw", "--entrypoint", "python3", image,
        "/morse/paper/adapters/run_vaxpress_shared.py", "--config", "/output/protocol.json", "--output", "/output/native64",
        "--seed", "7", "--alpha", "0", "--budget", "64", "--processes", "4"]
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".vaxpress-smoke-", dir=output.parent) as temporary:
        stage = Path(temporary) / "evidence"
        stage.mkdir()
        (stage / "protocol.json").write_text(json.dumps(config, indent=2) + "\n")
        (stage / "command.json").write_text(json.dumps(command, indent=2) + "\n")
        stage.rename(output)
    return command


def verify(output: Path):
    native = output / "native64"
    summary = json.loads((native / "summary.json").read_text())
    rows = [json.loads(line) for line in (native / "evaluations.jsonl").read_text().splitlines()]
    batches = [json.loads(line) for line in (native / "budget_batches.jsonl").read_text().splitlines()]
    if len(rows) != 64 or summary.get("events") != 64 or not summary.get("budget_reached"):
        raise ValueError("Native driver did not record the exact 64-event boundary")
    if sum(b["scored_events"] for b in batches) != 64 or batches[-1]["scored_events"] != 4 or batches[-1]["discarded_unscored_proposals"] != 16:
        raise ValueError("The last native batch was not truncated from 20 proposals to 4 scored events")
    if any(row["csi_codons"] != 117 or row["transcript_body_length"] != 470 for row in rows):
        raise ValueError("The frozen shared score context changed")
    return {"status": "pass", "mode": "diagnostic-native-64-events", "events": 64,
            "unique_cds": len({r["cds"] for r in rows}), "last_batch_scored": 4, "last_batch_discarded_unscored": 16,
            "csi_codons": 117, "transcript_body_length": 470, "native_summary": summary}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", default="morse/vaxpress:dfc366-shared")
    parser.add_argument("--output", type=Path, required=True, help="A new directory; plan and execution outputs are exclusive")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args(argv)
    image = args.image
    if args.execute:
        image = subprocess.run(["docker", "image", "inspect", image, "--format", "{{.Id}}"], capture_output=True, text=True, check=True).stdout.strip()
        if not image.startswith("sha256:"):
            raise ValueError("Cannot resolve the baseline image to an immutable local identity")
    command = prepare(args.output, image)
    print(json.dumps(command, indent=2), flush=True)
    if args.execute:
        output = args.output.expanduser().resolve()
        with (output / "stdout.log").open("w") as stdout, (output / "stderr.log").open("w") as stderr:
            result = subprocess.run(command, stdout=stdout, stderr=stderr)
        if result.returncode:
            (output / "execution.json").write_text(json.dumps({"status": "failed", "exit_code": result.returncode, "image": image}, indent=2) + "\n")
            raise SystemExit(result.returncode)
        receipt = {**verify(output), "image": image}
        (output / "execution.json").write_text(json.dumps(receipt, indent=2) + "\n")
        print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
