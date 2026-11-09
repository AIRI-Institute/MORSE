"""Benchmark distinct FASTA records with the native and persistent JSONL APIs."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import random
import selectors
import statistics
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

from .feasibility import SequenceFeasibilityScorer
from .sequence import normalize_rna, read_fasta


def statistics_ms(seconds):
    ordered = sorted(value * 1000 for value in seconds)
    return {
        "count": len(ordered),
        "mean_ms": statistics.mean(ordered),
        "median_ms": statistics.median(ordered),
        "p95_ms": ordered[math.ceil(0.95 * len(ordered)) - 1],
        "min_ms": ordered[0],
        "max_ms": ordered[-1],
        "total_ms": sum(ordered),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--count", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260918)
    parser.add_argument("--poly-a-length", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if args.count < 1 or args.poly_a_length < 0:
        parser.error("count must be positive and poly-a-length nonnegative")
    unique = {}
    for name, raw in read_fasta(args.fasta):
        sequence = normalize_rna(raw)
        if args.poly_a_length:
            if (
                len(sequence) <= args.poly_a_length
                or sequence[-args.poly_a_length :] != "A" * args.poly_a_length
            ):
                parser.error("Input has an invalid annotated poly(A) tail")
            sequence = sequence[: -args.poly_a_length]
        unique.setdefault(sequence, name)
    if len(unique) < args.count:
        parser.error(f"Need {args.count} distinct sequences; found {len(unique)}")
    records = random.Random(args.seed).sample(list(unique.items()), args.count)
    sequences = [seq for seq, name in records]
    scorer = SequenceFeasibilityScorer()
    # No memoization is used; warm-up excludes interpreter/container startup.
    for sequence in sequences[: min(10, args.count)]:
        scorer.score_sequence(sequence)
    native, expected = [], []
    for sequence in sequences:
        start = time.perf_counter_ns()
        output = scorer.score_sequence(sequence)
        native.append((time.perf_counter_ns() - start) / 1e9)
        expected.append(output)
    command = args.command or [
        sys.executable,
        "-m",
        "manufacturability_scorer.feasibility_server",
    ]
    process = None
    request_id = 0

    def exchange(batch):
        nonlocal request_id
        request_id += 1
        process.stdin.write(
            json.dumps({"request_id": request_id, "sequences": batch}) + "\n"
        )
        process.stdin.flush()
        with selectors.DefaultSelector() as selector:
            selector.register(process.stdout, selectors.EVENT_READ)
            if not selector.select(120):
                raise TimeoutError("Persistent scorer did not respond within 120 s")
        line = process.stdout.readline()
        if not line:
            raise RuntimeError("Persistent scorer closed stdout")
        response = json.loads(line)
        if response.get("request_id") != request_id or "error" in response:
            raise RuntimeError(f"Invalid scorer response: {response}")
        return response["outputs"]

    try:
        start = time.perf_counter_ns()
        process = subprocess.Popen(
            command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1
        )
        first = exchange([sequences[0]])
        startup = (time.perf_counter_ns() - start) / 1e9
        if first != expected[:1]:
            raise RuntimeError("Runtime and native outputs differ at startup")
        for sequence in sequences[: min(10, args.count)]:
            exchange([sequence])
        roundtrips = []
        for sequence, output in zip(sequences, expected):
            start = time.perf_counter_ns()
            actual = exchange([sequence])
            roundtrips.append((time.perf_counter_ns() - start) / 1e9)
            if actual != [output]:
                raise RuntimeError("Runtime and native outputs differ")
        start = time.perf_counter_ns()
        batch_outputs = exchange(sequences)
        batch_seconds = (time.perf_counter_ns() - start) / 1e9
        if batch_outputs != expected:
            raise RuntimeError("Batch and singleton outputs differ")
    finally:
        if process is not None:
            process.stdin.close()
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.terminate()
                process.wait(timeout=10)
            process.stdout.close()
    lengths = [len(s) for s in sequences]
    cpu = platform.processor()
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        cpu = next(
            (
                line.split(":", 1)[1].strip()
                for line in cpuinfo.read_text().splitlines()
                if line.startswith("model name")
            ),
            cpu,
        )
    failures = Counter(
        rule
        for result in expected
        for rule in json.loads(result["feasibility_failed_checks_json"])
    )
    report = {
        "schema_version": 1,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "hostname": platform.node(),
        "platform": platform.platform(),
        "cpu": cpu,
        "python": sys.version,
        "runtime_command": command,
        "input_fasta": str(args.fasta.resolve()),
        "input_sha256": hashlib.sha256(args.fasta.read_bytes()).hexdigest(),
        "sample_seed": args.seed,
        "unique_available": len(unique),
        "count": args.count,
        "excluded_annotated_poly_a_nt": args.poly_a_length,
        "length_nt": {
            "min": min(lengths),
            "mean": statistics.mean(lengths),
            "max": max(lengths),
        },
        "configuration": json.loads(scorer.configuration_json),
        "configuration_sha256": scorer.configuration_sha256,
        "native": statistics_ms(native),
        "persistent_single_sequence": statistics_ms(roundtrips),
        "persistent_batch": {
            "count": args.count,
            "total_ms": batch_seconds * 1000,
            "amortized_ms_per_sequence": batch_seconds * 1000 / args.count,
        },
        "startup_including_first_score_ms": startup * 1000,
        "warmup_sequences": min(10, args.count),
        "outputs_verified_equal": True,
        "passing_sequences": sum(o["feasibility_pass"] for o in expected),
        "failed_sequences_by_check": dict(failures),
        "per_sequence": [
            {
                "id": name,
                "length_nt": len(seq),
                "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest(),
                "native_ms": a * 1000,
                "persistent_ms": b * 1000,
                "outputs": out,
            }
            for (seq, name), a, b, out in zip(records, native, roundtrips, expected)
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    args.output.with_suffix(".fasta").write_text(
        "".join(f">{name}\n{seq}\n" for seq, name in records)
    )
    print(
        json.dumps(
            {
                k: v
                for k, v in report.items()
                if k not in ("per_sequence", "configuration")
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
