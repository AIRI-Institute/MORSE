"""Native VaxPress search with an experiment-owned exact event-budget boundary.

All search and selection are inherited. Only the final proposed batch is
truncated, in native proposal order, to enforce the declared scoring budget.
"""
from __future__ import annotations

import argparse
import json
import platform
import time
from pathlib import Path

import RNA
from vaxpress.evolution_chamber import CDSEvolutionChamber, ExecutionOptions
from vaxpress.log import initialize_logging
from shared_scores import SharedFitness, sha256


class BudgetedChamber(CDSEvolutionChamber):
    def __init__(self, *args, event_budget, budget_log, **kwargs):
        self.event_budget, self.budget_log = event_budget, Path(budget_log)
        self.requested_events = 0
        super().__init__(*args, **kwargs)

    def mutate_population(self, iteration):
        remaining = self.event_budget - self.requested_events
        if remaining == 0:
            raise StopIteration
        super().mutate_population(iteration)
        proposed = len(self.flatten_seqs)
        keep = min(proposed, remaining)
        self.population[:] = self.population[:keep]
        self.population_sources[:] = self.population_sources[:keep]
        self.flatten_seqs[:] = self.flatten_seqs[:keep]
        self.requested_events += keep
        with self.budget_log.open("a") as handle:
            handle.write(json.dumps({"iteration": iteration + 1,
                "proposed_events": proposed, "scored_events": keep,
                "discarded_unscored_proposals": proposed - keep,
                "cumulative_events": self.requested_events}) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--alpha", type=float, required=True)
    parser.add_argument("--budget", type=int, default=2560)
    parser.add_argument("--processes", type=int, default=4)
    args = parser.parse_args()
    if args.budget < 1 or not 1 <= args.processes <= 4:
        raise ValueError("Invalid budget or process count")
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    log_path, budget_log = out / "evaluations.jsonl", out / "budget_batches.jsonl"
    log_path.touch(exist_ok=False)
    budget_log.touch(exist_ok=False)
    cfg = json.loads(Path(args.config).read_text())
    if args.alpha not in cfg["alphas"] or args.seed not in cfg["seeds"]:
        raise ValueError("Seed or weight differs from frozen protocol")
    initialize_logging(str(out / "log.txt"), True)
    opts = ExecutionOptions(
        n_iterations=100000, n_population=20, n_survivors=2,
        initial_mutation_rate=0.1, winddown_trigger=15, winddown_rate=0.9,
        output=str(out), command_line="experiment run_vaxpress_shared.py",
        overwrite=False, seed=args.seed, processes=args.processes,
        random_initialization=False, conservative_start=None,
        boost_loop_mutations="1.5:15", full_scan_interval=300,
        species="Homo sapiens", codon_table="standard", protein=False,
        quiet=True, seq_description="MLANA fixed UTR shared MFE/nt and CSI",
        print_top_mutants=0, addons=[str(Path(__file__).with_name("shared_scores.py"))],
        lineardesign_dir=None, lineardesign_lambda=None, lineardesign_omit_start=5,
        folding_engine="vienna")
    scoreopts = {"shared": {"weight": 1.0, "alpha": args.alpha,
        "config": str(Path(args.config).resolve()), "log_path": str(log_path)}}
    with (out / "experiment_provenance.json").open("x") as handle:
        json.dump({"arguments": vars(args), "protocol": cfg,
            "protocol_sha256": sha256(args.config), "python": platform.python_version(),
            "viennarna": RNA.__version__, "shared_scores_sha256": sha256(
                Path(__file__).with_name("shared_scores.py")),
            "driver_sha256": sha256(__file__), "active_score_names": ["shared"],
            "budget_rule": "all scored CDS events including duplicates; final batch prefix"},
            handle, indent=2)
    started = time.monotonic()
    chamber = BudgetedChamber(cfg["initial_cds"], {"shared": SharedFitness},
        scoreopts, opts, event_budget=args.budget, budget_log=budget_log)
    last_status = None
    for last_status in chamber.run():
        pass
    if last_status is None or last_status["error"]:
        raise RuntimeError("VaxPress scoring failed; retained output contains diagnostics")
    chamber.save_results()
    rows = [json.loads(line) for line in log_path.read_text().splitlines()]
    if len(rows) != chamber.requested_events:
        raise AssertionError("Evaluation logging/count mismatch")
    if any(row["csi_codons"] != 117 or row["transcript_body_length"] != 470
           for row in rows):
        raise AssertionError("Scoring scope changed")
    summary = {"seed": args.seed, "alpha": args.alpha,
        "events": len(rows), "unique_cds": len({row["cds"] for row in rows}),
        "budget": args.budget, "budget_reached": len(rows) == args.budget,
        "stop_reason": "event_budget" if len(rows) == args.budget else "native_winddown",
        "elapsed_seconds": time.monotonic() - started,
        "best_fitness": max(row["fitness"] for row in rows),
        "best_mfe_per_nt": min(row["mfe_per_nt"] for row in rows),
        "best_csi_percentile": max(row["csi_percentile"] for row in rows)}
    with (out / "summary.json").open("x") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
