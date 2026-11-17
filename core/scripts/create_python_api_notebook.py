#!/usr/bin/env python3
"""Build the CPU quickstart notebook for the public Python API."""
from __future__ import annotations

import json
from pathlib import Path
import textwrap


def cell(kind: str, name: str, source: str) -> dict:
    item = {
        "cell_type": kind,
        "id": name,
        "metadata": {},
        "source": (textwrap.dedent(source).strip() + "\n").splitlines(keepends=True),
    }
    if kind == "code":
        item.update(execution_count=None, outputs=[])
    return item


def build_notebook() -> dict:
    return {
        "cells": [
            cell("markdown", "intro", """
                # Run your first MORSE experiment

                Use the same CPU quickstart YAML as the command-line example.
                The search preserves the protein and minimizes two sequence metrics:
                CpG per kilobase and distance from the target GC fraction.
                These metrics do not measure protein expression.

                Open this notebook inside your MORSE checkout and select the
                environment installed with the notebook extra. No GPU, model
                downloads or scorer services are required for this example.
                """),
            cell("code", "imports", """
                from pathlib import Path
                from tempfile import mkdtemp
                from IPython.display import Markdown, display
                import morse

                here = Path.cwd()
                repo = next((folder for folder in (here, *here.parents)
                             if (folder / "examples/quickstart/experiment.yaml").is_file()), None)
                if repo is None:
                    raise FileNotFoundError("Open this notebook from your MORSE checkout.")

                config = repo / "examples/quickstart/experiment.yaml"
                params = morse.load_config(config)
                outputs = repo / "outputs"
                outputs.mkdir(exist_ok=True)
                experiment = Path(mkdtemp(prefix="python_quickstart_", dir=outputs))
                print("Saving this experiment to:", experiment)
                """),
            cell("markdown", "validate-help", """
                ## Check the input

                `validate` checks the FASTA and coordinates without calling a scorer.
                In this config the CDS interval contains the body only; MORSE adds
                the configured start and stop codons during assembly.
                """),
            cell("code", "validate", """
                check = morse.validate(params=params)
                print("Assembled RNA:", check["sequence_length"], "nt")
                print("Allowed changes:", check["mutable_parts"])
                """),
            cell("markdown", "pilot-help", """
                ## Try the scorers before optimizing

                `pilot` scores the starting construct and a small random sample.
                Check that candidates can pass the constraints and that the objectives
                vary. To score just the starting construct, use `morse.score(params=params)`.
                """),
            cell("code", "pilot", """
                pilot = morse.pilot(params=params, samples=8, output_dir=experiment / "pilot")
                print("Starting scores:", pilot["reference"]["outputs"])
                print(f"Feasible fraction in pilot: {pilot['hard_feasible_fraction']:.0%}")
                print("Warnings:", pilot["warnings"])
                """),
            cell("markdown", "optimize-help", """
                ## Run the search

                This example runs eight generations. Use `verbose=2` for live plots,
                `1` for progress only or `0` for a quiet run.
                """),
            cell("code", "optimize", """
                result = morse.optimize(params=params, output_dir=experiment / "run",
                                        verbose=2, plots=["metrics", "hypervolume"])
                print("Pareto candidates:", len(result.pareto))
                """),
            cell("markdown", "results-help", """
                ## Compare candidates with the starting construct

                Select a candidate by its exclusive hypervolume contribution (HVC)
                for a compromise between objectives, or a specialist by one metric.
                The same candidate can win both selections. A negative change is
                favourable for the two minimized objectives; GC fraction is only logged.
                """),
            cell("code", "compare", """
                compromise = result.select("hvc")
                specialist = result.select(metric="sequence_metrics.cpg_per_kb")
                print("HVC:", compromise.id, "| CpG specialist:", specialist.id)
                chosen = list({item.id: item for item in (compromise, specialist)}.values())
                rows = result.compare(chosen, baseline=result.reference)

                table = ["| Candidate | Metric | Start | Final | Change |",
                         "|---|---|---:|---:|---:|"]
                for row in rows:
                    table.append(f"| {row['candidate_id']} | {row['metric']} | "
                                 f"{row['baseline']:.3g} | {row['value']:.3g} | {row['delta']:+.3g} |")
                display(Markdown("\\n".join(table)))
                """),
            cell("markdown", "reopen-help", """
                ## Reopen the report and save the selected RNA

                Loading a completed run does not recalculate scores. Retrieve a
                candidate by the ID shown in the report or Pareto table, then
                export its complete assembled RNA to FASTA.
                """),
            cell("code", "report", """
                reopened = morse.load_run(result.output_dir)
                reopened.report()
                """),
            cell("code", "export", """
                selected = reopened.get_candidate(compromise.id)
                fasta = experiment / "selected.fasta"
                selected.to_fasta(fasta)
                print(f"Saved {selected.id} ({len(selected.sequence)} nt):", fasta)
                """),
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.13"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def main() -> None:
    path = Path(__file__).resolve().parents[1] / "notebooks" / "MORSE_Python_API.ipynb"
    path.write_text(json.dumps(build_notebook(), indent=1, ensure_ascii=False) + "\n")
    print(path)


if __name__ == "__main__":
    main()
