#!/usr/bin/env python3
"""Build the self-contained CPU example for the public Python API."""
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
                # MORSE from Python

                Run a small CPU experiment, inspect its Pareto candidates and reopen
                the saved results. The bundled scorer measures CpG content and GC
                distance; these demonstration metrics do not predict expression.

                Install MORSE with the notebook extra following the repository README,
                then select that environment as the notebook kernel. No additional
                models or services are needed for this example.
                """),
            cell("code", "imports", """
                from datetime import datetime, timezone
                from pathlib import Path
                import shlex
                import sys

                import yaml
                import morse

                # A new folder for each independent run.
                stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
                experiment = Path.cwd() / "results" / f"python_api_demo_{stamp}"
                experiment.mkdir(parents=True, exist_ok=False)
                print(experiment)
                """),
            cell("markdown", "input-help", """
                ## Input and configuration

                This short synthetic RNA contains a 5′UTR, CDS, 3′UTR and poly(A).
                Coordinates cover the RNA consecutively and use 0-based, end-exclusive
                ranges. The UTR libraries provide two choices each.
                """),
            cell("code", "input-config", """
                utr5 = "GGGAAACCC"
                cds = "".join(["AUG", "ACU", "UGU", "GGU", "ACU", "AAG",
                               "CUU", "GGU", "AGU", "GGU", "AGU", "UAA"])
                utr3 = "GGAUCUCC"
                rna_seq = utr5 + cds + utr3 + "A" * 12

                (experiment / "sequence.fasta").write_text(f">demo_rna\\n{rna_seq}\\n")
                (experiment / "utr5.fasta").write_text(
                    ">leader_a\\nGGGAAACCC\\n>leader_b\\nGGAACUCC\\n")
                (experiment / "utr3.fasta").write_text(
                    ">trailer_a\\nGGAUCUCC\\n>trailer_b\\nCCUAAGGA\\n")

                cds_start = len(utr5)
                cds_end = cds_start + len(cds)
                utr3_end = cds_end + len(utr3)
                params = {
                    "schema_version": 4,
                    "name": "python_api_demo",
                    "output_dir": "run",
                    "construct": {
                        "fasta": "sequence.fasta",
                        "coordinates": {
                            "five_prime_utr": [0, cds_start],
                            "cds": [cds_start, cds_end],
                            "three_prime_utr": [cds_end, utr3_end],
                            "poly_a": [utr3_end, len(rna_seq)],
                        },
                        "changes": {
                            "cds": "synonymous",
                            "five_prime_utr": {"library": "utr5.fasta"},
                            "three_prime_utr": {"library": "utr3.fasta"},
                        },
                    },
                    "optimizer": {
                        "algorithm": "nsga2", "population_size": 12,
                        "generations": 3, "seed": 43,
                    },
                    "scorers": [{
                        "name": "sequence_metrics",
                        "command": f"{shlex.quote(sys.executable)} -m morse.resources.simple_scorer",
                        "scope": "cds",
                        "metrics": {
                            "cpg_per_kb": {"role": "optimize", "direction": "min"},
                            "gc_distance": {"role": "optimize", "direction": "min"},
                            "gc_fraction": {"role": "log"},
                        },
                    }],
                }
                config_path = experiment / "experiment.yaml"
                config_path.write_text(yaml.safe_dump(params, sort_keys=False))

                # Resolve paths relative to the YAML; the result is an ordinary dict.
                params = morse.load_config(config_path)
                print(f"RNA: {len(rna_seq)} nt; configuration: {config_path}")
                """),
            cell("markdown", "score-help", """
                ## Check the input and initial scores

                Passing `rna_seq` uses this string instead of the configured FASTA.
                The configuration still supplies coordinates and scorer settings.
                """),
            cell("code", "validate", """
                check = morse.validate(rna_seq, params)
                check
                """),
            cell("code", "score", """
                scores = morse.score(rna_seq, params)
                scores["outputs"]
                """),
            cell("markdown", "optimize-help", """
                ## Optimize

                `verbose=2` updates progress and plots inside the cell. Use `1` for
                progress only or `0` for a quiet run. Results remain in the experiment
                folder after the notebook session ends.
                """),
            cell("code", "optimize", """
                result = morse.optimize(
                    rna_seq, params, verbose=2,
                    plots=["metrics", "hypervolume", "utrs"],
                )
                print("Results:", result.output_dir)
                """),
            cell("markdown", "results-help", """
                ## Inspect and select

                The Pareto set contains trade-offs between the objectives. Select one
                candidate by exclusive hypervolume contribution, or by a named metric.
                """),
            cell("code", "pareto", """
                print("Status:", result.summary["status"])
                print("Pareto candidates:", len(result.pareto))
                result.pareto[:3]
                """),
            cell("code", "select", """
                candidate = result.select("hvc")
                print(candidate.id)
                print(candidate.sequence)
                candidate.metrics
                """),
            cell("code", "history", """
                # Lists of dictionaries can also be passed to pandas.DataFrame().
                result.history[-1]
                """),
            cell("code", "utrs", """
                import matplotlib.pyplot as plt
                from IPython.display import display

                figure = result.plot("utrs")
                display(figure)
                plt.close(figure)
                """),
            cell("code", "report", """
                result.report()
                """),
            cell("markdown", "reopen-help", """
                ## Reopen saved results

                Loading results does not rerun scorers. To continue an interrupted run,
                use `morse.resume(result.output_dir, verbose=2)` in the same output folder.
                """),
            cell("code", "reopen", """
                reopened = morse.load_run(result.output_dir)
                print(reopened.summary["status"])
                reopened.pareto[:2]
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
