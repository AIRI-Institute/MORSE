#!/usr/bin/env python3
"""Build the graphical Colab frontend for the MORSE Python API."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import textwrap
from typing import Any


def _cell(kind: str, cell_id: str, text: str) -> dict[str, Any]:
    metadata = {"id": cell_id}
    if kind == "code":
        metadata["cellView"] = "form"
    result = {"id": cell_id, "cell_type": kind, "metadata": metadata,
              "source": (textwrap.dedent(text).strip() + "\n").splitlines(keepends=True)}
    if kind == "code":
        result.update(execution_count=None, outputs=[])
    return result


def build_notebook() -> dict[str, Any]:
    cells = [
        _cell("markdown", "introduction", """
            # MORSE

            Design an mRNA construct, choose prepared scorer commands and generate
            one complete experiment YAML. Run it with the MORSE Python API in this
            notebook's kernel, follow live progress and inspect saved results.
            Start with the bundled CPU sequence metrics, or connect scorers
            prepared in the kernel's environment.
            """),
        _cell("markdown", "environment_help", """
            ## Environment

            Use the complete MORSE checkout. A fresh VM can clone the configured
            repository's `main` branch, or an exact release commit when supplied.
            Set `MONOREPO_ROOT` to use an existing checkout. Locally install the
            notebook environment with
            `uv sync --all-packages --extra paper --extra mfe --extra notebook --frozen`
            and select it as the kernel. The optional setup cell uses Git and pip
            only to install MORSE. Optimization uses Python calls directly.

            In Colab, mount Google Drive to keep configurations and results after
            the session ends. To use a compute server, open this notebook with a
            Jupyter kernel on that server.
            """),
        _cell("code", "bootstrap_environment", """
            # @title Locate or install this MORSE revision
            MONOREPO_ROOT = ""  # @param {type:"string"}
            MONOREPO_REPOSITORY_URL = "https://github.com/AIRI-Institute/MORSE.git"  # @param {type:"string"}
            MONOREPO_REF = "main"  # @param {type:"string"}
            MONOREPO_COMMIT = ""  # @param {type:"string"}
            INSTALL_NOTEBOOK_DEPENDENCIES = False  # @param {type:"boolean"}
            MOUNT_GOOGLE_DRIVE = False  # @param {type:"boolean"}

            import os
            import re
            import subprocess
            import sys
            from pathlib import Path

            IS_COLAB = "google.colab" in sys.modules
            candidates = [Path(MONOREPO_ROOT).expanduser()] if MONOREPO_ROOT.strip() else []
            candidates += [Path.cwd(), *Path.cwd().parents, Path("/content/morse")]
            ROOT = next((p.resolve() for p in candidates if (p / "core/morse/notebook_tools.py").is_file()), None)
            CLONED_CHECKOUT = False
            if ROOT is None:
                commit = MONOREPO_COMMIT.strip()
                if commit and not re.fullmatch(r"[0-9a-fA-F]{40}", commit):
                    raise ValueError("MONOREPO_COMMIT must be empty or an exact 40-character commit")
                revision = commit or MONOREPO_REF.strip()
                if not MONOREPO_REPOSITORY_URL.strip() or not revision or revision.startswith("-"):
                    raise ValueError("Set the checkout path, or a repository URL and branch/tag or exact commit")
                ROOT = Path("/content/morse") if IS_COLAB else Path.home() / "morse"
                try:
                    subprocess.run(["git", "clone", "--filter=blob:none", "--", MONOREPO_REPOSITORY_URL.strip(), str(ROOT)], check=True)
                except subprocess.CalledProcessError as error:
                    raise RuntimeError("Clone failed. Check repository access or provide an existing checkout through MONOREPO_ROOT.") from error
                subprocess.run(["git", "checkout", "--detach", revision], cwd=ROOT, check=True)
                CLONED_CHECKOUT = True
                INSTALL_NOTEBOOK_DEPENDENCIES = True
            try:
                RESOLVED_MONOREPO_COMMIT = subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True,
                    stderr=subprocess.DEVNULL).strip()
            except (OSError, subprocess.CalledProcessError):
                if CLONED_CHECKOUT:
                    raise RuntimeError("Could not resolve the cloned MORSE commit")
                RESOLVED_MONOREPO_COMMIT = None
            if INSTALL_NOTEBOOK_DEPENDENCIES:
                subprocess.run([
                    sys.executable, "-m", "pip", "install", "-q",
                    "-e", str(ROOT / "core") + "[notebook]",
                    "-e", str(ROOT / "warm-starter"),
                    "-e", str(ROOT / "scorers/manufacturability"),
                    "-e", str(ROOT) + "[mfe]",
                ], check=True)
            STATE_ROOT = Path.home() / "MORSE"
            if IS_COLAB and MOUNT_GOOGLE_DRIVE:
                from google.colab import drive
                drive.mount("/content/drive")
                STATE_ROOT = Path("/content/drive/MyDrive/MORSE")
            STATE_ROOT.mkdir(parents=True, exist_ok=True)
            os.environ.setdefault("MORSE_ASSETS_DIR", str(ROOT / ".cache/morse/assets"))
            print("Checkout:", ROOT, "· commit:", RESOLVED_MONOREPO_COMMIT or "unversioned checkout", "· state:", STATE_ROOT)
            """),
        _cell("code", "load_module", """
            # @title Load the notebook frontend
            import importlib
            import sys
            from IPython.display import display, SVG, FileLink

            core_path = str(ROOT / "core")
            if core_path not in sys.path:
                sys.path.insert(0, core_path)
            loaded = sys.modules.get("morse")
            if loaded is not None and Path(loaded.__file__).resolve().parent != (ROOT / "core/morse").resolve():
                raise RuntimeError("Another MORSE checkout is already imported; restart this kernel")
            import morse
            required_api = ("load_config", "validate", "score", "optimize", "load_run", "resume")
            if not all(callable(getattr(morse, name, None)) for name in required_api):
                raise RuntimeError("This MORSE revision has no Python API. Install the current checkout and restart the kernel.")
            import morse.notebook_tools as notebook_tools
            notebook_tools = importlib.reload(notebook_tools)
            from morse.notebook_tools import (
                ExperimentEditor, render_generated_scope_preview,
                export_configuration_bundle, download_file_to_browser,
            )
            print("MORSE Python API loaded:", morse.__version__)
            """),
        _cell("markdown", "editor_help", """
            ## 1–3. Design, scorers and search

            Enter one full CDS including its ATG start and optional stop, or upload
            a complete RNA and enter its part coordinates. UTR fields accept named
            FASTA choices. Leave absent UTRs and poly(A) blank. Select changes and search
            settings. Add scorer metadata from the local catalog, then enter the
            Bash command prepared on the machine that will execute the run.
            Enter calculation settings in **Parameters** as a YAML mapping, for
            example `max_bp_span: 300` and `compute_ensemble: false` for native MFE.
            Leave it blank to keep the scorer's startup settings. Select a human
            or mouse profile with `profile_path`; bundled defaults are human.

            For each metric choose `optimize` with `min`/`max`, `threshold` with
            bounds or `equals`, or `log`. A command reads one JSON batch from stdin
            followed by EOF and writes one JSON response. Parameters are included
            in that request. Use `curl --data-binary @-`
            for a prepared HTTP service or `docker exec -i` for a running container.

            **Generate YAML** updates the editable full configuration below the
            forms. Schema 4 stores the RNA in FASTA, part coordinates and permitted
            changes in YAML, and scorer outputs as `scorer.metric`. Scorer `cds`
            means the entire CDS; `cds_body` excludes its start and stop.
            Edits to the YAML are used directly by Save/Validate and Launch.
            Advanced YAML supports initialization, including existing warm-start
            manifests, without a second configuration or optimizer implementation.
            """),
        _cell("code", "experiment_editor", """
            # @title Open the experiment editor
            EDITOR = ExperimentEditor(ROOT)
            display(EDITOR.widget)
            """),
        _cell("code", "save_validate_preview", """
            # @title Save the exact YAML, validate it, and preview scopes
            YAML_FILENAME = "experiment.yaml"  # @param {type:"string"}
            import yaml
            run_name = str(yaml.safe_load(EDITOR.yaml.value).get("name", "notebook-run"))
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", run_name):
                raise ValueError("Use letters, numbers, dot, underscore or hyphen for the run name")
            CONFIG_PATH = STATE_ROOT / "experiments" / run_name / Path(YAML_FILENAME).name
            SPEC = EDITOR.save(CONFIG_PATH)
            print("Validated:", CONFIG_PATH, "·", SPEC.name)
            display(FileLink(str(CONFIG_PATH)))
            CONFIG_BUNDLE = export_configuration_bundle(CONFIG_PATH)
            display(FileLink(str(CONFIG_BUNDLE)))
            display(SVG(render_generated_scope_preview(CONFIG_PATH)))
            # Download YAML and FASTA together: download_file_to_browser(CONFIG_BUNDLE)
            """),
        _cell("markdown", "run_help", """
            ## 4. Run with Python

            Save the YAML before running. Validation, scoring and optimization
            read that saved configuration. Scorers execute in this kernel's
            environment using their configured commands.

            The optimization cell remains active and shows progress. Use the
            notebook's **Interrupt** button to stop it. The last completed
            generation can be continued with **Resume** below.
            """),
        _cell("code", "run_settings", """
            # @title Choose live output
            VERBOSE = 2  # @param [0, 1, 2]
            SHOW_METRICS = True  # @param {type:"boolean"}
            SHOW_HYPERVOLUME = True  # @param {type:"boolean"}
            SHOW_UTRS = True  # @param {type:"boolean"}
            PLOTS = [name for name, enabled in (
                ("metrics", SHOW_METRICS), ("hypervolume", SHOW_HYPERVOLUME),
                ("utrs", SHOW_UTRS),
            ) if enabled]
            # verbose: 0 quiet, 1 progress, 2 progress and plots.
            """),
        _cell("code", "preflight", """
            # @title Validate the saved configuration
            PARAMS = morse.load_config(CONFIG_PATH)
            CHECK = morse.validate(params=PARAMS)
            display(CHECK)
            """),
        _cell("code", "score_reference", """
            # @title Optionally score the starting construct
            SCORE_START = False  # @param {type:"boolean"}
            START_SCORES = None
            if SCORE_START:
                PARAMS = morse.load_config(CONFIG_PATH)
                START_SCORES = morse.score(params=PARAMS)
                display(START_SCORES)
            """),
        _cell("code", "launch", """
            # @title Optimize the saved configuration
            from datetime import datetime, timezone

            PARAMS = morse.load_config(CONFIG_PATH)
            stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
            RUN_DIR = STATE_ROOT / "runs" / f"{PARAMS['name']}_{stamp}"
            RESULT = None
            CANDIDATE = None
            print("Results:", RUN_DIR)
            RESULT = morse.optimize(
                params=PARAMS, output_dir=RUN_DIR, verbose=VERBOSE, plots=PLOTS,
            )
            """),
        _cell("code", "resume", """
            # @title Continue an interrupted run
            RESUME_RUN = False  # @param {type:"boolean"}
            RESUME_DIRECTORY = ""  # @param {type:"string"}
            if RESUME_RUN:
                chosen = RESUME_DIRECTORY.strip() or globals().get("RUN_DIR")
                if not chosen:
                    raise ValueError("Set the folder of the run to resume")
                RUN_DIR = Path(chosen).expanduser().resolve()
                RESULT = None
                CANDIDATE = None
                RESULT = morse.resume(RUN_DIR, verbose=VERBOSE, plots=PLOTS)
            """),
        _cell("markdown", "results_help", """
            ## 5. Results

            Reopen a saved run by entering its folder below. Loading and plotting
            use recorded results; they do not rerun scoring. Completed runs have
            a report, Pareto candidates and full RNA sequences.
            """),
        _cell("code", "load_results", """
            # @title Open the current run or a saved results folder
            EXISTING_RUN_DIR = ""  # @param {type:"string"}
            RESULT = None
            CANDIDATE = None
            chosen = EXISTING_RUN_DIR.strip() or globals().get("RUN_DIR")
            if not chosen:
                raise ValueError("Run an experiment or enter its results folder")
            RUN_DIR = Path(chosen).expanduser().resolve()
            RESULT = morse.load_run(RUN_DIR)
            print("Status:", RESULT.summary["status"], "· results:", RESULT.output_dir)
            """),
        _cell("code", "results", """
            # @title Show the saved report and figures
            if (RESULT.output_dir / "report.md").is_file():
                REPORT = RESULT.report()
            else:
                print("This run has no completed report yet. Its saved history is available below.")
            """),
        _cell("code", "plot_results", """
            # @title Plot recorded history
            SHOW_HISTORY_PLOT = False  # @param {type:"boolean"}
            PLOT_KIND = "metrics"  # @param ["metrics", "hypervolume", "utrs"]
            if SHOW_HISTORY_PLOT:
                import matplotlib.pyplot as plt
                figure = RESULT.plot(PLOT_KIND)
                display(figure)
                plt.close(figure)
            # Tables for further analysis: RESULT.pareto, RESULT.history.
            """),
        _cell("code", "select_candidate", """
            # @title Select a Pareto candidate by hypervolume contribution
            CANDIDATE = None
            if RESULT.pareto:
                CANDIDATE = RESULT.select("hvc")
                print("Selected:", CANDIDATE.id, "· RNA:", len(CANDIDATE.sequence), "nt")
                display(CANDIDATE.metrics)
                # Inspect the RNA and its parts: CANDIDATE.sequence, CANDIDATE.parts.
            else:
                print("No feasible Pareto candidate is available in this run.")
            """),
        _cell("code", "export_results", """
            # @title Download the saved results
            import shutil

            RESULT_ZIP = Path(shutil.make_archive(str(RUN_DIR), "zip", root_dir=RUN_DIR))
            display(FileLink(str(RESULT_ZIP)))
            # For a browser download: download_file_to_browser(RESULT_ZIP)
            """),
    ]
    return {"cells": cells, "metadata": {
        "colab": {"name": "MORSE_Server_Colab.ipynb", "provenance": []},
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.10"},
    }, "nbformat": 4, "nbformat_minor": 5}


def rebuild_notebook(source_path: Path, destination: Path) -> None:
    # source_path remains accepted for existing local automation. The builder is
    # the only source of cells; stale frontend cells cannot survive regeneration.
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(build_notebook(), ensure_ascii=False, indent=1) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    default = Path(__file__).resolve().parents[1] / "notebooks/MORSE_Server_Colab.ipynb"
    parser.add_argument("--output", type=Path, default=default)
    args = parser.parse_args()
    rebuild_notebook(args.output, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
