#!/usr/bin/env python3
"""Build the small YAML-first Colab frontend deterministically."""
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
            one complete experiment YAML. The same YAML runs through Core locally
            or through your authenticated server backend. This notebook shows
            optimizer progress and the plots produced by Core.

            Scorer scripts, services, containers and checkpoints must already be
            prepared in the execution environment. The notebook does not install
            models. Start with the built-in CPU sequence metrics to check the flow.
            """),
        _cell("markdown", "environment_help", """
            ## Environment

            Use the complete MORSE checkout. A fresh VM can clone the configured
            repository's `main` branch, or an exact release commit when supplied.
            For this private repository, use an authenticated Git environment or
            provide an existing checkout through `MONOREPO_ROOT`. Keep credentials
            out of notebook fields. Locally install the optional
            widgets with `uv sync --all-packages --extra paper --extra mfe --extra notebook --frozen`.
            The optional installer prepares MORSE, notebook widgets and native
            CPU scorers from this checkout.
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
                    raise RuntimeError("Clone failed. Use an authenticated Git environment for the private repository, or provide an existing checkout through MONOREPO_ROOT.") from error
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
            import morse.notebook_tools as notebook_tools
            notebook_tools = importlib.reload(notebook_tools)
            from morse.notebook_tools import (
                ExperimentEditor, LocalMorseClient, MorseBackendClient,
                colab_secret_or_prompt, require_backend_compatibility,
                read_configuration, render_generated_scope_preview, monitor_job,
                export_configuration_bundle,
                extract_artifacts, display_result_artifacts, download_file_to_browser,
            )
            print("Frontend loaded. Scorers are inspected as metadata; no model is loaded.")
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
            ## 4. Run and results

            Local mode uses the existing Core job manager and runs `python -m morse
            run --config ... --output-dir ...`. Server mode sends the complete YAML
            mapping and its small sequence/library inputs to backend API 6. Scorer commands resolve on the execution
            host, so use explicit paths or installed commands for remote services.
            The backend token authorizes execution of your configured Bash commands.
            Keep it in a Colab secret or the hidden prompt, never in a saved cell.
            """),
        _cell("code", "connect_executor", """
            # @title Connect to the execution environment
            EXECUTION_MODE = "Local"  # @param ["Local", "Server backend"]
            BACKEND_URL = "http://127.0.0.1:8765"  # @param {type:"string"}
            if EXECUTION_MODE == "Server backend":
                CLIENT = MorseBackendClient(BACKEND_URL, colab_secret_or_prompt())
            else:
                # Reusing the manager keeps existing local jobs running when this
                # cell is rerun, or after switching back from a remote executor.
                if "LOCAL_CLIENT" not in globals():
                    LOCAL_CLIENT = LocalMorseClient(ROOT / "core", STATE_ROOT)
                CLIENT = LOCAL_CLIENT
            require_backend_compatibility(CLIENT.health())
            print("Connected:", EXECUTION_MODE)
            """),
        _cell("code", "preflight", """
            # @title Preflight the saved YAML on the execution host
            CONFIGURATION = read_configuration(CONFIG_PATH)
            PREFLIGHT = CLIENT.preflight_job(configuration=CONFIGURATION)
            print({key: PREFLIGHT.get(key) for key in (
                "valid", "run_name", "algorithm", "population_size", "generations",
                "objective_count", "threshold_count", "configuration_fingerprint")})
            display(SVG(PREFLIGHT["sequence_scope_svg"]))
            """),
        _cell("code", "launch", """
            # @title Launch the saved YAML
            # Read the saved snapshot again, so unsaved form changes cannot alter it.
            CONFIGURATION = read_configuration(CONFIG_PATH)
            PREFLIGHT = CLIENT.preflight_job(configuration=CONFIGURATION)
            JOB = CLIENT.submit_job(configuration=CONFIGURATION)
            JOB_ID = JOB["id"]
            print("Job:", JOB_ID, "·", JOB["status"])
            """),
        _cell("code", "progress", """
            # @title Track progress, a live curve and recent logs
            WATCH_SECONDS = 300  # @param {type:"integer"}
            POLL_SECONDS = 3  # @param {type:"number"}
            JOB = monitor_job(CLIENT, JOB_ID, poll_seconds=POLL_SECONDS, max_seconds=WATCH_SECONDS)
            print("Job:", JOB_ID, "·", JOB["status"])
            # Rerun this cell to keep watching. Interrupting it only stops observation.
            """),
        _cell("code", "job_controls", """
            # @title Reconnect, cancel, or resume an existing job
            ACTION = "List jobs"  # @param ["List jobs", "Reconnect", "Cancel", "Resume"]
            EXISTING_JOB_ID = ""  # @param {type:"string"}
            if ACTION == "List jobs":
                display(CLIENT.jobs())
            else:
                chosen = EXISTING_JOB_ID.strip() or globals().get("JOB_ID", "")
                if not chosen:
                    raise ValueError("Set an existing job ID")
                JOB_ID = chosen
                if ACTION == "Cancel":
                    display(CLIENT.cancel(JOB_ID))
                elif ACTION == "Resume":
                    display(CLIENT.resume_job(JOB_ID))
                else:
                    display(CLIENT.job(JOB_ID))
            """),
        _cell("code", "results", """
            # @title Show Core plots and report, and download results
            RESULT_ROOT = STATE_ROOT / "downloads" / JOB_ID
            RESULT_ZIP = CLIENT.download_artifacts(JOB_ID, RESULT_ROOT.with_suffix(".zip"))
            ARTIFACTS = extract_artifacts(RESULT_ZIP, RESULT_ROOT)
            display_result_artifacts(ARTIFACTS)
            display(FileLink(str(RESULT_ZIP)))
            # For a browser download, run: download_file_to_browser(RESULT_ZIP)
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
