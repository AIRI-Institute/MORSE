# Research configuration examples

These schema-5 examples declare prepared Bash scorer commands, calculation
`parameters` and `metrics`.
They preserve the example objectives, thresholds, scopes, seeds and search
settings. They are independent research templates, not the frozen paper runs.
See `paper/` for the publication's exact protocols.

The dependency-free process example can run immediately after installation:

```sh
morse run --config core/configs/mlana_cds_single_cpg_ga.yaml
```

The remaining examples use native MFE and prepared containers for other models.
Install native MFE with `uv sync --all-packages --extra paper --extra mfe --frozen`.
Its YAML `parameters` specifies `max_bp_span: 300` and `compute_ensemble: false`.
Container commands use `docker exec -i`; they never create containers, install
dependencies or download checkpoints. Their setup settings remain in `docker exec -e`.

| Container variable | Default name | Script already present inside it |
|---|---|---|
| `MORSE_CAI_CONTAINER` | `morse-cai` | `/app/server.py` |
| `MORSE_TRANSLATION_EFFICIENCY_CONTAINER` | `morse-translation-efficiency` | `/app/server.py` |
| `MORSE_HALF_LIFE_CONTAINER` | `morse-half-life` | `/app/server.py` |
| `MORSE_IMMUNOGENICITY_CONTAINER` | `morse-immunogenicity` | `/app/server.py` |

Use the scorer README/manifest and `docs/scorers.md` to prepare images, model
mounts, dependency versions and GPU visibility. Checkpoint paths in these
commands are container paths. A prepared container must remain running and
accept `docker exec`; changing the command to your own HTTP client is also
supported. Invoking a Python model script per batch reloads that script; a
resident model service can keep the model loaded.

`MORSE_PYTHON` selects the Python environment for local process commands and
defaults to MORSE's interpreter. Shell variable references inside `command`
are evaluated by Bash, with normal quoting, when the command runs.

Inspect an interface or generate a new section without loading a model:

```sh
morse scorers describe cai
morse scorers template cai --command 'docker exec -i morse-cai python3 /app/server.py'
```

No `mode`, image lifecycle or command argument list is needed in an active run
configuration. Scorer setup manifests retain their own installation metadata.

Each config uses `construct.fasta` and `coordinates` for the input before codon
insertion. `input.fasta` omits the start and stop present in the reference
`sequence.fasta`; the codon lists restore the same assembled reference RNA.
`construct.changes` defines other permitted variants and `output_dir` sets results.
Names such as `cds_body` in scorer scopes explicitly exclude start and stop.
