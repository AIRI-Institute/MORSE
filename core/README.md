# MORSE Core

Core optimizes non-replicating RNA constructs, including CDS-only inputs, while
preserving the encoded protein. GA, NSGA-II and NSGA-III explore UTR alternatives,
synonymous CDS choices, start/stop choices and poly(A) lengths. Independent inference
processes receive only their declared sequence scopes. Objectives, feasibility
thresholds and logged measurements stay separate.

## Install and run

From the monorepo root, install with the root instructions or:

```sh
python -m pip install -e ./core -e ./warm-starter
morse-core run --config examples/quickstart/experiment.yaml --output-dir runs/quickstart
```

The quickstart requires CPU only, no checkpoints and no Docker. Its transparent
CpG/GC metrics demonstrate the interface, not a biological performance model.
See [the researcher guide](../docs/researcher-guide.md) for an experiment using
your own FASTA, UTRs and scorer. Models run inference from ready checkpoints;
model training is excluded from this public repository.

The engine package provides `morse-core` and `morse-backend`. The complete
monorepo installation also provides `morse`, which adds the paper and asset
commands and delegates experiment execution to the same Core engine.

```sh
morse-core init --cds my_antigen.fasta --output-dir experiments/my_antigen
morse-core validate --config experiments/my_antigen/experiment.yaml
morse-core preview --config experiments/my_antigen/experiment.yaml --output runs/preview.svg
morse-core pilot --config experiments/my_antigen/experiment.yaml --samples 24 --output-dir runs/pilot
morse-core run --config experiments/my_antigen/experiment.yaml --output-dir runs/my_antigen
morse-core report runs/my_antigen
```

`init` accepts a full single-record CDS beginning with ATG, with an optional
terminal stop. It writes schema-4 YAML, an RNA FASTA and a short README.
Missing UTRs, stop and poly(A) are not added. The CLI also reads schema-3 configs
with their original part and scope semantics.

## Experiment configuration

YAML owns all scientific choices. Relative FASTA/profile/initial-population paths
and explicit scorer `workdir` resolve from its directory (or `--input-root`).
Use the same configuration on a workstation, SLURM node or your own backend.
`construct.fasta` and `construct.coordinates` define the input with 0-based,
end-exclusive ranges. Only the full `cds` is required; UTRs and poly(A) are
optional. `construct.changes` lists what may change; unlisted parts stay fixed.
See the [commented template](../examples/experiment.template.yaml) and
[researcher guide](../docs/researcher-guide.md).
See `examples/quickstart/experiment.yaml` and the retained configs under
`core/configs` for complete examples. Model examples need their documented
inference dependencies and model directories. Prepare scripts, containers and
services before a run. `morse-core validate` checks the config and construct inputs.

The recommended scorer interface is a Bash command string:

```yaml
scorers:
  - name: sequence_metrics
    command: '"${MORSE_PYTHON}" scorer.py'
    scope: [start_codon, cds_body]
    metrics:
      cpg_per_kb: {role: optimize, direction: min}
      gc_fraction: {role: threshold, min: 0.35, max: 0.65}
      gc_distance: {role: log}
```

A metric has role `optimize` with direction `min` or `max`, `threshold`
with `min`/`max`/`equals`, or `log`. Scope `cds` selects the full CDS;
`cds_body` excludes start and stop. Lists, slices and anchored windows select
other regions. Poly(A) enters a scorer only when its scope includes it.
Outputs are named `scorer.metric` automatically. Numeric metrics must be finite.
Floats are the default; `value_type` can declare other types.

Each invocation receives one JSON object through stdin and EOF, then returns
one JSON object with the same `request_id` and ordered `outputs`. Optional
`sequence_contexts` retain physical-part ranges. Logs belong on stderr.
`MORSE_PYTHON` defaults to Core's interpreter; command variables and pipelines
are interpreted by Bash. The command may call a script, `docker exec -i` or an
HTTP client using `curl --data-binary @-`. The external scorer runtime is
prepared separately; no runtime mode or image is required in the run YAML.

Manufacturability is configurable as a local hard validator, a versioned
sequence-feasibility scorer, or explicit thresholds. Temporary adaptive epsilon
is a search mechanism; final hard feasibility is recorded separately. No generic
profile guarantees synthesis, IVT performance, expression or clinical utility.

## Provenance and resume

Runs save original/resolved YAML, reference and Pareto FASTA/CSV, the feasible
archive, structural/threshold observations, reports, plots, generation history,
compressed JSONL logs and an atomic optimizer checkpoint. `run_provenance.json`
records Core source content, Python/package versions and configured scorer
commands. Declare every local source dependency, lookup table and model file
in `scorers[].provenance_paths`; these inputs are hashed by content before
inference and included in continuation compatibility checks. Arbitrary shell
commands and external services are not parsed to discover hidden dependencies.
Record their image/model identities and preparation receipts separately. Legacy
managed runtimes additionally record their configured image/mount identities.
Content checks can take time for large checkpoints.

```sh
morse-core run --config experiment.yaml --output-dir runs/interrupted --resume
morse-core run --config continuation.yaml --output-dir runs/continued --resume-from runs/completed
```

An interrupted run resumes its exact optimizer/random state in its original
output directory. A completed run seeds a continuation from its final population;
that is a distinct operation. Only load locally trusted checkpoint files: Python
pickle is an executable serialization format. Exact continuation requires the
same compatible runtime and model content; historical run records without this
manifest retain their original provenance limits.

## Scorer extensions and optional notebook

[The minimal JSONL scorer](../examples/custom_scorer) is the shortest extension
path. Each full scorer under `scorers/` supplies a schema-2 `morse-scorer.yaml`
with setup metadata, parameters, scope, metric policy and a
self-test. In the complete CLI, `morse scorers describe ID` inspects that
interface and `morse scorers template ID --command 'YOUR_COMMAND'` generates a
section for your prepared client. Custom interfaces use `--manifest PATH`.
Learned predictors need ready weights;
see [assets](../assets) and [scorers](../scorers).

`notebooks/MORSE_Server_Colab.ipynb` is an optional interface to the same CLI.
It generates and exports schema-4 YAML, launches and tracks jobs, and displays
their plots and artifacts. Scorer commands use an already prepared environment;
the notebook does not install scorer runtimes or start model services. The
notebook builder is `scripts/update_colab_notebook.py`; edit the builder and
regenerate the notebook. Your own backend is optional; see
[jobs backend](docs/colab-backend.md).

Local mode uses the same job manager as the optional remote jobs API. Remote
mode submits the complete YAML configuration with a trusted execution token.

## Tests

Unit-test sources are local maintainer material and are excluded from the public
checkout. When those sources are present, run from `core/` after installing both packages:

```sh
python -m unittest discover -v
```

Core tests cover representation/translation, scopes, constraints, objectives,
plugin contracts, backend/notebook behavior and exact interrupted NSGA-III
continuation. Root CI also runs independent experiment and scorer smoke checks.
