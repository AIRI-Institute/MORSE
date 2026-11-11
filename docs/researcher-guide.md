# Independent experiments

Install from the monorepo root as described in the root README. MORSE searches
sequence designs using prepared inference scorers and fixed checkpoints. Model
training is outside this repository.

Use the [commented YAML template](../examples/experiment.template.yaml) or the
[working MAGE-A1 config](../examples/vaccines/magea1/experiment.yaml).
Schema 4 is the public run format; schema 3 remains readable for older studies.

## A CPU example

```sh
morse run --config examples/quickstart/experiment.yaml
```

The bundled CPU scorer needs no additional model or service. This real search
preserves the MLANA protein and minimizes CpG per kb and distance from 50% GC in
start codon + CDS. These demonstration metrics do not predict vaccine efficacy
or expression. Inspect `report.md`, `pareto.csv` and `pareto.fasta` in the output.

Colours follow scorer order in the configuration; outputs of one scorer share
its colour. Charts are saved as PNG and vector PDF. With at least three optimized
numeric metrics, the report includes a shared half-spider plot with direction-
aware min–max scaling over the final Pareto set and initial sequence.

```sh
morse report outputs/quickstart
```

## A complete RNA or CDS-only input

Provide a single-record RNA or DNA FASTA and define part coordinates in the YAML:

```yaml
construct:
  fasta: sequence.fasta
  coordinates:
    five_prime_utr: [0, 17]
    cds: [17, 947]
    three_prime_utr: [947, 1011]
    poly_a: [1011, 1081]
  changes:
    cds: synonymous
```

Coordinates are **0-based, end exclusive**: `[0, 17]` selects the first 17
nucleotides. The file path is relative to the config. Only `cds` is required;
omit other parts when they are absent. For a 930-nt CDS-only FASTA, use:

```yaml
construct:
  fasta: cds.fasta
  coordinates:
    cds: [0, 930]
  changes:
    cds: synonymous
```

The provided ranges must cover the FASTA consecutively, in biological order,
without gaps or overlaps. `cds` here includes the ATG and an optional terminal
stop; MORSE splits these internally. It does not add missing UTRs, a stop or
poly(A). RNA U is accepted and normalized to T.

`construct.changes` lists what may change. Unlisted parts stay fixed:

```yaml
changes:
  cds: synonymous
  five_prime_utr: {library: utr5.fasta}
  three_prime_utr: {library: utr3.tsv, set: curated}
  poly_a: {lengths: [70, 100, 120]}
```

UTR libraries accept multi-record FASTA or the existing parts TSV format.
The input variant stays available as the reference. Coordinates describe that
input; candidate coordinates are recalculated when UTR or tail lengths change.

Scope `cds` selects the full CDS, including its start and optional stop.
`cds_body` excludes both; `[start_codon, cds_body]` excludes only the stop.
Use `full_sequence` for the whole input, including poly(A) if present.

See [the vaccine examples](../examples/vaccines/) for complete run configs.
The old `construct.parts` format remains available through schema-3 configs.

## Your antigen and UTR candidates

Provide a single-record FASTA with the full coding sequence, beginning with ATG
and optionally ending with TAA/TAG/TGA. RNA U is accepted and normalized to T.
Internal stops, malformed frame and multiple FASTA records are rejected. The
helper writes the RNA FASTA and a schema-4 config into a new experiment directory:

```sh
morse init --cds antigen.fasta --output-dir experiments/antigen \
  --five-prime-utr utr5_a.fasta --five-prime-utr utr5_b.fasta \
  --three-prime-utr utr3.fasta --seed 43 --population-size 32 --generations 100
morse validate --config experiments/antigen/experiment.yaml
morse preview --config experiments/antigen/experiment.yaml --output runs/antigen_scope.svg
morse pilot --config experiments/antigen/experiment.yaml --samples 24 --output-dir runs/antigen_pilot
morse run --config experiments/antigen/experiment.yaml --output-dir runs/antigen_seed43
```

Repeat UTR flags to supply alternatives. Omitted parts are not added. Edit the YAML
to select permitted changes, poly(A) alternatives, scorer scopes and metrics. A
metric's `role` is `optimize`, `threshold` or `log`. Optimization metrics specify
`direction: min|max`; thresholds specify `min`, `max` or `equals`. Numeric float
values are the default; other values can declare `value_type`.

Use a new output directory and a recorded seed for each independent search.
The initializer refuses to overwrite an existing experiment directory.

`output_dir` in YAML is relative to the config, as are input and library paths.
`--output-dir` overrides it for an invocation. Run the saved experiment with
`morse run --config experiment.yaml`; no additional arguments are needed.

To run a prepared collection of configurations, give each YAML its intended
seed and use a short Bash loop. This continues after failures and returns a
nonzero status if any run failed:

```bash
failed=0
for config in experiments/*.yaml; do
  filename=${config##*/}
  morse run --config "$config" --output-dir "outputs/${filename%.yaml}" || failed=1
done
test "$failed" -eq 0
```

## Prepare a scorer, then configure its command

Install scripts, checkpoints, dependencies and containers, or start a model
service, before `morse run`. The optimizer invokes the Bash command you provide
for each batch. It does not prepare or manage the external container or service.
A command can call a local script, `docker exec -i`, `curl`, or a pipeline.

Inspect bundled interfaces without loading a model:

```sh
morse scorers
morse scorers describe cai
morse scorers describe cai --json
morse scorers template cai --command 'docker exec -i morse-cai python3 /app/server.py'
morse scorers describe --manifest /path/to/morse-scorer.yaml
morse scorers template --manifest /path/to/morse-scorer.yaml --command './my-scorer-client.sh'
```

Descriptions show metric names, types, scopes and calculation parameters. Templates provide a run YAML
section for your command; review its scope and metric roles. Scorer manifests
are installation/interface metadata, not a request to launch an image during
optimization. Use each scorer's README and the asset registry to prepare its
real checkpoint and runtime. Never substitute a demonstration metric for an
unavailable predictor.

## Bash command and JSON contract

Each invocation receives one JSON object on stdin, followed by EOF. Return one
JSON object on stdout with the same `request_id` and one output dictionary per
input sequence, in the same order:

```json
{"request_id":1,"sequences":["ATGGCTCGT"],"parameters":{"target_gc":0.5}}
```

```json
{"request_id":1,"outputs":[{"cpg_per_kb":111.11111111111111,"gc_fraction":0.5555555555555556,"gc_distance":0.05555555555555558}]}
```

MORSE routes the requested sequence scope before invoking the command. Optional
`parameters` comes from the scorer's YAML section and holds JSON-compatible
calculation settings. The adapter validates its supported names, types and ranges.
An omitted or empty mapping sends no `parameters` field. The request can also
contain `sequence_contexts` with physical-part ranges; forward both fields intact
when calling a model that uses them. Send logs to stderr. An error
object such as `{"request_id":1,"error":"message"}`, a nonzero exit, an invalid
response or a timeout stops scoring.

A local script, relative to the configuration directory:

```yaml
scorers:
  - name: sequence_metrics
    command: '"${MORSE_PYTHON}" scorer.py'
    parameters:
      target_gc: 0.5
    provenance_paths: [scorer.py]
    scope: [start_codon, cds_body]
    metrics:
      cpg_per_kb: {role: optimize, direction: min}
      gc_fraction: {role: threshold, min: 0.35, max: 0.65}
      gc_distance: {role: log}
```

`MORSE_PYTHON` defaults to the interpreter running MORSE; set it explicitly for a
separate scorer environment. `workdir` is optional and defaults to the config
folder. A command remains a Bash string: variables, quotes and pipelines are
interpreted by Bash at invocation time, rather than split into configuration
arguments.

Each scorer name is unique. Metrics are stored as `scorer.metric`; two scorers
may both return `score` without collision. Missing required outputs, invalid
numeric values, command errors and timeouts stop the run. Numeric metrics and
threshold bounds must be finite; NaN and infinity are rejected.

For an already running container with its model and files prepared:

```yaml
command: 'docker exec -i "${MY_SCORER_CONTAINER:-my-scorer}" python3 /app/server.py'
```

For a prepared HTTP service that accepts and returns the same JSON protocol:

```yaml
command: >-
  curl --fail --silent --show-error
  --header 'Content-Type: application/json'
  --data-binary @-
  "${MY_SCORER_URL:?Set MY_SCORER_URL before the run}"
```

Use `-i` for Docker stdin and `--data-binary @-` for the whole HTTP request body.
The response must retain `request_id` and output order. A local executable client
can also handle SSH, authentication or another transport. Copy
`examples/custom_scorer/scorer.py` for a minimal working implementation.

A script is invoked anew for each batch. A resident model service can retain its
loaded model across requests. Choose the client command to match your prepared
runtime; record its model/checkpoint identity and precision.

## Provenance, inspection and continuation

`provenance_paths` list the local source modules, weights and auxiliary files
that your command reads, relative to the experiment/input root. MORSE records
the literal command and declared local dependencies. Arbitrary shell commands
and remote services are not introspected to discover hidden dependencies: keep
external image/model versions and preparation receipts with the experiment.

Every run records its configuration, reference construct, selected sequences,
native scores, validation/threshold results, score history and report. Numeric
CSV/JSON artifacts are the basis for downstream analysis. A seed alone does not
specify an experiment; preserve the exact inputs, model versions and hardware.

```sh
morse report runs/antigen_seed43
morse run --config experiments/antigen/experiment.yaml --output-dir runs/antigen_seed43 --resume
```

Interrupted continuation restores optimizer and random state and checks recorded
runtime compatibility. Start a new run if a scorer or external model changed.
A completed run can initialize a new continuation with `--resume-from` and a
separate output directory. Results remain computational design candidates.
