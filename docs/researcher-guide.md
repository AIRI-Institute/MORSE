# Independent experiments

Install from the monorepo root as described in the root README. MORSE searches
sequence designs using prepared inference scorers and fixed checkpoints. Model
training is outside this repository.

Use the [commented YAML template](../examples/experiment.template.yaml) or the
[working MAGE-A1 config](../examples/vaccines/magea1/experiment.yaml).
Schema 5 is the public run format; schemas 3/4 retain their original input
semantics for older studies, including the paper runs.

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

## Python and notebooks

The Python API accepts the same YAML configuration as the CLI:

```python
import morse

params = morse.load_config("experiment.yaml")
morse.validate(params=params)
scores = morse.score(params=params)
result = morse.optimize(params=params, verbose=2)
```

To supply RNA directly, pass it as the first argument to `validate`, `score` or
`optimize`; the configuration still defines the part coordinates. `verbose=1`
shows progress, and `verbose=2` also updates plots inside the notebook.
Results are saved in the configured output directory.

See the [Python API guide](python-api.md) and the
[self-contained CPU notebook](../core/notebooks/MORSE_Python_API.ipynb).
The [graphical notebook](../core/notebooks/MORSE_Server_Colab.ipynb) builds the
YAML in forms and uses the same Python API to validate, score and optimize it
in the active kernel. Interrupt the launch cell to stop; use its Resume cell
to continue from the saved generation.

## Input sequence and assembly

Provide a single-record RNA or DNA FASTA containing the CDS body and any UTRs or
poly(A). The CDS body excludes its initiation and terminal stop codons:

```yaml
schema_version: 5
construct:
  fasta: input.fasta
  coordinates:  # 0-based; end excluded.
    five_prime_utr: [0, 63]
    cds: [63, 414]  # Body only.
    three_prime_utr: [414, 452]
    poly_a: [452, 552]
  start_codon: [AUG]
  stop_codon: [UAA]
  changes:
    cds: synonymous
```

Coordinates describe the supplied FASTA and must cover it consecutively, in
biological order, without gaps or overlaps. Only `cds` is required. For a body-only
FASTA of 351 nt, use `coordinates: {cds: [0, 351]}`. RNA U and DNA T are accepted.
The body must have a length divisible by three and contain no in-frame stops.
An initial AUG in the body is retained: it encodes a methionine after initiation.

MORSE assembles `5′UTR + start + CDS body + stop + 3′UTR + poly(A)` and recalculates
coordinates for each assembled candidate. The default codon lists are `[AUG]` and
`[UAA]`. A single entry fixes the codon; multiple entries allow optimization, with
the first entry defining the reference:

```yaml
start_codon: [AUG, CUG, GUG]
stop_codon: [UAA, UAG, UGA]
```

Put these lists directly under `construct`; `changes` lists the other permitted
changes. Unlisted parts stay fixed:

```yaml
changes:
  cds: synonymous
  five_prime_utr: {library: utr5.fasta}
  three_prime_utr: {library: utr3.tsv, set: curated}
  poly_a: {lengths: [70, 100, 120]}
```

UTR libraries accept multi-record FASTA or the existing parts TSV format.
The input variant stays available as the reference. No UTR or poly(A) is added
unless provided or configured.

Scorer scope `cds` selects the assembled start + body + stop. `cds_body` selects
only the body; `[start_codon, cds_body]` excludes the stop. `full_sequence` selects
the assembled RNA, including poly(A) if present. Scorers receive assembled
coordinates, rather than the input FASTA coordinates.

See [the vaccine examples](../examples/vaccines/) for complete run configs.
Schemas 3/4 keep their previous full-CDS interpretation; do not change only the
version number of an existing configuration.

## Your antigen and UTR candidates

Provide a single-record FASTA with the CDS body, without initiation or stop
codons. If starting from a full annotated CDS, remove its first codon and terminal
stop after checking the annotation. The helper writes the input FASTA and a
schema-5 config with explicit `[AUG]` and `[UAA]` lists:

```sh
morse init --cds antigen.fasta --output-dir experiments/antigen \
  --five-prime-utr utr5_a.fasta --five-prime-utr utr5_b.fasta \
  --three-prime-utr utr3.fasta --seed 43 --population-size 32 --generations 100
morse validate --config experiments/antigen/experiment.yaml
morse preview --config experiments/antigen/experiment.yaml --output runs/antigen_scope.svg
morse pilot --config experiments/antigen/experiment.yaml --samples 24 --output-dir runs/antigen_pilot
morse run --config experiments/antigen/experiment.yaml --output-dir runs/antigen_seed43
```

Repeat UTR flags to supply alternatives. Start and stop are inserted during
assembly. Edit the YAML to select permitted changes, poly(A) alternatives, scorer scopes and metrics. A
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

MORSE routes the requested sequence scope before invoking the command. The
optional `parameters` mapping holds JSON-compatible calculation settings. The adapter validates its supported names, types and ranges.
An omitted or empty mapping sends no `parameters` field. The request can also
contain `sequence_contexts` with physical-part ranges; forward both fields intact
when calling a model that uses them. Send logs to stderr. An error
object such as `{"request_id":1,"error":"message"}`, a nonzero exit, an invalid
response or a timeout stops scoring.

Only explicitly supplied values override the scorer's startup settings;
omitted values retain the selected runtime's defaults. The adapters return their
effective settings in the response's optional `parameters` object; MORSE saves
them to `scorer_parameters.json` during optimization.
See the [parameter reference](scorer-parameters.md) for supported settings and profiles.

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
