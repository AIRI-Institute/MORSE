# MORSE Warm Starter

Standalone CLI for generating a scorer-guided first population for
[MORSE Core](../core).

Core owns validation, thresholds and optimization. Warm Starter owns only anchor
search and writes the versioned JSON manifest consumed by
`optimizer.initial_population`. A schema-5 run can also declare `initialization`
directly; `morse run` invokes this package before the first population, so the
same YAML works from Colab, a workstation or a cluster CLI. See the commented
[run template](../examples/experiment.template.yaml).

## Installation

Use the shared monorepo environment from the [root setup](../README.md).
Warm Starter is included; it needs the same prepared scorers as the run.

## CLI

The CLI reads the existing Core run config to reuse its construct design space,
prepared Bash scorer commands, input scopes and optimization directions:

From the monorepo root, this lightweight example uses the quickstart scorer:

```sh
uv run morse-warm-start generate \
  --config examples/quickstart/experiment.yaml \
  --output outputs/initial-population.json \
  --greedy sequence_metrics:cpg_per_kb:3 \
  --population-size 16 --seed 43
```

`population-size` and `seed` default to the values in the Core run config. The
reference construct is included by default and the remaining slots are filled
with reproducible random design vectors. Use `--no-original` to omit the
reference, or `--anchors-only` to let Core fill all missing/rejected slots.

The generated file is then the only warm-start setting in Core:

```yaml
optimizer:
  initial_population: /absolute/path/to/initial-population.json
```

The file must have been generated for the same construct and design space.

## LinearDesign

LinearDesign warm starts optimize synonymous codons against a configurable CAI
profile and folding energy. The upstream source is pinned to commit
`f0126ca89a8b853088b4bccfd2cc8c378d3678be`; review its academic/non-profit
license before setting the explicit acceptance flag.

    initialization:
      mode: linear_design
      profile_path: profiles/mouse.json
      lambdas: [0.0, 1.0, 2.0]
      lineardesign_dir: .runtime/LinearDesign
      accept_lineardesign_license: true

`profile_path` accepts a MORSE schema-1 codon-count JSON or a LinearDesign
CSV/TSV. Lambda is repeatable in the standalone command:

    morse-warm-start generate \
      --config run.yaml \
      --output initial-population.json \
      --anchors-only \
      --lineardesign-profile profiles/mouse.json \
      --lineardesign-lambda 0 \
      --lineardesign-lambda 1 \
      --accept-lineardesign-license

## Anchor modes

Warm Starter has two scorer-guided modes plus the reference and random fillers:

- `--exhaustive SCORER:OBJECTIVE[:COUNT]` enumerates every mutable design
  variable intersecting the selected scorer's input scope, then keeps the best
  candidates in the objective's configured `min`/`max` direction. Use it for
  small scopes or libraries whose full combination count is practical.
- `--greedy SCORER:OBJECTIVE[:COUNT]` performs coordinate descent over mutable
  variables in the scorer scope. This is suitable for maximum
  `start_codon + cds_body` CAI and
  minimum full-sequence `mfe_per_nt`, and also works for other configured
  objectives such as translation efficiency or half-life.
- The original construct is an anchor unless `--no-original` is passed.
- Random candidates fill the population unless `--anchors-only` is passed.

The main controls are `--max-combinations` and `--batch-size` for exhaustive
search, `--passes` for greedy search, and `--candidate-pool` plus
`--min-hamming-distance` for diversity. Run `morse-warm-start generate --help`
for defaults.

Only the selected scorer command is invoked for each anchor request. Prepare
the script, model files, container or service before generation. Warm Starter
accepts the same command strings and metric declarations as Core:

```yaml
scorers:
  - name: cai_human
    command: 'docker exec -i "${MORSE_CAI_CONTAINER:-morse-cai}" python3 /app/server.py'
    scope: [start_codon, cds_body]
    metrics:
      cai: {role: optimize, direction: max}
```

Each Bash invocation receives one JSON batch and EOF, and must return one JSON
response with the same `request_id` and ordered `outputs`. Logs belong on stderr.
The command may call a local script, `docker exec -i` or a prepared HTTP service
using `curl --data-binary @-`. It does not prepare or manage the external runtime.
`MORSE_PYTHON` defaults to the interpreter running MORSE when a local Python
command uses it. Legacy process/container configurations remain readable.

Inspect an interface with `morse scorers describe cai`, or generate a section
with `morse scorers template cai --command 'YOUR_COMMAND'`. Custom manifests use
`morse scorers describe --manifest /path/to/morse-scorer.yaml`.

Warm Starter uses
Core's `ScorerInputRouter`, so the runtime receives batches already cut to the
configured scope together with protocol-v3 physical-part ranges. `parts`,
`window`, and `full_sequence` scopes all map back to the intersecting design
variables, including the independent start-codon variable. Routed sequences
and their part contexts are cached together, so repeated local contexts are not
evaluated again.

## Validation boundary

Warm Starter produces in-bounds design vectors, sequences and provenance. It
does not duplicate Core's run policy. When Core reads the manifest it:

- rejects a file with the wrong design-space fingerprint or too many candidates;
- verifies optional sequences against vectors;
- applies structural validation and all configured thresholds;
- warns with the rejected manifest IDs and records detailed
  `initial_population` events in `run_log.jsonl.gz`;
- replaces rejected or missing candidates with reproducible random candidates.

This boundary keeps scorer-guided initialization independently releasable while
making Core the single authority for feasibility.
