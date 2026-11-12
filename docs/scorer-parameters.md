# Scorer parameters

Put calculation settings in `parameters`, beside `command`, `scope` and `metrics`:

```yaml
scorers:
  - name: mfe
    command: '"${MORSE_PYTHON}" -m scorers.viennarna.scorer'
    parameters:
      max_bp_span: 300  # nt; 0 removes the limit.
      temperature: 37.0  # Celsius.
      compute_ensemble: false
    scope: full_sequence
    metrics:
      mfe: {role: optimize, direction: min}
```

Only supplied values override the prepared scorer's settings. Leaving out
`parameters`, using `{}`, or omitting a particular key retains the previous defaults,
including settings supplied through environment variables.
Unknown names, invalid types and unsupported values are rejected.

Inspect settings or generate a YAML section without loading the model:

```sh
morse scorers describe cai
morse scorers template cai --command '"${MORSE_PYTHON}" -m scorers.cai.server'
```

Generated templates contain explicit defaults. Remove a setting to inherit it
from your prepared runtime instead. Adapters report effective settings in the
response's `parameters` field;
optimization runs save them to `scorer_parameters.json`.

## Available settings

These are defaults without environment overrides. All batch sizes, window sizes
and strides are positive integers. Boolean settings use YAML `true` or `false`.

| Scorer | Settings and defaults |
|---|---|
| ViennaRNA MFE | `max_bp_span: 300`, `compute_ensemble: false`, `temperature: 37.0`, `dangles: 2`, `no_lonely_pairs: false` |
| CAI | `profile_path`: bundled human profile; `low_weight_threshold: 0.1` |
| CSI | `profile_path`: bundled human profile; `window_codons: 15`, `destabilizing_threshold: 0.0` |
| Orthrus TE full-fit | `batch_size: 8`, `max_tokens: 32000` |
| Orthrus TE held-out | `batch_size: 4`, `max_tokens: null` |
| Orthrus half-life | `batch_size: 12`, `max_tokens: 36000` |
| RiboNN | `batch_size: 16`, `cpu_threads: 2` |
| RiboDecode | `cell_context: HEK293T`, `context_file: null`, `batch_size: 256` |
| RiNALMo | `window_size: 27`, `window_stride: 3`, `window_offsets: []`, `batch_size: 256` |
| scanMiR | `profile_path`: bundled human panel; `species: human`, `profiles: null`, `max_log_kd: -1.0`, `min_site_distance: 7`, `canonical_only: true`, `workers: 4` |
| Sequence feasibility | `preset: sequence_feasibility_v1` and individual checks below |

Orthrus portable runtimes retain `batch_size: 2` and no token budget by default.
`max_tokens` limits padded tokens per model batch; it does not truncate sequences.
Use `null` to disable that budget. A single longer sequence is processed on its
own. Runtime device, checkpoints and model architecture remain part of scorer
setup.

CAI's `low_weight_threshold` and CSI's `window_codons` and
`destabilizing_threshold` affect diagnostic outputs, not their main score.
RiboDecode supports `HEK293T`, `A549` and `HeLa`; a `context_file` selects a custom
CSV with the official 10,552-gene order and one RPKM column. RiNALMo uses 27-nt
windows by default; `window_offsets` is a YAML list of nonnegative offsets.
scanMiR's `profiles: null` uses all tissue/cell profiles in the selected panel;
a list selects a subset.

For MFE, `dangles` accepts `0`, `1`, `2` or `3`; with `compute_ensemble: true`,
use `0` or `2`. Ensemble metrics are returned only when that calculation is enabled.

## Human and mouse profiles

Both files are included under each scorer's directory. Human is the default.
Select a file through `profile_path`:

| Scorer | Human | Mouse |
|---|---|---|
| CAI | `profiles/homo_sapiens_kazusa_93487_cds.json` | `profiles/mus_musculus_kazusa_53036_cds.json` |
| CSI | `profiles/homo_sapiens_consensus_wu2019_v1.json` | `profiles/mus_musculus_mesc_slamseq_wu2019_v1.json` |
| scanMiR | `panels/2026-08-21/human/core_975.tsv` | `panels/2026-08-21/mouse/core_975.tsv` |

These bundled paths resolve inside the scorer directory, including `/app` in a
container. Custom file paths must be accessible to the scorer process. For
Docker, mount a custom file and use its container path. The mouse CSI profile is
the bundled mES-cell profile. scanMiR infers species from the bundled panel path;
there is no need to repeat `species` when selecting one of these files.

When caching or resuming experiments with custom profiles or contexts, list the
host source files in `provenance_paths` so changes to their contents are recorded.

```yaml
parameters:
  profile_path: profiles/mus_musculus_kazusa_53036_cds.json
```

## Sequence feasibility

Missing settings inherit the selected preset. The default preset uses:

```yaml
parameters:
  preset: sequence_feasibility_v1
  # enabled_checks: [global_gc, local_gc, homopolymers, direct_repeats, tandem_repeats, forbidden_motifs]
  global_gc_min: 0.4
  global_gc_max: 0.6
  local_gc_min: 0.3
  local_gc_max: 0.7
  local_gc_window: 50
  local_gc_stride: 1
  homopolymer_max_a: 9
  homopolymer_max_u: 9
  homopolymer_max_g: 5
  homopolymer_max_c: 5
  direct_repeat_length: 21
  direct_repeat_min_copies: 2
  tandem_min_period: 2
  tandem_max_period: 6
  tandem_min_length: 18
  tandem_min_copies: 3
  forbidden_motifs: []
  exempt_regions: []
```

Lengths are in nucleotides; GC bounds are fractions. `exempt_regions` contains
objects with `start`, `end`, `checks` and `reason`, using zero-based, end-excluded
coordinates of the sequence passed to this scorer.

`manufacturing_strategy_v1` additionally enables restriction-site, capping-start
and inverted-repeat checks. Its settings include `restriction_enzymes: [BsaI,
BsmBI]`, `capping_method: cleancap_ag`, `custom_5prime_prefix: ""`,
`inverted_repeat_stem: 20`, `inverted_repeat_window: 200` and
`inverted_repeat_min_gap: 0`. These checks remain disabled in the ordinary preset
unless explicitly enabled.

## Using this checkout with containers

The request parameters above are implemented in this checkout. Existing registry
images retain their bundled adapter version. To use the updated adapter with an
existing image, mount the matching scorer source directory when preparing the
container, for example on a CUDA host:

```sh
docker run --rm -i --gpus all \
  --mount type=bind,src="$PWD/scorers/translation-efficiency",dst=/app,readonly \
  ghcr.io/shpakb/morse-translation-efficiency:2.2.1
```

Run this command from the repository root. The image still supplies its ready
weights and dependencies. Once your runtime is prepared, use its ordinary Bash
command in the run configuration.
