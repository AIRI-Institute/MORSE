# MORSE ViennaRNA scorer

Setup and schema-4 integration: [scorer guide](../../docs/scorers.md).

## What it scores

The scorer calculates thermodynamic and compositional RNA properties: minimum
free energy (MFE), MFE per nucleotide, GC fraction, minimum-free-energy
secondary structure, ensemble free energy, and ensemble diversity.

## Method and validation

This is a JSONL adapter around
[ViennaRNA Package 2.7.2](https://github.com/ViennaRNA/ViennaRNA), not a trained
model. It uses `RNA.fold_compound().mfe()` and, unless disabled, the partition
function. Primary method publication:
[Lorenz et al., Algorithms for Molecular Biology, 2011](https://doi.org/10.1186/1748-7188-6-26).

There is no dataset, set of weights, or ML validation metric. Unit tests verify
adapter formatting, normalization, and agreement with the installed library;
physical accuracy is determined by the ViennaRNA model and its parameters.

## Interface

Input: sequence-batch JSONL request `{"request_id": ..., "sequences": [...]}`.
Core resolves the scope before invoking the runtime. Optional request
`parameters` contains calculation settings:

- `max_bp_span`: maximum paired-base span in nucleotides; default 300, or 0 for unrestricted pairing.
- `compute_ensemble`: calculate partition-function metrics; default `false`.

The default scope is `full_sequence`; the manifest allows any parts,
part-relative slices, and windows with explicit upstream/downstream context.
The runtime receives only the sequence already resolved by Core.

Outputs: `mfe` (kcal/mol), `mfe_per_nt`, `gc_fraction`, `mfe_structure`,
`ensemble_free_energy` (kcal/mol), and `ensemble_diversity`.

## Native CPU setup

From the monorepo root, install MORSE with the `mfe` extra. This includes
`ViennaRNA==2.7.2`:

```sh
uv sync --python 3.13 --all-packages --extra paper --extra mfe --frozen
```

For a pip-based virtual environment, install the local packages together:

```sh
python -m pip install ./core ./warm-starter ./scorers/manufacturability '.[mfe]'
```

Use the installed module in the experiment YAML:

```yaml
scorers:
  - name: viennarna
    command: '"${MORSE_PYTHON}" -m scorers.viennarna.scorer'
    parameters:
      max_bp_span: 300
      compute_ensemble: false
    timeout_seconds: 1800
    scope: full_sequence
    metrics:
      mfe: {role: optimize, direction: min}
      mfe_per_nt: {role: log}
      gc_fraction: {role: log}
      mfe_structure: {role: log, value_type: structure}
```

`MORSE_PYTHON` defaults to MORSE's interpreter. MORSE sends one JSON request and
EOF per batch; the adapter returns one JSON response. The command works outside
the checkout as well. `morse-mfe` is the equivalent command in an activated
environment. A direct check is:

```sh
printf '%s\n' '{"request_id":1,"sequences":["AUGGCCGCCGCCGCCGCCUAA"],"parameters":{"max_bp_span":300,"compute_ensemble":false}}' \
  | python -m scorers.viennarna.scorer
```

Set `max_bp_span: 0` for unrestricted pairing. Set `compute_ensemble: true`
when requesting `ensemble_free_energy` or `ensemble_diversity`.

## Reproducibility

The paper uses a pinned container. Native CPU execution uses the same ViennaRNA
API; retain the configured version, span limit and ensemble setting.
The published 2.2.0 image uses the older `VIENNARNA_MAX_BP_SPAN` and
`VIENNARNA_ENABLE_PF` environment variables. Use the updated native adapter or
a new source build for JSON `parameters`; frozen paper configs retain their
original container commands.
Container preparation and optional archiving are separate from optimization:

```bash
# Run from the monorepo root.
docker build --platform linux/amd64 -t morse-viennarna:local scorers/viennarna
python scorers/viennarna/prepare_assets.py --docker-image morse-viennarna:local
```

The image preserves the adapter's MIT notice and ViennaRNA's exact upstream
copyright notice under `/app/licenses/`.

The built image archive and its SHA-256 digest are stored in
`${MORSE_ASSETS_DIR}/viennarna_files` (default `.cache/morse/assets/viennarna_files`).
This scorer requires no external model weights or data files.
