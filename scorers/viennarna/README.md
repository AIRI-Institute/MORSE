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
Core resolves the scope before invoking the runtime. `VIENNARNA_ENABLE_PF=0`
disables the expensive ensemble calculation; `VIENNARNA_MAX_BP_SPAN` limits
the maximum paired-base span, and `0` leaves it unlimited. The default scope is
`full_sequence`; the manifest allows any parts, part-relative slices, and
windows with explicit upstream/downstream context. The runtime always receives
only the sequence already resolved by Core.

Outputs: `mfe` (kcal/mol), `mfe_per_nt`, `gc_fraction`, `mfe_structure`,
`ensemble_free_energy` (kcal/mol), and `ensemble_diversity`.

## Native CPU setup

Install the official `ViennaRNA==2.7.2` Python wheel into the environment used
for the scorer. From the monorepo root:

```sh
uv pip install --python .venv/bin/python -r scorers/viennarna/requirements-portable.txt
```

For an experiment, use a Bash command and explicit metric policies. Replace
`/path/to/morse` with the checkout path:

```yaml
scorers:
  - name: viennarna
    command: 'VIENNARNA_ENABLE_PF=0 VIENNARNA_MAX_BP_SPAN=200 "${MORSE_PYTHON}" /path/to/morse/scorers/viennarna/scorer.py'
    timeout_seconds: 1800
    scope: {type: full_sequence}
    metrics:
      mfe_per_nt: {role: optimize, direction: min}
      mfe: {role: log}
      gc_fraction: {role: log}
      mfe_structure: {role: log, value_type: structure}
```

`MORSE_PYTHON` defaults to MORSE's interpreter. MORSE sends one JSON request and
EOF per batch; the adapter returns one JSON response. Add `scorer.py` to
`provenance_paths` when using caching or resume. Partition-function outputs are
disabled in this example; enable them explicitly for ensemble metrics.

## Reproducibility

The paper uses a pinned container. Native CPU execution uses the same ViennaRNA
API; retain the configured version, span limit and partition-function setting.
Container preparation and optional archiving are separate from optimization:

```bash
# Run from the monorepo root.
docker build --platform linux/amd64 -t ghcr.io/shpakb/morse-viennarna:2.2.0 scorers/viennarna
python scorers/viennarna/prepare_assets.py --docker-image ghcr.io/shpakb/morse-viennarna:2.2.0
```

The image preserves the adapter's MIT notice and ViennaRNA's exact upstream
copyright notice under `/app/licenses/`.

The built image archive and its SHA-256 digest are stored in
`${MORSE_ASSETS_DIR}/viennarna_files` (default `.cache/morse/assets/viennarna_files`).
This scorer requires no external model weights or data files.
