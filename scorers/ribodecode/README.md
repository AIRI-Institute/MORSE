# MORSE RiboDecode scorer

Setup and MORSE integration: [scorer guide](../../docs/scorers.md).

## What it scores

The scorer predicts coding-sequence translation level in a selected cellular
context. The
model was trained on `ln(Ribo-seq RPKM × 5 + 1)`; this is a computational
translation proxy, not absolute protein yield or translation efficiency.

## Method, weights, and validation

This is a thin JSONL adapter around the official
[`TranslationModel 1.1.0`](https://drive.google.com/file/d/19Yyl8uUbjQn0TAA4E52ZUGBd_pb_RyW8/view)
from the [RiboDecode repository](https://github.com/wangfanfff/RiboDecode). We
did not retrain the model. The architecture receives a one-hot CDS, adapter-fixed
abundance `4.5`, and an RNA-seq vector describing the cellular environment.
Publication and data description:
[Li et al., Nature Communications, 2025](https://doi.org/10.1038/s41467-025-64894-x).

The authors trained the model on 320 paired Ribo-seq/RNA-seq datasets from 24
tissues and cell lines, with more than 10,000 mRNAs per dataset. They report R²
values of 0.81, 0.89, and 0.81 for unseen genes, unseen environments, and both
factors together. These metrics belong to the upstream model; we have not
performed a separate biological validation of this adapter.

## Interface

Core sends only a batch of already extracted coding sequences, each assembled
as `start_codon + cds_body` without a terminal stop. Each coding sequence must
be at most 4,500 nt long and divisible by three. YAML `parameters.cell_context`
selects `HEK293T` (default), `A549`, or `HeLa`. Set `context_file` to use a custom
CSV: exactly 10,552 Gene IDs in the official template order and one column of
nonnegative RPKM values. `context_file: null` selects the built-in context.
`batch_size` defaults to `256`. Missing values retain the prepared runtime
settings, including existing `RIBODECODE_*` overrides. See the [parameter reference](../../docs/scorer-parameters.md).

Alternative start codons can be routed, but non-AUG starts may be outside the
training distribution and require independent validation.

Outputs: `ribodecode_translation_log` is the direct model output;
`ribodecode_translation` is the inverse transformation `(exp(log) - 1) / 5`.

For a prepared container named `morse-ribodecode` with the current source adapter
([container setup](../../docs/scorer-parameters.md#using-this-checkout-with-containers)), use:

```yaml
scorers:
  - name: ribodecode
    command: 'docker exec -i morse-ribodecode python3.8 /app/server.py'
    parameters:
      cell_context: HEK293T
      context_file: null
      batch_size: 256
    timeout_seconds: 1800
    scope: [start_codon, cds_body]
    metrics:
      ribodecode_translation_log: {role: optimize, direction: max}
      ribodecode_translation: {role: log}
```

## Reproducibility

The image includes the MORSE adapter, Python/CUDA dependencies, exact author
TranslationModel 1.1.0 package, ready checkpoint, built-in condition arrays, and
official environment CSV. It starts without asset mounts or model downloads.

```bash
docker run --rm -i --gpus all ghcr.io/shpakb/morse-ribodecode:2.1.1
```

The installed package and its original notices are under `/opt/ribodecode`,
which is already on `PYTHONPATH`. The adapter defaults to its bundled
`model_config.json`, `best_model.p`, and cellular-condition NPZs.
`RIBODECODE_ENV_TEMPLATE` defaults to `/assets/env_file.csv`. Explicit
`RIBODECODE_BEST_MODEL`, `RIBODECODE_MODEL_CONFIG`, or `RIBODECODE_ENV_CSV`
settings continue to support user-supplied alternatives.

To rebuild from this scorer directory, supply a dedicated named context with
only the registered author wheel and `env_file.csv`. Their exact hashes and
every installed source/weight member hash are in `notices/image-assets.json`.

```bash
docker build --platform linux/amd64 \
  --build-context "assets=$MORSE_ASSETS_DIR/container-image-inputs/ribodecode-2.1.1" \
  -t ghcr.io/shpakb/morse-ribodecode:2.1.1 .
```

The wheel is installed directly from a read-only build mount. Its temporary
distribution archive is not duplicated in a separate image layer. The ready
checkpoint is at `/opt/ribodecode/TranslationModel/best_model.p`. No training
pipeline or fitting dataset is added to the monorepo or image.

## Upstream license and source

TranslationModel 1.1.0 retains SSPL-1.0 and copyright 2024 Xie Lab. The full
original notice remains in its installed dist-info and under
`/usr/share/morse/notices/TranslationModel-SSPL-1.0.txt`; the MORSE adapter has its
separate MIT notice. All author-supplied Python/configuration/condition and
packaging sources are accessible under `/usr/share/morse/sources/TranslationModel`,
with asset locations/hashes and the exact installation recipe under
`/usr/share/morse/sources`. Adapter source is `/app/server.py`.

Retain those source materials and notices when redistributing the covered
program. SSPL section 13 additionally requires Service Source Code when offering
the program's functionality to third parties as a hosted service; the license
defines the scope of that service source. Installed dependencies retain their
own terms.
