# MORSE Half-Life scorer

Orthrus `large-6-track` regressor for ranking the expected stability of complete
synthetic vaccine mRNAs. The deployable score is `half_life_score`; larger is
better. It predicts the normalized target used by the RNA half-life dataset,
not hours of physical half-life.

## Validation

The selected end-to-end fine-tuning recipe uses Huber loss, a LayerNorm
regression head, backbone LR `3e-5`, head LR `3e-4`, and 75% splice-track
dropout. Four-fold homology-grouped OOF validation on 10,998 development
transcripts gave Spearman `0.6350`, Pearson `0.6421`, RMSE `0.7799`, and MAE
`0.6079`. The final checkpoint is trained for five epochs on all 12,939 labelled
transcripts. The deployment checkpoint has been stripped of optimizer state and
validated by loading it into the persistent scorer on an NVIDIA A100.

Synthetic IVT input uses four sequence channels, a codon-start CDS channel, and
an all-zero splice channel. The poly(A) tail should not be sent to this scorer.

## Interface

The scorer implements MORSE's sequence-batch JSONL protocol. Use schema-4 scope
`[five_prime_utr, cds, three_prime_utr]`, excluding poly(A). The scorer uses
protocol-v3 physical-part ranges to build its
codon-start coding track from the explicit start codon through the CDS body.
Requests without those ranges are rejected; the model never infers an ORF from
sequence motifs.

Alternative start codons are represented correctly in the coding track, but
non-AUG starts may be outside the training distribution and require independent
validation.

```yaml
- name: half_life
  command: 'docker exec -i morse-half-life python /app/server.py'
  parameters:
    batch_size: 12
    max_tokens: 36000
  timeout_seconds: 1800
  scope: [five_prime_utr, cds, three_prime_utr]
  metrics:
    half_life_score: {role: optimize, direction: max}
```

This example assumes a prepared container named `morse-half-life`, built from
`ghcr.io/shpakb/morse-half-life:2.2.1`. Its exact ready checkpoint and pinned
backbone are included at their default paths; model mounts are unnecessary.
See the [scorer guide](../../docs/scorers.md) for runtime setup.

Set `batch_size` and `max_tokens` in YAML `parameters`. The token budget limits
padded tokens per batch without truncating sequences; `null` removes the budget.
Missing values retain the prepared runtime defaults. Existing `ORTHRUS_*`
environment settings remain supported. Device and model paths remain runtime
setup options. See the [parameter reference](../../docs/scorer-parameters.md).

The YAML uses the current source adapter; see [container setup](../../docs/scorer-parameters.md#using-this-checkout-with-containers)
when preparing an existing registry image.

For Google-hosted Colab and local CPU/GPU environments, the manifest exposes a
portable venv profile. It loads the same deployment checkpoint with a
state-dict-compatible, pure-PyTorch Mamba selective scan and requires no Docker,
CUDA extension, Transformers, or backend server. `ORTHRUS_DEVICE=auto` selects
CUDA when available. This reference path is reproducible but slower than the
optimized server container. Portable defaults are `batch_size: 2` and
`max_tokens: null`.

## Ready inference assets

This public package uses ready checkpoints. Model fitting, training datasets,
and cross-validation generators are outside this release. Historical model
metadata is retained under `model_metadata/`.

From the monorepo root:

```bash
export MORSE_ASSETS_DIR="$PWD/.cache/morse/assets"
python assets/manage.py install --asset orthrus-half-life --source /path/to/ready/model_final.pt
python assets/manage.py install --asset orthrus-backbone
python assets/manage.py verify --asset orthrus-half-life
python assets/manage.py paths --asset orthrus-half-life --json
```

The self-contained image includes the ready checkpoint. Standalone installation
for native inference or rebuilding accepts an authorized local provider, checks
the exact SHA-256, and never substitutes the pretrained backbone for the fine-tuned
checkpoint. See the [scorer guide](../../docs/scorers.md) for the
portable CPU inference runtime and container usage. Published analysis must retain
the model ID, input scope and precision in addition to the scorer version.

## Self-contained GPU image

The Dockerfile shares the pinned PyTorch 2.2.0/CUDA 12.1 runtime and exact
official Mamba/causal-conv1d release wheels with both Orthrus TE images. The
build verifies the wheel SHA-256 digests and checks CPython 3.10, PyTorch 2.2.0,
CUDA major version 12, and the expected C++ ABI. Both native extensions must
import successfully. No compiler or GPU is required during the build. Downloads
are installed and removed within one layer; the retained Dockerfile records
their exact URLs and hashes. Upstream releases compile `sm_70` and `sm_80` for
V100/A100 and select these CUDA 12.2 wheels for PyTorch CUDA 12.x. Actual GPU
execution still requires a compatible host driver.

```bash
docker pull ghcr.io/shpakb/morse-half-life:2.2.1
docker run --rm -i --gpus all ghcr.io/shpakb/morse-half-life:2.2.1
```

This JSONL command reads requests from stdin. The image defaults to
`/models/model_final.pt` and `/backbone`; `/models/orthrus-backbone` is a
compatibility alias. The matching ready checkpoint, backbone configuration,
inference source and pretrained backbone weights are included and hash-checked.
No runtime download or model mount is required.

To rebuild from the monorepo root:

```bash
docker buildx build --platform linux/amd64 \
  --build-context "checkpoint=$MORSE_ASSETS_DIR/orthrus-half-life" \
  --build-context "backbone=$MORSE_ASSETS_DIR/orthrus-backbone" \
  -t ghcr.io/shpakb/morse-half-life:2.2.1 scorers/half-life
```

The narrow named asset contexts admit only the matching ready checkpoint and
three registered backbone files, alongside inference code, notices and the build
recipe. Model weights remain outside Git. First-party and upstream MIT notices
are retained under `/usr/share/morse/notices`, and the exact Dockerfile under
`/usr/share/morse/sources/Dockerfile`. Installed dependencies retain their own
notices and terms.
