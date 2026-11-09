# MORSE Translation Efficiency scorer

Orthrus `large-6-track` regressor for ranking mean translation efficiency of
complete synthetic vaccine mRNAs. The deployable output is
`translation_efficiency_score`; larger is better.

## Validation

The selected end-to-end fine-tuning recipe uses Huber loss, a LayerNorm
regression head, backbone LR `1e-4`, and head LR `5e-4`. Nine-fold OOF
validation on 9,324 development transcripts gave Spearman `0.7904`, Pearson
`0.7894`, RMSE `0.4255`, and MAE `0.3296`. The final checkpoint is trained for
four epochs on all 10,348 labelled transcripts. The deployment checkpoint has
been stripped of optimizer state and validated by loading it into the
persistent scorer on an NVIDIA A100.

The input representation reproduces the validated training pipeline: four
sequence channels, a continuous CDS mask, and an all-zero splice channel. The
poly(A) tail should not be sent to this scorer.

## Interface

The scorer implements MORSE's sequence-batch JSONL protocol. Use schema-4 scope
`[five_prime_utr, cds, three_prime_utr]`, excluding poly(A). The scorer uses
protocol-v3 physical-part ranges to build its
continuous coding track from the explicit start codon through the CDS body.
Requests without those ranges are rejected; the model never infers an ORF from
sequence motifs.

Alternative start codons are represented correctly in the coding track, but
non-AUG starts may be outside the training distribution and require independent
validation.

```yaml
- name: translation_efficiency
  command: 'docker exec -i morse-translation-efficiency python /app/server.py'
  timeout_seconds: 1800
  scope: [five_prime_utr, cds, three_prime_utr]
  metrics:
    translation_efficiency_score: {role: optimize, direction: max}
```

This example assumes a prepared container named `morse-translation-efficiency`,
built from `ghcr.io/shpakb/morse-translation-efficiency:2.2.1`. The image includes
the exact ready checkpoint and pinned backbone at its default paths; model
mounts are unnecessary. See the
[scorer guide](../../docs/scorers.md) for runtime setup.

Runtime settings: `ORTHRUS_CHECKPOINT`, `ORTHRUS_MODEL_SOURCE`,
`ORTHRUS_DEVICE`, `ORTHRUS_BATCH_SIZE`, and `ORTHRUS_MAX_TOKENS`.

For Google-hosted Colab and local CPU/GPU environments, the manifest exposes a
portable venv profile. It loads the same deployment checkpoint with a
state-dict-compatible, pure-PyTorch Mamba selective scan and requires no Docker,
CUDA extension, Transformers, or backend server. `ORTHRUS_DEVICE=auto` selects
CUDA when available. This reference path is reproducible but slower than the
optimized server container.

## Ready inference assets

This public package uses ready checkpoints. Model fitting, training datasets,
and cross-validation generators are outside this release. Historical model
metadata is retained under `model_metadata/`.

From the monorepo root:

```bash
export MORSE_ASSETS_DIR="$PWD/.cache/morse/assets"
python assets/manage.py install --asset orthrus-te-full-fit --source /path/to/ready/model_final.pt
python assets/manage.py install --asset orthrus-backbone
python assets/manage.py verify --asset orthrus-te-full-fit
python assets/manage.py paths --asset orthrus-te-full-fit --json
```

The self-contained images include these ready checkpoints. Standalone asset
installation for native inference or rebuilding accepts an authorized local
provider, checks the exact SHA-256, and never substitutes the pretrained backbone
for the fine-tuned checkpoint. See the [scorer guide](../../docs/scorers.md) for the
portable CPU inference runtime and container usage. Published analysis must retain
the model ID, input scope and precision in addition to the scorer version.

The paper also uses `orthrus-te-held-out`, SHA-256
`70913fc29ba9c073100913b56f59542520ec18a473fae358ff19bd57ecddd604`.
Run it with `server_heldout.py` or `morse-scorer-heldout.yaml`.
Its continuous CDS mask includes the terminal stop; the full-fit adapter excludes
that stop from the mask. The checkpoints and adapters are not interchangeable.

## Self-contained GPU images

Both Dockerfiles share the pinned PyTorch 2.2.0/CUDA 12.1 runtime and exact
official Mamba/causal-conv1d release wheels. The build verifies each wheel's
SHA-256 and checks CPython 3.10, PyTorch 2.2.0, CUDA major version 12, and the
expected C++ ABI before installing it. Both native extensions must import
successfully. No CUDA compiler or GPU is required during the build. Downloads
are installed and removed within one layer; their exact URLs and hashes remain
in the retained Dockerfile. The upstream releases compile `sm_70` and `sm_80`
for V100/A100; their CUDA 12.2 wheels are the authors' selected binaries for
PyTorch CUDA 12.x. Actual GPU execution still requires a compatible host driver.

```bash
docker pull ghcr.io/shpakb/morse-translation-efficiency:2.2.1
docker run --rm -i --gpus all ghcr.io/shpakb/morse-translation-efficiency:2.2.1

docker pull ghcr.io/shpakb/morse-translation-efficiency-held-out:2.2.1-heldout
docker run --rm -i --gpus all ghcr.io/shpakb/morse-translation-efficiency-held-out:2.2.1-heldout
```

These JSONL commands read requests from stdin. Full-fit defaults to
`/models/model_final.pt`; held-out defaults to `/models/best.pt` and its distinct
`server_heldout.py` adapter. Both use `/backbone`, also accessible through the
compatibility path `/models/orthrus-backbone`. Checkpoints, backbone configuration,
inference source and pretrained backbone weights are included and hash-checked.
No runtime download or model mount is required.

To rebuild from the monorepo root, supply two narrow named asset contexts:

```bash
docker buildx build --platform linux/amd64 \
  --build-context "checkpoint=$MORSE_ASSETS_DIR/orthrus-te-full-fit" \
  --build-context "backbone=$MORSE_ASSETS_DIR/orthrus-backbone" \
  -t ghcr.io/shpakb/morse-translation-efficiency:2.2.1 scorers/translation-efficiency
docker buildx build --platform linux/amd64 \
  --build-context "checkpoint=$MORSE_ASSETS_DIR/orthrus-te-held-out" \
  --build-context "backbone=$MORSE_ASSETS_DIR/orthrus-backbone" \
  -f scorers/translation-efficiency/Dockerfile.heldout \
  -t ghcr.io/shpakb/morse-translation-efficiency-held-out:2.2.1-heldout scorers/translation-efficiency
```

Explicit copies admit only the matching ready checkpoint and three registered
backbone files, alongside inference code, notices and the build recipe. Model
weights remain outside Git. First-party and upstream MIT notices are retained
under `/usr/share/morse/notices`, and the exact build recipe is in
`/usr/share/morse/sources/Dockerfile`. Installed dependencies retain their own
notices and terms.
