# MORSE immunogenicity scorer

Setup and schema-4 integration: [scorer guide](../../docs/scorers.md).

## What it scores

The scorer estimates the probability that a short RNA fragment belongs to the
immunomodulatory IMrnaDS class. In MORSE, a full construct is divided into
overlapping windows and `rinalmo_window_logit_p95` is minimized. This is a
MORSE application-level aggregation, not RiNALMo's original task.

## Method and data provenance

We fine-tuned [RiNALMo-mega](https://multimolecule.danling.org/models/rinalmo/)
with mean pooling and a binary head on 1,122
[IMrnaDS](https://webs.iiitd.edu.in/raghava/imrna/download.php) sequences: 602
positive and 520 negative examples approximately 17–27 nt long. The original
dataset publication is
[Chaudhary et al., Scientific Reports, 2016](https://doi.org/10.1038/srep20678).
Dataset preparation is outside this inference-only release. The shared asset
manager installs and verifies the ready checkpoint.

The reference run used a stratified random 80/20 split (`seed=42`) with 898
training and 224 validation examples. Best epoch: ROC-AUC 0.9815, average
precision 0.9879, accuracy 0.9464, balanced accuracy 0.9468, and F1 0.9496.
This is internal validation on short oligonucleotides; applying windowed scoring
to a full vaccine construct has not yet been experimentally validated.

## Interface

Input: sequence-batch JSONL request `{"request_id": ..., "sequences": [...]}`.
Core resolves the scope before invoking the runtime. `RINALMO_CHECKPOINT`,
`RINALMO_DEVICE`, `RINALMO_WINDOW_SIZE` (27), `RINALMO_WINDOW_STRIDE` (3),
`RINALMO_WINDOW_OFFSETS`, and `RINALMO_BATCH_SIZE` are fixed once at startup.
The default scope is `full_sequence`; the manifest allows any parts,
part-relative slices, and windows with explicit upstream/downstream context.
The runtime always receives only the sequence already resolved by Core.

Outputs: `rinalmo_window_logit_mean`, `rinalmo_window_logit_median`,
`rinalmo_window_logit_max`, `rinalmo_window_logit_p95`,
`rinalmo_window_probability_mean`, `rinalmo_window_probability_max`,
`rinalmo_window_probability_p95`, and `rinalmo_window_count`.

For a prepared container named `morse-immunogenicity`, add this schema-4 section:

```yaml
scorers:
  - name: immunogenicity
    command: 'docker exec -i morse-immunogenicity python /app/server.py'
    timeout_seconds: 1800
    scope: full_sequence
    metrics:
      rinalmo_window_logit_p95: {role: optimize, direction: min}
```

## Ready inference assets

The published image includes the ready checkpoint, configuration, and tokenizer;
it requires no external model files or mounts. To rebuild it from an authorized
local checkpoint, prepare the verified build asset from the monorepo root:

```bash
python assets/manage.py install --asset rinalmo-imrnads --source /path/to/ready/best_model.pt
python assets/manage.py verify --asset rinalmo-imrnads
```

Only ready-model inference is included. `runtime_model.py` reconstructs the
checkpoint-compatible network from the bundled configuration and loads every
parameter strictly; tokenizer files are bundled. No fitting code, fitting data,
losses, split builders, or evaluation loops are included. The ready checkpoint
has no registered public asset-manager download; source builds use an authorized
local copy. MultiMolecule and its converted backbone retain their AGPL terms; see
`THIRD_PARTY_NOTICES.md` before redistributing code, images, or derived weights.

## Public runtime image

The pinned PyTorch 2.12.1/CUDA 13.0 image is intended for A100 and compatible
newer GPUs with a CUDA-13-compatible host driver (NVIDIA branch R580 or newer)
and NVIDIA Container Toolkit. This image does not support V100/Volta: CUDA 13
removed support for pre-Turing GPUs. See the
[CUDA 13 release notes](https://docs.nvidia.com/cuda/archive/13.0.0/cuda-toolkit-release-notes/index.html)
and [driver compatibility table](https://docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html).
CPU inference remains available with `RINALMO_DEVICE=cpu`.

The self-contained `ghcr.io/shpakb/morse-immunogenicity:2.1.1` image includes
`/models/best_model.pt` with SHA-256
`923d1e9f843d12fe29a9515a95b6a32eae286017790dfeef1c82215446958f0a`.
Both the build and runtime verify this digest. `RINALMO_CHECKPOINT` and
`MORSE_EXPECTED_CHECKPOINT_SHA256` select this bundled checkpoint by default.
Run the published image directly:

```bash
docker run --rm -i --gpus all -e RINALMO_DEVICE=cuda \
  ghcr.io/shpakb/morse-immunogenicity:2.1.1
```

For a source build, the ordinary context allowlists the adapter,
`runtime_model.py`, configuration/tokenizer files, and notices. A separate
named context supplies only the prepared `rinalmo-imrnads` asset directory;
the Dockerfile copies its exact `best_model.pt`. From the monorepo root:

```bash
export MORSE_ASSETS_DIR="$PWD/.cache/morse/assets"
docker build --platform linux/amd64 \
  --build-context assets="$MORSE_ASSETS_DIR/rinalmo-imrnads" \
  -t ghcr.io/shpakb/morse-immunogenicity:2.1.1 scorers/immunogenicity
```

No model-fitting data enters either build context or the runtime image.

The image retains MultiMolecule 0.2.0's original `license.md` and
`license.header` from the hash-verified source archive, plus its matching wheel
metadata and notice provenance, under `/usr/share/morse/notices/multimolecule`.
Its unmodified wheel and source
archive, their SHA-256 digests, all installed Python package versions, and the
Dockerfile under `/usr/share/morse/sources`. The adapter source stays in `/app`.
Preserve these corresponding-source materials and the applicable AGPL terms in
redistributed copies. The image build reads the ready checkpoint from the named
asset context and records its identity in
`/usr/share/morse/sources/checkpoint-manifest.json`. Configuration and tokenizer
files contain no trained parameters. Dependency notices and the first-party MIT
notice remain separate.
