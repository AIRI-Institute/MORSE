# RiboNN scorer for MORSE

Setup and MORSE integration: [scorer guide](../../docs/scorers.md).

JSONL adapter for the official **human RiboNN ensemble**. The maximized output
`ribonn_mean_te` is the equal mean of predictions over 78 human cell types,
averaging the five highest-validation-R2 models within each of ten outer folds
(50 models total). No model is selected using the publication's reporter data.
All weights are frozen; no downstream fine-tuning is performed.

Source: <https://github.com/Sanofi-Public/RiboNN>, revision
`dc896e1d1f78bf6b4e9bb7525c7852fcf25532e8`.
Weights: <https://zenodo.org/records/17258709>, `weights.zip`, MD5
`dd088ee2feffcda6ce26124a471a140f`.
Paper: <https://doi.org/10.1038/s41587-025-02712-x>.

## Inputs and inference

Full transcript body: 5′UTR + AUG + CDS body + terminal stop + 3′UTR, excluding
poly(A). Explicit `sequence_contexts[].part_ranges` are required; the CDS is
never guessed by searching for the first AUG. Sequences are aligned at the CDS
start, with 1,381 positions to its left and 11,937 to its right. Nucleotide
channels follow upstream A/T/C/G order. A fifth channel marks the first base of
every CDS codon, including the stop codon. Out-of-range sequences and ambiguous
bases fail explicitly; there is no truncation.

The bundled upstream model is loaded without architecture changes, in FP32
evaluation mode. TF32 and stochastic shifts are disabled. No chemical-modification
input exists: U/m¹Ψ/mo⁵U variants of the same sequence receive the same prediction.
`ribonn_model_sd` is model dispersion, not a calibrated confidence interval.

## Run

```bash
docker run --rm -i --gpus device=0 ghcr.io/shpakb/morse-ribonn:1.0.1
```

The image contains the MORSE adapter, checkpoint-compatible upstream inference
source, and all 50 selected model weights. It starts with the bundled ensemble
at `/models` and needs no asset mounts or startup downloads. The source is at
`/app/third_party/ribonn`. An explicit `RIBONN_MODEL_DIR` override can select
another compatible ensemble; the manifest's host model-directory override is
optional.

Requests: `{"request_id":"...","sequences":["..."],"sequence_contexts":[{"part_ranges":[...]}]}`.
Responses: `{"request_id":"...","outputs":[{"ribonn_mean_te":0.1,"ribonn_model_sd":0.02,"ribonn_model_count":50}]}`.
Set `batch_size` (default `16`) and `cpu_threads` (default `2`) in YAML
`parameters`. Missing values retain the prepared runtime settings. `RIBONN_MODEL_DIR`
and `RIBONN_DEVICE` remain setup options; existing batch/thread environment
overrides are supported. See the [parameter reference](../../docs/scorer-parameters.md).

For a prepared container named `morse-ribonn` with the current source adapter
([container setup](../../docs/scorer-parameters.md#using-this-checkout-with-containers)), use:

```yaml
scorers:
  - name: ribonn
    command: 'docker exec -i morse-ribonn python /app/server.py'
    parameters:
      batch_size: 16
      cpu_threads: 2
    timeout_seconds: 1800
    scope: [five_prime_utr, cds, three_prime_utr]
    metrics:
      ribonn_mean_te: {role: optimize, direction: max}
      ribonn_model_sd: {role: log}
```

The target is normalized translation efficiency, not absolute protein yield.
Our Orthrus TE model uses labels from the same atlas: these are different
architectures, not independent training-label sources.

## Upstream license

Sanofi's software and weights notices cover academic research and non-commercial
use by academic or non-profit organizations. The software notice also describes
commercial-benefit and private-company-collaboration restrictions. Both original
notices are included at `/app/third_party/ribonn/LICENSE` and
`/app/third_party/ribonn/MODEL WEIGHTS LICENSE.txt`, together with the inference
extraction's modification and source-revision records. The MORSE adapter retains
its separate MIT notice under `/usr/share/morse/notices`. Installed dependencies
retain their own terms.

## Rebuilding the image

From this scorer directory, use a dedicated named asset context containing only
the 52 files listed in `notices/image-assets.json`: the selected 50 weights,
`ensemble.json`, and `human/runs.csv`, preserving relative paths. The build checks
every file's SHA-256 and the exact inference-source/notice hashes. The complete
upstream archive and unused ensemble members are unnecessary build inputs.

```bash
docker build --platform linux/amd64 \
  --build-context "assets=$MORSE_ASSETS_DIR/container-image-inputs/ribonn-1.0.1" \
  -t ghcr.io/shpakb/morse-ribonn:1.0.1 .
```

The ordinary source context remains small and allowlisted. The exact asset
manifest and build recipe are retained under `/usr/share/morse/notices` and
`/usr/share/morse/sources`. This image contains ready inference assets and no
model-fitting pipeline or training dataset.
