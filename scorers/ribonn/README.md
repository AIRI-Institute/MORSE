# RiboNN: human translation efficiency

`ribonn_mean_te` estimates normalized translation efficiency, averaged over 78 human cell types. Higher values indicate higher predicted TE. `ribonn_model_sd` measures disagreement among the 50 models; it is not a confidence interval for protein yield.

**Input:** the transcript from the 5′UTR through the 3′UTR, excluding poly(A). RiboNN requires an AUG start and a canonical in-frame stop. The 5′UTR can be at most 1,381 nt; the start through the end of the 3′UTR can be at most 11,937 nt. Longer inputs are rejected. MORSE supplies the part coordinates.

## Run with MORSE

The published Linux/amd64 image includes the human ensemble. On a Docker host with an NVIDIA GPU, run this from an activated MORSE environment and leave the terminal open:

```sh
morse scorers serve --port 8806 \
  --command 'docker run --rm -i --gpus all ghcr.io/shpakb/morse-ribonn:1.0.1'
```

On CPU, omit `--gpus all`. The service keeps the ensemble loaded between batches. In a second terminal on the same host, add this section to your experiment YAML:

```yaml
scorers:
  - name: ribonn
    command: 'curl -fsS -H "Content-Type: application/json" --data-binary @- http://127.0.0.1:8806/score'
    scope: [five_prime_utr, cds, three_prime_utr]
    timeout_seconds: 1800
    metrics:
      ribonn_mean_te: {role: optimize, direction: max}
      ribonn_model_sd: {role: log}
```

For a first check from the repository root, copy [the quickstart YAML](../../examples/quickstart/experiment.yaml) to `examples/quickstart/ribonn.yaml` (beside its input FASTA), replace its `scorers:` section with the one above, and run:

```sh
morse score --config examples/quickstart/ribonn.yaml
```

Look for `ribonn.ribonn_mean_te` under `outputs`. Stop the service with Ctrl-C when finished.

## Settings and interpretation

The example uses the image defaults. For a run-specific `parameters.batch_size` (default 16) or `parameters.cpu_threads` (default 2), prepare the image with the [current adapter](../../docs/scorer-parameters.md#using-this-checkout-with-containers). These settings control computation, not the target cell types. To inspect all outputs, run `morse scorers describe ribonn`.

The ensemble averages the five highest-validation-R² models from each of ten outer folds. Predictions use FP32. RiboNN and MORSE's Orthrus TE model use labels from the same atlas, so agreement is not independent evidence. The input does not encode chemical modifications.

[Paper](https://doi.org/10.1038/s41587-025-02712-x) · [Upstream source at the bundled revision](https://github.com/Sanofi-Public/RiboNN/tree/dc896e1d1f78bf6b4e9bb7525c7852fcf25532e8) · [Weights](https://zenodo.org/records/17258709)

## License and source builds

The upstream software and weights have separate Sanofi academic/non-commercial terms; see the [software license](third_party/ribonn/LICENSE) and [weights license](third_party/ribonn/MODEL%20WEIGHTS%20LICENSE.txt). Both notices accompany the image. The MORSE adapter has its own MIT notice.

To rebuild from this directory, prepare the files listed in [the asset manifest](notices/image-assets.json) and use the named build context:

```sh
docker build --platform linux/amd64 \
  --build-context "assets=$MORSE_ASSETS_DIR/container-image-inputs/ribonn-1.0.1" \
  -t morse-ribonn:local .
```
