# RiboDecode: translation in a cell context

RiboDecode predicts a coding sequence's translation level for a selected cellular RNA-seq context. Use `ribodecode_translation_log` as the objective: higher values indicate higher predicted translation. `ribodecode_translation` expresses the same prediction on the inverse-transformed scale. Neither is a measurement of protein yield or translation efficiency.

**Input:** start codon + CDS body, without the terminal stop or UTRs; complete codons only, up to 4,500 nt. Non-AUG starts can be scored but may be outside the training distribution. The default cellular context is **HEK293T**.

## Run with MORSE

The Linux/amd64 image includes the model and its built-in contexts. On a Docker host with an NVIDIA GPU, leave this running in a terminal with the MORSE environment active:

```sh
morse scorers serve --port 8807 \
  --command 'docker run --rm -i --gpus all ghcr.io/shpakb/morse-ribodecode:2.1.1'
```

The service reuses the loaded model for each batch. In a second terminal on the same host, use this scorer section:

```yaml
scorers:
  - name: ribodecode
    command: 'curl -fsS -H "Content-Type: application/json" --data-binary @- http://127.0.0.1:8807/score'
    scope: [start_codon, cds_body]
    timeout_seconds: 1800
    metrics:
      ribodecode_translation_log: {role: optimize, direction: max}
      ribodecode_translation: {role: log}
```

For a first check from the repository root, copy [the quickstart YAML](../../examples/quickstart/experiment.yaml) to `examples/quickstart/ribodecode.yaml` beside its input FASTA. Replace its `scorers:` section with the one above, then run:

```sh
morse score --config examples/quickstart/ribodecode.yaml
```

Look for `ribodecode.ribodecode_translation_log` under `outputs`. Stop the service with Ctrl-C when finished.

## Choose a context

To use a different built-in context with the published image, add, for example, `-e RIBODECODE_ENV=A549` before the image name. Supported contexts are `HEK293T`, `A549` and `HeLa`.

For per-run YAML settings, use the [current adapter](../../docs/scorer-parameters.md#using-this-checkout-with-containers). Its `parameters.cell_context` selects the context; `batch_size` defaults to 256. To use a custom context, mount the CSV read-only into the container and set `parameters.context_file` to its **container path**. It must contain the 10,552 gene IDs in the template order and one column of nonnegative RPKM values. The template is at `/assets/env_file.csv` in the image. See `morse scorers describe ribodecode` for the full interface.

## Method, license and source builds

The upstream [TranslationModel](https://github.com/wangfanfff/RiboDecode) was trained on `ln(Ribo-seq RPKM × 5 + 1)`. The adapter fixes RNA abundance at 4.5. If the model returns `z`, the second output is `(exp(z) − 1) / 5`. See [Li et al.](https://doi.org/10.1038/s41467-025-64894-x) for training datasets and validation.

TranslationModel 1.1.0 retains its SSPL-1.0 license and Xie Lab copyright. Its notice and corresponding sources are included in the image under `/usr/share/morse/notices` and `/usr/share/morse/sources/TranslationModel`. See the [original license](../../licenses/TranslationModel-1.1.0-SSPL-1.0.txt) for redistribution and hosted-service terms. The MORSE adapter has its own MIT notice.

To rebuild from this directory, supply the registered wheel and context CSV listed in [the asset manifest](notices/image-assets.json):

```sh
docker build --platform linux/amd64 \
  --build-context "assets=$MORSE_ASSETS_DIR/container-image-inputs/ribodecode-2.1.1" \
  -t morse-ribodecode:local .
```
