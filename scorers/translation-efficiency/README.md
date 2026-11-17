# Orthrus translation efficiency

`translation_efficiency_score` ranks transcripts by predicted mean translation
efficiency. Higher is better. The score does not predict an absolute protein yield.

Both variants take **5′UTR + start + CDS body + stop + 3′UTR**, without poly(A).
MORSE supplies the CDS coordinates. Choose one variant for a comparison:

| Variant | When to use it | Coding track |
|---|---|---|
| Full-fit (`translation-efficiency`) | New experiments; fitted to all labelled training transcripts | Marks start and CDS body, excluding stop. Accepts alternative starts; predictions for non-AUG starts are less well supported. |
| Held-out (`translation-efficiency-held-out`) | Comparisons that specify the paper's held-out model, including Fig. 5 | Marks the stop as well. Requires AUG and a canonical in-frame terminal stop. |

The two models return the same metric name on different scales. Do not compare
a full-fit score for one sequence with a held-out score for another.

## Run with MORSE

On a Linux GPU host (e.g. A100) with Docker, NVIDIA Container Toolkit and MORSE
installed, start **one** of the following services. Leave its terminal open.
The model stays loaded between batches.

Full-fit:

```sh
morse scorers serve --port 8811 \
  --command 'docker run --rm -i --gpus all ghcr.io/shpakb/morse-translation-efficiency:2.2.1'
```

Held-out, **instead** of full-fit:

```sh
morse scorers serve --port 8811 \
  --command 'docker run --rm -i --gpus all ghcr.io/shpakb/morse-translation-efficiency-held-out:2.2.1-heldout'
```

In another terminal on the same host, copy the [quickstart YAML](../../examples/quickstart/experiment.yaml)
to `examples/quickstart/te.yaml` beside its input FASTA. Replace its `scorers:`
section with:

```yaml
scorers:
  - name: orthrus_te_full_fit
    command: 'curl -fsS -H "Content-Type: application/json" --data-binary @- http://127.0.0.1:8811/score'
    timeout_seconds: 1800
    scope: [five_prime_utr, cds, three_prime_utr]
    metrics:
      translation_efficiency_score: {role: optimize, direction: max}
```

For the held-out service use `name: orthrus_te_held_out` so the results identify
the model. From the repository root, score the reference once:

```sh
morse score --config examples/quickstart/te.yaml
```

Look for `orthrus_te_full_fit.translation_efficiency_score` (or its held-out
counterpart) under `outputs`. Stop the service with Ctrl-C after the experiment.

## Settings

The example uses the published image's built-in settings. To set per-request
YAML `parameters`, [start the image with the current scorer adapter](../../docs/scorer-parameters.md#using-this-checkout-with-containers)
and select the adapter matching the checkpoint. `batch_size` limits sequences
per batch; `max_tokens` limits padded tokens per batch, without truncating
sequences. The current GPU defaults are 8 / 32000 for full-fit and 4 / null
(no token limit) for held-out. Omit these keys to inherit runtime settings.

For CPU inference with the full-fit checkpoint, see the
[portable PyTorch runtime](../../docs/scorers.md#runtime-choices).
The held-out model requires its own adapter and environment.

## Model and sources

Nine-fold out-of-fold validation of the full-fit recipe on 9,324 development
transcripts gave Spearman 0.7904 (Pearson 0.7894). The final model was fitted
on all 10,348 labelled transcripts. See the [validation record](model_metadata/reference_run.json),
[model cards](../../docs/model-cards.md) and [license notices](notices/README.md).
The [full-fit](Dockerfile) and [held-out](Dockerfile.heldout) build recipes pin
their respective checkpoints and the Orthrus backbone.
