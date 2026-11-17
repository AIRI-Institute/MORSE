# Orthrus RNA half-life

`half_life_score` ranks transcripts by predicted RNA stability. A higher value
means greater predicted stability on the model's normalized target scale.
A difference of one score unit is **not** an extra hour of half-life.

**Input:** 5′UTR + start + CDS body + stop + 3′UTR, without poly(A). MORSE
passes explicit CDS coordinates. The model marks the first base of each coding
codon except the terminal stop. It accepts alternative start codons, but
predictions for non-AUG starts are less well supported.

## Run with MORSE

On a Linux GPU host (e.g. A100) with Docker, NVIDIA Container Toolkit and MORSE
installed, start a service in a separate terminal. Leave it running to keep the
model loaded between batches:

```sh
morse scorers serve --port 8812 \
  --command 'docker run --rm -i --gpus all ghcr.io/shpakb/morse-half-life:2.2.1'
```

On the same host, copy the [quickstart YAML](../../examples/quickstart/experiment.yaml)
to `examples/quickstart/half-life.yaml` beside its input FASTA. Replace the
`scorers:` section with:

```yaml
scorers:
  - name: orthrus_half_life
    command: 'curl -fsS -H "Content-Type: application/json" --data-binary @- http://127.0.0.1:8812/score'
    timeout_seconds: 1800
    scope: [five_prime_utr, cds, three_prime_utr]
    metrics:
      half_life_score: {role: optimize, direction: max}
```

From the repository root, score the starting construct:

```sh
morse score --config examples/quickstart/half-life.yaml
```

Read `orthrus_half_life.half_life_score` under `outputs`. To use it only as a
readout alongside other objectives, set its metric to `{role: log}` as in the
paper's comparisons. Stop the service with Ctrl-C when finished.

## Settings

The example uses the published image as supplied. To pass YAML `parameters`,
[start it with the current scorer adapter](../../docs/scorer-parameters.md#using-this-checkout-with-containers).
`batch_size` sets the number of sequences per batch (GPU default 12);
`max_tokens` limits padded tokens per batch (GPU default 36000) without
truncating sequences. Omit these keys to inherit runtime settings.

For CPU inference with the same checkpoint, see the
[portable PyTorch runtime](../../docs/scorers.md#runtime-choices). Its defaults
are 2 sequences per batch and no token limit. Record the runtime and precision
when comparing scores from different machines.

## Model and sources

Four-fold homology-grouped out-of-fold validation on 10,998 development
transcripts gave Spearman 0.6350 (Pearson 0.6421). The final model was fitted on
12,939 labelled transcripts. See the [validation record](model_metadata/reference_run.json),
[model card](../../docs/model-cards.md#orthrus-half-life), [build recipe](Dockerfile)
and [license notices](notices/README.md).
