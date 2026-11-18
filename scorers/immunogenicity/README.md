# RiNALMo RNA window classifier

RiNALMo scores short RNA fragments against the immunomodulatory class in
IMrnaDS. MORSE reports `rinalmo_window_logit_p95`: the 95th percentile of the
classification scores (*logits*) of the windows in the selected region.
This score can rank synonymous RNA designs encoding the same protein
by the immunomodulatory signal predicted for their RNA windows.

**Input:** any routed RNA sequence. The manifest defaults to `full_sequence`,
including poly(A). The example below uses the paper's scope, transcript body
without poly(A), and records the result without using it to select candidates.
Use the same scope and window settings when comparing transcripts.

**Hardware:** the released Linux CUDA 13 image targets A100 and compatible
newer GPUs with a compatible NVIDIA driver (R580+) and NVIDIA Container Toolkit;
it does not support V100 on GPU. CPU inference is also available.

## Run with MORSE

On the machine running MORSE, start the model once and leave this terminal open:

```sh
morse scorers serve --port 8813 \
  --command 'docker run --rm -i --gpus all -e RINALMO_DEVICE=cuda ghcr.io/shpakb/morse-immunogenicity:2.1.1'
```

For CPU, replace `--gpus all -e RINALMO_DEVICE=cuda` with `-e RINALMO_DEVICE=cpu`.
In a second terminal, copy the [quickstart YAML](../../examples/quickstart/experiment.yaml)
to `examples/quickstart/rinalmo.yaml` beside its input FASTA. **Append** this list
entry under its existing `scorers:` key, aligned with the other `- name:` entry.
Keep the existing objectives:

```yaml
  - name: rinalmo
    command: 'curl -fsS -H "Content-Type: application/json" --data-binary @- http://127.0.0.1:8813/score'
    timeout_seconds: 1800
    scope: [five_prime_utr, cds, three_prime_utr]
    metrics:
      rinalmo_window_logit_p95: {role: log}
      rinalmo_window_count: {role: log, value_type: int}
```

From the repository root, score the starting construct:

```sh
morse score --config examples/quickstart/rinalmo.yaml
```

Look for `rinalmo.rinalmo_window_logit_p95` and `rinalmo.rinalmo_window_count`
under `outputs`. The service keeps the model loaded until you stop it with Ctrl-C.

## Settings and interpretation

The default scan uses 27-nt windows at a stride of 3 nt. Window count depends on
sequence length and window placement. The scorer also returns mean and maximum
window logits and probability summaries; `morse scorers describe immunogenicity`
lists the metrics. For per-request window settings and batch size,
[use the current adapter with the image](../../docs/scorer-parameters.md#using-this-checkout-with-containers).
Omitting `parameters` uses the prepared runtime's defaults.

## Model and sources

The classifier uses [RiNALMo-mega](https://multimolecule.danling.org/models/rinalmo/)
and 1,122 labelled short sequences from [IMrnaDS](https://doi.org/10.1038/srep20678).
Internal validation on 224 held-out fragments gave ROC-AUC 0.9815. That evaluation
does not establish performance on whole transcripts.
See the [model card](../../docs/model-cards.md#rinalmo-imrnads-window-classifier),
[image build recipe](Dockerfile) and [license and source notices](notices/README.md)
for the checkpoint, tokenizer and MultiMolecule terms.
