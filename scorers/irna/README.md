# scanMiR: predicted miRNA repression

`irna_max_profile_risk` is the largest predicted repression index across the tissue and cell profiles in the selected panel. Lower values indicate less predicted repression. Use the index to compare sequences scored with the same scope and panel; it is not a calibrated fold change in protein output.

**Input:** by default, the entire RNA including poly(A). scanMiR treats every routed base as UTR-like, even within the CDS. **Runtime:** Linux/amd64 CPU container with R, scanMiR, scanMiRData and both species panels included. Human is the default.

## Run with MORSE

In an activated MORSE environment, leave the local scoring service running in one terminal:

```sh
morse scorers serve --port 8803 \
  --command 'docker run --platform linux/amd64 --rm -i ghcr.io/shpakb/morse-irna:1.0.0'
```

For **mouse** instead, add `-e IRNA_SPECIES=mouse` before the image name when starting the service. This selects the mouse panel and affinity models in the published image. In a second terminal on the same host, use this scorer section for either species:

```yaml
scorers:
  - name: irna
    command: 'curl -fsS -H "Content-Type: application/json" --data-binary @- http://127.0.0.1:8803/score'
    scope: full_sequence
    timeout_seconds: 1800
    metrics:
      irna_max_profile_risk: {role: log}
      irna_species: {role: log, value_type: string}
```

For a first check from the repository root, copy [the quickstart YAML](../../examples/quickstart/experiment.yaml) to `examples/quickstart/irna.yaml` beside its input FASTA. Append the list item above to its existing `scorers:` list. This records the scanMiR results without changing the example's objectives. Then run:

```sh
morse score --config examples/quickstart/irna.yaml
```

Check `irna.irna_species` and `irna.irna_max_profile_risk` under `outputs`. The service keeps the models loaded between requests. Stop it with Ctrl-C when finished.

## Choose a panel

| Panel | miRNA models | Tissue and cell profiles |
|---|---:|---|
| [Human](panels/2026-08-21/human/core_975.tsv) | 157 | PBMC, dendritic cells, macrophage, mononuclear cells, limb muscle |
| [Mouse](panels/2026-08-21/mouse/core_975.tsv) | 170 | PBMC, macrophage, inflammatory monocyte, limb muscle, spleen, lymph node |

To select a panel in YAML instead of at service startup, use the [current adapter](../../docs/scorer-parameters.md#using-this-checkout-with-containers). Add this under the scorer:

```yaml
parameters:
  profile_path: panels/2026-08-21/mouse/core_975.tsv
```

The bundled path also selects the species. Set `profiles` to a subset of the panel's profile names, or omit it to use all. The defaults scan canonical sites, use `max_log_kd: -1.0` and a minimum same-miRNA site spacing of 7 nt. See `morse scorers describe irna` for all settings and outputs.

## How to read the score

scanMiR aggregates retained binding sites into a repression estimate per miRNA. Within each profile, the scorer weights each miRNA's nonnegative repression magnitude by its expression fraction and sums the contributions. The main score is the **maximum** of these profile sums. `irna_mean_profile_risk` reports their mean. Tissue expression is used as a proxy for the available miRNA pool.

The frozen [panel snapshot](panels/2026-08-21/sources.json) records sources from scanMiRData 1.16.0, MirGeneDB 3.0 and miRNATissueAtlas 2025. The runtime uses scanMiR 1.16.0. See [third-party notices](THIRD_PARTY_NOTICES.md) for attribution. The published image includes the MIT adapter notice, GPL-3 notices and the scanMiR/scanMiRData source archives in `/usr/share/licenses/morse-irna`.

For a native R installation with those packages, BiocParallel and jsonlite, the scorer command is `Rscript /absolute/path/to/scorers/irna/server.R`. To build a CPU image from the current checkout instead, run this from the repository root:

```sh
docker build --platform linux/amd64 -t morse-irna:local scorers/irna
```
