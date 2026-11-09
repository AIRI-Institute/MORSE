# MORSE interfering-RNA scorer

Setup and schema-4 integration: [scorer guide](../../docs/scorers.md).

This scorer identifies canonical miRNA binding sites that could reduce antigen
production from a vaccine RNA. It uses `scanMiR 1.16.0` affinity models and
expression-aware human and mouse panels built from `scanMiRData 1.16.0`,
MirGeneDB 3.0, and miRNATissueAtlas 2025.

## Score

`scanMiR` predicts binding-site affinity and aggregates all retained sites for
each miRNA into a predicted repression value. For candidate `m`, miRNA `g`, and
biological profile `p`, this scorer defines:

```text
miRNA risk(m,g)   = max(0, -scanMiR repression(m,g))
profile risk(m,p) = sum_g expression_fraction(g,p) * miRNA risk(m,g)
final risk(m)     = max_p profile risk(m,p)
```

`irna_max_profile_risk` is minimized by MORSE. The maximum across profiles is
deliberately conservative: a sequence should not acquire a strong site for a
highly expressed miRNA in any target immune or injection-site compartment.

The index is useful for ranking otherwise comparable vaccine candidates. It is
not a calibrated absolute fold-change in antigen output because tissue-average
RPMM is not a direct measurement of AGO-loaded miRNA concentration in the
experimental system.

## Fixed scan policy

- canonical sites only;
- `maxLogKd = c(-1.0, -1.5)`;
- minimum same-miRNA site distance: 7 nt;
- match sequences and extended 3-prime alignment details are not retained;
- full routed input is treated as UTR-like sequence by `scanMiR`, so the first
  `maxLogKd` threshold applies throughout the current implementation;
- 3-prime pairing still contributes through the standard `scanMiR`
  aggregation model.

The default scope is the full vaccine sequence. Core may route a configured
part, slice, or window instead.

## Panels

The bundled `2026-08-21` snapshot contains:

| Species | Models | Target profiles | Minimum covered total expression |
|---|---:|---:|---:|
| Human | 157 | PBMC, dendritic cells, macrophage, mononuclear cells, limb muscle | 96.8% |
| Mouse | 170 | PBMC, macrophage, inflammatory monocyte, limb muscle, spleen, lymph node | 96.3% |

The panel TSVs retain source identifiers, mature sequences, matching audit
flags, per-profile RPMM, and normalized expression fractions. Exact source
checksums are recorded in `panels/2026-08-21/sources.json`.

## JSONL interface

```json
{"request_id":1,"sequences":["AAAATAAGCTAAAAAAACGTACGT"]}
```

The process returns one output object per input sequence and remains alive for
subsequent batches. `U` is normalized to `T`; whitespace is ignored; other
bases are rejected.

## Build and verification

```bash
docker build -t morse/irna:1.0.0 .
printf '%s\n' '{"request_id":1,"sequences":["AAAATAAGCTAAAAAAACGTACGT"]}' \
  | docker run --rm -i morse/irna:1.0.0
python prepare_assets.py --root /tmp/morse-irna-assets \
  --docker-image morse/irna:1.0.0
```

The build context is restricted to the two R entry points, the bundled
panels, the Dockerfile, and license notices. The default base image is pinned
by digest; an existing compatible Bioconductor runtime can be selected with
`--build-arg BASE_IMAGE=IMAGE`. The build installs from source and requires
scanMiR and scanMiRData to resolve to exactly 1.16.0.

The image includes the adapter's MIT license, the GPL-3 license for scanMiR
and scanMiRData, the exact source archives used to install both packages,
and dependency and source-checksum manifests in
`/usr/share/licenses/morse-irna`. To retrieve these materials without
starting the scorer:

```sh
container=$(docker create morse/irna:1.0.0)
docker cp "$container":/usr/share/licenses/morse-irna ./morse-irna-licenses
docker rm "$container"
```

See [third-party notices](THIRD_PARTY_NOTICES.md) for package attribution
and panel-source provenance. The complete panel snapshot, including its
source checksums and coverage audit, remains in `/app/panels/2026-08-21`.

## MORSE integration

Select human or mouse models with `IRNA_SPECIES=human|mouse` in the scorer's
startup environment. A native CPU runtime needs R, scanMiR 1.16.0,
scanMiRData 1.16.0, BiocParallel and jsonlite. The bundled panels contain
expression weights; the affinity models come from scanMiRData's `hsa` and `mmu`
collections.

For example, prepare one native process in a separate terminal, replacing
`/path/to/morse` with the checkout path:

```sh
morse scorers serve --command 'IRNA_SPECIES=mouse IRNA_WORKERS=1 Rscript /path/to/morse/scorers/irna/server.R' --port 8803
```

The optional adapter keeps this R process and its models loaded across batches.
For an already prepared container, use
`docker exec -i -e IRNA_SPECIES=mouse -e IRNA_WORKERS=1 CONTAINER Rscript /app/server.R`
as the adapter's command instead. Prepare and start the container before the
experiment; the optimizer only calls the endpoint:

```yaml
scorers:
  - name: irna_mouse
    command: 'curl -fsS -H "Content-Type: application/json" --data-binary @- http://127.0.0.1:8803/score'
    timeout_seconds: 300
    scope: {type: full_sequence}
    metrics:
      irna_max_profile_risk: {role: optimize, direction: min}
      irna_mean_profile_risk: {role: log}
      irna_site_count: {role: log, value_type: int}
      irna_species: {role: log, value_type: string}
```

Stop the adapter with Ctrl-C after the experiment. Record the R sources, selected
panel, package versions and prepared runtime identity in experiment provenance.
