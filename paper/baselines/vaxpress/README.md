# Build the native VaxPress baseline

This recipe installs upstream [VaxPress at commit
dfc36638978e4c8f0925045018fd65678d3dcb3e](https://github.com/ChangLabSNU/VaxPress/tree/dfc36638978e4c8f0925045018fd65678d3dcb3e)
for the paper's shared MFE/nt and CSI comparison. The exact upstream Git archive
and six executed source files are SHA256-verified. All twenty Python/build
packages have fixed versions and wheel hashes. The Python base uses the immutable
digest already recorded in the monorepo runtime lock. The build context is only
this recipe directory; caches, environments, weights and archived runs are not
sent to Docker. From the monorepo root, build for Linux amd64:

```sh
docker build --platform linux/amd64 \
  -t morse/vaxpress:dfc366-shared paper/baselines/vaxpress
docker run --rm morse/vaxpress:dfc366-shared
```

The default container command verifies installed package versions, ViennaRNA
**2.7.2**, source checksums, and the active native `CDSEvolutionChamber` module.
Build-time installation also runs `pip check`. Network access is needed to build;
search execution is offline. No OS package installation or unpinned dependency
resolution is used. The source's MIT license stays in `/opt/VaxPress`.

The recorded historical image used Python 3.10.12; this rebuild uses the pinned
Python 3.12 base. Its image ID will differ from the historical local image ID.
The native optimizer/addon sources, NumPy 2.2.6 and folding version are pinned;
validate numerical/search replay using the generated evidence rather than
assuming image-byte or trajectory identity. The original [setup requirements](https://github.com/ChangLabSNU/VaxPress/blob/dfc36638978e4c8f0925045018fd65678d3dcb3e/setup.py)
and [native search](https://github.com/ChangLabSNU/VaxPress/blob/dfc36638978e4c8f0925045018fd65678d3dcb3e/vaxpress/evolution_chamber.py)
remain unchanged. This experiment supplies its sole shared-fitness addon directly;
iCodon/R, LinearFold and LinearDesign are not invoked or included.

## Check exactly 64 scoring events

From the full monorepo, run one bounded native search at seed 7 and alpha 0:

```sh
python paper/baselines/vaxpress/smoke.py \
  --image morse/vaxpress:dfc366-shared \
  --output outputs/vaxpress-native64 --execute
```

The helper checks the frozen CSI/ViennaRNA source and codon-profile hashes,
resolves the actual local image ID, and invokes the **same**
`paper/adapters/run_vaxpress_shared.py` used by the full baseline. It verifies 64
recorded events, the 470-nt transcript body and 117 scored codons. The last native
batch has 20 proposals; only its first 4 are scored, leaving 16 unscored proposals.
It records command/config, stdout/stderr, native evaluation and batch logs,
provenance, summary and `execution.json`. Existing evidence directories are never
replaced. Omit `--execute` to prepare a plan in a separate new directory.

## Use the rebuilt image for paper runs

Create an image map pointing at the local build:

```json
{"vaxpress": "morse/vaxpress:dfc366-shared"}
```

Save it as `outputs/baseline-image-map.json`, then run:

```sh
morse paper baseline vaxpress --seed 7 --runtime container \
  --image-map outputs/baseline-image-map.json \
  --output outputs/vaxpress-replay --execute
```

That frozen sweep uses five weights and 2,560 scoring events per weight, totaling
12,800 events for the seed. For all five independent baseline seeds, execute the
five seeds listed by `morse paper list`, each with a new output root, or omit
`--seed` to run all five native sweeps. `--scoring-events 320` is an explicitly
labelled diagnostic sweep of 64 events per weight.

The private local source/recipe suite is retained in the maintainer's checkout
and excluded from the public Git snapshot. In that local checkout it runs
without Docker or baseline dependencies:

```sh
python -m unittest discover -s paper/baselines/vaxpress/tests -v
```
