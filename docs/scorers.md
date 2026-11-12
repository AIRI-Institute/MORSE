# Scorers and ready inference assets

The paper uses CAI, CSI, ViennaRNA, sequence feasibility, RiboDecode, RiboNN,
Orthrus TE, Orthrus half-life, scanMiR and RiNALMo. Orthrus TE has two distinct
ready checkpoints and adapters. Additional predictors can be connected through
custom commands and manifests.
Scorers communicate through JSONL and remain isolated processes/environments.
Heavy dependencies are not installed with the lightweight optimizer.

## Container images

The runtime images use `linux/amd64`. Pin the immutable registry reference from
[`assets/runtime-lock.json`](../assets/runtime-lock.json) when recording an
experiment; its `published_images` entries contain the image and config digests.
The six neural images include their exact ready checkpoints and required inference
assets. They run without a host model cache or asset mounts. The five CPU image
versions are unchanged.

| Interface | Image | Runtime assets |
|---|---|---|
| CAI | `ghcr.io/shpakb/morse-cai:2.2.0` | Optional custom profile |
| CSI | `ghcr.io/shpakb/morse-csi:1.0.0` | Optional custom profile |
| ViennaRNA | `ghcr.io/shpakb/morse-viennarna:2.2.0` | None |
| Sequence feasibility | `ghcr.io/shpakb/morse-manufacturability:4.1.0` | None; frozen policies are included |
| scanMiR | `ghcr.io/shpakb/morse-irna:1.0.0` | None; frozen panels are included |
| RiboNN | `ghcr.io/shpakb/morse-ribonn:1.0.1` | Bundled ensemble and inference source |
| RiboDecode | `ghcr.io/shpakb/morse-ribodecode:2.1.1` | Bundled author package, checkpoint, conditions and CSV |
| Orthrus TE, full-fit | `ghcr.io/shpakb/morse-translation-efficiency:2.2.1` | Bundled full-fit checkpoint and backbone |
| Orthrus TE, held-out | `ghcr.io/shpakb/morse-translation-efficiency-held-out:2.2.1-heldout` | Bundled held-out checkpoint and backbone |
| Orthrus half-life | `ghcr.io/shpakb/morse-half-life:2.2.1` | Bundled half-life checkpoint and backbone |
| RiNALMo | `ghcr.io/shpakb/morse-immunogenicity:2.1.1` | Bundled task checkpoint, configuration and tokenizer |

For example, `docker pull ghcr.io/shpakb/morse-cai:2.2.0` prepares the CAI image.
Pull the selected neural image and start it directly; its default model paths
refer to bundled files. A compatible GPU host and container runtime are still
needed for CUDA inference.
The container names in research commands, such as `morse-cai`, are chosen when
starting containers; they are independent of the image references in this table.
CAI and CSI retain their dependency-free native Python defaults in the catalog.

## Inspect and connect an interface

```sh
morse scorers
morse scorers describe cai
morse scorers describe cai --json
morse scorers template cai --command '"${MORSE_PYTHON}" -m scorers.cai.server'
morse scorers describe --manifest /path/to/morse-scorer.yaml
morse scorers template --manifest /path/to/morse-scorer.yaml --command './my-client.sh'
```

These commands read interface metadata without loading a model. Prepare the
scorer process or container first, then place its Bash `command`, `scope`, optional `parameters` and `metrics`
in the run YAML. Each command receives one JSON request and EOF, and returns one
JSON response with matching `request_id` and ordered `outputs`. A local script,
`docker exec -i` or an HTTP client with `curl --data-binary @-` can implement this
contract. MORSE does not prepare or manage external services and containers.
See `docs/researcher-guide.md` for complete examples.

The optional `parameters` mapping holds calculation settings.
Values can be JSON scalars, lists or objects; the adapter validates its supported
names and ranges. Missing settings retain the prepared scorer's defaults,
including existing environment settings. See the [parameter reference](scorer-parameters.md) for all scorers,
human/mouse profile files and using updated adapters with existing images.

## Native Python scorers

CAI, CSI and sequence feasibility are installed with MORSE. For native MFE,
include the `mfe` extra when setting up the environment:

```sh
uv sync --python 3.13 --all-packages --extra paper --extra mfe --frozen
```

The installed modules work from any directory:

| Scorer | JSONL command |
|---|---|
| CAI | `python -m scorers.cai.server` or `morse-cai` |
| CSI | `python -m scorers.csi.server` or `morse-csi` |
| MFE | `python -m scorers.viennarna.scorer` or `morse-mfe` |
| Sequence feasibility | `python -m manufacturability_scorer.feasibility_server` |

In YAML, use `"${MORSE_PYTHON}"` in place of `python` to select MORSE's environment.
CAI and CSI include human and mouse profiles in their `profiles/` directories.
Select one through `parameters.profile_path`; human is the default. scanMiR bundles
both species' panels under `panels/2026-08-21/` and also accepts `profile_path`.
See the [CPU example](../examples/scorer-parameters/experiment.yaml).
When caching or resuming, list the adapter and selected profile files in
`provenance_paths`.
GPU scorer runtimes keep their separate dependencies and containers.
For native MFE, configure `max_bp_span` and `compute_ensemble` in `parameters`;
see the linked YAML example. The published ViennaRNA 2.2.0 container uses the
older environment-variable interface; use the updated native adapter for JSON parameters.

## Optional native assets and source builds

Published containers do not require these steps. For native processes, explicit
custom asset roots, or rebuilding an image from source, prepare ready files from
the monorepo root:

```bash
export MORSE_ASSETS_DIR="$PWD/.cache/morse/assets"
python assets/manage.py list
python assets/manage.py install --asset orthrus-backbone
python assets/manage.py install --asset orthrus-te-full-fit --source /path/to/ready/model_final.pt
python assets/manage.py install --asset orthrus-te-held-out --source /path/to/ready/best.pt
python assets/manage.py verify --asset orthrus-te-held-out
python assets/manage.py paths --asset orthrus-te-held-out --json
```

Ready checkpoint assets have an exact SHA-256. The installer validates all
registered files before publishing an asset directory, refuses an incorrect
checkpoint, and verifies an existing installation instead of overwriting it.
The published neural images provide the registered checkpoints as part of their
immutable container contents. The file-based installer uses an authorized local
checkpoint/provider for the task checkpoints; it does not extract them from an
image or substitute an unregistered model.
`paths` reports paths without loading weights; use `verify` to check contents.
Set `MORSE_EXPECTED_CHECKPOINT_SHA256` to the registered digest when invoking
an adapter directly. Paper manifests set it explicitly.

| Registry ID | File / purpose | Availability |
|---|---|---|
| `orthrus-te-full-fit` | `model_final.pt`; continuous CDS mask excludes stop | Local ready provider for source/native setup; bundled in image |
| `orthrus-te-held-out` | `best.pt`; continuous CDS mask includes stop | Local ready provider for source/native setup; bundled in image |
| `orthrus-half-life` | `model_final.pt`; codon-start CDS channel | Local ready provider for source/native setup; bundled in image |
| `orthrus-backbone` | Pinned HF config, model and patched inference code | Immutable public revision |
| `rinalmo-imrnads` | `best_model.pt`; bundled tokenizer/config | Local ready provider for source/native setup; bundled in image |
| `ribonn-human-ensemble` | 50 selected human models, official archive | Zenodo; Sanofi academic/non-commercial terms |
| `ribodecode-translationmodel` | Author wheel 1.1.0 and context CSV | Public author download; SSPL-1.0, obligations require review |
| `cai-profiles`, `csi-profiles` | Frozen empirical profiles | Bundled |
| `sequence-feasibility-profiles`, `scanmir-panels` | Frozen policies/panels | Bundled |

Only ready-model inference is included. Model-fitting code, fitting datasets,
split generation, loss functions and fitting/evaluation loops are excluded.
Historical metadata describing a checkpoint is retained as provenance.

## Runtime choices

CAI, CSI and sequence feasibility are lightweight Python processes. ViennaRNA
uses `ViennaRNA==2.7.2` on CPU. scanMiR uses a separate R/Bioconductor container.
RiboNN, RiboDecode and learned models use their own pinned dependency sets.

Full-fit Orthrus TE/half-life can run using the existing pure-PyTorch portable
runtime, selected by `ORTHRUS_PORTABLE=1`. It needs only the ready checkpoint
and PyTorch/NumPy. It is slower than the optimized CUDA/Mamba runtime; compare
scores within declared numerical tolerances and keep precision/device provenance.
The held-out TE adapter uses the pinned HF backbone and its matching environment.

```bash
# Use the corresponding scorer environment's Python.
ORTHRUS_PORTABLE=1 \
ORTHRUS_CHECKPOINT="$MORSE_ASSETS_DIR/orthrus-te-full-fit/model_final.pt" \
ORTHRUS_DEVICE=cpu python scorers/translation-efficiency/server.py
```

The published neural containers include everything needed by their inference
adapters. Orthrus reads its task checkpoint at `/models` and the pinned backbone
at `/backbone`; RiboNN includes its exact inference source and 50-model ensemble;
RiboDecode includes TranslationModel 1.1.0, its checkpoint, condition arrays and
environment CSV; RiNALMo includes its task checkpoint, configuration and tokenizer.
The builds verify registered asset digests and retain licenses and source records.

```bash
docker pull ghcr.io/shpakb/morse-translation-efficiency:2.2.1
docker run --rm -i --gpus all ghcr.io/shpakb/morse-translation-efficiency:2.2.1
docker run --rm -i --gpus all ghcr.io/shpakb/morse-ribonn:1.0.1
docker run --rm -i --gpus all ghcr.io/shpakb/morse-ribodecode:2.1.1
```

For held-out TE, use
`ghcr.io/shpakb/morse-translation-efficiency-held-out:2.2.1-heldout`. Its
checkpoint and preprocessing are distinct from the full-fit image. Source
manifests are `morse-scorer.yaml` and `morse-scorer-heldout.yaml` in the TE
directory. Optional manifest path parameters allow explicit read-only model
or backbone overrides; they have no required default mounts. Keep the matching
checkpoint digest and preprocessing when comparing with frozen paper results.
The scorer READMEs document each named-context source build and bundled asset.

The asset manager remains useful for native or custom setups. RiboDecode asset
installation verifies the wheel and CSV, then prepares its exact pure-Python
package without importing it or installing dependencies. RiboNN asset installation
also prepares the pinned inference extraction and original Sanofi notices.
`verify --runtime` checks these prepared trees; `paths` reports optional host
mount locations. These preparations are unnecessary for the published images.

For auxiliary paper replay, select immutable images with
`--image-map assets/runtime-lock.json`. A `self_contained: true` catalog entry
uses its bundled assets when no explicit `--asset-root` is supplied. Optional
`--source-overlay` mounts only the current adapter files for these images;
their bundled packages, weights and backbone remain in place. Supplying
`--asset-root /path/to/assets` explicitly selects host assets and their
verification. Omitting the image map retains historical runtime choices;
legacy maps and native processes keep their existing asset setup. Frozen paper
YAMLs and model batch/device presets remain unchanged; see
[the paper instructions](../paper/README.md).

RiboNN source and weights retain their Sanofi academic/non-commercial notices.
TranslationModel retains SSPL-1.0 and its source/service terms. MultiMolecule
retains its AGPL terms and corresponding source materials. These notices are
included with the relevant images; the root MIT notice does not replace them.
See [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).

## Independent experiments

Use your FASTA and UTR libraries with scorer instances declaring explicit scopes
and metric roles (`optimize`, `threshold`, `log`). CAI/CSI accept custom profile
JSON; miRNA scorers accept explicit profiles; sequence-feasibility rules remain
declared policies. For an additional predictor, implement the JSON contract and
provide its Bash command. An optional manifest enables interface inspection
and template generation. Do not combine unrelated
model outputs into an implicit biological quality score.

Record checkpoint/profile digest, adapter source identity, sequence scope,
precision, device, container digest, resolved configuration and seed with each
experiment. This information is required to compare an independent result
with the frozen paper release.

## Model cards and runtime identities

See `docs/model-cards.md` for input contracts, score scales, ready-checkpoint
identity, paper usage and limitations. `assets/runtime-lock.json` records
Dockerfile content hashes, resolved base-image digests where available and
unavailable upstream tags. Historical local image IDs are recorded separately;
they are not public pull URLs. Transitive apt/R/Python layers still require an
archived runtime or stronger dependency locks for byte-identical rebuilds.
