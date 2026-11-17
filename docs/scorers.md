# Choose a scorer

Start with the property you want to measure. Each scorer page explains what the
output means, how to connect it to MORSE, and how to score the starting RNA before
running an optimization. Here **transcript body** means 5′UTR + start + CDS body
+ stop + 3′UTR, without poly(A). **Coding input** means start + CDS body, without
the terminal stop.

| What do you want to measure? | Scorer | Input in the example | How it runs |
|---|---|---|---|
| Codon usage relative to a human or mouse reference | [CAI](../scorers/cai/README.md) | Coding input | Native CPU |
| Codon-associated stability relative to a reference profile | [CSI](../scorers/csi/README.md) | Coding input | Native CPU |
| RNA folding energy | [ViennaRNA](../scorers/viennarna/README.md) | Whole RNA or a chosen region | Native CPU, `mfe` extra |
| Violations of a set of sequence rules | [Sequence feasibility](../scorers/manufacturability/README.md) | Transcript body | Native CPU |
| Predicted translation efficiency | [Orthrus TE, full-fit](../scorers/translation-efficiency/README.md) | Transcript body | GPU image or portable CPU runtime |
| Translation efficiency with the paper’s held-out checkpoint | [Orthrus TE, held-out](../scorers/translation-efficiency/README.md) | Transcript body, AUG required | Separate GPU image |
| Predicted RNA stability | [Orthrus half-life](../scorers/half-life/README.md) | Transcript body | GPU image or portable CPU runtime |
| Translation efficiency averaged across 78 human cell types | [RiboNN](../scorers/ribonn/README.md) | Transcript body, AUG required | GPU or CPU image |
| Predicted translation in a selected cell context | [RiboDecode](../scorers/ribodecode/README.md) | Coding input, up to 4,500 nt | GPU or CPU image |
| Predicted repression by a miRNA panel | [scanMiR](../scorers/irna/README.md) | Whole RNA in the example | CPU image |
| Distribution of IMrnaDS classifier scores across short RNA windows | [RiNALMo](../scorers/immunogenicity/README.md) | Configurable region | GPU or CPU image; see hardware requirements |

## First score

After [installing MORSE](../README.md#install) with the `mfe` extra, try the
[four-scorer CPU example](../examples/scorer-parameters/experiment.yaml) from
the repository root:

```sh
morse score --config examples/scorer-parameters/experiment.yaml
```

The JSON `outputs` contains the scores for **one starting construct**, including
`mfe.mfe` and `cai.cai`. The first part of each key is the name assigned in the
YAML. `morse validate` alone checks the configuration but does not call the scorers.

To use another scorer, follow its page above. For a native scorer, the YAML
calls the installed module via `"${MORSE_PYTHON}"`. For a model used across
many batches, start its service in one terminal and leave it running. The
YAML in its guide calls the same local port in a second terminal. MORSE starts
neither Docker nor a model service during optimization.

## Inspect a scorer

```sh
morse scorers
morse scorers describe cai
morse scorers template cai --command '"${MORSE_PYTHON}" -m scorers.cai.server'
```

`describe` lists outputs, their roles and the available parameters without
loading a model. `template` prints a YAML section for the supplied command.
Use `scope` to choose **what** the scorer reads, `parameters` to change **how**
it computes a score, and `metrics` to optimize, constrain or log an output.
For someone else's scorer, pass `--manifest path/to/morse-scorer.yaml` instead
of a catalog ID. See [custom scorer example](../examples/custom_scorer) for the
request and response format.

## Container images

The published images are for Linux/amd64. The examples in the model guides
use the image's built-in settings; these images include the weights.
For **per-run YAML parameters**, use the native scorer or the
[current adapter with the container](scorer-parameters.md#using-this-checkout-with-containers).
For exact versions, see the image digests in the [runtime lock](../assets/runtime-lock.json).

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

## Runtime choices

For native CAI, CSI, ViennaRNA and sequence feasibility, use the commands on
their scorer pages. Orthrus full-fit TE and half-life also have portable CPU
implementations for environments with PyTorch, NumPy and the matching task
checkpoint. For example, from the repository root in that environment:

```sh
ORTHRUS_PORTABLE=1 ORTHRUS_DEVICE=cpu \
ORTHRUS_CHECKPOINT="$MORSE_ASSETS_DIR/orthrus-te-full-fit/model_final.pt" \
python scorers/translation-efficiency/server.py
```

This process reads JSON requests on stdin. For half-life use its own checkpoint
and `scorers/half-life/server.py`. Held-out TE uses a different checkpoint and
adapter. Record the device and precision when comparing model outputs.

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

For paper reproduction, use the [paper guide](../paper/README.md); the saved runs
specify their own scopes and runtime settings. See the [model cards](model-cards.md)
for checkpoint identity and validation context, and the [third-party notices](../THIRD_PARTY_NOTICES.md)
for the terms covering each model and runtime.
