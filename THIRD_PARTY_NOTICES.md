# Third-party materials and release boundaries

The root MIT license covers first-party MORSE code supplied by 2A2I. It does
not relicense third-party software, model weights, source datasets, biological
sequence libraries, or materials carrying their own notices. The existing
copyright convention is preserved; rights-holder confirmation remains a release
task. Scientific outputs and ready checkpoints do not acquire a new license
merely by being used with MORSE.

Component terms and distribution boundaries are recorded in
`licenses/inventory.json` and `licenses/redistribution-exclusions.json`.
The published neural images retain the notices for their bundled inference
code and ready weights; optional excluded components remain separate.

| Component | License / status | Boundary in this repository |
|---|---|---|
| First-party optimizer, adapters and sequence-feasibility code | MIT | Root [LICENSE](LICENSE) |
| RiboNN upstream code | Sanofi academic/non-commercial notice | Pinned inference source is bundled in the RiboNN image with its original notices; the root MIT license does not cover it |
| RiboNN weights | Separate Sanofi academic/non-commercial weights notice | The selected 50-model ensemble is bundled in the RiboNN image with its weights notice; no weights in Git |
| MultiMolecule / converted RiNALMo implementation | AGPL-3.0-or-later in the model card | Separate optional runtime; original dependency license must accompany images |
| Fine-tuned RiNALMo checkpoint | Separate checkpoint and upstream dependency terms | The exact ready checkpoint is bundled in the RiNALMo image alongside MultiMolecule notices and source archives |
| scanMiR and scanMiRData | GPL-3 | Optional R process/container; retain package notices and corresponding sources when required |
| ViennaRNA | ViennaRNA COPYING | Retain authors' notices and the COPYING conditions in distributions |
| Orthrus upstream backbone | MIT declared by the pinned model's model card | Preserve the upstream notice; pin the exact model/code revision |
| Our Orthrus TE/half-life checkpoints | First-party ready checkpoints; upstream backbone retains its MIT notice | Exact task checkpoints and the pinned backbone are bundled in their matching images; no weights in Git |
| RiboDecode TranslationModel wheel/model | SSPL-1.0 declared in the exact author wheel; copyright 2024 Xie Lab | The author package, checkpoint and conditions are bundled with the exact SSPL notice and package source materials; root MIT does not apply |
| LinearDesign benchmark tool | Author license restricts code redistribution | No vendored checkout/binary; external authorized installation only |
| VaxPress | MIT; optional LinearFold/LinearDesign have separate terms | External pinned benchmark dependency |

The license status of input data is recorded separately from code: GENCODE UTRs,
Kazusa/LinearDesign codon tables, Wu/eLife/iCodon CSI profiles, miRNATissueAtlas,
MirGeneDB and scanMiR panels all require their original provenance/attribution.
Bundled profile JSON and panel source manifests preserve the available source
URLs, revisions and hashes. A public URL alone is not a redistribution license.

Only ready-model inference is included. This release deliberately excludes
model-fitting code and fitting datasets. Bundled checkpoints and their inference
dependencies retain their original notices and applicable source/service terms.
Including a ready checkpoint does not relicense it under the root MIT notice.

## Primary sources

- RiboNN software: <https://github.com/Sanofi-Public/RiboNN/blob/main/LICENSE>
- RiboNN weights: <https://github.com/Sanofi-Public/RiboNN/blob/main/MODEL%20WEIGHTS%20LICENSE.txt>
- Converted RiNALMo: <https://huggingface.co/multimolecule/rinalmo-mega>
- MultiMolecule license FAQ: <https://multimolecule.danling.org/license/faq/>
- scanMiR: <https://bioconductor.org/packages/release/bioc/html/scanMiR.html>
- scanMiRData: <https://bioconductor.org/packages/release/data/experiment/html/scanMiRData.html>
- ViennaRNA: <https://github.com/ViennaRNA/ViennaRNA/blob/master/COPYING>
- Orthrus: <https://huggingface.co/quietflamingo/orthrus-large-6-track/tree/c17ac0db4611c95248001d6fb4bcd95e3293cec4>
- RiboDecode: <https://github.com/wangfanfff/RiboDecode>
- LinearDesign: <https://github.com/LinearDesignSoftware/LinearDesign/blob/main/license.txt>
- VaxPress: <https://github.com/ChangLabSNU/VaxPress/blob/main/LICENSE.txt>

Review date: 2026-11-05. Upstream default branches can change; preserve the
notices shipped with the actual pinned source, wheel, package, model and image.

## TranslationModel package license evidence

The registered `TranslationModel-1.1.0-py3-none-any.whl` (SHA-256
`22d4b398a399878f7573074fd1c381b2f3894bd3a62ede9cd27989dc47ac276a`)
contains `TranslationModel-1.1.0.dist-info/LICENSE.txt`, headed Server Side
Public License, version 1, with copyright 2024 Xie Lab. That exact notice is
preserved byte-for-byte; its SHA-256 and package provenance are recorded in
`licenses/TranslationModel-1.1.0-license.json`. This identifies the actual
package license; it does not clear redistribution, hosted services, or the
source obligations of a proposed distribution. Keep the notice and evaluate
the applicable SSPL terms before publishing a containing image or service.
