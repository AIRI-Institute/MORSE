# Ready-model inference cards

These cards describe the frozen inference artifacts used by MORSE, their input
contracts and paper roles. No model-fitting recipes or fitting datasets are
included. Exact file digests and ready-provider installation instructions are
in `assets/registry.json` and `docs/scorers.md`. The cards do not replace the
upstream models' licenses or establish experimental performance of a designed
construct.

## Orthrus TE, full-fit

- Asset: `orthrus-te-full-fit`, `model_final.pt`.
- SHA-256: `c380316fe282640d563b9254e9ab8031e129f9e768c9ce04c0c7955bb0ac6bab`.
- Input: complete transcript body (5′UTR, explicit start, CDS body, stop, 3′UTR),
  excluding poly(A). Four A/C/G/T channels, a continuous CDS mask **excluding
  terminal stop**, and a zero splice channel. Explicit part coordinates required.
- Output: `translation_efficiency_score`; inverse-standardized normalized mean
  TE score. Higher is favoured. It is not absolute protein yield.
- Runtime: `scorers/translation-efficiency/server.py`; optimized CUDA/Mamba or
  the separate pure-PyTorch portable implementation. CUDA inference uses BF16
  autocast; record the chosen implementation, precision and device.
- Paper: objectives in Fig. 3a–c, Fig. 4b and Extended Data Fig. 2; full-fit TE
  axes in all Fig. 4 radars; additional recorded scores of selected candidates.
- Availability: verified local ready-provider installation. Public distribution
  of our checkpoint is not yet registered. The Orthrus backbone is pinned to
  `c17ac0db4611c95248001d6fb4bcd95e3293cec4`.

## Orthrus TE, common held-out

- Asset: `orthrus-te-held-out`, `best.pt`.
- SHA-256: `70913fc29ba9c073100913b56f59542520ec18a473fae358ff19bd57ecddd604`.
- Input: same transcript-body scope, but its continuous CDS mask **includes
  terminal stop**. Four A/C/G/T channels and zero splice channel. The frozen
  adapter requires AUG and a canonical in-frame stop, with complete part ranges.
- Output: `translation_efficiency_score`, in its checkpoint's normalized score
  scale. Higher is favoured. Its scores are not interchangeable with full-fit TE.
- Runtime: `server_heldout.py` and `heldout_encoding.py`; frozen checkpoint task
  and SHA-256 are checked before loading. Optimized CUDA uses BF16 autocast.
- Paper: objective in Fig. 5; recorded in Fig. 3d and additional selected-candidate
  readouts; association analysis in Fig. 2.
- Availability: local ready provider; uses the same separately pinned backbone.

## Orthrus half-life

- Asset: `orthrus-half-life`, `model_final.pt`.
- SHA-256: `fe4f816572626d127f54e64d4a724923bd5a080b4898f225898f69f91f28a40e`.
- Input: transcript body without poly(A), A/C/G/T channels, codon-start CDS
  channel, zero splice channel, explicit physical-part coordinates. The current
  adapter excludes terminal stop from its coding interval.
- Output: `half_life_score`, inverse-standardized RNA-stability target score.
  Higher is favoured. The output is not a physical half-life in hours.
- Runtime: optimized CUDA/Mamba or portable PyTorch. Record precision/device.
- Paper: recorded in Fig. 3d and selected-candidate comparisons; association
  analysis in Fig. 2. It did not guide those searches.
- Availability: local ready provider; checkpoint distribution terms pending.

## RiboNN human ensemble

- Asset: `ribonn-human-ensemble`; official Zenodo record `17258709`.
- Archive SHA-256: `2ce3bc1adb38ec6752c79147241a4fb2664e9804e806b3eab1287489a15d8b64`.
- Trusted ensemble manifest SHA-256:
  `8e0d62ca61964ec660a43389ab309361f0c388a906c68c48a9d5b928c5a10cee`.
- Input: transcript body without poly(A), A/T/C/G channel order, codon-start
  channel including stop. The upstream alignment allows at most 1,381 bases
  before CDS start and 11,937 from CDS start onward; out-of-range inputs fail.
- Output: `ribonn_mean_te`, mean over five selected models per fold, ten folds
  and 78 human cellular outputs; `ribonn_model_sd` is dispersion among the 50
  model mean scores, not a calibrated uncertainty interval.
- Runtime: FP32, TF32 and stochastic shifts disabled. Each model digest is
  checked. The inference-only extraction retains checkpoint parameter names
  and numerical operations; see the preserved modification notice and
  `assets/ribonn-static-parity.json` for the limited static audit.
- Paper: objective in Fig. 5; recorded Fig. 3d and selected-candidate readouts;
  association analysis in Fig. 2.
- License: separate Sanofi academic/non-commercial software and weights terms.
  A publicly downloadable archive does not remove those restrictions.

## RiboDecode TranslationModel

- Asset: `ribodecode-translationmodel`, author wheel version `1.1.0`.
- Wheel SHA-256: `22d4b398a399878f7573074fd1c381b2f3894bd3a62ede9cd27989dc47ac276a`.
- Input: start codon + synonymous CDS body, **no terminal stop or UTRs**;
  maximum 4,500 nt. Paper context is HEK293T with fixed RNA abundance 4.5.
- Outputs: `ribodecode_translation_log` is the direct log-scale model output;
  `ribodecode_translation = expm1(z)/5` is the inverse-transformed proxy.
  Neither is a measured protein yield.
- Runtime: upstream inference implementation in its isolated Python 3.8/PyTorch
  environment. Cellular context and abundance must remain recorded.
- Paper: objectives Fig. 3a–c, Fig. 4c, Fig. 5 and Extended Data Fig. 3;
  association analysis Fig. 2; other selected-candidate readouts.
- License: the registered author wheel declares **SSPL-1.0**, copyright 2024
  Xie Lab. The exact notice and its SHA-256 are preserved under `licenses/`.
  This package is excluded from the root MIT license. Source, hosted-service
  and redistribution obligations require review; no clearance is claimed.
- Availability: author public wheel download with pinned SHA-256. The package
  notice is `licenses/TranslationModel-1.1.0-SSPL-1.0.txt`, with provenance in
  `licenses/TranslationModel-1.1.0-license.json`.

## RiNALMo IMrnaDS window classifier

- Asset: `rinalmo-imrnads`, `best_model.pt`.
- SHA-256: `923d1e9f843d12fe29a9515a95b6a32eae286017790dfeef1c82215446958f0a`.
- Input: routed RNA sequence, 27-nt windows at stride 3, including the terminal
  window. Paper input is transcript body without poly(A). Tokenizer and model
  configuration are bundled; no online model lookup is needed for loading.
- Outputs: window logit/probability summaries and count. The paper uses the
  95th percentile of positive-class window logits (`rinalmo_window_logit_p95`),
  with lower values favoured for display.
- Runtime: `runtime_model.py` contains only the ready-checkpoint network and
  masked nucleotide mean pooling; `server.py` performs windowed inference.
- Paper: recorded Fig. 3d and association analysis Fig. 2; did not guide search.
- Limits: short-oligonucleotide classification applied across longer constructs;
  it is not a direct measurement of whole-mRNA innate sensing or immunogenicity.
- License: MultiMolecule/converted backbone AGPL terms and derived-checkpoint
  source/distribution obligations require review. Ready provider is local.

## Profiles and sequence policies

CAI/CSI are empirical profile calculations, ViennaRNA is thermodynamic folding,
scanMiR is a separately versioned affinity/repression model with expression
panels, and sequence feasibility is a configurable engineering policy. Their
profile/source hashes and explicit scopes remain part of each experiment.
