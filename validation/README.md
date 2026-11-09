# Validation of the local assembly

Validation date: 2026-11-05. Real inference and optimization checks ran in an
isolated checkout on olig (NVIDIA A100). The CPU environment was also exercised
locally on macOS/Python 3.13.3. No results were published.

The 2026-11-08 repository cleanup is recorded in
[repository-cleanup.json](repository-cleanup.json): 361 CPU tests passed and four
optional IRNA container tests were skipped. All 47 configs and 14 README YAML
examples validated. Five CPU examples, the notebook's full local workflow and
an isolated installation from public files passed. Existing initial-generation
logs were rechecked and all eight paper figures regenerated from frozen data;
paper searches and neural inference were not rerun.

The 2026-11-08 schema-4 migration is recorded in [schema4.json](schema4.json):
47 YAML files validated, 363 CPU tests passed and four IRNA container tests
skipped. The commented template and working vaccine config completed local CPU
searches. Paper inputs and scientific settings were checked offline across all
35 configs; no paper searches or model inference were rerun. Colab uses backend
API 6 with YAML and sequence inputs exported together.

The 2026-11-06 prepared-command and Colab simplification is recorded in
`command-yaml-colab.json`: 35 ordinary paper YAML configs, 15 scorer interfaces,
323 current CPU tests (318 passed, 5 optional-runtime skips), and a fresh locked
installation with a real notebook-generated search, progress, PNG plots and ZIP.
The only new article search was a two-generation F3A diagnostic on olig. Its
128 evaluations matched the prior native control exactly; see
`direct-run-22e273e/`. No complete article campaigns were rerun for this change.

Native local Mac checks for CAI, CSI, sequence feasibility, ViennaRNA, iRNA and
miRNA-context are recorded in `local-light-scorers.json`. All six ran together
through ordinary YAML commands without Docker or GPU; reports and PNG/PDF plots
were generated. TIE and heavy scorers are outside this final validation scope.

The later MORSE namespace and native-figure update is validated separately in
`namespace-and-presentation.json`. The local test suite remains on the
maintainer's machine and is excluded from the public Git snapshot. Five real
ready-model checks were repeated on A100 after the rename, with zero score
differences (`branding-inference/`). The updated scope diagrams and half-spider
reports preserve numeric inputs, selection and normalization; report-only
regeneration is recorded in `native-report-rendering.json`.

| Check | Measured result | Evidence |
|---|---|---|
| CPU suites | On macOS/Python 3.13.3 and olig/Python 3.10.12: 335 discovered, 330 passed, 5 optional-runtime skips on each | `reproduction.json` |
| CPU researcher quick start | Real optimization, Pareto exports and report completed | `examples/quickstart/experiment.yaml` |
| Paper source integrity | 194 compact bundle files and 155 raw files verified | `analysis/raw-release.json` |
| Archived run reanalysis | All selections matched in all 40 runs/sweeps; Figure 5 CDF contains exactly 104,239 unique transcripts | `analysis/archive-statistics.json` |
| Figure 2 statistics | All 14 correlations and 99,999-permutation p-values matched; 200,000 paired bootstrap intervals matched | `analysis/archive-statistics.json` |
| Figures | Eight PDF/PNG pairs generated; scientific layouts visually inspected | `analysis/figures.json` |
| Full fresh Figure 4a search, seed 7 | All first 12,800 events exactly match archived candidate IDs, coding sequences, generations, MFE/nt, CSI percentile and raw CSI; all selected sequences and scores match | `searches/fig4a-seed7-full.json` |
| Fresh Figure 3 diagnostic, seed 7 | Eight generations; all 4,288 initial native score/validation observations match exactly | `replay/F2_filtered_s7_initial_scores.json` |
| Fresh Figure 5 diagnostic, seed 7 | Two generations; all 1,472 initial native observations match exactly | `replay/F4_morse_s7_initial_scores.json` |
| Ready-model adapter parity | RiboNN, full-fit TE, held-out TE, half-life and RiNALMo checked on three fixture sequences each; maximum absolute score difference zero | `inference/` |
| Auxiliary panel inference | One target per panel, original context siblings retained: Mauger 12/12 and Figure 3d 7/7 native measurements exact | `assays/` |
| Native VaxPress | Pinned image built; real 64-event execution and strict budget accounting passed | `baselines/` |
| Ready assets | All 11 assets marked `paper=true`, including dependencies and all 50 RiboNN members, verified | `assets/paper-assets.json` |
| Clean clone | Bundle integrity, eight figure outputs with warnings as errors, and independent EGFP initialization/validation/search passed; tracked tree stayed clean | `reproduction.json` |

The full Figure 4a search produced 16,369 events over 256 generations. Only its
first 12,800 are part of the paper comparison. Events after that budget do not
change the paper selection. Its source revision is recorded separately from
later launcher corrections.

The initialization checks preserve the original optimizer sampling and request
order: native studies use the reference plus 63 new candidates; the 15 supplied
continuations retain their original members. Historical CPU/GPU placement and
neural batch context are preserved. Initial attempts with injected populations
or misplaced RiboDecode execution were superseded; their failed diagnostics
remain in ignored private evidence and are not reported as successful replays.

Archived analysis is independent of new inference. The two short neural-model
searches do not establish equality of their full evolutionary trajectories.
The limited auxiliary checks do not reevaluate all 39 Mauger or all 5,960
Figure 3d constructs. All 40 fresh searches have executable recipes, but were
not all rerun during this assembly. VaxPress's Python 3.12 rebuild is a tested
native baseline runtime, without a claim of binary identity to the original
Python 3.10 image.

Figure generation reconstructs the scientific plots from numerical records;
pixel identity is not claimed. The exact eight submitted PDFs are preserved
under `paper/reference/`. The conceptual Figure 1 is reused as a static vector
source. Predictor fitting code and fitting datasets are absent.

Technical logs, full search outputs and machine-specific maps stay under ignored
`validation/private/` or server runtime directories. Public receipts preserve
content hashes without private host paths. Optional non-paper model assets and
their inference runtimes are outside the validated paper module set.
