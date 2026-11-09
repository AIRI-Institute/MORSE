# Paper experiments and figures

This directory packages the numerical evidence and inference recipes used for
the final MORSE paper. It contains no manuscript text, manuscript PDF, predictor
training implementation, or training dataset. The published neural images bundle ready checkpoints; `morse assets` supports
separate native runtimes, custom asset roots and source builds.

Run commands from the full repository checkout after installing its Python
environment. A wheel installation can point to the checkout with `MORSE_ROOT`.

```bash
morse paper list
morse paper verify
morse paper figures --output results/paper
morse paper analyze --output results/paper/analysis
```

`figures` rebuilds Figures 2–5 and Extended Data Figures 1–3 from frozen numerical
observations. It emits all seven scientific figures, the static vector Figure 1,
PDF/PNG files, the observations actually plotted, generated data tables and a
receipt. The scientific panel content and selections are preserved; the plotting
layout is reconstructed and is not claimed to be pixel-identical to the paper.
`reference/` holds the exact eight submitted figures for comparison. Figure 1 is
a conceptual vector diagram and is explicitly reused as a static diagram source.
No previous numerical PDF is used as an input to the numerical plots.

The compact tracked evidence is approximately 30 MB. It covers **40 independent
runs or sweeps**: eight groups with five independent seeds, 7, 43, 101, 137 and 211.
A VaxPress repeat consists of one five-weight sweep; its five weights are not
independent replicates. The registry retains the historical run IDs so every
selected construction can be traced to an archive.

| Command study | Final paper | Conditions |
|---|---|---|
| `fig3` | Figure 3 and Extended Data 1 | Filtered and logging-only, five matched pairs |
| `fig4a` | Figure 4a | MORSE and VaxPress, five runs/sweeps each |
| `fig4b` | Figure 4b and Extended Data 2 | Fixed and editable UTRs, five matched pairs |
| `fig4c` | Figure 4c and Extended Data 3 | Five native-start continuations |
| `fig5` | Figure 5 | Five three-predictor searches |

Figure 2 uses 39 archived luciferase constructs and fixed predictions. Its exact
rank composites, permutation protocol and paired bootstrap inputs are in
`data/SupplementaryFigure3_Mauger_U/`. Recompute all 200,000 paired bootstrap
resamples with:

```bash
morse paper analyze --bootstrap --output results/paper/statistics
```

## Replay the auxiliary ready-scorer inference

The two retrospective panels have executable inference recipes:

```bash
morse paper rescore mauger --output results/mauger-plan \
  --image-map assets/runtime-lock.json
morse paper rescore mauger --output results/mauger-inference --execute \
  --device cuda:0 --image-map assets/runtime-lock.json
morse paper rescore fig3-logged --output results/figure3-logged-inference --execute \
  --device cuda:0 --image-map assets/runtime-lock.json
```

Pull the [published scorer images](../docs/scorers.md#container-images) first.
The explicit runtime-lock image map selects immutable GHCR references. Its six
`self_contained: true` neural images use their bundled checkpoints, backbone,
packages and inference source without a host asset cache or setup mounts.
Optional `--source-overlay` replaces only adapter source files and keeps those
bundled assets. To use an explicitly prepared external asset tree, supply
`--asset-root /path/to/assets`; this retains file verification and matching
runtime mounts. Native processes, legacy maps and historical defaults keep
their existing asset behavior. Omitting the image map preserves the archived
runtime choices and never substitutes a current image for a historical ID.

`mauger` reevaluates all 39 constructs with the 12 measurements in Figure 2:
RiboNN, held-out Orthrus TE, RiboDecode, half-life, CSI, LinearDesign-profile CAI,
body MFE/nt, CDS MFE, scanMiR, RiNALMo, feasibility violations and GC fraction.
The complete inference output regenerates both fixed rank composites and all
14 correlation/permutation statistics. Native scores and their directions are
retained separately in `comparison.csv`.

`native_scores.csv` always retains native model scales; `oriented_scores.csv`
separately applies the prespecified Figure 2 directions. A completed replay with
score mismatches retains all evidence and returns exit status 1 so shell/CI
cannot treat a numerical discrepancy as a parity pass.

`fig3-logged` reevaluates the 5,960 unique saved constructs with the seven
Figure 3d outputs. The compact panel is derived from the raw
`F2_filtered_s211` archive; its source hash, exact original input order and all
22,784 population memberships are preserved. The command rebuilds means over
the same 64 members at each generation 1–356, including infeasible members,
and the direction-aware retrospective display normalization. It writes
`Figure3_modular_design.d_series.csv`, generation-level comparisons and the
new normalization. This does not modify search objectives or select candidates.

For a short diagnostic add `--limit 1`. A limited diagnostic still evaluates the
original sibling rows needed to preserve model batch context, and marks them
`context_only` in its raw responses. In particular, the held-out TE model can
depend on padded batch length; canonical row order, four-row TE requests and
its original inner batching are preserved. An explicit `--batch-size` override
is labelled diagnostic and records the altered protocol. No tolerances are
increased to hide a batch or numerical discrepancy.

The image map additionally accepts `half-life`, `immunogenicity` and `irna`.
The latter requires a real R/scanMiR 1.16.0 runtime and the frozen human
`core_975` panel. It is not replaced with a Python approximation. Learned
checkpoints in external asset trees and bundled profiles are checked before
inference. Self-contained images are selected by immutable digest and retain
their build-verified asset identities. Missing images, required external
checkpoints, runtime packages or panel mismatches produce a failed receipt.
`progress.json` records current scoring progress; native responses, request IDs,
context siblings, timings and observed differences remain reviewable.

`--runtime process` uses installed real Python/R runtimes. To exercise selected
light scorers, use repeatable `--scorer` selections such as `cai`, `csi`,
`manufacturability`, `viennarna-body` or `viennarna-cds`; a subset is always
labelled diagnostic. Full assay recipes are plan-only until `--execute` is set.

## Complete archived evaluations

The optional raw-source release bundle contains 155 files, approximately 609 MB
before archive compression. These files are excluded from Git and live in
`paper/artifacts/raw/` after extracting the release bundle. The manifest is
`provenance/raw_sources.json` and verifies every source with SHA256.

```bash
morse paper verify --raw
morse paper analyze --raw --output results/paper/archive-analysis
```

The full archive audit recomputes all specialist, HV, balanced and maximin
selections from the complete in-budget evaluations. It requires every selected
candidate ID to match the frozen selection and rebuilds the Figure 5 empirical
reference of exactly 104,239 unique full transcripts. This reanalysis makes no
model calls and does not rerun evolutionary searches.

Create or extract the independent release artifact with manifest-checked paths:

```bash
morse paper archive pack --file dist/morse-paper-raw.tar.gz
morse paper archive extract --file /path/to/morse-paper-raw.tar.gz
```

## Run the paper configurations directly

Prepare your scorer environments independently, then use the same command as
for any other experiment:

```bash
morse validate --config paper/configs/f2_filtered_s7.yaml
morse run --config paper/configs/f2_filtered_s7.yaml
```

Each schema-4 YAML specifies `output_dir`, a reference RNA FASTA, its part
coordinates and the allowed `changes`. UTR libraries retain their original order.
Each scorer has a plain Bash `command`, sequence `scope` and `metrics`. The
command receives one JSON batch on stdin through EOF and must return one JSON
response on stdout. It can call your existing API, container or executable;
MORSE does not create its runtime. The example paper commands use `curl` against
localhost endpoints. Change their addresses to match your prepared services.
[SCORERS.md](SCORERS.md) lists the frozen model/device/batch presets and example
ports needed for reference reproduction.

For a scorer that already speaks HTTP, no adapter is needed. A ready JSONL
server can optionally be exposed with:

```bash
morse scorers serve --command 'docker exec -i EXISTING_CONTAINER python -u /app/server.py' --port 8703
```

This explicit utility retains one backend process and its loaded model across
requests. It constructs no Docker commands, creates no containers and never
restarts a stopped scorer. Stop the utility with Ctrl-C when finished.

All 35 MORSE YAMLs load directly. The 20 native searches use ordinary sampling;
their frozen initial files only audit the first generation. The 15 supplied
Figure 4b/4c continuations retain their complete structurally valid populations,
including threshold failures. Portable fingerprints and exact vector-to-sequence
checks permit relocation without weakening this audit or the hard constraints.
A short diagnostic changes the YAML generation count; it is not a full paper
rerun. No-feasible outcomes still retain their complete evidence.

[The schema migration receipt](provenance/schema4_migration.json) records the
offline equivalence checks for all 35 configurations and the mapping from new
`scorer.metric` outputs to historical observation columns. Reference sequences,
search bounds, mutable codon/UTR choices, scorer inputs and all 2,240 initial
candidate vectors and sequences were checked without inference. Original initial
population files remain unchanged; 20 schema-4 copies bind the same vectors to
the new fixed-part labels after these checks. Frozen numerical observations and
submitted reference figures are unchanged.

VaxPress is a separate native baseline:

```bash
morse paper baseline vaxpress --seed 7 --runtime container --image-map images.json \
  --output outputs/paper-runs --execute
```

Its default sweep has five weights × 2,560 scoring events = 12,800 events.
`--scoring-events 320` explicitly selects a short 64-event-per-weight diagnostic.
Figure 4b uses the frozen LinearDesign lambda=3 start; Figure 4c uses five frozen
RiboDecode starts. [UPSTREAM.md](UPSTREAM.md) describes separate native regeneration.

## Frozen scientific distinctions

* Full-fit and held-out Orthrus TE have different weights and CDS masks. Their
  complete mapping is in `provenance/checkpoint_usage_0108.json`.
* Online stopping HV, candidate-selection HV, archive reporting HV, empirical
  percentiles and radar display scales are separate quantities. Current numerical
  display records and normalizations are bundled; SD uses denominator n−1.
* Early-stopped trajectories are summarized only over common observed support;
  missing correlations remain missing.
* Figure 4a uses 12,800 actual scoring events per MORSE run or VaxPress sweep,
  including repeats. A VaxPress search uses 2,560 events per weight.
* Figure 5 historical constructs contain A120; the final schematic standardizes
  the displayed tail to A100. All scorer inputs and sequence-feasibility checks
  exclude poly(A). Historical sequences, fingerprints and selection IDs are
  preserved rather than silently changed.

The code records the runtime Python/library versions and the exact scorer/asset
identity of new experiments. Identical evolutionary trajectories across different
hardware/library versions are not assumed. Frozen archive reanalysis is exact;
new inference and searches require comparison with the archived observations.

Compare fresh Core outputs with the archived specialist/HV/balanced/maximin
selection rules and numeric observations:

```bash
morse paper analyze --runs-root outputs/paper-runs --output outputs/fresh-analysis
```

Store each run under its historical run ID, such as `F2_filtered_s7`.
The reader accepts ordinary Core result directories as well as the older
nested result layout. It reports budgets, sequences, selected IDs and numeric
differences independently; a seed alone is not a reproducibility guarantee.

When all five fresh Figure 5 archives are available, selection rebuilds their
own pooled empirical CDF. A partial run uses the frozen historical CDF for an
explicitly labelled partial/diagnostic comparison. `paper figures` always draws
the frozen paper evidence; it does not silently substitute new experiments into
the published figures.
