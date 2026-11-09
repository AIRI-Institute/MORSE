# Audit the initial score replay

Use this helper after a paper run or a diagnostic generation override to compare
inference on the **same first-generation 64 transcripts**. It matches complete
transcript sequences, allowing candidate numbering and row order to differ. It
compares score fields, feasibility flags, native counts, structured observations,
profile IDs/hashes and structural checks that both logs contain.

```sh
uv run python scripts/compare_initial_scores.py \
  --run-root outputs/replay/F2_filtered_s7 \
  --raw-root paper/artifacts/raw \
  --output outputs/audits/F2_filtered_s7.initial.json
```

The run ID is inferred from `plan.json` or the directory name; use `--run-id ID`
when supplying a campaign root. The output file must be new. Run
`morse paper verify --raw` separately for complete archived-file checksum
verification: this helper stops after generation one and reports the registered
expected file hash without claiming to have rehashed a potentially large archive.

Schema-4 `scorer.metric` and threshold names are matched to historical columns
using the frozen mapping in `paper/provenance/schema4_migration.json`. Each report
records that mapping. Scorer prefixes are never guessed or merged; unmapped and
missing observations prevent a pass. Transcript `gc_fraction` and a scorer's
`output_gc_fraction` are checked separately.

Each metric reports exact matches, paired/missing observations, maximum absolute
and relative differences, tolerance, and up to five mismatch examples. Relative
differences with a zero archived value are counted separately. The fixed defaults
are explicit in `compare_initial_scores.py`: deterministic ratios use
`atol=rtol=1e-12`; native counts/lengths and flags are exact; MFE uses `atol=1e-5`
and MFE per nucleotide uses `atol=1e-7`; learned score tolerances are declared per
metric. These thresholds never adapt to observed results. To apply a separately
justified policy, supply a JSON file with exact names:

```json
{"orthrus_translation_efficiency.translation_efficiency_score": {"atol": 0.00001, "rtol": 0.00001}}
```

Use `--tolerances FILE`; the selected override and its file checksum are recorded.
Both current qualified names and the historical names shown in the report work.
Unknown numeric metrics and missing observations prevent a pass. Known numeric
fields must contain finite numbers; nulls, strings and booleans cannot confirm
numeric agreement. Unset threshold bounds may remain null. Exit codes are
`0` pass, `1` mismatch, `2` invalid input, and `3` partial evidence. Candidate IDs
alone never establish equality.

The older F3A MORSE benchmark prefixes stored coding sequences and three native
scores. Their full transcripts are reconstructed and structurally validated using
the frozen fixed UTR/start/stop/poly(A) configuration. Additional fields absent
from that archive remain explicitly unavailable, even when the three recorded
scores match. Native VaxPress has a different initialization protocol and is not
accepted by this 64-member Core check.

For an older interrupted/failed run without a JSONL export, `--allow-checkpoint`
can read its locally trusted `optimizer-checkpoint.pkl`. This deserializes Python
pickle; use it only for your own trusted checkpoints. It may read the full
checkpoint, while the archived JSONL remains bounded to generation one.

This checks numerical inference replay on identical inputs. A full optimization
trajectory, selection, budget and stopping-rule comparison uses
`morse paper analyze --runs-root DIR` separately.
