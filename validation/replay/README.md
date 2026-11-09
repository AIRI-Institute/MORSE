# Initial-generation numerical replay

These receipts independently compare the same 64 full transcripts against the
canonical archived experiment logs. They validate real model inference and the
initialization protocol for source commit
`8763d174066a00c1ad0d0d19e9c2ef7f60ebcbf5`.

| Run | Diagnostic generations | Fields × transcripts | Exact comparisons | Result |
| --- | ---: | ---: | ---: | --- |
| `F2_filtered_s7` (Figure 3) | 8 | 67 × 64 | 4,288 / 4,288 | Pass |
| `F4_morse_s7` (Figure 5) | 2 | 23 × 64 | 1,472 / 1,472 | Pass |

Every numeric value matches exactly; the additional check of parsed float64 bit
patterns also passes (960 F2 values and 384 F4 values). Structured count/flag
fields match semantically. No field is missing or unknown, and the audit uses
its original fixed per-metric tolerances with no overrides.

Both runs use native `ConstructSampling`, with the frozen initial files used
only for auditing. Full transcript order and the original scorer request sizes
of 1 followed by 63 are verified. Figure 5 RiboDecode runs on CPU with inner
batch size 256; Figure 3 RiboDecode runs on CUDA with inner batch size 8.
Other learned scorers retain their historical CUDA device policy.

The receipts include hashes of the archived files, immutable image IDs through
hashed private provenance, the launcher protocol, actual evidence files, and
source overlays verified against the commit. Absolute private paths are omitted.
Full original copied run evidence remains in the maintainer's excluded
`validation/private/canonical-runs/` directory.

Run the checked-in audit against a newly copied experiment directory:

```bash
uv run python scripts/compare_initial_scores.py \
  --run-root PATH_TO_COPIED_RUN/F4_morse_s7 \
  --raw-root paper/artifacts/raw \
  --output NEW_REPORT.json
```

The current helper also accepts schema-4 logs. It uses the frozen per-run metric
mapping, checks renamed threshold fields, and reconstructs F3A transcripts from
the current run configuration. It reads existing evidence without running
scorers. The receipts above retain their original source and script hashes.

These are initial-generation inference checks on diagnostic runs. They do not
establish reproduction of all later search trajectories or all 40 paper series.
