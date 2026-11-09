These receipts record real inference checks on one target construct per assay,
with the canonical batch siblings retained to preserve model input shapes.

| Diagnostic | Scorer adapters | Compared measurements | Exact matches |
|---|---:|---:|---:|
| `rescore-mauger-limit1` | 11 | 12 | 12 |
| `rescore-fig3-limit1` | 7 | 7 | 7 |

All compared native values have zero absolute difference from the archived paper
observations. ViennaRNA body scoring produces both MFE and GC, hence 12 Mauger
measurements from 11 adapters. Context siblings were scored but are not counted
as independently validated targets.

Each directory contains the native outputs, explicit reference comparisons and
a receipt with image IDs, scorer source hashes, verified checkpoint/profile
hashes, devices and request sizes. Hashes of the original technical artifacts
link these compact records to the locally retained raw evidence. Private host
paths and full technical logs are excluded from these public records.

These are limited execution checks. They do not establish fresh agreement for
all 39 Mauger constructs, all 5960 retrospective Fig3 constructs, their complete
statistics or the full optimization campaigns. Full replay commands are
documented under `docs/` and `paper/`.
