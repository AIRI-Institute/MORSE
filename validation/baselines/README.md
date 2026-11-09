`vaxpress-native64` records a real native VaxPress execution using the shared
MFE/CSI adapter, seed 7 and alpha 0, stopped after 64 scored events. Its 64 event
rows contain 58 distinct CDS hashes. Three batches scored 20 events each; the
last scored four and discarded 16 unscored proposals. The native stop reason is
`event_budget`.

The receipt independently checks event counts, budget accounting, sequence
lengths, the best fitness against the event log and alpha-zero fitness against
CSI/100. Its scorer/adapter hashes match the runtime provenance. `events.csv`
contains measured scores and CDS hashes; `budget_batches.json` preserves the
budget records. Runtime image ID and source hashes are recorded in the receipt.

This is a native execution and budget check, not a completed 2560-event paper
baseline or the full 25-run sweep. It does not claim matching the full published
baseline trajectory. Private technical logs are retained outside the public
validation records.
