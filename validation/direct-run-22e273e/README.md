# Plain MORSE run with prepared scorers

Exactly one article diagnostic ran on olig on 2026-11-06: F3A seed 7, limited
to **two generations and 128 evaluations**. The original public YAML was
validated first. Its diagnostic copy changed only the generation budget,
disabled early stopping and resolved its existing input paths before relocation.

The ordinary `morse run --config … --output-dir …` command called the YAML's
literal curl commands. ViennaRNA and CSI backends were prepared independently
before the optimizer. Both retained their process IDs and received batches
**1, 63, 64**. The frozen first-generation population and evaluation-order audit
passed. Container identity, start time and restart count remained unchanged.
Docker CLI calls were denied in the optimizer's environment; none occurred.

All 128 candidate IDs, full and coding sequences, generations, evaluation
indices, MFE/nt, CSI percentile and raw CSI matched the same prefix of the prior
full native control exactly. `first128-comparison.csv` exposes those comparisons
without duplicating the full transcripts.

Source commit `22e273eaf57df417fcb5c68e5d2be1ac63e60088` came from a verified Git
bundle and remained clean. A pre-existing immutable ViennaRNA image supplied
the runtime with the current read-only scorer source mounted. This was not a
new image build or a full paper campaign. Both HTTP utilities exited cleanly;
the one container created for this check was removed and both ports released.

`receipt.json` records source/config/log identities, model presets, exact-match
counts and lifecycle evidence. Full technical logs and result files remain in
ignored `validation/private/direct-run-22e273e/`.
