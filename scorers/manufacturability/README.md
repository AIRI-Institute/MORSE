# Sequence feasibility

This CPU scorer counts violations of sequence design rules. Zero
`feasibility_violation_count` means the input passes every enabled check
after explicit exemptions. This is a test of the selected policy, not a
guarantee of synthesis, expression or manufacturing success.

**Input:** the transcript body, including UTRs, start and stop, without poly(A).
Rules apply across part boundaries. The scorer ships with the
[MORSE environment](../../README.md#install); its catalog ID is `manufacturability`.

## Use in MORSE

Add this entry alongside your optimization objectives. It treats zero
violations as a hard constraint and logs which rules fail.

```yaml
scorers:
  - name: feasibility
    command: '"${MORSE_PYTHON}" -m manufacturability_scorer.feasibility_server'
    scope: [five_prime_utr, cds, three_prime_utr]
    metrics:
      feasibility_violation_count: {role: threshold, max: 0, value_type: int}
      feasibility_failed_checks_json: {role: log, value_type: string}
```

Try the [bundled four-scorer CPU example](../../examples/scorer-parameters/experiment.yaml)
from the repository root (MORSE installed with the `mfe` extra):

```sh
morse score --config examples/scorer-parameters/experiment.yaml
```

Look for `feasibility.feasibility_violation_count` under `outputs`. The bundled
example logs the count; the YAML above also uses it to filter candidates.

## Restriction sites

Add recognition sequences under the scorer's `parameters`. For example, to
include BspQI in the checks above:

```yaml
parameters:
  restriction_sites: [GCTCTTC]
```

Give each site once, in either orientation. The scorer checks both the motif
and its reverse complement across the selected `scope`. It counts overlapping
matches at different positions and counts a palindromic site once per position.
Use exact DNA or RNA sequences (T and U are equivalent).

With the default preset, supplying sites adds this check to the existing
checks. If you set `enabled_checks` explicitly, include `restriction_sites`
in that list to count its violations. The count is available separately as
`restriction_sites_violation_count` and contributes to `feasibility_violation_count`.

## Read a failure

`feasibility_failed_checks_json` lists the check names. For event coordinates,
observed values and limits, log `feasibility_violations_json` as a string.
Both are JSON-encoded diagnostics. Positions are zero-based, end-exclusive
offsets within the **scored sequence**. The violation count is the number of
events, not the number of distinct failed checks.

## Choose a preset

The default [`sequence_feasibility_v1`](manufacturability_scorer/profiles/sequence_feasibility_v1.json)
checks global and local GC, homopolymers, direct repeats, tandem repeats and
user-specified forbidden motifs (empty by default).
[`manufacturing_strategy_v1`](manufacturability_scorer/profiles/manufacturing_strategy_v1.json)
adds restriction sites, a capping-start requirement and inverted-repeat checks
for the specified manufacturing workflow. A capping-start check requires input
beginning at transcript nucleotide +1.

Set `parameters.preset` to choose a preset. For a record-only experiment keep
the same preset and change the count's role to `log`, so the same rules are
evaluated. For thresholds, units and per-check settings, see the
[parameter reference](../../docs/scorer-parameters.md#sequence-feasibility)
or run `morse scorers describe manufacturability`. Each result includes the
preset and effective configuration hashes.

For a standalone assembled FASTA, `morse-sequence-feasibility --help`
describes the CPU command and how to specify an annotated poly(A) tail.
The preset files above include the rule sources and assumptions.

[Adapter license](LICENSE) · [Container and parameter compatibility](../../docs/scorer-parameters.md#using-this-checkout-with-containers)
