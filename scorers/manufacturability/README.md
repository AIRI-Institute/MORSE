# MORSE sequence feasibility scorer

Deterministic CPU checks for mRNA sequence constraints. Choose the
`sequence_feasibility_v1` or `manufacturing_strategy_v1` preset. The manifest ID
is `manufacturability`; its display name is **Sequence feasibility**.

`feasibility_violation_count: {role: threshold, max: 0}` is the default policy.
`feasibility_pass` is true exactly when every enabled, non-exempt rule passes.
The preset defines engineering acceptance criteria, not biological viability or
an experimental guarantee of synthesis, expression, or manufacturing success.

## Frozen preset and parameter overrides

| Check | `sequence_feasibility_v1` |
| --- | --- |
| Global GC | Inclusive 0.40–0.60 |
| Local GC | Inclusive 0.30–0.70 in 50-nt windows, stride 1; terminal window always included |
| Homopolymers | Maximum A/U run 9 nt; G/C run 5 nt |
| Direct repeats | At least two pairwise nonoverlapping exact 21-mer copies fail |
| Tandem repeats | Primitive period 2–6 nt, tract at least 18 nt and at least 3 complete copies fails |
| Forbidden motifs | Zero occurrences of supplied exact uppercase RNA motifs; empty list by default |

The preset JSON records sources and distinguishes vendor-inspired thresholds
from our design heuristics. Do not edit a released preset: add a new version.
The output records both the preset-file SHA-256 and the resolved configuration
SHA-256, including all parameter overrides and exemptions.

### Manufacturing strategy preset

`manufacturing_strategy_v1` retains ALL original rules and bounds above and adds:

| Check | Default | Meaning |
| --- | --- | --- |
| `restriction_sites` | BsaI, BsmBI; zero internal sites | Exact recognition motifs in both orientations; RNA T/U normalization |
| `capping_start` | `cleancap_ag`: prefix AG | Positions +1/+2 of the transcript, NOT the CDS start codon |
| `inverted_repeats` | 20-nt arms, <=200-nt total span, minimum gap 0 | Exact nonoverlapping reverse-complement arms; no mismatches, bulges or GU wobble |

The sources and scope limits are recorded in the immutable profile JSON. These
are constraints for a selected BsaI/BsmBI assembly and CleanCap AG workflow, not
universal manufacturing requirements. Inverted-repeat screening uses DNA Chisel's
20/200 motif definition as a structural risk heuristic; it is NOT thermodynamic
folding or proof of a hairpin. Zero gap includes reverse-complement palindromes
of at least 40 nt. All qualifying arm pairs are recorded; a longer repeat may
produce several overlapping 20-nt matches. Matches at the last nucleotide count.

Restriction enzyme names are selected through `restriction_enzymes` (YAML list;
currently BsaI and BsmBI). Custom exact recognition words can instead use
`forbidden_motifs`, explicitly supplying both orientations. Intentional assembly
flanks can be exempted by check and coordinate when they are part of the routed
sequence. The default body scope does not include vector/promoter/flanks and
therefore cannot certify the complete manufacturing DNA template.

`capping_method` accepts `none`, `cleancap_ag`, or `custom`. The custom mode requires
`custom_5prime_prefix` as exact uppercase RNA. The prefix check requires input
starting at transcript nucleotide +1: do not route a CDS-only or truncated 5'UTR
scope when capping is active. The scorer checks sequence compatibility only,
not actual capping chemistry, promoter integrity, or yield.

Set `parameters.preset: manufacturing_strategy_v1` to select the extended checks.
Omit other settings to inherit the selected preset. Existing startup environment
overrides remain supported.
Explicit overrides are validated and hashed. For direct Python use, omit an
override to inherit it. The released v1 profile bytes and its default resolved
configuration hash remain unchanged; its three new count outputs are zero.

Before imposing thresholds on a fixed protein, establish a passing synonymous
construct using the actual available UTRs and fixed start/stop codons. A failing
reference nucleotide sequence alone is not proof that the protein is infeasible.
If a requirement truly cannot be satisfied, relax or disable it explicitly before
the optimization and record the reason and final configuration. Never change
thresholds mid-run or silently alter the protein. This scorer does not perform
automatic protein-dependent relaxation; Core remains responsible for identity.

Set calculation options in YAML `parameters`. `preset` selects the registered
policy. `enabled_checks`, `forbidden_motifs`, `restriction_enzymes`, and
`exempt_regions` are ordinary YAML lists. Missing values retain the selected
preset and existing runtime settings. See the
[parameter reference](../../docs/scorer-parameters.md#sequence-feasibility)
for all names and defaults.

For GC-only ablation, set `enabled_checks` to `["global_gc","local_gc"]`.
For a matched record-only comparison, retain the same enabled checks and
parameters and change only the aggregate output role to `log`. Setting checks to
`[]` suppresses active violations and is NOT equivalent for analysis.
Disabled checks still run for observational reporting; this is not a speed
optimization. Structural validation and protein identity remain Core's separate,
mandatory responsibility.

```yaml
scorers:
  - name: sequence_feasibility
    command: '"${MORSE_PYTHON}" -m manufacturability_scorer.feasibility_server'
    parameters:
      preset: sequence_feasibility_v1
      local_gc_stride: 1
      homopolymer_max_g: 5
      homopolymer_max_c: 5
      forbidden_motifs: [GGUCUC, GAGACC]
    scope: [five_prime_utr, cds, three_prime_utr]
    metrics:
      feasibility_violation_count: {role: threshold, max: 0, value_type: int}
      feasibility_pass: {role: log, value_type: bool}
```

The standard monorepo installation includes this package. In schema 4, `cds`
includes the start and stop codons; the scope above excludes the annotated poly(A). MORSE supplies
`MORSE_PYTHON` and sends one JSON request followed by EOF per batch. Include the
selected preset and scorer sources in `provenance_paths` when using caching or
resume. Parameters apply to each request without changing later requests.
For existing registry images, use the
[updated adapter](../../docs/scorer-parameters.md#using-this-checkout-with-containers)
or rebuild locally.

## Input scope and explicit exceptions

Core routes `5′UTR + start codon + CDS body + stop codon + 3′UTR` by default,
so cross-module motifs and repeats are checked. The manifest forbids poly(A)
in this scope. DNA is normalized to RNA (`T` → `U`); ambiguous bases and empty
sequences are rejected. Custom forbidden motifs must be exact RNA; supply both orientations explicitly.
Named restriction enzymes automatically check both orientations.

Short sequences use a single full-length GC window. Repeat events count
qualifying distinct k-mers, not all pairs or maximal merged repeat blocks.
Direct-repeat `copies` is the maximum nonoverlapping copy count, while `spans`
records all occurrences. Tandem events are maximal exact periodic intervals;
reducible motifs are assigned to their primitive period. A trailing partial
copy contributes to tract length, not to complete-copy count. Homopolymers are
handled separately, not counted again as period-2/3 tandem repeats.

An optional exemption must be explicit and justified:

```json
[{"start":0,"end":75,"checks":["tandem_repeats"],"reason":"Intentionally retained validated leader"}]
```

Coordinates are **zero-based, half-open offsets in the routed sequence**.
Every span of an event must be completely contained in a matching exemption;
partial overlaps do not suppress violations. Exemptions do not delete bases,
change GC calculations, or create artificial sequence junctions. Global GC
cannot be exempted. Out-of-range exemptions are errors. Pin the scope and
module lengths when using coordinates; precheck fixed UTRs before a search.

## Outputs

- GC fractions and window count; maximum homopolymer lengths overall and per base.
- `<check>_violation_count`: enabled, non-exempt events for each check family.
- `feasibility_violation_count`, `feasibility_failed_check_count`, `feasibility_pass`.
- `feasibility_failed_checks_json` and `feasibility_violations_json`: failed
  checks and precise event coordinates, observed values and limits.
- `feasibility_raw_counts_json`: event counts including disabled/exempt checks.
- `feasibility_exemptions_json`: suppressed events and exemption reasons.
- `feasibility_profile_id`, preset/config SHA-256, `feasibility_config_json`.

Structured diagnostics are JSON strings to preserve the scalar-output contract.
All diagnostics are retained; detailed outputs can be large for highly repetitive
sequences. No result cache or external model is used.

## Running and measuring

From the monorepo root after installation:

```bash
python -m manufacturability_scorer.feasibility_server
morse-sequence-feasibility candidates.fasta --poly-a-length 100
morse-sequence-feasibility candidates.fasta --parameters overrides.json
docker build -t morse/manufacturability:4.1.0 scorers/manufacturability
```

`overrides.json` contains a JSON object with the same parameter names as the
manifest. CLI poly(A) exclusion is explicit and verifies the tail is A-only.
The process protocol is the usual ordered sequence-batch JSONL contract.

The benchmark samples 100 distinct sequences, records their hashes and lengths,
warms the scorer, and measures 100 individual calls using `perf_counter_ns`:

```bash
python -m manufacturability_scorer.benchmark \
  --fasta candidates.fasta --count 100 --poly-a-length 100 \
  --output /path/outside/source/benchmark.json \
  --command docker run --rm -i morse/manufacturability:4.1.0
```

It reports native Python computation, singleton JSONL round trips (including
serialization), a single 100-sequence batch, and startup separately. Native,
singleton and batched outputs must agree exactly. The selected body sequences
are written beside the report. Times are specific to sequence lengths, rule
settings, diagnostics volume, Python version, CPU, and current machine load.
