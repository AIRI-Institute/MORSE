# MORSE Codon Stability Index scorer

Setup and schema-4 integration: [scorer guide](../../docs/scorers.md).

This scorer calculates a profile-specific Codon Stability Index (CSI) for the
CDS body. It expects MORSE to route `start_codon + cds_body`, excludes the first
codon, rejects stop codons, and averages the Codon Stability Coefficient (CSC)
assigned to every remaining codon.

The scorer is deliberately organism-agnostic. One process loads exactly one
profile through `CSI_PROFILE`. Run two or more instances with different
profiles and distinct scorer names to optimize the same candidate in several
contexts:

```yaml
scorers:
  - name: csi_human
    command: 'CSI_PROFILE=/absolute/path/homo_sapiens_consensus_wu2019_v1.json "${MORSE_PYTHON}" -m scorers.csi.server'
    scope: [start_codon, cds_body]
    metrics:
      csi_percentile: {role: optimize, direction: max}
      csi_raw: {role: log}

  - name: csi_mouse
    command: 'CSI_PROFILE=/absolute/path/mus_musculus_mesc_slamseq_wu2019_v1.json "${MORSE_PYTHON}" -m scorers.csi.server'
    scope: [start_codon, cds_body]
    metrics:
      csi_percentile: {role: optimize, direction: max}
      csi_raw: {role: log}
```

CSI and its default human profile are installed with MORSE and need no extra
Python dependencies. The module command works from any directory. Replace the
profile paths with your JSON files. To locate the bundled profiles directory:

```sh
python -c 'from importlib.resources import files; print(files("scorers.csi") / "profiles")'
```

MORSE supplies `MORSE_PYTHON`, sends one JSON request and EOF per batch, and
reads one response. Startup environment settings belong in the command.

Core exposes these as `csi_human.csi_percentile`, `csi_human.csi_raw`,
`csi_mouse.csi_percentile`, and `csi_mouse.csi_raw`. The scorer itself always returns
the same local output names.

## Outputs

- `csi_raw`: arithmetic mean CSC; this is the direct codon-stability signal.
- `csi_percentile`: the raw CSI mapped to the profile's reference CDS CDF. This
  is the default objective because its scale is easier to interpret across
  separately calibrated profiles.
- `csi_window_p05` and `csi_window_min`: diagnostics for locally poor regions.
- `csi_destabilizing_fraction`: fraction of codons with CSC below the configured
  threshold (zero by default).
- `csi_profile_id` and `csi_profile_sha256`: exact profile provenance.

The bundled human consensus coefficients are the median of the six CSC columns
published as [Wu et al. Figure 1 source data 2](https://doi.org/10.7554/eLife.45396.006).
The bundled mouse mES profile is derived from the mouse SLAM-seq records bundled
with [iCodon](https://github.com/santiago1234/iCodon) at the revision recorded
inside the profile. The corresponding papers and sequence counts are also
recorded in each JSON file.

## Ready profiles for independent experiments

Select a ready human, mouse, tissue or cell-state profile using `CSI_PROFILE`.
Custom profile JSON must follow schema version 1, include all 61 sense-codon
coefficients and retain its source identity and reference CDF. The scorer validates
and hashes that exact profile at startup. This inference-only release excludes
utilities that estimate coefficients from labelled stability measurements.

CSI models one codon-mediated component of stability. It does not model UTRs,
RNA structure, innate sensing, translation initiation, formulation, dose, or
clinical efficacy, so those remain separate objectives and validation layers.

## Direct inference

```bash
printf '%s\n' '{"request_id":1,"sequences":["ATGTTCTTT"]}' | python -m scorers.csi.server
```

`morse-csi` is the equivalent command in an activated MORSE environment.

## CPU container

From the monorepo root:

```bash
docker build --platform linux/amd64 -t ghcr.io/shpakb/morse-csi:1.0.0 scorers/csi
printf '%s\n' '{"request_id":1,"sequences":["ATGTTCTTT"]}' \
  | docker run --platform linux/amd64 --rm -i ghcr.io/shpakb/morse-csi:1.0.0
```

The image includes the bundled profiles and the adapter's MIT notice at
`/app/LICENSE`. Select another bundled profile through `CSI_PROFILE`, or mount
a custom JSON profile read-only and set `CSI_PROFILE` to its container path.
