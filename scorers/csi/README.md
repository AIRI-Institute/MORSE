# CSI: codon stability index

`csi_raw` is the mean codon stability coefficient for the CDS body.
`csi_percentile` places that value on the selected profile's reference
distribution (0–100). A score of 90 means approximately the 90th percentile
among reference CDSs, not a 90% increase in RNA half-life.

**Input:** start codon + CDS body, without UTRs or stop codons. The first
codon is skipped. Runs on CPU in the [MORSE environment](../../README.md#install).

## Use in MORSE

```yaml
scorers:
  - name: csi
    command: '"${MORSE_PYTHON}" -m scorers.csi.server'
    scope: [start_codon, cds_body]
    parameters:
      profile_path: profiles/homo_sapiens_consensus_wu2019_v1.json
    metrics:
      csi_percentile: {role: optimize, direction: max}
      csi_raw: {role: log}
      csi_profile_id: {role: log, value_type: string}
```

From the repository root, try the
[bundled four-scorer CPU example](../../examples/scorer-parameters/experiment.yaml)
(MORSE installed with the `mfe` extra):

```sh
morse score --config examples/scorer-parameters/experiment.yaml
```

This example logs `csi.csi_raw` and `csi.csi_profile_id`; the YAML above also
uses the percentile as an optimization objective.

## Choose a profile

Two profiles are included:

- [Human consensus](profiles/homo_sapiens_consensus_wu2019_v1.json), the default:
  the median of six published coefficient columns ([Wu et al., source data](https://doi.org/10.7554/eLife.45396.006)).
- [Mouse embryonic stem cells](profiles/mus_musculus_mesc_slamseq_wu2019_v1.json):
  use `profile_path: profiles/mus_musculus_mesc_slamseq_wu2019_v1.json`.
  The coefficients are derived from mouse SLAM-seq records distributed with
  [iCodon](https://github.com/santiago1234/iCodon).

Bundled paths resolve inside the scorer. For another tissue or cell type,
provide your own profile with the 61 sense-codon coefficients, source identity
and reference distribution in the same JSON format. The score describes only
the codon contribution captured by that profile.

`window_codons` changes the local-window diagnostics; `destabilizing_threshold`
changes the fraction of codons flagged. Neither changes the overall CSI.
Use `morse scorers describe csi` for all settings and outputs.

[Adapter license](LICENSE) · [Container and parameter compatibility](../../docs/scorer-parameters.md#using-this-checkout-with-containers)
