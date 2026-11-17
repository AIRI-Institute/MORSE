# CAI: codon adaptation

The Codon Adaptation Index (`cai`, 0 < CAI ≤ 1) measures how closely a CDS
uses the preferred synonymous codons in a reference profile. Higher means
closer agreement with that profile; CAI does not predict protein yield.

**Input:** start codon + CDS body, without the stop codon or UTRs. The scorer
skips the first codon and excludes AUG and UGG from the calculation. At least
one informative codon must remain. Runs on CPU in the [MORSE environment](../../README.md#install).

## Use in MORSE

```yaml
scorers:
  - name: cai
    command: '"${MORSE_PYTHON}" -m scorers.cai.server'
    scope: [start_codon, cds_body]
    parameters:
      profile_path: profiles/homo_sapiens_kazusa_93487_cds.json
    metrics:
      cai: {role: optimize, direction: max}
      cai_profile_id: {role: log, value_type: string}
```

To try it from the repository root, score the starting construct in the
[bundled four-scorer CPU example](../../examples/scorer-parameters/experiment.yaml)
(MORSE installed with the `mfe` extra):

```sh
morse score --config examples/scorer-parameters/experiment.yaml
```

Look for `cai.cai` and `cai.cai_profile_id` under `outputs`.

## Choose a profile

The bundled profiles describe genome-wide codon usage, not a specific tissue.
Human is the default. For mouse, replace the profile path with
[`profiles/mus_musculus_kazusa_53036_cds.json`](profiles/mus_musculus_kazusa_53036_cds.json).
Bundled paths resolve inside the scorer; for a custom profile use an absolute
path accessible to the scorer. The MORSE LinearDesign wrapper accepts the same
JSON and converts it to the CSV used by the LinearDesign binary.

`low_weight_threshold` changes the diagnostic `cai_low_weight_fraction`, not
the CAI itself. Use `morse scorers describe cai` to see all outputs and settings.

## Method and source

CAI is the geometric mean of the informative codons' relative frequencies
within their synonymous families ([Sharp & Li, 1987](https://doi.org/10.1093/nar/15.3.1281)).
The [human profile](profiles/homo_sapiens_kazusa_93487_cds.json) records the
source [Kazusa table](https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606);
both profile files record their source data. To build a profile from your own
reference CDS set, see `python -m scorers.cai.build_profile --help`.

[Adapter license](LICENSE) · [Container and parameter compatibility](../../docs/scorer-parameters.md#using-this-checkout-with-containers)
