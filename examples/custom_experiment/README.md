# Start an independent experiment

The included `input_cds.fasta` is the MLANA CDS body, without initiation and
terminal stop codons. Run from the monorepo root:

```sh
morse init --cds examples/custom_experiment/input_cds.fasta --output-dir experiments/my_antigen
morse run --config experiments/my_antigen/experiment.yaml
```

Replace the FASTA with your own CDS body. The generated config inserts AUG and
UAA; edit `construct.start_codon` and `construct.stop_codon` to allow alternatives. Add repeated
`--five-prime-utr FILE` and `--three-prime-utr FILE` flags for UTR alternatives.
See [the researcher guide](../../docs/researcher-guide.md) for real inference
scorers, explicit scope/objective policies, provenance, multiple seeds and resume.

The generated configuration uses a Bash `command` string, optional calculation
`parameters` and a `metrics`
mapping with `optimize`, `threshold` or `log` roles. Prepare any replacement
scorer environment, model, container or HTTP service before `morse run`.
Inspect its interface with `morse scorers describe ID` and generate a section
with `morse scorers template ID --command 'YOUR_COMMAND'`.

The [commented template](../experiment.template.yaml) shows the main sections.
Change `output_dir` in the YAML to select the results folder.
