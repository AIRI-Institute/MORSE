# ViennaRNA: RNA folding energy

`mfe` is the predicted minimum free energy of the scored RNA (kcal/mol).
More negative means a lower-energy predicted fold. Use `mfe` for the
total energy of the region, or `mfe_per_nt` (kcal/mol per nucleotide) to
compare energy per base when candidate lengths differ. At fixed length,
minimizing either metric gives the same ranking.

**Input:** the RNA region selected by `scope`. The default, `full_sequence`,
includes poly(A) if present. Runs on CPU with the `mfe` extra; no model weights
are required.

## Use in MORSE

From the repository root:

```sh
uv sync --python 3.13 --all-packages --extra mfe --frozen
source .venv/bin/activate
```

For pip installation, see [Install MORSE](../../README.md#install).
This example scores the transcript body without poly(A):

```yaml
scorers:
  - name: mfe
    command: '"${MORSE_PYTHON}" -m scorers.viennarna.scorer'
    scope: [five_prime_utr, cds, three_prime_utr]
    parameters:
      max_bp_span: 300
    metrics:
      mfe: {role: optimize, direction: min}
      mfe_per_nt: {role: log}
```

Try the [bundled four-scorer CPU example](../../examples/scorer-parameters/experiment.yaml):

```sh
morse score --config examples/scorer-parameters/experiment.yaml
```

Under `outputs`, `mfe.mfe` is the starting construct's energy for this scope.

## Choose folding settings

- `max_bp_span` limits the distance between paired bases in nucleotides
  (default 300; 0 allows unrestricted pairing).
- `temperature` sets the folding temperature in °C (default 37).
- Set `compute_ensemble: true` to obtain `ensemble_free_energy` and
  `ensemble_diversity`. Their calculation requires `dangles: 0` or `2`
  (default 2). Both outputs are absent when ensemble calculation is off.

Use the same scope, temperature and pairing settings when comparing scores.
`morse scorers describe viennarna` lists the remaining settings and outputs,
including the dot-bracket structure of the minimum-energy fold.

## Method and runtime

The adapter uses the [ViennaRNA Package](https://github.com/ViennaRNA/ViennaRNA)
(v2.7.2; [Lorenz et al., 2011](https://doi.org/10.1186/1748-7188-6-26)).
It predicts folding under a thermodynamic model; the predicted energy is not
a measurement of transcript half-life.

The native command above accepts YAML `parameters`. The published 2.2.0
container uses environment settings such as `VIENNARNA_MAX_BP_SPAN` and
`VIENNARNA_ENABLE_PF`. See [using the current adapter with an existing image](../../docs/scorer-parameters.md#using-this-checkout-with-containers)
for per-request parameters.

[Adapter license](LICENSE) · [ViennaRNA license](ViennaRNA-COPYING.txt)
