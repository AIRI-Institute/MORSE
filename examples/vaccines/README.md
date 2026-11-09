# Vaccine examples

| Example | RNA | Source | Run config |
|---|---|---|---|
| MAGE-A1 | [FASTA](magea1/sequence.fasta) | [README](magea1/README.md) | [YAML](magea1/experiment.yaml) |
| MAGE-A6 | [FASTA](magea6/sequence.fasta) | [README](magea6/README.md) | [YAML](magea6/experiment.yaml) |

Run from the monorepo root:

```sh
uv run morse run --config examples/vaccines/magea1/experiment.yaml
```
