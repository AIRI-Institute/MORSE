# MORSE CAI scorer

Setup and schema-4 integration: [scorer guide](../../docs/scorers.md).

This is a deterministic Codon Adaptation Index (CAI) scorer for explicit
start-codon-plus-CDS-body sequences using the standard genetic code. The
project has no Python dependencies and supports the MORSE sequence-batch
JSONL contract.

`morse-scorer.yaml` documents the scope, outputs, startup parameters and
self-test. Inspect it with `morse scorers describe cai`. Prepare the runtime
before starting an experiment, then put its Bash command in the experiment YAML.
CAI runs with MORSE's Python interpreter and needs no additional packages.
Its installed entry points are `python -m scorers.cai.server` and `morse-cai`.

## Method

For each codon, relative adaptiveness is calculated within its synonymous-codon
family:

```text
w(c) = count(c) / max(count(synonymous codons))
CAI  = exp(mean(log(w(c))))
```

The implementation follows Sharp and Li: a zero count is replaced by a
profile-defined value (`0.5` by default), `ATG` and `TGG` are excluded from the
geometric mean, and stop codons are not scored. Method reference:
[Sharp & Li, Nucleic Acids Research, 1987](https://doi.org/10.1093/nar/15.3.1281).

CAI measures agreement with a particular reference codon set; it is not a
direct estimate of translation. It does not account for rare-codon positions,
local RNA structure, start-site accessibility, codon-pair bias, or the
cell-dependent tRNA pool.

## Profiles

The repository includes `homo_sapiens_kazusa_93487_cds`, built from the
[Homo sapiens, 93,487 CDS table](https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606).
The JSON stores the original counts, sample description, and retrieval date.
The profile is labeled `genome_wide_cds`: it is a general human profile, not a
classic highly expressed gene set or a tissue-specific profile.

It also includes `mus_musculus_kazusa_53036_cds`, a mouse genome-wide CDS
profile for mouse experiments. Select either schema-1 JSON through `CAI_PROFILE`;
the same profiles can be used directly by LinearDesign:

- `profiles/homo_sapiens_kazusa_93487_cds.json`
- `profiles/mus_musculus_kazusa_53036_cds.json`

For an experiment, prefer building a profile from a predefined set of highly
expressed CDS sequences in the target cell or tissue:

```bash
python -m scorers.cai.build_profile reference_cds.fasta \
  --output profiles/human_apc_high_expression_v1.json \
  --profile-id human_apc_high_expression_v1 \
  --taxon "Homo sapiens" \
  --taxon-id 9606 \
  --reference-kind highly_expressed_gene_set \
  --source-name "dataset name and selection rule" \
  --source-version "release-or-analysis-id" \
  --source-url "https://example.org/dataset"
```

`CAI_PROFILE` accepts a path to a JSON profile. Without that variable, the
built-in genome-wide human profile is used. Every result includes the SHA-256
digest of the file that was actually loaded.

## Interface

The input is the explicit start codon followed by the CDS body, without a
terminal stop codon. It must contain complete codons and use the standard
genetic code. The first codon is always treated as the initiation codon and is
excluded from the geometric mean, including when an experimentally justified
non-`ATG` start is supplied. `U` is normalized to `T` automatically.

Outputs:

- `cai`: CAI in `(0, 1]`;
- `cai_log`: natural logarithm of CAI;
- `cai_min_weight`: minimum informative codon weight in the CDS body;
- `cai_low_weight_fraction`: fraction of codons with `w < CAI_LOW_WEIGHT_THRESHOLD` (`0.1` by default);
- `cai_informative_codons`: number of codons in the mean after excluding the
  initiation codon, `ATG`, and `TGG`;
- `cai_profile_id` and `cai_profile_sha256`: reference identifiers.

Core resolves the schema-4 `start_codon + cds_body` scope before scoring. The runtime receives
only a batch of sequences; profile parameters are fixed once at startup:

```json
{"request_id":1,"sequences":["ATGTTCTTT"]}
```

## MORSE integration

```yaml
- name: cai_human
  command: '"${MORSE_PYTHON}" -m scorers.cai.server'
  timeout_seconds: 120
  scope: [start_codon, cds_body]
  metrics:
    cai: {role: optimize, direction: max}
    cai_log: {role: log}
    cai_min_weight: {role: log}
    cai_low_weight_fraction: {role: log}
    cai_informative_codons: {value_type: int, role: log}
    cai_profile_id: {value_type: string, role: log}
    cai_profile_sha256: {value_type: string, role: log}
```

`MORSE_PYTHON` is supplied by MORSE. The module and default profile are installed
with MORSE, so this command works from any directory. Each batch is one JSON
request followed by EOF; the scorer returns one JSON response and exits.

`cds_body` excludes the start and stop; schema-4 `cds` includes both. This scope
sends the start codon and body without a terminal stop to the scorer.

Choose a custom profile in the command:

```yaml
  command: 'CAI_PROFILE=/absolute/path/profile.json "${MORSE_PYTHON}" -m scorers.cai.server'
```

For an already prepared container, use `docker exec -i CONTAINER python
/app/server.py` as the command. Configure profile mounts when preparing that
container, before `morse run`.

## Verification and reproducibility

From the monorepo root, build the CPU runtime with its scorer-directory context:

```bash
docker build --platform linux/amd64 -t ghcr.io/shpakb/morse-cai:2.2.0 scorers/cai
printf '%s\n' '{"request_id":1,"sequences":["ATGTTCTTT"]}' \
  | docker run --platform linux/amd64 --rm -i ghcr.io/shpakb/morse-cai:2.2.0
python scorers/cai/prepare_assets.py --docker-image ghcr.io/shpakb/morse-cai:2.2.0
```

The image includes the adapter's MIT notice at `/app/LICENSE`.

`prepare_assets.py` verifies every built-in profile and creates a manifest with
their SHA-256 digests. When a Docker image is supplied, its archive is stored in
`${MORSE_ASSETS_DIR}/cai_files/docker`.
