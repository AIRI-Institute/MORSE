# MORSE UTR library pipeline

This pipeline builds reproducible 5′ UTR and 3′ UTR libraries for MORSE
enumeration. Its source is
[GENCODE Human release 50](https://www.gencodegenes.org/human/release_50.html):
the GFF3 annotation and transcript FASTA. By default, the pipeline retains
protein-coding MANE Select transcripts, extracts and deduplicates UTRs, and can
down-select them with farthest-point sampling over length, GC content, and
simple sequence features.

`prepare_assets.py` downloads the two source files once, verifies pinned
SHA-256 digests, and stores them under
`MORSE_ASSETS_DIR/utr-library-pipeline/data/raw/gencode_v50`. When the environment
variable is unset, the cache root is `.cache/morse/assets` in the monorepo.
`build_library.py` reads that same cache and writes a table plus manifest under
`utr-library-pipeline/data/libraries` in the cache. Set `MORSE_ASSETS_DIR` to your
preferred cache root, or `MORSE_UTR_LIBRARY_PIPELINE_FILES` to an explicit UTR
asset directory. `prepare_assets.py --root DIR` overrides the cache root for
that invocation; use the same `MORSE_ASSETS_DIR=DIR` when building afterward.

```bash
python prepare_assets.py

python build_library.py --utr-type five_prime_utr \
  --min-length 30 --max-length 500 --max-candidates 100
python build_library.py --utr-type three_prime_utr \
  --min-length 30 --max-length 500 --max-candidates 100

python -m unittest discover -s tests
```

Without `--max-candidates`, the complete deduplicated library is retained. Each
result records selection parameters, input URLs and hashes, candidate count,
and the result SHA-256 digest. The committed `releases/` directory contains
only a small verification example; downloaded inputs and generated libraries
are stored in the user-selected cache outside version control.
