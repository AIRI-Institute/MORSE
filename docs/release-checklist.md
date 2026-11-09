# Source and container distribution

The source repository is [AIRI-Institute/MORSE](https://github.com/AIRI-Institute/MORSE),
currently private and accessible to authorized GitHub accounts. Runtime container
images are publicly available under `ghcr.io/shpakb/morse-*`; their immutable
identities are recorded in [`assets/runtime-lock.json`](../assets/runtime-lock.json).
The six neural images include ready checkpoints and required inference assets.
The source tree excludes fitting workflows, fitting datasets and checkpoint binaries.

The reproducible package contains a locked CPU environment, independent
experiment initialization, inference adapters, ready-asset verification,
paper study launchers, frozen numerical observations, archive analysis and
figure generation. Model-fitting workflows and fitting datasets are excluded.
The measured validation scope is recorded in [../validation/README.md](../validation/README.md).

Source publication uses a fresh root snapshot of the reviewed current tree
on `main`. The local assembly branches retain earlier history separately; that
history contains removed components and is not part of the published snapshot.
Do not upload those branches or their history. Subsequent source development
can use ordinary commits on the fresh `main` history.

The public release excludes the `perplexity`, `perplexity_cds`, `tie` and
`mirna-context` scorers and their container images, under any registry or tag.
The source paths and image names are recorded in
[../licenses/redistribution-exclusions.json](../licenses/redistribution-exclusions.json).

For source and model-asset distribution, retain these release records:

1. Confirm the first-party rights holder and MIT license already preserved from
   the source components. Review the exceptions in
   [../THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md) and
   [../licenses/inventory.json](../licenses/inventory.json). RiboNN, RiNALMo,
   scanMiR and the SSPL TranslationModel package have separate terms; LinearDesign
   is acquired externally. The root MIT license does not cover them.
2. Preserve the exact ready TE, half-life and RiNALMo checkpoint identities
   bundled in the neural images, together with their upstream license and
   corresponding-source materials. RiboNN and TranslationModel retain their
   separate Sanofi and SSPL notices. File-based native/source builds continue
   to accept authorized ready assets with fixed SHA256 checks; no fitting
   pipeline or fitting dataset is included.
3. Attach the optional raw evidence artifact to the eventual release and publish
   its checksum. The local canonical archive is
   `outputs/artifacts/morse-paper-raw.tar.gz`, 389,985,380 bytes, SHA256
   `7f768c397b26129657c72de12a0648d72299107f0a7c28aefbbf925c0eb893ab`.
   It contains the 155 frozen raw files and stays out of Git. The 30 MB compact
   tracked evidence is sufficient for figure reconstruction.
4. Keep the canonical repository URL and source revision synchronized with the
   notebook bootstrap and citation metadata. Access to the source currently
   requires an authorized GitHub account. Add the final article citation/DOI
   when available.
5. Preserve the published registry digests, build recipes and required notices
   recorded with the runtime images. Historical paper image IDs remain separate
   provenance records; published runtime identities are listed under
   `published_images` in the runtime lock.
