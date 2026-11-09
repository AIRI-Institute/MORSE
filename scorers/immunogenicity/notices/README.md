# Public runtime image notices

The MORSE inference adapter and classification head are provided under the
accompanying MIT license (`LICENSE`). MultiMolecule and its RiNALMo
implementation retain their own AGPL terms; the MIT notice does not relicense
that dependency or a combined covered work.

During the build, the original MultiMolecule 0.2.0 `license.md` and
`license.header` are extracted from the exact source distribution with SHA-256
verification and retained under `/usr/share/morse/notices/multimolecule`.
The matching wheel's complete `METADATA`, which also embeds the full AGPL
license, and notice provenance are retained there. Its unmodified wheel and source
distribution, archive SHA-256 digests, and complete installed package-version
list are retained under `/usr/share/morse/sources`. The image also contains the
adapter's Python source in `/app` and the exact Dockerfile under that source
folder. Preserve these materials in redistributed copies and comply with the
applicable source and network-service terms.

`/app/model_assets` contains model configuration and tokenizer files. The
self-contained runtime also includes the ready `rinalmo-imrnads` fine-tuned
checkpoint at `/models/best_model.pt`, SHA-256
`923d1e9f843d12fe29a9515a95b6a32eae286017790dfeef1c82215446958f0a`.
The Docker build verifies that digest and retains its identity and size in
`/usr/share/morse/sources/checkpoint-manifest.json`; the runtime verifies the
same digest before loading parameters. No model-fitting dataset is included.
The checkpoint and converted backbone retain their applicable upstream terms;
the first-party MIT notice does not replace them. Installed dependencies
retain their own notices and terms.
