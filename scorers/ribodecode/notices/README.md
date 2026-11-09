# Public runtime image notices

The MORSE adapter is provided under the accompanying MIT license (`LICENSE`).
The installed Python/PyTorch dependencies retain their own license notices.

This image includes the exact author TranslationModel 1.1.0 package at
`/opt/ribodecode`, including its ready checkpoint and built-in condition data,
plus the official environment CSV at `/assets/env_file.csv`. It retains
SSPL-1.0, copyright 2024 Xie Lab. The original notice is preserved both in the
installed dist-info and as `TranslationModel-SSPL-1.0.txt` here. MORSE's separate
MIT notice does not relicense the upstream program or a combined covered work.

The author-supplied preferred Python/configuration/condition and packaging
source files are also retained under `/usr/share/morse/sources/TranslationModel`.
That directory includes an index of installed source/weight locations and
their exact hashes. The native tensor checkpoint remains in the runnable
package; it is not duplicated in the source index. Adapter source is in
`/app/server.py`. Preserve the source materials and notices under the applicable
SSPL conveying terms. A hosted service offering the program's functionality
must additionally comply with section 13's Service Source Code requirements.

The exact Dockerfile is retained at `/usr/share/morse/sources/Dockerfile`.
