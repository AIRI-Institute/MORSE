# Self-contained inference image notices

The MORSE inference adapter and ready task checkpoint are supplied with the
accompanying first-party MIT notice (`LICENSE`). The image also contains the
frozen Orthrus backbone configuration, inference source and pretrained weights.
Its exact upstream MIT notice is `ORTHRUS-LICENSE.md`; revision, attribution and
file identities are recorded in `ORTHRUS-PROVENANCE.md`.

The image contains ready-model inference assets, without fitting code or fitting
datasets. Installed PyTorch, Mamba, causal-conv1d, Transformers, NumPy and other
dependencies retain their own original license notices.

The exact Dockerfile, including checkpoint/backbone SHA-256 checks, is retained
at `/usr/share/morse/sources/Dockerfile`. All notices above are retained under
`/usr/share/morse/notices`.
