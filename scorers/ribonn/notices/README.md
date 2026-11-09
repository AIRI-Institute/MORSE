# Public runtime image notices

The MORSE adapter is provided under the accompanying MIT license (`LICENSE`).
The installed Python/PyTorch dependencies retain their own license notices.

This image includes the pinned Sanofi RiboNN inference extraction under
`/app/third_party/ribonn` and the selected 50-model ensemble under `/models`.
The exact Sanofi software notice (`LICENSE`) and `MODEL WEIGHTS LICENSE.txt`
remain with that source. They cover academic research and non-commercial use
by academic/non-profit organizations and retain their stated restrictions.
The adapter's separate MIT license does not relicense Sanofi materials.

`image-assets.json` records every included checkpoint/metadata hash and the
upstream inference-source/notice hashes. No model-fitting code or training
dataset is included.

The exact Dockerfile is retained at `/usr/share/morse/sources/Dockerfile`.
