# MORSE interfering-RNA runtime notices

The MORSE scorer adapter (`/app/scorer.R` and `/app/server.R`) is distributed
under the MIT license reproduced in `LICENSE` in this directory.

## scanMiR and scanMiRData

This runtime includes scanMiR 1.16.0 and scanMiRData 1.16.0 by Pierre-Luc
Germain, Michael Soutschek, and Fridolin Gross. Their package DESCRIPTION
files declare `GPL-3`. The full GNU General Public License, version 3, is
reproduced in `GPL-3.0.txt`. The packages are provided without warranty as
specified by that license. Recipients may modify and redistribute them under
its terms.

The unmodified source package archives used to install these exact versions
are included in `source/scanMiR_1.16.0.tar.gz` and
`source/scanMiRData_1.16.0.tar.gz`, alongside the source archives downloaded
for additional R dependencies during the build. `source/SHA256SUMS` records
their SHA-256 checksums. Each archive retains the upstream source, package
metadata, documentation, data, and notices. `Dockerfile` records the build
and installation procedure; `installed-packages.tsv` records the installed
R package versions and declared licenses. The source archives can be
extracted from the image without starting the scorer.

Upstream locations:

- https://bioconductor.org/packages/3.22/bioc/html/scanMiR.html
- https://bioconductor.org/packages/3.22/data/experiment/html/scanMiRData.html
- https://git.bioconductor.org/packages/scanMiR
- https://git.bioconductor.org/packages/scanMiRData

## Expression panels

The complete bundled panel snapshot is `/app/panels/2026-08-21`. It contains
derived human and mouse miRNA panels assembled from scanMiRData 1.16.0,
MirGeneDB 3.0, and miRNATissueAtlas 2025. The original identifiers, source
URLs, checksums, matching audit flags, and expression profiles are retained
in the panel TSVs and `sources.json`; `coverage.tsv` and `summary.json`
retain the coverage audit. The MIT notice for the MORSE adapter does not
relicense these third-party sources.

- MirGeneDB 3.0: https://www.mirgenedb.org/
- miRNATissueAtlas 2025: https://web.ccb.uni-saarland.de/mirnatissueatlas_2025/

## Base runtime and other dependencies

The image also contains R, Bioconductor, and their dependencies from the
pinned Bioconductor base image. Their existing license and copyright
notices are retained in their installed package and system documentation
locations. The MIT license for MORSE does not replace any third-party
license.
