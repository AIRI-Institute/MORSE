ARG BASE_IMAGE=bioconductor/bioconductor_docker:RELEASE_3_22@sha256:e34177896c1579f66d48b756e4229b58237499183a089da084af72d7aba1d40a
FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.source="https://github.com/AIRI-Institute/MORSE" \
      org.opencontainers.image.version="1.0.0"

# Keep the source packages used by this build with the distributed runtime.
RUN mkdir -p /usr/share/licenses/morse-irna/source \
    && Rscript -e 'options(pkgType="source", timeout=600L); source_dir <- "/usr/share/licenses/morse-irna/source"; repositories <- c(CRAN="https://cloud.r-project.org", BioCsoft="https://bioconductor.org/packages/3.22/bioc", BioCexperiment="https://bioconductor.org/packages/3.22/data/experiment"); required <- c("scanMiR", "scanMiRData", "jsonlite"); archives <- file.path(source_dir, c("scanMiR_1.16.0.tar.gz", "scanMiRData_1.16.0.tar.gz")); available <- available.packages(repos=repositories, type="source"); stopifnot(all(required %in% rownames(available)), all(available[c("scanMiR", "scanMiRData"), "Version"] == "1.16.0")); for (attempt in seq_len(5L)) { pending <- required[!vapply(required, requireNamespace, logical(1L), quietly=TRUE) | c(!file.exists(archives), FALSE)]; if (length(pending)) try(install.packages(pending, repos=repositories, available=available, dependencies=c("Depends", "Imports", "LinkingTo"), Ncpus=4L, type="source", destdir=source_dir), silent=FALSE); if (all(vapply(required, requireNamespace, logical(1L), quietly=TRUE)) && all(file.exists(archives))) break; if (attempt < 5L) message("Retrying incomplete R package installation (attempt ", attempt + 1L, "/5)") }; stopifnot(all(file.exists(archives)), system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", archives)) == 0L, as.character(packageVersion("scanMiR")) == "1.16.0", as.character(packageVersion("scanMiRData")) == "1.16.0"); write.table(installed.packages()[, c("Package", "Version", "License"), drop=FALSE], file="/usr/share/licenses/morse-irna/installed-packages.tsv", sep="\t", row.names=FALSE, quote=TRUE)' \
    && cd /usr/share/licenses/morse-irna/source \
    && printf '%s\n' \
      '28bf8ad4b0bf1cecfd2a17d6278d2896b2786629a44b48a4e2a33a4ca376df79  scanMiR_1.16.0.tar.gz' \
      '268b0ceb69b2220266103032746e3dfe8dc219914afc2586354195c3ee021a96  scanMiRData_1.16.0.tar.gz' \
      | sha256sum -c - \
    && sha256sum *.tar.gz > SHA256SUMS

COPY LICENSE THIRD_PARTY_NOTICES.md /usr/share/licenses/morse-irna/
COPY licenses/GPL-3.0.txt /usr/share/licenses/morse-irna/GPL-3.0.txt
COPY Dockerfile /usr/share/licenses/morse-irna/Dockerfile

ENV IRNA_SPECIES=human \
    IRNA_WORKERS=4 \
    IRNA_PANEL_ROOT=/app/panels/2026-08-21 \
    IRNA_MAX_LOG_KD=-1.0 \
    IRNA_ORF_MAX_LOG_KD=-1.5

WORKDIR /app
COPY scorer.R server.R /app/
COPY panels /app/panels

ENTRYPOINT ["Rscript", "/app/server.R"]
