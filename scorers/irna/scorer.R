suppressPackageStartupMessages({
  library(BiocParallel)
  library(jsonlite)
  library(scanMiR)
})


IRNA_SPECIES <- list(
  human = list(
    collection = "hsa",
    profiles = c(
      "pbmc",
      "dendritic_cells",
      "macrophage",
      "mononuclear_cells",
      "limb_muscle"
    )
  ),
  mouse = list(
    collection = "mmu",
    profiles = c(
      "pbmc",
      "macrophage",
      "inflammatory_monocyte",
      "limb_muscle",
      "spleen",
      "lymph_node"
    )
  )
)


load_irna_scorer <- function(
  species = "human",
  panel_root = file.path("panels", "2026-08-21"),
  workers = 4L,
  max_log_kd = -1.0,
  orf_max_log_kd = -1.5
) {
  species <- tolower(trimws(species))
  if (!species %in% names(IRNA_SPECIES)) {
    stop("IRNA species must be one of: ", paste(names(IRNA_SPECIES), collapse = ", "))
  }
  workers <- as.integer(workers)
  if (is.na(workers) || workers < 1L) {
    stop("IRNA workers must be a positive integer")
  }
  max_log_kd <- as.numeric(max_log_kd)
  orf_max_log_kd <- as.numeric(orf_max_log_kd)
  if (!is.finite(max_log_kd) || !is.finite(orf_max_log_kd)) {
    stop("IRNA maxLogKd thresholds must be finite")
  }

  species_config <- IRNA_SPECIES[[species]]
  panel_path <- file.path(panel_root, species, "core_975.tsv")
  if (!file.exists(panel_path)) {
    stop("IRNA panel does not exist: ", panel_path)
  }
  panel <- read.delim(
    panel_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  required_columns <- c(
    "scanmir_name",
    paste0("fraction_total__", species_config$profiles)
  )
  missing_columns <- setdiff(required_columns, colnames(panel))
  if (length(missing_columns) > 0L) {
    stop("IRNA panel is missing columns: ", paste(missing_columns, collapse = ", "))
  }
  if (anyDuplicated(panel$scanmir_name)) {
    stop("IRNA panel contains duplicate scanMiR model names")
  }

  data(
    list = species_config$collection,
    package = "scanMiRData",
    envir = environment()
  )
  collection <- get(species_config$collection, envir = environment())
  missing_models <- setdiff(panel$scanmir_name, names(collection))
  if (length(missing_models) > 0L) {
    stop(
      "IRNA panel references missing scanMiR models: ",
      paste(head(missing_models, 10L), collapse = ", ")
    )
  }
  models <- collection[panel$scanmir_name]
  if (!identical(names(models), panel$scanmir_name)) {
    stop("IRNA scanMiR model order does not match the panel")
  }

  fraction_columns <- paste0("fraction_total__", species_config$profiles)
  expression_weights <- as.matrix(panel[, fraction_columns, drop = FALSE])
  storage.mode(expression_weights) <- "double"
  rownames(expression_weights) <- panel$scanmir_name
  colnames(expression_weights) <- species_config$profiles
  if (any(!is.finite(expression_weights)) || any(expression_weights < 0)) {
    stop("IRNA panel expression weights must be finite and non-negative")
  }

  bp <- if (workers == 1L) {
    SerialParam(progressbar = FALSE)
  } else {
    MulticoreParam(workers, progressbar = FALSE)
  }
  list(
    species = species,
    panel_id = paste(species, "core_975", "2026-08-21", sep = "-"),
    panel_path = normalizePath(panel_path),
    profiles = species_config$profiles,
    models = models,
    expression_weights = expression_weights,
    workers = workers,
    bp = bp,
    max_log_kd = max_log_kd,
    orf_max_log_kd = orf_max_log_kd
  )
}

normalize_irna_sequence <- function(sequence) {
  if (!is.character(sequence) || length(sequence) != 1L || is.na(sequence)) {
    stop("Every IRNA input sequence must be a string")
  }
  normalized <- chartr("U", "T", toupper(gsub("[[:space:]]+", "", sequence)))
  if (!nzchar(normalized)) {
    stop("IRNA input sequence must not be empty")
  }
  if (grepl("[^ACGT]", normalized)) {
    invalid <- unique(strsplit(gsub("[ACGT]", "", normalized), "")[[1L]])
    stop("IRNA input sequence contains unsupported bases: ", paste(invalid, collapse = ", "))
  }
  normalized
}


score_irna_batch <- function(scorer, sequences) {
  if (!is.list(sequences)) {
    stop("Request must define sequences as a JSON list")
  }
  if (length(sequences) == 0L) {
    return(list())
  }
  normalized <- vapply(sequences, normalize_irna_sequence, character(1L))
  candidate_ids <- sprintf("candidate_%06d", seq_along(normalized))
  names(normalized) <- candidate_ids

  hits <- findSeedMatches(
    normalized,
    scorer$models,
    onlyCanonical = TRUE,
    maxLogKd = c(scorer$max_log_kd, scorer$orf_max_log_kd),
    minDist = 7L,
    keepMatchSeq = FALSE,
    p3.extra = FALSE,
    ret = "data.frame",
    BP = scorer$bp,
    verbose = FALSE
  )
  if (nrow(hits) > 0L) {
    hits$miRNA <- as.character(hits$miRNA)
    hits$transcript <- as.character(hits$transcript)
    aggregated <- aggregateMatches(
      hits,
      keepSiteInfo = TRUE,
      toInt = FALSE,
      BP = scorer$bp
    )
    aggregated$miRNA <- as.character(aggregated$miRNA)
    aggregated$transcript <- as.character(aggregated$transcript)
  } else {
    aggregated <- data.frame()
  }

  lapply(candidate_ids, function(candidate_id) {
    candidate_hits <- if (nrow(hits) == 0L) {
      hits
    } else {
      hits[hits$transcript == candidate_id, , drop = FALSE]
    }
    candidate_aggregated <- if (nrow(aggregated) == 0L) {
      aggregated
    } else {
      aggregated[aggregated$transcript == candidate_id, , drop = FALSE]
    }
    score_irna_candidate(scorer, candidate_hits, candidate_aggregated)
  })
}


score_irna_candidate <- function(scorer, hits, aggregated) {
  profile_risks <- setNames(rep(0.0, length(scorer$profiles)), scorer$profiles)
  top_mirna <- ""

  if (nrow(aggregated) > 0L) {
    model_rows <- match(aggregated$miRNA, rownames(scorer$expression_weights))
    if (anyNA(model_rows)) {
      stop("scanMiR returned a model absent from the IRNA expression panel")
    }
    repression_risk <- pmax(0.0, -as.numeric(aggregated$repression))
    profile_weights <- scorer$expression_weights[model_rows, , drop = FALSE]
    profile_risks <- colSums(profile_weights * repression_risk)
    max_profile_index <- which.max(profile_risks)
    contributions <- repression_risk * profile_weights[, max_profile_index]
    if (length(contributions) > 0L && max(contributions) > 0.0) {
      top_mirna <- aggregated$miRNA[[which.max(contributions)]]
    }
  }

  max_profile_index <- which.max(profile_risks)
  max_profile <- names(profile_risks)[[max_profile_index]]
  strongest_log_kd <- if (nrow(hits) == 0L) {
    0.0
  } else {
    min(as.numeric(hits$log_kd)) / 1000.0
  }
  site_count <- function(name) {
    if (nrow(aggregated) == 0L || !name %in% colnames(aggregated)) {
      return(0L)
    }
    as.integer(sum(aggregated[[name]]))
  }

  list(
    irna_max_profile_risk = as.numeric(profile_risks[[max_profile_index]]),
    irna_mean_profile_risk = as.numeric(mean(profile_risks)),
    irna_site_count = as.integer(nrow(hits)),
    irna_8mer_site_count = site_count("8mer"),
    irna_7mer_site_count = site_count("7mer"),
    irna_6mer_site_count = site_count("6mer"),
    irna_strongest_log_kd = as.numeric(strongest_log_kd),
    irna_max_profile = max_profile,
    irna_top_mirna = top_mirna,
    irna_species = scorer$species,
    irna_panel_id = scorer$panel_id,
    irna_model_count = as.integer(length(scorer$models))
  )
}


score_irna_request <- function(scorer, request) {
  if (!is.list(request) || is.null(names(request))) {
    stop("Request must be a JSON object")
  }
  unknown <- setdiff(names(request), c("request_id", "sequences", "sequence_contexts"))
  if (length(unknown) > 0L) {
    stop("Unknown request fields: ", paste(sort(unknown), collapse = ", "))
  }
  if (is.null(request$sequences)) {
    stop("Request must define a list of sequences")
  }
  list(
    request_id = request$request_id,
    outputs = score_irna_batch(scorer, request$sequences)
  )
}
