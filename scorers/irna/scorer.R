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
  orf_max_log_kd = -1.5,
  panel_path = NULL,
  source_root = getwd()
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
  if (is.null(panel_path)) {
    panel_path <- file.path(panel_root, species, "core_975.tsv")
  }
  if (!file.exists(panel_path)) {
    stop("IRNA panel does not exist: ", panel_path)
  }
  panel_path <- normalizePath(panel_path)
  profile_md5 <- unname(tools::md5sum(panel_path))
  bundled_path <- normalizePath(file.path(source_root, "panels", "2026-08-21", species, "core_975.tsv"), mustWork = FALSE)
  panel_id <- if (identical(panel_path, bundled_path)) {
    paste(species, "core_975", "2026-08-21", sep = "-")
  } else {
    paste(species, basename(panel_path), profile_md5, sep = "-")
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
    panel_id = panel_id,
    panel_path = panel_path,
    profile_md5 = profile_md5,
    panel_root = normalizePath(panel_root),
    source_root = normalizePath(source_root),
    profiles = species_config$profiles,
    models = models,
    expression_weights = expression_weights,
    workers = workers,
    bp = bp,
    max_log_kd = max_log_kd,
    orf_max_log_kd = orf_max_log_kd,
    min_site_distance = 7L,
    canonical_only = TRUE,
    asset_cache = new.env(parent = emptyenv())
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
    onlyCanonical = scorer$canonical_only,
    maxLogKd = c(scorer$max_log_kd, scorer$orf_max_log_kd),
    minDist = scorer$min_site_distance,
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
  unknown <- setdiff(names(request), c("request_id", "sequences", "sequence_contexts", "parameters", "params"))
  if (length(unknown) > 0L) {
    stop("Unknown request fields: ", paste(sort(unknown), collapse = ", "))
  }
  if (is.null(request$sequences)) {
    stop("Request must define a list of sequences")
  }
  if (all(c("parameters", "params") %in% names(request))) {
    stop("Use only one of parameters and params")
  }
  parameters <- if ("parameters" %in% names(request)) {
    request$parameters
  } else if ("params" %in% names(request)) {
    request$params
  } else {
    list()
  }
  if (any(c("parameters", "params") %in% names(request)) && (!is.list(parameters) || is.null(names(parameters)))) {
    stop("parameters must be a JSON object")
  }
  effective <- resolve_irna_parameters(scorer, parameters)
  list(
    request_id = request$request_id,
    outputs = score_irna_batch(effective$scorer, request$sequences),
    parameters = effective$parameters
  )
}


resolve_irna_parameters <- function(scorer, parameters) {
  if (!is.list(parameters) || (length(parameters) > 0L && is.null(names(parameters)))) {
    stop("parameters must be a JSON object")
  }
  allowed <- c("profile_path", "species", "profiles", "max_log_kd", "min_site_distance", "canonical_only", "workers")
  unknown <- setdiff(names(parameters), allowed)
  if (length(unknown) > 0L) {
    stop("Unknown scanMiR parameters: ", paste(sort(unknown), collapse = ", "))
  }
  species <- if ("species" %in% names(parameters)) parameters$species else scorer$species
  if (!is.character(species) || length(species) != 1L || is.na(species) || !species %in% names(IRNA_SPECIES)) {
    stop("species must be human or mouse")
  }
  profile_path <- if ("profile_path" %in% names(parameters)) parameters$profile_path else NULL
  if (!is.null(profile_path)) {
    if (!is.character(profile_path) || length(profile_path) != 1L || is.na(profile_path) || !nzchar(profile_path)) {
      stop("profile_path must be a non-empty path")
    }
    candidates <- c(profile_path, file.path(scorer$source_root, profile_path), file.path(scorer$panel_root, profile_path))
    matches <- candidates[file.exists(candidates)]
    if (length(matches) == 0L) {
      stop("scanMiR profile file does not exist: ", profile_path)
    }
    profile_path <- normalizePath(matches[[1L]])
    inferred_species <- basename(dirname(profile_path))
    if (inferred_species %in% names(IRNA_SPECIES)) {
      if ("species" %in% names(parameters) && !identical(species, inferred_species)) {
        stop("species conflicts with profile_path")
      }
      species <- inferred_species
    }
  } else {
    if ("profile_path" %in% names(parameters)) {
      stop("profile_path must be a non-empty path")
    }
    profile_path <- if (identical(species, scorer$species)) scorer$panel_path else normalizePath(file.path(scorer$panel_root, species, "core_975.tsv"), mustWork = TRUE)
  }
  profile_md5 <- if (!"profile_path" %in% names(parameters) && identical(species, scorer$species)) {
    scorer$profile_md5
  } else {
    unname(tools::md5sum(profile_path))
  }
  positive_integer <- function(name, default, minimum = 1L) {
    value <- if (name %in% names(parameters)) parameters[[name]] else default
    if (!is.numeric(value) || length(value) != 1L || is.na(value) || !is.finite(value) || value != floor(value) || value < minimum || value > .Machine$integer.max) {
      stop(name, " must be an integer >= ", minimum)
    }
    as.integer(value)
  }
  workers <- positive_integer("workers", scorer$workers)
  min_site_distance <- positive_integer("min_site_distance", scorer$min_site_distance)
  max_log_kd <- if ("max_log_kd" %in% names(parameters)) parameters$max_log_kd else scorer$max_log_kd
  if (!is.numeric(max_log_kd) || length(max_log_kd) != 1L || is.na(max_log_kd) || !is.finite(max_log_kd)) {
    stop("max_log_kd must be a finite number")
  }
  canonical_only <- if ("canonical_only" %in% names(parameters)) parameters$canonical_only else scorer$canonical_only
  if (!is.logical(canonical_only) || length(canonical_only) != 1L || is.na(canonical_only)) {
    stop("canonical_only must be true or false")
  }
  profiles <- if ("profiles" %in% names(parameters)) parameters$profiles else NULL
  if (is.null(profiles)) {
    profiles <- IRNA_SPECIES[[species]]$profiles
  } else {
    if (!is.list(profiles) || !is.null(names(profiles)) || length(profiles) == 0L || any(!vapply(profiles, function(value) is.character(value) && length(value) == 1L && !is.na(value), logical(1L)))) {
      stop("profiles must be a non-empty list of profile names or null")
    }
    profiles <- vapply(profiles, identity, character(1L))
    if (anyDuplicated(profiles) || any(!profiles %in% IRNA_SPECIES[[species]]$profiles)) {
      stop("profiles must contain unique names from the selected species")
    }
  }
  # Cache only loaded species assets; every request builds its own settings view.
  effective <- scorer
  if (!identical(species, scorer$species) || !identical(profile_path, scorer$panel_path) || !identical(profile_md5, scorer$profile_md5)) {
    cache_key <- paste(species, profile_path, profile_md5, sep = "|")
    if (!exists(cache_key, envir = scorer$asset_cache, inherits = FALSE)) {
      loaded <- load_irna_scorer(
        species = species, panel_root = scorer$panel_root,
        workers = scorer$workers, max_log_kd = scorer$max_log_kd,
        orf_max_log_kd = scorer$orf_max_log_kd,
        panel_path = profile_path, source_root = scorer$source_root
      )
      assign(cache_key, loaded, envir = scorer$asset_cache)
    }
    effective <- get(cache_key, envir = scorer$asset_cache, inherits = FALSE)
  }
  effective$profiles <- profiles
  effective$expression_weights <- effective$expression_weights[, profiles, drop = FALSE]
  effective$max_log_kd <- as.numeric(max_log_kd)
  effective$min_site_distance <- min_site_distance
  effective$canonical_only <- canonical_only
  if (workers != effective$workers) {
    effective$bp <- if (workers == 1L) SerialParam(progressbar = FALSE) else MulticoreParam(workers, progressbar = FALSE)
    effective$workers <- workers
  }
  list(
    scorer = effective,
    parameters = list(
      profile_path = profile_path, profile_md5 = profile_md5,
      species = species, profiles = unname(as.list(profiles)),
      max_log_kd = as.numeric(max_log_kd), min_site_distance = min_site_distance,
      canonical_only = canonical_only, workers = workers
    )
  )
}
