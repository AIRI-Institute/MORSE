args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}
source(file.path(script_dir, "scorer.R"), local = FALSE)


integer_environment <- function(name, default) {
  value <- Sys.getenv(name, unset = default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed)) {
    stop(name, " must be an integer")
  }
  parsed
}


numeric_environment <- function(name, default) {
  value <- Sys.getenv(name, unset = default)
  parsed <- suppressWarnings(as.numeric(value))
  if (!is.finite(parsed)) {
    stop(name, " must be a finite number")
  }
  parsed
}


scorer <- load_irna_scorer(
  species = Sys.getenv("IRNA_SPECIES", unset = "human"),
  panel_root = Sys.getenv(
    "IRNA_PANEL_ROOT",
    unset = file.path(script_dir, "panels", "2026-08-21")
  ),
  workers = integer_environment("IRNA_WORKERS", "4"),
  max_log_kd = numeric_environment("IRNA_MAX_LOG_KD", "-1.0"),
  orf_max_log_kd = numeric_environment("IRNA_ORF_MAX_LOG_KD", "-1.5")
)

message(
  "IRNA scorer ready species=", scorer$species,
  " panel=", scorer$panel_id,
  " models=", length(scorer$models),
  " profiles=", length(scorer$profiles),
  " workers=", scorer$workers
)

input <- file("stdin", open = "r")
on.exit(close(input), add = TRUE)
repeat {
  line <- readLines(input, n = 1L, warn = FALSE)
  if (length(line) == 0L) {
    break
  }
  if (!nzchar(trimws(line))) {
    next
  }
  request_id <- NULL
  response <- tryCatch(
    {
      request <- fromJSON(line, simplifyVector = FALSE)
      request_id <- request$request_id
      score_irna_request(scorer, request)
    },
    error = function(error) {
      message("IRNA request failed: ", conditionMessage(error))
      list(request_id = request_id, error = conditionMessage(error))
    }
  )
  cat(
    toJSON(response, auto_unbox = TRUE, null = "null", digits = NA),
    "\n",
    sep = ""
  )
  flush.console()
}
