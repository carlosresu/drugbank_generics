# _shared.R — Common setup for DrugBank R scripts
# Sourced by each script to handle package loading, parallel setup, and shared utilities.
# Guard prevents double-loading if scripts are sourced together.

# Guard to prevent double-loading
if (exists("DRUGBANK_SHARED_LOADED") && isTRUE(DRUGBANK_SHARED_LOADED)) {

  # Already loaded, skip
} else {
  DRUGBANK_SHARED_LOADED <- TRUE
  
  # ============================================================================
  # PACKAGE INSTALLATION & LOADING
  # ============================================================================
  
  suppressWarnings({
    suppressPackageStartupMessages({
      ensure_installed <- function(pkg) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          install.packages(pkg, repos = "https://cloud.r-project.org")
        }
      }
      ensure_dbdataset <- function() {
        if (requireNamespace("dbdataset", quietly = TRUE)) return(invisible(TRUE))
        ensure_installed("remotes")
        installer <- NULL
        if (requireNamespace("remotes", quietly = TRUE)) {
          installer <- remotes::install_github
        } else if (requireNamespace("devtools", quietly = TRUE)) {
          installer <- devtools::install_github
        }
        if (is.null(installer)) {
          stop("dbdataset package is required; install remotes or devtools to proceed.")
        }
        installer("interstellar-Consultation-Services/dbdataset", quiet = TRUE, upgrade = "never")
        invisible(TRUE)
      }
      ensure_installed("data.table")
      ensure_dbdataset()
    })
  })
  
  tryCatch({
    ensure_installed("arrow")
  }, error = function(e) {})
  
  library(data.table)
  library(dbdataset)
  
  # ============================================================================
  # GLOBAL CONFIGURATION
  # ============================================================================
  
  quiet_mode <- identical(tolower(Sys.getenv("ESOA_DRUGBANK_QUIET", "0")), "1")
  
  resolve_workers <- function() {
    env_val <- suppressWarnings(as.integer(Sys.getenv("ESOA_DRUGBANK_WORKERS", "")))
    if (!is.na(env_val) && env_val > 0) return(env_val)
    cores <- NA_integer_
    try(cores <- parallel::detectCores(logical = TRUE), silent = TRUE)
    if (is.na(cores) && requireNamespace("future", quietly = TRUE)) {
      try(cores <- future::availableCores(), silent = TRUE)
    }
    if (is.na(cores)) cores <- 1L
    # Default to min(8, cores) per AGENTS.md #6
    min(8L, max(1L, cores))
  }
  
  worker_count <- resolve_workers()
  data.table::setDTthreads(worker_count)
  
  # ============================================================================
  # PARALLEL SETUP (done once, reused by all scripts)
  # ============================================================================
  
  argv <- commandArgs(trailingOnly = TRUE)
  parallel_enabled <- !("--no-parallel" %in% argv)
  
  init_parallel_plan <- function(enabled_flag, workers) {
    plan_state <- NULL
    if (!enabled_flag) return(plan_state)
    tryCatch({
      ensure_installed("future")
      ensure_installed("future.apply")
      library(future)
      library(future.apply)
      os_name <- tolower(Sys.info()[["sysname"]])
      if (os_name %in% c("windows")) {
        plan(future::multisession, workers = workers)
      } else if (future::supportsMulticore()) {
        plan(future::multicore, workers = workers)
      } else {
        plan(future::multisession, workers = workers)
      }
      plan_state <- TRUE
    }, error = function(e) {
      parallel_enabled <<- FALSE
    })
    plan_state
  }
  
  plan_reset <- init_parallel_plan(parallel_enabled, worker_count)
  
  if (!quiet_mode) {
    cat(sprintf("[shared] data.table threads: %s | parallel workers: %s\n", 
                data.table::getDTthreads(), worker_count))
  }
  
  # ============================================================================
  # SHARED UTILITY FUNCTIONS
  # ============================================================================
  
  get_script_dir <- function() {
    # Try to get from --file= argument
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
    }
    # Fall back to sourced file location or working directory
    if (exists("DRUGBANK_BASE_DIR") && !is.null(DRUGBANK_BASE_DIR)) {
      return(DRUGBANK_BASE_DIR)
    }
    getwd()
  }
  
  # Set base directory (can be overridden by drugbank_all.R)
  if (!exists("DRUGBANK_BASE_DIR")) {
    DRUGBANK_BASE_DIR <- get_script_dir()
  }
  
  collapse_ws <- function(x) {
    ifelse(is.na(x), NA_character_, trimws(gsub("\\s+", " ", as.character(x))))
  }
  
  empty_to_na <- function(x) {
    val <- collapse_ws(x)
    ifelse(!nzchar(val), NA_character_, val)
  }
  
  unique_canonical <- function(values) {
    vals <- values[!is.na(values) & nzchar(values)]
    if (!length(vals)) return(character())
    vals[order(tolower(vals), vals)] |> unique()
  }
  
  collapse_pipe <- function(values) {
    vals <- unique_canonical(values)
    if (!length(vals)) return(NA_character_)
    paste(vals, collapse = "|")
  }
  
  write_csv_and_parquet <- function(dt, csv_path) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      data.table::fwrite(dt, csv_path)
      return(csv_path)
    }
    arrow_table <- arrow::Table$create(dt)
    if ("write_csv_arrow" %in% getNamespaceExports("arrow")) {
      arrow::write_csv_arrow(arrow_table, csv_path)
    } else if ("write_delim_arrow" %in% getNamespaceExports("arrow")) {
      arrow::write_delim_arrow(arrow_table, csv_path, delim = ",")
    } else {
      data.table::fwrite(dt, csv_path)
    }
    parquet_path <- sub("\\.csv$", ".parquet", csv_path)
    arrow::write_parquet(arrow_table, parquet_path)
    parquet_path
  }
  
  # Parallel lapply with chunking
  resolve_chunk_size <- function(total_rows) {
    env_val <- suppressWarnings(as.integer(Sys.getenv("ESOA_DRUGBANK_CHUNK", "")))
    if (!is.na(env_val) && env_val > 0) return(env_val)
    max(1000L, ceiling(total_rows / max(1L, worker_count * 3L)))
  }
  
  parallel_lapply <- function(x, fun) {
    n <- length(x)
    if (n == 0L) return(lapply(x, fun))
    if (!parallel_enabled || worker_count <= 1L) return(lapply(x, fun))
    if (!requireNamespace("future.apply", quietly = TRUE)) return(lapply(x, fun))
    chunk_size <- max(1L, min(resolve_chunk_size(n), ceiling(n / max(1L, worker_count))))
    idx <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
    chunks <- lapply(idx, function(ix) x[ix])
    worker <- function(chunk) lapply(chunk, fun)
    tryCatch({
      res <- future.apply::future_lapply(
        chunks,
        worker,
        future.seed = FALSE,
        future.scheduling = 1,
        future.chunk.size = 1
      )
      out <- vector("list", n)
      for (i in seq_along(idx)) out[idx[[i]]] <- res[[i]]
      out
    }, error = function(err) {
      lapply(x, fun)
    })
  }
  
  # Load drugbank dataset (cached)
  load_drugbank_dataset <- function() {
    if (exists("drugbank", inherits = FALSE)) {
      return(get("drugbank", inherits = FALSE))
    }
    if (exists("drugbank", where = "package:dbdataset")) {
      return(get("drugbank", envir = as.environment("package:dbdataset")))
    }
    data("drugbank", package = "dbdataset", envir = environment())
    if (!exists("drugbank", inherits = FALSE)) {
      stop("dbdataset::drugbank dataset is unavailable; reinstall dbdataset.")
    }
    get("drugbank", inherits = FALSE)
  }
  
  # Pre-load drugbank dataset so it's cached for all scripts
  drugbank <- load_drugbank_dataset()
}
