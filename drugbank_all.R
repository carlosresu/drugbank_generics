#!/usr/bin/env Rscript

cat("Running DrugBank generics and mixtures exporters...\n")

ensure_packages <- function() {
  if (!is.null(opt_lib) && dir.exists(opt_lib)) {
    .libPaths(c(opt_lib, .libPaths()))
  }
  options(
    repos = c(CRAN = "https://cloud.r-project.org"),
    R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true"
  )
  base_pkgs <- c(
    "arrow", "data.table", "future", "future.apply", "httr2", "memoise",
    "pacman", "purrr", "readr", "rvest", "stringr", "tibble", "xml2", "remotes"
  )
  missing <- base_pkgs[!vapply(base_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    install.packages(missing, repos = "https://cloud.r-project.org", lib = opt_lib)
  }
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", lib = opt_lib)
  }
  if (!requireNamespace("dbdataset", quietly = TRUE)) {
    remotes::install_github(
      "interstellar-Consultation-Services/dbdataset",
      quiet = TRUE,
      upgrade = "never",
      lib = opt_lib
    )
  }
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  getwd()
}

find_rscript <- function() {
  env_path <- Sys.getenv("RSCRIPT_PATH", unset = "")
  if (nzchar(env_path) && file.exists(env_path)) {
    return(normalizePath(env_path, winslash = "/", mustWork = FALSE))
  }
  which_path <- Sys.which("Rscript")
  if (nzchar(which_path)) {
    return(normalizePath(which_path, winslash = "/", mustWork = FALSE))
  }
  candidates <- c(
    file.path(R.home("bin"), "Rscript"),
    file.path(R.home("bin"), "Rscript.exe"),
    file.path(R.home("bin"), "x64", "Rscript.exe")
  )
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
    }
  }
  stop("Rscript executable not found; set RSCRIPT_PATH or add Rscript to PATH.")
}

paths_equal <- function(a, b) {
  normalizePath(a, mustWork = FALSE) == normalizePath(b, mustWork = FALSE)
}

safe_copy <- function(src, dest) {
  if (!file.exists(src)) return(invisible(FALSE))
  if (paths_equal(src, dest)) return(invisible(TRUE))
  dest_dir <- dirname(dest)
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(src, dest, overwrite = TRUE, copy.mode = TRUE)
  invisible(TRUE)
}

propagate_outputs <- function(script_dir) {
  output_dir <- file.path(script_dir, "output")
  mixtures_master <- file.path(output_dir, "drugbank_mixtures_master.csv")
  if (file.exists(mixtures_master)) {
    safe_copy(mixtures_master, file.path(output_dir, "drugbank_mixtures_master.csv"))
  }
}

run_subscript <- function(script_path, env_vars = character()) {
  if (!file.exists(script_path)) {
    stop(sprintf("Required script not found: %s", script_path))
  }
  env_entries <- c(
    env_vars,
    sprintf("RSCRIPT_PATH=%s", rscript_bin),
    sprintf("R_LIBS_USER=%s", opt_lib),
    sprintf("R_LIBS=%s", lib_env)
  )
  entry_names <- names(env_entries)
  if (is.null(entry_names)) entry_names <- rep("", length(env_entries))
  split_entries <- strsplit(env_entries, "=", fixed = TRUE)
  fallback_names <- vapply(split_entries, function(x) if (length(x)) x[[1]] else "", character(1))
  values <- vapply(
    seq_along(split_entries),
    function(i) {
      parts <- split_entries[[i]]
      if (length(parts) <= 1) return(env_entries[[i]])
      paste(parts[-1], collapse = "=")
    },
    character(1)
  )
  final_names <- ifelse(nzchar(entry_names), entry_names, fallback_names)
  valid <- !is.na(final_names) & nzchar(final_names)
  env_all <- structure(values[valid], names = final_names[valid])
  env_names <- names(env_all)
  if (is.null(env_names)) env_names <- character(length(env_all))
  # Apply env vars for the subprocess without polluting the parent session.
  prior <- if (length(env_names)) Sys.getenv(env_names, unset = NA_character_) else character()
  on.exit({
    if (!length(env_names)) return(invisible())
    for (i in seq_along(prior)) {
      nm <- env_names[i]
      if (is.na(nm) || !nzchar(nm)) next
      if (is.na(prior[i])) {
        Sys.unsetenv(nm)
      } else {
        do.call(Sys.setenv, stats::setNames(list(prior[i]), nm))
      }
    }
  }, add = TRUE)
  if (length(env_all)) do.call(Sys.setenv, as.list(env_all))

  tmp_log <- tempfile("drugbank_subscript_", fileext = ".log")
  status <- system2(rscript_bin, c(script_path), stdout = tmp_log, stderr = tmp_log)
  if (!is.null(status) && status != 0) {
    output <- tryCatch(paste(readLines(tmp_log, warn = FALSE), collapse = "\n"), error = function(...) "")
    dest_log <- file.path(script_dir, "output", "drugbank_all_error.log")
    dir.create(dirname(dest_log), recursive = TRUE, showWarnings = FALSE)
    try(file.copy(tmp_log, dest_log, overwrite = TRUE), silent = TRUE)
    if (nzchar(output)) {
      message("\n--- Subscript output start ---\n", output, "\n--- Subscript output end ---\n")
    }
    if (file.exists(dest_log)) {
      message(sprintf("Subscript log saved to %s", dest_log))
    }
    stop(sprintf("Script %s exited with status %s", script_path, status))
  }
  unlink(tmp_log, force = TRUE)
}

script_dir <- get_script_dir()
opt_lib <- normalizePath(file.path(script_dir, "r_libs"), winslash = "/", mustWork = FALSE)
if (!dir.exists(opt_lib)) dir.create(opt_lib, recursive = TRUE, showWarnings = FALSE)
rscript_bin <- find_rscript()
orig_libs <- .libPaths()
lib_env <- paste(unique(c(opt_lib, orig_libs)), collapse = .Platform$path.sep)
Sys.setenv(R_LIBS_USER = opt_lib, R_LIBS = lib_env)
ensure_packages()
workers <- Sys.getenv("ESOA_DRUGBANK_WORKERS", unset = "")
if (nzchar(workers)) {
  message(sprintf("[drugbank_all] Using %s workers (ESOA_DRUGBANK_WORKERS)", workers))
}
subscripts <- c("drugbank_generics.R", "drugbank_mixtures.R")
env_forward <- sprintf("ESOA_DRUGBANK_QUIET=%s", Sys.getenv("ESOA_DRUGBANK_QUIET", unset = "0"))

for (script in subscripts) {
  full_path <- file.path(script_dir, script)
  message(sprintf("[drugbank_all] Running %s", full_path))
  run_subscript(full_path, env_vars = env_forward)
}

propagate_outputs(script_dir)

cat("DrugBank exports completed successfully.\n")
