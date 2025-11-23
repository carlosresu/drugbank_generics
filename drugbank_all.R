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
  env_all <- c(env_vars, sprintf("RSCRIPT_PATH=%s", rscript_bin), sprintf("R_LIBS_USER=%s", opt_lib), sprintf("R_LIBS=%s", lib_env))
  # Apply env vars for the subprocess without polluting the parent session.
  names(env_all) <- sub("=.*$", "", env_all)
  vals <- sub("^[^=]+=", "", env_all)
  prior <- Sys.getenv(names(env_all), unset = NA_character_)
  on.exit({for (i in seq_along(prior)) { if (is.na(prior[i])) Sys.unsetenv(names(env_all)[i]) else Sys.setenv(structure(prior[i], names = names(env_all)[i])) }}, add = TRUE)
  do.call(Sys.setenv, as.list(structure(vals, names = names(env_all))))

  tmp_log <- tempfile("drugbank_subscript_", fileext = ".log")
  status <- system2(rscript_bin, c(script_path), stdout = tmp_log, stderr = tmp_log)
  if (!is.null(status) && status != 0) {
    output <- tryCatch(paste(readLines(tmp_log, warn = FALSE), collapse = "\n"), error = function(...) "")
    if (nzchar(output)) {
      message("\n--- Subscript output start ---\n", output, "\n--- Subscript output end ---\n")
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
