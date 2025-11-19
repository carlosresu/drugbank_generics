#!/usr/bin/env Rscript

cat("Running DrugBank generics and mixtures exporters...\n")

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
  generics_master <- file.path(output_dir, "drugbank_generics_master.csv")
  if (file.exists(generics_master)) {
    dests <- c(
      file.path(output_dir, "drugbank_generics.csv"),
      generics_master
    )
    for (dest in dests) {
      safe_copy(generics_master, dest)
    }
  }

  mixtures_master <- file.path(output_dir, "drugbank_mixtures_master.csv")
  if (file.exists(mixtures_master)) {
    safe_copy(mixtures_master, file.path(output_dir, "drugbank_mixtures_master.csv"))
  }
}

run_subscript <- function(script_path, env_vars = character()) {
  if (!file.exists(script_path)) {
    stop(sprintf("Required script not found: %s", script_path))
  }
  status <- system2("Rscript", c(script_path), env = env_vars)
  if (!is.null(status) && status != 0) {
    stop(sprintf("Script %s exited with status %s", script_path, status))
  }
}

script_dir <- get_script_dir()
subscripts <- c("drugbank_generics.R", "drugbank_mixtures.R")
env_forward <- sprintf("ESOA_DRUGBANK_QUIET=%s", Sys.getenv("ESOA_DRUGBANK_QUIET", unset = "0"))

for (script in subscripts) {
  full_path <- file.path(script_dir, script)
  message(sprintf("[drugbank_all] Running %s", full_path))
  run_subscript(full_path, env_vars = env_forward)
}

propagate_outputs(script_dir)

cat("DrugBank exports completed successfully.\n")
