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

run_subscript <- function(script_path) {
  if (!file.exists(script_path)) {
    stop(sprintf("Required script not found: %s", script_path))
  }
  status <- system2("Rscript", c(script_path))
  if (!is.null(status) && status != 0) {
    stop(sprintf("Script %s exited with status %s", script_path, status))
  }
}

script_dir <- get_script_dir()
subscripts <- c("drugbank_generics.R", "drugbank_mixtures.R")

for (script in subscripts) {
  full_path <- file.path(script_dir, script)
  message(sprintf("[drugbank_all] Running %s", full_path))
  run_subscript(full_path)
}

cat("DrugBank exports completed successfully.\n")
