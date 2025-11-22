#!/usr/bin/env Rscript
# Legacy entrypoint kept for compatibility; delegates to polars-based pipeline.

message("drugbank.R is deprecated; running drugbank_all.R (polars-first).")

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
}

script_dir <- get_script_dir()
delegate <- file.path(script_dir, "drugbank_all.R")
if (!file.exists(delegate)) {
  stop("drugbank_all.R not found; cannot continue.")
}

status <- system2("Rscript", delegate)
if (!is.null(status) && status != 0) {
  stop(sprintf("drugbank_all.R exited with status %s", status))
}
