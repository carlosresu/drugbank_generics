#!/usr/bin/env Rscript
# Simple profiling of DrugBank R scripts using base R Rprof
# Usage: Rscript profile_simple.R [script_name]

args <- commandArgs(trailingOnly = TRUE)
script_name <- if (length(args) > 0) args[1] else "drugbank_generics"

script_dir <- tryCatch({
  file_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg[1])))
  else getwd()
}, error = function(e) getwd())

script_path <- file.path(script_dir, paste0(script_name, ".R"))
if (!file.exists(script_path)) {
  stop(sprintf("Script not found: %s", script_path))
}

output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
prof_file <- file.path(output_dir, paste0("Rprof_", script_name, ".out"))

cat(sprintf("Profiling: %s\n", script_path))
cat(sprintf("Profile output: %s\n\n", prof_file))

# Start profiling with 10ms interval
Rprof(prof_file, interval = 0.01, memory.profiling = FALSE)

# Source the script
tryCatch({
  source(script_path, local = new.env())
}, error = function(e) {
  cat(sprintf("Error: %s\n", e$message))
})

# Stop profiling
Rprof(NULL)

# Summarize
cat("\n=== PROFILE SUMMARY ===\n\n")
if (file.exists(prof_file) && file.info(prof_file)$size > 0) {
  summary <- summaryRprof(prof_file)
  
  cat("By self time (time spent in function itself):\n")
  by_self <- head(summary$by.self, 20)
  print(by_self)
  
  cat("\n\nBy total time (including calls to other functions):\n")
  by_total <- head(summary$by.total, 20)
  print(by_total)
  
  cat(sprintf("\n\nTotal profiling time: %.1f seconds\n", summary$sampling.time))
} else {
  cat("No profiling data collected (script may have run too fast)\n")
}
