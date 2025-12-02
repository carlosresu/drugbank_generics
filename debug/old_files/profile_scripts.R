#!/usr/bin/env Rscript
# Profile DrugBank R scripts to identify bottlenecks
# Usage: Rscript profile_scripts.R [script_name]
# Example: Rscript profile_scripts.R drugbank_generics

suppressPackageStartupMessages({
  if (!require("profvis", quietly = TRUE)) {
    install.packages("profvis", repos = "https://cloud.r-project.org")
  }
  library(profvis)
})

args <- commandArgs(trailingOnly = TRUE)
script_name <- if (length(args) > 0) args[1] else "drugbank_generics"

script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
if (length(script_dir) == 0) script_dir <- getwd()

script_path <- file.path(script_dir, paste0(script_name, ".R"))
if (!file.exists(script_path)) {
  stop(sprintf("Script not found: %s", script_path))
}

cat(sprintf("Profiling: %s\n", script_path))
cat("This may take a while...\n\n")

# Profile the script
prof_result <- profvis({
  source(script_path, local = new.env())
}, interval = 0.01)  # 10ms sampling

# Save HTML report
output_path <- file.path(script_dir, "output", paste0("profile_", script_name, ".html"))
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

htmlwidgets::saveWidget(prof_result, output_path, selfcontained = TRUE)
cat(sprintf("\nProfile saved to: %s\n", output_path))

# Print summary
cat("\n=== PROFILE SUMMARY ===\n")
prof_data <- prof_result$x$message$prof
if (!is.null(prof_data)) {
  # Aggregate by function
  funcs <- prof_data$ref
  times <- table(funcs)
  times_sorted <- sort(times, decreasing = TRUE)
  
  cat("\nTop 20 functions by sample count:\n")
  for (i in seq_len(min(20, length(times_sorted)))) {
    cat(sprintf("  %4d samples: %s\n", times_sorted[i], names(times_sorted)[i]))
  }
}
