#!/usr/bin/env Rscript
# drugbank_mixtures.R — build a mixtures-focused DrugBank master dataset

suppressWarnings({
  suppressPackageStartupMessages({
    ensure_installed <- function(pkg) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }
    ensure_installed("data.table")
    ensure_installed("dbdataset")
  })
})

tryCatch({
  ensure_installed("arrow")
}, error = function(e) {})

library(data.table)
library(dbdataset)

argv <- commandArgs(trailingOnly = TRUE)
parallel_enabled <- !("--no-parallel" %in% argv)
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
  max(1L, cores)
}

detect_os_name <- function() {
  os <- Sys.info()[["sysname"]]
  if (is.null(os)) os <- .Platform$OS.type
  tolower(os)
}

worker_count <- resolve_workers()
data.table::setDTthreads(worker_count)

resolve_chunk_size <- function(total_rows) {
  env_val <- suppressWarnings(as.integer(Sys.getenv("ESOA_DRUGBANK_CHUNK", "")))
  if (!is.na(env_val) && env_val > 0) return(env_val)
  # Aim for multiple chunks per worker while avoiding excessive overhead.
  max(1000L, ceiling(total_rows / max(1L, worker_count * 3L)))
}

init_parallel_plan <- function(enabled_flag, workers) {
  plan_state <- NULL
  if (!enabled_flag) return(plan_state)
  tryCatch({
    ensure_installed("future")
    ensure_installed("future.apply")
    library(future)
    library(future.apply)
    os_name <- detect_os_name()
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
cat(sprintf("[drugbank_mixtures] data.table threads: %s | parallel workers: %s\n", data.table::getDTthreads(), worker_count))

choose_backend <- function() {
  backend <- tolower(Sys.getenv("ESOA_DRUGBANK_BACKEND", "future"))
  if (backend == "auto") {
    if (.Platform$OS.type == "unix") return("mclapply")
    return("future")
  }
  if (backend %chin% c("mclapply", "future")) return(backend)
  "future"
}

parallel_lapply <- function(x, fun) {
  backend <- choose_backend()
  if (parallel_enabled && backend == "mclapply" && .Platform$OS.type == "unix" && requireNamespace("parallel", quietly = TRUE)) {
    cores <- max(1L, worker_count)
    res <- tryCatch(
      parallel::mclapply(
        x,
        fun,
        mc.cores = cores,
        mc.preschedule = FALSE
      ),
      error = function(err) {
        msg <- sprintf(
          "[parallel:mclapply] failed: %s | length(x)=%s | workers=%s",
          conditionMessage(err),
          length(x),
          cores
        )
        cat(msg, "\n")
        NULL
      }
    )
    if (!is.null(res)) return(res)
    cat("[parallel] Falling back to future backend after mclapply failure\n")
  }
  if (parallel_enabled) {
    return(
      tryCatch(
        future.apply::future_lapply(
          x,
          fun,
          future.seed = FALSE,
          future.scheduling = worker_count,
          future.chunk.size = NULL
        ),
        error = function(err) {
          msg <- sprintf(
            "[parallel:future] failed: %s | length(x)=%s | workers=%s",
            conditionMessage(err),
            length(x),
            worker_count
          )
          cat(msg, "\n")
          stop(msg, call. = FALSE)
        }
      )
    )
  }
  stop("[parallel] No supported parallel backend available on this platform.")
}

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmd_args)
  if (length(match) > 0) {
    return(normalizePath(dirname(sub(needle, "", cmd_args[match[1]]))))
  }
  normalizePath(getwd())
}

paths_equal <- function(a, b) {
  normalizePath(a, mustWork = FALSE) == normalizePath(b, mustWork = FALSE)
}

safe_copy <- function(src, dest) {
  if (is.null(src) || is.null(dest) || is.na(src) || is.na(dest)) {
    stop("safe_copy: invalid path")
  }
  if (paths_equal(src, dest)) return(invisible(dest))
  dest_dir <- dirname(dest)
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(src, dest, overwrite = TRUE, copy.mode = TRUE)
  invisible(dest)
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

combine_values <- function(...) {
  inputs <- list(...)
  combined <- unlist(inputs, use.names = FALSE)
  unique_canonical(combined)
}

split_ingredients <- function(value) {
  val <- collapse_ws(value)
  if (is.na(val) || !nzchar(val)) return(character())
  repeat {
    new_val <- gsub("(?<!\\S)\\([^()]*\\)(?=\\s|$)", " ", val, perl = TRUE)
    if (identical(new_val, val)) break
    val <- new_val
  }
  val <- gsub("\\s+", " ", val, perl = TRUE)
  if (grepl("\\+", val)) {
    val <- gsub(",\\s[^+]+(?=\\s*\\+)", "", val, perl = TRUE)
  }
  val <- gsub(",\\s[^+]+$", "", val, perl = TRUE)
  parts <- if (grepl("\\+", val)) {
    unlist(strsplit(val, "\\+", perl = TRUE), use.names = FALSE)
  } else {
    val
  }
  parts <- collapse_ws(parts)
  parts <- parts[nzchar(parts)]
  unique_canonical(parts)
}

split_raw_components <- function(value) {
  val <- collapse_ws(value)
  if (is.na(val) || !nzchar(val)) return(character())
  parts <- if (grepl("\\+", val)) {
    unlist(strsplit(val, "\\+", perl = TRUE), use.names = FALSE)
  } else {
    val
  }
  parts <- collapse_ws(parts)
  parts[nzchar(parts)]
}

collapse_pipe <- function(values) {
  vals <- unique_canonical(values)
  if (!length(vals)) return(NA_character_)
  paste(vals, collapse = "|")
}

SALT_SYNONYM_LOOKUP <- list(
  "hydrochloride" = c("hydrochloride", "hydrochlorid", "hcl"),
  "sodium" = c("sodium", "na"),
  "potassium" = c("potassium", "k"),
  "calcium" = c("calcium", "ca"),
  "sulfate" = c("sulfate", "sulphate"),
  "sulphate" = c("sulphate", "sulfate")
)

expand_salt_set <- function(values) {
  expanded <- values
  for (val in values) {
    key <- tolower(trimws(val))
    if (!nzchar(key)) next
    if (!is.null(SALT_SYNONYM_LOOKUP[[key]])) {
      expanded <- c(expanded, SALT_SYNONYM_LOOKUP[[key]])
    }
  }
  unique_canonical(expanded)
}

normalize_lexeme_key_scalar <- function(value) {
  val <- collapse_ws(value)
  if (!length(val)) return(NA_character_)
  if (is.na(val) || !nzchar(val)) return(NA_character_)
  val <- tolower(val)
  val <- gsub("\\s+", " ", val, perl = TRUE)
  val <- trimws(val)
  val <- gsub("^[[:punct:]]+", "", val)
  val <- gsub("[[:punct:]]+$", "", val)
  val <- trimws(val)
  if (!nzchar(val)) return(NA_character_)
  val
}

normalize_lexeme_key <- function(value) {
  if (length(value) <= 1) {
    return(normalize_lexeme_key_scalar(value))
  }
  vapply(value, normalize_lexeme_key_scalar, character(1), USE.NAMES = FALSE)
}

write_csv_and_parquet <- function(dt, path) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("arrow package is required for Parquet export.")
  }
  arrow_table <- arrow::Table$create(dt)
  if ("write_csv_arrow" %chin% getNamespaceExports("arrow")) {
    arrow::write_csv_arrow(arrow_table, path)
  } else if ("write_delim_arrow" %chin% getNamespaceExports("arrow")) {
    arrow::write_delim_arrow(arrow_table, path, delim = ",")
  } else {
    data.table::fwrite(dt, path)
  }
  parquet_path <- sub("\\.csv$", ".parquet", path)
  arrow::write_parquet(arrow_table, parquet_path)
  parquet_path
}

find_generics_master_path <- function(script_dir, filename) {
  cwd <- normalizePath(getwd())
  base_dirs <- c(
    script_dir,
    cwd,
    dirname(script_dir),
    dirname(cwd),
    file.path(script_dir, ".."),
    file.path(script_dir, "..", ".."),
    file.path(cwd, ".."),
    file.path(cwd, "..", "..")
  )
  extra_dirs <- c()
  for (base in base_dirs) {
    if (is.na(base) || !nzchar(base)) next
    extra_dirs <- c(extra_dirs,
      file.path(base, "github_repos", "drugbank_generics"),
      file.path(base, "github_repos", "esoa"),
      file.path(base, "github_repos", "esoa", "dependencies", "drugbank_generics")
    )
  }
  base_dirs <- unique(vapply(c(base_dirs, extra_dirs), function(p) normalizePath(p, mustWork = FALSE), character(1)))
  candidate_paths <- character()
  for (base in base_dirs) {
    if (is.na(base) || !nzchar(base)) next
    candidate_paths <- c(candidate_paths,
      file.path(base, "output", filename),
      file.path(base, "dependencies", "drugbank_generics", "output", filename),
      file.path(base, "inputs", "drugs", filename)
    )
  }
  candidate_paths <- unique(vapply(candidate_paths, function(p) normalizePath(p, mustWork = FALSE), character(1)))
  for (path in candidate_paths) {
    if (!is.na(path) && file.exists(path)) {
      return(path)
    }
  }
  stop("Generics master not found in candidate paths:\n", paste(candidate_paths, collapse = "\n"))
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
mixtures_output_path <- file.path(output_dir, "drugbank_mixtures_master.csv")
generics_master_path <- find_generics_master_path(script_dir, "drugbank_generics_master.csv")

dataset <- drugbank

groups_dt <- as.data.table(dataset$drugs$groups)
groups_dt[, drugbank_id := as.character(drugbank_id)]
groups_dt[, group_clean := tolower(trimws(group))]
excluded_ids <- unique(groups_dt[group_clean %chin% c("vet"), drugbank_id])

filter_excluded <- function(dt, id_col = "drugbank_id") {
  dt[!(get(id_col) %chin% excluded_ids)]
}

groups_clean_dt <- groups_dt[!(drugbank_id %chin% excluded_ids) & !is.na(group_clean) & nzchar(group_clean)]
groups_clean_dt <- groups_clean_dt[, .(groups_list = list(unique_canonical(group_clean))), by = drugbank_id]
groups_lookup <- setNames(groups_clean_dt$groups_list, groups_clean_dt$drugbank_id)

salts_dt <- as.data.table(dataset$salts)[
  , .(drugbank_id = as.character(drugbank_id), salt_name = collapse_ws(name))
]
salts_dt <- filter_excluded(salts_dt)
salts_dt <- salts_dt[!is.na(salt_name) & nzchar(salt_name)]
salts_dt <- salts_dt[, .(salt_names_list = list(unique_canonical(salt_name))), by = drugbank_id]
salts_lookup <- setNames(salts_dt$salt_names_list, salts_dt$drugbank_id)

if (!file.exists(generics_master_path)) {
  stop("Generics master not found at ", generics_master_path)
}
generics_dt <- fread(generics_master_path, na.strings = c("", "NA"))
if (!all(c("drugbank_id", "lexeme", "generic_components_key") %chin% names(generics_dt))) {
  stop("Generics master missing required columns.")
}
generics_dt[, drugbank_id := as.character(drugbank_id)]
generics_dt[, lexeme_key := normalize_lexeme_key(lexeme)]
generics_lexeme_dt <- generics_dt[!is.na(lexeme_key) & nzchar(lexeme_key)]
lexeme_map <- split(generics_lexeme_dt$drugbank_id, generics_lexeme_dt$lexeme_key)
lexeme_map <- lapply(lexeme_map, unique)

generic_key_by_id <- generics_dt[
  !is.na(generic_components_key) & nzchar(generic_components_key),
  .(generic_components_key = generic_components_key[1L]),
  by = drugbank_id
]
generic_key_lookup <- setNames(generic_key_by_id$generic_components_key, generic_key_by_id$drugbank_id)

mixtures_dt <- as.data.table(dataset$drugs$mixtures)[
  , .(
    mixture_drugbank_id = as.character(drugbank_id),
    mixture_name = collapse_ws(name),
    ingredients_raw = collapse_ws(ingredients)
  )
]
mixtures_dt <- filter_excluded(mixtures_dt, "mixture_drugbank_id")
mixtures_dt <- mixtures_dt[!(is.na(mixture_name) & is.na(ingredients_raw))]
mixtures_dt[, mixture_name_key := normalize_lexeme_key(mixture_name)]
raw_components_list <- parallel_lapply(as.list(mixtures_dt$ingredients_raw), split_raw_components)
raw_components_char <- unlist(parallel_lapply(raw_components_list, function(vec) {
  if (!length(vec)) return(NA_character_)
  paste(vec, collapse = " ; ")
}), use.names = FALSE)
mixtures_dt[, component_raw_segments := raw_components_char]
ingredients_list <- parallel_lapply(as.list(mixtures_dt$ingredients_raw), split_ingredients)
mixtures_dt[, ingredient_components_vec := ingredients_list]
ingredient_components_char <- unlist(parallel_lapply(ingredients_list, function(vec) {
  if (!length(vec)) return(NA_character_)
  paste(vec, collapse = "; ")
}), use.names = FALSE)
mixtures_dt[, ingredient_components := ingredient_components_char]
ingredient_components_key_char <- unlist(parallel_lapply(ingredients_list, function(vec) {
  if (!length(vec)) return(NA_character_)
  keys <- unique(normalize_lexeme_key(vec))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (!length(keys)) return(NA_character_)
  paste(sort(keys), collapse = "||")
}), use.names = FALSE)
mixtures_dt[, ingredient_components_key := ingredient_components_key_char]

resolve_component_ids <- function(vec) {
  if (is.null(vec) || !length(vec)) return(character())
  ids <- unique(unlist(lapply(vec, function(comp) {
    key <- normalize_lexeme_key(comp)
    if (is.null(key) || is.na(key) || !nzchar(key)) return(character())
    vals <- lexeme_map[[key]]
    if (is.null(vals)) character() else vals
  }), use.names = FALSE))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  unique(ids)
}

component_ids_list <- parallel_lapply(mixtures_dt$ingredient_components_vec, resolve_component_ids)
mixtures_dt[, component_ids_list := component_ids_list]
component_generic_keys_list <- parallel_lapply(component_ids_list, function(ids) {
  if (!length(ids)) return(character())
  vals <- unique(generic_key_lookup[ids])
  vals <- vals[!is.na(vals) & nzchar(vals)]
  unique(vals)
})
mixtures_dt[, component_generic_keys_list := component_generic_keys_list]
component_lexemes_char <- unlist(parallel_lapply(mixtures_dt$ingredient_components_vec, collapse_pipe), use.names = FALSE)
mixtures_dt[, component_lexemes := component_lexemes_char]
component_drugbank_ids_char <- unlist(parallel_lapply(component_ids_list, collapse_pipe), use.names = FALSE)
mixtures_dt[, component_drugbank_ids := component_drugbank_ids_char]
component_generic_keys_char <- unlist(parallel_lapply(component_generic_keys_list, collapse_pipe), use.names = FALSE)
mixtures_dt[, component_generic_keys := component_generic_keys_char]

groups_char <- unlist(parallel_lapply(as.list(mixtures_dt$mixture_drugbank_id), function(id) {
  vals <- groups_lookup[[id]]
  if (is.null(vals)) return(NA_character_)
  collapse_pipe(vals)
}), use.names = FALSE)
mixtures_dt[, groups := groups_char]

salt_names_char <- unlist(parallel_lapply(as.list(mixtures_dt$mixture_drugbank_id), function(id) {
  vals <- salts_lookup[[id]]
  if (is.null(vals)) return(NA_character_)
  collapse_pipe(expand_salt_set(vals))
}), use.names = FALSE)
mixtures_dt[, salt_names := salt_names_char]
mixtures_dt[, c("ingredient_components_vec", "component_ids_list", "component_generic_keys_list") := NULL]

mixtures_dt[, mixture_id := .I]
setcolorder(mixtures_dt, c(
  "mixture_id",
  "mixture_drugbank_id",
  "mixture_name",
  "mixture_name_key",
  "ingredients_raw",
  "component_raw_segments",
  "ingredient_components",
  "ingredient_components_key",
  "component_lexemes",
  "component_drugbank_ids",
  "component_generic_keys",
  "groups",
  "salt_names"
))

setorder(mixtures_dt, mixture_name_key, mixture_drugbank_id, mixture_id)

write_csv_and_parquet(mixtures_dt, mixtures_output_path)

cat(sprintf("Wrote %d rows to %s\n", nrow(mixtures_dt), mixtures_output_path))
if (!quiet_mode) {
  cat("Sample rows:\n")
  print(head(mixtures_dt, 5))
}

if (!is.null(plan_reset)) {
  try(future::plan(future::sequential), silent = TRUE)
}
