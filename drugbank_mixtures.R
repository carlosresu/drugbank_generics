#!/usr/bin/env Rscript
# drugbank_mixtures.R — build a mixtures-focused DrugBank master dataset

suppressWarnings({
  suppressPackageStartupMessages({
    ensure_installed <- function(pkg) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }
    ensure_installed("dbdataset")
    ensure_installed("polars")
  })
})

tryCatch({
  ensure_installed("arrow")
}, error = function(e) {})

library(dbdataset)
library(polars)

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
cat(sprintf("[drugbank_mixtures] parallel workers: %s\n", worker_count))

choose_backend <- function() {
  backend <- tolower(Sys.getenv("ESOA_DRUGBANK_BACKEND", "future"))
  if (backend == "auto") {
    if (.Platform$OS.type == "unix") return("mclapply")
    return("future")
  }
  if (backend %in% c("mclapply", "future")) return(backend)
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

write_csv_and_parquet <- function(df, path) {
  parquet_path <- sub("\\.csv$", ".parquet", path)
  polars_df <- pl$DataFrame(as.data.frame(df))
  polars_df$write_parquet(parquet_path)
  polars_df$write_csv(path)
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

groups_df <- as.data.frame(dataset$drugs$groups, stringsAsFactors = FALSE)
groups_df$drugbank_id <- as.character(groups_df$drugbank_id)
groups_df$group_clean <- tolower(trimws(groups_df$group))
excluded_ids <- unique(groups_df$drugbank_id[groups_df$group_clean %in% c("vet")])

filter_excluded <- function(df, id_col = "drugbank_id") {
  df[!(df[[id_col]] %in% excluded_ids), , drop = FALSE]
}

groups_clean_df <- groups_df[!(groups_df$drugbank_id %in% excluded_ids) & !is.na(groups_df$group_clean) & nzchar(groups_df$group_clean), , drop = FALSE]
groups_split <- split(groups_clean_df$group_clean, groups_clean_df$drugbank_id)
groups_lookup <- lapply(groups_split, unique_canonical)

salts_df <- as.data.frame(dataset$salts, stringsAsFactors = FALSE)
salts_df <- salts_df[, c("drugbank_id", "name"), drop = FALSE]
names(salts_df) <- c("drugbank_id", "salt_name")
salts_df$drugbank_id <- as.character(salts_df$drugbank_id)
salts_df$salt_name <- collapse_ws(salts_df$salt_name)
salts_df <- filter_excluded(salts_df)
salts_df <- salts_df[!is.na(salts_df$salt_name) & nzchar(salts_df$salt_name), , drop = FALSE]
salts_split <- split(salts_df$salt_name, salts_df$drugbank_id)
salts_lookup <- lapply(salts_split, function(vals) unique_canonical(vals))

if (!file.exists(generics_master_path)) {
  stop("Generics master not found at ", generics_master_path)
}
generics_df <- pl$read_parquet(generics_master_path)$to_r()
required_generics_cols <- c("drugbank_id", "lexeme", "generic_components_key")
if (!all(required_generics_cols %in% names(generics_df))) {
  stop("Generics master missing required columns.")
}
generics_df$drugbank_id <- as.character(generics_df$drugbank_id)
generics_df$lexeme_key <- normalize_lexeme_key(generics_df$lexeme)
generics_lexeme_df <- generics_df[!is.na(generics_df$lexeme_key) & nzchar(generics_df$lexeme_key), , drop = FALSE]
lexeme_map <- split(generics_lexeme_df$drugbank_id, generics_lexeme_df$lexeme_key)
lexeme_map <- lapply(lexeme_map, unique)

valid_generic_keys <- !is.na(generics_df$generic_components_key) & nzchar(generics_df$generic_components_key)
generic_key_by_id <- generics_df[valid_generic_keys, c("drugbank_id", "generic_components_key"), drop = FALSE]
generic_key_lookup <- setNames(generic_key_by_id$generic_components_key, generic_key_by_id$drugbank_id)

mixtures_df <- data.frame(
  mixture_drugbank_id = as.character(dataset$drugs$mixtures$drugbank_id),
  mixture_name = collapse_ws(dataset$drugs$mixtures$name),
  ingredients_raw = collapse_ws(dataset$drugs$mixtures$ingredients),
  stringsAsFactors = FALSE
)
mixtures_df <- filter_excluded(mixtures_df, "mixture_drugbank_id")
mixtures_df <- mixtures_df[!(is.na(mixtures_df$mixture_name) & is.na(mixtures_df$ingredients_raw)), , drop = FALSE]
mixtures_df$mixture_name_key <- normalize_lexeme_key(mixtures_df$mixture_name)
raw_components_list <- parallel_lapply(as.list(mixtures_df$ingredients_raw), split_raw_components)
raw_components_char <- unlist(parallel_lapply(raw_components_list, function(vec) {
  if (!length(vec)) return(NA_character_)
  paste(vec, collapse = " ; ")
}), use.names = FALSE)
mixtures_df$component_raw_segments <- raw_components_char
ingredients_list <- parallel_lapply(as.list(mixtures_df$ingredients_raw), split_ingredients)
mixtures_df$ingredient_components_vec <- ingredients_list
ingredient_components_char <- unlist(parallel_lapply(ingredients_list, function(vec) {
  if (!length(vec)) return(NA_character_)
  paste(vec, collapse = "; ")
}), use.names = FALSE)
mixtures_df$ingredient_components <- ingredient_components_char
ingredient_components_key_char <- unlist(parallel_lapply(ingredients_list, function(vec) {
  if (!length(vec)) return(NA_character_)
  keys <- unique(normalize_lexeme_key(vec))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (!length(keys)) return(NA_character_)
  paste(sort(keys), collapse = "||")
}), use.names = FALSE)
mixtures_df$ingredient_components_key <- ingredient_components_key_char

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

component_ids_list <- parallel_lapply(mixtures_df$ingredient_components_vec, resolve_component_ids)
mixtures_df$component_ids_list <- component_ids_list
component_generic_keys_list <- parallel_lapply(component_ids_list, function(ids) {
  if (!length(ids)) return(character())
  vals <- unique(generic_key_lookup[ids])
  vals <- vals[!is.na(vals) & nzchar(vals)]
  unique(vals)
})
mixtures_df$component_generic_keys_list <- component_generic_keys_list
component_lexemes_char <- unlist(parallel_lapply(mixtures_df$ingredient_components_vec, collapse_pipe), use.names = FALSE)
mixtures_df$component_lexemes <- component_lexemes_char
component_drugbank_ids_char <- unlist(parallel_lapply(component_ids_list, collapse_pipe), use.names = FALSE)
mixtures_df$component_drugbank_ids <- component_drugbank_ids_char
component_generic_keys_char <- unlist(parallel_lapply(component_generic_keys_list, collapse_pipe), use.names = FALSE)
mixtures_df$component_generic_keys <- component_generic_keys_char

groups_char <- unlist(parallel_lapply(as.list(mixtures_df$mixture_drugbank_id), function(id) {
  vals <- groups_lookup[[id]]
  if (is.null(vals)) return(NA_character_)
  collapse_pipe(vals)
}), use.names = FALSE)
mixtures_df$groups <- groups_char

salt_names_char <- unlist(parallel_lapply(as.list(mixtures_df$mixture_drugbank_id), function(id) {
  vals <- salts_lookup[[id]]
  if (is.null(vals)) return(NA_character_)
  collapse_pipe(expand_salt_set(vals))
}), use.names = FALSE)
mixtures_df$salt_names <- salt_names_char
mixtures_df$mixture_id <- seq_len(nrow(mixtures_df))
mixtures_df <- mixtures_df[, !(names(mixtures_df) %in% c("ingredient_components_vec", "component_ids_list", "component_generic_keys_list")), drop = FALSE]

mixtures_df <- mixtures_df[, c(
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
)]

mixtures_df <- mixtures_df[order(mixtures_df$mixture_name_key, mixtures_df$mixture_drugbank_id, mixtures_df$mixture_id), , drop = FALSE]

write_csv_and_parquet(mixtures_df, mixtures_output_path)

  cat(sprintf("Wrote %d rows to %s\n", nrow(mixtures_df), mixtures_output_path))
  if (!quiet_mode) {
    cat("Sample rows:\n")
    print(head(mixtures_df, 5))
  }

if (!is.null(plan_reset)) {
  try(future::plan(future::sequential), silent = TRUE)
}
