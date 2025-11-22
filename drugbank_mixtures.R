#!/usr/bin/env Rscript
# drugbank_mixtures.R — build a mixtures-focused DrugBank master dataset

suppressWarnings({
  suppressPackageStartupMessages({
    ensure_installed <- function(pkg, repos = "https://cloud.r-project.org") {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = repos)
      }
    }
    install_tidypolars <- function() {
      if (!requireNamespace("tidypolars", quietly = TRUE)) {
        Sys.setenv(NOT_CRAN = "true")
        install.packages(
          "tidypolars",
          repos = c("https://community.r-multiverse.org", "https://cloud.r-project.org")
        )
      }
    }
    install_tidypolars()
    ensure_installed("dbdataset")
  })
})

suppressPackageStartupMessages({
  library(tidypolars)
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr, warn.conflicts = FALSE)
  library(dbdataset)
})

pl <- polars::pl

argv <- commandArgs(trailingOnly = TRUE)
quiet_mode <- identical(tolower(Sys.getenv("ESOA_DRUGBANK_QUIET", "0")), "1")

worker_env <- suppressWarnings(as.integer(Sys.getenv("ESOA_DRUGBANK_WORKERS", "")))
if (!is.na(worker_env) && worker_env > 0) {
  Sys.setenv(POLARS_MAX_THREADS = worker_env)
}

collapse_ws <- function(x) {
  ifelse(is.na(x), NA_character_, trimws(gsub("\\s+", " ", as.character(x))))
}

unique_canonical <- function(values) {
  vals <- collapse_ws(values)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (!length(vals)) return(character())
  vals <- vals[order(tolower(vals), vals)]
  keys <- tolower(vals)
  vals[!duplicated(keys)]
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

write_csv_and_parquet <- function(df, csv_path) {
  parquet_path <- sub("\\.csv$", ".parquet", csv_path)
  df$write_parquet(parquet_path)
  df$write_csv(csv_path)
  parquet_path
}

find_generics_master_path <- function(script_dir, stem) {
  cwd <- normalizePath(getwd())
  base_dirs <- unique(vapply(
    c(
      script_dir,
      cwd,
      dirname(script_dir),
      dirname(cwd),
      file.path(script_dir, ".."),
      file.path(script_dir, "..", ".."),
      file.path(cwd, ".."),
      file.path(cwd, "..", "..")
    ),
    function(p) normalizePath(p, mustWork = FALSE),
    character(1)
  ))
  candidate_paths <- character()
  for (base in base_dirs) {
    if (is.na(base) || !nzchar(base)) next
    candidate_paths <- c(candidate_paths,
      file.path(base, "output", paste0(stem, ".parquet")),
      file.path(base, "dependencies", "drugbank_generics", "output", paste0(stem, ".parquet")),
      file.path(base, "inputs", "drugs", paste0(stem, ".parquet"))
    )
  }
  candidate_paths <- unique(vapply(candidate_paths, function(p) normalizePath(p, mustWork = FALSE), character(1)))
  for (path in candidate_paths) {
    if (!is.na(path) && file.exists(path)) {
      return(path)
    }
  }
  stop("Generics master (Parquet) not found in candidate paths:\n", paste(candidate_paths, collapse = "\n"))
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

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
mixtures_output_csv <- file.path(output_dir, "drugbank_mixtures_master.csv")
generics_master_path <- find_generics_master_path(script_dir, "drugbank_generics_master")

dataset <- drugbank

groups_df <- as_polars_df(dataset$drugs$groups) |>
  mutate(
    drugbank_id = pl$col("drugbank_id")$cast(pl$String),
    group_clean = pl$col("group")$map_elements(function(x) tolower(trimws(x)), return_dtype = pl$String)
  )
excluded_ids <- groups_df |>
  filter(pl$col("group_clean") == "vet") |>
  pull("drugbank_id") |>
  unique() |>
  to_r()

filter_excluded <- function(df, id_col = "drugbank_id") {
  if (!length(excluded_ids)) return(df)
  dplyr::filter(df, !pl$col(id_col)$is_in(excluded_ids))
}

groups_clean <- groups_df |>
  filter(
    !pl$col("drugbank_id")$is_in(excluded_ids) &
      pl$col("group_clean")$is_not_null() &
      pl$col("group_clean") != ""
  ) |>
  group_by(drugbank_id) |>
  summarise(groups_list = pl$col("group_clean")) |>
  mutate(
    groups_list = pl$col("groups_list")$map_elements(unique_canonical, return_dtype = pl$List(pl$String))
  )
groups_lookup <- groups_clean$to_r()
groups_lookup <- setNames(groups_lookup$groups_list, groups_lookup$drugbank_id)

salts_df <- as_polars_df(dataset$salts) |>
  mutate(
    drugbank_id = pl$col("drugbank_id")$cast(pl$String),
    salt_name = pl$col("name")$map_elements(collapse_ws, return_dtype = pl$String)
  ) |>
  select(drugbank_id, salt_name) |>
  filter_excluded() |>
  filter(pl$col("salt_name")$is_not_null() & pl$col("salt_name") != "") |>
  group_by(drugbank_id) |>
  summarise(salt_names_list = pl$col("salt_name")) |>
  mutate(
    salt_names_list = pl$col("salt_names_list")$map_elements(unique_canonical, return_dtype = pl$List(pl$String))
  )
salts_lookup <- salts_df$to_r()
salts_lookup <- setNames(salts_lookup$salt_names_list, salts_lookup$drugbank_id)

generics_df <- read_parquet_polars(generics_master_path)
required_cols <- c("drugbank_id", "lexeme", "generic_components_key")
missing_cols <- setdiff(required_cols, generics_df$names())
if (length(missing_cols)) {
  stop("Generics master missing required columns: ", paste(missing_cols, collapse = ", "))
}
generics_df <- generics_df |>
  mutate(
    lexeme_key = pl$col("lexeme")$map_elements(normalize_lexeme_key_scalar, return_dtype = pl$String)
  )
generics_lexeme_df <- generics_df |>
  select("drugbank_id", "lexeme_key") |>
  to_r()
generics_lexeme_df <- generics_lexeme_df[!is.na(generics_lexeme_df$lexeme_key) & nzchar(generics_lexeme_df$lexeme_key), ]
lexeme_map <- split(generics_lexeme_df$drugbank_id, generics_lexeme_df$lexeme_key)
lexeme_map <- lapply(lexeme_map, unique)

generic_key_lookup <- generics_df |>
  select("drugbank_id", "generic_components_key") |>
  to_r()
generic_key_lookup <- generic_key_lookup[!is.na(generic_key_lookup$generic_components_key) & nzchar(generic_key_lookup$generic_components_key), ]
generic_key_lookup <- generic_key_lookup[!duplicated(generic_key_lookup$drugbank_id), ]
generic_key_lookup <- setNames(generic_key_lookup$generic_components_key, generic_key_lookup$drugbank_id)

mixtures_df <- as_polars_df(dataset$drugs$mixtures) |>
  mutate(
    mixture_drugbank_id = pl$col("drugbank_id")$cast(pl$String),
    mixture_name = pl$col("name")$map_elements(collapse_ws, return_dtype = pl$String),
    ingredients_raw = pl$col("ingredients")$map_elements(collapse_ws, return_dtype = pl$String)
  ) |>
  select(mixture_drugbank_id, mixture_name, ingredients_raw) |>
  filter_excluded(id_col = "mixture_drugbank_id") |>
  filter(!(pl$col("mixture_name")$is_null() & pl$col("ingredients_raw")$is_null())) |>
  mutate(
    mixture_name_key = pl$col("mixture_name")$map_elements(normalize_lexeme_key_scalar, return_dtype = pl$String),
    raw_components_list = pl$col("ingredients_raw")$map_elements(split_raw_components, return_dtype = pl$List(pl$String)),
    ingredient_components_vec = pl$col("ingredients_raw")$map_elements(split_ingredients, return_dtype = pl$List(pl$String))
  ) |>
  mutate(
    component_raw_segments = pl$col("raw_components_list")$map_elements(
      function(vec) {
        if (!length(vec)) return(NA_character_)
        paste(vec, collapse = " ; ")
      },
      return_dtype = pl$String
    ),
    ingredient_components = pl$col("ingredient_components_vec")$map_elements(
      function(vec) {
        if (!length(vec)) return(NA_character_)
        paste(vec, collapse = "; ")
      },
      return_dtype = pl$String
    ),
    ingredient_components_key = pl$col("ingredient_components_vec")$map_elements(
      function(vec) {
        if (!length(vec)) return(NA_character_)
        keys <- unique(normalize_lexeme_key(vec))
        keys <- keys[!is.na(keys) & nzchar(keys)]
        if (!length(keys)) return(NA_character_)
        paste(sort(keys), collapse = "||")
      },
      return_dtype = pl$String
    )
  )

resolve_component_ids <- function(vec, lex_map) {
  if (is.null(vec) || !length(vec)) return(character())
  ids <- unique(unlist(lapply(vec, function(comp) {
    key <- normalize_lexeme_key_scalar(comp)
    if (is.null(key) || is.na(key) || !nzchar(key)) return(character())
    vals <- lex_map[[key]]
    if (is.null(vals)) character() else vals
  }), use.names = FALSE))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  unique(ids)
}

mixtures_r <- mixtures_df$to_r()
mixtures_r$component_ids_list <- lapply(mixtures_r$ingredient_components_vec, resolve_component_ids, lex_map = lexeme_map)
mixtures_r$component_generic_keys_list <- lapply(mixtures_r$component_ids_list, function(ids) {
  if (!length(ids)) return(character())
  vals <- unique(generic_key_lookup[ids])
  vals <- vals[!is.na(vals) & nzchar(vals)]
  unique(vals)
})
mixtures_r$component_lexemes <- vapply(mixtures_r$ingredient_components_vec, collapse_pipe, character(1), USE.NAMES = FALSE)
mixtures_r$component_drugbank_ids <- vapply(mixtures_r$component_ids_list, collapse_pipe, character(1), USE.NAMES = FALSE)
mixtures_r$component_generic_keys <- vapply(mixtures_r$component_generic_keys_list, collapse_pipe, character(1), USE.NAMES = FALSE)
mixtures_r$groups <- vapply(mixtures_r$mixture_drugbank_id, function(id) {
  vals <- groups_lookup[[id]]
  if (is.null(vals)) return(NA_character_)
  collapse_pipe(vals)
}, character(1), USE.NAMES = FALSE)
mixtures_r$salt_names <- vapply(mixtures_r$mixture_drugbank_id, function(id) {
  vals <- salts_lookup[[id]]
  if (is.null(vals)) return(NA_character_)
  collapse_pipe(expand_salt_set(vals))
}, character(1), USE.NAMES = FALSE)
mixtures_r$mixture_id <- seq_len(nrow(mixtures_r))

mixtures_df <- as_polars_df(mixtures_r) |>
  mutate(
    component_ids_list = pl$col("component_ids_list")$cast(pl$List(pl$String)),
    component_generic_keys_list = pl$col("component_generic_keys_list")$cast(pl$List(pl$String))
  ) |>
  select(
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
  ) |>
  arrange(mixture_name_key, mixture_drugbank_id, mixture_id)

write_csv_and_parquet(mixtures_df, mixtures_output_csv)

cat(sprintf("Wrote %d rows to %s (Parquet + CSV)\n", mixtures_df$height(), mixtures_output_csv))
if (!quiet_mode) {
  cat("Sample rows:\n")
  print(mixtures_df$head(5)$to_r())
}
