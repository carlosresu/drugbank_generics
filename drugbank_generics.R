#!/usr/bin/env Rscript
# drugbank_generics.R — build a generics-focused DrugBank master dataset
# Columns: drugbank_id, generic_name, generic_components, generic_components_key,
#          synonyms, normalized_doses, raw_doses, normalized_forms, raw_forms,
#          normalized_routes, raw_routes, salt_names, atc_codes, groups
# Only veterinary-only drugs are excluded; all other groups remain.

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

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmd_args)
  if (length(match) > 0) {
    return(normalizePath(dirname(sub(needle, "", cmd_args[match[1]]))))
  }
  return(normalizePath(getwd()))
}

paths_equal <- function(a, b) {
  normalizePath(a, mustWork = FALSE) == normalizePath(b, mustWork = FALSE)
}

safe_copy <- function(src, dest) {
  if (is.null(src) || is.null(dest) || is.na(src) || is.na(dest)) {
    stop("safe_copy: invalid path")
  }
  if (paths_equal(src, dest)) {
    return(invisible(dest))
  }
  dest_dir <- dirname(dest)
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(src, dest, overwrite = TRUE, copy.mode = TRUE)
  invisible(dest)
}

copy_outputs_to_superproject <- function(output_path) {
  script_dir <- get_script_dir()
  super_root <- normalizePath(file.path(script_dir, "..", ".."))
  dest1 <- file.path(super_root, "dependencies", "drugbank_generics", "output", basename(output_path))
  dest2 <- file.path(super_root, "inputs", "drugs", "drugbank_generics_master.csv")
  safe_copy(output_path, dest1)
  safe_copy(output_path, dest2)
}

collapse_ws <- function(x) {
  ifelse(is.na(x), NA_character_, trimws(gsub("\\s+", " ", x)))
}

empty_to_na <- function(x) {
  val <- collapse_ws(x)
  ifelse(val == "" | is.na(val), NA_character_, val)
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
  parts <- unlist(strsplit(val, "(?i)\\sand\\s|\\swith\\s|\\splus\\s|\"?\\+\"?|/|,|;", perl = TRUE))
  parts <- collapse_ws(parts)
  parts <- parts[nzchar(parts)]
  unique_canonical(parts)
}

split_semicolon <- function(value) {
  val <- collapse_ws(value)
  if (is.na(val) || !nzchar(val)) return(character())
  collapse_ws(unlist(strsplit(val, ";")))
}

collapse_pipe <- function(values) {
  vals <- unique_canonical(values)
  if (!length(vals)) return(NA_character_)
  paste(vals, collapse = "|")
}

clean_numeric_string <- function(x) {
  tolower(gsub(",", "", collapse_ws(x)))
}

extract_numeric <- function(x) {
  num <- suppressWarnings(as.numeric(x))
  ifelse(is.na(num), NA_real_, num)
}

mass_to_mg <- function(value, unit) {
  if (is.na(value) || is.na(unit)) return(NA_real_)
  switch(unit,
         "mg" = value,
         "g" = value * 1000,
         "mcg" = value / 1000,
         "µg" = value / 1000,
         "ug" = value / 1000,
         "kg" = value * 1e6,
         value)
}

format_numeric <- function(x) {
  if (is.na(x)) return(NA_character_)
  if (abs(x - round(x)) < 1e-9) {
    sprintf("%d", as.integer(round(x)))
  } else {
    format(round(x, 3), trim = TRUE, scientific = FALSE)
  }
}

normalize_dose_value <- function(value) {
  val <- clean_numeric_string(value)
  if (is.na(val) || !nzchar(val)) return(NA_character_)
  pattern <- "^([0-9]+(?:\\.[0-9]+)?)\\s*(mg|g|mcg|ug|µg|kg)(?:\\s*/\\s*([0-9]+(?:\\.[0-9]+)?)\\s*(ml|l))?$"
  if (!grepl(pattern, val, perl = TRUE)) return(val)
  match <- regexec(pattern, val, perl = TRUE)
  pieces <- regmatches(val, match)[[1]]
  dose_val <- mass_to_mg(extract_numeric(pieces[2]), pieces[3])
  dose_fmt <- format_numeric(dose_val)
  if (length(pieces) >= 5 && nzchar(pieces[4])) {
    per_val <- extract_numeric(pieces[4])
    per_unit <- pieces[5]
    per_norm <- ifelse(per_unit == "l", per_val * 1000, per_val)
    per_fmt <- format_numeric(per_norm)
    paste0(dose_fmt, " mg/", per_fmt, " mL")
  } else {
    paste0(dose_fmt, " mg")
  }
}

normalize_dose_vector <- function(values) {
  unique_canonical(na.omit(sapply(values, normalize_dose_value, USE.NAMES = FALSE)))
}

FORM_CANONICAL_MAP <- c(
  "tabs" = "tablet",
  "tab" = "tablet",
  "tablets" = "tablet",
  "caps" = "capsule",
  "cap" = "capsule",
  "capsules" = "capsule",
  "sol" = "solution",
  "susp" = "suspension",
  "inj" = "injection",
  "liq" = "solution",
  "liquid" = "solution"
)

normalize_form_value <- function(value) {
  val <- tolower(collapse_ws(value))
  if (is.na(val) || !nzchar(val)) return(NA_character_)
  FORM_CANONICAL_MAP[[val]] %||% val
}

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b

normalize_form_vector <- function(values) {
  unique_canonical(na.omit(sapply(values, normalize_form_value, USE.NAMES = FALSE)))
}

ROUTE_ALIAS_MAP <- c(
  "po" = "oral",
  "per os" = "oral",
  "by mouth" = "oral",
  "iv" = "intravenous",
  "im" = "intramuscular",
  "sc" = "subcutaneous",
  "subcut" = "subcutaneous",
  "subling" = "sublingual",
  "sl" = "sublingual"
)

normalize_route_entry <- function(value) {
  val <- tolower(collapse_ws(value))
  if (is.na(val) || !nzchar(val)) return(NA_character_)
  ROUTE_ALIAS_MAP[[val]] %||% val
}

normalize_route_vector <- function(values) {
  unique_canonical(na.omit(sapply(values, normalize_route_entry, USE.NAMES = FALSE)))
}

normalize_list_column <- function(values) {
  if (is.null(values) || !length(values)) return(list(character()))
  list(unique_canonical(values))
}

combine_list_column <- function(...) {
  values <- combine_values(...)
  list(values)
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_master_path <- file.path(output_dir, "drugbank_generics_master.csv")

dataset <- drugbank

groups_dt <- as.data.table(dataset$drugs$groups)
groups_dt[, drugbank_id := as.character(drugbank_id)]
groups_dt[, group_clean := tolower(trimws(group))]
excluded_groups <- c("vet")
excluded_ids <- unique(groups_dt[group_clean %chin% excluded_groups, as.character(drugbank_id)])

filter_excluded <- function(dt, id_col = "drugbank_id") {
  dt <- dt[!(get(id_col) %chin% excluded_ids)]
  dt
}

general_dt <- as.data.table(dataset$drugs$general_information)[
  , .(drugbank_id = as.character(drugbank_id), generic_name = collapse_ws(name))
]
general_dt <- filter_excluded(general_dt)
general_dt <- general_dt[!is.na(generic_name) & nzchar(generic_name)]
general_dt[, generic_components := lapply(generic_name, split_ingredients)]

syn_dt <- as.data.table(dataset$drugs$synonyms)[
  , .(drugbank_id = as.character(drugbank_id), synonym = collapse_ws(synonym))
]
syn_dt <- filter_excluded(syn_dt)
syn_dt <- syn_dt[!is.na(synonym) & nzchar(synonym)]
syn_dt <- syn_dt[, .(synonyms_list = list(unique_canonical(synonym))), by = drugbank_id]

atc_dt <- as.data.table(dataset$drugs$atc_codes)[
  , .(drugbank_id = as.character(drugbank_id), atc_code = trimws(atc_code))
]
atc_dt <- filter_excluded(atc_dt)
atc_dt <- atc_dt[!is.na(atc_code) & nzchar(atc_code)]
atc_dt <- atc_dt[, .(atc_codes_list = list(unique_canonical(atc_code))), by = drugbank_id]

groups_clean_dt <- groups_dt[!(drugbank_id %chin% excluded_ids)]
groups_clean_dt <- groups_clean_dt[, .(groups_list = list(unique_canonical(group_clean))), by = drugbank_id]

dosages_dt <- as.data.table(dataset$drugs$dosages)[
  , .(drugbank_id = as.character(drugbank_id), route, form, strength)
]
dosages_dt <- filter_excluded(dosages_dt)
dosages_dt[, route := empty_to_na(route)]
dosages_dt[, form := empty_to_na(form)]
dosages_dt[, strength := empty_to_na(strength)]
dosages_dt <- dosages_dt[, .(
  dosage_routes_list = list(unique_canonical(route)),
  dosage_forms_list = list(unique_canonical(form)),
  dosage_doses_list = list(unique_canonical(strength))
), by = drugbank_id]

products_dt <- as.data.table(dataset$products)[
  , .(drugbank_id = as.character(drugbank_id), route, dosage_form, strength)
]
products_dt <- filter_excluded(products_dt)
products_dt[, route := empty_to_na(route)]
products_dt[, dosage_form := empty_to_na(dosage_form)]
products_dt[, strength := empty_to_na(strength)]
products_dt <- products_dt[, .(
  product_routes_list = list(unique_canonical(route)),
  product_forms_list = list(unique_canonical(dosage_form)),
  product_doses_list = list(unique_canonical(strength))
), by = drugbank_id]

dose_form_route_dt <- merge(products_dt, dosages_dt, by = "drugbank_id", all = TRUE)
dose_form_route_dt[is.na(product_routes_list), product_routes_list := list(character())]
dose_form_route_dt[is.na(dosage_routes_list), dosage_routes_list := list(character())]
dose_form_route_dt[is.na(product_forms_list), product_forms_list := list(character())]
dose_form_route_dt[is.na(dosage_forms_list), dosage_forms_list := list(character())]
dose_form_route_dt[is.na(product_doses_list), product_doses_list := list(character())]
dose_form_route_dt[is.na(dosage_doses_list), dosage_doses_list := list(character())]
dose_form_route_dt[, raw_routes_list := lapply(seq_len(.N), function(i) combine_values(product_routes_list[[i]], dosage_routes_list[[i]]))]
dose_form_route_dt[, raw_forms_list := lapply(seq_len(.N), function(i) combine_values(product_forms_list[[i]], dosage_forms_list[[i]]))]
dose_form_route_dt[, raw_doses_list := lapply(seq_len(.N), function(i) combine_values(product_doses_list[[i]], dosage_doses_list[[i]]))]
dose_form_route_dt[, normalized_routes_list := lapply(raw_routes_list, normalize_route_vector)]
dose_form_route_dt[, normalized_forms_list := lapply(raw_forms_list, normalize_form_vector)]
dose_form_route_dt[, normalized_doses_list := lapply(raw_doses_list, normalize_dose_vector)]
dose_form_route_dt <- dose_form_route_dt[, .(
  drugbank_id,
  raw_routes_list,
  normalized_routes_list,
  raw_forms_list,
  normalized_forms_list,
  raw_doses_list,
  normalized_doses_list
)]

salts_dt <- as.data.table(dataset$salts)[
  , .(drugbank_id = as.character(drugbank_id), salt_name = collapse_ws(name))
]
salts_dt <- filter_excluded(salts_dt)
salts_dt <- salts_dt[!is.na(salt_name) & nzchar(salt_name)]
salts_dt <- salts_dt[, .(salt_names_list = list(unique_canonical(salt_name))), by = drugbank_id]

general_dt[, generic_components_vec := lapply(generic_components, unique_canonical)]
general_dt[, generic_components_key := sapply(generic_components_vec, function(vec) {
  if (!length(vec)) return(NA_character_)
  canon <- tolower(trimws(gsub("\\s+", " ", vec)))
  canon <- canon[nzchar(canon)]
  if (!length(canon)) return(NA_character_)
  paste(sort(canon), collapse = "||")
})]

final_dt <- copy(general_dt[, .(drugbank_id, generic_name, generic_components_vec, generic_components_key)])
final_dt <- merge(final_dt, syn_dt, by = "drugbank_id", all.x = TRUE)
final_dt <- merge(final_dt, atc_dt, by = "drugbank_id", all.x = TRUE)
final_dt <- merge(final_dt, groups_clean_dt, by = "drugbank_id", all.x = TRUE)
final_dt <- merge(final_dt, dose_form_route_dt, by = "drugbank_id", all.x = TRUE)
final_dt <- merge(final_dt, salts_dt, by = "drugbank_id", all.x = TRUE)

final_dt[is.na(synonyms_list), synonyms_list := list(character())]
final_dt[is.na(atc_codes_list), atc_codes_list := list(character())]
final_dt[is.na(groups_list), groups_list := list(character())]
final_dt[is.na(raw_routes_list), raw_routes_list := list(character())]
final_dt[is.na(normalized_routes_list), normalized_routes_list := list(character())]
final_dt[is.na(raw_forms_list), raw_forms_list := list(character())]
final_dt[is.na(normalized_forms_list), normalized_forms_list := list(character())]
final_dt[is.na(raw_doses_list), raw_doses_list := list(character())]
final_dt[is.na(normalized_doses_list), normalized_doses_list := list(character())]
final_dt[is.na(salt_names_list), salt_names_list := list(character())]

final_dt[, generic_components := sapply(generic_components_vec, function(vec) {
  if (!length(vec)) return(NA_character_)
  paste(vec, collapse = "; ")
})]

final_dt[, synonyms := sapply(synonyms_list, collapse_pipe)]
final_dt[, normalized_doses := sapply(normalized_doses_list, collapse_pipe)]
final_dt[, raw_doses := sapply(raw_doses_list, collapse_pipe)]
final_dt[, normalized_forms := sapply(normalized_forms_list, collapse_pipe)]
final_dt[, raw_forms := sapply(raw_forms_list, collapse_pipe)]
final_dt[, normalized_routes := sapply(normalized_routes_list, collapse_pipe)]
final_dt[, raw_routes := sapply(raw_routes_list, collapse_pipe)]
final_dt[, salt_names := sapply(salt_names_list, collapse_pipe)]
final_dt[, atc_codes := sapply(atc_codes_list, collapse_pipe)]
final_dt[, groups := sapply(groups_list, collapse_pipe)]

final_dt <- final_dt[, .(
  drugbank_id,
  generic_name,
  generic_components,
  generic_components_key,
  synonyms,
  normalized_doses,
  raw_doses,
  normalized_forms,
  raw_forms,
  normalized_routes,
  raw_routes,
  salt_names,
  atc_codes,
  groups
)]

setorder(final_dt, generic_name, drugbank_id)

write_arrow_csv <- function(dt, path) {
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_csv_arrow(dt, path)
  } else {
    fwrite(dt, path)
  }
}

write_arrow_csv(final_dt, output_master_path)
copy_outputs_to_superproject(output_master_path)

cat(sprintf("Wrote %d rows to %s\n", nrow(final_dt), output_master_path))
cat("Sample rows:\n")
print(head(final_dt, 5))
