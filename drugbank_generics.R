#!/usr/bin/env Rscript
# drugbank_generics.R — build a generics-focused DrugBank master dataset
# Columns: drugbank_id, generic_name, generic_components, generic_components_key,
#          synonyms, dose_synonyms, raw_doses, form_synonyms, raw_forms,
#          route_synonyms, raw_routes, salt_names, atc_codes, groups
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

clean_numeric_string <- function(value) {
  if (is.null(value) || is.na(value)) return(NA_character_)
  s <- as.character(value)
  s <- gsub("([0-9]),([0-9])", "\\1.\\2", s, perl = TRUE)
  s <- gsub(",", "", s, fixed = TRUE)
  tolower(trimws(s))
}

extract_numeric <- function(value) {
  s <- clean_numeric_string(value)
  if (is.na(s)) return(NA_real_)
  suppressWarnings(as.numeric(s))
}

mass_to_mg <- function(value, unit) {
  if (is.na(value) || is.na(unit)) return(NA_real_)
  u <- tolower(trimws(unit))
  switch(u,
         "mg" = value,
         "g" = value * 1000,
         "mcg" = value / 1000,
         "ug" = value / 1000,
         "µg" = value / 1000,
         "kg" = value * 1e6,
         value)
}

format_numeric <- function(x) {
  if (is.na(x)) return(NA_character_)
  out <- sprintf("%.6f", x)
  out <- sub("0+$", "", out, perl = TRUE)
  out <- sub("\\.$", "", out, perl = TRUE)
  if (identical(out, "-0")) out <- "0"
  out
}

PER_UNIT_MAP <- list(
  "ml" = "ml",
  "l" = "l",
  "tab" = "tablet",
  "tabs" = "tablet",
  "tablet" = "tablet",
  "tablets" = "tablet",
  "chewing gum" = "tablet",
  "cap" = "capsule",
  "caps" = "capsule",
  "capsule" = "capsule",
  "capsules" = "capsule",
  "sachet" = "sachet",
  "sachets" = "sachet",
  "drop" = "drop",
  "drops" = "drop",
  "gtt" = "drop",
  "actuation" = "actuation",
  "actuations" = "actuation",
  "spray" = "spray",
  "sprays" = "spray",
  "puff" = "puff",
  "puffs" = "puff",
  "dose" = "dose",
  "doses" = "dose",
  "application" = "application",
  "applications" = "application",
  "ampule" = "ampule",
  "ampules" = "ampule",
  "ampoule" = "ampule",
  "ampoules" = "ampule",
  "amp" = "ampule",
  "vial" = "vial",
  "vials" = "vial"
)

normalize_per_unit <- function(unit) {
  if (is.null(unit)) return(NA_character_)
  u <- tolower(trimws(unit))
  if (!nzchar(u)) return(NA_character_)
  u <- gsub("\\s+", " ", u, perl = TRUE)
  mapped <- PER_UNIT_MAP[[u]]
  if (!is.null(mapped)) return(mapped)
  mapped <- PER_UNIT_MAP[[sub("s$", "", u)]]
  if (!is.null(mapped)) return(mapped)
  u
}

normalize_dose_value <- function(value) {
  if (is.null(value) || is.na(value)) return(NA_character_)
  original <- collapse_ws(value)
  if (!nzchar(original)) return(NA_character_)
  s <- tolower(original)
  s <- gsub("µ|μ|\u00b5", "mcg", s)
  s <- gsub("microgram", "mcg", s)
  s <- gsub("milligram", "mg", s)
  s <- gsub("millilitre|milliliter", "ml", s)
  s <- gsub("litre|liter", "l", s)
  s <- gsub("cc", "ml", s, fixed = TRUE)
  s <- gsub("\\s+", " ", s, perl = TRUE)
  s <- trimws(s)
  if (!nzchar(s)) return(NA_character_)

  pct_match <- regexec("^([0-9]+(?:\\.[0-9]+)?)\\s*%$", s, perl = TRUE)
  pct_groups <- regmatches(s, pct_match)[[1]]
  if (length(pct_groups)) {
    pct_val <- extract_numeric(pct_groups[2])
    pct_str <- format_numeric(pct_val)
    if (!is.na(pct_str)) return(paste0(pct_str, "%"))
  }

  ratio_ml_match <- regexec("([0-9]+(?:\\.[0-9]+)?)\\s*(mg|g|mcg|ug)\\s*(?:/|\\s+per\\s+)\\s*([0-9]+(?:\\.[0-9]+)?)?\\s*(ml|l)\\b", s, perl = TRUE)
  ratio_ml_groups <- regmatches(s, ratio_ml_match)[[1]]
  if (length(ratio_ml_groups)) {
    strength_val <- extract_numeric(ratio_ml_groups[2])
    unit_val <- ratio_ml_groups[3]
    per_val <- ratio_ml_groups[4]
    per_val_num <- if (nzchar(per_val)) extract_numeric(per_val) else 1
    per_unit_val <- ratio_ml_groups[5]
    if (!is.na(per_val_num) && per_val_num > 0) {
      if (per_unit_val == "l") per_val_num <- per_val_num * 1000
      mg_val <- mass_to_mg(strength_val, unit_val)
      if (!is.na(mg_val)) {
        ratio_val <- mg_val / per_val_num
        ratio_str <- format_numeric(ratio_val)
        if (!is.na(ratio_str)) return(paste0(ratio_str, "mg/ml"))
      }
    }
  }

  per_unit_pattern <- "(tab(?:let)?s?|caps?(?:ule)?s?|sachets?|drops?|gtt|actuations?|sprays?|puffs?|doses?|applications?|ampoules?|ampules?|vials?)"
  ratio_unit_match <- regexec(
    paste0("([0-9]+(?:\\.[0-9]+)?)\\s*(mg|g|mcg|ug|iu)\\s*(?:/|\\s+per\\s+)\\s*([0-9]+(?:\\.[0-9]+)?)?\\s*(", per_unit_pattern, ")\\b"),
    s,
    perl = TRUE
  )
  ratio_unit_groups <- regmatches(s, ratio_unit_match)[[1]]
  if (length(ratio_unit_groups)) {
    strength_val <- extract_numeric(ratio_unit_groups[2])
    unit_val <- ratio_unit_groups[3]
    per_val <- ratio_unit_groups[4]
    per_val_num <- if (nzchar(per_val)) extract_numeric(per_val) else 1
    per_unit_val <- normalize_per_unit(ratio_unit_groups[5])
    if (!is.na(per_val_num) && !is.na(per_unit_val) && nzchar(per_unit_val)) {
      mg_val <- if (unit_val %chin% c("mg", "g", "mcg", "ug")) mass_to_mg(strength_val, unit_val) else NA_real_
      num_str <- if (!is.na(mg_val)) paste0(format_numeric(mg_val), "mg") else {
        amt_str <- format_numeric(strength_val)
        if (is.na(amt_str)) return(original)
        paste0(amt_str, unit_val)
      }
      per_part <- if (is.na(per_val_num) || per_val_num == 1) {
        per_unit_val
      } else {
        per_val_str <- format_numeric(per_val_num)
        if (is.na(per_val_str)) return(original)
        paste0(per_val_str, per_unit_val)
      }
      if (nzchar(num_str) && nzchar(per_part)) return(paste0(num_str, "/", per_part))
    }
  }

  amount_match <- regexec("([0-9]+(?:\\.[0-9]+)?)\\s*(mg|g|mcg|ug|iu)\\b", s, perl = TRUE)
  amount_groups <- regmatches(s, amount_match)[[1]]
  if (length(amount_groups)) {
    strength_val <- extract_numeric(amount_groups[2])
    unit_val <- amount_groups[3]
    if (unit_val %chin% c("mg", "g", "mcg", "ug")) {
      mg_val <- mass_to_mg(strength_val, unit_val)
      mg_str <- format_numeric(mg_val)
      if (!is.na(mg_str)) return(paste0(mg_str, "mg"))
    } else {
      amt_str <- format_numeric(strength_val)
      if (!is.na(amt_str)) return(paste0(amt_str, unit_val))
    }
  }

  volume_match <- regexec("([0-9]+(?:\\.[0-9]+)?)\\s*(ml|l)\\b", s, perl = TRUE)
  volume_groups <- regmatches(s, volume_match)[[1]]
  if (length(volume_groups)) {
    vol_val <- extract_numeric(volume_groups[2])
    vol_unit <- volume_groups[3]
    if (vol_unit == "l") vol_val <- vol_val * 1000
    vol_str <- format_numeric(vol_val)
    if (!is.na(vol_str)) return(paste0(vol_str, "ml"))
  }

  gsub("\\s+", " ", tolower(original), perl = TRUE)
}

normalize_dose_vector <- function(values) {
  if (is.null(values) || !length(values)) return(character())
  out <- unlist(lapply(values, function(v) {
    norm <- normalize_dose_value(v)
    if (is.null(norm) || is.na(norm) || !nzchar(norm)) return(character())
    norm
  }), use.names = FALSE)
  if (!length(out)) return(character())
  unique_canonical(out)
}

build_inverse_map <- function(map) {
  inv <- list()
  for (key in names(map)) {
    canon <- map[[key]]
    inv[[canon]] <- unique_canonical(c(inv[[canon]], canon, key))
  }
  inv
}

FORM_CANONICAL_MAP <- list(
  "tab" = "tablet",
  "tabs" = "tablet",
  "tablet" = "tablet",
  "tablets" = "tablet",
  "chewing gum" = "tablet",
  "cap" = "capsule",
  "caps" = "capsule",
  "capsule" = "capsule",
  "capsulee" = "capsule",
  "capsules" = "capsule",
  "susp" = "suspension",
  "suspension" = "suspension",
  "syr" = "syrup",
  "syrup" = "syrup",
  "sol" = "solution",
  "soln" = "solution",
  "solution" = "solution",
  "inhal.solution" = "solution",
  "instill.solution" = "solution",
  "lamella" = "solution",
  "ointment" = "ointment",
  "oint" = "ointment",
  "gel" = "gel",
  "cream" = "cream",
  "lotion" = "lotion",
  "patch" = "patch",
  "supp" = "suppository",
  "suppository" = "suppository",
  "dpi" = "dpi",
  "inhal.powder" = "dpi",
  "mdi" = "mdi",
  "inhal.aerosol" = "mdi",
  "oral aerosol" = "mdi",
  "ampu" = "ampule",
  "ampul" = "ampule",
  "ampule" = "ampule",
  "ampoule" = "ampule",
  "amp" = "ampule",
  "vial" = "vial",
  "inj" = "injection",
  "injection" = "injection",
  "implant" = "solution",
  "s.c. implant" = "solution",
  "metered dose inhaler" = "mdi",
  "dry powder inhaler" = "dpi",
  "spray" = "spray",
  "nasal spray" = "spray",
  "nebule" = "solution",
  "neb" = "solution",
  "inhaler" = "mdi"
)

FORM_CANONICAL_KEYS <- names(FORM_CANONICAL_MAP)[order(nchar(names(FORM_CANONICAL_MAP)), decreasing = TRUE)]
FORM_CANONICAL_VALUES <- unique(unlist(FORM_CANONICAL_MAP, use.names = FALSE))
RELEASE_PATTERN <- "(extended(?:-|\u0020)?release|immediate(?:-|\u0020)?release|delayed(?:-|\u0020)?release|sustained(?:-|\u0020)?release|controlled(?:-|\u0020)?release)"
FORM_SYNONYM_LOOKUP <- build_inverse_map(FORM_CANONICAL_MAP)

normalize_form_value <- function(value) {
  if (is.null(value) || is.na(value)) return(NA_character_)
  s <- tolower(collapse_ws(value))
  if (!nzchar(s)) return(NA_character_)
  release_match <- regexpr(RELEASE_PATTERN, s, perl = TRUE)
  release_suffix <- ""
  if (release_match[1] > -1) {
    release_suffix <- regmatches(s, list(release_match))[[1]]
    release_suffix <- gsub("-", " ", release_suffix, fixed = TRUE)
    s <- trimws(sub(RELEASE_PATTERN, "", s, perl = TRUE))
  }
  base_part <- strsplit(s, "[,;/]", perl = TRUE)[[1]]
  base_token <- if (length(base_part)) trimws(base_part[1]) else s
  canonical <- if (base_token %chin% names(FORM_CANONICAL_MAP)) FORM_CANONICAL_MAP[[base_token]] else NULL
  if (is.null(canonical)) {
    for (key in FORM_CANONICAL_KEYS) {
      if (grepl(paste0("\\b", key, "\\b"), base_token, perl = TRUE)) {
        canonical <- FORM_CANONICAL_MAP[[key]]
        break
      }
    }
  }
  if (is.null(canonical)) {
    for (value_candidate in FORM_CANONICAL_VALUES) {
      if (grepl(paste0("\\b", value_candidate, "\\b"), base_token, perl = TRUE)) {
        canonical <- value_candidate
        break
      }
    }
  }
  if (is.null(canonical) || !nzchar(canonical)) {
    canonical <- base_token
  }
  canonical <- trimws(gsub("\\s+", " ", canonical, perl = TRUE))
  if (nzchar(release_suffix)) {
    suffix <- trimws(gsub("\\s+", " ", release_suffix, perl = TRUE))
    return(paste0(canonical, ", ", suffix))
  }
  canonical
}

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b

normalize_form_vector <- function(values) {
  if (is.null(values) || !length(values)) return(character())
  out <- unlist(lapply(values, function(v) {
    norm <- normalize_form_value(v)
    if (is.null(norm) || is.na(norm) || !nzchar(norm)) return(character())
    norm
  }), use.names = FALSE)
  if (!length(out)) return(character())
  unique_canonical(out)
}

ROUTE_ALIAS_MAP <- list(
  "oral" = "oral",
  "po" = "oral",
  "per orem" = "oral",
  "per os" = "oral",
  "by mouth" = "oral",
  "iv" = "intravenous",
  "intravenous" = "intravenous",
  "im" = "intramuscular",
  "intramuscular" = "intramuscular",
  "sc" = "subcutaneous",
  "subcut" = "subcutaneous",
  "subcutaneous" = "subcutaneous",
  "subdermal" = "subcutaneous",
  "sl" = "sublingual",
  "sublingual" = "sublingual",
  "bucc" = "buccal",
  "buccal" = "buccal",
  "topical" = "topical",
  "cutaneous" = "topical",
  "dermal" = "transdermal",
  "oph" = "ophthalmic",
  "eye" = "ophthalmic",
  "ophthalmic" = "ophthalmic",
  "otic" = "otic",
  "ear" = "otic",
  "inh" = "inhalation",
  "neb" = "inhalation",
  "inhalation" = "inhalation",
  "inhaler" = "inhalation",
  "rectal" = "rectal",
  "per rectum" = "rectal",
  "pr" = "rectal",
  "vaginal" = "vaginal",
  "per vaginam" = "vaginal",
  "pv" = "vaginal",
  "intrathecal" = "intrathecal",
  "nasal" = "nasal",
  "per nasal" = "nasal",
  "intranasal" = "nasal",
  "td" = "transdermal",
  "transdermal" = "transdermal",
  "intradermal" = "intradermal",
  "id" = "intradermal",
  "urethral" = "urethral",
  "intravesical" = "intravesical",
  "s.c. implant" = "subcutaneous"
)
ROUTE_SYNONYM_LOOKUP <- build_inverse_map(ROUTE_ALIAS_MAP)

normalize_route_entry <- function(value) {
  if (is.null(value) || is.na(value)) return(character())
  s <- tolower(collapse_ws(value))
  if (!nzchar(s)) return(character())
  parts <- unique(trimws(unlist(strsplit(s, "[/|,;]", perl = TRUE), use.names = FALSE)))
  parts <- parts[nzchar(parts)]
  suggestions <- character()
  for (part in parts) {
    if (!is.null(ROUTE_ALIAS_MAP[[part]])) {
      suggestions <- c(suggestions, ROUTE_ALIAS_MAP[[part]])
    }
  }
  unique(suggestions)
}

normalize_route_vector <- function(values) {
  if (is.null(values) || !length(values)) return(character())
  out <- unlist(lapply(values, normalize_route_entry), use.names = FALSE)
  if (!length(out)) return(character())
  unique_canonical(out)
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
final_dt[, dose_synonyms := sapply(seq_len(.N), function(i) collapse_pipe(expand_dose_set(normalized_doses_list[[i]], raw_doses_list[[i]])))]
final_dt[, raw_doses := sapply(raw_doses_list, collapse_pipe)]
final_dt[, form_synonyms := sapply(seq_len(.N), function(i) collapse_pipe(expand_value_set(normalized_forms_list[[i]], raw_forms_list[[i]], FORM_SYNONYM_LOOKUP)))]
final_dt[, raw_forms := sapply(raw_forms_list, collapse_pipe)]
final_dt[, route_synonyms := sapply(seq_len(.N), function(i) collapse_pipe(expand_value_set(normalized_routes_list[[i]], raw_routes_list[[i]], ROUTE_SYNONYM_LOOKUP)))]
final_dt[, raw_routes := sapply(raw_routes_list, collapse_pipe)]
final_dt[, salt_names := sapply(seq_len(.N), function(i) collapse_pipe(expand_salt_set(salt_names_list[[i]])))]
final_dt[, atc_codes := sapply(atc_codes_list, collapse_pipe)]
final_dt[, groups := sapply(groups_list, collapse_pipe)]

final_dt <- final_dt[, .(
  drugbank_id,
  generic_name,
  generic_components,
  generic_components_key,
  synonyms,
  dose_synonyms,
  raw_doses,
  form_synonyms,
  raw_forms,
  route_synonyms,
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
PER_UNIT_SYNONYM_LOOKUP <- build_inverse_map(PER_UNIT_MAP)

expand_value_set <- function(primary, raw_values, lookup) {
  values <- combine_values(primary, raw_values)
  if (!length(primary)) return(values)
  for (val in primary) {
    if (!is.null(lookup[[val]])) {
      values <- c(values, lookup[[val]])
    }
  }
  unique_canonical(values)
}

expand_dose_set <- function(primary, raw_values) {
  values <- combine_values(primary, raw_values)
  for (val in primary) {
    if (!nzchar(val)) next
    if (grepl("/", val, fixed = TRUE)) {
      parts <- strsplit(val, "/", fixed = TRUE)[[1]]
      if (length(parts) == 2) {
        base <- parts[1]
        per_part <- trimws(parts[2])
        syns <- PER_UNIT_SYNONYM_LOOKUP[[per_part]]
        if (!is.null(syns)) {
          values <- c(values, paste0(base, "/", syns))
        }
      }
    }
  }
  unique_canonical(values)
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
