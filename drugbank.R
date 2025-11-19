# =================================
# DrugBank → aggregated ATC-level generics + brand mapping for ESOA ingestion
# Outputs:
# - drugbank_generics.csv: generic, dose, form, route, atc_code, drugbank_ids
# - drugbank_brands.csv: brand, generic
# =================================

ensure_installed <- function(packages, repos = "https://cloud.r-project.org") {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = repos)
    }
  }
}

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  getwd()
}

# --- deps ---
ensure_installed(c("data.table", "arrow"))
if (!requireNamespace("dbdataset", quietly = TRUE)) {
  ensure_installed("remotes")
  remotes::install_github("interstellar-Consultation-Services/dbdataset", quiet = TRUE, upgrade = "never")
}

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(dbdataset)
})

collapse_ws <- function(x) {
  x <- gsub("\\s+", " ", x, perl = TRUE)
  trimws(x)
}

empty_to_na <- function(x) {
  x <- trimws(x)
  x[!nzchar(x)] <- NA_character_
  x
}

unique_canonical <- function(x) {
  if (!length(x)) return(character())
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  if (!length(x)) return(character())
  ord <- order(tolower(x), x)
  x <- x[ord]
  x[!duplicated(tolower(x))]
}

combine_values <- function(...) {
  vals <- unlist(list(...), use.names = FALSE)
  unique_canonical(vals)
}

write_arrow_csv <- function(dt, path) {
  cat(sprintf("Writing %s rows to %s\n", format(nrow(dt), big.mark = ","), path))
  arrow_table <- arrow::Table$create(dt)
  if ("write_csv_arrow" %in% getNamespaceExports("arrow")) {
    arrow::write_csv_arrow(arrow_table, path)
  } else if ("write_delim_arrow" %in% getNamespaceExports("arrow")) {
    arrow::write_delim_arrow(arrow_table, path, delim = ",")
  } else {
    warning("Falling back to data.table::fwrite because Arrow CSV writers are unavailable.")
    fwrite(dt, path)
  }
  cat("Done.\n")
}

split_ingredients <- function(value) {
  if (is.null(value) || length(value) == 0L) return(character())
  parts <- unlist(lapply(value, function(v) {
    if (is.na(v) || !nzchar(v)) return(character())
    v <- collapse_ws(v)
    if (!nzchar(v)) return(character())
    v <- gsub("\\s+(and|with|plus)\\s+", "|", v, ignore.case = TRUE, perl = TRUE)
    v <- gsub("\\s*/\\s*|\\s*\\+\\s*|\\s*,\\s*|\\s*;\\s*", "|", v, perl = TRUE)
    res <- strsplit(v, "\\|", perl = TRUE)[[1]]
    res <- trimws(res)
    res[nzchar(res)]
  }), use.names = FALSE)
  unique_canonical(parts)
}

split_semicolon <- function(value) {
  if (is.null(value) || length(value) == 0L) return(character())
  parts <- unlist(lapply(value, function(v) {
    if (is.na(v) || !nzchar(v)) return(character())
    res <- strsplit(v, ";", fixed = TRUE)[[1]]
    res <- trimws(res)
    res[nzchar(res)]
  }), use.names = FALSE)
  unique_canonical(parts)
}

normalize_list_column <- function(col) {
  lapply(col, function(item) {
    if (is.null(item)) return(character())
    if (length(item) == 1L && is.na(item)) return(character())
    item <- unlist(item, use.names = FALSE)
    if (!length(item)) return(character())
    item <- trimws(as.character(item))
    item <- item[nzchar(item)]
    unique_canonical(item)
  })
}

ensure_list_column <- function(DT, col) {
  if (!col %in% names(DT)) {
    DT[, (col) := vector("list", .N)]
  }
  DT[, (col) := normalize_list_column(.SD[[1]]), .SDcols = col]
}

combine_list_column <- function(lst) {
  if (is.null(lst) || !length(lst)) return(character())
  combine_values(unlist(lst, use.names = FALSE))
}

collapse_pipe <- function(values) {
  values <- combine_values(values)
  if (!length(values)) {
    return(NA_character_)
  }
  paste(values, collapse = "|")
}

format_numeric <- function(x) {
  if (is.null(x) || is.na(x) || is.nan(x) || is.infinite(x)) return(NA_character_)
  out <- sprintf("%.6f", x)
  out <- sub("0+$", "", out, perl = TRUE)
  out <- sub("\\.$", "", out, perl = TRUE)
  if (identical(out, "-0")) out <- "0"
  out
}

mass_to_mg <- function(value, unit) {
  if (is.null(value) || is.na(value) || is.null(unit) || is.na(unit)) return(NA_real_)
  u <- tolower(trimws(unit))
  switch(u,
    "mg" = value,
    "g" = value * 1000,
    "mcg" = value / 1000,
    "ug" = value / 1000,
    NA_real_
  )
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
  "ampu" = "ampule",
  "ampul" = "ampule",
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

clean_numeric_string <- function(value) {
  if (is.null(value) || is.na(value)) return(NA_character_)
  s <- as.character(value)
  s <- gsub("([0-9]),([0-9])", "\\1.\\2", s, perl = TRUE)
  s <- gsub(",", "", s, fixed = TRUE)
  trimws(s)
}

extract_numeric <- function(value) {
  s <- clean_numeric_string(value)
  if (is.na(s)) return(NA_real_)
  suppressWarnings(as.numeric(s))
}

normalize_dose_value <- function(value) {
  if (is.null(value) || is.na(value)) return(NA_character_)
  original <- collapse_ws(value)
  if (!nzchar(original)) return(NA_character_)
  s <- tolower(original)
  s <- gsub("\u00b5|µ|μ", "mcg", s, perl = TRUE)
  s <- gsub("microgram", "mcg", s, fixed = TRUE)
  s <- gsub("milligram", "mg", s, fixed = TRUE)
  s <- gsub("millilitre", "ml", s, fixed = TRUE)
  s <- gsub("milliliter", "ml", s, fixed = TRUE)
  s <- gsub("litre", "l", s, fixed = TRUE)
  s <- gsub("liter", "l", s, fixed = TRUE)
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

FORM_TO_ROUTE_CANONICAL <- list(
  "tablet" = "oral", "tab" = "oral", "tabs" = "oral", "chewing gum" = "oral",
  "capsule" = "oral", "cap" = "oral", "caps" = "oral", "capsules" = "oral",
  "syrup" = "oral", "suspension" = "oral", "susp" = "oral", "solution" = "oral", "soln" = "oral", "sol" = "oral",
  "sachet" = "oral",
  "drop" = "ophthalmic", "eye drop" = "ophthalmic", "ear drop" = "otic",
  "cream" = "topical", "ointment" = "topical", "gel" = "topical", "lotion" = "topical",
  "patch" = "transdermal",
  "inhaler" = "inhalation", "nebule" = "inhalation", "neb" = "inhalation",
  "inhal.aerosol" = "inhalation", "inhal.powder" = "inhalation", "inhal.solution" = "inhalation", "oral aerosol" = "inhalation",
  "ampoule" = "intravenous", "amp" = "intravenous", "ampul" = "intravenous", "ampule" = "intravenous", "vial" = "intravenous",
  "inj" = "intravenous", "injection" = "intravenous",
  "suppository" = "rectal", "supp" = "rectal",
  "mdi" = "inhalation", "dpi" = "inhalation",
  "metered dose inhaler" = "inhalation",
  "dry powder inhaler" = "inhalation",
  "spray" = "nasal", "nasal spray" = "nasal",
  "td" = "transdermal",
  "instill.solution" = "ophthalmic", "lamella" = "ophthalmic",
  "implant" = "subcutaneous", "s.c. implant" = "subcutaneous"
)

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

ROUTE_ALIAS_KEYS <- names(ROUTE_ALIAS_MAP)[order(nchar(names(ROUTE_ALIAS_MAP)), decreasing = TRUE)]
ALLOWED_ROUTE_SET <- sort(unique(c(unlist(ROUTE_ALIAS_MAP, use.names = FALSE), unlist(FORM_TO_ROUTE_CANONICAL, use.names = FALSE))))

normalize_route_entry <- function(value) {
  if (is.null(value) || is.na(value)) return(character())
  s <- tolower(collapse_ws(value))
  if (!nzchar(s)) return(character())
  tokens <- character()
  for (alias in ROUTE_ALIAS_KEYS) {
    if (grepl(paste0("\\b", alias, "\\b"), s, perl = TRUE)) {
      tokens <- c(tokens, ROUTE_ALIAS_MAP[[alias]])
    }
  }
  parts <- unlist(strsplit(s, "[/|,;]", perl = TRUE), use.names = FALSE)
  if (!length(parts)) parts <- c(s)
  for (part in parts) {
    part_trim <- trimws(part)
    if (!nzchar(part_trim)) next
    if (!is.null(ROUTE_ALIAS_MAP[[part_trim]])) {
      tokens <- c(tokens, ROUTE_ALIAS_MAP[[part_trim]])
      next
    }
    if (!is.null(FORM_TO_ROUTE_CANONICAL[[part_trim]])) {
      tokens <- c(tokens, FORM_TO_ROUTE_CANONICAL[[part_trim]])
      next
    }
    words <- strsplit(part_trim, " ", fixed = TRUE)[[1]]
    matched <- FALSE
    for (word in words) {
      word_trim <- trimws(word)
      if (!nzchar(word_trim)) next
      if (!is.null(ROUTE_ALIAS_MAP[[word_trim]])) {
        tokens <- c(tokens, ROUTE_ALIAS_MAP[[word_trim]])
        matched <- TRUE
      } else if (!is.null(FORM_TO_ROUTE_CANONICAL[[word_trim]])) {
        tokens <- c(tokens, FORM_TO_ROUTE_CANONICAL[[word_trim]])
        matched <- TRUE
      } else if (word_trim %chin% ALLOWED_ROUTE_SET) {
        tokens <- c(tokens, word_trim)
        matched <- TRUE
      }
    }
    if (!matched && part_trim %chin% ALLOWED_ROUTE_SET) {
      tokens <- c(tokens, part_trim)
    }
  }
  tokens <- unique(tokens)
  if (!length(tokens)) return(c(s))
  tokens
}

normalize_route_vector <- function(values) {
  if (is.null(values) || !length(values)) return(character())
  out <- unlist(lapply(values, normalize_route_entry), use.names = FALSE)
  if (!length(out)) return(character())
  unique_canonical(out)
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_generics_path <- file.path(output_dir, "drugbank_generics.csv")
output_brands_path <- file.path(output_dir, "drugbank_brands.csv")

paths_equal <- function(path_a, path_b) {
  if (is.null(path_a) || is.null(path_b)) {
    return(FALSE)
  }
  a_norm <- tryCatch(
    normalizePath(path_a, winslash = "/", mustWork = FALSE),
    error = function(...) NA_character_
  )
  b_norm <- tryCatch(
    normalizePath(path_b, winslash = "/", mustWork = FALSE),
    error = function(...) NA_character_
  )
  !is.na(a_norm) && !is.na(b_norm) && identical(a_norm, b_norm)
}

safe_copy <- function(src, dest) {
  tryCatch({
    if (!file.exists(src) || paths_equal(src, dest)) {
      return(FALSE)
    }
    dest_dir <- dirname(dest)
    if (!dir.exists(dest_dir)) {
      dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(dest_dir)) {
      return(FALSE)
    }
    file.copy(src, dest, overwrite = TRUE, copy.mode = TRUE)
  }, error = function(...) FALSE)
}

copy_outputs_to_superproject <- function(src_file) {
  repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
  dependencies_dir <- file.path(repo_root, "dependencies")
  if (!dir.exists(dependencies_dir)) {
    return(invisible(FALSE))
  }

  super_output_dir <- file.path(repo_root, "dependencies", "drugbank_generics", "output")
  safe_copy(src_file, file.path(super_output_dir, basename(src_file)))

  inputs_dir <- file.path(repo_root, "inputs", "drugs")
  safe_copy(src_file, file.path(inputs_dir, basename(src_file)))
}

dataset <- drugbank

excluded_groups <- c("experimental", "withdrawn", "illicit", "vet")
groups_dt <- as.data.table(dataset$drugs$groups)
groups_dt[, group_clean := tolower(trimws(group))]
excluded_ids <- unique(groups_dt[group_clean %chin% excluded_groups, as.character(drugbank_id)])

general_dt <- as.data.table(dataset$drugs$general_information)[
  , .(
    drugbank_id = as.character(drugbank_id),
    generic_name = collapse_ws(name)
  )
]
general_dt <- general_dt[!is.na(generic_name) & nzchar(generic_name)]
if (length(excluded_ids)) general_dt <- general_dt[!(drugbank_id %chin% excluded_ids)]

generic_dt <- copy(general_dt)
generic_dt[, generic_parts := lapply(generic_name, split_ingredients)]

ensure_list_column(generic_dt, "generic_parts")
generic_dt[, generic := vapply(generic_parts, function(parts) {
  if (!length(parts)) return(NA_character_)
  paste(parts, collapse = "; ")
}, character(1))]

products_dt <- as.data.table(dataset$products)[
  , .(
    drugbank_id = as.character(drugbank_id),
    product_name = collapse_ws(name),
    route = collapse_ws(route),
    dosage_form = collapse_ws(dosage_form),
    strength = collapse_ws(strength)
  )
]
if (length(excluded_ids)) products_dt <- products_dt[!(drugbank_id %chin% excluded_ids)]
products_dt[, c("product_name", "route", "dosage_form", "strength") := lapply(
  .SD, empty_to_na
), .SDcols = c("product_name", "route", "dosage_form", "strength")]

product_agg <- products_dt[
  ,
  {
    names_vec <- unique_canonical(product_name[!is.na(product_name)])
    route_vec <- unique_canonical(unlist(lapply(route[!is.na(route)], split_semicolon), use.names = FALSE))
    form_vec <- unique_canonical(unlist(lapply(dosage_form[!is.na(dosage_form)], split_semicolon), use.names = FALSE))
    dose_vec <- unique_canonical(strength[!is.na(strength)])
    list(
      product_names = list(names_vec),
      product_routes = list(route_vec),
      product_forms = list(form_vec),
      product_doses = list(dose_vec)
    )
  },
  by = drugbank_id
]

dosages_dt <- as.data.table(dataset$drugs$dosages)[
  , .(
    drugbank_id = as.character(drugbank_id),
    route = collapse_ws(route),
    dosage_form = collapse_ws(form),
    strength = collapse_ws(strength)
  )
]
if (length(excluded_ids)) dosages_dt <- dosages_dt[!(drugbank_id %chin% excluded_ids)]
dosages_dt[, c("route", "dosage_form", "strength") := lapply(
  .SD, empty_to_na
), .SDcols = c("route", "dosage_form", "strength")]

dosage_agg <- dosages_dt[
  ,
  {
    route_vec <- unique_canonical(unlist(lapply(route[!is.na(route)], split_semicolon), use.names = FALSE))
    form_vec <- unique_canonical(unlist(lapply(dosage_form[!is.na(dosage_form)], split_semicolon), use.names = FALSE))
    dose_vec <- unique_canonical(strength[!is.na(strength)])
    list(
      dosage_routes = list(route_vec),
      dosage_forms = list(form_vec),
      dosage_doses = list(dose_vec)
    )
  },
  by = drugbank_id
]

brands_dt <- as.data.table(dataset$drugs$international_brands)[
  , .(
    drugbank_id = as.character(drugbank_id),
    brand_name = collapse_ws(brand)
  )
]
if (length(excluded_ids)) brands_dt <- brands_dt[!(drugbank_id %chin% excluded_ids)]
brands_dt[, brand_name := empty_to_na(brand_name)]

brand_agg <- brands_dt[
  !is.na(brand_name),
  .(brand_names = list(unique_canonical(brand_name))),
  by = drugbank_id
]

atc_dt <- as.data.table(dataset$drugs$atc_codes)[
  , .(
    drugbank_id = as.character(drugbank_id),
    atc_code = trimws(atc_code)
  )
]
atc_dt <- atc_dt[!is.na(atc_code) & nzchar(atc_code)]
if (length(excluded_ids)) atc_dt <- atc_dt[!(drugbank_id %chin% excluded_ids)]
atc_dt <- unique(atc_dt, by = c("drugbank_id", "atc_code"))

drug_ids <- unique(atc_dt$drugbank_id)
info_dt <- data.table(drugbank_id = drug_ids)

info_dt <- merge(info_dt, generic_dt[, .(drugbank_id, generic, generic_parts)], by = "drugbank_id", all.x = TRUE, sort = FALSE)
info_dt <- merge(info_dt, product_agg, by = "drugbank_id", all.x = TRUE, sort = FALSE)
info_dt <- merge(info_dt, dosage_agg, by = "drugbank_id", all.x = TRUE, sort = FALSE)
info_dt <- merge(info_dt, brand_agg, by = "drugbank_id", all.x = TRUE, sort = FALSE)

ensure_list_column(info_dt, "generic_parts")
ensure_list_column(info_dt, "product_names")
ensure_list_column(info_dt, "product_routes")
ensure_list_column(info_dt, "product_forms")
ensure_list_column(info_dt, "product_doses")
ensure_list_column(info_dt, "dosage_routes")
ensure_list_column(info_dt, "dosage_forms")
ensure_list_column(info_dt, "dosage_doses")
ensure_list_column(info_dt, "brand_names")

brands_output <- info_dt[
  ,
  {
    brand_vec <- brand_names[[1]]
    if (!length(brand_vec)) {
      NULL
    } else {
      generic_vec <- unique_canonical(generic_parts[[1]])
      if (!length(generic_vec)) {
        generic_vec <- combine_values(generic)
      }
      generic_str <- if (length(generic_vec)) paste(generic_vec, collapse = "; ") else NA_character_
      data.table(
        brand = brand_vec,
        generic = generic_str
      )
    }
  },
  by = drugbank_id
]
if (!is.null(brands_output) && nrow(brands_output)) {
  brands_output <- unique(brands_output, by = c("brand", "generic"))
  setorder(brands_output, brand, generic)
} else {
  brands_output <- data.table(brand = character(), generic = character())
}

info_dt[, route_values := Map(function(prod_route, dosage_route) {
  combine_values(prod_route, dosage_route)
}, product_routes, dosage_routes)]

info_dt[, form_values := Map(function(prod_form, dosage_form) {
  combine_values(prod_form, dosage_form)
}, product_forms, dosage_forms)]

info_dt[, dose_values := Map(function(prod_dose, dosage_dose) {
  combine_values(prod_dose, dosage_dose)
}, product_doses, dosage_doses)]

info_dt[, dose_values := lapply(dose_values, normalize_dose_vector)]
info_dt[, form_values := lapply(form_values, normalize_form_vector)]
info_dt[, route_values := lapply(route_values, normalize_route_vector)]

info_dt[, c(
  "product_names", "brand_names",
  "product_routes", "dosage_routes",
  "product_forms", "dosage_forms",
  "product_doses", "dosage_doses"
) := NULL]

drug_atc <- merge(atc_dt, info_dt, by = "drugbank_id", all.x = TRUE, sort = FALSE)

generics_dt <- drug_atc[
  ,
  {
    generic_vec <- combine_list_column(generic_parts)
    if (!length(generic_vec)) {
      generic_vec <- combine_values(generic)
    }
    generic_str <- if (length(generic_vec)) paste(generic_vec, collapse = "; ") else NA_character_

    doses <- combine_list_column(dose_values)
    forms <- combine_list_column(form_values)
    routes <- combine_list_column(route_values)
    ids <- combine_values(drugbank_id)

    list(
      generic = generic_str,
      dose = list(doses),
      form = list(forms),
      route = list(routes),
      drugbank_ids = list(ids)
    )
  },
  by = atc_code
]

generics_dt[, dose := vapply(dose, collapse_pipe, character(1))]
generics_dt[, form := vapply(form, collapse_pipe, character(1))]
generics_dt[, route := vapply(route, collapse_pipe, character(1))]
generics_dt[, drugbank_ids := vapply(drugbank_ids, collapse_pipe, character(1))]

setcolorder(generics_dt, c("generic", "dose", "form", "route", "atc_code", "drugbank_ids"))
setorder(generics_dt, atc_code)

if (nrow(generics_dt) > 100000) {
  message(sprintf("Generics output currently has %s rows (>100k)", format(nrow(generics_dt), big.mark = ",")))
}

write_arrow_csv(generics_dt, output_generics_path)
write_arrow_csv(brands_output, output_brands_path)
copy_outputs_to_superproject(output_generics_path)
copy_outputs_to_superproject(output_brands_path)
