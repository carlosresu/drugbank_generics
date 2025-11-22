#!/usr/bin/env Rscript
# drugbank_generics.R — build a generics-focused DrugBank master dataset
# Columns: drugbank_id, lexeme, canonical_generic_name, generic_components,
#          generic_components_key, dose_norm, raw_dose, dose_synonyms,
#          form_norm, raw_form, form_synonyms, route_norm, raw_route,
#          route_synonyms, atc_code, salt_names, groups

suppressWarnings({
  suppressPackageStartupMessages({
    ensure_installed <- function(pkg, repos = c("https://community.r-multiverse.org", "https://cloud.r-project.org"), min_version = NULL) {
      needs_install <- !requireNamespace(pkg, quietly = TRUE)
      if (!needs_install && !is.null(min_version)) {
        needs_install <- utils::packageVersion(pkg) < min_version
      }
      if (needs_install) {
        install.packages(pkg, repos = repos)
      }
    }
    install_polars <- function() {
      ensure_installed(
        "polars",
        repos = c("https://community.r-multiverse.org", "https://cloud.r-project.org"),
        min_version = "1.6.0"
      )
    }
    install_polars()
    ensure_installed("dbdataset")
  })
})

suppressPackageStartupMessages({
  library(polars)
  library(dbdataset)
})

pl <- polars::pl

to_polars_df <- function(x) {
  if (is_polars_df(x)) return(x)
  if (is_polars_lf(x)) return(x$collect())
  if (is.data.frame(x)) {
    return(as_polars_df(x))
  }
  if (is.list(x)) {
    list_cols <- vapply(x, is.list, logical(1))
    if (any(list_cols)) {
      x[list_cols] <- lapply(x[list_cols], I)
    }
    return(as_polars_df(x))
  }
  as_polars_df(list(value = x))
}

argv <- commandArgs(trailingOnly = TRUE)
keep_all_flag <- "--keep-all" %in% argv
quiet_mode <- identical(tolower(Sys.getenv("ESOA_DRUGBANK_QUIET", "0")), "1")

worker_env <- suppressWarnings(as.integer(Sys.getenv("ESOA_DRUGBANK_WORKERS", "")))
if (!is.na(worker_env) && worker_env > 0) {
  Sys.setenv(POLARS_MAX_THREADS = worker_env)
}

collapse_ws <- function(x) {
  ifelse(is.na(x), NA_character_, trimws(gsub("\\s+", " ", as.character(x))))
}

empty_to_na <- function(x) {
  val <- collapse_ws(x)
  ifelse(!nzchar(val), NA_character_, val)
}

split_coder_tokens <- function(value) {
  if (is.na(value)) return(character())
  cleaned <- trimws(value)
  if (!nzchar(cleaned)) return(character())
  tokens <- unlist(strsplit(cleaned, "[,\\/]", perl = TRUE), use.names = FALSE)
  tokens <- tolower(trimws(tokens))
  tokens[tokens != ""]
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

`%nin%` <- function(x, table) !(x %in% table)

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
  parts <- unlist(strsplit(val, "(?i)\\sand\\s|\\swith\\s|\\splus\\s|\\+|/|,(?=\\s)|;", perl = TRUE))
  parts <- collapse_ws(parts)
  parts <- parts[nzchar(parts)]
  unique_canonical(parts)
}

strip_parenthetical_segments <- function(value) {
  original <- collapse_ws(value)
  if (is.null(original) || is.na(original) || !nzchar(original)) {
    return(list(base = NA_character_, details = NA_character_))
  }
  details <- character()
  val <- original
  repeat {
    match <- regexpr("(?<!\\S)\\([^()]*\\)(?=\\s|$|[;,:/])", val, perl = TRUE)
    if (match[1] == -1) break
    seg <- substr(val, match[1], match[1] + attr(match, "match.length") - 1)
    detail <- trimws(substr(seg, 2, nchar(seg) - 1))
    if (nzchar(detail)) details <- c(details, detail)
    val <- paste0(
      substr(val, 1, match[1] - 1),
      substr(val, match[1] + attr(match, "match.length"), nchar(val))
    )
  }
  val <- gsub("\\s+", " ", val, perl = TRUE)
  val <- trimws(val)
  list(base = if (nzchar(val)) val else NA_character_, details = details)
}

strip_comma_segments <- function(value) {
  if (is.null(value) || is.na(value) || !nzchar(value)) return(list(base = value, detail = NA_character_))
  match <- regexpr(",\\s+.+$", value, perl = TRUE)
  if (match[1] == -1) return(list(base = value, detail = NA_character_))
  detail <- trimws(substr(value, match[1] + 1, nchar(value)))
  base <- trimws(substr(value, 1, match[1] - 1))
  list(
    base = if (nzchar(base)) base else NA_character_,
    detail = if (nzchar(detail)) detail else NA_character_
  )
}

clean_form_route_entry <- function(value) {
  parenthetical <- strip_parenthetical_segments(value)
  base <- parenthetical$base
  details <- parenthetical$details
  comma_clean <- strip_comma_segments(base)
  base <- comma_clean$base
  if (!is.na(comma_clean$detail)) details <- c(details, comma_clean$detail)
  detail_out <- if (length(details)) paste(details[details != ""], collapse = " | ") else NA_character_
  list(
    base = if (!is.null(base)) base else NA_character_,
    detail = if (!is.null(detail_out) && nzchar(detail_out)) detail_out else NA_character_
  )
}

clean_form_route_base <- function(value) {
  clean_form_route_entry(value)$base
}

clean_form_route_detail <- function(value) {
  clean_form_route_entry(value)$detail
}

collapse_pipe <- function(values) {
  vals <- unique_canonical(values)
  if (!length(vals)) return(NA_character_)
  paste(vals, collapse = "|")
}

canonical_case_key <- function(value) {
  val <- collapse_ws(value)
  if (is.na(val) || !nzchar(val)) return(NA_character_)
  tolower(val)
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
    "\\u00b5g" = value / 1000,
    "kg" = value * 1e6,
    value
  )
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
  "ml" = "ml", "l" = "l",
  "tab" = "tablet", "tabs" = "tablet", "tablet" = "tablet", "tablets" = "tablet",
  "chewing gum" = "tablet",
  "cap" = "capsule", "caps" = "capsule", "capsule" = "capsule", "capsules" = "capsule",
  "sachet" = "sachet", "sachets" = "sachet",
  "drop" = "drop", "drops" = "drop", "gtt" = "drop",
  "actuation" = "actuation", "actuations" = "actuation",
  "spray" = "spray", "sprays" = "spray",
  "puff" = "puff", "puffs" = "puff",
  "dose" = "dose", "doses" = "dose",
  "application" = "application", "applications" = "application",
  "ampule" = "ampule", "ampules" = "ampule", "ampoule" = "ampule", "ampoules" = "ampule",
  "amp" = "ampule", "vial" = "vial", "vials" = "vial"
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
  s <- gsub("µ|μ|\\u00b5", "mcg", s)
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
      mg_val <- if (unit_val %in% c("mg", "g", "mcg", "ug")) mass_to_mg(strength_val, unit_val) else NA_real_
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
    if (unit_val %in% c("mg", "g", "mcg", "ug")) {
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
  "tab" = "tablet", "tabs" = "tablet", "tablet" = "tablet", "tablets" = "tablet",
  "chewing gum" = "tablet",
  "cap" = "capsule", "caps" = "capsule", "capsule" = "capsule", "capsulee" = "capsule", "capsules" = "capsule",
  "susp" = "suspension", "suspension" = "suspension",
  "syr" = "syrup", "syrup" = "syrup",
  "sol" = "solution", "soln" = "solution", "solution" = "solution",
  "inhal.solution" = "solution", "instill.solution" = "solution", "lamella" = "solution",
  "ointment" = "ointment", "oint" = "ointment",
  "gel" = "gel", "cream" = "cream", "lotion" = "lotion",
  "patch" = "patch", "supp" = "suppository", "suppository" = "suppository",
  "dpi" = "dpi", "inhal.powder" = "dpi",
  "mdi" = "mdi", "inhal.aerosol" = "mdi", "oral aerosol" = "mdi",
  "ampu" = "ampule", "ampul" = "ampule", "ampule" = "ampule", "ampoule" = "ampule", "amp" = "ampule",
  "vial" = "vial", "inj" = "injection", "injection" = "injection",
  "implant" = "solution", "s.c. implant" = "solution",
  "metered dose inhaler" = "mdi", "dry powder inhaler" = "dpi",
  "spray" = "spray", "nasal spray" = "spray",
  "nebule" = "solution", "neb" = "solution", "inhaler" = "mdi"
)

FORM_CANONICAL_KEYS <- names(FORM_CANONICAL_MAP)[order(nchar(names(FORM_CANONICAL_MAP)), decreasing = TRUE)]
FORM_CANONICAL_VALUES <- unique(unlist(FORM_CANONICAL_MAP, use.names = FALSE))
RELEASE_PATTERN <- "(extended(?:-|\\s)?release|immediate(?:-|\\s)?release|delayed(?:-|\\s)?release|sustained(?:-|\\s)?release|controlled(?:-|\\s)?release)"

FORM_TO_ROUTE_CANONICAL <- list(
  "tablet" = "oral", "tab" = "oral", "tabs" = "oral", "chewing gum" = "oral",
  "capsule" = "oral", "cap" = "oral", "caps" = "oral", "capsulee" = "oral",
  "syrup" = "oral", "suspension" = "oral", "susp" = "oral", "solution" = "oral", "soln" = "oral", "sol" = "oral",
  "sachet" = "oral",
  "drop" = "ophthalmic", "eye drop" = "ophthalmic", "ear drop" = "otic",
  "cream" = "topical", "ointment" = "topical", "gel" = "topical", "lotion" = "topical",
  "patch" = "transdermal",
  "inhaler" = "inhalation", "nebule" = "inhalation", "neb" = "inhalation",
  "inhal.aerosol" = "inhalation", "inhal.powder" = "inhalation", "inhal.solution" = "inhalation", "oral aerosol" = "inhalation",
  "ampoule" = c("intravenous", "intravesical"), "amp" = c("intravenous", "intravesical"), "ampul" = c("intravenous", "intravesical"), "ampule" = c("intravenous", "intravesical"), "vial" = c("intravenous", "intravesical"),
  "inj" = c("intravenous", "intravesical"), "injection" = c("intravenous", "intravesical"),
  "suppository" = "rectal", "supp" = "rectal",
  "mdi" = "inhalation", "dpi" = "inhalation",
  "metered dose inhaler" = "inhalation",
  "dry powder inhaler" = "inhalation",
  "spray" = "nasal", "nasal spray" = "nasal",
  "td" = "transdermal",
  "instill.solution" = "ophthalmic", "lamella" = "ophthalmic",
  "implant" = "subcutaneous", "s.c. implant" = "subcutaneous"
)

normalize_form_value <- function(value) {
  if (is.null(value) || is.na(value)) return(character())
  s <- tolower(collapse_ws(value))
  if (!nzchar(s)) return(character())
  raw_parts <- strsplit(s, "[;|/]", perl = TRUE)[[1]]
  if (!length(raw_parts)) raw_parts <- c(s)
  canon_values <- character()
  for (part in raw_parts) {
    token <- trimws(part)
    if (!nzchar(token)) next
    release_match <- regexpr(RELEASE_PATTERN, token, perl = TRUE)
    release_suffix <- ""
    if (release_match[1] > -1) {
      release_suffix <- regmatches(token, list(release_match))[[1]]
      release_suffix <- gsub("-", " ", release_suffix, fixed = TRUE)
      token <- trimws(sub(RELEASE_PATTERN, "", token, perl = TRUE))
    }
    base_part <- strsplit(token, "[,]", perl = TRUE)[[1]]
    base_token <- if (length(base_part)) trimws(base_part[1]) else token
    canonical <- if (base_token %in% names(FORM_CANONICAL_MAP)) FORM_CANONICAL_MAP[[base_token]] else NULL
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
      canonical <- paste0(canonical, ", ", suffix)
    }
    if (nzchar(canonical)) canon_values <- c(canon_values, canonical)
  }
  unique_canonical(canon_values)
}

normalize_form_vector <- function(values) {
  if (is.null(values) || !length(values)) return(character())
  out <- unlist(lapply(values, function(v) {
    norm <- normalize_form_value(v)
    if (is.null(norm) || length(norm) == 0) return(character())
    norm
  }), use.names = FALSE)
  if (!length(out)) return(character())
  unique_canonical(out)
}

ROUTE_ALIAS_MAP <- list(
  "oral" = "oral", "po" = "oral", "per orem" = "oral", "per os" = "oral", "by mouth" = "oral",
  "iv" = "intravenous", "intravenous" = "intravenous",
  "im" = "intramuscular", "intramuscular" = "intramuscular",
  "sc" = "subcutaneous", "subcut" = "subcutaneous", "subcutaneous" = "subcutaneous", "subdermal" = "subcutaneous",
  "sl" = "sublingual", "sublingual" = "sublingual",
  "bucc" = "buccal", "buccal" = "buccal",
  "topical" = "topical", "cutaneous" = "topical",
  "dermal" = "transdermal", "td" = "transdermal", "transdermal" = "transdermal",
  "oph" = "ophthalmic", "eye" = "ophthalmic", "ophthalmic" = "ophthalmic",
  "otic" = "otic", "ear" = "otic",
  "inh" = "inhalation", "neb" = "inhalation", "inhalation" = "inhalation", "inhaler" = "inhalation",
  "rectal" = "rectal", "per rectum" = "rectal", "pr" = "rectal",
  "vaginal" = "vaginal", "per vaginam" = "vaginal", "pv" = "vaginal",
  "intrathecal" = "intrathecal", "intranasal" = "nasal", "nasal" = "nasal", "per nasal" = "nasal",
  "intradermal" = "intradermal", "id" = "intradermal",
  "urethral" = "urethral", "intravesical" = "intravesical",
  "endotracheal" = "endotracheal",
  "s.c. implant" = "subcutaneous"
)

ROUTE_ALIAS_KEYS <- names(ROUTE_ALIAS_MAP)[order(nchar(names(ROUTE_ALIAS_MAP)), decreasing = TRUE)]
ALLOWED_ROUTE_SET <- sort(unique(c(
  unlist(ROUTE_ALIAS_MAP, use.names = FALSE),
  unlist(FORM_TO_ROUTE_CANONICAL, use.names = FALSE)
)))

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
    } else if (!is.null(FORM_TO_ROUTE_CANONICAL[[part_trim]])) {
      tokens <- c(tokens, FORM_TO_ROUTE_CANONICAL[[part_trim]])
    } else if (part_trim %in% ALLOWED_ROUTE_SET) {
      tokens <- c(tokens, part_trim)
    } else {
      words <- strsplit(part_trim, " ", fixed = TRUE)[[1]]
      for (word in words) {
        word_trim <- trimws(word)
        if (!nzchar(word_trim)) next
        if (!is.null(ROUTE_ALIAS_MAP[[word_trim]])) {
          tokens <- c(tokens, ROUTE_ALIAS_MAP[[word_trim]])
        } else if (!is.null(FORM_TO_ROUTE_CANONICAL[[word_trim]])) {
          tokens <- c(tokens, FORM_TO_ROUTE_CANONICAL[[word_trim]])
        } else if (word_trim %in% ALLOWED_ROUTE_SET) {
          tokens <- c(tokens, word_trim)
        }
      }
    }
  }
  tokens <- unique(tokens)
  if (!length(tokens)) return(character())
  tokens
}

normalize_route_vector <- function(values) {
  if (is.null(values) || !length(values)) return(character())
  out <- unlist(lapply(values, normalize_route_entry), use.names = FALSE)
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

PER_UNIT_SYNONYM_LOOKUP <- build_inverse_map(PER_UNIT_MAP)
FORM_SYNONYM_LOOKUP <- build_inverse_map(FORM_CANONICAL_MAP)
ROUTE_SYNONYM_LOOKUP <- build_inverse_map(ROUTE_ALIAS_MAP)

expand_value_set <- function(primary, raw_values, lookup) {
  values <- combine_values(primary, raw_values)
  if (!length(primary)) return(values)
  for (val in primary) {
    if (!is.na(val) && nzchar(val) && !is.null(lookup[[val]])) {
      values <- c(values, lookup[[val]])
    }
  }
  unique_canonical(values)
}

expand_dose_set <- function(primary, raw_values) {
  values <- combine_values(primary, raw_values)
  for (val in primary) {
    if (is.na(val) || !nzchar(val)) next
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

is_route_compatible <- function(form_val, route_val) {
  if (is.null(form_val) || is.na(form_val) || !nzchar(form_val)) return(TRUE)
  if (is.null(route_val) || is.na(route_val) || !nzchar(route_val)) return(TRUE)
  expected <- FORM_TO_ROUTE_CANONICAL[[form_val]]
  if (is.null(expected)) return(TRUE)
  route_val %in% expected
}

collect_disallowed_pairs <- function(df) {
  if (!nrow(df)) return(data.frame())
  rows <- list()
  idx <- 1
  for (i in seq_len(nrow(df))) {
    routes <- unique_canonical(df$route_norm_list[[i]])
    forms <- unique_canonical(df$form_norm_list[[i]])
    if (!length(routes) || all(is.na(routes)) || !length(forms) || all(is.na(forms))) next
    for (r in routes) {
      for (f in forms) {
        if (is.na(r) || is.na(f)) next
        if (!is_route_compatible(f, r)) {
          rows[[idx]] <- data.frame(
            route_norm = r,
            form_norm = f,
            drugbank_id = df$drugbank_id[[i]],
            raw_route = df$route_raw[[i]],
            raw_form = df$form_raw[[i]],
            raw_route_original = df$raw_route_original[[i]],
            raw_form_original = df$raw_form_original[[i]],
            stringsAsFactors = FALSE
          )
          idx <- idx + 1
        }
      }
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
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

write_csv_and_parquet <- function(df, csv_path) {
  parquet_path <- sub("\\.csv$", ".parquet", csv_path)
  df$write_parquet(parquet_path)
  df$write_csv(csv_path)
  parquet_path
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_master_csv <- file.path(output_dir, "drugbank_generics_master.csv")

dataset <- drugbank

groups_r <- data.frame(
  drugbank_id = as.character(dataset$drugs$groups$drugbank_id),
  group_clean = tolower(trimws(dataset$drugs$groups$group)),
  stringsAsFactors = FALSE
)
groups_df <- to_polars_df(groups_r)
excluded_ids <- unique(groups_r$drugbank_id[groups_r$group_clean == "vet"])

filter_excluded <- function(df, id_col = "drugbank_id") {
  if (!length(excluded_ids)) return(df)
  df$filter(!pl$col(id_col)$is_in(excluded_ids))
}

general_r <- data.frame(
  drugbank_id = as.character(dataset$drugs$general_information$drugbank_id),
  canonical_generic_name = collapse_ws(dataset$drugs$general_information$name),
  stringsAsFactors = FALSE
)
general_r <- general_r[general_r$drugbank_id %nin% excluded_ids, ]
general_r <- general_r[!is.na(general_r$canonical_generic_name) & general_r$canonical_generic_name != "", ]
general_r$generic_components_raw <- lapply(general_r$canonical_generic_name, split_ingredients)
general_r$generic_components <- vapply(general_r$generic_components_raw, function(vec) {
  if (!length(vec)) return(NA_character_)
  paste(vec, collapse = "; ")
}, character(1), USE.NAMES = FALSE)
general_r$generic_components_key <- vapply(general_r$generic_components_raw, function(vec) {
  if (!length(vec)) return(NA_character_)
  canon <- tolower(trimws(gsub("\\s+", " ", vec)))
  canon <- canon[nzchar(canon)]
  if (!length(canon)) return(NA_character_)
  paste(sort(canon), collapse = "||")
}, character(1), USE.NAMES = FALSE)
general_r$generic_components_raw <- NULL
general_df <- to_polars_df(general_r)

allowed_synonym_coders <- c("inn", "usan", "ban", "jan", "dcj", "usp", "dcit")
syn_raw <- data.frame(
  drugbank_id = as.character(dataset$drugs$synonyms$drugbank_id),
  synonym = collapse_ws(dataset$drugs$synonyms$synonym),
  language = tolower(trimws(dataset$drugs$synonyms$language)),
  coder = tolower(trimws(dataset$drugs$synonyms$coder)),
  stringsAsFactors = FALSE
)
syn_raw <- syn_raw[!(is.na(syn_raw$synonym) | syn_raw$synonym == ""), ]
syn_raw <- syn_raw[!(is.na(syn_raw$language) | !grepl("english", syn_raw$language) | is.na(syn_raw$coder)), ]
syn_raw$coder_tokens <- lapply(syn_raw$coder, split_coder_tokens)
syn_raw$has_allowed <- vapply(syn_raw$coder_tokens, function(vals) length(vals) > 0 && any(vals %in% allowed_synonym_coders), logical(1))
syn_raw$only_iupac <- vapply(syn_raw$coder_tokens, function(vals) length(vals) > 0 && all(vals == "iupac"), logical(1))
syn_raw <- syn_raw[syn_raw$has_allowed & !syn_raw$only_iupac, c("drugbank_id", "synonym")]
if (nrow(syn_raw)) {
  syn_list <- tapply(syn_raw$synonym, syn_raw$drugbank_id, function(x) list(unique_canonical(x)))
  syn_df <- to_polars_df(list(
    drugbank_id = names(syn_list),
    synonyms_list = unname(syn_list)
  ))
} else {
  syn_df <- to_polars_df(list(drugbank_id = character(), synonyms_list = list()))
}
syn_df <- filter_excluded(syn_df)

atc_raw <- data.frame(
  drugbank_id = as.character(dataset$drugs$atc_codes$drugbank_id),
  atc_code = trimws(dataset$drugs$atc_codes$atc_code),
  stringsAsFactors = FALSE
)
atc_raw <- atc_raw[!(is.na(atc_raw$atc_code) | atc_raw$atc_code == ""), ]
if (nrow(atc_raw)) {
  atc_list <- tapply(atc_raw$atc_code, atc_raw$drugbank_id, function(x) list(unique_canonical(x)))
  atc_df <- to_polars_df(list(
    drugbank_id = names(atc_list),
    atc_codes_list = unname(atc_list)
  ))
} else {
  atc_df <- to_polars_df(list(drugbank_id = character(), atc_codes_list = list()))
}
atc_df <- filter_excluded(atc_df)

groups_r_nonvet <- groups_r[!(groups_r$drugbank_id %in% excluded_ids) & groups_r$group_clean != "vet", ]
if (nrow(groups_r_nonvet)) {
  groups_list <- tapply(groups_r_nonvet$group_clean, groups_r_nonvet$drugbank_id, function(x) list(unique_canonical(x)))
  groups_clean <- to_polars_df(list(
    drugbank_id = names(groups_list),
    groups_list = unname(groups_list)
  ))
} else {
  groups_clean <- to_polars_df(list(drugbank_id = character(), groups_list = list()))
}

salts_raw <- data.frame(
  drugbank_id = as.character(dataset$salts$drugbank_id),
  salt_name = collapse_ws(dataset$salts$name),
  stringsAsFactors = FALSE
)
salts_raw <- salts_raw[!(is.na(salts_raw$salt_name) | salts_raw$salt_name == ""), ]
if (nrow(salts_raw)) {
  salts_list <- tapply(salts_raw$salt_name, salts_raw$drugbank_id, function(x) list(unique_canonical(x)))
  salts_df <- to_polars_df(list(
    drugbank_id = names(salts_list),
    salt_names_list = unname(salts_list)
  ))
} else {
  salts_df <- to_polars_df(list(drugbank_id = character(), salt_names_list = list()))
}
salts_df <- filter_excluded(salts_df)

process_source <- function(df) {
  if (is_polars_lf(df)) {
    df <- df$collect()
  }
  if (is_polars_df(df)) {
    df <- as.data.frame(df)
  }
  if (is.null(df) || !nrow(df)) {
    return(to_polars_df(list(
      drugbank_id = character(),
      route_raw = character(),
      form_raw = character(),
      dose_raw = character(),
      raw_route_original = character(),
      raw_form_original = character(),
      raw_route_details = character(),
      raw_form_details = character(),
      route_norm_list = list(),
      form_norm_list = list(),
      dose_norm = character(),
      raw_dose = character()
    )))
  }
  df$drugbank_id <- as.character(df$drugbank_id)
  df$route_raw <- collapse_ws(df$route_raw)
  df$form_raw <- collapse_ws(df$form_raw)
  df$dose_raw <- collapse_ws(df$dose_raw)
  df <- df[!(df$drugbank_id %in% excluded_ids), ]
  if (!nrow(df)) {
    return(to_polars_df(list(
      drugbank_id = character(),
      route_raw = character(),
      form_raw = character(),
      dose_raw = character(),
      raw_route_original = character(),
      raw_form_original = character(),
      raw_route_details = character(),
      raw_form_details = character(),
      route_norm_list = list(),
      form_norm_list = list(),
      dose_norm = character(),
      raw_dose = character()
    )))
  }
  df$raw_route_original <- df$route_raw
  df$raw_form_original <- df$form_raw
  df$raw_route_details <- vapply(df$route_raw, clean_form_route_detail, character(1), USE.NAMES = FALSE)
  df$raw_form_details <- vapply(df$form_raw, clean_form_route_detail, character(1), USE.NAMES = FALSE)
  df$route_raw <- vapply(df$route_raw, clean_form_route_base, character(1), USE.NAMES = FALSE)
  df$form_raw <- vapply(df$form_raw, clean_form_route_base, character(1), USE.NAMES = FALSE)
  df$route_norm_list <- lapply(df$route_raw, normalize_route_entry)
  df$form_norm_list <- lapply(df$form_raw, normalize_form_value)
  df$dose_norm <- vapply(df$dose_raw, normalize_dose_value, character(1), USE.NAMES = FALSE)
  df$route_norm_list <- lapply(df$route_norm_list, function(vals) {
    vals <- unique_canonical(vals)
    if (!length(vals)) return(list(NA_character_))
    vals
  })
  df$form_norm_list <- lapply(df$form_norm_list, function(vals) {
    vals <- unique_canonical(vals)
    if (!length(vals)) return(list(NA_character_))
    vals
  })
  df$form_norm <- vapply(df$form_norm_list, function(vals) {
    vals <- vals[!is.na(vals)]
    if (!length(vals)) return(NA_character_)
    vals[[1]]
  }, character(1), USE.NAMES = FALSE)
  df$dose_norm <- ifelse(is.na(df$dose_norm) | df$dose_norm == "", NA_character_, df$dose_norm)
  to_polars_df(df)
}

if (length(dataset$drugs$dosages)) {
  dosages_raw <- to_polars_df(dataset$drugs$dosages)$select(
    drugbank_id = pl$col("drugbank_id")$cast(pl$String),
    route_raw = pl$col("route"),
    form_raw = pl$col("form"),
    dose_raw = pl$col("strength")
  )
} else {
  dosages_raw <- to_polars_df(list(
    drugbank_id = character(),
    route_raw = character(),
    form_raw = character(),
    dose_raw = character()
  ))
}
dosages_proc <- process_source(dosages_raw)

if (length(dataset$products)) {
  products_raw <- to_polars_df(dataset$products)$select(
    drugbank_id = pl$col("drugbank_id")$cast(pl$String),
    route_raw = pl$col("route"),
    form_raw = pl$col("dosage_form"),
    dose_raw = pl$col("strength")
  )
} else {
  products_raw <- to_polars_df(list(
    drugbank_id = character(),
    route_raw = character(),
    form_raw = character(),
    dose_raw = character()
  ))
}
products_proc <- process_source(products_raw)

combined <- pl$concat(list(dosages_proc, products_proc), how = "vertical")
combined_r <- as.data.frame(combined)

disallowed_pairs <- collect_disallowed_pairs(combined_r)

expand_combos <- function(df) {
  if (!nrow(df)) return(data.frame())
  rows <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    routes <- unique_canonical(df$route_norm_list[[i]])
    forms <- unique_canonical(df$form_norm_list[[i]])
    if (!length(routes) || all(is.na(routes))) routes <- NA_character_
    if (!length(forms) || all(is.na(forms))) forms <- NA_character_
    if (length(routes) == 1L && length(forms) == 1L) {
      pairs <- data.frame(route_norm = routes, form_norm = forms, stringsAsFactors = FALSE)
    } else {
      pairs <- expand.grid(
        route_norm = routes,
        form_norm = forms,
        stringsAsFactors = FALSE
      )
    }
    if (!nrow(pairs)) {
      pairs <- data.frame(route_norm = NA_character_, form_norm = NA_character_, stringsAsFactors = FALSE)
    }
    n <- nrow(pairs)
    rows[[i]] <- data.frame(
      drugbank_id = rep(df$drugbank_id[[i]], n),
      route_norm = pairs$route_norm,
      form_norm = pairs$form_norm,
      raw_route = rep(df$route_raw[[i]], n),
      raw_route_original = rep(df$raw_route_original[[i]], n),
      raw_form = rep(df$form_raw[[i]], n),
      raw_form_original = rep(df$raw_form_original[[i]], n),
      raw_route_details = rep(df$raw_route_details[[i]], n),
      raw_form_details = rep(df$raw_form_details[[i]], n),
      dose_norm = rep(df$dose_norm[[i]], n),
      raw_dose = rep(df$dose_raw[[i]], n),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

combo_base <- expand_combos(combined_r)
if (nrow(combo_base)) {
  dedup_cols <- c(
    "drugbank_id",
    "dose_norm",
    "raw_dose",
    "form_norm",
    "raw_form",
    "raw_form_details",
    "route_norm",
    "raw_route",
    "raw_route_details"
  )
  combo_base <- combo_base[!duplicated(combo_base[, dedup_cols, drop = FALSE]), ]
  combo_base <- combo_base[!(is.na(combo_base$dose_norm) & is.na(combo_base$form_norm) & is.na(combo_base$route_norm)), ]
}

if (!nrow(combo_base)) {
  stop("No real combinations found in DrugBank data.")
}

combo_df <- to_polars_df(combo_base)
distinct_ids <- as.vector(combo_df$select("drugbank_id")$unique()$get_column("drugbank_id"))

empty_atc_list <- pl$lit(list())$cast(pl$List(pl$String))
atc_map <- to_polars_df(list(drugbank_id = distinct_ids))$join(atc_df, on = "drugbank_id", how = "left")
atc_map <- atc_map$with_columns(
  atc_codes_list = pl$when(pl$col("atc_codes_list")$is_null())$then(empty_atc_list)$otherwise(pl$col("atc_codes_list"))
)
atc_long <- atc_map$explode("atc_codes_list")$rename(list(atc_codes_list = "atc_code"))

name_df <- general_df$join(syn_df, on = "drugbank_id", how = "left")$with_columns(
  synonyms_list = pl$when(pl$col("synonyms_list")$is_null())$then(pl$lit(list())$cast(pl$List(pl$String)))$otherwise(pl$col("synonyms_list"))
)
name_df <- name_df$with_columns(
  lexeme_list = pl$struct(c("canonical_generic_name", "synonyms_list"))$map_batches(
    function(s) {
      canonical <- s$canonical_generic_name
      syns <- s$synonyms_list
      vals <- unique_canonical(c(canonical, syns))
      if (!length(vals)) return(list(canonical))
      vals
    },
    return_dtype = pl$List(pl$String)
  )
)

combo_dt <- combo_df$join(atc_long, on = "drugbank_id", how = "left")
combo_dt <- combo_dt$join(
  name_df$select("drugbank_id", "canonical_generic_name", "generic_components", "generic_components_key", "lexeme_list"),
  on = "drugbank_id",
  how = "left"
)
combo_dt <- combo_dt$join(salts_df, on = "drugbank_id", how = "left")
combo_dt <- combo_dt$with_columns(
  salt_names_list = pl$when(pl$col("salt_names_list")$is_null())$then(pl$lit(list())$cast(pl$List(pl$String)))$otherwise(pl$col("salt_names_list")),
  salt_names = pl$col("salt_names_list")$map_batches(function(lst) collapse_pipe(expand_salt_set(lst)), return_dtype = pl$String)
)$select(-"salt_names_list")

combo_dt <- combo_dt$join(groups_clean, on = "drugbank_id", how = "left")
combo_dt <- combo_dt$with_columns(
  groups_list = pl$when(pl$col("groups_list")$is_null())$then(pl$lit(list())$cast(pl$List(pl$String)))$otherwise(pl$col("groups_list")),
  groups = pl$col("groups_list")$map_batches(collapse_pipe, return_dtype = pl$String)
)$select(-"groups_list")

combo_dt <- combo_dt$with_columns(
  lexeme_list = pl$when(pl$col("lexeme_list")$is_null())$then(
    pl$col("canonical_generic_name")$map_batches(function(x) list(x), return_dtype = pl$List(pl$String))
  )$otherwise(pl$col("lexeme_list"))
)
combo_dt <- combo_dt$explode("lexeme_list")$rename(list(lexeme_list = "lexeme"))

combo_dt <- combo_dt$with_columns(
  dose_synonyms = pl$struct(c("dose_norm", "raw_dose"))$map_batches(
    function(x) collapse_pipe(expand_dose_set(x$dose_norm, x$raw_dose)),
    return_dtype = pl$String
  ),
  form_synonyms = pl$struct(c("form_norm", "raw_form"))$map_batches(
    function(x) collapse_pipe(expand_value_set(x$form_norm, x$raw_form, FORM_SYNONYM_LOOKUP)),
    return_dtype = pl$String
  ),
  route_synonyms = pl$struct(c("route_norm", "raw_route"))$map_batches(
    function(x) collapse_pipe(expand_value_set(x$route_norm, x$raw_route, ROUTE_SYNONYM_LOOKUP)),
    return_dtype = pl$String
  )
)

final_df <- combo_dt$select(
  "drugbank_id",
  "lexeme",
  "canonical_generic_name",
  "generic_components",
  "generic_components_key",
  "dose_norm",
  "raw_dose",
  "dose_synonyms",
  "form_norm",
  "raw_form",
  "raw_form_details",
  "raw_form_original",
  "form_synonyms",
  "route_norm",
  "raw_route",
  "raw_route_details",
  "raw_route_original",
  "route_synonyms",
  "atc_code",
  "salt_names",
  "groups"
)

if (!keep_all_flag) {
  final_df <- final_df$with_columns(
    approved_flag = pl$col("groups")$str$contains("approved", literal = FALSE),
    atc_flag = pl$col("atc_code")$is_not_null() & pl$col("atc_code") != ""
  )$filter(pl$col("approved_flag") | pl$col("atc_flag"))$select(
    -"approved_flag",
    -"atc_flag"
  )
}

final_df <- final_df$sort(c(
  "canonical_generic_name",
  "lexeme",
  "atc_code",
  "dose_norm",
  "form_norm",
  "route_norm"
))

write_csv_and_parquet(final_df, output_master_csv)

cat(sprintf("Wrote %d rows to %s (Parquet + CSV)\n", final_df$height, output_master_csv))
if (!quiet_mode) {
  cat("Sample rows:\n")
  print(as.data.frame(final_df$head(5)))
}

if (exists("disallowed_pairs") && nrow(disallowed_pairs)) {
  disallowed_pairs <- merge(
    disallowed_pairs,
    as.data.frame(name_df$select("drugbank_id", "canonical_generic_name")),
    by = "drugbank_id",
    all.x = TRUE
  )
  disallowed_cols <- c(
    "drugbank_id",
    "canonical_generic_name",
    "form_norm",
    "route_norm",
    "raw_form",
    "raw_route",
    "raw_form_original",
    "raw_route_original"
  )
  disallowed_path <- file.path(output_dir, "drugbank_route_form_exclusions.csv")
  disallowed_df <- to_polars_df(disallowed_pairs[, disallowed_cols])
  write_csv_and_parquet(disallowed_df, disallowed_path)
}
