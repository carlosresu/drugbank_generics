#!/usr/bin/env Rscript
# drugbank_generics.R — build a generics-focused DrugBank master dataset
# Columns: drugbank_id, lexeme, canonical_generic_name, generic_components,
#          generic_components_key, dose_norm, raw_dose, dose_synonyms,
#          form_norm, raw_form, form_synonyms, route_norm, raw_route,
#          route_synonyms, atc_code, salt_names, groups

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
keep_all_flag <- "--keep-all" %in% argv
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
  if (is.na(cores)) cores <- 13L
  max(13L, cores)
}

detect_os_name <- function() {
  os <- Sys.info()[["sysname"]]
  if (is.null(os)) os <- .Platform$OS.type
  tolower(os)
}

worker_count <- resolve_workers()
data.table::setDTthreads(worker_count)

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
cat(sprintf("[drugbank_generics] data.table threads: %s | parallel workers: %s\n", data.table::getDTthreads(), worker_count))

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

resolve_chunk_size <- function(total_rows) {
  env_val <- suppressWarnings(as.integer(Sys.getenv("ESOA_DRUGBANK_CHUNK", "")))
  if (!is.na(env_val) && env_val > 0) return(env_val)
  # Aim for multiple chunks per worker while avoiding excessive overhead.
  max(1000L, ceiling(total_rows / max(1L, worker_count * 3L)))
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
  parts <- unlist(strsplit(val, "(?i)\\sand\\s|\\swith\\s|\\splus\\s|\\+|/|,(?=\\s)|;", perl = TRUE))
  parts <- collapse_ws(parts)
  parts <- parts[nzchar(parts)]
  unique_canonical(parts)
}

strip_parenthetical_segments <- function(value) {
  original <- collapse_ws(value)
  if (is.na(original) || !nzchar(original)) {
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
  if (is.na(value) || !nzchar(value)) return(list(base = value, detail = NA_character_))
  match <- regexpr(",\\s+.+$", value, perl = TRUE)
  if (match[1] == -1) return(list(base = value, detail = NA_character_))
  detail <- trimws(substr(value, match[1] + 1, nchar(value)))
  base <- trimws(substr(value, 1, match[1] - 1))
  list(base = if (nzchar(base)) base else NA_character_, detail = if (nzchar(detail)) detail else NA_character_)
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
    detail = if (nzchar(detail_out)) detail_out else NA_character_
  )
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
RELEASE_PATTERN <- "(extended(?:-|\u0020)?release|immediate(?:-|\u0020)?release|delayed(?:-|\u0020)?release|sustained(?:-|\u0020)?release|controlled(?:-|\u0020)?release)"

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
  # Split multi-form strings (e.g., "kit; tablet") before normalizing each token.
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
    } else if (part_trim %chin% ALLOWED_ROUTE_SET) {
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
        } else if (word_trim %chin% ALLOWED_ROUTE_SET) {
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
  route_val %chin% expected
}

collect_disallowed_pairs <- function(dt) {
  if (!nrow(dt)) return(data.table())
  chunk_size <- resolve_chunk_size(nrow(dt))
  worker <- function(ix) {
    rows <- lapply(ix, function(i) {
      routes <- unique_canonical(dt$route_norm_list[[i]])
      forms <- unique_canonical(dt$form_norm_list[[i]])
      if (!length(routes) || all(is.na(routes)) || !length(forms) || all(is.na(forms))) return(NULL)
      pairs <- data.table::data.table(
        route_norm = rep(routes, each = length(forms)),
        form_norm = rep(forms, times = length(routes))
      )
      pairs <- pairs[!is.na(route_norm) & !is.na(form_norm)]
      if (!nrow(pairs)) return(NULL)
      allowed_mask <- mapply(is_route_compatible, pairs$form_norm, pairs$route_norm)
      allowed_mask[is.na(allowed_mask)] <- TRUE
      disallowed <- pairs[!as.logical(allowed_mask)]
      if (!nrow(disallowed)) return(NULL)
      disallowed[, drugbank_id := dt$drugbank_id[[i]]]
      disallowed[, raw_route := dt$route_raw[[i]]]
      disallowed[, raw_form := dt$form_raw[[i]]]
      disallowed[, raw_route_original := dt$raw_route_original[[i]]]
      disallowed[, raw_form_original := dt$raw_form_original[[i]]]
      disallowed
    })
    data.table::rbindlist(rows, fill = TRUE)
  }
  idx <- split(seq_len(nrow(dt)), ceiling(seq_len(nrow(dt)) / chunk_size))
  parts <- parallel_lapply(idx, worker)
  data.table::rbindlist(parts, fill = TRUE, use.names = TRUE)
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_master_path <- file.path(output_dir, "drugbank_generics_master.csv")

dataset <- drugbank

groups_dt <- as.data.table(dataset$drugs$groups)
groups_dt[, drugbank_id := as.character(drugbank_id)]
groups_dt[, group_clean := tolower(trimws(group))]
excluded_ids <- unique(groups_dt[group_clean %chin% c("vet"), drugbank_id])

filter_excluded <- function(dt, id_col = "drugbank_id") {
  dt[!(get(id_col) %chin% excluded_ids)]
}

# general info
general_dt <- as.data.table(dataset$drugs$general_information)[
  , .(drugbank_id = as.character(drugbank_id), canonical_generic_name = collapse_ws(name))
]
general_dt <- filter_excluded(general_dt)
general_dt <- general_dt[!is.na(canonical_generic_name) & nzchar(canonical_generic_name)]
general_dt[, generic_components_raw := lapply(canonical_generic_name, split_ingredients)]
general_dt[, generic_components_vec := lapply(generic_components_raw, unique_canonical)]
general_dt[, generic_components := sapply(generic_components_vec, function(vec) {
  if (!length(vec)) return(NA_character_)
  paste(vec, collapse = "; ")
})]
general_dt[, generic_components_key := sapply(generic_components_vec, function(vec) {
  if (!length(vec)) return(NA_character_)
  canon <- tolower(trimws(gsub("\\s+", " ", vec)))
  canon <- canon[nzchar(canon)]
  if (!length(canon)) return(NA_character_)
  paste(sort(canon), collapse = "||")
})]
general_dt[, generic_components_vec := NULL]
general_dt[, generic_components_raw := NULL]

# synonyms
syn_dt <- as.data.table(dataset$drugs$synonyms)[
  , .(
    drugbank_id = as.character(drugbank_id),
    synonym = collapse_ws(synonym),
    language = tolower(trimws(language)),
    coder = tolower(trimws(coder))
  )
]
syn_dt <- filter_excluded(syn_dt)
allowed_synonym_coders <- c("inn", "usan", "ban", "jan", "dcj", "usp", "dcit")
syn_dt <- syn_dt[
  !is.na(synonym) & nzchar(synonym) &
    !is.na(language) & grepl("english", language, fixed = TRUE) &
    !is.na(coder) & coder %chin% allowed_synonym_coders
]
syn_dt <- syn_dt[, .(synonyms_list = list(unique_canonical(synonym))), by = drugbank_id]

# atc
atc_dt <- as.data.table(dataset$drugs$atc_codes)[
  , .(drugbank_id = as.character(drugbank_id), atc_code = trimws(atc_code))
]
atc_dt <- filter_excluded(atc_dt)
atc_dt <- atc_dt[!is.na(atc_code) & nzchar(atc_code)]
atc_dt <- atc_dt[, .(atc_codes_list = list(unique_canonical(atc_code))), by = drugbank_id]

# groups aggregation
groups_clean_dt <- groups_dt[!(drugbank_id %chin% excluded_ids)]
groups_clean_dt <- groups_clean_dt[, .(groups_list = list(unique_canonical(group_clean))), by = drugbank_id]

# salts
salts_dt <- as.data.table(dataset$salts)[
  , .(drugbank_id = as.character(drugbank_id), salt_name = collapse_ws(name))
]
salts_dt <- filter_excluded(salts_dt)
salts_dt <- salts_dt[!is.na(salt_name) & nzchar(salt_name)]
salts_dt <- salts_dt[, .(salt_names_list = list(unique_canonical(salt_name))), by = drugbank_id]

process_source <- function(dt) {
  if (!nrow(dt)) return(dt)
  chunk_size <- resolve_chunk_size(nrow(dt))
  process_block <- function(block) {
    source_dt <- copy(block)
    source_dt <- filter_excluded(source_dt)
    if (!"route_raw" %in% names(source_dt)) source_dt[, route_raw := NA_character_]
    if (!"form_raw" %in% names(source_dt)) source_dt[, form_raw := NA_character_]
    if (!"dose_raw" %in% names(source_dt)) source_dt[, dose_raw := NA_character_]
    source_dt[, route_raw := collapse_ws(route_raw)]
    source_dt[, form_raw := collapse_ws(form_raw)]
    source_dt[, dose_raw := collapse_ws(dose_raw)]
    source_dt[, raw_form_original := form_raw]
    source_dt[, raw_route_original := route_raw]
    form_vec <- source_dt$form_raw
    route_vec <- source_dt$route_raw
    form_clean_list <- lapply(form_vec, clean_form_route_entry)
    route_clean_list <- lapply(route_vec, clean_form_route_entry)
    source_dt[, raw_form_details := vapply(form_clean_list, function(x) x$detail, character(1), USE.NAMES = FALSE)]
    source_dt[, raw_route_details := vapply(route_clean_list, function(x) x$detail, character(1), USE.NAMES = FALSE)]
    source_dt[, raw_form_original_key := vapply(raw_form_original, canonical_case_key, character(1), USE.NAMES = FALSE)]
    source_dt[, raw_form_details_key := vapply(raw_form_details, canonical_case_key, character(1), USE.NAMES = FALSE)]
    source_dt[, raw_route_original_key := vapply(raw_route_original, canonical_case_key, character(1), USE.NAMES = FALSE)]
    source_dt[, raw_route_details_key := vapply(raw_route_details, canonical_case_key, character(1), USE.NAMES = FALSE)]
    source_dt[, form_raw := vapply(form_clean_list, function(x) collapse_ws(x$base), character(1), USE.NAMES = FALSE)]
    # Preserve full route string to allow multi-token splits downstream.
    source_dt[, route_raw := route_raw]
    source_dt[, route_norm_list := parallel_lapply(as.list(route_raw), normalize_route_entry)]
    source_dt[, form_norm_list := parallel_lapply(as.list(form_raw), normalize_form_value)]
    source_dt[, form_norm_list := lapply(form_norm_list, function(vals) {
      vals <- unique_canonical(vals)
      if (!length(vals)) vals <- NA_character_
      vals
    })]
    source_dt[, form_norm := vapply(form_norm_list, function(vals) {
      if (length(vals) && !all(is.na(vals))) vals[[1]] else NA_character_
    }, character(1), USE.NAMES = FALSE)]
    source_dt[, dose_norm := unlist(parallel_lapply(as.list(dose_raw), normalize_dose_value), use.names = FALSE)]
    source_dt[, raw_dose := dose_raw]
    source_dt[form_norm == "", form_norm := NA_character_]
    source_dt[dose_norm == "", dose_norm := NA_character_]
    source_dt
  }

  if (parallel_enabled && worker_count > 1L && nrow(dt) > 1L) {
    idx <- split(seq_len(nrow(dt)), ceiling(seq_len(nrow(dt)) / chunk_size))
    parts <- parallel_lapply(idx, function(ix) process_block(dt[ix]))
    return(rbindlist(parts, fill = TRUE, use.names = TRUE))
  }
  source_dt <- copy(dt)
  source_dt <- filter_excluded(source_dt)
  if (!"route_raw" %in% names(source_dt)) source_dt[, route_raw := NA_character_]
  if (!"form_raw" %in% names(source_dt)) source_dt[, form_raw := NA_character_]
  if (!"dose_raw" %in% names(source_dt)) source_dt[, dose_raw := NA_character_]
  source_dt[, route_raw := collapse_ws(route_raw)]
  source_dt[, form_raw := collapse_ws(form_raw)]
  source_dt[, dose_raw := collapse_ws(dose_raw)]
  source_dt[, raw_form_original := form_raw]
  source_dt[, raw_route_original := route_raw]
  form_vec <- source_dt$form_raw
  route_vec <- source_dt$route_raw
  form_clean_list <- lapply(form_vec, clean_form_route_entry)
  route_clean_list <- lapply(route_vec, clean_form_route_entry)
  source_dt[, raw_form_details := vapply(form_clean_list, function(x) x$detail, character(1), USE.NAMES = FALSE)]
  source_dt[, raw_route_details := vapply(route_clean_list, function(x) x$detail, character(1), USE.NAMES = FALSE)]
  source_dt[, raw_form_original_key := vapply(raw_form_original, canonical_case_key, character(1), USE.NAMES = FALSE)]
  source_dt[, raw_form_details_key := vapply(raw_form_details, canonical_case_key, character(1), USE.NAMES = FALSE)]
  source_dt[, raw_route_original_key := vapply(raw_route_original, canonical_case_key, character(1), USE.NAMES = FALSE)]
  source_dt[, raw_route_details_key := vapply(raw_route_details, canonical_case_key, character(1), USE.NAMES = FALSE)]
  source_dt[, form_raw := vapply(form_clean_list, function(x) x$base, character(1), USE.NAMES = FALSE)]
  source_dt[, route_raw := vapply(route_clean_list, function(x) x$base, character(1), USE.NAMES = FALSE)]
  source_dt[, route_norm_list := parallel_lapply(as.list(route_raw), normalize_route_entry)]
  source_dt[, form_norm_list := parallel_lapply(as.list(form_raw), normalize_form_value)]
  source_dt[, form_norm_list := lapply(form_norm_list, function(vals) {
    vals <- unique_canonical(vals)
    if (!length(vals)) vals <- NA_character_
    vals
  })]
  source_dt[, form_norm := vapply(form_norm_list, function(vals) {
    if (length(vals) && !all(is.na(vals))) vals[[1]] else NA_character_
  }, character(1), USE.NAMES = FALSE)]
  source_dt[, dose_norm := unlist(parallel_lapply(as.list(dose_raw), normalize_dose_value), use.names = FALSE)]
  source_dt[, raw_dose := dose_raw]
  source_dt[form_norm == "", form_norm := NA_character_]
  source_dt[dose_norm == "", dose_norm := NA_character_]
  source_dt
}

# dosages
if (length(dataset$drugs$dosages)) {
  dosages_raw <- as.data.table(dataset$drugs$dosages)[
    , .(drugbank_id = as.character(drugbank_id), route = route, form = form, strength = strength)
  ]
  setnames(dosages_raw, c("route", "form", "strength"), c("route_raw", "form_raw", "dose_raw"))
} else {
  dosages_raw <- data.table(drugbank_id = character(), route_raw = character(), form_raw = character(), dose_raw = character())
}
dosages_proc <- process_source(dosages_raw)

# products
if (length(dataset$products)) {
  products_raw <- as.data.table(dataset$products)[
    , .(drugbank_id = as.character(drugbank_id), route = route, form = dosage_form, strength = strength)
  ]
  setnames(products_raw, c("route", "form", "strength"), c("route_raw", "form_raw", "dose_raw"))
} else {
  products_raw <- data.table(drugbank_id = character(), route_raw = character(), form_raw = character(), dose_raw = character())
}
products_proc <- process_source(products_raw)

combined <- rbindlist(list(dosages_proc, products_proc), fill = TRUE)
disallowed_pairs <- collect_disallowed_pairs(combined)
combined[, route_norm_list := lapply(route_norm_list, function(x) {
  vals <- unique_canonical(x)
  if (!length(vals)) vals <- NA_character_
  vals
})]
combined[, form_norm_list := lapply(form_norm_list, function(x) {
  vals <- unique_canonical(x)
  if (!length(vals)) vals <- NA_character_
  vals
})]

expand_combos_chunk <- function(chunk_dt) {
  result <- vector("list", nrow(chunk_dt))
  for (i in seq_len(nrow(chunk_dt))) {
    routes <- unique_canonical(chunk_dt$route_norm_list[[i]])
    if (!length(routes) || all(is.na(routes))) routes <- NA_character_
    forms <- unique_canonical(chunk_dt$form_norm_list[[i]])
    if (!length(forms) || all(is.na(forms))) forms <- NA_character_
    if (length(routes) == 1L && length(forms) == 1L) {
      pairs <- data.table::data.table(route_norm = routes, form_norm = forms)
    } else {
      pairs <- data.table::data.table(
        route_norm = rep(routes, each = length(forms)),
        form_norm = rep(forms, times = length(routes))
      )
      if (!nrow(pairs)) pairs <- data.table::data.table(route_norm = NA_character_, form_norm = NA_character_)
      data.table::setkeyv(pairs, c("route_norm", "form_norm"))
      pairs <- unique(pairs)
    }
    n <- nrow(pairs)
    raw_route_val <- if (length(chunk_dt$route_raw)) chunk_dt$route_raw[[i]] else NA_character_
    raw_route_original_val <- if (length(chunk_dt$raw_route_original)) chunk_dt$raw_route_original[[i]] else NA_character_
    raw_route_details_val <- if (length(chunk_dt$raw_route_details)) chunk_dt$raw_route_details[[i]] else NA_character_
    raw_route_original_key_val <- if (length(chunk_dt$raw_route_original_key)) chunk_dt$raw_route_original_key[[i]] else NA_character_
    raw_route_details_key_val <- if (length(chunk_dt$raw_route_details_key)) chunk_dt$raw_route_details_key[[i]] else NA_character_
    raw_form_val <- if (length(chunk_dt$form_raw)) chunk_dt$form_raw[[i]] else NA_character_
    raw_form_original_val <- if (length(chunk_dt$raw_form_original)) chunk_dt$raw_form_original[[i]] else NA_character_
    raw_form_details_val <- if (length(chunk_dt$raw_form_details)) chunk_dt$raw_form_details[[i]] else NA_character_
    raw_form_original_key_val <- if (length(chunk_dt$raw_form_original_key)) chunk_dt$raw_form_original_key[[i]] else NA_character_
    raw_form_details_key_val <- if (length(chunk_dt$raw_form_details_key)) chunk_dt$raw_form_details_key[[i]] else NA_character_
    dose_norm_val <- if (length(chunk_dt$dose_norm)) chunk_dt$dose_norm[[i]] else NA_character_
    raw_dose_val <- if (length(chunk_dt$raw_dose)) chunk_dt$raw_dose[[i]] else NA_character_
    result[[i]] <- data.table(
      drugbank_id = chunk_dt$drugbank_id[[i]],
      route_norm = pairs$route_norm,
      form_norm = pairs$form_norm,
      raw_route = rep(raw_route_val, n),
      raw_route_original = rep(raw_route_original_val, n),
      raw_route_details = rep(raw_route_details_val, n),
      raw_route_original_key = rep(raw_route_original_key_val, n),
      raw_route_details_key = rep(raw_route_details_key_val, n),
      raw_form = rep(raw_form_val, n),
      raw_form_original = rep(raw_form_original_val, n),
      raw_form_details = rep(raw_form_details_val, n),
      raw_form_original_key = rep(raw_form_original_key_val, n),
      raw_form_details_key = rep(raw_form_details_key_val, n),
      dose_norm = rep(dose_norm_val, n),
      raw_dose = rep(raw_dose_val, n)
    )
  }
  data.table::rbindlist(result, fill = TRUE, use.names = TRUE)
}

combo_chunk_size <- resolve_chunk_size(nrow(combined))
combo_indices <- split(seq_len(nrow(combined)), ceiling(seq_len(nrow(combined)) / combo_chunk_size))
combo_parts <- parallel_lapply(combo_indices, function(ix) expand_combos_chunk(combined[ix]))
combo_base <- data.table::rbindlist(combo_parts, fill = TRUE, use.names = TRUE)
dedup_cols <- c(
  "drugbank_id",
  "dose_norm",
  "raw_dose",
  "form_norm",
  "raw_form",
  "raw_form_original_key",
  "raw_form_details_key",
  "route_norm",
  "raw_route",
  "raw_route_original_key",
  "raw_route_details_key"
)
combo_base <- unique(combo_base, by = dedup_cols)
combo_base[, c("raw_form_original_key", "raw_form_details_key", "raw_route_original_key", "raw_route_details_key") := NULL]
combo_base <- combo_base[!(is.na(dose_norm) & is.na(form_norm) & is.na(route_norm))]

if (!nrow(combo_base)) {
  stop("No real combinations found in DrugBank data.")
}

distinct_ids <- unique(combo_base$drugbank_id)

# ensure atc map covers all combos
atc_map <- data.table(drugbank_id = distinct_ids)
atc_map <- merge(atc_map, atc_dt, by = "drugbank_id", all.x = TRUE)
atc_map[is.na(atc_codes_list), atc_codes_list := list(list(character()))]
atc_long <- atc_map[, {
  codes <- atc_codes_list[[1]]
  if (!length(codes)) codes <- NA_character_
  .(atc_code = codes)
}, by = drugbank_id]

combo_dt <- merge(combo_base, atc_long, by = "drugbank_id", allow.cartesian = TRUE)

name_dt <- merge(general_dt[, .(drugbank_id, canonical_generic_name, generic_components, generic_components_key)],
                 syn_dt, by = "drugbank_id", all.x = TRUE)
name_dt[is.na(synonyms_list), synonyms_list := list(list(character()))]
name_dt[, lexeme_list := lapply(seq_len(.N), function(i) unique_canonical(c(canonical_generic_name[i], synonyms_list[[i]])))]

combo_dt <- merge(combo_dt, name_dt[, .(drugbank_id, canonical_generic_name, generic_components, generic_components_key, lexeme_list)],
                  by = "drugbank_id", all.x = TRUE)
combo_dt[is.na(lexeme_list), lexeme_list := list(list(canonical_generic_name))]

combo_dt <- merge(combo_dt, salts_dt, by = "drugbank_id", all.x = TRUE)
combo_dt[is.na(salt_names_list), salt_names_list := list(list(character()))]
combo_dt[, salt_names := unlist(parallel_lapply(salt_names_list, function(lst) collapse_pipe(expand_salt_set(lst))), use.names = FALSE)]
combo_dt[, salt_names_list := NULL]
combo_dt <- merge(combo_dt, groups_clean_dt, by = "drugbank_id", all.x = TRUE)
combo_dt[is.na(groups_list), groups_list := list(list(character()))]
combo_dt[, groups := unlist(parallel_lapply(groups_list, collapse_pipe), use.names = FALSE)]
combo_dt[, groups_list := NULL]

combo_dt <- combo_dt[, {
  lexes <- lexeme_list[[1]]
  if (!length(lexes)) lexes <- canonical_generic_name
  data.table(lexeme = lexes)
}, by = .(drugbank_id, canonical_generic_name, generic_components, generic_components_key,
          dose_norm, raw_dose, form_norm, raw_form, raw_form_details, raw_form_original,
          route_norm, raw_route, raw_route_details, raw_route_original, atc_code,
          salt_names, groups)]

dose_syn_vec <- unlist(parallel_lapply(seq_len(nrow(combo_dt)), function(i) {
  collapse_pipe(expand_dose_set(combo_dt$dose_norm[i], combo_dt$raw_dose[i]))
}), use.names = FALSE)
form_syn_vec <- unlist(parallel_lapply(seq_len(nrow(combo_dt)), function(i) {
  collapse_pipe(expand_value_set(combo_dt$form_norm[i], combo_dt$raw_form[i], FORM_SYNONYM_LOOKUP))
}), use.names = FALSE)
route_syn_vec <- unlist(parallel_lapply(seq_len(nrow(combo_dt)), function(i) {
  collapse_pipe(expand_value_set(combo_dt$route_norm[i], combo_dt$raw_route[i], ROUTE_SYNONYM_LOOKUP))
}), use.names = FALSE)

combo_dt[, dose_synonyms := dose_syn_vec]
combo_dt[, form_synonyms := form_syn_vec]
combo_dt[, route_synonyms := route_syn_vec]

final_dt <- combo_dt[, .(
  drugbank_id,
  lexeme,
  canonical_generic_name,
  generic_components,
  generic_components_key,
  dose_norm,
  raw_dose,
  dose_synonyms,
  form_norm,
  raw_form,
  raw_form_details,
  raw_form_original,
  form_synonyms,
  route_norm,
  raw_route,
  raw_route_details,
  raw_route_original,
  route_synonyms,
  atc_code,
  salt_names,
  groups
)]

write_csv_and_parquet <- function(dt, csv_path) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("arrow package is required for Parquet export.")
  }
  arrow_table <- arrow::Table$create(dt)
  if ("write_csv_arrow" %chin% getNamespaceExports("arrow")) {
    arrow::write_csv_arrow(arrow_table, csv_path)
  } else if ("write_delim_arrow" %chin% getNamespaceExports("arrow")) {
    arrow::write_delim_arrow(arrow_table, csv_path, delim = ",")
  } else {
    data.table::fwrite(dt, csv_path)
  }
  parquet_path <- sub("\\.csv$", ".parquet", csv_path)
  arrow::write_parquet(arrow_table, parquet_path)
  parquet_path
}

if (exists("disallowed_pairs") && nrow(disallowed_pairs)) {
  disallowed_pairs <- merge(
    disallowed_pairs,
    name_dt[, .(drugbank_id, canonical_generic_name)],
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
  write_csv_and_parquet(disallowed_pairs[, ..disallowed_cols], disallowed_path)
}

if (!keep_all_flag) {
  final_dt[, approved_flag := grepl("approved", groups, ignore.case = TRUE)]
  final_dt[, atc_flag := !is.na(atc_code) & nzchar(trimws(atc_code))]
  final_dt <- final_dt[approved_flag | atc_flag]
  final_dt[, c("approved_flag", "atc_flag") := NULL]
}

setorder(final_dt, canonical_generic_name, lexeme, atc_code, dose_norm, form_norm, route_norm)

write_csv_and_parquet(final_dt, output_master_path)

cat(sprintf("Wrote %d rows to %s\n", nrow(final_dt), output_master_path))
if (!quiet_mode) {
  cat("Sample rows:\n")
  print(head(final_dt, 5))
}

if (!is.null(plan_reset)) {
  try(future::plan(future::sequential), silent = TRUE)
}
