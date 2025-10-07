# =================================
# DrugBank → authoritative SOA text dataset (single CSV)
# + Authoritative text→ATC vocabulary for SOA matching
# Memory-aware (data.table + progressr)
#
# Rules (generic names):
# - Keep ALL names from general_information$name exactly as-is (preserve casing/Unicode).
# - Add synonyms only where:
#     coder contains any of {inn, ban, usan, jan, usp} (case-insensitive),
#     language contains "english",
#     and the synonym (case/whitespace-insensitively) is not already in the GI set for that drug.
# - Do NOT infer combinations from text/formatting. Do NOT turn names into ';' lists.
# - Add new column `is_combo` using ONLY ATC structure (e.g., suffix ≥ 50 or special combo groups).
# - Do NOT add salts at all.
# Output:
#   output/drugbank_authority.csv
#     (columns: text, text_type, text_subtype, generic_name, brand_name,
#      brand_company, drugbank_id, atc_code, atc_level_1_code, atc_level_1_desc,
#      atc_level_2_code, atc_level_2_desc, atc_level_3_code, atc_level_3_desc,
#      atc_level_4_code, atc_level_4_desc, is_combo, route, dosage_form,
#      strength, country, source, labeller, data_origin)
# =================================

ensure_installed <- function(packages, repos = "https://cloud.r-project.org") {
  for (pkg in packages) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = repos)
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
ensure_installed(c("data.table", "stringi", "progressr", "devtools"))
if (!requireNamespace("dbdataset", quietly = TRUE)) {
  devtools::install_github("interstellar-Consultation-Services/dbdataset", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(stringi)
  library(progressr)
  library(dbdataset)
})

data.table::setDTthreads(percent = 100)

handlers(global = TRUE)
handlers("cli") # or handlers("txtprogressbar")

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, "drugbank_authority.csv")

# ---------- helpers ----------
collapse_ws <- function(x) {
  x <- gsub("\\s+", " ", x, perl = TRUE)
  trimws(x)
}
dedupe_key <- function(x) {
  x <- collapse_ws(x)
  tolower(x)
}

simplify_form_one <- function(f) {
  if (is.null(f) || is.na(f) || !nzchar(f)) {
    return(NA_character_)
  }
  s <- tolower(gsub("\\s+", " ", trimws(f)))
  parts <- strsplit(s, ",", fixed = TRUE)[[1]]
  main <- trimws(parts[1])
  keep_suffix <- NA_character_
  if (length(parts) > 1) {
    suffixes <- trimws(parts[-1])
    i <- which(grepl("release", suffixes, fixed = TRUE))[1]
    if (!is.na(i)) keep_suffix <- suffixes[i]
  }
  if (is.na(keep_suffix)) main else paste0(main, ", ", keep_suffix)
}
simplify_form_vec <- function(v) vapply(v, simplify_form_one, character(1))

normalize_route_field <- function(v) {
  v <- trimws(v)
  v[is.na(v)] <- NA_character_
  tolower(v)
}
split_semicolon <- function(s) {
  if (is.na(s) || !nzchar(s)) {
    return(character(0))
  }
  parts <- strsplit(s, ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  tolower(parts)
}

split_ingredients <- function(x) {
  if (length(x) == 0L) {
    return(character(0))
  }
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) {
    return(character(0))
  }
  unlist(lapply(x, function(s) {
    y <- gsub("\\s+(and|with|plus)\\s+", "|", s, ignore.case = TRUE, perl = TRUE)
    y <- gsub("\\s*/\\s*|\\s*\\+\\s*|\\s*,\\s*|\\s*;\\s*", "|", y, perl = TRUE)
    parts <- strsplit(y, "\\|", perl = TRUE)[[1]]
    parts <- collapse_ws(parts)
    parts <- parts[nzchar(parts)]
    unique(parts)
  }), use.names = FALSE)
}

extract_product_brand <- function(names, brand_pool) {
  if (!length(names)) {
    return(character(0))
  }
  res <- rep(NA_character_, length(names))
  if (is.null(brand_pool) || length(brand_pool) == 0L) {
    return(res)
  }
  keep <- brand_pool[!is.na(brand_pool) & nzchar(brand_pool)]
  keep <- unique(keep)
  if (!length(keep)) {
    return(res)
  }
  ord <- order(-nchar(keep), tolower(keep))
  keep <- keep[ord]

  names_chr <- as.character(names)
  names_chr[is.na(names_chr)] <- ""
  names_lc <- stringi::stri_trans_tolower(names_chr)

  remaining <- seq_along(names_chr)
  for (brand in keep) {
    brand_lc <- stringi::stri_trans_tolower(brand)
    pattern <- paste0("^\\s*", stringi::stri_replace_all_regex(brand_lc, "([\\^$.*+?()\\[\\]{}|\\\\])", "\\\\$1"), "(\\b|[^[:alnum:]])")
    hits <- stringi::stri_detect_regex(names_lc[remaining], pattern)
    if (any(hits)) {
      res[remaining[hits]] <- brand
      remaining <- remaining[!hits]
    }
    if (!length(remaining)) break
  }
  if (length(remaining)) {
    for (brand in keep) {
      brand_lc <- stringi::stri_trans_tolower(brand)
      hits <- stringi::stri_detect_fixed(names_lc[remaining], brand_lc)
      if (any(hits)) {
        res[remaining[hits]] <- brand
        remaining <- remaining[!hits]
      }
      if (!length(remaining)) break
    }
  }
  if (length(remaining)) {
    fallback <- stringi::stri_trim_both(stringi::stri_extract_first_regex(names_chr[remaining], "^[[:alpha:]][^\\s,/;:]*"))
    fallback[is.na(fallback) | !nzchar(fallback)] <- NA_character_
    res[remaining] <- fallback
  }
  res
}

simplify_form_column <- function(DT, column, label) {
  n <- nrow(DT)
  if (!n) return(invisible(NULL))
  chunk_size <- 50000L
  starts <- seq.int(1L, n, by = chunk_size)
  progressr::with_progress({
    p <- progressr::progressor(steps = length(starts), label = label)
    for (start in starts) {
      end <- min(start + chunk_size - 1L, n)
      idx <- seq.int(start, end)
      data.table::set(DT, i = idx, j = column, value = simplify_form_vec(DT[[column]][idx]))
      p(message = sprintf("%s rows %d–%d / %d", label, start, end, n))
      if ((end %% (chunk_size * 4L)) == 0L) invisible(gc(FALSE))
    }
  })
  invisible(NULL)
}

expand_route_form <- function(DT) {
  if (!nrow(DT)) {
    return(DT)
  }
  tmp <- copy(DT)
  tmp[, route_list := lapply(route, split_semicolon)]
  tmp[, dosage_form_list := lapply(dosage_form, split_semicolon)]
  tmp[, rid := .I]
  expanded <- tmp[, {
    rl <- route_list[[1L]]
    if (length(rl) == 0L) rl <- NA_character_
    fl <- dosage_form_list[[1L]]
    if (length(fl) == 0L) fl <- NA_character_
    CJ(route = rl, dosage_form = fl, unique = TRUE)
  }, by = .(rid)]
  expanded <- merge(
    expanded,
    tmp[, setdiff(names(tmp), c("route", "dosage_form", "route_list", "dosage_form_list")), with = FALSE],
    by = "rid",
    all.x = TRUE,
    sort = FALSE,
    allow.cartesian = TRUE
  )
  expanded[, rid := NULL]
  expanded[route == "", route := NA_character_]
  expanded[dosage_form == "", dosage_form := NA_character_]
  expanded[, route := normalize_route_field(route)]
  expanded[, dosage_form := simplify_form_vec(dosage_form)]
  unique(expanded)
}

# ---------- load ----------
cat("Loading DrugBank slices…\n")
flush.console()
dataset <- drugbank

excluded_groups <- c("experimental", "withdrawn", "illicit", "vet")
grp <- as.data.table(dataset$drugs$groups)
grp[, drugbank_id := as.character(drugbank_id)]
grp[, group_lc := tolower(trimws(group))]
bad_ids <- unique(grp[group_lc %chin% excluded_groups, drugbank_id])
rm(grp)
invisible(gc())

# ---------- ATC (codes + hierarchy + combo flag) ----------
atc_tbl <- as.data.table(dataset$drugs$atc_codes)[
  !is.na(atc_code) & nzchar(trimws(atc_code)),
  .(
    drugbank_id = as.character(drugbank_id),
    atc_code = trimws(atc_code),
    atc_level_1_code = trimws(code_1),
    atc_level_1_desc = trimws(level_1),
    atc_level_2_code = trimws(code_2),
    atc_level_2_desc = trimws(level_2),
    atc_level_3_code = trimws(code_3),
    atc_level_3_desc = trimws(level_3),
    atc_level_4_code = trimws(code_4),
    atc_level_4_desc = trimws(level_4)
  )
]
atc_tbl[, suffix := suppressWarnings(as.integer(sub(".*(..)$", "\\1", atc_code)))]
atc_tbl[, atc_combo := grepl("^J05AR", atc_code) | (!is.na(suffix) & suffix >= 50L)]
combo_flag <- atc_tbl[, .(is_combo = any(atc_combo, na.rm = TRUE)), by = drugbank_id]
atc_tbl[, c("suffix", "atc_combo") := NULL]
atc_tbl <- atc_tbl[!(drugbank_id %chin% bad_ids)]
combo_flag <- combo_flag[!(drugbank_id %chin% bad_ids)]
atc_tbl <- unique(atc_tbl)
setkey(atc_tbl, drugbank_id)
setkey(combo_flag, drugbank_id)

# ---------- Name sources ----------
primary_gi <- as.data.table(dataset$drugs$general_information)[
  , .(drugbank_id = as.character(drugbank_id), generic_name = collapse_ws(name))
]
primary_gi <- primary_gi[!(drugbank_id %chin% bad_ids)]
primary_gi <- primary_gi[!is.na(generic_name) & nzchar(generic_name)]
setkey(primary_gi, drugbank_id)

gi_raw <- as.data.table(dataset$drugs$general_information)[
  , .(drugbank_id = as.character(drugbank_id), text = collapse_ws(name))
]
gi_raw <- gi_raw[!(drugbank_id %chin% bad_ids)]
gi_raw <- gi_raw[!is.na(text) & nzchar(text)]
gi_raw[, text_key := dedupe_key(text)]
generic_keys <- unique(gi_raw[, .(drugbank_id, text_key)])
generic_texts <- unique(gi_raw, by = c("drugbank_id", "text_key"))
generic_texts[, `:=`(text_type = "generic", text_subtype = NA_character_, brand_name = NA_character_, brand_company = NA_character_)]
generic_texts[, text_key := NULL]

keep_coder_re <- "(?i)\\b(inn|ban|usan|jan|usp)\\b"
synonym_texts <- as.data.table(dataset$drugs$synonyms)[
  !is.na(coder) & nzchar(trimws(coder)) &
    grepl(keep_coder_re, coder, perl = TRUE) &
    !is.na(language) & grepl("english", language, ignore.case = TRUE) &
    !is.na(synonym) & nzchar(trimws(synonym)),
  .(
    drugbank_id = as.character(drugbank_id),
    text = collapse_ws(synonym),
    text_subtype = tolower(trimws(coder))
  )
]
synonym_texts <- synonym_texts[!(drugbank_id %chin% bad_ids)]
synonym_texts[, text_key := dedupe_key(text)]
synonym_texts <- synonym_texts[!generic_keys, on = .(drugbank_id, text_key)]
synonym_texts <- unique(synonym_texts, by = c("drugbank_id", "text_key"))
synonym_texts[, `:=`(text_type = "synonym", brand_name = NA_character_, brand_company = NA_character_)]
synonym_texts[, text_key := NULL]

brand_international <- as.data.table(dataset$drugs$international_brands)[
  , .(
    drugbank_id = as.character(drugbank_id),
    text = collapse_ws(brand),
    brand_company = collapse_ws(company)
  )
]
brand_international <- brand_international[!(drugbank_id %chin% bad_ids)]
brand_international <- brand_international[!is.na(text) & nzchar(text)]
brand_international[!nzchar(brand_company), brand_company := NA_character_]
brand_international[, text_key := dedupe_key(text)]
brand_international <- unique(brand_international, by = c("drugbank_id", "text_key"))
brand_international[, `:=`(text_type = "brand", text_subtype = "international", brand_name = text)]
brand_international[, text_key := NULL]

brand_lookup <- brand_international[, .(brand_pool = list(unique(brand_name))), by = drugbank_id]

mix <- as.data.table(dataset$drugs$mixtures)[
  , .(
    drugbank_id = as.character(drugbank_id),
    name = collapse_ws(name),
    ingredients = ingredients
  )
]
mix <- mix[!(drugbank_id %chin% bad_ids)]

mix_name <- mix[
  !is.na(name) & nzchar(name),
  .(
    drugbank_id, text = name,
    text_type = "mixture",
    text_subtype = "mixture_name",
    brand_name = NA_character_,
    brand_company = NA_character_
  )
]

mix_ing <- mix[!is.na(ingredients) & nzchar(ingredients)]
if (nrow(mix_ing)) {
  mix_ing <- mix_ing[, .(text = split_ingredients(ingredients)), by = .(drugbank_id)]
  mix_ing <- mix_ing[!is.na(text) & nzchar(text)]
  mix_ing[, `:=`(
    text_type = "mixture",
    text_subtype = "mixture_ingredient",
    brand_name = NA_character_,
    brand_company = NA_character_
  )]
} else {
  mix_ing <- mix[0, .(
    drugbank_id = character(),
    text = character(),
    text_type = character(),
    text_subtype = character(),
    brand_name = character(),
    brand_company = character()
  )]
}

# ---------- Route/Form/Dose (dosages + products) ----------
dosages_tbl <- as.data.table(dataset$drugs$dosages)[
  , .(
    drugbank_id = as.character(drugbank_id),
    route = normalize_route_field(route),
    dosage_form = trimws(form),
    strength = tolower(trimws(strength)),
    country = NA_character_,
    source = "DOSAGES",
    labeller = NA_character_,
    data_origin = "dosage"
  )
]
dosages_tbl <- dosages_tbl[!(drugbank_id %chin% bad_ids)]
dosages_tbl[!nzchar(dosage_form), dosage_form := NA_character_]
dosages_tbl[!nzchar(strength), strength := NA_character_]

products_base <- as.data.table(dataset$products)[
  , .(
    drugbank_id = as.character(drugbank_id),
    product_name = collapse_ws(name),
    labeller = collapse_ws(labeller),
    route = normalize_route_field(route),
    dosage_form = trimws(dosage_form),
    strength = tolower(trimws(strength)),
    country = trimws(country),
    source = trimws(source)
  )
]
products_base <- products_base[!(drugbank_id %chin% bad_ids)]
for (col in c("product_name", "labeller", "dosage_form", "strength", "country", "source")) {
  products_base[, (col) := {
    val <- get(col)
    val[!nzchar(val)] <- NA_character_
    val
  }]
}
products_base <- merge(products_base, brand_lookup, by = "drugbank_id", all.x = TRUE)
products_base[, brand := extract_product_brand(product_name, brand_pool[[1L]]), by = .(drugbank_id)]
products_base[, brand_pool := NULL]
products_base[, data_origin := "product"]
simplify_form_column(dosages_tbl, "dosage_form", "Simplifying dosage forms (dosages)")
simplify_form_column(products_base, "dosage_form", "Simplifying dosage forms (products)")

route_form_all <- unique(rbindlist(
  list(
    dosages_tbl[, .(drugbank_id, route, dosage_form, strength, country, source, labeller, data_origin)],
    products_base[, .(drugbank_id, route, dosage_form, strength, country, source, labeller, data_origin)]
  ),
  use.names = TRUE,
  fill = TRUE
))

product_brand_admin <- products_base[
  !is.na(brand) & nzchar(brand),
  .(
    drugbank_id, brand_name = brand, route, dosage_form,
    strength, country, source, labeller, data_origin
  )
]

products_base[, product_name := NULL]

# ---------- Assemble non-product texts ----------
non_product_texts <- rbindlist(
  list(
    generic_texts,
    synonym_texts,
    brand_international,
    mix_name,
    mix_ing
  ),
  use.names = TRUE,
  fill = TRUE
)
non_product_texts <- merge(non_product_texts, primary_gi, by = "drugbank_id", all.x = TRUE)
desired_cols <- c(
  "text", "text_type", "text_subtype", "generic_name", "brand_name",
  "brand_company", "drugbank_id", "atc_code",
  "atc_level_1_code", "atc_level_1_desc",
  "atc_level_2_code", "atc_level_2_desc",
  "atc_level_3_code", "atc_level_3_desc",
  "atc_level_4_code", "atc_level_4_desc",
  "is_combo", "route", "dosage_form", "strength",
  "country", "source", "labeller", "data_origin"
)

# ---------- Product-derived brand rows ----------
brand_company_lookup <- unique(brand_international[, .(drugbank_id, brand_name, brand_company)])
product_brand_rows <- merge(
  product_brand_admin,
  brand_company_lookup,
  by = c("drugbank_id", "brand_name"),
  all.x = TRUE
)
product_brand_rows <- merge(
  product_brand_rows,
  primary_gi,
  by = "drugbank_id",
  all.x = TRUE
)
product_brand_rows[, `:=`(
  text = brand_name,
  text_type = "brand",
  text_subtype = "product"
)]

# ---------- Streaming combine + attach ATC ----------
atc_ids <- unique(atc_tbl$drugbank_id)
non_product_texts <- non_product_texts[drugbank_id %chin% atc_ids]
product_brand_rows <- product_brand_rows[drugbank_id %chin% atc_ids]
route_form_all <- route_form_all[drugbank_id %chin% atc_ids]
combo_flag <- combo_flag[drugbank_id %chin% atc_ids]

make_empty_output <- function(cols) {
  template <- setNames(
    lapply(cols, function(col) if (col == "is_combo") logical(0) else character(0)),
    cols
  )
  as.data.table(template)
}

rm(primary_gi)
invisible(gc(FALSE))

setkey(non_product_texts, drugbank_id)
setkey(product_brand_rows, drugbank_id)
setkey(route_form_all, drugbank_id)
setkey(atc_tbl, drugbank_id)
setkey(combo_flag, drugbank_id)

if (file.exists(output_path)) file.remove(output_path)

batch_size <- 200L
all_ids <- atc_ids
if (!length(all_ids)) {
  empty_out <- make_empty_output(desired_cols)
  fwrite(empty_out, output_path)
  message("Wrote: ", output_path, " (0 rows, 0 unique texts)")
} else {
  batches <- split(all_ids, ceiling(seq_along(all_ids) / batch_size))
  rows_written <- 0L
  first_batch <- TRUE
  unique_tracker <- local({
    seen <- new.env(parent = emptyenv(), hash = TRUE)
    count <- 0L
    add <- function(values) {
      vals <- values[!is.na(values)]
      if (!length(vals)) return(invisible(NULL))
      vals <- unique(vals)
      new_mask <- !vapply(vals, exists, logical(1), envir = seen, inherits = FALSE)
      if (!any(new_mask)) return(invisible(NULL))
      new_vals <- vals[new_mask]
      for (v in new_vals) assign(v, TRUE, envir = seen)
      count <<- count + length(new_vals)
      invisible(NULL)
    }
    get <- function() count
    list(add = add, get = get)
  })

  progressr::with_progress({
    p <- progressr::progressor(steps = length(batches), label = "Writing output batches")
    for (i in seq_along(batches)) {
      ids <- batches[[i]]
      non_local <- non_product_texts[J(ids), nomatch = 0L]
      route_local <- route_form_all[J(ids), nomatch = 0L]
      non_rows <- if (nrow(non_local)) {
        merge(
          non_local,
          route_local,
          by = "drugbank_id",
          all.x = TRUE,
          allow.cartesian = TRUE
        )
      } else {
        non_local
      }
      product_local <- product_brand_rows[J(ids), nomatch = 0L]

      chunk_components <- list()
      if (nrow(non_rows)) chunk_components[[length(chunk_components) + 1L]] <- non_rows
      if (nrow(product_local)) chunk_components[[length(chunk_components) + 1L]] <- product_local
      if (!length(chunk_components)) {
        p(message = sprintf("Batch %d/%d (0 input rows)", i, length(batches)))
        next
      }

      combined_local <- rbindlist(chunk_components, use.names = TRUE, fill = TRUE)
      combined_local <- combined_local[!is.na(text) & nzchar(text)]
      if (!nrow(combined_local)) {
        p(message = sprintf("Batch %d/%d (0 rows after text filter)", i, length(batches)))
        next
      }

      atc_local <- atc_tbl[J(ids), nomatch = 0L]
      combined_local <- merge(
        combined_local,
        atc_local,
        by = "drugbank_id",
        all.x = FALSE,
        allow.cartesian = TRUE
      )
      if (!nrow(combined_local)) {
        p(message = sprintf("Batch %d/%d (no ATC rows)", i, length(batches)))
        next
      }

      combined_local <- merge(
        combined_local,
        combo_flag,
        by = "drugbank_id",
        all.x = TRUE
      )
      combined_local[is.na(is_combo), is_combo := FALSE]

      expanded_local <- expand_route_form(combined_local)
      if (!nrow(expanded_local)) {
        p(message = sprintf("Batch %d/%d (0 rows after expand)", i, length(batches)))
        next
      }

      expanded_local[!nzchar(strength), strength := NA_character_]

      missing_cols <- setdiff(desired_cols, names(expanded_local))
      if (length(missing_cols)) {
        for (col in missing_cols) expanded_local[, (col) := NA_character_]
      }
      setcolorder(expanded_local, desired_cols)
      expanded_local <- unique(expanded_local)
      setorder(expanded_local, text, text_type, generic_name, atc_code, route, dosage_form, strength)

      fwrite(
        expanded_local,
        output_path,
        append = !first_batch,
        col.names = first_batch
      )
      rows_written <- rows_written + nrow(expanded_local)
      unique_tracker$add(as.character(expanded_local$text))
      first_batch <- FALSE
      p(message = sprintf("Batch %d/%d (%d rows)", i, length(batches), nrow(expanded_local)))
      if ((i %% 5L) == 0L) invisible(gc(FALSE))
    }
  })

  if (rows_written == 0L) {
    fwrite(make_empty_output(desired_cols), output_path)
  }
  message(
    "Wrote: ", output_path,
    " (", format(rows_written, big.mark = ","), " rows, ",
    format(unique_tracker$get(), big.mark = ","), " unique texts)"
  )
}
