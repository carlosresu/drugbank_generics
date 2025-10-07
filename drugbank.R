# ================================
# DrugBank → (generic[s], route, form, dosage) → ATC (all)
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
# Outputs:
#   1) output/generics.csv   (columns: generic, route, form, dosage, atc_code, is_combo)
#   2) output/authority.csv  (columns: text, type, generic, route, form, dosage, atc_code, is_combo, country, source, labeller)
# ================================

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
output_path_generics <- file.path(output_dir, "generics.csv")
output_path_authority <- file.path(output_dir, "authority.csv")

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

# Vector-safe splitter for mixture ingredients
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

# ---------- load ----------
cat("Loading DrugBank slices…\n")
flush.console()
dataset <- drugbank

# Optional: Filter out certain non-therapeutic groups
excluded_groups <- c("experimental", "withdrawn", "illicit", "vet")
grp <- as.data.table(dataset$drugs$groups)
grp[, drugbank_id := as.character(drugbank_id)]
grp[, group_lc := tolower(trimws(group))]
bad_ids <- unique(grp[group_lc %chin% excluded_groups, drugbank_id])
rm(grp)
invisible(gc())

# ---------- ATC (for inclusion + combo flag) ----------
atc_tbl <- unique(as.data.table(dataset$drugs$atc_codes)[
  !is.na(atc_code) & nzchar(atc_code),
  .(drugbank_id, atc_code)
])
atc_tbl[, drugbank_id := as.character(drugbank_id)]
atc_tbl[, last2 := suppressWarnings(as.integer(sub(".*(..)$", "\\1", atc_code)))]
atc_tbl[, is_combo := grepl("^J05AR", atc_code) | (!is.na(last2) & last2 >= 50L)]
combo_flag <- atc_tbl[, .(is_combo = any(is_combo, na.rm = TRUE)), by = drugbank_id]
setkey(atc_tbl, drugbank_id)
setkey(combo_flag, drugbank_id)

# ---------- Build primary name map (generic + clinical synonyms) ----------
gi_raw <- as.data.table(dataset$drugs$general_information)[
  , .(drugbank_id, name)
][!is.na(name) & nzchar(trimws(name))]
gi_raw[, drugbank_id := as.character(drugbank_id)]
gi_raw <- gi_raw[!(drugbank_id %chin% bad_ids)]
gi_raw[, generic := collapse_ws(name)]
gi <- unique(gi_raw[, .(drugbank_id, generic_key = dedupe_key(generic), generic)],
  by = c("drugbank_id", "generic_key")
)[, .(drugbank_id, generic)]
setkey(gi, drugbank_id)

keep_coder_re <- "(?i)\\b(inn|ban|usan|jan|usp)\\b"
syn <- as.data.table(dataset$drugs$synonyms)[
  , .(drugbank_id, synonym, language, coder)
][
  !is.na(coder) & nzchar(trimws(coder)) &
    grepl(keep_coder_re, coder, perl = TRUE) &
    !is.na(language) & grepl("english", language, ignore.case = TRUE) &
    !is.na(synonym) & nzchar(trimws(synonym))
]
syn[, drugbank_id := as.character(drugbank_id)]
syn <- syn[!(drugbank_id %chin% bad_ids)]
syn[, generic := collapse_ws(synonym)]
syn[, generic_key := dedupe_key(generic)]
gi_keys <- gi[, .(generic_key = dedupe_key(generic)), by = drugbank_id]
setkey(gi_keys, drugbank_id, generic_key)
setkey(syn, drugbank_id, generic_key)
syn <- syn[!gi_keys][, .(drugbank_id, generic)]
setkey(syn, drugbank_id)

generic_name_map <- unique(rbindlist(list(gi, syn), use.names = TRUE, fill = TRUE))
generic_name_map <- merge(generic_name_map, combo_flag, by = "drugbank_id", all.x = TRUE)
generic_name_map[is.na(is_combo), is_combo := FALSE]
setkey(generic_name_map, drugbank_id)

# ---------- Route/Form/Dose (dosages + products) ----------
dosages_tbl <- as.data.table(dataset$drugs$dosages)[
  , .(
    drugbank_id = as.character(drugbank_id),
    route = normalize_route_field(route),
    form = trimws(form),
    dosage = tolower(trimws(strength)),
    source = "DOSAGES",
    country = NA_character_,
    labeller = NA_character_,
    fda_application_number = NA_character_,
    ndc_product_code = NA_character_,
    dpd_id = NA_character_,
    ema_product_code = NA_character_,
    ema_ma_number = NA_character_
  )
]
dosages_tbl <- dosages_tbl[!(drugbank_id %chin% bad_ids)]

products_tbl <- as.data.table(dataset$products)[
  , .(
    drugbank_id = as.character(drugbank_id),
    product_name = collapse_ws(name),
    labeller = collapse_ws(labeller),
    route = normalize_route_field(route),
    form = trimws(dosage_form),
    dosage = tolower(trimws(strength)),
    source = trimws(source),
    country = trimws(country),
    fda_application_number = trimws(fda_application_number),
    ndc_product_code = trimws(ndc_product_code),
    dpd_id = trimws(dpd_id),
    ema_product_code = trimws(ema_product_code),
    ema_ma_number = trimws(ema_ma_number)
  )
]
products_tbl <- products_tbl[!(drugbank_id %chin% bad_ids)]

route_form <- unique(rbindlist(
  list(
    dosages_tbl[, .(
      drugbank_id, route, form, dosage,
      source, country, labeller,
      fda_application_number, ndc_product_code,
      dpd_id, ema_product_code, ema_ma_number
    )],
    products_tbl[, .(
      drugbank_id, route, form, dosage,
      source, country, labeller,
      fda_application_number, ndc_product_code,
      dpd_id, ema_product_code, ema_ma_number
    )]
  ),
  use.names = TRUE, fill = TRUE
))

# --- Fix progressr steps to match loop iterations exactly ---
{
  n <- nrow(route_form)
  if (n > 0) {
    chunk_size <- 50000L
    starts <- seq.int(1L, n, by = chunk_size)
    progressr::with_progress({
      p <- progressr::progressor(steps = length(starts), label = "Simplifying `form`")
      for (start in starts) {
        end <- min(start + chunk_size - 1L, n)
        route_form[start:end, form := simplify_form_vec(form)]
        p(message = sprintf("forms %d–%d / %d", start, end, n))
        if ((end %% (chunk_size * 4L)) == 0L) invisible(gc(FALSE))
      }
    })
  }
}
setkey(route_form, drugbank_id)

# ---------- Assemble GENERICS: keep ALL ATC codes (no ranking) ----------
ga <- merge(generic_name_map, atc_tbl[, .(drugbank_id, atc_code)], by = "drugbank_id", allow.cartesian = TRUE)
gaf <- merge(ga, route_form, by = "drugbank_id", all.x = TRUE, allow.cartesian = TRUE)
gaf <- gaf[, .(generic, is_combo, route, form, dosage, atc_code)]

# ---------- Explode route & form to true long format ----------
explode_long <- function(DT) {
  tmp <- copy(DT)
  tmp[, route_list := lapply(route, split_semicolon)]
  tmp[, form_list := lapply(form, split_semicolon)]
  tmp[, rid := .I]
  out <- tmp[,
    {
      rl <- route_list[[1L]]
      if (length(rl) == 0L) rl <- NA_character_
      fl <- form_list[[1L]]
      if (length(fl) == 0L) fl <- NA_character_
      CJ(route = rl, form = fl, unique = TRUE)
    },
    by = .(generic, is_combo, dosage, atc_code, rid)
  ]
  out[, rid := NULL]
  out[route == "", route := NA_character_]
  out[form == "", form := NA_character_]
  out[, route := ifelse(is.na(route), NA_character_, tolower(trimws(route)))]
  out[, form := ifelse(is.na(form), NA_character_, simplify_form_vec(form))]
  unique(out[, .(generic, route, form, dosage, atc_code, is_combo)])
}

# --- Time + explicit post-explode progress message ---
cat("Exploding route/form combinations…\n")
flush.console()
t0 <- proc.time()[["elapsed"]]
final_long <- explode_long(gaf)
t1 <- proc.time()[["elapsed"]]
cat(sprintf(
  "Explode complete: %s rows in %.1fs\n",
  format(nrow(final_long), big.mark = ","), (t1 - t0)
))
flush.console()

# ---------- Save #1 ----------
setorder(final_long, generic, route, form, atc_code, dosage)
fwrite(final_long, output_path_generics)
message("Wrote: ", output_path_generics)

# =====================================================================
#                Authoritative text → ATC vocabulary (SOA)
# =====================================================================

primary_gi <- as.data.table(dataset$drugs$general_information)[
  , .(drugbank_id = as.character(drugbank_id), primary_generic = collapse_ws(name))
][!(drugbank_id %chin% bad_ids)]
setkey(primary_gi, drugbank_id)

# 1) generics & synonyms as text entries
auth_generics <- unique(rbindlist(list(
  gi[, .(drugbank_id, text = generic, type = "generic")],
  syn[, .(drugbank_id, text = generic, type = "synonym")]
), use.names = TRUE))

# 2) brands
brands <- as.data.table(dataset$drugs$international_brands)[
  , .(
    drugbank_id = as.character(drugbank_id),
    text = collapse_ws(brand),
    type = "brand"
  )
][!is.na(text) & nzchar(text)]
brands <- brands[!(drugbank_id %chin% bad_ids)]

# 3) products (keep own route/form/dose; already memory-bounded)
prod_text <- products_tbl[, .(
  drugbank_id,
  text = product_name, type = "product",
  route, form, dosage, country, source, labeller
)][!is.na(text) & nzchar(text)]

# 4) mixtures (name + split ingredients; vector-safe)
mix <- as.data.table(dataset$drugs$mixtures)[
  , .(
    drugbank_id = as.character(drugbank_id),
    name = collapse_ws(name),
    ingredients = ingredients
  ) # keep raw; split later with vector-safe helper
]
mix <- mix[!(drugbank_id %chin% bad_ids)]

mix_name <- mix[
  !is.na(name) & nzchar(name),
  .(drugbank_id, text = name, type = "mixture")
]

mix_ing <- mix[!is.na(ingredients) & nzchar(ingredients)]
if (nrow(mix_ing)) {
  mix_ing <- mix_ing[, .(text = split_ingredients(ingredients)), by = .(drugbank_id)]
  mix_ing <- mix_ing[!is.na(text) & nzchar(text)]
  mix_ing[, type := "ingredient"]
} else {
  mix_ing <- mix[0, .(drugbank_id, text = character(), type = character())]
}

# -------------------- MEMORY-SAFE AUTHORITY BUILD --------------------
# Combine all text sources (no structure yet)
authority_texts <- rbindlist(list(auth_generics, brands, mix_name, mix_ing, prod_text),
  use.names = TRUE, fill = TRUE
)
# Attach primary GI for labeling
authority_texts <- merge(authority_texts, primary_gi, by = "drugbank_id", all.x = TRUE)

# Precompute a MEMORY-SAFE per-drug collapsed route/form/dose pool
#   - collapse unique values per drugbank_id into semicolon strings
#   - avoids the huge cartesian of (text × all route/form combos)
collapse_unique <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) {
    return(NA_character_)
  }
  paste(unique(x), collapse = "; ")
}
collapse_split_unique <- function(x) {
  # split on semicolons across the vector, normalize case/space, and collapse unique
  if (!length(x)) {
    return(NA_character_)
  }
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) {
    return(NA_character_)
  }
  pieces <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  pieces <- tolower(trimws(pieces))
  pieces <- pieces[nzchar(pieces)]
  if (!length(pieces)) {
    return(NA_character_)
  }
  paste(unique(pieces), collapse = "; ")
}

cat("Collapsing route/form/dose pools per drug…\n")
flush.console()
t_rf0 <- proc.time()[["elapsed"]]
rf_agg <- route_form[, .(
  route = collapse_split_unique(route),
  form = collapse_split_unique(form),
  dosage = collapse_unique(dosage),
  country = collapse_unique(country),
  source = collapse_unique(source),
  labeller = collapse_unique(labeller)
), by = .(drugbank_id)]
t_rf1 <- proc.time()[["elapsed"]]
cat(sprintf(
  "Collapsed pools for %s drugs in %.1fs\n",
  format(nrow(rf_agg), big.mark = ","), (t_rf1 - t_rf0)
))
setkey(rf_agg, drugbank_id)

# Split into product vs non-product rows to avoid cartesian bloat
is_prod <- authority_texts$type == "product"

# Non-product: attach collapsed route/form/dose (single row per drug), THEN attach ATC
nonprod <- merge(
  authority_texts[
    !is_prod,
    .(drugbank_id, text, type, primary_generic)
  ],
  rf_agg,
  by = "drugbank_id", all.x = TRUE, allow.cartesian = FALSE
)

# Attach ATC & is_combo (this may fan-out per atc_code, but not per route/form)
nonprod <- merge(nonprod,
  atc_tbl[, .(drugbank_id, atc_code, is_combo)],
  by = "drugbank_id", allow.cartesian = TRUE
)

# Product rows already have their own route/form/dose; just attach ATC
prod <- merge(
  authority_texts[
    is_prod,
    .(
      drugbank_id, text, type, primary_generic,
      route, form, dosage, country, source, labeller
    )
  ],
  atc_tbl[, .(drugbank_id, atc_code, is_combo)],
  by = "drugbank_id", allow.cartesian = TRUE
)

# Combine back; keep unified schema
authority_final <- rbindlist(list(
  nonprod[, .(text, type, generic = primary_generic, route, form, dosage, atc_code, is_combo, country, source, labeller)],
  prod[, .(text, type, generic = primary_generic, route, form, dosage, atc_code, is_combo, country, source, labeller)]
), use.names = TRUE, fill = TRUE)

# Light de-dup & ordering (memory-friendly)
setkey(authority_final, NULL)
setorder(authority_final, text, type, atc_code, route, form, dosage)
authority_final <- unique(authority_final)

cat(sprintf(
  "Authority table ready: %s rows, %s unique texts.\n",
  format(nrow(authority_final), big.mark = ","),
  format(length(unique(authority_final$text)), big.mark = ",")
))
flush.console()

# ---------- Save #2 ----------
fwrite(authority_final, output_path_authority)
message("Wrote: ", output_path_authority)
