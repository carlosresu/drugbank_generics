# ================================
# DrugBank → (generic[s], route, form, dosage) → clinically-preferred ATC
# Map-Reduce, memory-aware (data.table + future/furrr + progressr)
# Requirements implemented:
# 1) Use drugbank_id internally only; exclude from final output.
# 2) Keep ALL names from general_information$name.
# 3) Add synonyms only where coder present, language contains "english",
#    and synonym (normalized) not already in #2 set.
# 4) Normalize/simplify `form`: keep main term; keep comma-suffix only if it contains "release".
# 5) Normalize & explode `route` and `form` by ';' to a true long format.
# Final columns: generic, route, form, dosage, atc_code (all lowercase column names)
# ================================

ensure_installed <- function(packages, repos = "https://cloud.r-project.org") {
  for (pkg in packages) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = repos)
}

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  getwd()
}

# --- deps ---
ensure_installed(c("data.table", "future", "furrr", "stringi", "progressr", "devtools"))
if (!requireNamespace("dbdataset", quietly = TRUE)) {
  devtools::install_github("interstellar-Consultation-Services/dbdataset", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(future)
  library(furrr)
  library(stringi)
  library(progressr)
  library(dbdataset)
})

data.table::setDTthreads(percent = 100)

handlers(global = TRUE)
handlers("cli")   # or handlers("txtprogressbar")

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, "generics.csv")

# ---------- helpers ----------
normalize_name <- function(x) {
  x <- trimws(x)
  x <- stri_trans_general(x, "Latin-ASCII")
  x[is.na(x)] <- ""
  tolower(x)
}

# Normalize raw names into (possibly) combo generics split by “and/with/plus”, “/”, “+”, “;”
standardize_combo_generic_one <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(NA_character_)
  x <- gsub("\\s+(and|with|plus)\\s+", "|", x, ignore.case = TRUE)
  x <- gsub("\\s*/\\s*|\\s*,\\s*|\\s*\\+\\s*|\\s*;\\s*", "|", x, perl = TRUE)
  parts <- strsplit(x, "\\|", fixed = FALSE)[[1]]
  if (!length(parts)) parts <- x
  parts <- normalize_name(parts)
  parts <- stri_replace_all_regex(parts, "[^a-z ]", " ")
  parts <- stri_trim_both(parts)
  parts <- unique(parts[nzchar(parts)])
  if (!length(parts)) return(NA_character_)
  if (length(parts) == 1) return(parts)
  paste(sort(parts), collapse = "; ")
}
standardize_combo_generic <- function(x) vapply(x, standardize_combo_generic_one, character(1))

# Form simplification (vectorized):
# - lowercase everything
# - if comma suffix contains 'release', keep "main, that suffix"; otherwise keep main only
simplify_form_one <- function(f) {
  if (is.null(f) || is.na(f) || !nzchar(f)) return(NA_character_)
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

# Normalize route: lower-case + trim; explode later by ';'
normalize_route_field <- function(v) {
  v <- trimws(v)
  v[is.na(v)] <- NA_character_
  tolower(v)
}

# Fast cross join helper
CJ_row <- function(A, B) {
  A[, `__tmp__` := 1L]
  B[, `__tmp__` := 1L]
  on.exit({ A[, `__tmp__` := NULL]; B[, `__tmp__` := NULL] }, add = TRUE)
  merge(A, B, by = "__tmp__", allow.cartesian = TRUE)[, `__tmp__` := NULL][]
}

split_semicolon <- function(s) {
  if (is.na(s) || !nzchar(s)) return(character(0))
  parts <- strsplit(s, ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  tolower(parts)
}

ensure_list_nonempty <- function(lst) {
  lst[lengths(lst) == 0L] <- list(NA_character_)
  lst
}

cat("Loading DrugBank slices…\n")
dataset <- drugbank

# ---------- Build generic map per requirements (GI + filtered Synonyms only) ----------
gi <- as.data.table(dataset$drugs$general_information)[
  , .(drugbank_id, name)
]
gi[, generic := standardize_combo_generic(name)]
gi <- gi[!is.na(generic) & nzchar(generic), .(drugbank_id, generic)]
setkey(gi, drugbank_id)

gi_names_norm <- unique(standardize_combo_generic(as.character(dataset$drugs$general_information$name)))
gi_names_norm <- gi_names_norm[!is.na(gi_names_norm) & nzchar(gi_names_norm)]

syn <- as.data.table(dataset$drugs$synonyms)[
  , .(drugbank_id, synonym, language, coder)
][
  !is.na(coder) & trimws(coder) != "" &
    !is.na(language) & grepl("english", tolower(language)) &
    !is.na(synonym) & nzchar(trimws(synonym))
]
syn[, generic := standardize_combo_generic(synonym)]
syn <- syn[!is.na(generic) & nzchar(generic)]
syn <- syn[!(generic %chin% gi_names_norm), .(drugbank_id, generic)]

generic_name_map <- unique(rbindlist(list(gi, syn), use.names = TRUE, fill = TRUE))
rm(gi, syn, gi_names_norm); invisible(gc())

# ---------- ATC ----------
atc_tbl <- unique(as.data.table(dataset$drugs$atc_codes)[
  !is.na(atc_code) & nzchar(atc_code),
  .(drugbank_id, atc_code)
])

# ---------- Route/Form/Dose ----------
dosages_tbl <- as.data.table(dataset$drugs$dosages)[
  , .(
    drugbank_id,
    route  = normalize_route_field(route),
    form   = trimws(form),
    dosage = tolower(trimws(strength)),
    source = "DOSAGES",
    country = NA_character_,
    fda_application_number = NA_character_,
    ndc_product_code       = NA_character_,
    dpd_id                 = NA_character_,
    ema_product_code       = NA_character_,
    ema_ma_number          = NA_character_
  )
]

products_tbl <- as.data.table(dataset$products)[
  , .(
    drugbank_id,
    route  = normalize_route_field(route),
    form   = trimws(dosage_form),
    dosage = tolower(trimws(strength)),
    source, country,
    fda_application_number = trimws(fda_application_number),
    ndc_product_code       = trimws(ndc_product_code),
    dpd_id                 = trimws(dpd_id),
    ema_product_code       = trimws(ema_product_code),
    ema_ma_number          = trimws(ema_ma_number)
  )
]

route_form <- unique(rbindlist(list(dosages_tbl, products_tbl), use.names = TRUE, fill = TRUE))[
  , .(drugbank_id, route, form, dosage,
      source, country, fda_application_number, ndc_product_code,
      dpd_id, ema_product_code, ema_ma_number)
]
rm(dosages_tbl, products_tbl); invisible(gc())

# Simplify/normalize form BEFORE explode — with progress updates
{
  n <- nrow(route_form)
  chunk_size <- 50000L  # adjust if you want finer/coarser updates
  
  progressr::with_progress({
    p <- progressr::progressor(steps = ceiling(n / chunk_size), label = "Simplifying `form`")
    for (start in seq.int(1L, n, by = chunk_size)) {
      end <- min(start + chunk_size - 1L, n)
      # in-place update on this chunk
      route_form[start:end, form := simplify_form_vec(form)]
      p(message = sprintf("forms %d–%d / %d", start, end, n))
      if ((end %% (chunk_size * 4L)) == 0L) invisible(gc(FALSE))  # occasional GC to keep memory healthy
    }
  })
}

# ---------- Index ----------
setkey(generic_name_map, drugbank_id)
setkey(atc_tbl,          drugbank_id)
setkey(route_form,       drugbank_id)

ids <- unique(atc_tbl$drugbank_id)

# ---------- Worker ----------
compute_support_for_ids <- function(id_vec) {
  res_list <- vector("list", length(id_vec))
  j <- 0L
  for (id in id_vec) {
    gens <- generic_name_map[.(id), .(generic)]
    if (nrow(gens) == 0L) next
    atcs <- atc_tbl[.(id), .(atc_code)]
    if (nrow(atcs) == 0L) next
    rf <- route_form[.(id),
                     .(route, form, dosage, source, country,
                       fda_application_number, ndc_product_code, dpd_id,
                       ema_product_code, ema_ma_number)]
    if (nrow(rf) == 0L) {
      rf <- data.table(
        route = NA_character_, form = NA_character_, dosage = NA_character_,
        source = NA_character_, country = NA_character_,
        fda_application_number = NA_character_, ndc_product_code = NA_character_,
        dpd_id = NA_character_, ema_product_code = NA_character_, ema_ma_number = NA_character_
      )
    }
    ga  <- CJ_row(copy(gens), copy(atcs))
    gaf <- CJ_row(ga, rf)
    gaf[, last2 := suppressWarnings(as.integer(sub(".*(..)$", "\\1", atc_code)))]
    gaf[, is_atc_combo := grepl("^J05AR", atc_code) | (!is.na(last2) & last2 >= 50L)]
    gaf[, is_combo_generic := grepl(";", generic, fixed = TRUE)]
    gaf[, support_key := fifelse(
      !is.na(ndc_product_code) & nzchar(ndc_product_code), ndc_product_code,
      fifelse(!is.na(fda_application_number) & nzchar(fda_application_number), fda_application_number,
              fifelse(!is.na(dpd_id) & nzchar(dpd_id), dpd_id,
                      fifelse(!is.na(ema_product_code) & nzchar(ema_product_code), ema_product_code,
                              fifelse(!is.na(ema_ma_number) & nzchar(ema_ma_number), ema_ma_number,
                                      paste0(fifelse(is.na(source), "", source), "|",
                                             fifelse(is.na(country), "", country), "|",
                                             fifelse(is.na(dosage), "", dosage))
                              ))))
    )]
    agg <- gaf[, .(
      support_count = uniqueN(support_key, na.rm = TRUE),
      ids_union     = id
    ), by = .(generic, route, form, dosage, atc_code, is_atc_combo, is_combo_generic)]
    j <- j + 1L
    res_list[[j]] <- agg
  }
  if (j == 0L) return(data.table())
  rbindlist(res_list[seq_len(j)], use.names = TRUE, fill = TRUE)
}

# ---------- Chunking & Parallel plan ----------
chunk_size <- 400L
chunks <- split(ids, ceiling(seq_along(ids) / chunk_size))

n_workers <- max(1L, future::availableCores() - 1L)
future::plan(multisession, workers = n_workers)

cat("Computing per-chunk support in parallel…\n")
partial_supports <- with_progress({
  p_chunks <- progressor(steps = length(chunks), label = "Chunks")
  out <- furrr::future_map(
    chunks,
    function(id_vec) {
      res <- compute_support_for_ids(id_vec)
      p_chunks()
      res
    },
    .options  = furrr::furrr_options(seed = TRUE),
    .progress = FALSE
  )
  out
})

future::plan(sequential)
rm(chunks, ids); invisible(gc())

# ---------- Reduce across chunks ----------
cat("Reducing supports across chunks…\n")
support_tbl <- rbindlist(partial_supports, use.names = TRUE, fill = TRUE)
rm(partial_supports); invisible(gc())

union_ids <- function(x) paste(sort(unique(unlist(strsplit(x, ";", fixed = TRUE)))), collapse = ";")
support_tbl <- support_tbl[
  , .(
    support_count = sum(support_count, na.rm = TRUE),
    ids_union     = union_ids(ids_union)
  ),
  by = .(generic, route, form, dosage, atc_code, is_atc_combo, is_combo_generic)
]

# ---------- Rank & select best ATC ----------
support_tbl[, suffix := suppressWarnings(as.integer(sub(".*(..)$", "\\1", atc_code)))]
support_tbl[, pref_combo := fifelse(is_combo_generic == is_atc_combo, 1L, 2L)]
support_tbl[, pref_base_mono := fifelse(!is_combo_generic & !is_atc_combo & !is.na(suffix) & suffix == 1L, 1L,
                                        fifelse(!is_combo_generic & !is_atc_combo, 2L, 3L))]

setorder(support_tbl, generic, route, form, dosage,
         pref_combo, pref_base_mono, -support_count, suffix, atc_code)

best_atc <- support_tbl[, .SD[1L], by = .(generic, route, form, dosage)][
  , .(generic, route, form, dosage, atc_code)
]

# ---------- Explode route & form to true long format (robust cartesian) ----------
explode_long <- function(DT) {
  tmp <- copy(DT)
  
  # create list-cols of routes/forms per row
  tmp[, route_list := lapply(route, split_semicolon)]
  tmp[, form_list  := lapply(form,  split_semicolon)]
  tmp[, rid := .I]
  
  # build cartesian product per row (generic, dosage, atc_code, rid)
  out <- tmp[, {
    rl <- route_list[[1L]]; if (length(rl) == 0L) rl <- NA_character_
    fl <- form_list[[1L]];  if (length(fl) == 0L) fl <- NA_character_
    CJ(route = rl, form = fl, unique = TRUE)
  }, by = .(generic, dosage, atc_code, rid)]
  
  out[, rid := NULL]
  
  # clean & normalize
  out[route == "", route := NA_character_]
  out[form  == "", form  := NA_character_]
  
  out[, route := ifelse(is.na(route), NA_character_, tolower(trimws(route)))]
  out[, form  := ifelse(is.na(form),  NA_character_, simplify_form_vec(form))]
  
  unique(out[, .(generic, route, form, dosage, atc_code)])
}

final_long <- explode_long(best_atc)

# ---------- Save ----------
setorder(final_long, generic, route, form, atc_code, dosage)
fwrite(final_long, output_path)
message("Wrote: ", output_path)