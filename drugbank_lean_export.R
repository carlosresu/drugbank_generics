#!/usr/bin/env Rscript
# drugbank_lean_export.R - Export LEAN tables from dbdataset
# NO explosions - just raw valid combinations from source
# Applies normalization/standardization rules from existing scripts

library(dbdataset)
drugbank <- drugbank

output_dir <- file.path(getwd(), "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("LEAN DrugBank Export (with normalization)\n")
cat("============================================================\n\n")

# =============================================================================
# SHARED: Helper functions
# =============================================================================

# Collapse whitespace
collapse_ws <- function(x) {
  ifelse(is.na(x), NA_character_, trimws(gsub("\\s+", " ", as.character(x))))
}

# Normalize to uppercase, trim
normalize_name <- function(x) {
  toupper(trimws(x))
}

# Normalize key (lowercase, clean)
normalize_key <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("\\s+", " ", x)
  x <- gsub("[^a-z0-9 ]", "", x)
  trimws(x)
}

# Split ingredients on + (from drugbank_mixtures.R logic)
split_ingredients <- function(value) {
  val <- collapse_ws(value)
  if (is.na(val) || !nzchar(val)) return(character())
  # Remove parenthetical content
  repeat {
    new_val <- gsub("(?<!\\S)\\([^()]*\\)(?=\\s|$)", " ", val, perl = TRUE)
    if (identical(new_val, val)) break
    val <- new_val
  }
  val <- gsub("\\s+", " ", val, perl = TRUE)
  # Clean up commas before +
  if (grepl("\\+", val)) {
    val <- gsub(",\\s[^+]+(?=\\s*\\+)", "", val, perl = TRUE)
  }
  val <- gsub(",\\s[^+]+$", "", val, perl = TRUE)
  # Split on +
  parts <- if (grepl("\\+", val)) {
    unlist(strsplit(val, "\\+", perl = TRUE), use.names = FALSE)
  } else {
    val
  }
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  unique(parts)
}

# Collapse to pipe-separated
collapse_pipe <- function(values) {
  vals <- unique(values[!is.na(values) & nzchar(values)])
  if (!length(vals)) return(NA_character_)
  paste(vals, collapse = "|")
}

# =============================================================================
# SHARED: Build exclusion filter (vet-only drugs)
# =============================================================================
groups_raw <- drugbank$drugs$groups
groups_raw$group_clean <- tolower(trimws(groups_raw$group))

# Exclude vet-only drugs (vet_approved but NOT human approved)
vet_ids <- unique(groups_raw$drugbank_id[groups_raw$group_clean == "vet_approved"])
approved_ids <- unique(groups_raw$drugbank_id[groups_raw$group_clean == "approved"])
vet_only_ids <- setdiff(vet_ids, approved_ids)

cat(sprintf("Total drugs in DB: %d\n", length(unique(drugbank$drugs$general_information$drugbank_id))))
cat(sprintf("Vet-ONLY excluded: %d\n", length(vet_only_ids)))

filter_vet <- function(df, id_col = "drugbank_id") {
  df[!(df[[id_col]] %in% vet_only_ids), ]
}

# =============================================================================
# 1. GENERICS_LEAN - One row per drug (drugbank_id → name)
# =============================================================================
cat("\n[1] generics_lean...\n")
generics <- drugbank$drugs$general_information[, c("drugbank_id", "name", "type", "cas_number", "unii")]
generics$name <- normalize_name(generics$name)
generics$name_key <- normalize_key(generics$name)
generics <- filter_vet(generics)
generics <- unique(generics)
cat(sprintf("    %d rows\n", nrow(generics)))
write.csv(generics, file.path(output_dir, "generics_lean.csv"), row.names = FALSE)

# =============================================================================
# 2. SYNONYMS_LEAN - drugbank_id → synonym
#    Filters: english, has allowed coder, not ONLY iupac (iupac ok if others present)
# =============================================================================
cat("[2] synonyms_lean...\n")
synonyms <- drugbank$drugs$synonyms
synonyms$synonym <- normalize_name(synonyms$synonym)
synonyms$synonym_key <- normalize_key(synonyms$synonym)
synonyms$language <- tolower(trimws(synonyms$language))
synonyms$coder <- tolower(trimws(synonyms$coder))

# Filter: English only
synonyms <- synonyms[!is.na(synonyms$language) & grepl("english", synonyms$language, fixed = TRUE), ]

# Filter: coder must exist
synonyms <- synonyms[!is.na(synonyms$coder) & synonyms$coder != "", ]

# Allowed coders
allowed_coders <- c("inn", "usan", "ban", "jan", "dcj", "usp", "dcit")

# Split coder tokens
split_coder <- function(x) {
  tokens <- unlist(strsplit(x, "[,;/\\s]+"))
  tokens <- tolower(trimws(tokens))
  tokens[tokens != ""]
}
synonyms$coder_tokens <- lapply(synonyms$coder, split_coder)

# Keep if has at least one allowed coder
has_allowed <- sapply(synonyms$coder_tokens, function(vals) {
  length(vals) > 0 && any(vals %in% allowed_coders)
})

# Exclude ONLY if coder is solely "iupac" (iupac + others is OK)
only_iupac <- sapply(synonyms$coder_tokens, function(vals) {
  length(vals) > 0 && length(vals) == 1 && vals[1] == "iupac"
})

synonyms <- synonyms[has_allowed | !only_iupac, ]
synonyms$coder_tokens <- NULL

synonyms <- filter_vet(synonyms)
synonyms <- unique(synonyms[, c("drugbank_id", "synonym", "synonym_key", "coder")])
cat(sprintf("    %d rows\n", nrow(synonyms)))
write.csv(synonyms, file.path(output_dir, "synonyms_lean.csv"), row.names = FALSE)

# =============================================================================
# 3. DOSAGES_LEAN - drugbank_id × form × route × strength (VALID combos)
# =============================================================================
cat("[3] dosages_lean...\n")
dosages <- drugbank$drugs$dosages[, c("drugbank_id", "form", "route", "strength")]
dosages$form <- normalize_name(dosages$form)
dosages$route <- normalize_name(dosages$route)
dosages$strength <- normalize_name(dosages$strength)
dosages <- filter_vet(dosages)
dosages <- unique(dosages)
cat(sprintf("    %d rows\n", nrow(dosages)))
write.csv(dosages, file.path(output_dir, "dosages_lean.csv"), row.names = FALSE)

# =============================================================================
# 4. BRANDS_LEAN - brand → drugbank_id
# =============================================================================
cat("[4] brands_lean...\n")
brands <- drugbank$drugs$international_brands[, c("drugbank_id", "brand", "company")]
brands$brand <- normalize_name(brands$brand)
brands$brand_key <- normalize_key(brands$brand)
brands <- filter_vet(brands)
brands <- unique(brands)
cat(sprintf("    %d rows\n", nrow(brands)))
write.csv(brands, file.path(output_dir, "brands_lean.csv"), row.names = FALSE)

# =============================================================================
# 5. SALTS_LEAN - parent drugbank_id → salt info
# =============================================================================
cat("[5] salts_lean...\n")
salts <- drugbank$salts[, c("drugbank_id", "db_salt_id", "name", "cas_number", "unii", "inchikey")]
salts$name <- normalize_name(salts$name)
salts$name_key <- normalize_key(salts$name)
salts <- filter_vet(salts)
salts <- unique(salts)
cat(sprintf("    %d rows\n", nrow(salts)))
write.csv(salts, file.path(output_dir, "salts_lean.csv"), row.names = FALSE)

# =============================================================================
# 6. MIXTURES_LEAN - with split components
#    Splits ingredients on + and creates component keys
# =============================================================================
cat("[6] mixtures_lean...\n")
mixtures <- drugbank$drugs$mixtures[, c("drugbank_id", "name", "ingredients")]
mixtures$mixture_name <- normalize_name(mixtures$name)
mixtures$mixture_name_key <- normalize_key(mixtures$name)
mixtures$ingredients_raw <- normalize_name(mixtures$ingredients)

# Split ingredients into components
mixtures$component_list <- lapply(mixtures$ingredients_raw, split_ingredients)
mixtures$component_generics <- sapply(mixtures$component_list, function(parts) {
  collapse_pipe(normalize_name(parts))
})
mixtures$component_keys <- sapply(mixtures$component_list, function(parts) {
  collapse_pipe(normalize_key(parts))
})
mixtures$component_count <- sapply(mixtures$component_list, length)

# Create sorted component key for matching
mixtures$component_key_sorted <- sapply(mixtures$component_list, function(parts) {
  keys <- sort(normalize_key(parts))
  collapse_pipe(keys)
})

mixtures$component_list <- NULL
mixtures$name <- NULL
mixtures$ingredients <- NULL

mixtures <- filter_vet(mixtures)
mixtures <- unique(mixtures)
cat(sprintf("    %d rows\n", nrow(mixtures)))
write.csv(mixtures, file.path(output_dir, "mixtures_lean.csv"), row.names = FALSE)

# =============================================================================
# 7. PRODUCTS_LEAN - drugbank_id × dosage_form × strength × route
#    name_type = "generic" if generic==true, else "brand"
# =============================================================================
cat("[7] products_lean...\n")
products <- drugbank$products[, c("drugbank_id", "name", "labeller", "dosage_form", 
                                   "strength", "route", "generic", "approved", 
                                   "over_the_counter", "country", "source")]
products$dosage_form <- normalize_name(products$dosage_form)
products$strength <- normalize_name(products$strength)
products$route <- normalize_name(products$route)
products$name <- normalize_name(products$name)
products$name_key <- normalize_key(products$name)

# Add name_type column based on generic flag
products$name_type <- ifelse(tolower(products$generic) == "true", "generic", "brand")

products <- filter_vet(products)
products <- unique(products)
cat(sprintf("    %d rows\n", nrow(products)))
write.csv(products, file.path(output_dir, "products_lean.csv"), row.names = FALSE)

# =============================================================================
# 8. ATC_LEAN - drugbank_id → atc_code (with hierarchy)
# =============================================================================
cat("[8] atc_lean...\n")
atc <- drugbank$drugs$atc_codes[, c("drugbank_id", "atc_code", 
                                     "level_1", "code_1", 
                                     "level_2", "code_2",
                                     "level_3", "code_3",
                                     "level_4", "code_4")]
atc$atc_code <- trimws(atc$atc_code)
atc <- filter_vet(atc)
atc <- atc[!is.na(atc$atc_code) & atc$atc_code != "", ]
atc <- unique(atc)
cat(sprintf("    %d rows\n", nrow(atc)))
write.csv(atc, file.path(output_dir, "atc_lean.csv"), row.names = FALSE)

# =============================================================================
# Summary
# =============================================================================
cat("\n============================================================\n")
cat("Summary (vet-only excluded, normalized)\n")
cat("============================================================\n")
cat(sprintf("  generics_lean:  %6d (one per drug, with name_key)\n", nrow(generics)))
cat(sprintf("  synonyms_lean:  %6d (english, allowed coders, iupac ok if not alone)\n", nrow(synonyms)))
cat(sprintf("  dosages_lean:   %6d (valid form × route × strength)\n", nrow(dosages)))
cat(sprintf("  brands_lean:    %6d (with brand_key)\n", nrow(brands)))
cat(sprintf("  salts_lean:     %6d (with name_key)\n", nrow(salts)))
cat(sprintf("  mixtures_lean:  %6d (with split components, component_key_sorted)\n", nrow(mixtures)))
cat(sprintf("  products_lean:  %6d (name_type=generic/brand, with name_key)\n", nrow(products)))
cat(sprintf("  atc_lean:       %6d (with hierarchy)\n", nrow(atc)))
cat("============================================================\n")
cat("Files saved to:", output_dir, "\n")
