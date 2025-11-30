#!/usr/bin/env Rscript
# drugbank_lean_export.R - Export LEAN tables from dbdataset
# NO explosions - just raw valid combinations from source
# Filters applied based on existing rules in drugbank_generics.R

library(dbdataset)
drugbank <- drugbank

output_dir <- file.path(getwd(), "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("LEAN DrugBank Export (with proper filters)\n")
cat("============================================================\n\n")

# =============================================================================
# SHARED: Build exclusion filter (vet drugs)
# =============================================================================
groups_raw <- drugbank$drugs$groups
groups_raw$group_clean <- tolower(trimws(groups_raw$group))

# Exclude vet-only drugs
vet_ids <- unique(groups_raw$drugbank_id[groups_raw$group_clean == "vet_approved"])
approved_ids <- unique(groups_raw$drugbank_id[groups_raw$group_clean == "approved"])
# Drugs that are ONLY vet (no human approval)
vet_only_ids <- setdiff(vet_ids, approved_ids)

cat(sprintf("Total drugs in DB: %d\n", length(unique(drugbank$drugs$general_information$drugbank_id))))
cat(sprintf("Vet-approved: %d\n", length(vet_ids)))
cat(sprintf("Human-approved: %d\n", length(approved_ids)))
cat(sprintf("Vet-ONLY (excluded): %d\n", length(vet_only_ids)))

filter_vet <- function(df, id_col = "drugbank_id") {
  df[!(df[[id_col]] %in% vet_only_ids), ]
}

# =============================================================================
# 1. GENERICS_LEAN - One row per drug (drugbank_id → name)
#    Source: $drugs$general_information
# =============================================================================
cat("\n[1] generics_lean...\n")
generics <- drugbank$drugs$general_information[, c("drugbank_id", "name", "type", "cas_number", "unii")]
generics$name <- toupper(trimws(generics$name))
generics <- filter_vet(generics)
generics <- unique(generics)
cat(sprintf("    %d rows\n", nrow(generics)))
write.csv(generics, file.path(output_dir, "generics_lean.csv"), row.names = FALSE)

# =============================================================================
# 2. SYNONYMS_LEAN - drugbank_id → synonym
#    Source: $drugs$synonyms
#    Filters: language == "english", coder has allowed values, not solely iupac
# =============================================================================
cat("[2] synonyms_lean...\n")
synonyms <- drugbank$drugs$synonyms
synonyms$synonym <- toupper(trimws(synonyms$synonym))
synonyms$language <- tolower(trimws(synonyms$language))
synonyms$coder <- tolower(trimws(synonyms$coder))

# Filter: English only
synonyms <- synonyms[!is.na(synonyms$language) & grepl("english", synonyms$language, fixed = TRUE), ]

# Filter: coder must exist and not be blank
synonyms <- synonyms[!is.na(synonyms$coder) & synonyms$coder != "", ]

# Allowed coders (INN, USAN, BAN, JAN, etc.)
allowed_coders <- c("inn", "usan", "ban", "jan", "dcj", "usp", "dcit")

# Split coder tokens and filter
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

# Exclude if ONLY iupac
only_iupac <- sapply(synonyms$coder_tokens, function(vals) {
  length(vals) > 0 && all(vals == "iupac")
})

synonyms <- synonyms[has_allowed & !only_iupac, ]
synonyms$coder_tokens <- NULL

# Apply vet filter
synonyms <- filter_vet(synonyms)
synonyms <- unique(synonyms[, c("drugbank_id", "synonym", "coder")])
cat(sprintf("    %d rows\n", nrow(synonyms)))
write.csv(synonyms, file.path(output_dir, "synonyms_lean.csv"), row.names = FALSE)

# =============================================================================
# 3. DOSAGES_LEAN - drugbank_id × form × route × strength (VALID combos)
#    Source: $drugs$dosages
# =============================================================================
cat("[3] dosages_lean...\n")
dosages <- drugbank$drugs$dosages[, c("drugbank_id", "form", "route", "strength")]
dosages$form <- toupper(trimws(dosages$form))
dosages$route <- toupper(trimws(dosages$route))
dosages$strength <- toupper(trimws(dosages$strength))
dosages <- filter_vet(dosages)
dosages <- unique(dosages)
cat(sprintf("    %d rows\n", nrow(dosages)))
write.csv(dosages, file.path(output_dir, "dosages_lean.csv"), row.names = FALSE)

# =============================================================================
# 4. GROUPS_LEAN - drugbank_id → group (for filtering)
#    Source: $drugs$groups
# =============================================================================
cat("[4] groups_lean...\n")
groups <- drugbank$drugs$groups[, c("drugbank_id", "group")]
groups$group <- tolower(trimws(groups$group))
groups <- unique(groups)
cat(sprintf("    %d rows\n", nrow(groups)))
write.csv(groups, file.path(output_dir, "groups_lean.csv"), row.names = FALSE)

# =============================================================================
# 5. BRANDS_LEAN - brand → drugbank_id
#    Source: $drugs$international_brands
# =============================================================================
cat("[5] brands_lean...\n")
brands <- drugbank$drugs$international_brands[, c("drugbank_id", "brand", "company")]
brands$brand <- toupper(trimws(brands$brand))
brands <- filter_vet(brands)
brands <- unique(brands)
cat(sprintf("    %d rows\n", nrow(brands)))
write.csv(brands, file.path(output_dir, "brands_lean.csv"), row.names = FALSE)

# =============================================================================
# 6. SALTS_LEAN - parent drugbank_id → salt info
#    Source: $salts (top-level, not under $drugs)
# =============================================================================
cat("[6] salts_lean...\n")
salts <- drugbank$salts[, c("drugbank_id", "db_salt_id", "name", "cas_number", "unii", "inchikey")]
salts$name <- toupper(trimws(salts$name))
salts <- filter_vet(salts)
salts <- unique(salts)
cat(sprintf("    %d rows\n", nrow(salts)))
write.csv(salts, file.path(output_dir, "salts_lean.csv"), row.names = FALSE)

# =============================================================================
# 7. MIXTURES_LEAN - mixture name, ingredients, drugbank_id
#    Source: $drugs$mixtures
# =============================================================================
cat("[7] mixtures_lean...\n")
mixtures <- drugbank$drugs$mixtures[, c("drugbank_id", "name", "ingredients")]
mixtures$name <- toupper(trimws(mixtures$name))
mixtures$ingredients <- toupper(trimws(mixtures$ingredients))
mixtures <- filter_vet(mixtures)
mixtures <- unique(mixtures)
cat(sprintf("    %d rows\n", nrow(mixtures)))
write.csv(mixtures, file.path(output_dir, "mixtures_lean.csv"), row.names = FALSE)

# =============================================================================
# 8. PRODUCTS_LEAN - drugbank_id × dosage_form × strength × route
#    Source: $products (top-level)
#    Rule: name = generic name if generic==true, name = brand if generic==false
# =============================================================================
cat("[8] products_lean...\n")
products <- drugbank$products[, c("drugbank_id", "name", "labeller", "dosage_form", 
                                   "strength", "route", "generic", "approved", 
                                   "over_the_counter", "country", "source")]
products$dosage_form <- toupper(trimws(products$dosage_form))
products$strength <- toupper(trimws(products$strength))
products$route <- toupper(trimws(products$route))
products$name <- toupper(trimws(products$name))

# Add name_type column based on generic flag
products$name_type <- ifelse(tolower(products$generic) == "true", "generic", "brand")

products <- filter_vet(products)
products <- unique(products)
cat(sprintf("    %d rows\n", nrow(products)))
write.csv(products, file.path(output_dir, "products_lean.csv"), row.names = FALSE)

# =============================================================================
# 9. ATC_LEAN - drugbank_id → atc_code (with hierarchy)
#    Source: $drugs$atc_codes
# =============================================================================
cat("[9] atc_lean...\n")
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
cat("Summary (all vet-only drugs excluded)\n")
cat("============================================================\n")
cat(sprintf("  generics_lean:  %6d (one per drug)\n", nrow(generics)))
cat(sprintf("  synonyms_lean:  %6d (english, allowed coders, not only iupac)\n", nrow(synonyms)))
cat(sprintf("  dosages_lean:   %6d (valid form × route × strength)\n", nrow(dosages)))
cat(sprintf("  groups_lean:    %6d (for filtering)\n", nrow(groups)))
cat(sprintf("  brands_lean:    %6d (international brands)\n", nrow(brands)))
cat(sprintf("  salts_lean:     %6d (salt forms)\n", nrow(salts)))
cat(sprintf("  mixtures_lean:  %6d (mixture ingredients)\n", nrow(mixtures)))
cat(sprintf("  products_lean:  %6d (product details, name_type=generic/brand)\n", nrow(products)))
cat(sprintf("  atc_lean:       %6d (ATC codes with hierarchy)\n", nrow(atc)))
cat("============================================================\n")
cat("Files saved to:", output_dir, "\n")
