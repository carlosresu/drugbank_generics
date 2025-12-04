#!/usr/bin/env Rscript
# drugbank_lean_export.R - Export LEAN tables from dbdataset
# NO explosions - just raw valid combinations from source
# Applies normalization/standardization rules from existing scripts
# Exports CSV only (CSV-first policy per AGENTS.md)

library(dbdataset)

drugbank <- drugbank

output_dir <- file.path(getwd(), "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: write CSV only (canonical format)
write_csv_only <- function(df, basename) {
  csv_path <- file.path(output_dir, paste0(basename, ".csv"))
  write.csv(df, csv_path, row.names = FALSE)
}

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
  parts[order(tolower(parts), parts)] |> unique()  # consistent ordering
}

# Unique with consistent ordering (lowercase sort, then unique)
unique_canonical <- function(values) {
  vals <- values[!is.na(values) & nzchar(values)]
  if (!length(vals)) return(character())
  vals[order(tolower(vals), vals)] |> unique()
}

# Collapse to pipe-separated
collapse_pipe <- function(values) {
  vals <- unique_canonical(values)
  if (!length(vals)) return(NA_character_)
  paste(vals, collapse = "|")
}

# Normalize brand name (remove trademark symbols)
normalize_brand <- function(x) {
  val <- collapse_ws(x)
  if (is.na(val)) return(NA_character_)
  val <- gsub("[®™©]", "", val, perl = TRUE)
  trimws(gsub("\\s+", " ", val))
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
# Standardize: name -> generic_name
names(generics)[names(generics) == "name"] <- "generic_name"
generics$generic_name <- normalize_name(generics$generic_name)
generics$name_key <- normalize_key(generics$generic_name)
generics <- filter_vet(generics)
generics <- unique(generics)
cat(sprintf("    %d rows\n", nrow(generics)))
write_csv_only(generics, "generics_lean")

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
# Standardize: synonym -> synonyms (plural)
names(synonyms)[names(synonyms) == "synonym"] <- "synonyms"
names(synonyms)[names(synonyms) == "synonym_key"] <- "synonyms_key"
cat(sprintf("    %d rows\n", nrow(synonyms)))
write_csv_only(synonyms, "synonyms_lean")

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
write_csv_only(dosages, "dosages_lean")

# =============================================================================
# 4. BRANDS_LEAN - brand → drugbank_id (with trademark removal)
# =============================================================================
cat("[4] brands_lean...\n")
brands <- drugbank$drugs$international_brands[, c("drugbank_id", "brand", "company")]
# Standardize: brand -> brand_name
names(brands)[names(brands) == "brand"] <- "brand_name"
brands$brand_name <- sapply(brands$brand_name, function(x) normalize_name(normalize_brand(x)))
brands$brand_key <- normalize_key(brands$brand_name)
brands <- filter_vet(brands)
brands <- unique(brands)
cat(sprintf("    %d rows\n", nrow(brands)))
write_csv_only(brands, "brands_lean")

# =============================================================================
# 5. SALTS_LEAN - parent drugbank_id → salt info
# =============================================================================
cat("[5] salts_lean...\n")
salts <- drugbank$salts[, c("drugbank_id", "db_salt_id", "name", "cas_number", "unii", "inchikey")]
# Standardize: name -> salt_name (to distinguish from generic_name)
names(salts)[names(salts) == "name"] <- "salt_name"
salts$salt_name <- normalize_name(salts$salt_name)
salts$name_key <- normalize_key(salts$salt_name)
salts <- filter_vet(salts)
salts <- unique(salts)
cat(sprintf("    %d rows\n", nrow(salts)))
write_csv_only(salts, "salts_lean")

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
write_csv_only(mixtures, "mixtures_lean")

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
write_csv_only(products, "products_lean")

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
write_csv_only(atc, "atc_lean")

# =============================================================================
# 9. LOOKUP TABLES - Canonical mappings for normalization
# =============================================================================
cat("[9] Lookup tables...\n")

# Salt suffixes (for stripping from generic names)
salt_suffixes <- c(
  "HYDROCHLORIDE", "HCL", "DIHYDROCHLORIDE",
  "SODIUM", "POTASSIUM", "CALCIUM", "MAGNESIUM", "ZINC", "IRON", "ALUMINUM", "ALUMINIUM",
  "SULFATE", "SULPHATE", "BISULFATE",
  "PHOSPHATE", "DIPHOSPHATE",
  "ACETATE", "DIACETATE",
  "CITRATE", "MALATE", "MALEATE", "FUMARATE", "TARTRATE", "SUCCINATE", "LACTATE", "GLUCONATE",
  "NITRATE", "CHLORIDE", "BROMIDE", "IODIDE", "FLUORIDE",
  "CARBONATE", "BICARBONATE",
  "OXIDE", "HYDROXIDE",
  "MESYLATE", "MESILATE", "BESYLATE", "BESILATE", "TOSYLATE",
  "STEARATE", "PALMITATE", "OLEATE",
  "PROPIONATE", "BUTYRATE", "VALERATE",
  "BENZOATE", "SALICYLATE",
  "TROMETHAMINE", "TROMETAMOL", "MEGLUMINE",
  "HYDROBROMIDE", "HYDROIODIDE",
  "MONOHYDRATE", "DIHYDRATE", "TRIHYDRATE", "HEMIHYDRATE", "ANHYDROUS"
)
salt_suffixes_df <- data.frame(
  salt_suffix = salt_suffixes,
  salt_suffix_key = tolower(salt_suffixes),
  stringsAsFactors = FALSE
)
write_csv_only(salt_suffixes_df, "lookup_salt_suffixes")
cat(sprintf("    salt_suffixes: %d\n", nrow(salt_suffixes_df)))

# Pure salt compounds (should NOT have salt stripped - the whole thing IS the salt)
pure_salt_compounds <- c(
  "SODIUM CHLORIDE", "POTASSIUM CHLORIDE", "CALCIUM CHLORIDE", "MAGNESIUM CHLORIDE",
  "SODIUM BICARBONATE", "POTASSIUM BICARBONATE", "CALCIUM CARBONATE", "MAGNESIUM CARBONATE",
  "SODIUM SULFATE", "POTASSIUM SULFATE", "MAGNESIUM SULFATE", "CALCIUM SULFATE",
  "SODIUM PHOSPHATE", "POTASSIUM PHOSPHATE", "CALCIUM PHOSPHATE",
  "SODIUM ACETATE", "POTASSIUM ACETATE", "CALCIUM ACETATE",
  "SODIUM CITRATE", "POTASSIUM CITRATE", "CALCIUM CITRATE", "MAGNESIUM CITRATE",
  "SODIUM LACTATE", "CALCIUM LACTATE",
  "SODIUM GLUCONATE", "CALCIUM GLUCONATE", "MAGNESIUM GLUCONATE", "POTASSIUM GLUCONATE",
  "FERROUS SULFATE", "FERROUS FUMARATE", "FERROUS GLUCONATE", "FERRIC SULFATE",
  "ZINC SULFATE", "ZINC OXIDE", "ZINC GLUCONATE",
  "CALCIUM HYDROXIDE", "MAGNESIUM HYDROXIDE", "ALUMINUM HYDROXIDE",
  "SODIUM HYDROXIDE", "POTASSIUM HYDROXIDE",
  "SODIUM FLUORIDE", "POTASSIUM FLUORIDE",
  "SODIUM IODIDE", "POTASSIUM IODIDE",
  "SODIUM BROMIDE", "POTASSIUM BROMIDE",
  "SODIUM NITRATE", "POTASSIUM NITRATE",
  "BARIUM SULFATE",
  "LITHIUM CARBONATE", "LITHIUM CITRATE"
)
pure_salts_df <- data.frame(
  compound = pure_salt_compounds,
  compound_key = tolower(pure_salt_compounds),
  stringsAsFactors = FALSE
)
write_csv_only(pure_salts_df, "lookup_pure_salts")
cat(sprintf("    pure_salts: %d\n", nrow(pure_salts_df)))

# Form canonical map (normalize dosage form names)
form_canonical <- data.frame(
  alias = c("tab", "tabs", "tablet", "tablets", "chewing gum",
            "cap", "caps", "capsule", "capsulee", "capsules",
            "susp", "suspension", "syr", "syrup",
            "sol", "soln", "solution", "inhal.solution", "instill.solution", "lamella",
            "ointment", "oint", "gel", "cream", "lotion",
            "patch", "supp", "suppository",
            "dpi", "inhal.powder", "mdi", "inhal.aerosol", "oral aerosol",
            "ampu", "ampul", "ampule", "ampoule", "amp", "vial",
            "inj", "injection", "implant", "s.c. implant",
            "metered dose inhaler", "dry powder inhaler",
            "spray", "nasal spray", "nebule", "neb", "inhaler"),
  canonical = c("TABLET", "TABLET", "TABLET", "TABLET", "TABLET",
                "CAPSULE", "CAPSULE", "CAPSULE", "CAPSULE", "CAPSULE",
                "SUSPENSION", "SUSPENSION", "SYRUP", "SYRUP",
                "SOLUTION", "SOLUTION", "SOLUTION", "SOLUTION", "SOLUTION", "SOLUTION",
                "OINTMENT", "OINTMENT", "GEL", "CREAM", "LOTION",
                "PATCH", "SUPPOSITORY", "SUPPOSITORY",
                "DPI", "DPI", "MDI", "MDI", "MDI",
                "AMPULE", "AMPULE", "AMPULE", "AMPULE", "AMPULE", "VIAL",
                "INJECTION", "INJECTION", "IMPLANT", "IMPLANT",
                "MDI", "DPI",
                "SPRAY", "SPRAY", "SOLUTION", "SOLUTION", "MDI"),
  stringsAsFactors = FALSE
)
write_csv_only(form_canonical, "lookup_form_canonical")
cat(sprintf("    form_canonical: %d\n", nrow(form_canonical)))

# Route alias map (normalize route names)
route_canonical <- data.frame(
  alias = c("oral", "po", "per orem", "per os", "by mouth",
            "iv", "intravenous",
            "im", "intramuscular",
            "sc", "subcut", "subcutaneous", "subdermal",
            "sl", "sublingual", "bucc", "buccal",
            "topical", "cutaneous",
            "dermal", "td", "transdermal",
            "oph", "eye", "ophthalmic",
            "otic", "ear",
            "inh", "neb", "inhalation", "inhaler",
            "rectal", "per rectum", "pr",
            "vaginal", "per vaginam", "pv",
            "intrathecal", "intranasal", "nasal", "per nasal",
            "intradermal", "id",
            "urethral", "intravesical", "endotracheal", "s.c. implant"),
  canonical = c("ORAL", "ORAL", "ORAL", "ORAL", "ORAL",
                "INTRAVENOUS", "INTRAVENOUS",
                "INTRAMUSCULAR", "INTRAMUSCULAR",
                "SUBCUTANEOUS", "SUBCUTANEOUS", "SUBCUTANEOUS", "SUBCUTANEOUS",
                "SUBLINGUAL", "SUBLINGUAL", "BUCCAL", "BUCCAL",
                "TOPICAL", "TOPICAL",
                "TRANSDERMAL", "TRANSDERMAL", "TRANSDERMAL",
                "OPHTHALMIC", "OPHTHALMIC", "OPHTHALMIC",
                "OTIC", "OTIC",
                "INHALATION", "INHALATION", "INHALATION", "INHALATION",
                "RECTAL", "RECTAL", "RECTAL",
                "VAGINAL", "VAGINAL", "VAGINAL",
                "INTRATHECAL", "NASAL", "NASAL", "NASAL",
                "INTRADERMAL", "INTRADERMAL",
                "URETHRAL", "INTRAVESICAL", "ENDOTRACHEAL", "SUBCUTANEOUS"),
  stringsAsFactors = FALSE
)
write_csv_only(route_canonical, "lookup_route_canonical")
cat(sprintf("    route_canonical: %d\n", nrow(route_canonical)))

# Form to route inference (infer route from form when not specified)
form_to_route <- data.frame(
  form = c("tablet", "capsule", "syrup", "suspension", "solution", "sachet",
           "drop", "eye drop", "ear drop",
           "cream", "ointment", "gel", "lotion",
           "patch",
           "inhaler", "nebule", "neb", "mdi", "dpi",
           "ampoule", "amp", "ampul", "ampule", "vial", "inj", "injection",
           "suppository", "supp",
           "spray", "nasal spray",
           "implant", "s.c. implant"),
  route = c("ORAL", "ORAL", "ORAL", "ORAL", "ORAL", "ORAL",
            "OPHTHALMIC", "OPHTHALMIC", "OTIC",
            "TOPICAL", "TOPICAL", "TOPICAL", "TOPICAL",
            "TRANSDERMAL",
            "INHALATION", "INHALATION", "INHALATION", "INHALATION", "INHALATION",
            "INTRAVENOUS", "INTRAVENOUS", "INTRAVENOUS", "INTRAVENOUS", "INTRAVENOUS", "INTRAVENOUS", "INTRAVENOUS",
            "RECTAL", "RECTAL",
            "NASAL", "NASAL",
            "SUBCUTANEOUS", "SUBCUTANEOUS"),
  stringsAsFactors = FALSE
)
write_csv_only(form_to_route, "lookup_form_to_route")
cat(sprintf("    form_to_route: %d\n", nrow(form_to_route)))

# Per-unit map (normalize dosage units)
per_unit <- data.frame(
  alias = c("ml", "l",
            "tab", "tabs", "tablet", "tablets", "chewing gum",
            "cap", "caps", "capsule", "capsules",
            "sachet", "sachets",
            "drop", "drops", "gtt",
            "actuation", "actuations",
            "spray", "sprays",
            "puff", "puffs",
            "dose", "doses",
            "application", "applications",
            "ampule", "ampules", "ampoule", "ampoules", "amp",
            "vial", "vials"),
  canonical = c("ML", "L",
                "TABLET", "TABLET", "TABLET", "TABLET", "TABLET",
                "CAPSULE", "CAPSULE", "CAPSULE", "CAPSULE",
                "SACHET", "SACHET",
                "DROP", "DROP", "DROP",
                "ACTUATION", "ACTUATION",
                "SPRAY", "SPRAY",
                "PUFF", "PUFF",
                "DOSE", "DOSE",
                "APPLICATION", "APPLICATION",
                "AMPULE", "AMPULE", "AMPULE", "AMPULE", "AMPULE",
                "VIAL", "VIAL"),
  stringsAsFactors = FALSE
)
write_csv_only(per_unit, "lookup_per_unit")
cat(sprintf("    per_unit: %d\n", nrow(per_unit)))

# =============================================================================
# Summary
# =============================================================================
cat("\n============================================================\n")
cat("Summary (vet-only excluded, normalized)\n")
cat("============================================================\n")
cat("Data tables:\n")
cat(sprintf("  generics_lean:  %6d (one per drug, with name_key)\n", nrow(generics)))
cat(sprintf("  synonyms_lean:  %6d (english, allowed coders, iupac ok if not alone)\n", nrow(synonyms)))
cat(sprintf("  dosages_lean:   %6d (valid form × route × strength)\n", nrow(dosages)))
cat(sprintf("  brands_lean:    %6d (with trademark symbols removed)\n", nrow(brands)))
cat(sprintf("  salts_lean:     %6d (with name_key)\n", nrow(salts)))
cat(sprintf("  mixtures_lean:  %6d (with split components, component_key_sorted)\n", nrow(mixtures)))
cat(sprintf("  products_lean:  %6d (name_type=generic/brand)\n", nrow(products)))
cat(sprintf("  atc_lean:       %6d (with hierarchy)\n", nrow(atc)))
cat("\nLookup tables:\n")
cat(sprintf("  lookup_salt_suffixes:   %3d (for stripping salts)\n", nrow(salt_suffixes_df)))
cat(sprintf("  lookup_pure_salts:      %3d (compounds that ARE salts)\n", nrow(pure_salts_df)))
cat(sprintf("  lookup_form_canonical:  %3d (form normalization)\n", nrow(form_canonical)))
cat(sprintf("  lookup_route_canonical: %3d (route normalization)\n", nrow(route_canonical)))
cat(sprintf("  lookup_form_to_route:   %3d (infer route from form)\n", nrow(form_to_route)))
cat(sprintf("  lookup_per_unit:        %3d (dose unit normalization)\n", nrow(per_unit)))
cat("============================================================\n")
cat("Files saved to:", output_dir, "\n")
