#!/usr/bin/env Rscript
# drugbank_salts.R — extract salt forms from DrugBank
# Output: drugbank_salts_master.csv with salt names for detection/stripping
# Used to normalize drug names by removing salt suffixes
# Can be run standalone: Rscript drugbank_salts.R

# ============================================================================
# SHARED SETUP - Source _shared.R if not already loaded
# ============================================================================
if (!exists("DRUGBANK_SHARED_LOADED") || !isTRUE(DRUGBANK_SHARED_LOADED)) {
  script_dir <- if (exists("DRUGBANK_BASE_DIR")) {
    DRUGBANK_BASE_DIR
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      getwd()
    }
  }
  shared_path <- file.path(script_dir, "_shared.R")
  if (file.exists(shared_path)) {
    source(shared_path, local = FALSE)
  } else {
    stop("_shared.R not found. Run from drugbank_generics directory or via drugbank_all_v2.R")
  }
}

collapse_ws <- function(x) {
  ifelse(is.na(x), NA_character_, trimws(gsub("\\s+", " ", as.character(x))))
}

normalize_salt <- function(x) {
  val <- collapse_ws(x)
  if (is.na(val) || !nzchar(val)) return(NA_character_)
  # Uppercase and clean
  val <- toupper(val)
  # Remove extra whitespace
  val <- trimws(gsub("\\s+", " ", val))
  return(val)
}

# ============================================================================
# MAIN
# ============================================================================

if (!quiet_mode) cat("[drugbank_salts] Loading DrugBank salts data...\n")

# Get salts dataset from DrugBank
# The salts dataset contains salt forms of drugs
salts_raw <- tryCatch({
  as.data.table(drugbank$salts)
}, error = function(e) {
  if (!quiet_mode) cat("[drugbank_salts] Warning: Could not load drugbank$salts:", conditionMessage(e), "\n")
  data.table()
})

if (!quiet_mode) cat("[drugbank_salts] Raw salts rows:", nrow(salts_raw), "\n")

if (nrow(salts_raw) > 0) {
  # Inspect columns
  if (!quiet_mode) {
    cat("[drugbank_salts] Columns:", paste(names(salts_raw), collapse = ", "), "\n")
  }
  
  # Expected columns: drugbank_id, name, unii, cas_number, inchikey, average_mass, monoisotopic_mass
  # We need: drugbank_id (parent drug), name (salt name)
  
  salts_dt <- salts_raw[, .(
    parent_drugbank_id = if ("drugbank_id" %in% names(salts_raw)) drugbank_id else NA_character_,
    salt_name = if ("name" %in% names(salts_raw)) name else NA_character_,
    salt_unii = if ("unii" %in% names(salts_raw)) unii else NA_character_,
    salt_cas = if ("cas_number" %in% names(salts_raw)) cas_number else NA_character_,
    salt_inchikey = if ("inchikey" %in% names(salts_raw)) inchikey else NA_character_
  )]
  
  # Normalize salt names
  salts_dt[, salt_name_normalized := sapply(salt_name, normalize_salt)]
  
  # Remove rows with no salt name
  salts_dt <- salts_dt[!is.na(salt_name_normalized) & nzchar(salt_name_normalized)]
  
  if (!quiet_mode) cat("[drugbank_salts] Valid salt entries:", nrow(salts_dt), "\n")
  
  # Extract unique salt suffixes for matching
  # Common patterns: "HYDROCHLORIDE", "SODIUM", "POTASSIUM", "SULFATE", etc.
  # We'll extract the last word(s) that are common salt indicators
  
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
  
  # Create a lookup table of known salt suffixes
  salt_suffix_dt <- data.table(
    salt_suffix = salt_suffixes,
    salt_suffix_normalized = toupper(salt_suffixes)
  )
  
  # Also extract unique salt names from the dataset
  unique_salts <- unique(salts_dt[, .(salt_name_normalized)])
  unique_salts[, is_from_drugbank := TRUE]
  
  if (!quiet_mode) cat("[drugbank_salts] Unique salt names from DrugBank:", nrow(unique_salts), "\n")
  
} else {
  # Fallback: create from known salt suffixes
  if (!quiet_mode) cat("[drugbank_salts] No salts data found, using hardcoded list\n")
  
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
  
  salts_dt <- data.table(
    parent_drugbank_id = NA_character_,
    salt_name = salt_suffixes,
    salt_name_normalized = toupper(salt_suffixes),
    salt_unii = NA_character_,
    salt_cas = NA_character_,
    salt_inchikey = NA_character_
  )
  
  unique_salts <- data.table(
    salt_name_normalized = toupper(salt_suffixes),
    is_from_drugbank = FALSE
  )
}

# Also identify "pure salt" compounds - these should NOT have salt stripped
# These are compounds where the entire name IS the salt (e.g., SODIUM CHLORIDE, CALCIUM CARBONATE)
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

pure_salts_dt <- data.table(
  compound_name = pure_salt_compounds,
  compound_name_normalized = toupper(pure_salt_compounds),
  is_pure_salt = TRUE
)

if (!quiet_mode) cat("[drugbank_salts] Pure salt compounds:", nrow(pure_salts_dt), "\n")

# ============================================================================
# OUTPUT
# ============================================================================

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
}

output_dir <- file.path(get_script_dir(), "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Write salts master (CSV + Parquet)
salts_out <- file.path(output_dir, "drugbank_salts_master.csv")
write_csv_and_parquet(salts_dt, salts_out)
if (!quiet_mode) cat("[drugbank_salts] Wrote:", salts_out, "(+parquet,", nrow(salts_dt), "rows)\n")

# Write pure salts reference (CSV + Parquet)
pure_salts_out <- file.path(output_dir, "drugbank_pure_salts.csv")
write_csv_and_parquet(pure_salts_dt, pure_salts_out)
if (!quiet_mode) cat("[drugbank_salts] Wrote:", pure_salts_out, "(+parquet,", nrow(pure_salts_dt), "rows)\n")

# Write unique salt suffixes for quick lookup (CSV + Parquet)
salt_suffixes_dt <- data.table(
  salt_suffix = salt_suffixes,
  salt_suffix_normalized = toupper(salt_suffixes)
)
suffixes_out <- file.path(output_dir, "drugbank_salt_suffixes.csv")
write_csv_and_parquet(salt_suffixes_dt, suffixes_out)
if (!quiet_mode) cat("[drugbank_salts] Wrote:", suffixes_out, "(+parquet,", nrow(salt_suffixes_dt), "rows)\n")

if (!quiet_mode) cat("[drugbank_salts] Done.\n")
