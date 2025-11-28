#!/usr/bin/env Rscript
# drugbank_brands.R — build a brands-focused DrugBank master dataset
# Extracts brand names from products (generic=false) and mixture names
# Output: drugbank_brands_master.csv with brand→generic mapping
# Can be run standalone: Rscript drugbank_brands.R

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

normalize_brand <- function(x) {
  val <- collapse_ws(x)
  if (is.na(val) || !nzchar(val)) return(NA_character_)
  # Uppercase and clean
  val <- toupper(val)
  # Remove trademark symbols
  val <- gsub("[®™©]", "", val, perl = TRUE)
  # Remove extra whitespace
  val <- trimws(gsub("\\s+", " ", val))
  if (!nzchar(val)) return(NA_character_)
  val
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

load_drugbank_dataset <- function() {
  if (exists("drugbank", inherits = FALSE)) {
    return(get("drugbank", inherits = FALSE))
  }
  if (exists("drugbank", where = "package:dbdataset")) {
    return(get("drugbank", envir = as.environment("package:dbdataset")))
  }
  data("drugbank", package = "dbdataset", envir = environment())
  if (!exists("drugbank", inherits = FALSE)) {
    stop("dbdataset::drugbank dataset is unavailable; reinstall dbdataset.")
  }
  get("drugbank", inherits = FALSE)
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("[drugbank_brands] Loading DrugBank dataset...\n")
dataset <- load_drugbank_dataset()

# ============================================================================
# PART 1: Extract brand names from products where generic = FALSE
# ============================================================================
cat("[drugbank_brands] Extracting brand names from products...\n")

products_dt <- as.data.table(dataset$products)
cat(sprintf("[drugbank_brands] Total products: %d\n", nrow(products_dt)))

# Filter to non-generic products (these are branded products)
# The 'generic' column indicates if the product is a generic version
brand_products <- products_dt[
  !is.na(generic) & tolower(generic) == "false" & !is.na(name) & nzchar(trimws(name))
]
cat(sprintf("[drugbank_brands] Brand products (generic=false): %d\n", nrow(brand_products)))

# Also get generic products for reference (products where generic=true)
generic_products <- products_dt[
  !is.na(generic) & tolower(generic) == "true" & !is.na(name) & nzchar(trimws(name))
]
cat(sprintf("[drugbank_brands] Generic products (generic=true): %d\n", nrow(generic_products)))

# Get general info to map drugbank_id to canonical generic name
general_dt <- as.data.table(dataset$drugs$general_information)[
  , .(drugbank_id = as.character(drugbank_id), canonical_generic_name = collapse_ws(name))
]

# Build brand→generic mapping from products
brand_from_products <- brand_products[
  , .(
    brand_name = collapse_ws(name),
    drugbank_id = as.character(drugbank_id),
    dosage_form = collapse_ws(dosage_form),
    strength = collapse_ws(strength),
    route = collapse_ws(route),
    country = collapse_ws(country),
    source = "products"
  )
]

# Normalize brand names
brand_from_products[, brand_name_normalized := sapply(brand_name, normalize_brand)]

# Join with general info to get generic name
brand_from_products <- merge(
  brand_from_products,
  general_dt,
  by = "drugbank_id",
  all.x = TRUE
)

# Remove rows where brand name equals generic name (not really a brand)
brand_from_products <- brand_from_products[
  is.na(canonical_generic_name) | 
  toupper(trimws(brand_name)) != toupper(trimws(canonical_generic_name))
]

cat(sprintf("[drugbank_brands] Brand entries from products: %d\n", nrow(brand_from_products)))

# ============================================================================
# PART 2: Extract brand names from mixtures dataset
# ============================================================================
cat("[drugbank_brands] Extracting brand names from mixtures...\n")

mixtures_path <- file.path(output_dir, "drugbank_mixtures_master.csv")
if (file.exists(mixtures_path)) {
  mixtures_dt <- fread(mixtures_path)
  
  # Get unique mixture_drugbank_id to ingredient_components mapping
  # (avoid cartesian join by deduplicating first)
  mixture_generics <- unique(mixtures_dt[
    !is.na(mixture_drugbank_id) & !is.na(ingredient_components),
    .(mixture_drugbank_id, ingredient_components)
  ], by = "mixture_drugbank_id")
  
  # mixture_name is often a brand name for the combination
  # Get unique brand names per mixture
  mixture_brands <- unique(mixtures_dt[
    !is.na(mixture_name) & nzchar(trimws(mixture_name)),
    .(
      brand_name = collapse_ws(mixture_name),
      drugbank_id = as.character(mixture_drugbank_id)
    )
  ], by = c("brand_name", "drugbank_id"))
  
  mixture_brands[, `:=`(
    dosage_form = NA_character_,
    strength = NA_character_,
    route = NA_character_,
    country = NA_character_,
    source = "mixtures"
  )]
  
  # Normalize brand names
  mixture_brands[, brand_name_normalized := sapply(brand_name, normalize_brand)]
  
  # Get generic name (ingredients) from mixtures - use deduplicated lookup
  mixture_brands <- merge(
    mixture_brands,
    mixture_generics,
    by.x = "drugbank_id",
    by.y = "mixture_drugbank_id",
    all.x = TRUE
  )
  setnames(mixture_brands, "ingredient_components", "canonical_generic_name")
  
  cat(sprintf("[drugbank_brands] Brand entries from mixtures: %d\n", nrow(mixture_brands)))
} else {
  cat("[drugbank_brands] Warning: mixtures master not found, skipping mixture brands\n")
  mixture_brands <- data.table(
    brand_name = character(),
    drugbank_id = character(),
    dosage_form = character(),
    strength = character(),
    route = character(),
    country = character(),
    source = character(),
    brand_name_normalized = character(),
    canonical_generic_name = character()
  )
}

# ============================================================================
# PART 3: Combine and deduplicate
# ============================================================================
cat("[drugbank_brands] Combining brand sources...\n")

# Ensure columns match
brand_from_products[, is_mixture := FALSE]
mixture_brands[, is_mixture := TRUE]

all_brands <- rbindlist(
  list(brand_from_products, mixture_brands),
  use.names = TRUE,
  fill = TRUE
)

# Remove empty/NA brand names
all_brands <- all_brands[!is.na(brand_name_normalized) & nzchar(brand_name_normalized)]

# Deduplicate by normalized brand name + drugbank_id
all_brands <- unique(all_brands, by = c("brand_name_normalized", "drugbank_id"))

cat(sprintf("[drugbank_brands] Total unique brand entries: %d\n", nrow(all_brands)))

# ============================================================================
# PART 4: Create summary statistics
# ============================================================================
unique_brands <- length(unique(all_brands$brand_name_normalized))
unique_generics <- length(unique(all_brands$canonical_generic_name[!is.na(all_brands$canonical_generic_name)]))
from_products <- sum(all_brands$source == "products")
from_mixtures <- sum(all_brands$source == "mixtures")

cat(sprintf("[drugbank_brands] Unique brand names: %d\n", unique_brands))
cat(sprintf("[drugbank_brands] Unique generic names: %d\n", unique_generics))
cat(sprintf("[drugbank_brands] From products: %d\n", from_products))
cat(sprintf("[drugbank_brands] From mixtures: %d\n", from_mixtures))

# ============================================================================
# PART 5: Write output
# ============================================================================
output_path <- file.path(output_dir, "drugbank_brands_master.csv")

# Select and order columns for output
output_cols <- c(
  "brand_name",
  "brand_name_normalized", 
  "drugbank_id",
  "canonical_generic_name",
  "is_mixture",
  "dosage_form",
  "strength",
  "route",
  "country",
  "source"
)

# Ensure all columns exist
for (col in output_cols) {
  if (!col %in% names(all_brands)) {
    all_brands[, (col) := NA_character_]
  }
}

fwrite(all_brands[, ..output_cols], output_path)
cat(sprintf("[drugbank_brands] Output written to %s\n", output_path))

# ============================================================================
# PART 6: Also export products info for enriching generics
# ============================================================================
cat("[drugbank_brands] Exporting products info for generics enrichment...\n")

# Export all products (both generic and brand) with dose/form/route info
# This can be used to enrich the generics master with real-world variants
products_export <- products_dt[
  !is.na(drugbank_id),
  .(
    drugbank_id = as.character(drugbank_id),
    product_name = collapse_ws(name),
    is_generic = tolower(generic) == "true",
    dosage_form = collapse_ws(dosage_form),
    strength = collapse_ws(strength),
    route = collapse_ws(route),
    country = collapse_ws(country),
    approved = collapse_ws(approved),
    over_the_counter = collapse_ws(over_the_counter)
  )
]

products_output_path <- file.path(output_dir, "drugbank_products_export.csv")
fwrite(products_export, products_output_path)
cat(sprintf("[drugbank_brands] Products export written to %s (%d rows)\n", 
            products_output_path, nrow(products_export)))

cat("[drugbank_brands] Done.\n")
