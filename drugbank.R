# =================================
# DrugBank → aggregated ATC-level generics + brand mapping for ESOA ingestion
# Outputs:
# - drugbank_generics.csv: generic, dose, form, route, atc_code, drugbank_ids
# - drugbank_brands.csv: brand, generic
# =================================

ensure_installed <- function(packages, repos = "https://cloud.r-project.org") {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = repos)
    }
  }
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
ensure_installed(c("data.table", "arrow"))
if (!requireNamespace("dbdataset", quietly = TRUE)) {
  ensure_installed("remotes")
  remotes::install_github("interstellar-Consultation-Services/dbdataset", quiet = TRUE, upgrade = "never")
}

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(dbdataset)
})

collapse_ws <- function(x) {
  x <- gsub("\\s+", " ", x, perl = TRUE)
  trimws(x)
}

empty_to_na <- function(x) {
  x <- trimws(x)
  x[!nzchar(x)] <- NA_character_
  x
}

unique_canonical <- function(x) {
  if (!length(x)) return(character())
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  if (!length(x)) return(character())
  ord <- order(tolower(x), x)
  x <- x[ord]
  x[!duplicated(tolower(x))]
}

combine_values <- function(...) {
  vals <- unlist(list(...), use.names = FALSE)
  unique_canonical(vals)
}

write_arrow_csv <- function(dt, path) {
  cat(sprintf("Writing %s rows to %s\n", format(nrow(dt), big.mark = ","), path))
  arrow_table <- arrow::Table$create(dt)
  if ("write_csv_arrow" %in% getNamespaceExports("arrow")) {
    arrow::write_csv_arrow(arrow_table, path)
  } else if ("write_delim_arrow" %in% getNamespaceExports("arrow")) {
    arrow::write_delim_arrow(arrow_table, path, delim = ",")
  } else {
    warning("Falling back to data.table::fwrite because Arrow CSV writers are unavailable.")
    fwrite(dt, path)
  }
  cat("Done.\n")
}

split_ingredients <- function(value) {
  if (is.null(value) || length(value) == 0L) return(character())
  parts <- unlist(lapply(value, function(v) {
    if (is.na(v) || !nzchar(v)) return(character())
    v <- collapse_ws(v)
    if (!nzchar(v)) return(character())
    v <- gsub("\\s+(and|with|plus)\\s+", "|", v, ignore.case = TRUE, perl = TRUE)
    v <- gsub("\\s*/\\s*|\\s*\\+\\s*|\\s*,\\s*|\\s*;\\s*", "|", v, perl = TRUE)
    res <- strsplit(v, "\\|", perl = TRUE)[[1]]
    res <- trimws(res)
    res[nzchar(res)]
  }), use.names = FALSE)
  unique_canonical(parts)
}

split_semicolon <- function(value) {
  if (is.null(value) || length(value) == 0L) return(character())
  parts <- unlist(lapply(value, function(v) {
    if (is.na(v) || !nzchar(v)) return(character())
    res <- strsplit(v, ";", fixed = TRUE)[[1]]
    res <- trimws(res)
    res[nzchar(res)]
  }), use.names = FALSE)
  unique_canonical(parts)
}

normalize_list_column <- function(col) {
  lapply(col, function(item) {
    if (is.null(item)) return(character())
    if (length(item) == 1L && is.na(item)) return(character())
    item <- unlist(item, use.names = FALSE)
    if (!length(item)) return(character())
    item <- trimws(as.character(item))
    item <- item[nzchar(item)]
    unique_canonical(item)
  })
}

ensure_list_column <- function(DT, col) {
  if (!col %in% names(DT)) {
    DT[, (col) := vector("list", .N)]
  }
  DT[, (col) := normalize_list_column(.SD[[1]]), .SDcols = col]
}

combine_list_column <- function(lst) {
  if (is.null(lst) || !length(lst)) return(character())
  combine_values(unlist(lst, use.names = FALSE))
}

collapse_pipe <- function(values) {
  values <- combine_values(values)
  if (!length(values)) {
    return(NA_character_)
  }
  paste(values, collapse = "|")
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_generics_path <- file.path(output_dir, "drugbank_generics.csv")
output_brands_path <- file.path(output_dir, "drugbank_brands.csv")

dataset <- drugbank

excluded_groups <- c("experimental", "withdrawn", "illicit", "vet")
groups_dt <- as.data.table(dataset$drugs$groups)
groups_dt[, group_clean := tolower(trimws(group))]
excluded_ids <- unique(groups_dt[group_clean %chin% excluded_groups, as.character(drugbank_id)])

general_dt <- as.data.table(dataset$drugs$general_information)[
  , .(
    drugbank_id = as.character(drugbank_id),
    generic_name = collapse_ws(name)
  )
]
general_dt <- general_dt[!is.na(generic_name) & nzchar(generic_name)]
if (length(excluded_ids)) general_dt <- general_dt[!(drugbank_id %chin% excluded_ids)]

generic_dt <- copy(general_dt)
generic_dt[, generic_parts := lapply(generic_name, split_ingredients)]

ensure_list_column(generic_dt, "generic_parts")
generic_dt[, generic := vapply(generic_parts, function(parts) {
  if (!length(parts)) return(NA_character_)
  paste(parts, collapse = "; ")
}, character(1))]

products_dt <- as.data.table(dataset$products)[
  , .(
    drugbank_id = as.character(drugbank_id),
    product_name = collapse_ws(name),
    route = collapse_ws(route),
    dosage_form = collapse_ws(dosage_form),
    strength = collapse_ws(strength)
  )
]
if (length(excluded_ids)) products_dt <- products_dt[!(drugbank_id %chin% excluded_ids)]
products_dt[, c("product_name", "route", "dosage_form", "strength") := lapply(
  .SD, empty_to_na
), .SDcols = c("product_name", "route", "dosage_form", "strength")]

product_agg <- products_dt[
  ,
  {
    names_vec <- unique_canonical(product_name[!is.na(product_name)])
    route_vec <- unique_canonical(unlist(lapply(route[!is.na(route)], split_semicolon), use.names = FALSE))
    form_vec <- unique_canonical(unlist(lapply(dosage_form[!is.na(dosage_form)], split_semicolon), use.names = FALSE))
    dose_vec <- unique_canonical(strength[!is.na(strength)])
    list(
      product_names = list(names_vec),
      product_routes = list(route_vec),
      product_forms = list(form_vec),
      product_doses = list(dose_vec)
    )
  },
  by = drugbank_id
]

dosages_dt <- as.data.table(dataset$drugs$dosages)[
  , .(
    drugbank_id = as.character(drugbank_id),
    route = collapse_ws(route),
    dosage_form = collapse_ws(form),
    strength = collapse_ws(strength)
  )
]
if (length(excluded_ids)) dosages_dt <- dosages_dt[!(drugbank_id %chin% excluded_ids)]
dosages_dt[, c("route", "dosage_form", "strength") := lapply(
  .SD, empty_to_na
), .SDcols = c("route", "dosage_form", "strength")]

dosage_agg <- dosages_dt[
  ,
  {
    route_vec <- unique_canonical(unlist(lapply(route[!is.na(route)], split_semicolon), use.names = FALSE))
    form_vec <- unique_canonical(unlist(lapply(dosage_form[!is.na(dosage_form)], split_semicolon), use.names = FALSE))
    dose_vec <- unique_canonical(strength[!is.na(strength)])
    list(
      dosage_routes = list(route_vec),
      dosage_forms = list(form_vec),
      dosage_doses = list(dose_vec)
    )
  },
  by = drugbank_id
]

brands_dt <- as.data.table(dataset$drugs$international_brands)[
  , .(
    drugbank_id = as.character(drugbank_id),
    brand_name = collapse_ws(brand)
  )
]
if (length(excluded_ids)) brands_dt <- brands_dt[!(drugbank_id %chin% excluded_ids)]
brands_dt[, brand_name := empty_to_na(brand_name)]

brand_agg <- brands_dt[
  !is.na(brand_name),
  .(brand_names = list(unique_canonical(brand_name))),
  by = drugbank_id
]

atc_dt <- as.data.table(dataset$drugs$atc_codes)[
  , .(
    drugbank_id = as.character(drugbank_id),
    atc_code = trimws(atc_code)
  )
]
atc_dt <- atc_dt[!is.na(atc_code) & nzchar(atc_code)]
if (length(excluded_ids)) atc_dt <- atc_dt[!(drugbank_id %chin% excluded_ids)]
atc_dt <- unique(atc_dt, by = c("drugbank_id", "atc_code"))

drug_ids <- unique(atc_dt$drugbank_id)
info_dt <- data.table(drugbank_id = drug_ids)

info_dt <- merge(info_dt, generic_dt[, .(drugbank_id, generic, generic_parts)], by = "drugbank_id", all.x = TRUE, sort = FALSE)
info_dt <- merge(info_dt, product_agg, by = "drugbank_id", all.x = TRUE, sort = FALSE)
info_dt <- merge(info_dt, dosage_agg, by = "drugbank_id", all.x = TRUE, sort = FALSE)
info_dt <- merge(info_dt, brand_agg, by = "drugbank_id", all.x = TRUE, sort = FALSE)

ensure_list_column(info_dt, "generic_parts")
ensure_list_column(info_dt, "product_names")
ensure_list_column(info_dt, "product_routes")
ensure_list_column(info_dt, "product_forms")
ensure_list_column(info_dt, "product_doses")
ensure_list_column(info_dt, "dosage_routes")
ensure_list_column(info_dt, "dosage_forms")
ensure_list_column(info_dt, "dosage_doses")
ensure_list_column(info_dt, "brand_names")

brands_output <- info_dt[
  ,
  {
    brand_vec <- brand_names[[1]]
    if (!length(brand_vec)) {
      NULL
    } else {
      generic_vec <- unique_canonical(generic_parts[[1]])
      if (!length(generic_vec)) {
        generic_vec <- combine_values(generic)
      }
      generic_str <- if (length(generic_vec)) paste(generic_vec, collapse = "; ") else NA_character_
      data.table(
        brand = brand_vec,
        generic = generic_str
      )
    }
  },
  by = drugbank_id
]
if (!is.null(brands_output) && nrow(brands_output)) {
  brands_output <- unique(brands_output, by = c("brand", "generic"))
  setorder(brands_output, brand, generic)
} else {
  brands_output <- data.table(brand = character(), generic = character())
}

info_dt[, route_values := Map(function(prod_route, dosage_route) {
  combine_values(prod_route, dosage_route)
}, product_routes, dosage_routes)]

info_dt[, form_values := Map(function(prod_form, dosage_form) {
  combine_values(prod_form, dosage_form)
}, product_forms, dosage_forms)]

info_dt[, dose_values := Map(function(prod_dose, dosage_dose) {
  combine_values(prod_dose, dosage_dose)
}, product_doses, dosage_doses)]

info_dt[, c(
  "product_names", "brand_names",
  "product_routes", "dosage_routes",
  "product_forms", "dosage_forms",
  "product_doses", "dosage_doses"
) := NULL]

drug_atc <- merge(atc_dt, info_dt, by = "drugbank_id", all.x = TRUE, sort = FALSE)

generics_dt <- drug_atc[
  ,
  {
    generic_vec <- combine_list_column(generic_parts)
    if (!length(generic_vec)) {
      generic_vec <- combine_values(generic)
    }
    generic_str <- if (length(generic_vec)) paste(generic_vec, collapse = "; ") else NA_character_

    doses <- combine_list_column(dose_values)
    forms <- combine_list_column(form_values)
    routes <- combine_list_column(route_values)
    ids <- combine_values(drugbank_id)

    list(
      generic = generic_str,
      dose = list(doses),
      form = list(forms),
      route = list(routes),
      drugbank_ids = list(ids)
    )
  },
  by = atc_code
]

generics_dt[, dose := vapply(dose, collapse_pipe, character(1))]
generics_dt[, form := vapply(form, collapse_pipe, character(1))]
generics_dt[, route := vapply(route, collapse_pipe, character(1))]
generics_dt[, drugbank_ids := vapply(drugbank_ids, collapse_pipe, character(1))]

setcolorder(generics_dt, c("generic", "dose", "form", "route", "atc_code", "drugbank_ids"))
setorder(generics_dt, atc_code)

if (nrow(generics_dt) > 100000) {
  message(sprintf("Generics output currently has %s rows (>100k)", format(nrow(generics_dt), big.mark = ",")))
}

write_arrow_csv(generics_dt, output_generics_path)
write_arrow_csv(brands_output, output_brands_path)
