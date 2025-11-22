# DrugBank Pipeline (r-polars)

This repository builds DrugBank exports with `polars` (r-polars). Parquet is the canonical output; CSVs are sidecars for inspection.

## Entrypoints & configuration
- `drugbank_all.R` runs both exporters (`drugbank_generics.R` then `drugbank_mixtures.R`) from this directory. It forwards `ESOA_DRUGBANK_QUIET` and respects `ESOA_DRUGBANK_WORKERS`/`POLARS_MAX_THREADS` for threading.
- `drugbank_generics.R` accepts `--keep-all` to skip the default “approved-or-ATC” filter.
- Outputs live under `output/` beside the scripts:
  - `drugbank_generics_master.parquet` (+ CSV sidecar)
  - `drugbank_mixtures_master.parquet` (+ CSV sidecar)

## Common setup
- Guard rails install `polars (>= 1.6.0)` and `dbdataset` if missing.
- `dbdataset` exposes the in-memory `drugbank` object; all wrangling is done with `polars::pl` DataFrames and `map_batches`/lazy-style expressions.
- Vet-only drugs are excluded early so they never re-enter downstream joins.

## Generics exporter (`drugbank_generics.R`)
1. Build base tables from `drugbank`:
   - Generic names (`general_information`), English synonyms filtered to allowed coders (inn/usan/ban/jan/dcj/usp/dcit), ATC codes, group memberships, and salts.
2. Normalize dosage/product routes, forms, and strengths:
   - Collapse whitespace, strip parenthetical/commas into base + detail fields, normalize dose units, and canonicalize route/form values.
   - Deduplicate and expand route/form cartesian pairs per drug, carrying raw/original text for traceability.
3. Join enriched attributes:
   - Attach ATC codes, salts, group flags, canonical names + synonyms (lexeme list), and derive synonym expansions for dose/form/route.
4. Filter and sort:
   - By default keep drugs that are approved or have any ATC code (use `--keep-all` to disable).
   - Sort by canonical name/lexeme/ATC/dose/form/route, then write Parquet+CSV.
5. Emit diagnostics:
   - Writes `drugbank_route_form_exclusions.*` capturing disallowed route/form combinations spotted during processing.

## Mixtures exporter (`drugbank_mixtures.R`)
1. Locate the generics Parquet (prefers `output/drugbank_generics_master.parquet` in this project) and build lookups:
   - `lexeme -> drugbank_id` and `drugbank_id -> generic_components_key`.
2. Prepare group/salt lookup tables (with vet exclusions) from `drugbank`.
3. Parse mixtures from `dataset$drugs$mixtures`:
   - Normalize names and ingredient strings, split components, derive lexeme keys, and map to component DrugBank IDs + generic component keys.
4. Attach group/salt names to each mixture, order rows deterministically, and write Parquet+CSV.

## Run notes
- Scripts are designed to run from the `dependencies/drugbank_generics` root; outputs stay under `output/`.
- The legacy `drugbank.R` simply delegates to `drugbank_all.R` and is kept for compatibility with older callers.
