# DrugBank lean exports (PIDS DRG)

Generates lean, normalized CSV extracts from the DrugBank dataset for downstream PIDS DRG matching pipelines (FDA scrapers, authority builders, etc.). Outputs live under `./output/` and are meant to be copied into consumer repos (e.g., `../pids-drg-fda-scraper/input/`).

## Requirements
- R ≥ 4.1
- `dbdataset` package available in the R library (provides the DrugBank object). Install from its GitHub source if not already present.

## Quick start
```bash
Rscript drugbank_lean_export.R
```
Creates `./output` and writes lean tables (CSV only). No arguments are needed; rerun whenever the upstream DrugBank dataset updates.

## Outputs
- `generics_lean.csv` – one row per drugbank_id with canonical generic name and key
- `synonyms_lean.csv` – English synonyms with allowed coders (INN/USAN/BAN/JAN/etc.) and keys
- `dosages_lean.csv` – form × route × strength combos
- `brands_lean.csv` – international brands with trademark symbols removed
- `salts_lean.csv` – salt variants with normalized keys
- `mixtures_lean.csv` – mixtures with split component lists and keys
- `products_lean.csv` – product rows with name_type flag (generic/brand)
- `atc_lean.csv` – ATC hierarchy per drug
- Lookup tables: `lookup_salt_suffixes`, `lookup_pure_salts`, `lookup_form_canonical`, `lookup_route_canonical`, `lookup_form_to_route`, `lookup_per_unit`

## Behavior
- Filters out vet-only drugs (keeps human-approved entries).
- Normalizes names, removes trademark symbols, builds deterministic keys, and dedupes rows.
- Collapses whitespace, strips salts where appropriate, and produces sorted, reproducible outputs.

## Downstream use
Copy `generics_lean.csv` and `synonyms_lean.csv` into `../pids-drg-fda-scraper/input/` to improve brand/generic detection, and use the other tables in matching/authority pipelines as needed.
