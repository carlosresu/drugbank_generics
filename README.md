# pids-drg-drugbank-generics

Generates lean, normalized CSV extracts from the DrugBank dataset for downstream PIDS DRG matching pipelines (FDA scrapers, authority builders, etc.). Outputs live under `./output/` and are meant to be copied into consumer repos (e.g., `../pids-drg-fda-scraper/input/`).

**Last Updated:** December 2025

---

## Progress Since June 2025

| Date | Milestone |
|------|-----------|
| **Jun-Jul 2025** | Initial R script with basic generics export |
| **Aug 2025** | Added synonyms with language/coder filtering |
| **Sep 2025** | Added dosages, brands, salts, mixtures lean tables |
| **Oct 2025** | Added lookup tables for normalization (salt suffixes, form/route canonical) |
| **Nov 2025** | Created `_shared.R` for common setup, optimized parallelization |
| **Dec 2025** | Default workers set to min(8, cores), cross-platform support |

### Current Outputs

- **8 Data Tables:** generics, synonyms, dosages, brands, salts, mixtures, products, atc_lean
- **6 Lookup Tables:** salt_suffixes, pure_salts, form_canonical, route_canonical, form_to_route, per_unit
- **Runtime:** ~433s total with 8 workers

---

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
