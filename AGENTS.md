# AGENT INSTRUCTIONS

Apply these rules when editing `dependencies/drugbank_generics`:

1. **r-polars first.** Default to `polars` for all tabular work. Avoid introducing `data.table`/tidyverse/arrow flows; if a temporary fallback is unavoidable, add a `TODO(polars): convert` next to it.
2. **Parquet as the canonical artifact.** Write Parquet outputs as the primary deliverable and emit CSV only as a debugging convenience. Do not add CSV-based reads in this folder.
3. **Keep scripts standalone.** `drugbank_generics.R`, `drugbank_mixtures.R`, and helpers must run from this directory with outputs staying under `output/`.
4. **Stay in sync.** When schema or column expectations change, propagate updates across the generics and mixtures scripts so downstream consumers (e.g., `drugbank_all.R`) keep working.
5. **Minimize dependencies.** Prefer base R helpers and `polars`; avoid adding new non-Polars tabular dependencies unless essential.
