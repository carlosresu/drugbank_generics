# DrugBank Authority Pipeline Walkthrough

This document breaks down `drugbank.R` into the exact sequence of actions it performs to build `output/drugbank_authority.csv`.

## 1. Environment Setup
1. **Package guard rails** – `ensure_installed()` installs `data.table`, `stringi`, `progressr`, and `devtools` if they are missing. If the custom `dbdataset` package is unavailable, it is installed from GitHub.
2. **Library loading** – the script suppresses startup noise and attaches `data.table`, `stringi`, `progressr`, and `dbdataset`. All data wrangling uses `data.table`.
3. **Threading & progress handlers** – `setDTthreads(percent = 100)` allows `data.table` to use all available cores. `progressr::handlers(global = TRUE)` and `handlers("cli")` switch on CLI-style progress bars.
4. **Paths** – `get_script_dir()` resolves the script’s own directory (supports interactive/source usage). The script ensures `output/` exists beside the script and defines `output/drugbank_authority.csv` as the target.

## 2. Helper Functions
5. **Whitespace normalization** – `collapse_ws()` squashes repeated whitespace and trims the ends; `dedupe_key()` lowercases these results to build dedupe keys.
6. **Dosage form simplification** – `simplify_form_one()` keeps the leading portion of a dosage form value and, if present, the first “*release*” qualifier (e.g., “extended release”). `simplify_form_vec()` maps the helper across a vector.
7. **Route parsing** – `normalize_route_field()` lowercases trimmed routes; `split_semicolon()` splits `;`-separated routes/forms into normalized vectors.
8. **Mixture ingredient parsing** – `split_ingredients()` breaks ingredient text on connectors (`and`, `/`, `+`, commas, semicolons) while preserving multi-word names.
9. **Brand extraction** – `extract_product_brand()` tries to pull a known brand from free-text product names, using case-insensitive pattern matches and falling back to a heuristic first word.
10. **Dosage form column simplifier** – `simplify_form_column()` processes massive columns in ~50k-row chunks, using `progressr` to report progress and occasionally calling `gc()` to keep memory down.
11. **Route/form expansion** – `expand_route_form()` splits `route` and `dosage_form` strings into vectors, cartesian-joins them with `CJ()`, recombines the rest of the columns, normalizes, and dedupes the result.

## 3. Load & Filter Base Dataset
12. **Dataset import** – `dataset <- drugbank` pulls the in-memory object exposed by `dbdataset`.
13. **Drug exclusions** – the script builds `bad_ids` for any drug whose groups include “experimental”, “withdrawn”, “illicit”, or “vet” (case-insensitive) and drops them from all downstream tables.

## 4. ATC Mapping and Combo Flag
14. **ATC table** – reads `dataset$drugs$atc_codes`, keeping non-empty codes. Each row contributes `drugbank_id` and its `atc_code`.
15. **Combination detection** – the last two digits of each ATC code (parsed into `suffix`) mark combination therapies if `suffix >= 50` or the code begins with `J05AR`. The flag is collapsed per drug in `combo_flag`.
16. **Cleanup** – temporary columns are removed, duplicates collapsed, and both `atc_tbl` and `combo_flag` are filtered against `bad_ids` and keyed by `drugbank_id`.

## 5. Derive Name/Text Sources
17. **Primary generic names** – from `dataset$drugs$general_information`, extract `(drugbank_id, generic_name)` with whitespace collapsed.
18. **Generic text rows** – duplicate the same table to capture the canonical name text, dedupe by a case-folded key so variations map to one entry, and tag rows as `text_type = "generic"`.
19. **Synonyms** – pull English synonyms whose coder mentions `inn`, `ban`, `usan`, `jan`, or `usp`; strip whitespace, dedupe against the generic names, and tag `text_type = "synonym"` with the lowercase coder as `text_subtype`.
20. **International brands** – gather brand names and companies, dedupe case-insensitively, and tag as `text_type = "brand"`, `text_subtype = "international"`. A `brand_lookup` list per drug is built here for later brand extraction.
21. **Mixtures** – create two sets: one for mixture names (`text_subtype = "mixture_name"`) and one for extracted mixture ingredients (`text_subtype = "mixture_ingredient"`), skipping empty entries.

## 6. Route / Form / Strength Inputs
22. **Dosage records** – from `dataset$drugs$dosages` build `(drugbank_id, route, dosage_form, strength, source="DOSAGES", data_origin="dosage")`, normalize fields, and simplify dosage forms in chunks.
23. **Product records** – from `dataset$products`, collect `(drugbank_id, product_name, route, dosage_form, strength, source)`; normalize blank strings to `NA`, attach the `brand_lookup`, infer a `brand` name via `extract_product_brand()`, set `data_origin="product"`, and simplify dosage forms in chunks.
24. **Unified route/form table** – stack dosage and product records into `route_form_all`, keeping unique combinations.
25. **Product brand admin rows** – capture `(drugbank_id, brand_name, route, dosage_form, strength, source, data_origin)` where inferred `brand` exists; drop the raw `product_name`.

## 7. Assemble Text Rows
26. **Non-product texts** – bind generic names, synonyms, international brands, mixture names, and mixture ingredients, then attach each drug’s `generic_name` from `primary_gi`.
27. **Desired schema** – define the final column order (`desired_cols`) ahead of the streaming writer.
28. **Product brand texts** – join `product_brand_admin` with company details and `primary_gi`, then convert each into textual rows tagged as `text_type = "brand"`, `text_subtype = "product"`.
29. **Lookup for progress labels** – cache `generic_name` by `drugbank_id` (used later while iterating).

## 8. Filter to ATC-Covered Drugs
30. **Limit scope** – compute `atc_ids` (drugs with ATC codes) and filter `non_product_texts`, `product_brand_rows`, `route_form_all`, and `combo_flag` down to those IDs only.
31. **Housekeeping** – drop the now-unused `primary_gi`, call `gc()`, and key all working tables by `drugbank_id` for fast joins.
32. **Reset output** – delete any existing `output/drugbank_authority.csv`.

## 9. Streaming Writer Primitives
33. **Empty output template** – `make_empty_output()` produces a zero-row table with the target schema, used if no data survives filtering.
34. **Unique text tracker** – builds an environment-backed set to count distinct text strings as they are emitted.
35. **Flush function** – `flush_buffer()` binds buffered chunks, fills missing columns, enforces column order, dedupes again, sorts, and appends to the CSV (writing headers only for the first batch).
36. **Buffer configuration** – initialize an empty buffer, track total rows in the buffer (`buffer_rows`), cap each flush at ~5 k rows (`target_rows`), and cap per-merge chunk sizes at ~4 k rows (`max_merge_rows`).
37. **Chunking utility** – `split_indices()` yields index lists of bounded size, used to break oversized tables into manageable chunks.

## 10. Per-Drug Processing Loop
38. **Progress handler** – wrap the loop in `progressr::with_progress()`, preparing a step for every DrugBank ID that has ATC coverage.
39. **Per-drug bookkeeping** – for each `drugbank_id`:
    - Derive a display label `"<generic name> [<id>]"` for progress messages.
    - Pull the drug’s rows from `atc_tbl`, `combo_flag`, `non_product_texts`, `route_form_all`, and `product_brand_rows`.
    - Default `combo_flag` to `FALSE` when no record exists.
40. **Chunk handler** – define `handle_chunk()` which:
    - Filters out empty text rows.
    - Joins ATC codes (allowing cartesian expansion when multiple codes exist).
    - Joins the combo flag, filling missing flags with `FALSE`.
    - Expands routes and dosage forms via `expand_route_form()`.
    - Normalizes strengths, backfills missing columns, enforces column order, dedupes, sorts, and pushes the chunk into the buffer.
    - Updates per-drug row counts, the running buffer size, and the unique text tracker.
    - Flushes the buffer to disk once the 5 k row threshold is hit.
41. **Non-product text processing** – if a drug has non-product text rows:
    - When route/form data exists, split it into chunks sized so that the merged table stays below `max_merge_rows` and send each merge through `handle_chunk()`.
    - If no route/form data exists, chunk the non-product rows themselves to stay below the cap before calling `handle_chunk()`.
42. **Product text processing** – similarly chunk product-derived rows to stay under `max_merge_rows` and feed them to `handle_chunk()`.
43. **Progress reporting** – after handling both sources, update the progress bar with the drug name plus the per-drug row count and current buffer size. If no rows were emitted, display a “no output rows” note for that drug.

## 11. Finalization
44. **Flush remainder** – after the loop, flush any residual chunks still in memory.
45. **Empty result guard** – if nothing was written, output an empty CSV with the correct header via `make_empty_output()`.
46. **Summary message** – print a final message reporting the number of rows written and the count of unique `text` values.

The resulting `drugbank_authority.csv` now contains one row per (text, drug, route/form/strength, ATC code) combination, ready for downstream authoritative SOA matching.
