# Book Chapter Benchmark Integration Audit

Date: 2026-08-18

## Scope

Converted the Chapter 6, Chapter 7 and Chapter 11 example data from the local book example folder into JSDM Studio Universal Benchmark-compatible data sets, using the book PDF chapter text only as source documentation.

## Added Examples

- `examples/book_chapter_benchmarks/ch06_plant_traits_whittaker`
- `examples/book_chapter_benchmarks/ch07_deadwood_fungi`
- `examples/book_chapter_benchmarks/ch11_finnish_birds`
- `examples/book_chapter_benchmarks/raw_source`
- `workflow_scripts/prepare_jsdm_book_chapter_cases.R`

Each converted case has 36 sites and 6 species, no all-zero/all-one species, full-rank benchmark design matrices and the complete Universal Benchmark input contract.

## Fixes Made During Real Runs

- Fixed book long-table conversion so repeated site/species labels are preserved and only final matrix row/column names are uniquified.
- Fixed Chapter 6 `xtabs` response conversion by explicitly rebuilding a numeric matrix from table dimensions and dimnames.
- Adjusted derived predictors to avoid exact linear dependence in GJAM formula matrices.
- Hardened CSV readers in generated jSDM, sjSDM and boral scripts so numeric-looking site IDs and empty first-column names are treated as row identifiers.
- Added a final Hmsc `studyDesign` factor guard immediately before constructing the Hmsc model.

## Real Run Results

All three book-derived cases were run through Hmsc, Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM and boral via `workflow_scripts/universal_benchmark_runner.R --action=all --input_dir=...`.

| Case | Engines fitted | Master ZIP |
| --- | ---: | --- |
| Chapter 6 plant traits / Whittaker gradient | 7/7 | `output/universal_benchmark_20260818_153015/JSDMStudio_universal_benchmark_ALL_MODELS.zip` |
| Chapter 7 dead wood fungi | 7/7 | `output/universal_benchmark_20260818_153132/JSDMStudio_universal_benchmark_ALL_MODELS.zip` |
| Chapter 11 Finnish birds | 7/7 | `output/universal_benchmark_20260818_153244/JSDMStudio_universal_benchmark_ALL_MODELS.zip` |

The output-contract audit found zero missing required engine files and zero missing required master comparison files. The machine-readable audit table is `output/book_chapter_benchmark_run_audit_20260818.csv`.

## Notes and Residual Risk

- The Chapter 7 spOccupancy detection data are pseudo-replicates for software validation only.
- The `truth/` tables are empirical reference summaries from observed occurrence, not hidden simulation truth.
- These are quick benchmark runs with very short MCMC or iteration settings; they verify execution and output contracts, not publication-quality inference.
