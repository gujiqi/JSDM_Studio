# Universal Benchmark real-input and image display fix

Date: 2026-06-03

## What changed

- Workflow and engine images are now embedded as SVG data URIs from `www/` assets during Shiny UI construction. This avoids broken images in browser mode, WebView2 mode, installed mode, and ZIP-expanded paths.
- Universal Benchmark now supports four data sources:
  - generate benchmark data from one latent ecological truth;
  - select a built-in real-style benchmark case;
  - upload real benchmark files one by one;
  - advanced ZIP import for already prepared benchmark folders.
- The recommended real-data path is now one file per upload control. Each control has an adjacent annotation explaining role, format, engine use and a concrete example. This avoids hiding missing or mismatched files inside a ZIP.
- Uploaded single files are assembled into `input/universal_benchmark_uploads/<timestamp>/` and validated against the full benchmark input contract before fitting.
- Uploaded ZIPs remain available only as an advanced import mode. They are validated against the full benchmark input contract before fitting and may contain files at the ZIP root or inside one top-level folder.
- Universal Benchmark runner now accepts `--input_dir=` so real or built-in benchmark inputs are used directly instead of being overwritten by generated data.
- Master benchmark input copying now preserves subdirectories such as `truth/`.

## Required benchmark input contract

The input directory must contain `Y_occurrence.csv`, `Y_count.csv`, `Y_normal.csv`, `XData.csv`, `traits.csv`, `studyDesign.csv`, `coordinates.csv`, `phylogeny.nwk`, `phylo_cov.csv`, `newdata.csv`, `newcoords.csv`, `folds.csv`, `trial.size.csv`, `offset.csv`, `row.ids.csv`, `ranef.ids.csv`, `distmat.csv`, `spOccupancy_y_detection.csv`, `occ.covs.csv`, `det.covs.csv`, and the `truth/` files used for effect and prediction comparison.

## Built-in real-style cases

The built-in cases are stored in `examples/universal_benchmark_cases/`.

| Case | Sites | Species | Visits | Validation |
|---|---:|---:|---:|---|
| case_01_forest_gradient_presence | 36 | 6 | 3 | all seven engines fitted |
| case_02_wetland_detection_replicates | 34 | 6 | 4 | all seven engines fitted |
| case_03_alpine_categorical_substrate | 38 | 7 | 3 | all seven engines fitted |
| case_04_coastal_phylogeny_traits | 32 | 5 | 3 | all seven engines fitted |
| case_05_grassland_count_offset | 40 | 8 | 3 | all seven engines fitted |
| case_06_restoration_newsite_prediction | 36 | 7 | 5 | all seven engines fitted |
| case_07_small_publication_smoke | 30 | 5 | 3 | all seven engines fitted |

Full run evidence is in `examples/universal_benchmark_cases/real_style_case_run_results.csv`.

## Verification

- `app.R` parse: passed.
- `app.R` source: passed.
- Strict Shiny input/output audit: passed.
- HTTP image smoke test: returned HTML with 20 embedded SVG data URIs.
- Seven built-in real-style benchmark cases: all engines fitted for each case.
