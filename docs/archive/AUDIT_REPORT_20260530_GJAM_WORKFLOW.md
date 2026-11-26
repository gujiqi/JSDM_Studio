# GJAM Workflow Audit Report - 2026-05-30

## Scope

This audit focused on the GJAM workflow in JSDM Studio. The workflow was checked against the installed `gjam` package and local reference files under `C:/Users/Google/Downloads/GJAM参考`.

Runtime used for verification:

- R: 4.5.3
- Rscript: `C:/Program Files/R/R-4.5.3/bin/Rscript.exe`
- gjam: 2.7

## Main Fixes

- Added a real GJAM adapter in `R/gjam_adapter.R` that writes and runs a standalone executable script: `reproducible_script/run_this_GJAM_analysis.R`.
- Connected the GJAM GUI workflow to real `gjam::gjam()` fitting when `Run real GJAM package fit` is enabled.
- Fixed `~ .` formula handling by expanding it to explicit predictors before calling `gjam()`.
- Fixed CSV row-name handling so exported or uploaded files with a leading row-id column do not create extra response columns.
- Added typed response handling so `CAT` columns remain categorical while all non-CAT response columns are converted to numeric.
- Fixed `gjamCensorY()` integration by replacing the censored columns back into the complete Y matrix instead of dropping all other responses.
- Added safe fallback for `gjamPredict()` with `newdata`: if gjam 2.7 fails for a mixed response model, the script writes a warning and falls back to fitted-data prediction.
- Added safe handling for `gjamPriorTemplate()`:
  - sign-constrained priors are supported for non-CAT response models;
  - CAT models skip betaPrior injection with a warning because gjam 2.7 expands categorical responses internally.
- Added safe handling for `gjamTrimY()`:
  - automatic trimming is only applied to single-type `CC` composition-count data;
  - other types run without trimming and write a warning.
- Added full machine-readable outputs under `standard/`, `tables/`, `results/`, `models/`, `chains/`, `predictions/`, `diagnostics/`, `report/`, `workflow_scripts/`, and `reproducible_script/`.
- Fixed GJAM ZIP behavior: each real run now produces a real ZIP containing the output tree and `standard/output_manifest.csv`.
- Audited GJAM UI/server wiring:
  - 61 GJAM controls checked;
  - no missing `input$gjam_*` server reads;
  - no server-only `input$gjam_*` references without UI controls;
  - all GJAM outputs and download controls have server definitions.

## Real Test Suite

Reproducible suite:

`examples/GJAM/run_real_example_suite.R`

Latest passing summary:

`output/GJAM_real_example_suite_summary_20260530_004109.csv`

Portable packaged summary:

`examples/GJAM/last_real_example_suite_summary.csv`

All 7 synthetic cases completed with `status = fitted`, `exit_code = 0`, non-empty required files, and real ZIP outputs.

| Case | Coverage | Status |
|---|---|---|
| 01_mixed_scales_default_formula_predict_sensitivity | Mixed `CON, PA, CA, DA, OC`, `~ .`, newdata fallback, sensitivity, plots | fitted |
| 02_DA_effort_random_FULL_notStandard | DA counts, effort, random factor, FULL, PREDICTX off, notStandard | fitted |
| 03_CA_censor_holdout_conditional | CA abundance, censor file/text settings, holdoutN, sign-constrained prior, conditional parameters | fitted |
| 04_FC_composition_groups_ordination | FC composition, FCgroups, holdoutIndex, ordination | fitted |
| 05_CC_composition_groups_iie | CC composition counts, CCgroups, IIE | fitted |
| 06_CAT_categorical_safe_prior_skip | CAT categorical responses, safe prior skip, prediction, sensitivity | fitted |
| 07_traits_missing_trimY | traits/specByTrait/traitTypes, missing X/Y, trimY safety guard | fitted |

ZIP validation for the latest suite:

- Case ZIPs existed for all 7 cases.
- ZIP sizes ranged from about 170 KB to 515 KB.
- ZIP entry counts ranged from 78 to 111 files.

## Important Output Files

Every passing case writes the following key files:

- `used_config.yml`
- `inputs/*.csv`
- `data/*.csv`
- `models/gjam_model.rds`
- `models/gjam_modelList.rds`
- `chains/chain_manifest.csv`
- `tables/typeNames_used.csv`
- `tables/betaMu.csv`
- `tables/corMu.csv`
- `tables/fit_DIC_rmspe_xscore_yscore.csv`
- `standard/run_summary.csv`
- `standard/effects_long.csv`
- `standard/predictions_long.csv`
- `standard/associations_long.csv`
- `standard/fit_metrics.csv`
- `results/run_summary.csv`
- `results/fit_metrics.csv`
- `results/README_GJAM_results.txt`
- `predictions/gjam_prediction_from_fit.rds`
- `diagnostics/engine_status.json`
- `diagnostics/data_check_messages.csv`
- `diagnostics/session_info.txt`
- `workflow_scripts/run_GJAM_workflow.R`
- `reproducible_script/run_this_GJAM_analysis.R`
- `report/GJAM_report.html`
- `RUN_COMPLETE.txt`

## Remaining Risks

- `gjamPredict(output, newdata=...)` can fail in gjam 2.7 for some mixed response models. The adapter now falls back to fitted-data prediction and records the warning in `diagnostics/engine_status.json`.
- `gjamPriorTemplate()` betaPrior is not injected for CAT response models because gjam expands categorical columns internally. Use sign-constrained priors with non-CAT models.
- `gjamTrimY()` is only applied automatically for single-type `CC` composition-count data. Other response types are safer without trimming.
- REDUCT with very small response matrices, especially CAT responses, can be unstable in gjam 2.7. The checker warns users before fitting.
- The quick synthetic cases use small `ng` and `burnin` values for software testing only. Publication runs need larger MCMC settings and convergence review.

## Diagnostics Guide

If a run fails, inspect these files in this order:

1. `diagnostics/engine_status.json`
2. `diagnostics/GJAM_reproducible_error.txt`
3. `diagnostics/GJAM_real_fit_stdout_stderr.txt`
4. `diagnostics/data_check_messages.csv`
5. `diagnostics/session_info.txt`
6. `standard/run_summary.csv`

Successful runs should contain `RUN_COMPLETE.txt`. Failed fits should contain `RUN_FAILED.txt` and `status = fit_failed`.
