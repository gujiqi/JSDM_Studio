# spOccupancy Workflow Audit Report - 2026-05-30

## Scope

This audit focused on the spOccupancy workflow in JSDM Studio. I checked the GUI wiring, current server logic, installed package API, local reference material under `C:/Users/Google/Downloads/spOccupancy参考`, real model execution, output completeness and ZIP generation.

Runtime used for verification:

- R: 4.5.3
- Rscript: `C:/Program Files/R/R-4.5.3/bin/Rscript.exe`
- spOccupancy: 0.8.0

## Main Fixes

- Added `R/spoccupancy_adapter.R` with a real spOccupancy fitting adapter.
- Connected the GUI workflow to real package execution through `Run real spOccupancy package fit`.
- Added export of `reproducible_script/run_this_spOccupancy_analysis.R`.
- Added `workflow_scripts/run_spOccupancy_workflow.R`.
- Replaced scaffold-only behavior for supported models with real calls to spOccupancy functions.
- Added CSV reshaping for:
  - single-species `sites x replicates` detection histories;
  - multi-species wide CSV columns like `sp1_rep1`, `sp1_rep2`, ... into `species x sites x replicates`;
  - integrated source CSV columns like `src1_rep1`, `src2_rep1` into spOccupancy data-source lists;
  - detection covariates from site-level, site-by-replicate and flattened formats.
- Added robust parameter validation for model/data conflicts, spatial coordinates, latent factors, MCMC settings and SVC columns.
- Added real outputs to `models/`, `samples/`, `tables/`, `predictions/`, `spatial/`, `model_assessment/`, `standard/`, `diagnostics/`, `results/`, `report/`.
- Fixed S3 prediction usage by calling `predict(output, ...)`, not `spOccupancy::predict`.
- Added failure discipline: real fit errors write `diagnostics/spOccupancy_reproducible_error.txt`, `diagnostics/engine_status.json`, standard fit-failed tables and `RUN_FAILED.txt`.
- Ensured real fitting copies uploaded files into both `inputs/` and `data/` even if the copy-inputs checkbox is disabled.
- Verified download ZIPs are real non-empty archives.

## Supported Real-Fit Functions

The adapter currently runs these spOccupancy functions for real:

- `PGOcc`
- `spPGOcc`
- `msPGOcc`
- `spMsPGOcc`
- `lfMsPGOcc`
- `intPGOcc`
- `svcPGOcc`

The UI still exposes additional advanced spOccupancy functions such as temporal, integrated multi-species, spatial-factor and binomial SVC variants. These are not faked: if selected before their adapter branch is implemented, the run fails with `fit_failed` and a diagnostic message.

## Real Test Suite

Reproducible suite:

`examples/spOccupancy/run_real_example_suite.R`

Latest passing summary:

`output/spOccupancy_real_example_suite_summary_20260530_010555.csv`

Portable packaged summary:

`examples/spOccupancy/last_real_example_suite_summary.csv`

All 7 synthetic cases completed with `status = fitted`, `exit_code = 0`, no missing required files and real ZIP outputs.

| Case | Coverage | Status |
|---|---|---|
| 01_PGOcc_single_predict_ppc_waic | single-species occupancy, detection covariate, prediction, PPC, WAIC | fitted |
| 02_spPGOcc_spatial_NNGP | spatial single-species occupancy, coords, NNGP, spatial parameters, prediction | fitted |
| 03_msPGOcc_multispecies | multi-species detection histories, species table, community priors | fitted |
| 04_spMsPGOcc_spatial_multispecies | spatial multi-species occupancy, coords, NNGP | fitted |
| 05_lfMsPGOcc_latent_factor | latent-factor multi-species occupancy, coords, newdata/newcoords prediction | fitted |
| 06_intPGOcc_integrated_sources | integrated single-species occupancy with two data sources | fitted |
| 07_svcPGOcc_spatial_varying_coefficients | spatially varying coefficient occupancy model | fitted |

ZIP validation for the latest suite:

- All 7 case ZIPs exist.
- ZIP sizes ranged from about 119 KB to 249 KB.
- ZIP entry counts ranged from 69 to 82 files.

## UI/Server Audit

- spOccupancy controls checked: 61.
- Missing server reads: none.
- Server reads without UI controls: none.
- Output widgets checked: `spocc_check_messages`, `spocc_data_table`, `spocc_files`, `spocc_log`, `spocc_run_summary`.
- Missing output render functions: none.
- Download button checked: `spocc_download`.
- Missing download handler: none.

Shiny page HTTP check:

- Status code: 200.
- Page contains `JSDM Studio`.
- Page contains `spOccupancy workflow`.
- Page contains `Run real spOccupancy package fit`.
- Page contains `Download spOccupancy ZIP`.

## Important Output Files

Successful runs write:

- `used_config.yml`
- `inputs/*.csv`
- `data/*.csv`
- `models/spOccupancy_model.rds`
- `models/spOccupancy_call_args.rds`
- `samples/posterior_sample_manifest.csv`
- `tables/model_settings_used.csv`
- `tables/summary_beta.csv`
- `tables/summary_alpha.csv`
- `tables/waicOcc_results.csv` when enabled
- `model_assessment/waicOcc_results.csv` when enabled
- `model_assessment/ppcOcc.rds` when enabled
- `results/fitted_values.rds`
- `results/spOccupancy_summary.txt`
- `standard/run_summary.csv`
- `standard/effects_long.csv`
- `standard/predictions_long.csv`
- `standard/associations_long.csv`
- `standard/fit_metrics.csv`
- `diagnostics/engine_status.json`
- `diagnostics/data_check_messages.csv`
- `diagnostics/session_info.txt`
- `reproducible_script/run_this_spOccupancy_analysis.R`
- `workflow_scripts/run_spOccupancy_workflow.R`
- `report/spOccupancy_report.html`
- `RUN_COMPLETE.txt`

## Remaining Risks

- Temporal models (`tPGOcc`, `stPGOcc`, `tMsPGOcc`, `stMsPGOcc`) require season/primary-period data structures that are not yet represented by a dedicated GUI upload format.
- `sfMsPGOcc`, `lfJSDM`, `sfJSDM`, `svcPGBinom`, temporal SVC and integrated multi-species variants still need dedicated real-fit branches.
- The synthetic suite uses intentionally tiny MCMC settings for software verification only. Publication analyses need much larger `n.batch`, `batch.length`, convergence checks and ecological model diagnostics.
- Prediction support differs by spOccupancy class. For latent-factor models, genuine new `newdata.csv` and `newcoords.csv` are required for formal prediction.
- Random effects with lme4 syntax are validated lightly; production use should verify grouping columns are numeric when spOccupancy requires numeric random-effect IDs.

## Diagnostics Guide

If a run fails, inspect these files first:

1. `diagnostics/engine_status.json`
2. `diagnostics/spOccupancy_reproducible_error.txt`
3. `diagnostics/spOccupancy_real_fit_stdout_stderr.txt`
4. `diagnostics/data_check_messages.csv`
5. `diagnostics/session_info.txt`
6. `standard/run_summary.csv`

Successful real fits should contain `RUN_COMPLETE.txt`. Failed fits should contain `RUN_FAILED.txt` and `status = fit_failed`.

