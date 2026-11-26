# s-jSDM / sjSDM Workflow Audit Report

Audit date: 2026-05-30

Scope: `app.R`, the s-jSDM Shiny workflow, exported reproducible scripts, output ZIP logic, diagnostics, standard comparison outputs, and real sjSDM 1.0.7 execution.

Local references checked:

- `C:/Users/Google/Downloads/s-jSDM参考/sjSDM.pdf`
- `C:/Users/Google/Downloads/s-jSDM参考/Help for package sjSDM.mhtml`
- `C:/Users/Google/Downloads/s-jSDM参考/sjSDM_ Getting started with sjSDM - a scalable joint Species Distribution Model.mhtml`
- Installed package API: `sjSDM 1.0.7`

## Fixed

- Verified the GUI input mapping for all `sjsdm_*` controls. No server-side `input$sjsdm_*` reference is missing from the UI; the only UI id not read as `input$` is the download button, as expected.
- Extended the exported `reproducible_script/run_this_sjSDM_analysis.R` so it now writes complete real-run outputs when executed outside Shiny:
  - `results/README_sjSDM_results.txt`
  - `results/sjSDM_run_summary.csv`
  - `report/sjSDM_report.html`
  - `diagnostics/session_info.txt`
  - `diagnostics/torch_diagnostic.txt`
  - `diagnostics/data_check_messages.csv`
  - `tables/sjSDM_input_manifest.csv`
- Fixed empty-folder behavior. Direct script runs now fill intentionally unused output directories with `README.txt`; successful example ZIPs have no empty folders.
- Added real persisted outputs for previously under-saved branches:
  - `tables/Rsquared_total.csv`, `tables/Rsquared_species.csv`, `tables/Rsquared_sites.csv`
  - `residuals/residuals.csv`
  - `anova/sjSDM_anova_results.csv`, `anova/anova_species.csv`, `anova/anova_sites.csv`, `anova/anova_plot.pdf`
  - `internal_structure/internal_structure_species.csv`, `internal_structure/internal_structure_sites.csv`, `internal_structure/internal_structure_long.csv`, `internal_structure/internal_structure_plot.pdf`
  - `internal_structure/assembly_effects.csv`, `internal_structure/assembly_effects.pdf`
  - `importance/importance_summary.csv`, `importance/importance_plot.pdf`
  - `tables/traits_metadata.csv`, `tables/species_groups.csv`, `tables/folds_uploaded.csv`
  - `tables/standard_errors.csv`, `tables/p_values.csv`
- Fixed the recursive object flattener so it skips R functions, environments and pointers inside sjSDM model/internalStructure objects. This removes the previous `cannot coerce type 'closure' to vector` failure.
- Replaced fragile sjSDM plotting-method calls with stable base-R summary plots for ANOVA/internal-structure exports.
- Added optional support/diagnostics for `sjSDM_cv` tuning and `setWeights`. The main fitted model remains reproducible even if optional package-level tuning is unavailable in a user environment.
- Updated the real example suite to 9 cases and added a stable summary file:
  - `examples/sjSDM/run_real_example_suite.R`
  - `examples/sjSDM/last_real_example_suite_summary.csv`

## Real Fit Coverage

All cases below ran with the installed `sjSDM 1.0.7` package and produced `models/sjSDM_model.rds`, standard comparison tables, diagnostics, predictions and real ZIPs.

| Case | Family | Env | Spatial | Extra coverage |
|---|---|---|---|---|
| `binomial_linear_categorical` | binomial probit | linear | none | categorical env, newdata |
| `poisson_spatial_linear` | poisson log | linear | linear | spatial predictors, ANOVA, internalStructure |
| `gaussian_dnn_env` | gaussian identity | DNN | none | DNN env, SGD |
| `binomial_logit_dnn_scheduler` | binomial logit | DNN | none | dropout, AdaBound, scheduler, early stopping |
| `nbinom_intercept_only` | nbinom | intercept-only | none | no env upload required by model design |
| `binomial_spatial_ev` | binomial probit | linear | eigenvectors | `generateSpatialEV`, retained EVs |
| `binomial_spatial_dnn_madgrad` | binomial probit | linear | DNN | spatial DNN, madgrad |
| `binomial_spatial_traits_assembly_se` | binomial probit | linear | linear | traits, groups, folds, SE, assembly effects |
| `binomial_cv_tuning` | binomial probit | linear | none | `sjSDM_cv`, uploaded folds, tuning output |

Latest suite summary:

`output/sjSDM_real_example_suite_summary_20260530_020453.csv`

Result: 9/9 fitted, 0 missing required outputs, 0 empty directories, 0 warnings.

## Verification

- `app.R` parse check: passed.
- `examples/sjSDM/run_real_example_suite.R`: passed.
- Shiny smoke check: `http://127.0.0.1:6894` returned HTTP 200 and contained the sjSDM workflow text.
- ZIP checks: all 9 example ZIPs are real files with non-zero sizes.

## Residual Risks

- `sjSDM_cv` is part of sjSDM 1.0.7, but it is sensitive to torch/PyTorch and to cross-validation folds where a species column has no 0/1 variation. The adapter records tuning diagnostics and keeps the main model fit reproducible.
- GPU mode is exposed only as `device = "gpu"` and still depends on a correctly installed CUDA-enabled torch stack. CPU mode is the tested default.
- `p_values.csv` is written as a documented compatibility table because sjSDM 1.0.7 exposes standard-error objects through `getSe`; stable p-value extraction is not a public API contract.
- Traits are not part of the core `sjSDM()` likelihood. They are used for metadata, reporting and assembly-effect interpretation where supported by sjSDM outputs.

## Minimal Test

From the app root:

```powershell
& 'C:\Program Files\R\R-4.5.3\bin\Rscript.exe' examples\sjSDM\run_real_example_suite.R
```

If a run fails, inspect:

- `diagnostics/engine_status.json`
- `diagnostics/session_info.txt`
- `diagnostics/torch_diagnostic.txt`
- `diagnostics/sjSDM_reproducible_error.txt`
- `diagnostics/sjSDM_cv_error.txt` for optional tuning issues
