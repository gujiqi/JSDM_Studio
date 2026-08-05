# sjSDM Multi-Case Real Fitting Audit - 2026-05-29

## Scope

This audit expanded sjSDM testing from one smoke example to seven small real-fitting cases. The goal was to cover the main response families, environmental modules, spatial modules, DNN settings, categorical predictors, prediction paths, output writing and ZIP generation without using large species matrices.

## Cases Run

All cases were run with:

```bat
"C:\Program Files\R\R-4.5.3\bin\Rscript.exe" examples\sjSDM\run_real_example_suite.R
```

Final suite summary:

```text
output/sjSDM_real_example_suite_summary_20260529_203736.csv
```

| Case | Family | Env module | Spatial module | Extra coverage |
|---|---|---|---|---|
| binomial_linear_categorical | binomial(probit) | linear | none | categorical predictors, newdata prediction |
| poisson_spatial_linear | poisson(log) | linear | linear | spatial predictors, new_spatial prediction, spatial coefficients |
| gaussian_dnn_env | gaussian(identity) | DNN | none | DNN env module, continuous responses |
| binomial_logit_dnn_scheduler | binomial(logit) | DNN | none | dropout, scheduler, early stopping, AdaBound |
| nbinom_intercept_only | nbinom | intercept-only | none | count response, no environmental predictors |
| binomial_spatial_ev | binomial(probit) | linear | eigenvectors | generateSpatialEV, spatial eigenvector output |
| binomial_spatial_dnn_madgrad | binomial(probit) | linear | DNN | spatial DNN, madgrad, DNN spatial weights |

## Final Result

All seven cases passed:

- exit status: `0`
- engine status: `fitted`
- errors: `0`
- warnings: `0`
- missing required outputs: none
- real ZIP files were created for every case

## Required Outputs Verified Per Case

Each case was checked for:

- `models/sjSDM_model.rds`
- `tables/coef_environment.csv`
- `tables/covariance_matrix.csv`
- `tables/correlation_matrix.csv`
- `predictions/predictions.csv`
- `tables/residuals.csv`
- `standard/run_summary.csv`
- `standard/effects_long.csv`
- `standard/predictions_long.csv`
- `standard/associations_long.csv`
- `plots/sjSDM_plot.pdf`
- case ZIP file

Spatial linear / eigenvector cases additionally save:

- `tables/coef_spatial.csv`
- `importance/spatial_importance.csv`

Spatial DNN cases save:

- `weights/spatial_dnn_coef_weights.rds`
- `importance/importance_skipped.csv`

The skip is intentional because sjSDM's linear `getImportance` formula expects linear spatial coefficients, not hidden-layer spatial DNN weights.

## Errors Found and Fixed

- Prediction with categorical `newdata` failed when the prediction subset contained only one factor level. Fixed by making prediction data inherit factor levels from training data.
- Spatial `coef(model)` returned nested `coef(model)$env[[1]]` and `coef(model)$spatial[[1]]` lists. Fixed matrix extraction and row/column naming.
- Covariance/correlation matrices sometimes lacked species names. Fixed by adding species names before writing tables and standard association outputs.
- Single DNN hidden-layer input such as `hidden = "5"` caused a PyTorch-side `TypeError: object of type 'int' has no len()`. Fixed by expanding one value to two hidden layers in the exported script.
- Spatial DNN importance was incorrectly sent to linear `getImportance`. Fixed by saving spatial DNN weights and writing an explicit skipped-importance file.
- Non-spatial DNN no longer creates a misleading `coef_spatial.csv`.

## Remaining Interpretation Note

These are small software-validation examples, not ecological inference examples. They prove that the package bridge, exported reproducible scripts, fitting calls, output tables, plots and ZIP generation work across important sjSDM parameter combinations. Formal analyses still need larger sample sizes, longer optimization, sensitivity checks and ecological validation.
