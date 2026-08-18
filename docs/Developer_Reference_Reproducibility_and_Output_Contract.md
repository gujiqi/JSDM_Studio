# JSDMWorkbench Developer Reference, Reproducibility and Output Contract

This document is intended for software reviewers, package maintainers and developers who need more detail than a user manual. It describes the GUI-to-script contract, expected output files and the included completed example.

## Reproducibility Contract

Every GUI run must export:

- `used_config.yml`: exact GUI settings and selected model options.
- `inputs/` and `data/`: copied input files used by the run.
- `diagnostics/engine_status.json`: machine-readable status. Failed fits must use `fit_failed` or `check_failed`, not `Completed`.
- `diagnostics/data_check_messages.csv`: validation messages.
- `reproducible_script/run_this_<engine>_analysis.R`: executable script for regenerating the exported scaffold or production outputs.
- `workflow_scripts/`: engine-specific step scripts or workflow wrapper.
- `standard/`: comparison tables used by Compare Models.
- `report/`: browser-readable run report.

For GUI-executed software, reviewers should be able to run the exported script from the output folder or by calling `Rscript reproducible_script/run_this_<engine>_analysis.R`.

## Completed Example

The package includes a completed example at:

`examples/completed_HMSC_spatial_full/`

This example was generated from the bundled HMSC example data with a minimal two-sample MCMC smoke run. It includes fitted HMSC workflow artifacts, S1-S7 scripts, standard tables, diagnostics and a reproducible script. It is intentionally small so reviewers can inspect the full output contract quickly.

## Example Data by Engine

- `examples/HMSC/`: `Y.csv`, `XData.csv`, `traits.csv`, `studyDesign.csv`, `coordinates.csv`.
- `examples/jSDM/`: `Y.csv`, `XData.csv`, `trait_data.csv`, `trials.csv`, `newdata.csv`, `prediction_ids.csv`.
- `examples/GJAM/`: mixed-scale `Y.csv`, `XData.csv`, `typeNames.csv`, trait type files and prediction/holdout examples.
- `examples/spOccupancy/`: detection/nondetection `y.csv`, occurrence and detection covariates, coordinates, species metadata and folds.
- `examples/sjSDM/`: response matrix, environmental predictors, spatial predictors, traits, prediction data and folds.
- `examples/boral/`: response matrix, covariates, traits, row/random-effect IDs, distance matrix, offset, prediction data and trial sizes.

## Engine Output Contract

### HMSC

HMSC is the production-connected engine. It writes real S1-S7 outputs when the Hmsc package is available and validation passes:

- S1: model definition and random-level summary.
- S2: fitted model objects from `sampleMcmc`.
- S3: MCMC convergence diagnostics.
- S4: model fit and cross-validation objects.
- S5: model fit plots.
- S6: parameter estimates, variance partitioning and association outputs.
- S7: gradient predictions.

### jSDM

The GUI currently exports a diagnostic scaffold and output contract. Production fitting should populate `models/`, `mcmc/`, `tables/`, `predictions/` and `plots/` with jSDM-specific posterior samples, coefficients, latent factors, residual/environmental correlations and predictions.

### GJAM

The GUI exports GJAM-specific output manifests for mixed-scale response models. Production fitting should populate chains, observation-scale parameters, covariance/correlation matrices, sensitivity, inverse prediction, ordination, missing-data imputation and plotting outputs.

### spOccupancy

The GUI exports occupancy-specific output manifests. Production fitting should populate occurrence and detection parameter summaries, posterior samples, fitted values, prediction tables, PPC/WAIC/k-fold assessment, spatial random effects, latent factors and SVC outputs.

### sjSDM

The GUI exports sjSDM-specific output manifests. Production fitting should populate coefficients, covariance/correlation matrices, standard errors, predictions, R-squared tables, ANOVA/variation partitioning, internal metacommunity structure, importance, weights and residual outputs.

### boral

The GUI exports boral-specific output manifests. Production fitting should populate JAGS model code, fitted boral object, MCMC samples, posterior summaries, ordination outputs, Dunn-Smyth residuals, coefficient/correlation tables, random effects, SSVS results, predictions and model comparison metrics.

## Scaffold Policy

When an engine is not production-connected or a dependency is missing, the run must:

- use `model_defined`, `check_failed` or `fit_failed`;
- never report `Completed`;
- write non-empty diagnostic files;
- write engine-specific result manifests;
- write an executable reproducibility script;
- preserve enough metadata for developers to connect a production adapter without changing the GUI contract.
