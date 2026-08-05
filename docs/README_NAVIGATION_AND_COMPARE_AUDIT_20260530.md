# Navigation and Compare Models audit

Date: 2026-05-30

This pass audited the non-fitting navigation pages: Home, Project, Engine Guide, Compare Models and Guides.

## Fixes

- Home now frames JSDM Studio as a developer-reference workflow system with explicit input, parameter, output, status and reproducible-script contracts.
- Project guidance now avoids implying that zeros or large species matrices have one universal solution. It states that each engine handles zeros, high-dimensional response matrices and traits differently.
- Engine Guide now describes Hmsc-HPC as a CPU pyhmsc/HDF5 workflow with Hmsc-style outputs, not as full R Hmsc, GPU Hmsc-HPC or Slurm submission.
- Compare Models now reports status classes:
  - `fitted_result`
  - `model_defined_no_posterior`
  - `failed_diagnostic_output`
  - `missing_folder`
  - `not_ready_or_unknown`
- Compare Models now checks more of the developer output contract:
  - `diagnostics/engine_status.json`
  - `used_config.yml`
  - `standard/run_summary.csv`
  - `standard/fit_metrics.csv`
  - `standard/effects_long.csv`
  - `standard/predictions_long.csv`
  - `standard/associations_long.csv`
  - `standard/diagnostics_long.csv`
  - `standard/output_manifest.csv`
  - workflow step maps or engine-specific workflow maps
  - report HTML
  - adjacent ZIP file
- Compare Models now distinguishes fitted runs from model-defined or failed diagnostic folders before metric comparison.
- Guides now include a dedicated status and comparison guide explaining why fitted, model_defined, check_failed and fit_failed must not be interpreted the same way.
- `R/engine_registry.R` now registers all seven engines, not only Hmsc and jSDM.
- Added `workflow_scripts/test_compare_models.R` to exercise the Compare Models server logic with a real Hmsc-HPC output folder.

## Validation

The following checks passed after the audit:

- `app.R` parse: passed.
- `app.R` source: passed.
- Static Shiny binding audit: 512 UI inputs / 512 server input references; 49 UI outputs / 49 output handlers.
- Engine registry audit: 7 registered engines.
- Compare Models test: passed with a real fitted Hmsc-HPC output folder.

## Scientific cautions now stated in the UI

- Hmsc Omega, Hmsc-HPC Lambda-derived associations, jSDM residual correlations, GJAM corMu/sigMu, spOccupancy latent factors, sjSDM covariance and boral residual correlations are not numerically interchangeable.
- `model_defined` is a valid reproducibility artifact but not a fitted posterior model.
- `check_failed` and `fit_failed` folders are diagnostic artifacts, not successful analyses.
- Prediction metrics should only be compared when validation design, response family and response scale are compatible.
