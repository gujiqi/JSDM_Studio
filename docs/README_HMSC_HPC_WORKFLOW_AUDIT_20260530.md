# Hmsc-HPC workflow audit and integration report

Date: 2026-05-30

## Scope

Hmsc-HPC was added as a separate JSDM Studio workflow placed beside the R Hmsc workflow. The GUI exposes the CPU-supported pyhmsc/Hmsc-HPC path only; GPU execution is intentionally not exposed.

## Implemented workflow

- New tab: `Hmsc-HPC Workflow`.
- New engine guide card and project recommendation scoring.
- New icon: `www/engine_hmschpc.svg`.
- Bundled source: `external_packages/hmsc-hpc-main`.
- New adapter: `R/hmschpc_adapter.R`.
- New Python runner: `workflow_scripts/hmschpc_runner.py`.
- New install helper: `install_hmschpc_python_packages.bat`.
- New reproducible test suite: `examples/HmscHPC/run_real_example_suite.R`.
- New uploadable example data: `examples/HmscHPC/minimal/`.

## Parameter coverage

The workflow exposes and validates:

- Response distribution: `poisson`, `probit`, `normal` / `gaussian`.
- X formula, including categorical predictors through patsy syntax such as `C(substrate)`.
- Traits and trait formula.
- Phylogeny by covariance matrix or Newick tree.
- Random-level mode: none, `iid`, `spatial_full`, guarded `random_slope_iid` compile-only.
- Random-level name, grouping column, coordinate columns, random-slope formula, nf/nfMin/nfMax and spatial alpha.
- CPU sampler settings: samples, transient, thin, chains, chains-to-run, verbose, seed, precision, truncated-normal backend, HMC leapfrog/thin, update Beta/Eta, save Eta, eager/profile toggles.
- Output settings: predictions, diagnostics, plots, ZIP. The predictions flag now controls `predictions/` export while still preserving an empty-standard `standard/predictions_long.csv` for comparison compatibility.
- Python executable and Hmsc-HPC source directory.

## Real run results

The suite at `examples/HmscHPC/run_real_example_suite.R` was run with R 4.5.3 and Python 3.10 on CPU.

Latest summary CSV:

`output/HmscHPC_real_example_suite_summary_20260530_141510.csv`

Cases:

- `01_fixed_poisson_categorical`: fitted.
- `02_fixed_probit_binary`: fitted.
- `03_fixed_normal_hmc_toggle`: fitted.
- `04_traits_phylo_cov`: fitted.
- `05_iid_random_intercept`: fitted.
- `06_spatial_full_random_intercept`: fitted.
- `07_phylo_newick_fixed`: fitted.
- `08_random_slope_compile_only`: model_defined, as expected, because Hmsc-HPC native sampling guards random slopes.
- `09_gaussian_alias_tfd_predictions_off`: fitted. This case covers `gaussian` -> `normal` normalization, `tnlib = tfd`, eager/profile flags, `save_eta = FALSE` safety handling and `predictions = FALSE`.

All cases produced real ZIP files, no missing required output files, and no empty output directories. The fitted cases now produce 82-100 files each, including Hmsc-style S1-S7 result tables.

## Important fixes made during audit

- Fixed YAML relative paths. pyhmsc resolves paths relative to the YAML location, so workflow YAML now uses `../data/...` because it is stored in `workflow_scripts/`.
- Guarded `verbose = 0`; current Hmsc-HPC TensorFlow sampler fails on modulo by zero, so the GUI and validator require verbose >= 1.
- Fixed reproducible script Python path handling on Windows.
- Added newdata prediction export to `predictions/newdata_predicted_mean.csv`.
- Expanded fitted Hmsc-HPC outputs to match the classic Hmsc result style: `results/S1_model_definition.csv`, `S2_fit_models.csv`, `S3_convergence_summary.csv`, `S4_model_fit_summary.csv`, `S5_model_fit_prediction_summary.csv`, `S6_parameter_estimates_Beta.csv` and `S7_predictions_training.csv`.
- Added richer posterior exports: Beta support tables, Beta/Gamma mean matrices, sigma/rho summaries, random-level Eta/Lambda summaries, lambda-derived association matrices and `tables/HmscHPC_S1S7_result_index.csv`.
- Added `diagnostics/data_check_messages.csv`, `diagnostics/session_info.txt`, `diagnostics/posterior_hdf5_summary.csv` and `standard/diagnostics_long.csv`.
- Added exported `workflow_scripts/S1_define_models.py` through `S7_make_predictions.py` as auditable Hmsc-HPC step scripts.
- Made `Predictions` a real output switch. When disabled, `predictions/` contains `predictions_disabled.txt`, while `standard/predictions_long.csv` remains present with the standard schema.
- Standardized Hmsc-HPC `predictions_long.csv` to include `engine`, `site_id`, `response_id`, `observed`, `predicted_mean`, `predicted_lower`, `predicted_upper` and `prediction_set`.
- Guarded the HDF5 `save Eta = FALSE` edge case. The current Hmsc-HPC HDF5 writer expects Eta entries, so the runner forces `--fse 1` and writes a warning instead of failing after sampling.
- Fixed ArviZ/xarray diagnostic export to CSV for Beta Rhat and ESS.
- Fixed Hmsc-HPC Python runtime text-input layout so `help_text()` is not passed as a `textInput()` width argument.
- Restored valid Shiny UI closure after the guide expansion; `app.R` now parses and sources cleanly.
- Removed obsolete global preview outputs from `app.R`; each workflow now owns its data preview block.
- Updated Home, Project, Engine Guide, Compare Models and Guides so Hmsc-HPC is described as a separate CPU pyhmsc/HDF5 workflow beside, not inside, classic R Hmsc.
- Added static Shiny binding audit script and verified:
  - input IDs referenced by server: 512
  - UI input IDs: 512
  - missing input IDs: 0
  - output UI IDs: 49
  - output handlers: 49
  - missing output handlers: 0

## Output contract

Each successful Hmsc-HPC run writes:

- `used_config.yml`
- `inputs/`
- `data/`
- `models/compiled_model/init.json`
- `models/compiled_model/init_arrays.h5`
- `samples/posterior.h5` for fitted sampler runs
- `tables/Beta_summary.csv`
- `tables/Beta_support.csv`
- `tables/Beta_mean_matrix.csv`
- `tables/Gamma_mean_matrix.csv` when Gamma is available
- `tables/random_level_*_Eta_mean.csv` and `tables/random_level_*_Lambda_mean.csv` when random levels are fitted
- `tables/random_level_*_association_correlation.csv` when Lambda samples are available
- `tables/HmscHPC_S1S7_result_index.csv`
- `results/S1_model_definition.csv`
- `results/S2_fit_models.csv`
- `results/S3_convergence_summary.csv`
- `results/S4_model_fit_summary.csv`
- `results/S5_model_fit_prediction_summary.csv`
- `results/S6_parameter_estimates_Beta.csv`
- `results/S7_predictions_training.csv`
- `diagnostics/engine_status.json`
- `diagnostics/data_check_messages.csv`
- `diagnostics/session_info.txt`
- `diagnostics/posterior_hdf5_summary.csv`
- `diagnostics/HmscHPC_compile.log`
- `diagnostics/HmscHPC_validate_init.log`
- `diagnostics/HmscHPC_sample.log` for fitted sampler runs
- `predictions/predicted_mean.csv` when prediction export is enabled
- `predictions/training_predictions_long.csv` when prediction export is enabled
- `predictions/newdata_predicted_mean.csv` when prediction export is enabled and newdata is uploaded
- `predictions/predictions_disabled.txt` when prediction export is disabled
- `standard/run_summary.csv`
- `standard/effects_long.csv`
- `standard/predictions_long.csv`
- `standard/fit_metrics.csv`
- `standard/diagnostics_long.csv`
- `standard/associations_long.csv`
- `standard/output_manifest.csv`
- `workflow_scripts/hmschpc_model.yaml`
- `workflow_scripts/S1_define_models.py` through `workflow_scripts/S7_make_predictions.py`
- `reproducible_script/run_this_HmscHPC_analysis.R`
- `reproducible_script/run_this_HmscHPC_analysis.py`
- `report/Hmsc-HPC_report.html`

## Remaining scientific/software risks

- Hmsc-HPC is a Python/TensorFlow research code path, not a CRAN R package.
- GPU execution is not exposed in this workflow by design.
- Random slopes compile but are not sampled by the native sampler; status is `model_defined` when sampler is disabled.
- GPP/NNGP Hmsc spatial modes are not native Hmsc-HPC modes here; use R Hmsc for those.
- Quick-test MCMC settings only prove execution and output integrity, not convergence.

## Failure diagnostics

If a Hmsc-HPC run fails, inspect in this order:

1. `diagnostics/engine_status.json`
2. `diagnostics/HmscHPC_error.txt`
3. `diagnostics/HmscHPC_compile.log`
4. `diagnostics/HmscHPC_validate_init.log`
5. `diagnostics/HmscHPC_sample.log`
6. `diagnostics/HmscHPC_Rscript_stdout.log`
7. `diagnostics/HmscHPC_Rscript_stderr.log`
