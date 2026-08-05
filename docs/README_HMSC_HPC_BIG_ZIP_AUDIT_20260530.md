# JSDM Studio Hmsc-HPC Big-Zip Integration Audit

Date: 2026-05-30

This build starts from `JSDMStudio_AUDITED_FIXED_20260530_ALL_WORKFLOWS_RECHECK_OK.zip` and adds the Hmsc-HPC workflow as a separate workflow next to Hmsc.

## Added

- Hmsc-HPC workflow tab in `app.R`.
- CPU-oriented Hmsc-HPC runner: `workflow_scripts/hmschpc_runner.py`.
- R wrapper/helper: `R/hmschpc_adapter.R`.
- Hmsc-HPC external reference package: `external_packages/hmsc-hpc-main`.
- Hmsc-HPC examples and multi-case test suite: `examples/HmscHPC`.
- Hmsc-HPC engine icon: `www/engine_hmschpc.svg`.
- Python dependency installer: `install_hmschpc_python_packages.bat`.
- Hmsc-HPC documentation: `docs/README_HMSC_HPC_WORKFLOW_AUDIT_20260530.md`.

## Verification Run

The following checks were run on this big source tree after integration.

- `app.R` parse/source: passed.
- Static Shiny binding audit: `ui_ids=512 input_refs=512 refs_not_in_ui=0`; `output_ui_ids=49 output_handlers=49 outputs_without_handler=0 handlers_without_ui=0`.
- Python runner compile: passed with Python 3.10.
- Installer build: passed with the staged Inno Setup build script.
- Hmsc-HPC real example suite: passed for all supported fitting cases, plus the guarded compile-only random-slope case.

## Hmsc-HPC Example Results

Eight fitted cases completed with real outputs and ZIP files:

- `01_fixed_poisson_categorical`: fitted.
- `02_fixed_probit_binary`: fitted.
- `03_fixed_normal_hmc_toggle`: fitted.
- `04_traits_phylo_cov`: fitted.
- `05_iid_random_intercept`: fitted.
- `06_spatial_full_random_intercept`: fitted.
- `07_phylo_newick_fixed`: fitted.
- `09_gaussian_alias_tfd_predictions_off`: fitted.

One intentionally guarded case completed at model-definition level:

- `08_random_slope_compile_only`: `model_defined`.

The guarded status is deliberate because the current pyhmsc/Hmsc-HPC CPU workflow does not expose complete random-slope sampling equivalent to the full Hmsc R random-level machinery. The app records this honestly instead of reporting a fake fitted model.

Additional audit fixes after the first big build:

- `Predictions = FALSE` is now respected; the predictions folder receives a disabled note and the standard comparison schema remains present.
- `standard/predictions_long.csv` now uses the same comparison columns as the other engines, including lower/upper placeholder columns and `prediction_set`.
- `save Eta = FALSE` is handled safely for the current Hmsc-HPC HDF5 exporter, which requires Eta entries. The runner forces `--fse 1` and records the warning in `diagnostics/engine_status.json`.
- Home, Project, Engine Guide, Compare Models and Guides were reread and corrected to describe Hmsc-HPC as a separate CPU pyhmsc workflow, not as full Hmsc-R or a GPU/Slurm launcher.

## Runtime Dependencies

Hmsc-HPC requires Python plus TensorFlow/TensorFlow Probability/pyhmsc. The included installer script is `install_hmschpc_python_packages.bat`. GPU mode is not exposed in this workflow; the workflow is written for CPU runs.
