# HMSC S1-S7 real code fix

This build fixes the previous issue where the log showed `HMSC S1-S7 results workflow status: not_started`.

Cause:
The error handler inside `run_hmsc_s1s7_pipeline()` used local assignment, so if S1/S2 failed, the returned status stayed `not_started`.

Fix:
- The error handler now correctly returns `fit_failed`.
- The exact error is written to `diagnostics/HMSC_S1S7_error.txt`.
- The pipeline now generates executable scripts:
  - workflow_scripts/S1_define_models.R
  - workflow_scripts/S2_fit_models.R
  - workflow_scripts/S3_evaluate_convergence.R
  - workflow_scripts/S4_compute_model_fit.R
  - workflow_scripts/S5_show_model_fit.R
  - workflow_scripts/S6_show_parameter_estimates.R
  - workflow_scripts/S7_make_predictions.R
  - reproducible_script/run_all_HMSC_S1_to_S7.R
- Input data are copied to `data/` so these scripts can be run from the downloaded ZIP.

Expected successful output:
- models/unfitted_models.RData
- models/models_thin_<thin>_samples_<samples>_chains_<chains>.Rdata
- results/MCMC_convergence.pdf
- results/model_fit_nfolds_<nfolds>.pdf
- results/parameter_estimates.pdf
- results/predictions.pdf
