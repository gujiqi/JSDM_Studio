# HMSC S1-S7 real results workflow

This build changes HMSC from scaffold-only output to a real S1-S7 results workflow.

The HMSC panel now runs, by default, the same result sequence as the uploaded scripts:

S1_define_models_template:
- creates `models/unfitted_models.RData`

S2_fit_models:
- fits HMSC models with `sampleMcmc`
- saves `models/models_thin_<thin>_samples_<samples>_chains_<chains>.Rdata`

S3_evaluate_convergence:
- computes coda/Gelman PSRF diagnostics where possible
- saves `results/MCMC_convergence.txt`, `results/MCMC_convergence.pdf`, and CSV diagnostics

S4_compute_model_fit:
- computes `computePredictedValues`, `evaluateModelFit`, CV fit and WAIC where possible
- saves `models/MF_thin_<thin>_samples_<samples>_chains_<chains>_nfolds_<nfolds>.Rdata`
- saves `tables/S4_model_fit_summary.csv`

S5_show_model_fit:
- saves `results/model_fit_nfolds_<nfolds>.pdf`

S6_show_parameter_estimates:
- saves `results/parameter_estimates.pdf`
- saves Beta/Gamma/Omega/variance-partitioning CSV outputs where available

S7_make_predictions:
- saves `results/predictions.pdf`
- saves gradient prediction RDS files in `predictions/`

The original uploaded S1-S7 scripts are copied into:
`workflow_templates/HMSC_user_S1S7_original/`

Each run also copies those original scripts into the downloaded ZIP:
`workflow_scripts/`
