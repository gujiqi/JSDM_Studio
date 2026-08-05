# Hmsc-HPC Hmsc-style output expansion

Date: 2026-05-30

This pass specifically addressed the question: can Hmsc-HPC produce real result folders comparable to the classic Hmsc workflow?

## Answer

Yes. Hmsc-HPC fitted runs now produce a real posterior file and a richer Hmsc-style result package. The output is not a fake scaffold: fitted cases run the pyhmsc compiler, validate the compiled model boundary, run the CPU TensorFlow sampler, read `samples/posterior.h5`, then export posterior summaries, predictions, diagnostics, plots and standard comparison tables.

## Added output files

Each fitted Hmsc-HPC run now writes:

- `results/S1_model_definition.csv`
- `results/S2_fit_models.csv`
- `results/S3_convergence_summary.csv`
- `results/S4_model_fit_summary.csv`
- `results/S5_model_fit_prediction_summary.csv`
- `results/S6_parameter_estimates_Beta.csv`
- `results/S6_Beta_support.csv`
- `results/S7_predictions_training.csv`
- `tables/Beta_summary.csv`
- `tables/Beta_support.csv`
- `tables/Beta_mean_matrix.csv`
- `tables/Gamma_mean_matrix.csv` when available
- `tables/sigma_mean.csv` when available
- `tables/rho_mean.csv` when available
- `tables/random_level_*_Eta_mean.csv` when a random level is fitted
- `tables/random_level_*_Lambda_mean.csv` when a random level is fitted
- `tables/random_level_*_association_correlation.csv` when Lambda samples are available
- `tables/HmscHPC_S1S7_result_index.csv`
- `diagnostics/data_check_messages.csv`
- `diagnostics/session_info.txt`
- `diagnostics/posterior_hdf5_summary.csv`
- `diagnostics/S3_convergence_summary.csv`
- `standard/diagnostics_long.csv`

The existing reproducible files remain:

- `used_config.yml`
- `workflow_scripts/hmschpc_model.yaml`
- `workflow_scripts/S1_define_models.py` through `workflow_scripts/S7_make_predictions.py`
- `reproducible_script/run_this_HmscHPC_analysis.R`
- `reproducible_script/run_this_HmscHPC_analysis.py`
- `report/Hmsc-HPC_report.html`

## Validation

The stricter Hmsc-HPC example suite was rerun after the expansion.

- Summary CSV: `output/HmscHPC_real_example_suite_summary_20260530_141510.csv`
- Cases run: 9
- Fitted cases: 8
- Guarded compile-only case: 1 (`random_slope_iid`, status `model_defined`)
- Missing required files: none
- Empty required directories: none
- Fitted case file counts: 82-100 files per output folder

## Scientific limits

The expanded outputs are Hmsc-style in organization, but they do not claim that Hmsc-HPC implements every classic R Hmsc diagnostic. Classic R Hmsc remains the engine for GPP/NNGP, full variance partitioning and full Hmsc Omega interpretation. Hmsc-HPC exports Eta/Lambda random-level summaries and lambda-derived association tables where the CPU pyhmsc posterior contains those samples.
