# jSDM complete-parameters and results update

This version expands the jSDM Workflow to better reflect the uploaded jSDM 0.2.7 material.

## jSDM functions covered

- `jSDM_binomial_probit()`
- `jSDM_binomial_logit()`
- `jSDM_binomial_probit_long_format()`
- `jSDM_binomial_probit_sp_constrained()`
- `jSDM_poisson_log()`
- `jSDM_gaussian()`

## Added jSDM parameter groups

### Data inputs

- Y.csv response matrix
- XData.csv / site_data
- trait_data.csv
- long_format.csv
- trials.csv
- newdata.csv
- prediction_ids.csv

### Function-specific settings

- response argument mapping: presence_data / count_data / response_data / long_format_data
- constrained latent probit option
- constrained preliminary chains
- long-format site/species/response columns
- scale site_data
- include intercept
- allow traits

### Existing and enhanced model settings

- model_type
- site_formula
- trait_formula
- n_latent
- site_effect
- trials
- burnin
- mcmc
- thin
- seed
- verbose
- ropt
- beta_start
- gamma_start
- lambda_start
- W_start
- alpha_start
- V_alpha
- V_start
- shape_Valpha
- rate_Valpha
- shape_V
- rate_V
- mu_beta
- V_beta
- mu_gamma
- V_gamma
- mu_lambda
- V_lambda

### Prediction settings

- Run predict.jSDM
- prediction type: mean / quantile / posterior
- probs for quantiles
- max prediction sites
- Id_sites
- Id_species
- prediction histograms

### Correlation and diagnostics

- HPD probability for correlations
- mean/median correlation estimate
- plot_residual_cor
- plot_associations
- beta/lambda/W/alpha/V_alpha/V/deviance diagnostics
- coda summaries

### Output tables

- model_spec table
- beta summary
- lambda summary
- gamma summary
- alpha summary
- prediction table

## jSDM downloaded ZIP structure

- `used_config.yml`: exact jSDM run configuration.
- `inputs/`: copied input files.
- `models/jsdm_model.rds`: fitted jSDM model object when production fitting is connected.
- `mcmc/`: raw posterior samples such as mcmc.sp, mcmc.gamma, mcmc.latent, mcmc.alpha, mcmc.V_alpha, mcmc.V and mcmc.Deviance.
- `tables/`: CSV summaries for model_spec, beta, lambda, gamma, alpha, variance, deviance, correlations and predictions.
- `plots/`: traceplots, density plots, correlation figures, association plots and prediction histograms.
- `diagnostics/`: engine_status.json and warning/error diagnostics.
- `results/README_jSDM_results.txt`: user-facing explanation of the output workflow.
- `report/jSDM_report.html`: HTML report.
- `tables/jSDM_result_workflow_map.csv`: map of output files and meanings.

## Hmsc result presentation

The Hmsc workflow remains separate and writes its own ZIP structure. See the Guides page for a detailed explanation of `inputs/`, `models/`, `results/`, `tables/`, `plots/`, `diagnostics/`, `report/` and `used_config.yml`.
