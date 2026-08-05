
# jSDM engine guide

This package includes a new **jSDM engine** for the workbench.

## What jSDM does

`jSDM` fits joint species distribution models in a hierarchical Bayesian framework. Its Gibbs sampler is written in C++ and uses Rcpp, Armadillo and GSL for computational efficiency.

## Supported jSDM model types in this GUI

```text
binomial_probit
binomial_logit
binomial_probit_sp_constrained
poisson_log
gaussian
```

## Main jSDM inputs

```text
Y.csv      response matrix: rows = sites, columns = species or responses
XData.csv  site-level predictors
TrData.csv optional species traits / response attributes
```

## Main jSDM parameters shown in the GUI

```text
model_type
site_formula
use_traits
trait_formula
n_latent
site_effect
burnin
mcmc
thin
beta_start
gamma_start
lambda_start
W_start
alpha_start
V_alpha
V_start
shape_Valpha
rate_Valpha
shape_V
rate_V
mu_beta
V_beta
mu_gamma
V_gamma
mu_lambda
V_lambda
ropt
seed
verbose
```

## Main jSDM outputs

```text
models/jsdm_model.rds
tables/jsdm_beta_estimates.csv
tables/jsdm_gamma_estimates.csv
results/residual_correlation.rds
tables/residual_correlation_mean.csv
results/environmental_correlation.rds
tables/environmental_correlation_mean.csv
diagnostics/traceplots_first_species.pdf
report.html
```

## Difference from Hmsc

Hmsc is the main engine for the full interpretive workflow with traits, phylogeny, random effects, variance partitioning, predictions and many Hmsc-specific outputs. jSDM is provided as a simpler Bayesian JSDM engine focused on latent variables, residual correlations and basic multivariate regression.
