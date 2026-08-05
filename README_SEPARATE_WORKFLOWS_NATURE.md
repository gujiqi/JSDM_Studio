# JSDMWorkbench separate workflows redesign

This version follows the principle requested by the user:

**Separate workflows, shared comparison.**

## Navigation

- Home
- Project
- Engine Guide
- Hmsc Workflow
- jSDM Workflow
- Compare Models
- Guides

## Why this design

Hmsc and jSDM have different data requirements, parameter systems, assumptions and outputs. Therefore:

- Hmsc has its own data upload.
- Hmsc has its own data check.
- Hmsc has its own model settings.
- Hmsc has its own output settings.
- Hmsc has its own run and results.

The same is true for jSDM.

Model comparison happens only after both models have separately created their own output folders.

## Hmsc parameters preserved

The Hmsc Workflow includes detailed settings for:

- distr
- XFormula
- TrFormula
- traits
- phylogeny/taxonomy
- random effect mode
- random effect column
- spatial method
- nNeighbours
- coordinate columns
- seed
- samples
- transient
- thin
- nChains
- nParallel
- verbose
- output modules: predicted values, model fit, CV, WAIC, diagnostics, Beta/Gamma, variance partitioning, Omega, gradient predictions, support thresholds

## jSDM parameters added

The jSDM Workflow includes detailed settings for:

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
- residual correlations
- environmental correlations
- traceplots
- predictions
- model object
- report

## Runtime safety

This build uses a safe workflow scaffold: it creates output folders, copies inputs, saves configs, writes diagnostic status files, creates reports and ZIP files. Production model fitting can be connected through engine adapters.


## Spatial method help correction

The Hmsc spatial-method help now explains Full, NNGP and GPP instead of only mentioning NNGP.
