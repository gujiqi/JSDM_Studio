# boral complete workflow update

This build adds a separate boral Workflow to JSDMWorkbench_Hmsc_jSDM_GJAM_spOccupancy_sjSDM.zip.

## What boral is for

boral means Bayesian Ordination and Regression AnaLysis. It fits Bayesian models for multivariate ecological data through JAGS. It can fit:

1. independent response GLMs with explanatory variables only,
2. pure latent-variable models for model-based unconstrained ordination,
3. correlated response GLMs with explanatory variables and latent variables for residual correlation.

## Added UI sections

- boral data upload
- boral data check
- boral model settings
- traits, random effects and variable selection
- MCMC and prior settings
- diagnostics, ordination and prediction settings
- output settings
- run and results

## Key parameters exposed

- family and family vector
- num.lv
- lv.control$type
- distmat
- formula.X
- X.ind
- trial.size
- row.eff
- row.ids
- ranef.ids
- offset
- do.fit
- traits
- which.traits
- SSVS: ssvs.index, ssvs.traitsindex, ssvs.g
- save.model
- mcmc.control: n.burnin, n.iteration, n.thin, seed
- prior.control: type, hypparams
- calc.ics
- summary.boral
- lvsplot
- plot.boral
- coefsplot
- ranefsplot
- get.enviro.cor
- get.residual.cor
- calc.varpart
- predict.boral
- fitted.boral
- tidyboral

## boral output ZIP structure

- `used_config.yml`
- `inputs/`
- `models/boral_model.rds`
- `jags/jagsboralmodel.txt`
- `mcmc/`
- `tables/`
- `ordination/`
- `residuals/`
- `random_effects/`
- `variable_selection/`
- `predictions/`
- `plots/`
- `diagnostics/`
- `results/README_boral_results.txt`
- `report/boral_report.html`

## External requirement

boral uses JAGS through R2jags. Users must install JAGS separately from R before real production fitting.
