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

## 2026-05-30 audit update

The boral workflow is no longer scaffold-only. It now writes and runs an executable reproducible script:

- `reproducible_script/run_this_boral_analysis.R`
- `workflow_scripts/run_boral_workflow.R`

The adapter attempts real `boral::boral()` fitting when JAGS, rjags, R2jags and boral are available. If those dependencies are missing or cannot load, the workflow reports `fit_failed` and writes complete diagnostic files instead of pretending that fitting completed.

Key repaired details:

- `lognormal` is normalized to boral's `lnormal`.
- `power.exponential` is corrected to `powered.exponential`.
- `trial.size` is checked correctly for binomial successes.
- `distmat.csv` and `offset.csv` matrix dimensions are checked after removing CSV row-name columns.
- categorical XData and traits are encoded before fitting.
- dependency failures write `diagnostics/boral_dependency_error.txt`, `diagnostics/JAGS_status.txt`, `diagnostics/engine_status.json`, `diagnostics/session_info.txt`, standard tables and a real ZIP.

Run the automated boral suite:

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" examples\boral\run_real_example_suite.R
```

The suite covers seven small cases: Poisson latent-variable regression, binomial trials and X.ind, normal pure ordination, negative-binomial row effects with offset, structured spatial latent variables with distmat, traits/fourth-corner SSVS, and mixed-family full-upload model definition.

Latest rerun on 2026-05-30 loaded JAGS 4.3.2 through `rjags`: six cases completed as real `fitted` boral runs and the mixed-family `do.fit=FALSE` case completed as `model_defined`. No required files were missing, no output directories were empty, and all seven ZIP files were non-empty.
