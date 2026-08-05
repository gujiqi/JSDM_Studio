# Hmsc complete-parameters update

This version expands the Hmsc Workflow so that more official Hmsc settings are visible and explained.

## Sources used

- Uploaded Hmsc package manual
- Uploaded Hmsc workflow scripts S1-S7
- Uploaded Hmsc vignettes and JSDM book excerpts

## What was added to Hmsc Workflow

### Hmsc constructor-level options

- XScale
- TrScale
- YScale
- truncateNumberOfFactors
- Loff offset file
- ranLevelsUsed
- C phylogenetic correlation matrix file
- XRRR reduced-rank regression
- ncRRR
- XRRRFormula
- XRRRScale
- XRRR file

### HmscRandomLevel details

- random_level_type
- sMethod: Full / NNGP / GPP
- N
- longlat
- units column
- distMat file
- xData file
- sKnot file
- nfMin
- nfMax

### Priors and advanced fitting

- setDefault priors
- a1, b1, a2, b2
- alphapw grid
- sample_prior
- pool_chains

### sampleMcmc/updater settings

- initPar
- alignPost
- updater$GammaEta
- updater$Beta
- updater$Gamma
- updater$Omega

### Output and diagnostic settings

- convergence for Beta, Gamma, Omega, rho, alpha
- maxOmega
- effective sample size
- Gelman PSRF
- var.part.order.explained
- var.part.order.raw
- show.sp.names.beta
- plotTree
- omega.order
- show.sp.names.omega
- plotBeta
- plotGamma
- species.list
- trait.list
- env.list
- nfolds
- partition column
- computeSAIR

## Hmsc downloaded result ZIP

The Hmsc result ZIP is designed to follow the script workflow:

- S1_define_models_template.R
- S2_fit_models.R
- S3_evaluate_convergence.R
- S4_compute_model_fit.R
- S5_show_model_fit.R
- S6_show_parameter_estimates.R
- S7_make_predictions.R

Expected production outputs:

- `used_config.yml`
- `inputs/`
- `models/unfitted_models.RData`
- `models/models_thin_[thin]_samples_[samples]_chains_[chains].Rdata`
- `results/MCMC_convergence.pdf`
- `results/MCMC_convergence.txt`
- `models/MF_thin_[thin]_samples_[samples]_chains_[chains]_nfolds_[nfolds].Rdata`
- `results/model_fit.pdf`
- `results/parameter_estimates.pdf`
- `results/parameter_estimates.txt`
- `results/parameter_estimates_[parameter].csv`
- `results/predictions.pdf`
- `tables/Hmsc_result_workflow_map.csv`
- `diagnostics/engine_status.json`
- `report/Hmsc_report.html`

## Safety

This is still a SAFE build. It exposes settings and writes auditable output scaffolds. Full production fitting should be connected through the Hmsc engine adapter.
