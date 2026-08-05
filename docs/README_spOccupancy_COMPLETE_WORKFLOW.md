# spOccupancy complete workflow update

This build adds a separate spOccupancy Workflow to JSDMWorkbench_Project_Compare_Fix.zip.

## What spOccupancy is for

spOccupancy is for single-species, multi-species, integrated, spatial, temporal, latent-factor and spatially varying coefficient occupancy models. Its key distinction is explicit modelling of imperfect detection using replicated detection-nondetection data.

## Added UI sections

- spOccupancy data upload
- spOccupancy data check
- spOccupancy model settings
- spatial / latent-factor / SVC settings
- MCMC settings
- priors, inits and tuning
- validation, prediction and outputs
- run and results

## Key parameters exposed

- model function: PGOcc, spPGOcc, msPGOcc, spMsPGOcc, lfJSDM, sfJSDM, lfMsPGOcc, sfMsPGOcc, intPGOcc, spIntPGOcc, tPGOcc, stPGOcc, tMsPGOcc, stMsPGOcc, SVC models
- occ.formula
- det.formula
- formula
- range.ind
- cov.model
- NNGP
- n.neighbors
- search.type
- n.factors
- svc.cols
- ar1
- x.positive
- n.batch
- batch.length
- n.burn
- n.thin
- n.chains
- accept.rate
- n.report
- n.omp.threads
- verbose
- seed
- updateMCMC
- beta.normal
- alpha.normal
- community priors
- sigma.sq.ig
- phi.unif
- nu.unif
- inits
- tuning
- fix
- ppcOcc
- waicOcc
- k.fold
- predict
- fitted values

## spOccupancy output ZIP structure

- `used_config.yml`
- `inputs/`
- `models/spOccupancy_model.rds`
- `samples/`
- `tables/`
- `predictions/`
- `plots/`
- `diagnostics/`
- `results/README_spOccupancy_results.txt`
- `report/spOccupancy_report.html`

## Safety

This is a SAFE build. It creates the output scaffold and records parameters. Full model fitting should be connected through a production spOccupancy engine adapter.
