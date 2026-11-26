# GJAM complete workflow update

This build adds an independent GJAM Workflow to JSDMWorkbench.

## GJAM purpose

GJAM is used for generalized joint attribute modelling of multivariate ecological responses that can mix data types:
presence-absence, continuous, continuous abundance with zero censoring, discrete abundance/count, fractional composition, count composition, ordinal counts and categorical classes.

## Added UI sections

- GJAM data upload
- GJAM data check
- GJAM model settings
- GJAM response type and composition settings
- GJAM prior and censoring settings
- GJAM prediction, sensitivity and ordination settings
- GJAM output settings
- GJAM run and results

## Key GJAM parameters exposed

- formula
- typeNames
- ng
- burnin
- holdoutN
- holdoutIndex
- censor
- effort
- FULL
- notStandard
- reductList$N
- reductList$r
- random
- REDUCT
- FCgroups
- CCgroups
- PREDICTX
- ematAlpha
- traitList inputs
- priorTemplate inputs
- gjamPredict options
- gjamSensitivity
- gjamOrdination
- gjamConditionalParameters
- gjamIIE
- missing X/Y prediction
- inverse prediction

## GJAM output ZIP structure

- `used_config.yml`: exact settings used.
- `inputs/`: copied inputs.
- `models/gjam_model.rds`: fitted GJAM object when production fitting is connected.
- `chains/`: MCMC chain components.
- `tables/`: parameter summaries, fit diagnostics, sensitivity, correlations and output maps.
- `predictions/`: predicted Y, richness and inverse-predicted X.
- `plots/`: gjamPlot, sensitivity, ordination and IIE figures.
- `diagnostics/`: data checks and engine status.
- `results/README_GJAM_results.txt`: human-readable output explanation.
- `report/GJAM_report.html`: HTML report.

## Safety

This is a SAFE build. It creates the output structure and records all parameters. Full production GJAM fitting should be connected through a GJAM engine adapter.
