# HMSC complete workflow added

This build adds a new HMSC Workflow panel as the first model workflow.

## Design

The panel is designed to be safe by default:
- distr = probit
- XFormula = ~ .
- traits off
- phylogeny off
- random effect mode = none
- quick-test MCMC
- reproducible scripts exported
- standard comparison tables exported
- real fitting off by default for interface safety

## HMSC output workflow

The downloaded HMSC ZIP follows the S1-S7 result workflow:

1. S1 Define models
2. S2 Fit models
3. S3 Evaluate convergence
4. S4 Compute model fit
5. S5 Show model fit
6. S6 Show parameter estimates
7. S7 Make predictions

## Output folders

- used_config.yml
- inputs/
- workflow_scripts/
- reproducible_script/
- models/
- tables/
- results/
- plots/
- diagnostics/
- predictions/
- standard/
- report/

## Important default

For a first test, upload only:
- Y.csv
- XData.csv

Then keep:
- distr = probit for 0/1 data
- Random effect mode = none
- Use traits = false
- Use phylogeny = false
