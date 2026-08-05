# Project and Compare Models fix

This build makes the Project page functional and updates Compare Models to include GJAM.

## Project page changes

The Project page now records:

- n observations / sites
- S responses / species
- Q predictors
- response structure
- traits
- phylogeny/taxonomy
- spatial coordinates
- median-zero / many zeros
- prediction/inverse prediction importance

These fields are saved into each engine's `used_config.yml` and are used for a simple engine recommendation table.

## Compare Models changes

Compare Models now supports:

- Hmsc output folder
- jSDM output folder
- GJAM output folder
- selectable engine list
- status/config/report/workflow-map checks
- table/plot/diagnostic file counts
- downloadable comparison CSV
- comparable-output matrix for Hmsc vs jSDM vs GJAM

## Important interpretation

Comparison is at the workflow/output level. Do not equate Hmsc Omega, jSDM residual correlations, and GJAM corMu/sigMu directly.
