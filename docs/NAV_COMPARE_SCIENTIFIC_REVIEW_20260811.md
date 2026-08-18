# Navigation And Compare Models Scientific Review

Date: 2026-08-11

Scope: Home, Project, Engine Guide, Compare Models and Guides text/logic.

## Reviewer Standard Applied

JSDM Studio should describe engines as separate statistical workflows, not interchangeable menu items. A successful software run is not the same as an ecologically adequate analysis. The interface must distinguish:

- `fitted`: fitted outputs and standard summaries exist, but convergence, identifiability and ecological adequacy still require review.
- `model_defined`: a model boundary, compiled representation or executable script exists, but fitted posterior/performance evidence is not claimed.
- `check_failed`: data or settings failed before fitting.
- `fit_failed`: fitting or post-processing failed after the run started.

## Changes Made

- Project recommendation logic now treats replicated detection-nondetection occupancy surveys as a first-class response structure and makes spOccupancy primary only when the occupancy design is present.
- Project metadata text now distinguishes spOccupancy, GJAM and single-family JSDM-family use cases.
- Engine Guide now avoids treating `fitted` as proof of scientific adequacy.
- sjSDM text now avoids implying direct equivalence with Hmsc-style random-level Omega.
- Compare Models now exports explicit columns for:
  - effect-direction comparison eligibility,
  - prediction comparison eligibility,
  - association numeric comparison eligibility,
  - association pattern comparison eligibility,
  - comparable outputs,
  - non-comparable outputs,
  - primary diagnostics.
- Compare Models matrix now separates directly comparable, conditionally comparable and not directly comparable outputs.
- Guides now explicitly separate direct comparisons from conditional comparisons and non-comparable raw parameter comparisons.

## Comparability Rules

Directly comparable:

- workflow status,
- output contract completeness,
- ZIP existence,
- diagnostic file presence,
- report/script presence,
- runtime as practical cost on the same machine.

Conditionally comparable:

- prediction metrics, only under matched validation units, response scale and metric definition;
- predictor-to-species effect direction, only after checking predictor names, scaling, contrasts, link functions and response family;
- broad association patterns, only with engine-specific `association_type`, scale and notes.

Not directly comparable:

- Hmsc Omega as a raw numeric equivalent to Hmsc-HPC Eta/Lambda,
- jSDM residual correlations as raw equivalents to Hmsc or boral associations,
- GJAM corMu/sigMu as raw equivalents to link-scale JSDM coefficients,
- spOccupancy detection effects as occurrence/environment effects,
- sjSDM bioticStruct covariance/correlation as a raw equivalent to MCMC latent covariance,
- boral latent residual correlations as raw equivalents to Hmsc Omega.

## Remaining Reviewer Cautions

- The Project recommendation remains a guide, not a statistical decision rule.
- Prediction metric comparability still depends on the actual validation design exported by each workflow.
- Raw model evidence criteria such as WAIC, DIC, AUC and RMSE should not be compared across incompatible response scales.
