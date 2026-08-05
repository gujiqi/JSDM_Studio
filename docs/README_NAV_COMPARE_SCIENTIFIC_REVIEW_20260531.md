# Navigation and Compare Models Scientific Review - 2026-05-31

Scope: Home, Project, Engine Guide, Compare Models and Guides in `app.R`.

## Review Standard

The interface should read like ecological-statistical software documentation, not marketing copy. It must not imply that different JSDM engines share identical parameters, assumptions, residual association objects, validation metrics or output scales.

## Changes Made

- Home now describes JSDM Studio as a reviewer-oriented workflow system instead of making a publication-quality claim.
- Project recommendation text now states that project settings are guidance only and never replace data checks, dependency checks or convergence diagnostics.
- Project scoring no longer promotes spOccupancy strongly from spatial coordinates or generic prediction alone. spOccupancy is now strongly recommended only for occupancy, detection, false-absence, imperfect-detection or spatial-occupancy questions.
- Project recommendation output now includes a `Reviewer_caution` column for each engine.
- Engine Guide wording was tightened:
  - Hmsc is described as the full Hmsc-R workflow in this app, not a universal primary model.
  - Hmsc-HPC is described as a CPU pyhmsc/HDF5 subset, not a full Hmsc-R replacement.
  - GJAM is described as a mixed-scale/median-zero observation-model engine with typeNames dependence.
  - spOccupancy is described as an occupancy/imperfect-detection engine, not a generic community matrix model.
  - sjSDM dependency and scope limits are explicit.
  - boral system JAGS and MCMC diagnostics are explicit.
- Compare Models now exports columns that clarify:
  - `Metric_comparison_allowed`
  - `Metric_comparison_note`
  - `Association_comparison_note`
  - `Failure_diagnostic_file`
  - `Data_check_messages`
  - `Session_info`
- `Status_class` now distinguishes `check_failed_before_fit` from `fit_failed_after_start`.
- Compare Models download filename now uses `JSDMStudio_model_comparison_*`.

## Status Definitions

- `fitted`: model/posterior object and standard summaries were produced. This does not by itself prove scientific adequacy.
- `model_defined`: model boundary, scripts or compiled representation exists, but no posterior or performance evidence is claimed.
- `check_failed`: data/settings failed before fitting; inspect `diagnostics/data_check_messages.csv` and `used_config.yml`.
- `fit_failed`: fitting or post-processing failed after the run started; inspect `diagnostics/engine_status.json`, `diagnostics/session_info.txt` and engine-specific error logs.

## Directly Comparable

- Workflow status and output completeness.
- ZIP, script, report and diagnostic file presence.
- Runtime and output-file counts as software/reproducibility indicators.
- Prediction metrics only under matched validation units, response scale, response family and metric definition.

## Not Directly Comparable

- Raw coefficients across engines when scaling, links, priors or encodings differ.
- Hmsc Omega, Hmsc-HPC Eta/Lambda summaries, jSDM residual correlations, GJAM corMu/sigMu, spOccupancy latent factors, sjSDM covariance/correlation and boral residual correlations.
- WAIC/DIC/AUC/RMSE-like metrics unless they target the same response scale and validation design.

