# JSDM Studio Semantic Parameter Icon Audit - 2026-05-29

## Fixed in this round

- Added a global semantic parameter icon decorator for all Shiny form labels.
- The decorator uses each control's `inputId` plus visible label text to assign a meaningful inline SVG icon.
- It covers ordinary labels, checkbox labels and synchronized slider/number parameter labels.
- It re-runs after Shiny binding, tab/collapse changes and dynamic DOM updates, so late-rendered preview selectors also receive icons.

## Semantic icon mapping

- Response data: `Y.csv`, response matrix, presence/absence, count and response inputs use a matrix icon.
- Predictors: `XData.csv`, environmental predictors, covariates, site data and new-data inputs use an environment/predictor icon.
- Traits: `TrData.csv`, trait formulas and Gamma/trait-mediated settings use a trait icon.
- Grouping/random effects: `studyDesign.csv`, grouping columns, partitions and random-effect units use a grouped-node icon.
- Spatial settings: coordinates, longitude/latitude, NNGP/GPP/knots/distances and spatial covariance parameters use a map icon.
- Phylogeny/taxonomy: phylogeny files, trees, rho and tree plots use a tree icon.
- Formulas/families: formulas, distributions, link functions, optimizers and model modules use a formula curve icon.
- MCMC controls: samples, burn-in, transient, thin, iterations and epochs use an MCMC clock icon.
- Chains and parallelism: chain counts use a link-chain icon; parallel/thread controls use a parallel-arrow icon.
- Priors and supports: priors, hyperparameters, support summaries and regularization controls use a prior-curve icon.
- Predictions, diagnostics and outputs: prediction controls, convergence/model-fit controls and save/export/download/script controls each use their own icon family.
- Engine-specific terms such as detection, latent factors, type names, mixed response types and associations are mapped to dedicated icons instead of a generic marker.

## Verification

- R 4.5.3 parse check: passed.
- R 4.5.3 package availability check: `shiny` is not installed in `C:/Program Files/R/R-4.5.3/library`, so full app startup under R 4.5.3 still requires installing the app dependencies into that R version.
- R 4.0.5 fallback source check: passed.
- Static UI/server audit:
  - UI inputs detected: 462
  - Server inputs detected: 455
  - Missing UI for server inputs: 0
  - Static UI outputs detected: 43
  - Server outputs detected: 46
  - Missing render/download handler for static UI outputs: 0
  - Unused non-action UI inputs: 0
- `app.R` and local `www/` image assets contain only ASCII text.

## Notes

- This change intentionally avoids manually attaching one-off icons to hundreds of controls. The central mapper makes future parameters inherit the correct icon behavior automatically when their `inputId` or label text is semantically named.
- Local engine guide SVG images remain under `www/`; per-parameter icons are inline SVG so they do not add network dependencies.
