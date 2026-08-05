# JSDM Studio Semantic 3D-lite Icon Correction

Date: 2026-05-29

## Scope

This pass focuses on the parameter-icon audit requested after the previous transparent-icon build.

## Fixed

- Replaced the core parameter icons with transparent-background, colored 3D-lite SVGs.
- Added `www/param_family.svg` for response family / distribution / likelihood parameters.
- Enlarged label and command-button icons in `app.R` so file inputs and parameter controls are easier to read.
- Updated automatic icon matching so important parameter meanings map correctly:
  - `Y.csv`, response matrix -> species-response matrix icon
  - `XData.csv`, predictors, environmental covariates -> environmental thermometer icon
  - species count/list/metadata -> plant/species icon
  - family, distr, distribution, likelihood -> DNA/family icon
  - spatial coordinates and random effects -> map/pin icon
  - MCMC samples, transient, thin, iterations -> stopwatch/trace icon
  - diagnostics, convergence, residual checks -> chart icon
  - model fit / assessment -> laboratory-fit icon
  - output, ZIP, download, export -> package/ZIP icon
- Confirmed no byte-identical duplicate SVG icons remain in `www/`.
- Added a visual contact sheet:
  - `icon_previews/semantic_3dlite_contact_sheet.html`
  - `icon_previews/semantic_3dlite_contact_sheet.png`

## Verification

- Parsed all SVG files as XML: passed.
- Checked byte-identical SVG duplicates: none found.
- Parsed `app.R` with `C:/Program Files/R/R-4.5.3/bin/Rscript.exe`: passed.
- Launched Shiny with R 4.5.3 on `http://127.0.0.1:6894`.
- Captured browser screenshots for the app and icon contact sheet.

## Notes

The engine-level banner images are preserved because they already provide workflow-level context. The corrected parameter-label icons now carry the stricter semantic mapping requested for individual inputs and controls.
