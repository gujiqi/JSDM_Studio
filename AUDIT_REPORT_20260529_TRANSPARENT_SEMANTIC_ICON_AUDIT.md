# JSDM Studio Transparent Semantic Icon Audit - 2026-05-29

## Scope

This pass audits and corrects the icon layer used by parameter labels, command buttons, engine banners, and the workflow illustration.

## Fixed

- Replaced repeated emoji/entity icons in `param_icon_script` with local SVG assets loaded as `<img>` tags.
- Made parameter and command icon containers transparent. The visual icon is now the SVG itself, not a colored tile or emoji badge.
- Reworked all 33 `www/param_*.svg` files with transparent backgrounds and distinct scientific meanings.
- Corrected the response/Y icon from a generic dish-like symbol to an actual species-by-site response matrix.
- Separated previously duplicated meanings:
  - response matrix vs species count
  - traits vs species
  - preview vs detection
  - study grouping vs file list
  - model settings vs latent factors
  - diagnostics vs fit metrics
- Replaced the six engine images with transparent-background vector illustrations for HMSC, jSDM, GJAM, spOccupancy, sjSDM, and boral.
- Replaced the workflow banner emoji artwork with vector workflow artwork.
- Tightened automatic icon assignment for `nParallel`, threads/cores, saved weights, and latent-factor settings.

## Semantic Mapping

- `param_response.svg`: Y/response matrix, site-by-species table.
- `param_predictor.svg`: XData/environmental predictors.
- `param_trait.svg`: species trait table and trait measurements.
- `param_species.svg`: response/species count and species lists.
- `param_group.svg`: studyDesign, grouping, random-effect units, folds.
- `param_spatial.svg`: coordinates, longitude/latitude, NNGP/GPP spatial settings.
- `param_phylogeny.svg`: phylogeny, taxonomy, tree, rho.
- `param_detection.svg`: occupancy/detection survey settings.
- `param_latent.svg`: latent factors, ordination, hidden layers, covariance.
- `param_fit.svg`: model-fit and validation metrics.
- `param_diagnostic.svg`: convergence and diagnostic traces.
- `param_association.svg`: residual/species association network.
- `param_partition.svg`: variance partitioning.

## Verification

- SVG XML check: 40/40 SVG files are well formed.
- SVG raster render check with `rsvg`: 40/40 rendered successfully.
- R 4.5.3 parse: `app.R` parses successfully.
- R 4.0.5 source: `app.R` sources successfully.
- Static UI/server check:
  - UI inputs: 455
  - Server inputs: 455
  - Missing UI for server inputs: 0
  - Unused non-action UI inputs: 0
  - UI outputs/downloads: 43
  - Missing render/download handlers: 0

## Notes

- The app still depends on installed Shiny ecosystem packages for full runtime launch. On this machine, R 4.5.3 can parse the app but does not have the full package set installed; R 4.0.5 was used for full `source()` verification.
- Generated preview PNGs are in `icon_previews/` and are not required by the app at runtime.
