# HMSC Random-Level Follow-Up Fix - 2026-05-25

## Fixed

- Reworked HMSC random-level selection into explicit designs: `none`, `sample`, `spatial_full`, `spatial_nngp`, and `spatial_gpp`.
- Fixed the real HMSC S1 failure `studyDesign columns must be factors` by converting every `studyDesign.csv` column to factor before calling `Hmsc()`.
- Stopped silently falling back from a selected spatial method to a different method. Full, NNGP and GPP now call `HmscRandomLevel()` with method-specific arguments.
- Added coordinate auto-detection for `longitude/latitude`, `lon/lat`, `x/y`, `easting/northing`, or the first two numeric columns.
- Added GPP knot handling: use `sKnot_file` when present, otherwise generate knots with `Hmsc::constructKnots()`.
- Updated generated `workflow_scripts/S1_define_models.R` and `reproducible_script/run_this_HMSC_analysis.R` so downloaded scripts follow the same logic as the GUI.
- Removed mojibake / non-ASCII UI symbols from `app.R` and fixed R parsing under R 4.0.5.
- Updated `R/helpers.R` and `R/helpers_hmsc_fixed.R` so CLI/helper workflows use the same factor and spatial random-level logic.

## Verification

- `app.R` parses successfully with R 4.0.5 using UTF-8.
- `R/helpers.R`, `R/helpers_hmsc_fixed.R`, and `main.R` parse successfully.
- Shiny wiring scan: missing UI inputs = 0; unused non-download inputs = 0; missing rendered outputs = 0.
- Real HMSC S1 construction was tested with the previously failing uploaded data for all three spatial modes:
  - `spatial_full`
  - `spatial_nngp`
  - `spatial_gpp`
- A complete S1-S7 smoke test with minimal MCMC settings finished with `pipeline_status=fitted`.

## Notes

- The user's earlier failure was caused by passing a `studyDesign` data frame where some columns were still character. Hmsc requires all `studyDesign` columns to be factors, even columns that are not the active random level.
- For first tests, use `No random level` or `Spatial Full Gaussian process` on small datasets. Use `Spatial NNGP` for larger spatial data. Use `Spatial GPP / knots` when a knot approximation is desired.
