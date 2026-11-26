# JSDM Studio Reference ZIP Icon Mapping - 2026-05-29

## Reference source

- Reviewed `D:/Otso_Ovaskainen/HmscWorkbench/HmscGUI_progress_eta_version_FIXED.zip`.
- The reference package does not include standalone icon image files.
- Its icons are embedded in `app.R` as UI emoji/pictogram choices, so the replacement was implemented by matching those reference icon meanings and recreating them as local SVG assets.

## Direct semantic replacements from the reference UI

- `Y.csv` / response matrix: reference `petri dish` meaning -> `param_response.svg`.
- `XData.csv` / environmental predictors: reference `thermometer` meaning -> `param_predictor.svg`.
- species count / species list: reference `sprout` meaning -> `param_species.svg`.
- model settings, distribution and model type: reference `DNA` meaning -> `param_model.svg`.
- random effects and spatial settings: reference `map` meaning -> `param_spatial.svg`.
- MCMC / progress / sampling controls: reference `timer` meaning -> `param_mcmc.svg`.
- output ZIP / download / export: reference `package` meaning -> `param_output.svg`.
- save model: reference `save` meaning -> `param_save.svg`.
- MCMC diagnostics and plots: reference `chart` meaning -> `param_diagnostic.svg`.
- model fit / evaluation / WAIC / CV: reference `test tube` meaning -> `param_fit.svg`.
- parameter estimates / tables / HTML report: reference `document` meaning -> `param_report.svg`.
- variance partitioning: reference `puzzle` meaning -> `param_partition.svg`.
- species associations / Omega / correlations: reference `network` meaning -> `param_association.svg`.
- predictions / gradients: reference `prediction` meaning -> `param_prediction.svg`.
- data preview: reference `eye` meaning -> `param_preview.svg`.
- data check buttons: reference `magnifier` meaning -> `param_check.svg`.
- run workflow buttons: reference `launch/run` meaning -> `param_run.svg`.
- logs and file lists: reference `log/file` meanings -> `param_log.svg`, `param_filelist.svg`.

## Added or changed behavior

- Expanded the JavaScript icon registry beyond broad categories, so more JSDM Studio-specific parameters get a reference-aligned icon instead of a generic setting icon.
- Added command-button icon decoration for Shiny buttons and download links:
  - Check buttons get `param_check.svg`.
  - Run workflow buttons get `param_run.svg`.
  - Download ZIP buttons get `param_output.svg`.
- Kept fallback icons for parameters that are not present in the reference package, but designed them in the same soft scientific mint/teal style.

## Verification

- R 4.5.3 parse check: passed.
- R 4.0.5 fallback source check: passed.
- SVG XML validation:
  - 6 engine SVGs: valid.
  - 33 parameter/command SVGs: valid.
  - workflow SVG: valid.
- Rendered PNG previews with `rsvg` for:
  - `param_response.svg`
  - `param_predictor.svg`
  - `param_model.svg`
  - `param_output.svg`
  - `engine_hmsc.svg`
- Static UI/server audit:
  - UI inputs detected: 462
  - Server inputs detected: 455
  - Missing UI for server inputs: 0
  - Static UI outputs detected: 43
  - Server outputs detected: 46
  - Missing render/download handler for static UI outputs: 0
  - Unused non-action UI inputs: 0
- `app.R` and `www/` SVG assets contain only ASCII text.

## Remaining note

- Full R 4.5.3 app launch still requires installing `shiny` and the Shiny dependency stack into the R 4.5.3 library. Syntax parsing is clean under R 4.5.3.
