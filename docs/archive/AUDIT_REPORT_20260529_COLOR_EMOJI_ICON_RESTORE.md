# JSDM Studio Color Emoji Icon Restore - 2026-05-29

## Why this round was needed

- The previous redesign converted the reference-style icons into local mint/teal SVG pictograms.
- That made the UI look flatter and more monochrome than the supplied HmscGUI reference.
- The reference package uses colorful system emoji/pictograms directly in `app.R`, so this round restores that behavior.

## Fixed in this round

- Changed visible parameter-label icons from SVG image files to color emoji/pictograms rendered by the browser.
- Kept all emoji in source as HTML numeric entities, so R 4.0.5 on Windows can parse and source `app.R` safely.
- Mapped same-reference parameters to the same reference icon meanings:
  - `Y.csv` / response matrix -> petri dish.
  - `XData.csv` / environmental predictors -> thermometer.
  - species count/list -> sprout.
  - model settings / distribution / family -> DNA.
  - random effects / spatial settings -> map.
  - MCMC / progress controls -> timer.
  - save model -> save disk.
  - output / ZIP / download -> package.
  - diagnostics -> chart.
  - model fit / WAIC / CV -> test tube.
  - reports / parameter tables -> document.
  - variance partitioning -> puzzle.
  - associations / Omega / correlations -> network.
  - predictions / gradients -> crystal ball.
  - checks -> magnifier.
  - run workflow -> rocket.
- Added color emoji icons to command buttons and download links, not only to form labels.
- Rebuilt the six engine SVGs and the workflow SVG so they use the same colorful emoji/pictogram style instead of monochrome line icons.
- Removed stale PNG preview files from the project to avoid showing old flat-icon previews.

## Verification

- R 4.5.3 parse check: passed.
- R 4.0.5 fallback source check: passed.
- SVG XML validation:
  - 6 engine SVGs: valid.
  - workflow SVG: valid.
  - 33 `param_*.svg` fallback assets: valid.
- Static UI/server audit:
  - UI inputs detected: 462
  - Server inputs detected: 455
  - Missing UI for server inputs: 0
  - Static UI outputs detected: 43
  - Server outputs detected: 46
  - Missing render/download handler for static UI outputs: 0
  - Unused non-action UI inputs: 0
- `app.R` and `www/` SVG assets are ASCII-safe because emoji are encoded as HTML numeric entities.

## Remaining note

- Full R 4.5.3 app launch still requires installing `shiny` and the Shiny dependency stack into the R 4.5.3 library. Syntax parsing is clean under R 4.5.3.
