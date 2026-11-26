# JSDM Studio Reference-Based Icon Redraw - 2026-05-29

## Reference reviewed

- Reviewed `D:/Otso_Ovaskainen/HmscWorkbench/HmscGUI_progress_eta_version_FIXED.zip`.
- The reference package does not contain separate bitmap or SVG image assets.
- Its useful visual direction is the UI language: large immediately recognizable icons, soft rounded cards, restrained teal/blue scientific colors, and simple pictograms instead of dense technical line art.

## Fixed in this round

- Replaced inline parameter SVG line drawings with real local SVG image assets under `www/`.
- Added 19 independent parameter icon images:
  - `param_response.svg`
  - `param_predictor.svg`
  - `param_trait.svg`
  - `param_group.svg`
  - `param_spatial.svg`
  - `param_phylogeny.svg`
  - `param_formula.svg`
  - `param_mcmc.svg`
  - `param_chain.svg`
  - `param_parallel.svg`
  - `param_seed.svg`
  - `param_prior.svg`
  - `param_prediction.svg`
  - `param_diagnostic.svg`
  - `param_output.svg`
  - `param_mixed.svg`
  - `param_detection.svg`
  - `param_latent.svg`
  - `param_setting.svg`
- Updated the parameter icon decorator so labels now render `<img>` icons from `www/param_*.svg`.
- Increased parameter icon size to 44 px for normal controls, 38 px for checkboxes and 34 px for synchronized number fields.
- Redrew all six engine images in the same larger, cleaner app-icon style:
  - HMSC: hierarchy, traits, phylogeny and spatial/random-effect cues.
  - jSDM: response matrix plus latent-variable network.
  - GJAM: mixed response types and response-scale curve.
  - spOccupancy: detection eye plus spatial survey map.
  - sjSDM: neural/covariance compute motif.
  - boral: ordination axes, latent curve and export/MCMC cue.

## Verification

- R 4.5.3 parse check: passed.
- R 4.0.5 fallback source check: passed.
- SVG XML validation:
  - 6 engine SVGs: valid.
  - 19 parameter SVGs: valid.
  - workflow SVG: valid.
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

- Full R 4.5.3 app launch still requires installing `shiny` and the rest of the Shiny dependency stack into the R 4.5.3 library. Syntax parsing is clean under R 4.5.3.
