# JSDM Studio Apple-Style Parameter and Engine Icon Restyle - 2026-05-29

## Fixed in this round

- Replaced the small line-style parameter icons with larger app-icon style semantic icons.
- Increased main parameter icon display size from 22 px to 36 px.
- Increased checkbox parameter icons to 32 px and synchronized number-field icons to 28 px.
- Added rounded gradient tiles, inner highlights, depth shadows and higher-detail SVG drawings for each parameter family.
- Preserved semantic mapping so each parameter icon still matches the control meaning instead of using a generic marker.
- Enlarged Engine Guide and workflow images:
  - Engine banner images increased from 132 px to 172 px.
  - Engine card images increased from 104 px to 132 px.
- Rebuilt the six engine image assets as polished 512-style SVG app icons:
  - `engine_hmsc.svg`: HMSC hierarchy, traits, phylogeny, spatial and association/random-effect motifs.
  - `engine_jsdm.svg`: response matrix, predictors and latent-variable network.
  - `engine_gjam.svg`: mixed response types and inverse-prediction curve.
  - `engine_spoccupancy.svg`: imperfect detection, occupancy eye and spatial survey map.
  - `engine_sjsdm.svg`: scalable neural/covariance network with compute-chip motif.
  - `engine_boral.svg`: ordination axes, latent variables and JAGS/MCMC export motif.

## Verification

- R 4.5.3 parse check: passed.
- R 4.0.5 fallback source check: passed.
- Static UI/server audit:
  - UI inputs detected: 462
  - Server inputs detected: 455
  - Missing UI for server inputs: 0
  - Static UI outputs detected: 43
  - Server outputs detected: 46
  - Missing render/download handler for static UI outputs: 0
  - Unused non-action UI inputs: 0
- `app.R` and `www/` image assets contain only ASCII text.

## Remaining note

- Full app launch under R 4.5.3 still requires installing Shiny and the app dependencies into the R 4.5.3 library. The R 4.5.3 syntax parse is clean, but `shiny` is not available in `C:/Program Files/R/R-4.5.3/library`.
