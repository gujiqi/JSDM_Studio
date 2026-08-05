# JSDM Studio Reference Icon Restyle - 2026-05-29

## Reference check

- Inspected `D:/Otso_Ovaskainen/HmscWorkbench/HmscGUI_progress_eta_version_FIXED.zip`.
- The ZIP does not contain standalone image assets such as PNG/JPG/SVG.
- The useful reference element is the small icon-badge / metric-card visual language in `app.R`.
- No workflow logic, layout, text or code from the reference app was copied into JSDM Studio.

## Updated in this round

- Replaced the previous engine SVG pictures with small polished scientific icon illustrations:
  - `www/engine_hmsc.svg`
  - `www/engine_jsdm.svg`
  - `www/engine_gjam.svg`
  - `www/engine_spoccupancy.svg`
  - `www/engine_sjsdm.svg`
  - `www/engine_boral.svg`
- Adjusted Engine Guide card images so they render as compact icon tiles instead of large rough diagrams.
- Adjusted workflow banners to use the same compact local icon assets.

## Verification

- `app.R` parses successfully.
- `app.R` sources successfully.
- All updated icon assets are local SVG files under `www/`.
- `app.R` and `www/` remain ASCII-only.
