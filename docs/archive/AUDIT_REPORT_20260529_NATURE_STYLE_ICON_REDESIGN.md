# JSDM Studio Nature-Style Icon Redesign - 2026-05-29

## Design target

- Reworked the icon system to match a small soft scientific app icon direction:
  - pastel 3D-lite utility icon style
  - rounded square containers
  - mint and teal palette
  - subtle academic dashboard aesthetic
  - minimal scientific symbols that remain legible at small size
  - no mascots, no cartoon scene, no realistic 3D render and no clutter

## Fixed in this round

- Redrew all 19 parameter icons in a unified low-saturation mint/teal system.
- Redrew all 6 engine icons in the same soft scientific app-icon language.
- Redrew the main workflow schematic to match the same visual system.
- Removed SVG filter dependencies after testing showed some renderers could drop filtered content.
- Normalized all SVG files to UTF-8 without BOM.

## Verification

- R 4.5.3 parse check: passed.
- R 4.0.5 fallback source check: passed.
- SVG XML validation:
  - 6 engine SVGs: valid.
  - 19 parameter SVGs: valid.
  - workflow SVG: valid.
- Rendered PNG previews with `rsvg` for:
  - `engine_hmsc.svg`
  - `param_response.svg`
  - `jsdm_workflow.svg`
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
