# JSDMWorkbench Workflow Preview, HMSC Random Effects and Visual Audit - 2026-05-26

## Fixed in this round

- Updated HMSC MCMC slider limits:
  - `hmsc_transient` maximum is now `1000000`.
  - `hmsc_thin` maximum is now `100000`.
- Reworked the HMSC random-effect UI into a stepwise design:
  - Step 1: choose one random-level design.
  - Step 2: only settings relevant to the selected design are displayed.
  - Spatial Full, Spatial NNGP and Spatial GPP are now visually separated.
  - Advanced compatibility settings are grouped separately so they do not interfere with the default workflow.
- Removed the top-level `Data Preview` tab.
- Added per-workflow upload previews inside each engine's data-check step:
  - HMSC
  - jSDM
  - GJAM
  - spOccupancy
  - sjSDM
  - boral
- Added one local workflow image per engine in `www/`:
  - `engine_hmsc.svg`
  - `engine_jsdm.svg`
  - `engine_gjam.svg`
  - `engine_spoccupancy.svg`
  - `engine_sjsdm.svg`
  - `engine_boral.svg`
- Added small visual icons to parameter labels and checkbox labels so each parameter area has a clear visual cue.

## Verification

- `app.R` parses successfully.
- `app.R` sources successfully.
- Static UI/server audit after excluding dynamic per-engine preview selectors:
  - Missing UI for server inputs: 0
  - Missing render/download handlers for static UI outputs: 0
  - Unused non-action UI inputs: 0
- `app.R` and all `www/*.svg` assets contain only ASCII text.

## Notes

- The per-engine preview controls are registered dynamically by `register_engine_file_preview()`.
- CSV/TSV previews show up to 100 rows. Text/tree/YAML/JSON previews show the first 100 lines. RDS and other binary files are intentionally not parsed in the browser preview.
- Advanced HMSC random-level controls remain available for custom or legacy configurations, but the ordinary workflow should use the main mode selector.
