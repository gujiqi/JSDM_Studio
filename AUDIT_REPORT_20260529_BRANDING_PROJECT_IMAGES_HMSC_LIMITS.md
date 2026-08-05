# JSDM Studio Branding, Project, Images and HMSC Limit Audit - 2026-05-29

## Fixed in this round

- Changed the top-left application title from `JSDMWorkbench` to `JSDM Studio`.
- Updated visible home/guide wording to use `JSDM Studio`.
- Simplified the Project page:
  - Removed visible `Project name` and `Main scientific question` fields.
  - Kept hidden defaults so output folder naming and exported configs still work.
  - Expanded `Question template` choices so model recommendation is better matched to real use cases.
- Increased Project study-size slider limits:
  - `n observations / sites` max: `100000000`
  - `S responses / species` max: `100000000`
- Added Engine Guide images to each engine card:
  - Hmsc, jSDM, GJAM, spOccupancy, sjSDM and boral.
- Added workflow upload images inside each model's Step 1 data-upload panel.
- Kept per-parameter visual cues by adding small SVG-style icons to labels and checkbox labels.
- Increased HMSC parameter maxima:
  - `nfMin` max: `1000`
  - `nfMax` max: `10000`
  - `nChains` max: `100`
  - `nParallel` max: `256`
  - `verbose` max: `100000`
  - `transient` max remains `1000000`
  - `thin` max remains `100000`
- Changed synchronized slider min/max display to plain numbers with commas, avoiding scientific notation such as `1e+06`.

## Verification

- `app.R` parses successfully.
- `app.R` sources successfully.
- Static UI/server audit:
  - Missing UI for server inputs: 0
  - Missing render/download handler for static UI outputs: 0
  - Unused non-action UI inputs: 0
- `app.R` and `www/` SVG assets contain only ASCII text.

## Notes

- Hidden `project_name` and `project_question` are retained intentionally because existing run/export logic writes them into configs and output folders.
- The image assets are local SVG figures under `www/`, so the app does not depend on internet image loading.
