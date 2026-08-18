# JSDM Studio Final Integration Audit - 2026-08-12

## Scope

Final integration validation and packaging for the full JSDM Studio developer bundle.

## Fresh Validation Evidence

- `app.R` parse/source: passed with `APP_PARSE_SOURCE_OK`.
- Strict Shiny binding audit: passed with 533 UI input ids, 533 server input references, 87 UI outputs, 87 output handlers, 10 download buttons and 10 download handlers.
- Download handler static audit: passed with no static fake/empty download handlers detected.
- Compare Models `testServer`: passed with a synthetic Hmsc fitted fixture and a real non-empty ZIP.
- HTTP smoke test: passed on localhost port 7846 with HTTP 200 and an HTML payload containing `JSDM Studio`.
- Key workflow example suite: passed. 28 synthetic workflow cases were run across Hmsc, Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM and boral. All 28 cases finished with `status = fitted`, `missing_count = 0` and non-empty ZIP files.

## Latest Workflow Suite Output

- Output root: `output/workflow_synthetic_case_audit_20260812_112359`
- Summary CSV: `output/workflow_synthetic_case_audit_20260812_112359/workflow_synthetic_case_summary.csv`
- Audit report: `output/workflow_synthetic_case_audit_20260812_112359/workflow_synthetic_case_audit_report.md`

Status by engine:

- Hmsc: 5 fitted cases.
- Hmsc-HPC: 3 fitted cases.
- jSDM: 5 fitted cases.
- GJAM: 5 fitted cases.
- spOccupancy: 4 fitted cases.
- sjSDM: 3 fitted cases.
- boral: 3 fitted cases.

## Build Evidence

- WebView2 launcher build: passed. `JSDMStudioLauncher.exe` and 8 runtime/dependency files were copied from `launcher_build`.
- Inno Setup build: passed.
- Installer path: `installer/output/JSDMStudio_Setup.exe`
- Installer SHA256: `5866CEC80C000F315CEF33DA58FBFB8E4AC3441E732FE1D76065A98EC07045F6`

## Important Fixes Covered By This Validation

- Hmsc GPP/NNGP quick fitting disables unsupported or fragile `GammaEta`/posterior alignment paths where Hmsc can otherwise return NA latent-factor summaries on tiny smoke-test posterior samples.
- GJAM synthetic FC and mixed-scale cases now use scientifically valid GJAM response semantics. `CA` is numeric semi-continuous abundance, and FC is generated with the package's own fractional-composition simulator.

## Remaining Risks

- The workflow suite uses quick synthetic data and very short MCMC/neural settings. It validates software execution, output contracts, ZIP creation and diagnostics, not publication-level convergence.
- Hmsc-HPC is validated in the CPU-supported subset exposed by this GUI; it is not a full GPU/cluster validation.
- GJAM composition and mixed-scale real data can still fail when users provide invalid `typeNames`, singular compositions, too few observations or inconsistent censoring/effort metadata.
- sjSDM, Hmsc-HPC and boral remain sensitive to external Python/PyTorch/JAGS/R package availability on other machines.
