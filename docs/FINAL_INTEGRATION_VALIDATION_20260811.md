# Final Integration Validation Report

Date: 2026-08-11

Project root:

`C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio`

## Validation Commands

All commands were run with:

`C:/Program Files/R/R-4.5.3/bin/Rscript.exe`

## Results

| Check | Result | Evidence |
|---|---:|---|
| app.R parse | PASS | `APP_PARSE_OK` |
| app.R source | PASS | `APP_SOURCE_OK` |
| Shiny strict input/output audit | PASS | 533 UI inputs matched 533 server refs; 87 outputs matched 87 handlers; 10 download buttons matched 10 download handlers; `STRICT_AUDIT_OK` |
| Compare Models testServer | PASS | `COMPARE_TEST_OK` |
| HTTP smoke test | PASS | Local Shiny server responded at `http://127.0.0.1:7860` with page content containing `JSDM Studio`; no R/Rscript process remained after cleanup |
| WebView2 launcher rebuild | PASS WITH WARNING | `dotnet publish` succeeded; .NET emitted `NETSDK1138` because `net5.0-windows` is out of support |
| Inno Setup installer build | PASS | Inno Setup 6.7.1 successful compile; installer copied to `installer/output/JSDMStudio_Setup.exe` |
| Key workflow synthetic suite | PASS WITH EXPECTED MODEL FAILURES | 27 cases ran; 25 `fitted`, 2 `fit_failed`; all cases had `missing_count = 0` and real ZIP outputs |

## Latest Workflow Suite

Summary file:

`output/workflow_synthetic_case_audit_20260811_231232/workflow_synthetic_case_summary.csv`

Status counts:

- Hmsc: 4 fitted, 1 fit_failed
- Hmsc-HPC: 3 fitted
- jSDM: 5 fitted
- GJAM: 3 fitted, 1 fit_failed
- spOccupancy: 4 fitted
- sjSDM: 3 fitted
- boral: 3 fitted

All 27 cases had:

- `missing_count = 0`
- non-empty engine ZIP
- `used_config.yml`
- `diagnostics/engine_status.json`
- `diagnostics/data_check_messages.csv`
- `diagnostics/session_info.txt`
- `standard/run_summary.csv`
- `standard/effects_long.csv`
- `standard/predictions_long.csv`
- `standard/associations_long.csv`
- `standard/fit_metrics.csv`
- `standard/diagnostics_long.csv`
- `standard/output_manifest.csv`

## Expected/Remaining Workflow Risks

- Hmsc GPP quick synthetic case remains `fit_failed` with a recorded diagnostic error: `missing value where TRUE/FALSE needed`. The output contract and ZIP still pass.
- GJAM FC/OC mixed-scale synthetic case remains `fit_failed` with a recorded diagnostic error related to a computationally singular system. The output contract and ZIP still pass.
- These two failures are not hidden as `Completed`; they are intentionally represented as diagnostic `fit_failed` outputs.

## Installer

Installer:

`installer/output/JSDMStudio_Setup.exe`

Size at build time: 136.4 MB.

Important build warning:

- The WebView2 launcher targets `net5.0-windows`, which is no longer supported by Microsoft. The launcher currently builds and is packaged, but future release hardening should retarget a supported Windows framework such as `net8.0-windows`.

## Packaging Note

The final maximum ZIP is generated outside the project directory so that the archive does not include itself. It is expected to include:

- source code,
- `R/`,
- `workflow_scripts/`,
- `examples/`,
- `docs/`,
- `external_packages/`,
- `installer/output/JSDMStudio_Setup.exe`,
- latest workflow output cases under `output/`,
- this validation report.
