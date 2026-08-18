# JSDM Studio Maximum Upload Book-Case Audit

Date: 2026-08-18

## Scope

This audit executed the maximum-upload book chapter cases in `examples/00_UPLOAD_READY_book_chapter_cases` for:

- `ch06_plant_traits_whittaker`
- `ch07_deadwood_fungi`
- `ch11_finnish_birds`

Each case was tested against Hmsc, Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM and boral. Mutually exclusive settings were run as separate branches rather than being enabled simultaneously. All engine-specific upload folders were copied into the corresponding run under both `inputs/maximum_upload_files/` and `data/maximum_upload_files/`.

## Code Changes

- Added `workflow_scripts/book_maximum_upload_case_suite.R`.
- Updated `workflow_scripts/workflow_synthetic_case_suite.R` so every completed branch re-runs the output-contract repair after standard synonym tables are created. This prevents optional prediction skips from leaving missing standard tables.
- Updated `workflow_scripts/hmschpc_runner.py` to export Hmsc-HPC `sigma` and `rho` summaries with pandas-compatible index handling.

## Final Maximum Upload Run

- Audit root: `output/book_maximum_upload_audit_20260818_162749`
- Summary CSV: `output/book_maximum_upload_audit_20260818_162749/maximum_upload_case_summary.csv`
- Audit report: `output/book_maximum_upload_audit_20260818_162749/maximum_upload_audit_report.md`
- Audit ZIP: `output/JSDMStudio_book_maximum_upload_audit_20260818_162749.zip`
- Audit ZIP size: 241,567,696 bytes (230.38 MB)
- Audit ZIP SHA256: `4687630DFF7836567D5F2519D50C1CC0E1219047C889D879CF1438046D86719D`

## Status Summary

| Engine | Branches | Status |
|---|---:|---|
| Hmsc | 15 | fitted |
| Hmsc-HPC | 9 | fitted |
| jSDM | 15 | fitted |
| GJAM | 9 | fitted |
| spOccupancy | 12 | fitted |
| sjSDM | 9 | fitted |
| boral | 9 | fitted |

Total branches: 78

Non-fitted branches: 0

Branches with missing required output-contract items: 0

Independent per-branch ZIP files: 78

## Output Contract Checked

Every branch was checked for:

- `used_config.yml`
- `inputs/`
- `data/`
- `models/`
- `tables/`
- `results/`
- `plots/`
- `predictions/`
- `diagnostics/`
- `workflow_scripts/`
- `reproducible_script/`
- `standard/`
- `report/`
- real non-empty ZIP file

## File Count And ZIP Ranges

| Engine | Branches | File count range | ZIP size range (KB) | Branches with non-fatal warnings |
|---|---:|---:|---:|---:|
| Hmsc | 15 | 105-110 | 5918.1-6602.4 | 3 |
| Hmsc-HPC | 9 | 135-147 | 234.3-266.0 | 9 |
| jSDM | 15 | 104-114 | 127.7-345.5 | 12 |
| GJAM | 9 | 135-140 | 327.5-842.8 | 6 |
| spOccupancy | 12 | 103-111 | 180.9-660.0 | 0 |
| sjSDM | 9 | 99-106 | 283.1-797.3 | 5 |
| boral | 9 | 105-105 | 159.0-259.4 | 6 |

## Integration Checks

- `app.R` parse: passed (`APP_PARSE_OK`)
- `app.R` source: passed (`APP_SOURCE_OK`)
- Strict Shiny binding audit: passed
  - UI input IDs: 533
  - server input references: 533
  - refs not in UI: 0
  - UI inputs not read: 0
  - output UI IDs: 87
  - output handlers: 87
  - download buttons: 10
  - download handlers: 10
- Compare Models test: passed (`COMPARE_TEST_OK`)
- HTTP smoke test: passed (`HTTP_SMOKE_OK`)
- WebView2 launcher rebuild: passed
- Inno Setup installer rebuild: passed
  - Installer: `installer/output/JSDMStudio_Setup.exe`
  - Installer size: 143,234,971 bytes (136.60 MB)
  - Installer SHA256: `68F78606472741EAFE3D169B759EA89D31386A7625708378F675A3019230C1EA`

## Remaining Non-Fatal Warnings

These warnings were retained in diagnostics because they describe package/API limitations or quick-test edge cases, not hidden fitting failures:

- Hmsc GPP quick/small posterior runs disabled `alignPost`; one ch06 GPP branch could not compute one model-fit control statistic.
- Hmsc-HPC fitted all CPU branches, but ArviZ diagnostics were skipped for the current xarray object shape. pyhmsc package metadata is also unavailable in the local Python environment even though the import and workflow run succeeded. Some newdata predictions with categorical factors were skipped when the design matrix lacked training-only dummy columns.
- jSDM fitted all branches. Residual correlations are skipped when `n_latent <= 1`, which is expected. Some optional `predict.jSDM` newdata calls report matrix/species mismatches; fitted-value standard predictions remain exported.
- GJAM fitted all branches. Some optional `gjamPredict(newdata=...)` calls fall back to fitted-data prediction when factor levels cannot be represented safely.
- sjSDM fitted all branches. `internalStructure` is skipped for non-spatial models, and ch11 newdata contains an unseen `substrate` level for some optional prediction calls.
- boral fitted all branches. Some optional `fitted.boral` or `predict.boral(newX=...)` calls report package-level newX/newrow ID restrictions; standard outputs and ZIPs are complete.

## Conclusion

The maximum-upload book-case suite is now executable, reproducible and fully packaged. All 78 tested branches reached `fitted` status with complete output-contract folders and real ZIP files. Remaining warnings should be interpreted as engine-specific optional post-processing limits and are recorded in each branch's diagnostics.
