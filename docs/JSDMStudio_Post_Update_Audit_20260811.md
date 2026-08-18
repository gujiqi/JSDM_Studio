# JSDM Studio Post-Update Audit

Date: 2026-08-11

## Scope

This audit checked the updated JSDM Studio package after the explanatory-interface update inspired by HMSC-Studio style strengths.

## Checks Passed

- `app.R` parse check: passed.
- `app.R` source/startup construction check: passed.
- All R files under the project tree parse successfully.
- Shiny UI/server static binding:
  - UI input IDs: 535
  - server `input$...` reads: 526
  - missing input IDs: none
  - UI output IDs: 55
  - missing render functions: none
  - download buttons: 9
  - missing `downloadHandler`: none
- Workflow/engine image assets referenced by `app.R`: all present in `www/`.
- Non-ASCII filenames: none found.
- WebView2 launcher rebuild: passed.
- Inno Setup rebuild: passed.

## Issue Found And Fixed

The source tree had been updated, but `installer/output/JSDMStudio_Setup.exe` was still the older 2026-06-03 installer. This meant a user installing from Setup could miss the newest explanatory UI layer. The installer was rebuilt successfully on 2026-08-11.

Updated installer:

- `installer/output/JSDMStudio_Setup.exe`
- SHA256: `D0E357E7BBDBBF8CB8F637B0B08C255F36E380E120A3790C3F8ECBD11C69B383`

## Remaining Risks

- Several legacy documentation files still contain Chinese text. This does not affect Shiny startup or workflow execution, but an international GitHub/software-paper release would benefit from translating or moving those files into a clearly labelled legacy folder.
- The installer still includes historical backup files and audit notes. This keeps the large developer package complete, but it makes the installed tree larger and less tidy than a minimal end-user installer.
- This pass did not rerun all heavy model-fitting examples. It focused on startup, static wiring, assets, packaging and installer correctness after the explanatory-interface update.
