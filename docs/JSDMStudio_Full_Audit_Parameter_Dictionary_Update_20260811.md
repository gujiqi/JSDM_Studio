# JSDM Studio Full Audit and Parameter Dictionary Update

Date: 2026-08-11

## Scope

This audit focused on release-blocking Shiny errors and reviewer-facing parameter documentation in the current full JSDM Studio package.

Checked areas:

- `app.R` parse/source startup
- Shiny `inputId`, `input$...`, `output$...`, `render*` and `downloadHandler` bindings
- Runtime-facing helper messages in `R/helpers.R` and `R/helpers_hmsc_fixed.R`
- Existence and UI integration of a complete parameter dictionary
- Real downloadable parameter documentation for review and reproducibility

## Fixes Applied

- Added `docs/JSDMStudio_Parameter_Dictionary.csv`, indexing 536 exposed Shiny controls by workflow, `input_id`, control type, semantic role, pre-run check, common mistake, output connection and reviewer note.
- Added a new Guides accordion panel, **Complete parameter dictionary**, with searchable table preview and a real CSV download handler.
- Cleaned legacy HMSC helper messages that previously contained broken or machine-translated English.
- Improved HMSC progress and report text so exported diagnostics/reports use clear status language.
- Preserved all model workflows and output contracts; no fitting adapter was replaced or bypassed.

## Current Verification Targets

Expected after this update:

- `app.R` parses and sources without error.
- Every UI output has a corresponding render function.
- Every download button has a corresponding `downloadHandler`.
- The parameter dictionary panel appears in Guides.
- `parameter_dictionary_download` writes a real CSV, not a placeholder.

## Remaining Risks

- Some archived legacy backup files and historical audit reports are still included for provenance and may contain old wording. They are not used by the running Shiny app.
- Full ecological validation still depends on optional engine dependencies such as Hmsc, pyhmsc/Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM/PyTorch and boral/JAGS.
- Quick-test model outputs are software smoke tests, not publication-grade inference.
