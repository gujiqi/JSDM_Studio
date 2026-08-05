# JSDMWorkbench UI Sync, Preview and Visual Audit - 2026-05-25

## Fixed in this round

- Added a reusable synchronized numeric control in `app.R`.
  - The original Shiny `inputId` is preserved.
  - The left slider and right numeric value box update each other.
  - Server-side code continues reading the same `input$...` values.
- Replaced core numeric controls in all major workflows with synchronized slider/value controls:
  - Project
  - Hmsc
  - jSDM
  - GJAM
  - spOccupancy
  - sjSDM
  - boral
- Added a global `Data Preview` tab.
  - It detects uploaded files from all engine workflows.
  - CSV/TSV files preview the first 100 rows.
  - Tree/text/YAML/JSON files preview the first 100 lines.
  - Binary files such as RDS are listed and protected from unsafe parsing.
- Added a local workflow schematic at `www/jsdm_workflow.svg`.
- Added visual cards and stronger icon badges for a cleaner journal-style overview.
- Fixed a Project recommendation UI bug where the sjSDM score card was incorrectly shown as a duplicate jSDM card.

## Verification

- `app.R` parses successfully.
- `app.R` sources successfully.
- Static UI/server audit:
  - Missing UI for server inputs: 0
  - Missing render/download handler for UI outputs: 0
  - Unused non-action UI inputs: 0
- `app.R` and `www/jsdm_workflow.svg` contain only ASCII text.

## Remaining risks

- Advanced prior fields that deliberately allow `NA` remain plain numeric inputs. This is intentional because sliders cannot represent `NA` safely.
- The synchronized slider ranges are broad defaults. Users can still type exact values in the numeric box within each control's allowed range.
- File preview uses the Shiny temporary upload path, so previews are only available during the active session.
