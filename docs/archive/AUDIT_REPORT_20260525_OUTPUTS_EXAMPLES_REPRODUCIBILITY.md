# Output, Example Data and Reproducibility Follow-Up - 2026-05-25

## Fixed

- Reduced empty output folders by changing directory creation from one global all-engine directory set to core directories plus engine-specific directories.
- Added engine-specific scaffold outputs for jSDM, GJAM, spOccupancy, sjSDM and boral.
- Added non-empty model, prediction, plot, result, table and standard manifests for scaffold-only engines.
- Changed standard comparison CSV files from zero-row placeholders to one-row status-aware scaffold rows.
- Exported executable reproducibility scripts for scaffold-only engines:
  - `reproducible_script/run_this_jSDM_analysis.R`
  - `reproducible_script/run_this_GJAM_analysis.R`
  - `reproducible_script/run_this_spOccupancy_analysis.R`
  - `reproducible_script/run_this_sjSDM_analysis.R`
  - `reproducible_script/run_this_boral_analysis.R`
- Added workflow wrappers in `workflow_scripts/` for non-HMSC engines.
- Added per-engine example datasets:
  - `examples/HMSC/`
  - `examples/jSDM/`
  - `examples/GJAM/`
  - `examples/spOccupancy/`
  - `examples/sjSDM/`
  - `examples/boral/`
- Added a completed executable example:
  - `examples/completed_HMSC_spatial_full/`
- Added developer/reviewer reference documentation:
  - `docs/Developer_Reference_Reproducibility_and_Output_Contract.md`

## Verification

- R parse passed for `app.R`, `R/helpers.R`, `R/helpers_hmsc_fixed.R` and `main.R`.
- Shiny wiring scan passed:
  - missing UI inputs: 0
  - unused non-download UI inputs: 0
  - hidden rendered outputs: 0
  - missing render functions: 0
- Scaffold helper smoke test passed for:
  - jSDM
  - GJAM
  - spOccupancy
  - sjSDM
  - boral
- The scaffold smoke test reported `empty_dirs=0` for all five non-HMSC engines.
- `examples/completed_HMSC_spatial_full/reproducible_script/run_this_HMSC_analysis.R` was executed successfully with Rscript.
- Completed example directory scan reported `empty_completed_dirs=0`.

## Remaining Scope

- HMSC has a real fitting path.
- jSDM, GJAM, spOccupancy, sjSDM and boral still require production adapters for real statistical fitting. Their current runs are explicit `scaffold_only` runs with complete output contracts, examples, diagnostics and executable scaffold reproduction scripts.
