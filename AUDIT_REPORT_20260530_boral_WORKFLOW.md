# boral Workflow Audit Report - 2026-05-30

## Scope

Audited and repaired the current JSDM Studio boral workflow against boral 2.0.3 source/manual behavior. The audit covered UI inputs, server references, parameter validation, generated reproducible scripts, output folders, ZIP creation, dependency failure handling, and seven synthetic workflow cases.

## Main Fixes

- Replaced the previous scaffold-only boral run path with a real boral adapter in `R/boral_adapter.R`.
- The GUI now exports and runs `reproducible_script/run_this_boral_analysis.R` for every boral run.
- Added `tables/boral_result_workflow_map.csv` listing expected boral outputs and their meaning.
- Fixed boral API names:
  - `lognormal` is normalized to boral's `lnormal`.
  - `power.exponential` is corrected to `powered.exponential`.
- Strengthened validators for:
  - binomial successes with `trial.size`,
  - zero-truncated count families,
  - positive-only `lnormal`, `gamma`, and `exponential`,
  - traits row count,
  - offset matrix dimensions,
  - structured latent variable distance matrices,
  - categorical X/traits conversion.
- Fixed CSV row-name handling for matrix uploads such as `distmat.csv` and `offset.csv`; numeric row indices from `write.csv(row.names=TRUE)` are no longer misread as data columns.
- Fixed a potential `arguments imply differing number of rows` error in dependency diagnostics.
- Added robust dependency failure diagnostics instead of reporting Completed/scaffold:
  - `diagnostics/engine_status.json`
  - `diagnostics/JAGS_status.txt`
  - `diagnostics/boral_dependency_error.txt`
  - `diagnostics/boral_fit_error.txt`
  - `diagnostics/session_info.txt`
- Ensured failed runs still write non-empty output directories, standard comparison tables, a report, executable script, and a real ZIP.
- Added `examples/boral/run_real_example_suite.R`, a seven-case reproducibility and workflow test suite.
- After rerunning with JAGS available, fixed additional real-run issues:
  - skipped self-copying the generated JAGS model file,
  - saved MCMC samples by automatically enabling boral `save.model=TRUE` when `Save MCMC samples` is requested,
  - wrote fitted values from boral list outputs,
  - skipped incompatible newdata predictions with a manifest instead of warnings,
  - wrote unavailable correlation outputs as manifest rows when a model branch lacks the needed coefficients.

## UI / Server Audit

- boral UI IDs found: 68.
- boral server `input$boral_*` references found: 62.
- Server references without UI controls: none after accounting for `synced_numeric_slider()`.
- boral outputs without render/download handlers: none.
- UI-only output IDs are expected display/download outputs.

## Seven Synthetic Cases

All cases were executed with `Rscript examples/boral/run_real_example_suite.R`.

Current rerun result: all seven cases completed with exit code 0, no missing required files, no empty output directories, and real non-empty ZIP files. Six cases performed real boral fitting; the full-upload `do.fit=FALSE` case correctly stopped at `model_defined`.

| Case | Branches covered | Status on this machine |
|---|---|---|
| 01_poisson_lv_covariates | Poisson, XData numeric+factor, latent variables, newdata | `fitted` |
| 02_binomial_trials_xind | Binomial successes, trial.size.csv, X.ind, no latent variables | `fitted` |
| 03_normal_pure_ordination | Normal response, pure latent-variable ordination, no XData | `fitted` |
| 04_negative_binomial_row_offset | Negative binomial, row random effects, row.ids, offset | `fitted` |
| 05_spatial_exponential_distmat | Structured latent variables, exponential distance matrix | `fitted` |
| 06_traits_fourth_corner_ssvs | Traits/fourth-corner setup, categorical traits, SSVS traits branch | `fitted` |
| 07_mixed_family_full_upload_define | Mixed family vector, ranef.ids, distmat, offset, do.fit=FALSE model-definition branch | `model_defined` |

Latest suite summary: `examples/boral/last_real_example_suite_summary.csv`.

Latest verification rerun: 2026-05-30 09:28 local session.

| Check | Result |
|---|---|
| `app.R` parse/source | Passed |
| `R/boral_adapter.R` parse/source | Passed |
| R version | R 4.5.3 |
| JAGS/rjags | `rjags` loaded and linked to JAGS 4.3.2 |
| R2jags | Loaded |
| boral | Loaded, version 2.0.3 |
| Required files missing across 7 cases | None |
| Empty output directories across 7 cases | None |
| ZIP files | All seven ZIPs created and non-empty |
| Engine warnings/errors | None |

## Dependency Status

This machine now loads the required boral stack:

- `rjags`: loads and links to JAGS 4.3.2.
- `R2jags`: available.
- `boral`: available.

Rerun command:

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" examples\boral\run_real_example_suite.R
```

Observed result after dependencies are available: ordinary fit cases became `fitted`; the full-upload `do.fit=FALSE` case became `model_defined`.

## Important Output Files

For any boral run:

- Start with `diagnostics/engine_status.json`.
- If dependency failure: read `diagnostics/boral_dependency_error.txt` and `diagnostics/JAGS_status.txt`.
- If model construction/fitting failure: read `diagnostics/boral_fit_error.txt`.
- Reproducible script: `reproducible_script/run_this_boral_analysis.R`.
- Standard comparison tables: `standard/run_summary.csv`, `standard/effects_long.csv`, `standard/predictions_long.csv`, `standard/associations_long.csv`, `standard/fit_metrics.csv`.
- Human report: `report/boral_report.html`.

## Remaining Risks

- boral 2.0.3 runs one JAGS chain internally for latent-variable identifiability; convergence interpretation differs from engines with multiple chains.
- Some optional plotting/extraction functions are wrapped as safe optional steps because boral output object structure varies by family/model branch.
- Mixed-family plus traits/random/spatial branches should still be treated as advanced after basic single-family tests pass, especially for publication-scale MCMC.
