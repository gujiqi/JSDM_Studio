# Hmsc Audit Report - 2026-05-29

## Scope

Audited and repaired the Hmsc workflow against the supplied S1-S7 reference workflow and the installed Hmsc 3.3.7 API under R 4.5.3 on Windows.

Primary files changed:

- `app.R`
- `examples/Hmsc/run_real_example_suite.R`

## Fixed

- Connected Hmsc random-level priors to real `setPriors.HmscRandomLevel()` calls, including `nfMin`, `nfMax`, `a1`, `b1`, `a2`, `b2`, `alphapw`, and `setDefault`.
- Fixed sample random effects so `studyDesign` columns are forced to ordered factors before `Hmsc()`.
- Split and tested the three spatial random-effect modes: `spatial_full`, `spatial_nngp`, and `spatial_gpp`.
- Added NNGP neighbour clamping so `nNeighbours` cannot exceed the number of spatial units minus one.
- Connected real `sampleMcmc()` options: `initPar`, `alignPost`, `sample_prior/fromPrior`, `GammaEta`, `adaptNf`, `verbose`, `nChains`, and `nParallel`.
- Normalized unsupported `initPar` UI labels so only Hmsc-valid `"fixed effects"` is passed; other labels fall back to `NULL`.
- Added optional Hmsc constructor adapters for `Loff`, `C`, `phyloTree`, `XRRRData`, `XRRRFormula`, `XRRRScale`, `ncRRR`, and `ranLevelsUsed`.
- Added taxonomy-table support: a CSV/TSV with `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, and `species` columns is converted into a positive-definite taxonomy correlation matrix `C` for Hmsc.
- Removed unsupported `XSelect` from the active GUI/config/test path after review, so the app no longer exposes an unimplemented Hmsc selector.
- Made S3/S4 robust for quick-test edge cases: non-finite PSRF values and metric functions such as AUC/Tjur with no control observations are recorded as diagnostics instead of aborting a fitted model.
- Respected the WAIC output switch; WAIC is no longer computed when `computeWAIC` is off.
- Added real `standard/` outputs from fitted Hmsc results:
  - `standard/run_summary.csv`
  - `standard/fit_metrics.csv`
  - `standard/effects_long.csv`
  - `standard/predictions_long.csv`
  - `standard/associations_long.csv`
  - `standard/output_manifest.csv`
- Filled `plots/` with copied figure PDFs and `plots/plot_manifest.csv`.
- Filled `report/` with `report/Hmsc_report.html`.
- Added required diagnostics:
  - `diagnostics/engine_status.json`
  - `diagnostics/session_info.txt`
  - `diagnostics/HMSC_S1S7_pipeline_status.json`
  - `diagnostics/HMSC_S1S7_error.txt` on failure.
- Exported executable scripts now include real S1-S7 steps plus `S8_standardize_outputs.R`; the master runner executes all scripts and regenerates plots, standard tables, report, and manifest.
- Added a real multi-case Hmsc example suite covering categorical predictors, traits, sample random effects, Full spatial, NNGP spatial, GPP spatial, probit, poisson, normal responses, a time-calibrated Newick phylogeny, and a taxonomy-only classification table.

## Verified Cases

Final suite summary:

`output/Hmsc_real_example_suite_summary_20260529_232240.csv`

| case_id | family | random mode | traits | phylogeny/taxonomy | exported script | status |
|---|---:|---:|---:|---:|---:|---:|
| probit_linear_categorical | probit | none | no | none | ok | fitted |
| sample_random_trait | probit | sample | yes | none | ok | fitted |
| poisson_spatial_full | poisson | spatial_full | no | none | ok | fitted |
| poisson_spatial_nngp | poisson | spatial_nngp | no | none | ok | fitted |
| normal_spatial_gpp | normal | spatial_gpp | no | none | ok | fitted |
| normal_time_phylogeny_traits | normal | spatial_full | yes | time-calibrated Newick tree plus derived `C` matrix | ok | fitted |
| normal_taxonomy_tree_spatial_traits | normal | spatial_full | yes | taxonomy CSV converted to `C` | ok | fitted |
| normal_linear | normal | none | no | none | ok | fitted |

All eight cases produced non-empty ZIP files, no missing required files, non-empty standard predictions, and no empty required output directories.

## Reproducibility Check

The exported script was run successfully for all eight cases, including the phylogeny and taxonomy cases:

```r
Rscript output/JSDMStudio_HMSC_SUITE_normal_time_phylogeny_traits_Hmsc_20260529_232226/reproducible_script/run_this_HMSC_analysis.R
Rscript output/JSDMStudio_HMSC_SUITE_normal_taxonomy_tree_spatial_traits_Hmsc_20260529_232231/reproducible_script/run_this_HMSC_analysis.R
```

They refit the models and regenerated the S1-S7 outputs plus standardized outputs.

## Remaining Risks

- The quick-test examples use very small MCMC settings. Hmsc correctly warns about unstable toy fits; publication analyses must increase samples, transient, thin, and chains, then inspect convergence.
- Hmsc itself disables `GammaEta` for no-random-level models and for GPP/NNGP, which is expected Hmsc package behavior.
- `XSelect` is no longer exposed. It should only be reintroduced if the GUI gains a structured Hmsc `XSelect` object upload/constructor and matching tests.
- Phylogeny support now reads `C` matrices, Newick trees, and taxonomy tables, but scientific validity still depends on matching species names and using a biologically meaningful tree or classification.

## Minimal Tests

Run all Hmsc smoke tests:

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" "examples\Hmsc\run_real_example_suite.R"
```

If a case fails, inspect:

- `diagnostics/engine_status.json`
- `diagnostics/HMSC_S1S7_pipeline_status.json`
- `diagnostics/HMSC_S1S7_error.txt`
- `diagnostics/data_check_messages.csv`
- `diagnostics/session_info.txt`
