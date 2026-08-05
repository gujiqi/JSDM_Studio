# HMSC Parameter Sweep Audit - 2026-05-29

## Scope

This pass specifically audited whether the visible HMSC GUI parameters are actually read by the server, mapped to Hmsc/HmscRandomLevel/sampleMcmc/result functions, and runnable in small real examples.

## Code changes

- Added advanced `HmscRandomLevel` support for:
  - unstructured `units`
  - spatial `sData` with Full / NNGP / GPP
  - distance matrix `distMat`
  - covariate-dependent `xData`
  - `N`-only random levels
- Fixed advanced random-level key normalization so mixed-case UI labels such as `distMat`, `xData`, and `N only` are not misread.
- Added real use of `partition_column` in `createPartition()`.
- Added real use of `env.list` for gradient prediction filtering, and saved `species.list`, `trait.list`, `env.list`, and `computeSAIR` metadata to `tables/S7_prediction_settings_*.csv`.
- Added `pool_chains` handling through `poolMcmcChains()` when enabled.
- Hardened `sample_prior = TRUE` output handling so missing ordinary posterior diagnostics do not crash the workflow.
- Hardened generated S3/S5 scripts so quick-test runs with missing or all-NA diagnostics do not crash exported reproducible scripts.
- Wired the HMSC MCMC preset selector so it updates the numeric sample/transient/thin/nChains/nParallel/verbose controls.
- Made `save_model` control the optional `models/hmsc_model_main.rds` export while still preserving the reproducibility RData used by generated scripts.
- Made diagnostics, convergence, parameter-estimate, variance-partitioning, Omega, plotBeta, plotGamma, plotTree, maxOmega, species-name and standard-output switches affect both the direct run and exported scripts.
- Fixed Omega order mapping: UI `alphabetical` is normalized to corrplot's valid `alphabet`, and `FPC` is now available.
- Added `tables/HMSC_parameter_audit.csv` to every run.
- Added `examples/Hmsc/run_parameter_sweep.R`, a small-data sweep covering constructor, random-level, prior, MCMC, phylogeny/trait, prediction, and output parameters.
- Added `HMSC_PARAMETER_COVERAGE_20260529.md`, a visible-parameter coverage map for the HMSC panel.
- Removed the unsupported HMSC `XSelect` UI/config/test path after review; the app no longer exposes a feature that cannot be safely wired to the installed Hmsc interface.
- Quieted optional post-fit Hmsc plotting/gradient helpers so recoverable package-level messages do not masquerade as exported-script failures.
- Fixed S7 prediction scripts to honor the `Environmental gradients` switch and to write `predictions/prediction_manifest.csv` when prediction files are intentionally disabled.
- Updated the time-calibrated phylogeny example to include a time-tree-derived `C` matrix and a small spatial random level, avoiding all-NA quick-test predictions while still exercising time phylogeny, traits, rho/Gamma/Beta, spatial and reproducible-script outputs.

## Parameter Sweep Result

Command:

```r
Rscript examples/Hmsc/run_parameter_sweep.R
```

R executable used:

```text
C:/Program Files/R/R-4.5.3/bin/Rscript.exe
```

Latest summary file:

```text
output/Hmsc_parameter_sweep_summary_20260529_231224.csv
```

Cases that fitted and exported runnable scripts:

- `constructor_offsets`
- `xrrr`
- `sample_priors_partition`
- `spatial_full_alphapw`
- `spatial_nngp`
- `spatial_gpp_knots`
- `advanced_units`
- `advanced_distmat`
- `advanced_xdata`
- `advanced_N`
- `phylogeny_C_traits`
- `phylogeny_newick_outputs`
- `sample_prior_pool`
- `output_and_plot_switches`
- `omega_plot_controls`

Result: 15 HMSC parameter cases fitted real models and exported runnable scripts successfully. All cases had `check_ok = TRUE` and no missing required output files.

## Regression Example Suite

Command:

```r
Rscript examples/Hmsc/run_real_example_suite.R
```

Latest summary file:

```text
output/Hmsc_real_example_suite_summary_20260529_232240.csv
```

All 8 HMSC real examples fitted and exported runnable scripts:

- `probit_linear_categorical`
- `sample_random_trait`
- `poisson_spatial_full`
- `poisson_spatial_nngp`
- `normal_spatial_gpp`
- `normal_time_phylogeny_traits`
- `normal_taxonomy_tree_spatial_traits`
- `normal_linear`

All 8 real examples now have non-empty standard predictions, non-empty required output directories, `engine_status = fitted`, and `exported_script = ok`.

## Remaining HMSC Risk

- Quick-test MCMC settings are only software checks. They are not publication settings and convergence values can be non-finite with tiny chains.
- `sample_prior = TRUE` is a prior predictive / teaching mode; ordinary posterior fit metrics may be unavailable, so the workflow now records warnings and still exports reproducible scripts.
- Hmsc itself may force or ignore some updater settings depending on model structure, especially no-random-effect models and NNGP/GPP spatial levels. These package-level messages are expected and are preserved in logs.
- `updater$Beta`, `updater$Gamma`, and `updater$Omega` are recorded in `used_config.yml` and `tables/HMSC_parameter_audit.csv`; the installed Hmsc 3.3.7 `sampleMcmc()` exposes a single `updater` object rather than separate formal arguments for these three switches, so they are not blindly passed as fake unsupported formals.
- Very short MCMC examples can print package warnings or plotting messages; the workflow status is determined from `diagnostics/engine_status.json`, exported script exit status, and required output-file checks.

## Diagnostics To Check If A Run Fails

- `diagnostics/engine_status.json`
- `diagnostics/HMSC_S1S7_pipeline_status.json`
- `diagnostics/HMSC_S1S7_error.txt`
- `diagnostics/data_check_messages.csv`
- `diagnostics/session_info.txt`
- `tables/HMSC_parameter_audit.csv`
