# JSDM Studio all-workflow full recheck

Date: 2026-05-30

Project root:

`C:/Users/Google/Documents/CodeX/JSDMStudio_AUDITED_FIXED_20260529_SEMANTIC_3DLITE_ICONS/JSDMStudio`

This pass rechecked the Shiny application wiring and executed the real example suites for all modelling engines exposed in the application.

## Environment

- R: 4.5.3 from `C:/Program Files/R/R-4.5.3/bin/Rscript.exe`
- Hmsc: 3.3.7
- jSDM: 0.2.7
- gjam: 2.7
- spOccupancy: 0.8.0
- sjSDM: 1.0.7
- boral: 2.0.3
- rjags: 4.17
- R2jags: 0.8.9
- sjSDM Python stack reported available through reticulate: torch, torch_optimizer, pyro, madgrad

## Code change made in this pass

- `install_packages.R` now installs the HMSC result-script helper packages `colorspace`, `writexl`, and `vioplot`.
- These packages were installed in the local R 4.5 library and HMSC was rerun after installation.
- Reason: HMSC was already fitting successfully, but the generated S1-S7 result scripts can use these packages for full convergence/result outputs.

## Shiny wiring audit

`app.R` parsed and sourced successfully after the installer update.

Global static wiring counts:

| Item | Count | Result |
|---|---:|---|
| UI input/control ids | 472 | OK |
| `input$...` server references | 465 | OK |
| Server references missing from UI | 0 | OK |
| `output$...` ids | 46 | OK |
| render/download handlers | 46 | OK |
| Outputs without render/download handler | 0 | OK |

The only UI ids not directly read as `input$...` are download button ids such as `hmsc_download`, `jsdm_download`, `gjam_download`, `spocc_download`, `sjsdm_download`, `boral_download`, and `compare_download`. This is expected in Shiny because their logic is attached through `output$... <- downloadHandler(...)`.

Per-panel wiring counts:

| Prefix | UI controls | Server input refs | Outputs |
|---|---:|---:|---:|
| hmsc | 103 | 102 | 6 |
| jsdm | 87 | 86 | 6 |
| gjam | 62 | 61 | 6 |
| spocc | 62 | 61 | 6 |
| sjsdm | 73 | 72 | 6 |
| boral | 63 | 62 | 6 |
| compare | 9 | 8 | 4 |
| project | 11 | 11 | 2 |

## Executed workflow suites

All suites were executed with small synthetic or bundled example data so that real fitting and output generation were exercised.

| Engine/suite | Script | Cases | Result | Latest summary |
|---|---|---:|---|---|
| HMSC real examples | `examples/HMSC/run_real_example_suite.R` | 8 | 8 fitted, bad=0 | `output/Hmsc_real_example_suite_summary_20260530_100058.csv` |
| HMSC parameter sweep | `examples/HMSC/run_parameter_sweep.R` | 15 | 15 fitted, bad=0 | `output/Hmsc_parameter_sweep_summary_20260530_100210.csv` |
| jSDM | `examples/jSDM/run_real_example_suite.R` | 7 | 7 fitted, bad=0 | `output/jSDM_real_example_suite_summary_20260530_095028.csv` |
| GJAM | `examples/GJAM/run_real_example_suite.R` | 7 | 7 fitted, bad=0 | `output/GJAM_real_example_suite_summary_20260530_095042.csv` |
| spOccupancy | `examples/spOccupancy/run_real_example_suite.R` | 7 | 7 fitted, bad=0 | `output/spOccupancy_real_example_suite_summary_20260530_095106.csv` |
| sjSDM | `examples/sjSDM/run_real_example_suite.R` | 9 | 9 fitted, bad=0 | `output/sjSDM_real_example_suite_summary_20260530_095350.csv` |
| boral | `examples/boral/run_real_example_suite.R` | 7 | 6 fitted + 1 model_defined, bad=0 | `examples/boral/last_real_example_suite_summary.csv` |

Total executed checks in this all-workflow pass: 60 cases.

## Parameter and branch coverage

HMSC:

- Covered probit, poisson, and normal responses.
- Covered no random level, sample/group random level, spatial Full, spatial NNGP, and spatial GPP.
- Covered traits, Newick phylogeny, taxonomic/correlation phylogeny, offsets, XRRR, advanced random-level units, distance matrix, random-level xData, N-only level, nfMin/nfMax, spatial alpha prior grid, partition column, MCMC controls, output toggles, omega plot controls, and executable S1-S7 output.
- Latest HMSC summaries show `engine_status=fitted`, `exported_script=ok`, and no missing required outputs.

jSDM:

- Covered binomial probit, binomial logit with trial size, poisson log, gaussian, latent/random branches, traits, prediction ids, constrained case, and output switches.

GJAM:

- Covered mixed `typeNames`, CON, PA, CA, DA, FC, CC, OC/CAT style branches, effort, random column, holdout/censoring, conditional prediction, sensitivity, and trait/missing-data trimming.

spOccupancy:

- Covered PGOcc, spPGOcc, msPGOcc, spMsPGOcc, lfMsPGOcc, intPGOcc, and svcPGOcc branches with detection/occupancy inputs, spatial coordinates, prediction, PPC/WAIC-style outputs where configured.

sjSDM:

- Covered binomial, poisson, gaussian, negative-binomial branch, linear environment model, DNN environment model, spatial model, generated spatial eigenvectors, Madgrad optimizer, trait/assembly outputs, CV/tuning branch, predictions, weights, importance, ANOVA/internal-structure outputs.

boral:

- Covered poisson latent-variable covariates, binomial with trial size and X.ind, normal pure ordination, negative-binomial with row offset, spatial exponential distance matrix, traits/fourth-corner/SSVS, and mixed-family define-only branch.
- The define-only branch intentionally returns `model_defined` when `do.fit=FALSE`; it still writes the real output structure and ZIP.

## Output contract check

The suite summaries verify that required outputs were present and non-empty:

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
- real ZIP files for each completed suite case

The application also creates diagnostic files on failure paths, including:

- `diagnostics/engine_status.json`
- `diagnostics/data_check_messages.csv`
- engine-specific `diagnostics/*_error.txt`
- `diagnostics/session_info.txt`

## Remaining risks

- These are small synthetic/example-data runs with quick MCMC or quick fitting settings. They prove software execution, parameter wiring, output generation, and ZIP integrity; they are not publication-level convergence evidence.
- Numeric sliders were not exhaustively enumerated over every possible value. Instead, every visible parameter was checked for UI/server/config wiring, and representative valid values were run through the relevant engine branches.
- HMSC prints internal messages such as `setting updater$GammaEta=FALSE` for no-random-effect or NNGP/GPP cases and `setting updater$Gamma2=FALSE` with phylogeny matrices. These are Hmsc package compatibility messages, not failed JSDM Studio statuses.
- boral production fitting depends on external JAGS. JAGS/rjags/R2jags were available in this run.
- sjSDM depends on the reticulate/PyTorch stack. It was available in this run; GPU/CUDA is optional and environment-dependent.

## Reproduction commands

Run from the project root with R 4.5.3:

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" examples\HMSC\run_real_example_suite.R
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" examples\HMSC\run_parameter_sweep.R
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" examples\jSDM\run_real_example_suite.R
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" examples\GJAM\run_real_example_suite.R
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" examples\spOccupancy\run_real_example_suite.R
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" examples\sjSDM\run_real_example_suite.R
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" examples\boral\run_real_example_suite.R
```

