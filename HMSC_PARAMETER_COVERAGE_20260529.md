# HMSC Visible Parameter Coverage - 2026-05-29

This file maps the visible HMSC workflow controls to the server config, runtime path, and test coverage. "Covered" means the control is read into `hmsc_config()`, recorded in `used_config.yml` and `tables/HMSC_parameter_audit.csv`, and exercised by at least one small HMSC run unless marked as package-limited.

Latest sweep: `output/Hmsc_parameter_sweep_summary_20260529_231224.csv`.

## Result

- 15 HMSC parameter-sweep cases were run with `C:/Program Files/R/R-4.5.3/bin/Rscript.exe`.
- 15 cases fitted real HMSC models and exported runnable reproducible scripts.
- All cases have `check_ok = TRUE`, expected engine status, exported-script status, and no missing required output files.
- The 8-case real example suite also passed: `output/Hmsc_real_example_suite_summary_20260529_232240.csv`.
- All real-suite cases have non-empty standard predictions and non-empty required output directories.

## Data Upload And Preview

| Control | inputId | Coverage | Notes |
|---|---|---|---|
| Y.csv response matrix | `hmsc_Y_file` | Covered | Tested by all fitted cases; coerced to matrix and response-family checked. |
| XData.csv environmental predictors | `hmsc_X_file` | Covered | Tested by all fitted cases; character columns are factorized and numeric-looking columns are numeric. |
| TrData.csv traits | `hmsc_Tr_file` | Covered | Tested by `phylogeny_C_traits`, `sample_random_trait`, `normal_time_phylogeny_traits`, `normal_taxonomy_tree_spatial_traits`. |
| studyDesign.csv grouping design | `hmsc_study_file` | Covered | Tested by sample and advanced random-level cases; grouping columns are converted to factors. |
| coordinates.csv spatial coordinates | `hmsc_coord_file` | Covered | Tested by spatial Full, NNGP and GPP cases. |
| phylogeny / taxonomy file | `hmsc_phylo_file` | Covered | Tested with Newick and taxonomy/correlation workflows. |
| HMSC uploaded data preview | `hmsc_file_preview_choice`, `hmsc_file_preview_meta`, `hmsc_file_preview_table` | Covered | Registered by the engine-specific preview renderer; previews CSV/TSV/text/tree uploads in step 2. |
| Check HMSC data | `hmsc_check` | Covered | Uses `validate_hmsc()` before fitting and writes check messages in diagnostics. |

## Core Model

| Control | inputId | Coverage | Notes |
|---|---|---|---|
| Response distribution | `hmsc_distr` | Covered | Tested with `probit`, `poisson`, and `normal` in real examples; invalid family/data combinations are stopped. |
| Environmental formula | `hmsc_XFormula` | Covered | Tested with `~ .` and explicit formulas; formula variables must exist in XData. |
| Trait formula | `hmsc_TrFormula` | Covered | Tested when traits are enabled; formula variables must exist in TrData. |
| Use TrData / traits | `hmsc_use_traits` | Covered | Controls whether TrData enters the Hmsc constructor. |
| Use phylogeny / taxonomy | `hmsc_use_phylogeny` | Covered | Tested with `C_file`, Newick tree, and taxonomy-derived correlation. |
| Save model object | `hmsc_save_model` | Covered | `FALSE` tested by `output_and_plot_switches`; reproducibility RData is still exported. |
| MCMC preset | `hmsc_preset` | Covered | Now updates samples, transient, thin, nChains, nParallel and verbose; custom values remain editable. |

## Random Effects And Spatial Settings

| Control | inputId | Coverage | Notes |
|---|---|---|---|
| Random level design | `hmsc_random_mode` | Covered | Tested: none, sample, spatial_full, spatial_nngp, spatial_gpp. |
| Grouping column | `hmsc_random_effect_column` | Covered | Tested with `plot`; missing columns fail data check instead of silently using wrong grouping. |
| Advanced units mirror | `hmsc_units_column` | Covered | Tested by unstructured, distMat, xData and N-only cases. |
| Longitude / x column | `hmsc_lon_col` | Covered | Tested with `x`. |
| Latitude / y column | `hmsc_lat_col` | Covered | Tested with `y`. |
| longlat | `hmsc_longlat` | Covered | Tested in spatial Full setup; passed to spatial random-level construction. |
| NNGP nNeighbours | `hmsc_nNeighbours` | Covered | Tested by `spatial_nngp`; slider and numeric field are synced. |
| GPP knots file | `hmsc_sKnot_file` | Covered | Tested by `spatial_gpp_knots`. |
| Legacy spatial method override | `hmsc_spatial_method` | Covered | Kept for older configs; normalized with the selected random mode. |
| Advanced HmscRandomLevel type | `hmsc_random_level_type` | Covered | Tested: none, unstructured units, spatial sData through spatial modes, distance matrix, covariate-dependent xData, N only. |
| Advanced sMethod | `hmsc_sMethod` | Covered | Used for advanced spatial random-level construction. |
| N for N-only level | `hmsc_random_N` | Covered | Tested by `advanced_N`; slider and numeric field are synced. |
| distMat file name | `hmsc_distMat_file` | Covered | Tested by `advanced_distmat`; matrix is aligned to unit levels. |
| random-level xData file name | `hmsc_xData_file` | Covered | Tested by `advanced_xdata`; rows are aligned to unit levels. |
| nfMin | `hmsc_nfMin` | Covered | Tested with sample and xData random-level priors; slider and numeric field are synced. |
| nfMax | `hmsc_nfMax` | Covered | Tested with sample and xData random-level priors; slider and numeric field are synced. |

## Constructor Advanced Options

| Control | inputId | Coverage | Notes |
|---|---|---|---|
| XScale | `hmsc_XScale` | Covered | Tested by `constructor_offsets`. |
| TrScale | `hmsc_TrScale` | Covered | Tested by trait cases. |
| YScale | `hmsc_YScale` | Covered | Tested by `constructor_offsets`. |
| truncateNumberOfFactors | `hmsc_truncateNumberOfFactors` | Covered | Tested by `constructor_offsets`. |
| Loff offset file | `hmsc_Loff_file` | Covered | Tested by `constructor_offsets`; dimensions are checked. |
| ranLevelsUsed | `hmsc_ranLevelsUsed` | Covered | Tested by `sample_priors_partition`. |
| C correlation matrix file | `hmsc_C_file` | Covered | Tested by `phylogeny_C_traits`; matrix is repaired to positive definite when needed. |
| Use XRRR reduced-rank regression | `hmsc_use_XRRR` | Covered | Tested by `xrrr`. |
| ncRRR | `hmsc_ncRRR` | Covered | Tested by `xrrr`; slider and numeric field are synced. |
| XRRRFormula | `hmsc_XRRRFormula` | Covered | Tested by `xrrr`. |
| XRRRScale | `hmsc_XRRRScale` | Covered | Tested by `xrrr`. |
| XRRRData file name | `hmsc_XRRR_file` | Covered | Tested by `xrrr`; row count is checked against Y. |

## MCMC And Priors

| Control | inputId | Coverage | Notes |
|---|---|---|---|
| samples | `hmsc_samples` | Covered | Used by `sampleMcmc()`; slider and numeric field are synced. |
| transient | `hmsc_transient` | Covered | Used by `sampleMcmc()`; max is 1,000,000. |
| thin | `hmsc_thin` | Covered | Used by `sampleMcmc()`; max is 100,000. |
| nChains | `hmsc_nChains` | Covered | Used by `sampleMcmc()`; tested with two chains. |
| nParallel | `hmsc_nParallel` | Covered | Used by `sampleMcmc()`; Windows default remains 1. |
| verbose | `hmsc_verbose` | Covered | Used by `sampleMcmc()` progress reporting. |
| seed | `hmsc_seed` | Covered | Sets reproducible seed before model fitting. |
| initPar | `hmsc_initPar` | Covered | Normalized for Hmsc; unsupported non-default strategies fall back safely. |
| alignPost | `hmsc_alignPost` | Covered | Tested by `sample_prior_pool`. |
| sample_prior | `hmsc_sample_prior` | Covered | Tested by `sample_prior_pool`; missing posterior metrics are handled safely. |
| pool_chains | `hmsc_pool_chains` | Covered | Tested by `sample_prior_pool`; uses `poolMcmcChains()` when enabled. |
| Use HMSC default priors | `hmsc_use_default_priors` | Covered | Tested by `sample_priors_partition`. |
| updater$GammaEta | `hmsc_updater_GammaEta` | Covered | Passed through the supported `updater` object; Hmsc may force it off for no-random-effect or NNGP/GPP models. |
| updater$Beta | `hmsc_updater_Beta` | Package-limited | Recorded in config/audit. Hmsc 3.3.7 does not expose separate formal `sampleMcmc(..., Beta=)` arguments. |
| updater$Gamma | `hmsc_updater_Gamma` | Package-limited | Recorded in config/audit. Not blindly passed as unsupported fake formal. |
| updater$Omega | `hmsc_updater_Omega` | Package-limited | Recorded in config/audit. Not blindly passed as unsupported fake formal. |
| a1 | `hmsc_a1` | Covered | Tested by `sample_priors_partition`. |
| b1 | `hmsc_b1` | Covered | Tested by `sample_priors_partition`. |
| a2 | `hmsc_a2` | Covered | Tested by `sample_priors_partition`. |
| b2 | `hmsc_b2` | Covered | Tested by `sample_priors_partition`. |
| alphapw | `hmsc_alphapw` | Covered | Tested by `spatial_full_alphapw`; parsed as numeric distance/probability pairs. |

## Results, Diagnostics And Plots

| Control | inputId | Coverage | Notes |
|---|---|---|---|
| computePredictedValues | `hmsc_out_predicted` | Covered | `FALSE` tested by `output_and_plot_switches`; TRUE tested by real examples. |
| evaluateModelFit | `hmsc_out_fit` | Covered | `FALSE` tested by `output_and_plot_switches`; TRUE tested by real examples. |
| Cross-validation | `hmsc_out_cv` | Covered | TRUE tested in most examples; FALSE tested by output switch and sample-prior cases. |
| computeWAIC | `hmsc_out_waic` | Covered | TRUE tested by `phylogeny_newick_outputs` and `output_and_plot_switches`. |
| MCMC diagnostics | `hmsc_out_diag` | Covered | FALSE tested by `output_and_plot_switches`; TRUE tested elsewhere. |
| Parameter estimates | `hmsc_out_params` | Covered | FALSE tested by `output_and_plot_switches`; TRUE tested elsewhere. |
| Variance partitioning | `hmsc_out_vp` | Covered | FALSE tested by `output_and_plot_switches`; TRUE tested elsewhere. |
| Associations / Omega | `hmsc_out_omega` | Covered | TRUE tested by `omega_plot_controls`; FALSE tested by `output_and_plot_switches`. |
| Environmental gradients | `hmsc_out_gradients` | Covered | TRUE path tested through S7 metadata/prediction settings; FALSE tested by output switch. |
| Beta support | `hmsc_beta_support` | Covered | TRUE and FALSE paths tested. |
| Gamma support | `hmsc_gamma_support` | Covered | TRUE and FALSE paths tested. |
| Omega support | `hmsc_omega_support` | Covered | TRUE and FALSE paths tested. |
| Convergence Beta | `hmsc_conv_beta` | Covered | Controls S3 Beta diagnostics. |
| Convergence Gamma | `hmsc_conv_gamma` | Covered | Controls S3 Gamma diagnostics. |
| Convergence Omega | `hmsc_conv_omega` | Covered | Controls S3 Omega diagnostics. |
| maxOmega | `hmsc_maxOmega` | Covered | Tested by output switch and Omega plotting case; slider and numeric field are synced. |
| Convergence rho | `hmsc_conv_rho` | Covered | Controls S3 rho diagnostics when phylogeny is present. |
| Convergence alpha | `hmsc_conv_alpha` | Covered | Controls S3 alpha diagnostics when spatial random levels exist. |
| effectiveSize | `hmsc_effective_size` | Covered | TRUE/FALSE paths tested; unnamed vectors get safe fallback names. |
| Gelman PSRF | `hmsc_gelman_psrf` | Covered | TRUE/FALSE paths tested; unnamed vectors get safe fallback names. |
| plotBeta | `hmsc_plot_beta` | Covered | TRUE and FALSE paths tested. |
| plotGamma | `hmsc_plot_gamma` | Covered | TRUE and FALSE paths tested. |
| VP order explained | `hmsc_varpart_order_explained` | Covered | TRUE default and FALSE path tested. |
| VP order raw | `hmsc_varpart_order_raw` | Covered | TRUE path tested by `omega_plot_controls`. |
| Show species names in Beta plots | `hmsc_show_sp_names_beta` | Covered | TRUE path tested by `omega_plot_controls`. |
| plotTree | `hmsc_plotTree` | Covered | Tested by phylogeny/trait case; passed to Beta plotting. |
| Omega order | `hmsc_omega_order` | Covered | UI values now map to valid corrplot order values; `alphabetical` normalization tested. |
| Show species names in Omega | `hmsc_show_sp_names_omega` | Covered | FALSE path tested by `omega_plot_controls`. |
| species.list | `hmsc_species_list` | Covered | Tested by `phylogeny_newick_outputs`; saved in S7 prediction settings. |
| trait.list | `hmsc_trait_list` | Covered | Tested by `phylogeny_newick_outputs`; saved in S7 prediction settings. |
| env.list | `hmsc_env_list` | Covered | Tested by constructor and phylogeny output cases; filters S7 gradients/prediction settings. |
| nfolds | `hmsc_nfolds` | Covered | Tested with 2 and 3 folds; slider and numeric field are synced. |
| partition column | `hmsc_partition_column` | Covered | Tested by `sample_priors_partition`; used in `createPartition()`. |
| computeSAIR | `hmsc_compute_sair` | Metadata-limited | Recorded in S7 settings. A dedicated SAIR numerical adapter is still future work. |

## Output Workflow And Download

| Control | inputId | Coverage | Notes |
|---|---|---|---|
| Copy input files | `hmsc_out_copy_inputs` | Covered | Default TRUE path copies uploads to `inputs/`; real pipeline also writes clean data to `data/`. |
| Export executable R scripts | `hmsc_export_scripts` | Covered | Real pipeline always exports executable S1-S7/S8 scripts for review reproducibility. |
| Create standard comparison tables | `hmsc_create_standard` | Covered | Real pipeline writes standard tables; scaffold path also honors this switch. |
| Run real HMSC S1-S7 results workflow | `hmsc_real_fit` | Covered | TRUE path tested extensively; FALSE remains explicit scaffold mode. |
| Run HMSC workflow | `hmsc_run` | Covered | Runs validation, fitting/scaffold, status files and ZIP creation. |
| Download HMSC ZIP | `hmsc_download` | Covered | Uses the real generated ZIP path; check-failed runs also create a real diagnostics ZIP. |

## Required Diagnostics If A Run Fails

- `diagnostics/engine_status.json`
- `diagnostics/HMSC_S1S7_pipeline_status.json`
- `diagnostics/HMSC_S1S7_error.txt`
- `diagnostics/data_check_messages.csv`
- `diagnostics/session_info.txt`
- `tables/HMSC_parameter_audit.csv`
