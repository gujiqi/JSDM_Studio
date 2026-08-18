# Workflow Synthetic Audit 2026-08-11

This audit ran small synthetic cases for all JSDM Studio engines:
Hmsc, Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM and boral.

Audit command:

```r
Rscript workflow_scripts/workflow_synthetic_case_suite.R
```

Latest audit root:

```text
C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio/output/workflow_synthetic_case_audit_20260811_225215
```

## Fixes Applied

- Fixed sjSDM DNN coefficient handling in `app.R`. DNN weights are now saved as reproducible weight objects and manifests instead of being forced into a species-by-predictor coefficient matrix.
- Fixed formula validation in `app.R` so one-sided formulas using `.` are not falsely flagged as missing variable `.`.
- Fixed jSDM long-format adapter logic in `R/jsdm_adapter.R`. If a long-format formula lacks required `species:` interactions, the adapter rewrites it to a jSDM-compatible formula and records a warning.
- Added `workflow_scripts/workflow_synthetic_case_suite.R` to reproduce the multi-engine workflow audit.
- Added stable benchmark coverage cases for Hmsc Spatial Full and GJAM PA/CON/DA mixed-scale data.

## Output Contract Check

Every case was checked for:

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
- a real non-empty ZIP archive

All 27 cases had `missing_count = 0`.

## Case Results

| Engine | Case | Status | Files | ZIP KB | Note |
|---|---|---:|---:|---:|---|
| Hmsc | hmsc_probit_sample_traits_phylo | fitted | 93 | 5858.6 | Probit, sample random effect, traits and phylogeny. |
| Hmsc | hmsc_poisson_spatial_nngp | fitted | 92 | 6236.5 | Poisson counts, spatial NNGP. |
| Hmsc | hmsc_probit_spatial_full | fitted | 92 | 6226.0 | Probit, spatial Full Gaussian process. |
| Hmsc | hmsc_probit_spatial_gpp_traits | fit_failed | 68 | 2040.7 | Hmsc GPP branch failed with `missing value where TRUE/FALSE needed`; diagnostics and ZIP complete. |
| Hmsc | hmsc_normal_no_random_guardrail | fitted | 88 | 6259.2 | Normal family guardrail without random effect. |
| Hmsc-HPC | hmschpc_probit_iid_traits_phylo | fitted | 127 | 190.3 | CPU Hmsc-HPC probit, iid random level, traits and phylogeny. |
| Hmsc-HPC | hmschpc_poisson_no_random | fitted | 118 | 193.1 | CPU Hmsc-HPC poisson, no random level. |
| Hmsc-HPC | hmschpc_normal_spatial_full | fitted | 127 | 201.2 | CPU Hmsc-HPC normal, spatial_full. |
| jSDM | jsdm_binomial_probit_latent_random_traits | fitted | 98 | 290.9 | Probit matrix, latent variables, random site effect, traits. |
| jSDM | jsdm_binomial_logit_trials | fitted | 89 | 111.8 | Binomial logit with trials. |
| jSDM | jsdm_poisson_log_counts | fitted | 89 | 111.0 | Poisson log counts. |
| jSDM | jsdm_gaussian_continuous | fitted | 89 | 143.2 | Gaussian continuous response. |
| jSDM | jsdm_binomial_probit_long_format | fitted | 90 | 140.0 | Long-format probit branch after formula rewrite fix. |
| GJAM | gjam_pa_traits_prediction | fitted | 128 | 749.2 | PA data, traits and prediction. |
| GJAM | gjam_mixed_pa_con_da_stable | fitted | 120 | 290.2 | Stable mixed PA, CON and DA case. |
| GJAM | gjam_mixed_pa_con_ca_fc_oc | fit_failed | 69 | 38.4 | Tiny mixed FC/OC case remains numerically singular; diagnostics and ZIP complete. |
| GJAM | gjam_discrete_abundance_counts | fitted | 124 | 709.7 | DA count-style response branch. |
| spOccupancy | spocc_msPGOcc_replicated_multi | fitted | 88 | 517.5 | Multi-species replicated occupancy. |
| spOccupancy | spocc_PGOcc_single_species | fitted | 83 | 163.5 | Single-species non-spatial PGOcc. |
| spOccupancy | spocc_spPGOcc_single_spatial | fitted | 86 | 218.6 | Single-species spatial PGOcc. |
| spOccupancy | spocc_lfMsPGOcc_latent_factors | fitted | 91 | 554.9 | Multi-species latent-factor occupancy. |
| sjSDM | sjsdm_binomial_linear_spatial | fitted | 93 | 687.3 | Binomial, linear environment and spatial predictors. |
| sjSDM | sjsdm_gaussian_no_spatial | fitted | 87 | 257.3 | Gaussian response without spatial model. |
| sjSDM | sjsdm_binomial_dnn_environment | fitted | 90 | 681.8 | CPU DNN environment branch; weights saved as DNN artifacts. |
| boral | boral_binomial_latent | fitted | 85 | 124.2 | Binomial with latent variable. |
| boral | boral_poisson_row_effect | fitted | 85 | 188.4 | Poisson with row effect. |
| boral | boral_normal_offset | fitted | 85 | 130.5 | Normal response with offset. |

## Remaining Risks

- Hmsc GPP is still a risk in very small quick-test runs. Spatial Full and NNGP fitted; GPP produced complete `fit_failed` diagnostics instead of a false success.
- GJAM FC/OC mixed-scale examples are sensitive to tiny synthetic data and can become singular. A stable PA/CON/DA mixed case now fits; FC/OC should be tested with larger and carefully balanced data before claiming fitted status.
- All quick-test MCMC settings are software smoke tests only. They are not convergence evidence and should not be used for publication inference.

## Diagnostics

For failed cases, inspect:

- `Hmsc/hmsc_probit_spatial_gpp_traits/diagnostics/engine_status.json`
- `Hmsc/hmsc_probit_spatial_gpp_traits/diagnostics/HMSC_S1S7_error.txt`
- `GJAM/gjam_mixed_pa_con_ca_fc_oc/diagnostics/engine_status.json`
- `GJAM/gjam_mixed_pa_con_ca_fc_oc/diagnostics/GJAM_reproducible_error.txt`

