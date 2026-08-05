# JSDM Studio Universal Benchmark Integration Audit

Date: 2026-06-03

## Scope

This audit covers the new top-level **Universal Benchmark** workflow in JSDM Studio. The workflow derives engine-specific inputs from one latent ecological truth, runs engines sequentially on Windows, writes per-engine output contracts and ZIP files, then builds one all-model benchmark folder and ZIP.

## Modified files

- `app.R`
  - Added the `Universal Benchmark` tab after all single-engine workflows and before `Compare Models`.
  - Added UI controls for benchmark design, data generation, dependency preflight, sequential all-engine run, unified result tables, comparison summaries and all-model ZIP download.
  - Added server handlers for `univ_generate`, `univ_preflight`, `univ_run_all`, `univ_compare`, result tables and `univ_download`.
  - Added Hmsc-HPC Python package-level preflight checks.
- `workflow_scripts/universal_benchmark_runner.R`
  - Added command-line actions: `generate`, `preflight`, `all`, `compare`, `zip`.
  - Added deterministic synthetic benchmark data generation under `examples/universal_benchmark/`.
  - Added sequential runners for Hmsc, Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM and boral using existing adapters/output contracts.
  - Added standard synonym tables per engine:
    - `standard/effects_species_environment.csv`
    - `standard/predictions_site_species.csv`
    - `standard/associations_species_species.csv`
  - Added master result aggregation, truth comparison, comparability rules, reports and all-model ZIP generation.
  - Hardened aggregation for empty/missing result classes and different engine-specific table columns.

## Generated benchmark data

Generated under:

`examples/universal_benchmark/`

Required files are present, including occurrence/count/normal Y matrices, XData, traits, study design, coordinates, phylogeny, phylogenetic covariance, prediction data, folds, offset/trial/row/random-effect files, spOccupancy detection data and truth tables.

The data are small and stable for smoke testing while still covering repeated detection:

- `n_sites = 36`
- `n_species = 6`
- `n_visits = 3`
- seed `20260601`

## Final benchmark run

Master folder:

`output/universal_benchmark_20260603_155335/master`

Master ZIP:

`output/universal_benchmark_20260603_155335/JSDMStudio_universal_benchmark_ALL_MODELS.zip`

Master ZIP size: `16.36 MB`

SHA256:

`8d56f7a613ab44ec14b9114f5fb6abc14f2bb667c1fd3b7e2af9cd82f66fb024`

## Engine status

| Engine | Status | Runtime seconds | Files | ZIP size KB |
|---|---:|---:|---:|---:|
| Hmsc | fitted | 3.00 | 98 | 5630.3 |
| Hmsc-HPC | fitted | 17.49 | 124 | 110.9 |
| jSDM | fitted | 0.94 | 96 | 320.9 |
| GJAM | fitted | 1.40 | 128 | 815.6 |
| spOccupancy | fitted | 1.87 | 88 | 602.7 |
| sjSDM | fitted | 24.67 | 93 | 769.2 |
| boral | fitted | 1.75 | 85 | 136.7 |

All seven engines produced real per-engine ZIP files and the required output contract folders:

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

All seven engines include the required standard tables:

- `standard/run_summary.csv`
- `standard/effects_long.csv`
- `standard/predictions_long.csv`
- `standard/fit_metrics.csv`
- `standard/diagnostics_long.csv`
- `standard/output_manifest.csv`
- `standard/effects_species_environment.csv`
- `standard/predictions_site_species.csv`
- `standard/associations_species_species.csv`

## Master output contract

The master output contains:

- `01_benchmark_inputs/`
- `02_truth/`
- `03_engine_outputs/`
- `04_unified_standard_results/`
- `05_truth_comparison/`
- `06_compare_models/`
- `07_reports/`
- `08_reproducible_scripts/`
- `09_configs/`

Important master outputs:

- `04_unified_standard_results/model_status_matrix.csv`
- `04_unified_standard_results/effects_all_models_long.csv`
- `04_unified_standard_results/effects_species_environment_wide.csv`
- `04_unified_standard_results/predictions_all_models_long.csv`
- `04_unified_standard_results/predictions_site_species_wide.csv`
- `04_unified_standard_results/fit_metrics_all_models.csv`
- `04_unified_standard_results/associations_all_models_long.csv`
- `04_unified_standard_results/associations_species_species_wide.csv`
- `05_truth_comparison/truth_vs_estimated_effects.csv`
- `05_truth_comparison/truth_vs_predictions.csv`
- `05_truth_comparison/effect_direction_agreement.csv`
- `05_truth_comparison/prediction_metric_by_model.csv`
- `05_truth_comparison/association_recovery_summary.csv`
- `06_compare_models/comparable_results_matrix.csv`
- `06_compare_models/non_comparable_results_notes.csv`
- `06_compare_models/compare_models_summary.csv`
- `06_compare_models/compare_models_report.html`
- `07_reports/master_report.html`
- `07_reports/master_report.md`
- `07_reports/benchmark_methods_summary.md`
- `07_reports/dependency_preflight_report.csv`
- `07_reports/run_log.txt`
- `07_reports/session_info.txt`
- `07_reports/master_output_manifest.csv`

## Dependency preflight

The preflight report now checks both R engines and Hmsc-HPC Python packages. In the final run:

- Hmsc, jSDM, GJAM, spOccupancy, sjSDM, boral, rjags and R2jags were available.
- JAGS was loadable through `rjags`.
- Hmsc-HPC Python packages were available:
  - `numpy 1.26.4`
  - `pandas 2.3.3`
  - `patsy 1.0.2`
  - `PyYAML 6.0.2`
  - `h5py 3.16.0`
  - `scipy 1.15.3`
  - `tensorflow 2.16.1`
  - `tensorflow_probability 0.24.0`
  - `ujson 5.12.1`
  - `pyhmsc` from `external_packages/hmsc-hpc-main`

Note: standalone R package `torch` is still reported as not installed by `requireNamespace("torch")`, but the sjSDM package and the benchmark run completed successfully.

## Validation commands

Passed:

- `app.R` parse
- `app.R` source
- `workflow_scripts/universal_benchmark_runner.R` parse
- strict Shiny input/output/download audit
- dependency preflight
- Universal Benchmark data generation
- Universal Benchmark all-engine run
- Universal Benchmark compare rebuild
- Universal Benchmark all-model ZIP rebuild
- HTTP smoke test: `GET http://127.0.0.1:6939` returned status 200 and contained `JSDM Studio` and `Universal Benchmark`
- per-engine output contract and ZIP inspection
- master output contract and SHA256 inspection

Strict Shiny audit summary:

- UI input IDs: 528
- server input references: 528
- missing UI IDs for server refs: 0
- UI inputs not read: 0
- output UI IDs: 83
- output handlers: 83
- outputs without handler: 0
- handlers without UI: 0
- download buttons: 9
- download handlers: 9
- broken static download handlers: 0

## Comparability rules

Directly comparable:

- workflow status
- output-contract completeness
- runtime
- real ZIP existence and size
- prediction metrics in this controlled benchmark
- truth-vs-prediction residuals
- effect direction agreement

Not directly interchangeable as one statistical parameter:

- Hmsc Omega
- Hmsc-HPC Eta/Lambda-derived summaries
- jSDM residual correlations
- GJAM corMu/sigMu
- spOccupancy latent-factor/random-effect associations
- sjSDM covariance/correlation summaries
- boral residual correlations

The master reports and association tables include `association_type`, `scale`, `comparable` and `note` fields.

## Remaining risks

- The benchmark uses quick/smoke-test sampler settings. It verifies software execution and output contracts, not publication-grade convergence.
- Hmsc-HPC was validated on the CPU pyhmsc path. GPU/cluster/HPC scheduling was not added or claimed.
- Hmsc-HPC Python dependencies were installed in the currently selected Python environment during validation; a fresh machine must run `install_hmschpc_python_packages.bat` or provide an equivalent Python environment.
- `torch` is not installed as a standalone R package, although sjSDM itself ran successfully in the benchmark.
- The benchmark compares outputs derived from one synthetic latent truth. Real ecological inference still requires model-specific data checks, sufficient sample sizes and convergence diagnostics.
