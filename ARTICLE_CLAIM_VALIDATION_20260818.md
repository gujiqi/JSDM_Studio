# Article Claim Validation for JSDM Studio

Date: 2026-08-18

This report checks whether the main claims proposed for the MEE Application manuscript draft are implemented in the current JSDM Studio code and outputs. The validation uses fresh runtime evidence from the current machine.

## Runtime Commands

```text
Rscript workflow_scripts/workflow_synthetic_case_suite.R
Rscript workflow_scripts/universal_benchmark_runner.R --action=all
```

## Synthetic Workflow Case Suite

Output root:

```text
C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio/output/workflow_synthetic_case_audit_20260818_144112
```

Summary CSV:

```text
C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio/output/workflow_synthetic_case_audit_20260818_144112/workflow_synthetic_case_summary.csv
```

Result:

```text
WORKFLOW_SUITE_OK cases=28 missing_sum=0
```

| Engine | Cases | Status | Missing outputs |
| --- | ---: | --- | ---: |
| Hmsc | 5 | fitted | 0 |
| Hmsc-HPC | 3 | fitted | 0 |
| jSDM | 5 | fitted | 0 |
| GJAM | 5 | fitted | 0 |
| spOccupancy | 4 | fitted | 0 |
| sjSDM | 3 | fitted | 0 |
| boral | 3 | fitted | 0 |

All 28 cases produced non-empty ZIP files.

## Universal Benchmark

Output root:

```text
C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio/output/universal_benchmark_20260818_144351
```

Master directory:

```text
C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio/output/universal_benchmark_20260818_144351/master
```

Master ZIP:

```text
C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio/output/universal_benchmark_20260818_144351/JSDMStudio_universal_benchmark_ALL_MODELS.zip
```

Master ZIP SHA256:

```text
ebcbeecc1aecd91b5015f574b4ab2933d5e51fb95f1cbb453dae824b4eaffafc
```

| Engine | Status | Runtime seconds | File count | ZIP size KB |
| --- | --- | ---: | ---: | ---: |
| Hmsc | fitted | 2.43 | 98 | 5985.4 |
| Hmsc-HPC | fitted | 14.41 | 124 | 111.2 |
| jSDM | fitted | 0.79 | 96 | 320.9 |
| GJAM | fitted | 1.17 | 128 | 815.6 |
| spOccupancy | fitted | 1.50 | 88 | 602.7 |
| sjSDM | fitted | 20.35 | 93 | 769.5 |
| boral | fitted | 1.36 | 85 | 137.0 |

Each engine output contains:

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
- a non-empty engine ZIP

## Universal Benchmark Master Sections

| Section | Files |
| --- | ---: |
| `01_benchmark_inputs/` | 28 |
| `02_truth/` | 7 |
| `03_engine_outputs/` | 719 |
| `04_unified_standard_results/` | 10 |
| `05_truth_comparison/` | 6 |
| `06_compare_models/` | 4 |
| `07_reports/` | 7 |
| `08_reproducible_scripts/` | 4 |
| `09_configs/` | 8 |

## Manuscript Claim Matrix

| Manuscript claim | Implementation status | Evidence |
| --- | --- | --- |
| JSDM Studio is an application layer, not a new statistical model. | Implemented | Separate engine adapters and engine-specific workflows are present. |
| It supports Hmsc, Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM and boral. | Implemented and tested | Synthetic suite and Universal Benchmark fitted all seven engines. |
| It performs data checking and validation. | Implemented | Each workflow has check logic; outputs include `diagnostics/data_check_messages.csv`. |
| It captures parameters and settings. | Implemented | Each engine output includes `used_config.yml`; benchmark master includes `09_configs/`. |
| It exports executable scripts. | Implemented | Each engine output includes `reproducible_script/`; master includes `08_reproducible_scripts/`. |
| It writes diagnostics and honest status. | Implemented | Each engine output includes `diagnostics/engine_status.json`; all tested cases were `fitted`. |
| It writes standard tables. | Implemented | `standard/` exists for all engines; master `04_unified_standard_results/` has 10 tables. |
| It packages real ZIP archives. | Implemented | All 28 synthetic cases and all 7 Universal Benchmark engines wrote non-empty ZIP files. |
| It derives benchmark inputs from a shared latent ecological truth. | Implemented | `examples/universal_benchmark/truth/` and master `02_truth/` are populated. |
| It writes unified effects, predictions and associations. | Implemented | Master includes long and wide effects, predictions and associations tables. |
| It compares predictions to truth. | Implemented | `truth_vs_predictions.csv` and `prediction_metric_by_model.csv` are present. |
| It includes accuracy, calibration and rank metrics. | Implemented | `prediction_metric_by_model.csv` includes RMSE, calibration intercept/slope, Pearson correlation and Spearman rank correlation. |
| It distinguishes non-comparable association parameters. | Implemented | Compare outputs include `non_comparable_results_notes.csv`; associations include type/scale/comparability notes. |
| It supports teaching, peer review and collaboration. | Supported as use cases | Reports, reproducible scripts, diagnostics and shareable ZIPs support these use cases; they are not separate statistical engines. |

## Caution for Manuscript Wording

The validation supports software-execution claims, output-contract claims and reproducibility-layer claims. It does not prove publication-quality ecological inference, MCMC convergence for real datasets, or general superiority of one JSDM engine over another. The manuscript should state that small synthetic cases validate execution and output integrity, while real ecological interpretation still requires expert model checking.

## Verdict

The main Application-paper claims currently proposed for JSDM Studio are implemented and supported by fresh runtime evidence on this machine.
