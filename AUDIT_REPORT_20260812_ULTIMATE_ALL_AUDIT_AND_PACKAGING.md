# JSDM Studio Ultimate Final Audit and Packaging Report

Date: 2026-08-12

This report records the final source-to-installer audit for the current JSDM Studio build. The audit was run as an evidence-first release check: parsing, source loading, Shiny binding checks, HTTP smoke testing, Compare Models testing, workflow synthetic suites, Universal Benchmark execution, launcher build, installer build, and ZIP contract checks.

## Project Audited

Project root:

`C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio`

R executable:

`C:/Program Files/R/R-4.5.3/bin/Rscript.exe`

## Scope

Audited source and release areas:

- `app.R`
- `R/`
- `workflow_scripts/`
- `examples/`
- `docs/`
- `external_packages/`
- `installer/`
- `output/`
- WebView2 launcher files
- Inno Setup installer files
- `.bat` startup/build scripts
- workflow ZIP outputs
- `diagnostics/`
- `standard/`
- `report/`
- `reproducible_script/`

## Static Validation

### R Parse

Command class:

`Rscript -e "parse app.R, R/*.R, workflow_scripts/*.R, tests/*.R"`

Result:

`ALL_R_PARSE_OK files=16`

### app.R Parse and Source

Result:

`APP_PARSE_SOURCE_OK`

### Shiny UI/Server Binding Audit

Script:

`workflow_scripts/static_shiny_binding_audit.R`

Result:

```text
STRICT_SHINY_BINDING_AUDIT
ui_ids=533 input_refs=533 refs_not_in_ui=0 ui_not_read=0
output_ui_ids=87 output_handlers=87 outputs_without_handler=0 handlers_without_ui=0 non_download_without_render=0
download_buttons=10 download_handlers=10 downloads_without_handler=0 handlers_without_button=0
download_content_static_bad=0
STRICT_AUDIT_OK
```

Interpretation:

- Every referenced `input$...` is represented in the UI.
- Every UI input is read by server logic or documented binding logic.
- Every output placeholder has a render function or download handler.
- Every download button has a corresponding `downloadHandler`.
- No statically detectable fake or empty download handler remained.

### Compare Models Test

Script:

`workflow_scripts/test_compare_models.R`

Result:

```text
COMPARE_TEST_OK
COMPARE_TEST_FIXTURE=C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio/output/CompareModelsTest_CompareModels_test_Hmsc_20260812_114316
COMPARE_TEST_ZIP=C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio/output/CompareModelsTest_CompareModels_test_Hmsc_20260812_114316.zip
```

### HTTP Smoke Test

Script:

`tests/http_smoke.R`

Result:

`HTTP_SMOKE_OK port=7856 status=200 bytes=759844`

## Workflow Synthetic Suite

Script:

`workflow_scripts/workflow_synthetic_case_suite.R`

Output root:

`C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio/output/workflow_synthetic_case_audit_20260812_114447`

Summary table:

`output/workflow_synthetic_case_audit_20260812_114447/workflow_synthetic_case_summary.csv`

Result:

```text
WORKFLOW_SUITE_OK all_cases_fitted=28 missing_count_sum=0 zips_nonempty=28
```

Engine case coverage:

| Engine | Cases | Status | Missing outputs |
| --- | ---: | --- | ---: |
| Hmsc | 5 | fitted | 0 |
| Hmsc-HPC | 3 | fitted | 0 |
| jSDM | 5 | fitted | 0 |
| GJAM | 5 | fitted | 0 |
| spOccupancy | 4 | fitted | 0 |
| sjSDM | 3 | fitted | 0 |
| boral | 3 | fitted | 0 |

Coverage includes small synthetic cases for occurrence/binomial/probit-compatible workflows, count/poisson workflows, normal/gaussian workflows, traits, phylogeny, random or spatial effects where supported, detection data for spOccupancy, latent factors, MCMC controls, prediction outputs, diagnostics, plots, `standard/` tables, and real ZIP outputs.

## Universal Benchmark

Script:

`workflow_scripts/universal_benchmark_runner.R --action=all`

Benchmark input directory:

`examples/universal_benchmark`

Master output directory:

`C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio/output/universal_benchmark_20260812_114724/master`

Master benchmark ZIP:

`C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio/output/universal_benchmark_20260812_114724/JSDMStudio_universal_benchmark_ALL_MODELS.zip`

Master ZIP SHA256:

`b61d43fba42b78153fbbda2c41d076e0a597d3f025360e743bfffad83cc037f9`

Universal Benchmark statuses:

| Engine | Status | Runtime seconds | File count | Engine ZIP size KB |
| --- | --- | ---: | ---: | ---: |
| Hmsc | fitted | 2.54 | 98 | 6068.1 |
| Hmsc-HPC | fitted | 14.97 | 124 | 111.2 |
| jSDM | fitted | 0.82 | 96 | 320.9 |
| GJAM | fitted | 1.14 | 128 | 815.6 |
| spOccupancy | fitted | 1.48 | 88 | 602.7 |
| sjSDM | fitted | 20.26 | 93 | 770.0 |
| boral | fitted | 1.48 | 85 | 136.8 |

Contract check result:

`UNIVERSAL_MASTER_CONTRACT_OK`

Compare action:

`COMPARE_OK`

Universal Benchmark standard output includes engine-specific and unified tables for:

- `effects_species_environment.csv`
- `predictions_site_species.csv`
- `associations_species_species.csv`
- `run_summary.csv`
- `effects_long.csv`
- `predictions_long.csv`
- `fit_metrics.csv`
- `diagnostics_long.csv`
- `output_manifest.csv`

The master output also includes truth-comparison tables, model-comparison summaries, reproducible scripts, configs, reports, and independent engine ZIP files.

## Scientific Text and Status Semantics Audit

The main UI and scripts were scanned for Chinese characters and mojibake. Core UI and executable scripts returned no matches for Chinese or replacement characters.

The following scientific/status statements were confirmed in the app:

- `fitted` means the workflow completed and wrote outputs; it does not prove convergence or publication readiness.
- `model_defined` means inputs and scripts were generated, but fitting was intentionally skipped or unavailable.
- `check_failed` and `fit_failed` are explicit failure states with diagnostics.
- HMSC/Hmsc-HPC association outputs are not equated blindly with jSDM, GJAM, sjSDM, spOccupancy, or boral association parameters.
- spOccupancy detection effects are separated from occurrence/environment effects.
- GJAM observation-scale and typeNames semantics are documented.
- Compare Models explains which outputs are comparable and which are not directly comparable.

Legacy Chinese documentation files remain in `docs/` and root markdown files because this is a maximum archival bundle. They do not affect the English Shiny UI or executable workflow scripts.

## Installer and Launcher

WebView2 launcher build:

`Build complete`

Launcher:

`JSDMStudioLauncher.exe`

Installer:

`installer/output/JSDMStudio_Setup.exe`

Installer size:

`143034433 bytes`

Installer SHA256:

`0484823162080F554C38B10068038A9BF25DBE665950F4AAFDE9957A4BBBE16A`

## Output Contract

Validated workflow outputs include:

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
- real non-empty engine ZIP files

Failure-state logic was checked by static audit and workflow suite behavior. The final suite contained no remaining workflow failures.

## Remaining Risks

1. External dependency sensitivity remains unavoidable for real user machines. Hmsc, jSDM, GJAM, spOccupancy, sjSDM, boral, JAGS, PyTorch, reticulate, and Hmsc-HPC Python can still fail on systems with broken package installations. The app records these as `check_failed` or `fit_failed` with diagnostics.
2. Synthetic cases are intentionally small for release validation. Publication analyses still require larger MCMC settings, convergence review, model-specific diagnostics, and domain-specific study design.
3. Cross-engine association parameters are reported in a unified table for navigation, but the report explicitly marks their scale and comparability. They should not be interpreted as identical statistical estimands.
4. Legacy archival documentation includes Chinese-language files. This is retained for the maximum bundle, but the executable UI and core scripts were checked for English/mojibake cleanliness.

## Final Status

The audited build passed:

- R parse
- app source load
- strict Shiny binding audit
- Compare Models test
- HTTP smoke test
- complete workflow synthetic suite
- Universal Benchmark all-engine run
- Universal Benchmark compare action
- Universal Benchmark ZIP contract
- WebView2 launcher build
- Inno Setup installer build

Final release packaging is performed after this report is written so the final ZIP hash can be calculated on the immutable archive.
