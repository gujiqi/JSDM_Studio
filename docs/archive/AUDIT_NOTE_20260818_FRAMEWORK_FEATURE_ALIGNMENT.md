# Framework Figure Feature Alignment Note

Date: 2026-08-18

This note records a targeted follow-up audit against the manuscript conceptual framework figure.

## Checked

The following figure concepts were mapped to implementation:

- Ecological input data: implemented through workflow upload panels and `examples/universal_benchmark/`.
- JSDM Studio reproducibility layer: implemented through data checks, parameter capture in config files, engine adapters, executable scripts, diagnostics, standard tables and ZIP outputs.
- Model engines: implemented for Hmsc, Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM and boral.
- Auditable outputs: implemented through `models/`, `tables/`, `plots/`, `predictions/`, `diagnostics/`, `report/`, `standard/` and engine ZIPs.
- Scientific use: supported through reproducible scripts, reports, model-comparison outputs and shareable ZIP packages. Peer review and collaboration are use cases enabled by these outputs, not separate statistical engines.
- Universal Benchmark: implemented through shared latent ecological truth, engine-specific inputs, unified effects/predictions/associations and truth-comparison outputs.

## Follow-up fix

The figure used the phrase "accuracy, calibration and rank" for truth comparison. The existing Universal Benchmark already wrote prediction RMSE to truth probability, but the calibration and rank summaries were not explicit columns in `prediction_metric_by_model.csv`.

`workflow_scripts/universal_benchmark_runner.R` was updated so `05_truth_comparison/prediction_metric_by_model.csv` now includes:

- `rmse_to_truth_probability`
- `calibration_intercept`
- `calibration_slope`
- `pearson_correlation`
- `spearman_rank_correlation`
- `n_predictions`

The existing Universal Benchmark master output was regenerated with `--action=compare`.

## Verification

Commands run:

```text
Rscript -e "parse('workflow_scripts/universal_benchmark_runner.R')"
Rscript workflow_scripts/universal_benchmark_runner.R --action=compare --master_dir='output/universal_benchmark_20260812_114724/master'
```

Results:

```text
UNIVERSAL_RUNNER_PARSE_OK
COMPARE_OK
```

The updated `prediction_metric_by_model.csv` contains explicit calibration and rank-correlation fields.
