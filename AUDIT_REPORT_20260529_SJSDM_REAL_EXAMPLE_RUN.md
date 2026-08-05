# sjSDM Real Example Run Audit - 2026-05-29

## Scope

This audit specifically checked whether the bundled sjSDM example can run a real package fit and produce real, non-empty results rather than a scaffold.

## Fixes Applied

- Fixed `examples/sjSDM/Y.csv`: every response column now has both 0 and 1 values. The previous example had an all-zero species column, which is invalid for a useful binomial/probit smoke test.
- Added `examples/sjSDM/run_real_example.R`, an executable smoke test that creates an output folder, writes `used_config.yml`, exports `reproducible_script/run_this_sjSDM_analysis.R`, runs it with Rscript and creates a ZIP.
- Connected the sjSDM GUI run path to real fitting when `sjSDM` is installed and `Run real sjSDM package fit` is enabled.
- Updated the generated sjSDM reproducible script so `step_size` is clamped to the number of rows.
- Fixed `coef.sjSDM` output handling for sjSDM 1.0.7, where `coef(model)` returns an unnamed list.
- Fixed standard association-table generation when covariance/correlation matrices do not carry row or column names.
- Fixed `getImportance` use by transposing the environmental coefficient matrix to the shape expected by sjSDM's internal `getImportance`.
- Replaced fragile package plot output with robust PDF heatmaps for environmental coefficients and species correlations.
- For non-spatial examples, `internalStructure` is now skipped with a warning because sjSDM only supports it for spatial models.

## Verified Real Output

Command run:

```bat
"C:\Program Files\R\R-4.5.3\bin\Rscript.exe" examples\sjSDM\run_real_example.R
```

Verified output folder:

```text
output/JSDMStudio_REAL_EXAMPLE_sjSDM_20260529_200621
```

Verified status:

```json
{
  "engine": "sjSDM",
  "status": "fitted",
  "errors": []
}
```

Real non-empty files were produced:

- `models/sjSDM_model.rds`
- `tables/coef_environment.csv`
- `tables/covariance_matrix.csv`
- `tables/correlation_matrix.csv`
- `tables/Rsquared.csv`
- `tables/residuals.csv`
- `predictions/predictions.csv`
- `importance/env_importance.csv`
- `importance/biotic_importance.csv`
- `anova/sjSDM_anova_summary.csv`
- `plots/sjSDM_plot.pdf`
- `standard/run_summary.csv`
- `standard/effects_long.csv`
- `standard/predictions_long.csv`
- `standard/associations_long.csv`
- `standard/fit_metrics.csv`
- `output/JSDMStudio_REAL_EXAMPLE_sjSDM_20260529_200621.zip`

## Remaining Note

The bundled smoke test is intentionally tiny and non-spatial. It is appropriate for proving that the software, package bridge, exported script and result-writing logic work. Formal ecological inference still requires larger data, more iterations/sampling and convergence/validation checks.
