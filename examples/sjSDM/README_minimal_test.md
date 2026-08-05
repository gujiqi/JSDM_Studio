# sjSDM Minimal Test

Use these files for the first sjSDM check:

- `Y.csv`
- `env.csv`

Recommended GUI settings:

- `family`: `binomial(probit)`
- `env module`: `linear`
- `env formula`: `~ .`
- `spatial module`: `none`
- `iter`: `5` for the bundled smoke test, larger for real analysis
- `sampling`: `100` for the bundled smoke test, larger for real analysis
- `step_size / batch size`: `4` for the bundled 12-row example; JSDM Studio clamps this to the number of rows in the exported script
- `parallel data-loader workers`: `0` on Windows
- `device`: `cpu`
- `dtype`: `float32`
- `optimizer`: `Adamax`
- `generateSpatialEV`: off for the first test
- `Run real sjSDM package fit`: on

Expected workflow:

1. Upload `Y.csv` and `env.csv`.
2. Click `Check sjSDM data`.
3. Click `Run sjSDM workflow`.
4. Open `diagnostics/engine_status.json` in the output folder.
5. To run the real package fit outside the scaffold, run:

```r
source("reproducible_script/run_this_sjSDM_analysis.R")
```

The exported script is mapped to sjSDM 1.0.7 and uses `sjSDM()`, `linear()`, `DNN()`, `bioticStruct()`, `sjSDMControl()`, `getCov()`, `getCor()`, `predict()`, `Rsquared()`, `anova()` and `internalStructure()` where requested.

Command-line smoke test:

```bat
"C:\Program Files\R\R-4.5.3\bin\Rscript.exe" examples\sjSDM\run_real_example.R
```

The smoke test creates a timestamped `output/JSDMStudio_REAL_EXAMPLE_sjSDM_*` folder and a real ZIP. The expected real outputs include `models/sjSDM_model.rds`, `tables/coef_environment.csv`, `tables/covariance_matrix.csv`, `tables/correlation_matrix.csv`, `predictions/predictions.csv`, `tables/residuals.csv`, `importance/env_importance.csv`, `anova/sjSDM_anova_summary.csv`, `plots/sjSDM_plot.pdf` and populated `standard/*.csv` comparison tables.

For this first non-spatial example, `internalStructure` is skipped with a diagnostic warning because sjSDM only supports that step for spatial models.

Multi-case real fitting suite:

```bat
"C:\Program Files\R\R-4.5.3\bin\Rscript.exe" examples\sjSDM\run_real_example_suite.R
```

The suite generates and fits seven small cases: binomial probit with categorical predictors, poisson with linear spatial predictors, gaussian with DNN environmental response, binomial logit with DNN/dropout/scheduler, nbinom intercept-only, binomial with generated spatial eigenvectors, and binomial with spatial DNN/madgrad. Each case must finish with `engine_status = fitted`, `exit_status = 0`, no missing required outputs, and a real ZIP file.
