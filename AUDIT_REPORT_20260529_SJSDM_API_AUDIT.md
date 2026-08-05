# sjSDM API Audit - 2026-05-29

## Reference Basis

Checked against the local sjSDM 1.0.7 package and the supplied reference files:

- `C:/Users/Google/Downloads/s-jSDM参考/sjSDM_1.0.7.tar.gz`
- `C:/Users/Google/Downloads/s-jSDM参考/sjSDM_ Getting started with sjSDM - a scalable joint Species Distribution Model.mhtml`
- `C:/Users/Google/Downloads/s-jSDM参考/Help for package sjSDM.mhtml`
- Installed R package: `sjSDM 1.0.7` under R 4.5.3

The audited core function signature is:

```r
sjSDM(
  Y, env, biotic, spatial, family, iter, step_size,
  learning_rate, se, sampling, parallel, control,
  device, dtype, seed, verbose
)
```

## Fixed Issues

- Corrected the response-family UI mapping:
  - `negative_binomial_log` was replaced by sjSDM's actual `nbinom` option.
  - Exported script maps families to `binomial("probit")`, `binomial("logit")`, `poisson("log")`, `"nbinom"`, or `gaussian("identity")`.
- Corrected environmental module logic:
  - Previous `none` option was misleading because sjSDM requires an `env` design matrix.
  - It is now `intercept-only`, and the script creates `linear(data.frame(intercept=1), ~1)`.
- Corrected device options:
  - Removed `cuda` as a UI option because sjSDM 1.0.7 uses `device = "cpu"` or `device = "gpu"` / numeric GPU id.
- Removed invalid optimizer option:
  - `DiffGrad` exists in source but is not exported in sjSDM 1.0.7, so it is no longer selectable.
- Added missing real sjSDM controls:
  - `parallel`
  - `dtype`
  - `verbose`
  - `sjSDMControl$scheduler`
  - `sjSDMControl$lr_reduce_factor`
  - `sjSDMControl$early_stopping_training`
  - `sjSDMControl$mixed`
- Clarified `step_size`:
  - It is now labelled as `step_size / batch size`, because sjSDM passes it internally as the stochastic-gradient batch size.
  - The previous independent `batch_size` control was removed.
- Added `generateSpatialEV threshold` and stricter validation:
  - `generateSpatialEV` now requires a spatial module and at least two numeric coordinate columns.
- Improved validation:
  - Checks sjSDM family values against actual supported families.
  - Checks formula variables in env/spatial data.
  - Converts character predictors to factor for model-matrix encoding.
  - Checks GPU, dtype, optimizer, scheduler and mixed-precision consistency.
- Replaced the generic reproducibility stub with an executable sjSDM script:
  - `reproducible_script/run_this_sjSDM_analysis.R`
  - Builds `linear`, `DNN`, `bioticStruct`, `sjSDMControl`, `generateSpatialEV`, and `sjSDM()` calls from `used_config.yml`.
  - Writes model object, coefficients, covariance/correlation, predictions, R2, residuals, ANOVA/internal structure, importance and plots when the package environment supports them.
- Added output API mapping:
  - `tables/sjSDM_api_mapping.csv`
  - This explains how each GUI setting maps to sjSDM 1.0.7 API calls.
- Removed the duplicate sjSDM input-map image from the upload accordion; only the main workflow banner remains.

## Verified

- `app.R` parses successfully under R 4.5.3.
- The generated sjSDM reproducible script parses successfully.
- A minimal installed-package fit completed:

```r
sjSDM(
  Y = small_binary_matrix,
  env = linear(small_env, ~ .),
  family = binomial("probit"),
  iter = 2L,
  sampling = 100L,
  step_size = 2L,
  device = "cpu"
)
```

## Remaining Risks

- The GUI run still creates a safe scaffold by default; real sjSDM fitting is exported through the reproducible script.
- GPU mode depends on a working PyTorch/CUDA installation and should be checked outside Shiny before large runs.
- ANOVA/internal-structure and importance can be unstable for many species with few detections or small sample sizes, as noted in sjSDM documentation.
