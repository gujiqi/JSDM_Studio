# jSDM Workflow Audit and Real-Run Fix Report

Date: 2026-05-29

## Scope

Audited and repaired the JSDM Studio `jSDM Workflow` against the local reference files in `C:/Users/Google/Downloads/jSDM参考`, the installed `jSDM` package API, and real synthetic runs.

Installed runtime used:

- R: `C:/Program Files/R/R-4.5.3/bin/Rscript.exe`
- jSDM: `0.2.7`
- Exported jSDM engines: `jSDM_binomial_probit`, `jSDM_binomial_logit`, `jSDM_poisson_log`, `jSDM_gaussian`, `jSDM_binomial_probit_long_format`, `jSDM_binomial_probit_sp_constrained`

## Fixed

1. Connected the jSDM workflow to real package fitting instead of scaffold-only output.
   - Added `R/jsdm_adapter.R`.
   - Added executable export: `reproducible_script/run_this_jSDM_analysis.R`.
   - Added workflow wrapper: `workflow_scripts/run_jSDM_workflow.R`.

2. Replaced basic jSDM checks with full parameter/data validation.
   - Validates MCMC constraints required by jSDM: `burnin + mcmc >= 100`, divisible by 10, and `mcmc %% thin == 0`.
   - Validates binary, count, and Gaussian response families.
   - Validates formulas against uploaded data.
   - Converts character predictors to factors and numeric-looking columns to numeric.
   - Validates trait dimensions and trait formulas.
   - Validates `newdata.csv` against prediction formulas.

3. Fixed `binomial_logit` trials handling.
   - jSDM 0.2.7 expects one trial count per site, applied to all species.
   - The adapter now accepts a scalar or one column/vector per site.
   - Species-specific trial matrices are rejected during check instead of failing inside jSDM.

4. Added a safety block for an unstable jSDM 0.2.7 combination.
   - `binomial_logit` with both `n_latent > 0` and `site_effect = fixed/random` is blocked at check time.
   - The message tells users to use latent variables with `site_effect = none`, or set `n_latent = 0` when using site effects.

5. Fixed long-format probit behavior.
   - jSDM 0.2.7 can fail if long-format formulas omit species terms.
   - The exported script rewrites unsafe long-format formulas to species-specific terms, e.g. `~ species + species:x1 + species:x2`.

6. Completed output writing.
   - Real fitted model: `models/jsdm_model.rds`, `models/jsdm_model_primary.rds`.
   - MCMC: `mcmc/mcmc_sp.rds`, `mcmc/mcmc_latent.rds`, `mcmc/mcmc_alpha.rds`, `mcmc/mcmc_V_alpha.rds`, `mcmc/mcmc_V.rds`, `mcmc/mcmc_Deviance.rds`.
   - Tables: posterior summaries, model spec, residual/env correlations where available.
   - Standard comparison outputs: `standard/run_summary.csv`, `effects_long.csv`, `predictions_long.csv`, `associations_long.csv`, `fit_metrics.csv`.
   - Results folder now contains real summaries and a README.
   - Diagnostics: `diagnostics/engine_status.json`, `diagnostics/data_check_messages.csv`, `diagnostics/session_info.txt`, and real error files on failure.
   - ZIP creation remains real and non-empty through `make_zip()`.

7. Added automated real-run coverage.
   - New script: `examples/jSDM/run_real_example_suite.R`.
   - The script generates small synthetic datasets, runs real jSDM fits, checks required outputs, creates ZIPs, and fails if any case is incomplete.

## Real Run Results

Latest validation summary:

`output/jSDM_real_example_suite_summary_20260529_235800.csv`

All 7 cases passed with `status = fitted`:

| Case | Model branch | Main coverage |
|---|---|---|
| 01 | `binomial_probit` | traits, categorical XData, latent variables, random site effects, newdata, prediction IDs |
| 02 | `binomial_logit` | trials.csv, fixed site effects |
| 03 | `poisson_log` | count response, latent variables, random site effects |
| 04 | `gaussian` | continuous response, traits, fixed site effects |
| 05 | `binomial_probit_long_format` | long-format upload, formula rewrite, species-specific terms |
| 06 | `binomial_probit_sp_constrained` | constrained latent probit, 2 chains, Windows-safe `ncores = 1` |
| 07 | `binomial_probit` | no latent effects, output switch coverage, disabled prediction/figures still leaves non-empty folders |

Additional check:

- No required output file was missing.
- No required output file was empty.
- No core output directory was empty across the latest 7 cases.

## Still Real Risks

1. These are quick smoke tests, not publication MCMC settings. Real analyses need larger `burnin`, `mcmc`, convergence checks, and ecological review.
2. `binomial_logit + n_latent > 0 + fixed/random site_effect` is blocked because it fails inside jSDM 0.2.7 on the tested runtime.
3. `binomial_probit_sp_constrained` is slower and more fragile than ordinary probit. It requires at least `n_latent + 2` species and uses `ncores = 1` for Windows safety.
4. Prediction with latent variables is safest when prediction site IDs match training site IDs, as in the jSDM reference example.
5. jSDM does not require JAGS, Python, PyTorch, or CUDA, but it does require the R package `jSDM` and companion packages such as `coda`, `yaml`, `jsonlite`, and `zip`.

## Minimum Reproducible Test

From the JSDM Studio folder:

```powershell
& 'C:\Program Files\R\R-4.5.3\bin\Rscript.exe' examples\jSDM\run_real_example_suite.R
```

Expected result:

- 7 rows printed.
- Every row has `status = fitted` and `pass = TRUE`.
- A new `output/jSDM_real_example_suite_summary_*.csv` is written.

## If a jSDM Run Fails

Open these first:

1. `diagnostics/engine_status.json`
2. `diagnostics/data_check_messages.csv`
3. `diagnostics/jSDM_reproducible_error.txt`
4. `diagnostics/jSDM_real_fit_stdout_stderr.txt` or `diagnostics/suite_stdout_stderr.txt`
5. `diagnostics/session_info.txt`

