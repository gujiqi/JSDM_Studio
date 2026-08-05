
# JSDMWorkbench Hmsc + jSDM update

This ZIP is based on the original `JSDMStudio_DETAILED_GUIDE` package and adds a new jSDM engine.

Key changes:

1. The app is reorganized as an engine-based workflow.
2. All original Hmsc settings are kept in the `Hmsc Settings` tab.
3. A new `jSDM Settings` tab exposes jSDM model type, latent variables, site effects, MCMC settings, priors, starting values and correlation outputs.
4. Uploaded jSDM package materials are stored in `external_packages/jSDM/`.
5. The original Hmsc-only `app.R` is backed up as `app_Hmsc_original_backup.R`.
6. `R/engines/engine_jsdm.R` implements the jSDM adapter.
7. `R/engine_registry.R` defines available engines.

Recommended first test:

```text
1. Run install_packages.bat
2. Run run_app.bat or Start_WebView2_Window.bat
3. Upload examples/Y.csv and examples/XData.csv
4. Select only Hmsc first and run a quick test
5. Then select jSDM with binomial_probit, n_latent=2, site_effect=random, burnin=100, mcmc=100
```

jSDM installation note: the package uses C++ and GNU GSL through RcppGSL. If CRAN binary installation fails, users may need Rtools/GSL or the included source package in `external_packages/jSDM/`.
