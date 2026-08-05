# Hmsc-HPC parameter and guide re-audit

Date: 2026-05-30

This re-audit checked the Hmsc-HPC workflow from the GUI boundary through the R adapter, Python runner, output folders, comparison tables and user-facing navigation pages.

## Hmsc-HPC parameter audit

The visible Hmsc-HPC parameters are now either used directly, validated before fitting, or explicitly guarded:

- Input files: `Y.csv`, `XData.csv`, `traits.csv`, `newdata.csv`, `studyDesign.csv`, `coordinates.csv`, `phylo_cov.csv`, `phylo_tree.nwk`.
- Model settings: distribution, X formula, trait switch/formula, phylogeny mode, random-level design, random-level name, grouping column, coordinate columns, random-slope formula, `nf`, `nfMin`, `nfMax`, spatial `alpha`.
- Sampler settings: run sampler, samples, transient, thin, chains, chains-to-run, verbose, rngseed, precision, truncated-normal backend, HMC leapfrog/thin, updbe, save Eta, TensorFlow eager/profile, Python executable and Hmsc-HPC source directory.
- Output settings: predictions, diagnostics, plots and ZIP.

The audit found and fixed these concrete issues:

- `Predictions` was visible in the UI but did not actually control prediction export. It now does.
- Fitted Hmsc-HPC runs now write a Hmsc-style S1-S7 result chain in `results/`, so the downloaded ZIP contains model definition, fit summary, convergence diagnostics, model fit, model-fit display summaries, parameter estimates and predictions.
- Random-level fitted cases now export Eta/Lambda summaries and lambda-derived species association matrices instead of only a generic placeholder association table.
- Diagnostics now include `data_check_messages.csv`, `session_info.txt`, `posterior_hdf5_summary.csv` and `standard/diagnostics_long.csv`.
- Hmsc-HPC `standard/predictions_long.csv` used a narrower schema than other engines. It now includes `engine`, `site_id`, `response_id`, `observed`, `predicted_mean`, `predicted_lower`, `predicted_upper` and `prediction_set`.
- `save Eta = FALSE` caused the current Hmsc-HPC HDF5 exporter to fail after sampling. The runner now forces Eta export for HDF5 and records the reason in diagnostics.
- `app.R` was missing one UI-closing parenthesis after guide edits. The Shiny source now parses again.
- The Hmsc-HPC Python executable text input incorrectly placed `help_text()` inside `textInput()` as an unnamed width argument. It now sits beside the input inside the column, so the UI can be constructed cleanly.
- Validation now checks formula variables, response distribution/data compatibility, phylogenetic covariance shape/names/symmetry, random-level requirements, spatial coordinate numeric validity, `nfMin <= nf <= nfMax`, chain ids, sampler counts, Python path/source and output toggles.

## Case coverage

The current real example suite is `examples/HmscHPC/run_real_example_suite.R`.

Latest run summary:

- Summary CSV: `output/HmscHPC_real_example_suite_summary_20260530_141510.csv`
- Fitted cases: fixed poisson with categorical covariate, fixed probit, normal/HMC settings, traits + phylogenetic covariance, iid random intercept, spatial_full random intercept, Newick phylogeny, gaussian alias + tfd + eager/profile + predictions off.
- Guarded model-defined case: random_slope_iid compile-only.
- Missing required files: none.
- Empty output directories: none.
- Fitted Hmsc-HPC case ZIPs now contain 82-100 files each, including `results/S1_*` through `results/S7_*`, posterior HDF5 summaries, standard comparison tables and richer parameter outputs.

## Navigation tab audit

Home:

- Updated to list all seven engines and show Hmsc/Hmsc-HPC as adjacent but separate.
- Kept the core message: upload/preview, check assumptions, fit by engine, export scripts, compare later.

Project:

- Expanded question templates for occupancy, traits, phylogeny, mixed scales, large eDNA/OTU matrices, inverse prediction, sparse effects and fourth-corner/ordination use cases.
- Updated engine recommendation logic so Hmsc-HPC is recommended for traits, phylogeny, spatial_full random intercepts, predictions and larger species matrices where a CPU pyhmsc/HDF5 path is useful.

Engine Guide:

- Added a separate Hmsc-HPC engine card beside Hmsc.
- Clarified that the GUI exposes CPU pyhmsc/HDF5 workflow support, not GPU or Slurm submission.
- Clarified that GPP/NNGP and full Hmsc-R variance partitioning remain classic Hmsc territory.

Compare Models:

- Now includes Hmsc-HPC folder input and comparison status.
- Clarifies that Hmsc Omega, Hmsc-HPC Eta/Lambda samples, jSDM residual correlations, GJAM covariance, spOccupancy latent factors, sjSDM covariance and boral residual correlations are related but not identical.
- Reads `standard/` outputs from each engine rather than assuming one common model object.

Guides:

- Added a Hmsc-HPC guide with parameter meanings.
- Added a Hmsc-HPC result workflow / downloaded ZIP guide.
- Replaced an old Hmsc/jSDM-only separation statement with the complete seven-engine statement.
- Updated the output guide to include all engine output folder patterns.

## Remaining risks

- Hmsc-HPC in JSDM Studio is CPU-only. It is not a GPU/Slurm launcher.
- Random-slope native sampling is guarded; it compiles and validates but does not claim fitted posterior samples.
- Quick synthetic cases prove code paths and output integrity, not MCMC convergence.
- The Hmsc-HPC Python stack must be installed in the selected Python environment.
