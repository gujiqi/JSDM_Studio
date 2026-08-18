# Maximum Upload Guide for Book Chapter Cases

Use this file when you want to upload every available input file for one workflow.

First choose one case folder:

- `ch06_plant_traits_whittaker`
- `ch07_deadwood_fungi`
- `ch11_finnish_birds`

Then open the engine folder listed below and upload every file in that folder to the corresponding JSDM Studio workflow panel.

## Hmsc

Folder: `01_Hmsc`

Upload all files:

- `Y.csv` -> Y.csv response matrix
- `XData.csv` -> XData.csv environmental predictors
- `TrData.csv` -> TrData.csv traits / response attributes
- `studyDesign.csv` -> studyDesign.csv grouping / random-effect design
- `coordinates.csv` -> coordinates.csv spatial coordinates
- `phylogeny.nwk` -> phylogeny / taxonomy file
- `phylo_cov.csv` -> C correlation matrix file name, when the advanced C/correlation option is used

Recommended full-test settings:

- Enable traits.
- Enable phylogeny only when using `phylogeny.nwk` or `phylo_cov.csv`.
- For grouping random effect, choose Sample/grouping and grouping column `plot`.
- For spatial random effect, choose Spatial Full, Spatial NNGP or Spatial GPP and use coordinate columns `x` and `y` if those names are shown in the file.

## Hmsc-HPC

Folder: `02_Hmsc-HPC`

Upload all files:

- `Y.csv` -> response matrix
- `XData.csv` -> environmental predictors
- `traits.csv` -> traits
- `studyDesign.csv` -> grouping / random-effect design
- `coordinates.csv` -> spatial coordinates
- `phylo_cov.csv` -> phylogenetic / taxonomic covariance matrix
- `newdata.csv` -> prediction covariates

Recommended full-test settings:

- Use CPU mode.
- Enable traits if the panel provides a traits option.
- Enable random/spatial options only one at a time unless the panel explicitly supports both.
- Use quick iterations for a GUI smoke test.

## jSDM

Folder: `03_jSDM`

Upload all files:

- `Y.csv` -> response matrix
- `XData.csv` -> environmental predictors
- `trait_data.csv` -> species traits
- `long_format.csv` -> long-format input, only if the selected jSDM branch uses long format
- `trials.csv` -> binomial trial sizes
- `newdata.csv` -> prediction covariates
- `prediction_ids.csv` -> prediction site/species ID filter

Recommended full-test settings:

- Use binomial/probit occurrence mode for the main benchmark.
- Use `trials.csv` only for a binomial-with-trials branch.
- Do not enable long-format input unless that branch is selected.

## GJAM

Folder: `04_GJAM`

Upload all files:

- `Y.csv` -> response matrix
- `XData.csv` -> predictor matrix
- `typeNames.csv` -> response type names
- `specByTrait.csv` -> species-by-trait table
- `traitTypes.csv` -> trait type definitions
- `effort.csv` -> effort / offset file
- `newdata.csv` -> prediction covariates
- `holdoutIndex.csv` -> holdout rows
- `censor.csv` -> censoring file

Recommended full-test settings:

- Keep the main benchmark response type as PA.
- Use trait files only in trait-enabled GJAM branches.
- Use censoring only if the panel explicitly enables censoring.

## spOccupancy

Folder: `05_spOccupancy`

Upload all files:

- `y.csv` -> replicated detection-nondetection response
- `occ.covs.csv` -> occupancy covariates
- `det.covs.csv` -> detection covariates
- `coords.csv` -> site coordinates
- `species.csv` -> species metadata
- `integrated_sources.csv` -> integrated-data source metadata
- `newdata.csv` -> prediction occupancy covariates
- `newcoords.csv` -> prediction coordinates
- `folds.csv` -> cross-validation folds

Recommended full-test settings:

- Use multi-species occupancy for the main benchmark.
- Detection covariates should be mapped to detection effects, not occurrence effects.
- Use spatial coordinates only for spatial model branches.

## sjSDM

Folder: `06_sjSDM`

Upload all files:

- `Y.csv` -> response matrix
- `env.csv` -> environmental predictors
- `spatial.csv` -> spatial coordinates or spatial predictors
- `traits.csv` -> species traits
- `newdata.csv` -> prediction covariates

Recommended full-test settings:

- Use occurrence/binomial-compatible response.
- Enable environmental module with `env.csv`.
- Enable spatial module only when using `spatial.csv`.
- Enable traits only when the selected sjSDM branch supports trait input.

## boral

Folder: `07_boral`

Upload all files:

- `Y.csv` -> response matrix
- `XData.csv` -> environmental predictors
- `traits.csv` -> species traits
- `row.ids.csv` -> row IDs
- `ranef.ids.csv` -> random-effect IDs
- `distmat.csv` -> distance matrix
- `offset.csv` -> offset matrix
- `newdata.csv` -> prediction covariates
- `trial.size.csv` -> binomial trial size matrix

Recommended full-test settings:

- Use binomial/probit-compatible occurrence for the main benchmark.
- Use `trial.size.csv` with binomial response.
- Use random-effect IDs only when random effects are enabled.
- Use `distmat.csv` only for distance/spatial branches supported by the selected boral settings.

## Validation

The upload guide files were checked programmatically. Every file listed above exists for all three case folders.

The three book-derived cases were also run through the Universal Benchmark runner. See:

`output/book_chapter_benchmark_run_audit_20260818.csv`
