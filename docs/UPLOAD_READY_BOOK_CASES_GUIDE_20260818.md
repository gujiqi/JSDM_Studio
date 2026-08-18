# Upload-Ready Book Chapter Case Guide

Use this guide when you want to run the converted book chapter examples in the JSDM Studio GUI.

Clean entry folder:

`examples/00_UPLOAD_READY_book_chapter_cases`

Do not start from the larger `examples/book_chapter_benchmarks` folder unless you are inspecting the raw generated benchmark contract. The upload-ready folder is organized for manual GUI upload.

## Step 1. Choose One Case

- `ch06_plant_traits_whittaker`
- `ch07_deadwood_fungi`
- `ch11_finnish_birds`

## Step 2. Choose One Engine Folder

Inside each case folder:

- `01_Hmsc`
- `02_Hmsc-HPC`
- `03_jSDM`
- `04_GJAM`
- `05_spOccupancy`
- `06_sjSDM`
- `07_boral`

Upload files from only that engine folder into the matching workflow panel.

## Fast Default Uploads

| Workflow | Minimum files for quick run | Optional files for fuller test |
|---|---|---|
| Hmsc | `Y.csv`, `XData.csv` | `TrData.csv`, `studyDesign.csv`, `coordinates.csv`, `phylogeny.nwk`, `phylo_cov.csv` |
| Hmsc-HPC | `Y.csv`, `XData.csv` | `traits.csv`, `studyDesign.csv`, `coordinates.csv`, `phylo_cov.csv`, `newdata.csv` |
| jSDM | `Y.csv`, `XData.csv` | `trait_data.csv`, `long_format.csv`, `trials.csv`, `newdata.csv`, `prediction_ids.csv` |
| GJAM | `Y.csv`, `XData.csv`, `typeNames.csv` | `specByTrait.csv`, `traitTypes.csv`, `effort.csv`, `newdata.csv`, `holdoutIndex.csv`, `censor.csv` |
| spOccupancy | `y.csv`, `occ.covs.csv`, `det.covs.csv` | `coords.csv`, `species.csv`, `integrated_sources.csv`, `newdata.csv`, `newcoords.csv`, `folds.csv` |
| sjSDM | `Y.csv`, `env.csv` | `spatial.csv`, `traits.csv`, `newdata.csv` |
| boral | `Y.csv`, `XData.csv` | `traits.csv`, `row.ids.csv`, `ranef.ids.csv`, `distmat.csv`, `offset.csv`, `newdata.csv`, `trial.size.csv` |

## Full Slot Mapping

Each case folder contains:

- `README_UPLOAD_FIRST.md`
- `UPLOAD_GUIDE.csv`

Open `UPLOAD_GUIDE.csv` for exact GUI slot names and file names. The guide was validated against the actual files.

## Existing Verified Benchmark Runs

The converted chapter cases were run through all seven engines with the Universal Benchmark runner. The run audit is:

`output/book_chapter_benchmark_run_audit_20260818.csv`

That audit records the engine status, output ZIP existence and output contract completeness for each case.
