# Book Chapter Benchmarks for JSDM Studio

These examples are derived from the book example data for *Joint Species Distribution Modelling: With Applications in R* and converted into the JSDM Studio Universal Benchmark input contract.

## Included Cases

- `ch06_plant_traits_whittaker`: Chapter 6 plant abundance, leaf C:N trait and Whittaker topographic moisture gradient.
- `ch07_deadwood_fungi`: Chapter 7 dead wood-inhabiting fungal sequencing counts, decay class and sequencing read count.
- `ch11_finnish_birds`: Chapter 11 Finnish bird presence-absence data with habitat, climate, coordinates, traits and phylogeny.

## Conversion Contract

Each case includes the same JSDM Studio benchmark files: `Y_occurrence.csv`, `Y_count.csv`, `Y_normal.csv`, `XData.csv`, `traits.csv`, `studyDesign.csv`, `coordinates.csv`, `phylogeny.nwk`, `phylo_cov.csv`, `newdata.csv`, `newcoords.csv`, `folds.csv`, `trial.size.csv`, `offset.csv`, `row.ids.csv`, `ranef.ids.csv`, `distmat.csv`, `spOccupancy_y_detection.csv`, `occ.covs.csv`, `det.covs.csv` and the `truth/` reference summaries.

The spOccupancy detection-nondetection replicates are deterministic pseudo-replicates derived from observed occurrence because the original chapter data are not replicated occupancy surveys. These files are suitable for software execution and cross-engine output validation, not for drawing biological conclusions about detection probability.

## Reproducibility

Run this from the JSDM Studio root to regenerate the three cases:

```r
source("workflow_scripts/prepare_jsdm_book_chapter_cases.R")
```

The raw book data copied from the local source folder are kept under `raw_source/`.
