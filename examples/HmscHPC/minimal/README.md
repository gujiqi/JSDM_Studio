# Hmsc-HPC minimal example

Use these files in the Hmsc-HPC Workflow panel.

Quick fitted test:
- Upload `Y.csv` and `XData.csv`.
- Set `distribution = poisson`.
- Set `X formula = ~ forest_cover + elevation + moisture + C(substrate)`.
- Keep random level design as `none`.
- Use small sampler settings such as samples 3, transient 3, thin 1 and verbose 3.

Feature tests:
- Traits: upload `traits.csv`, enable Use traits, and use `~ body_size + forest_specialist + C(life_form)`.
- Phylogeny covariance: upload `phylo_cov.csv` and choose covariance mode.
- Newick: upload `phylo_tree.nwk` and choose Newick mode.
- IID random level: upload `studyDesign.csv`, choose iid and grouping column `plot`.
- Spatial full random level: upload `studyDesign.csv` and `coordinates.csv`, choose spatial_full and coordinate columns `xcoord`, `ycoord`.
- New site predictions: upload `newdata.csv`; fitted runs export `predictions/newdata_predicted_mean.csv`.
