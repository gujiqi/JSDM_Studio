# Chapter 7 dead wood-inhabiting fungi sequencing data

Open the engine folder that matches the JSDM Studio workflow panel.
Upload the files using the exact GUI slot names listed in UPLOAD_GUIDE.csv.

Fastest default tests:
- Hmsc: upload only Y.csv and XData.csv, keep probit, random level none.
- Hmsc-HPC: upload Y.csv and XData.csv, CPU quick settings.
- jSDM: upload Y.csv and XData.csv; trait_data.csv is optional.
- GJAM: upload Y.csv, XData.csv and typeNames.csv; keep PA response type.
- spOccupancy: upload y.csv, occ.covs.csv and det.covs.csv; coords optional for spatial models.
- sjSDM: upload Y.csv and env.csv; spatial.csv optional.
- boral: upload Y.csv and XData.csv; trial.size.csv optional for binomial.

For full benchmark coverage, upload every file in that engine folder to the matching optional slot.
