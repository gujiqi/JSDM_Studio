# HMSC categorical predictor fix

This build fixes the HMSC error:

`XData variables had bad types: substrate:character`

Cause:
Hmsc accepts categorical predictors as factors, but the CSV reader was importing text columns as character. The example variable `substrate` is a categorical predictor and must be converted to factor before calling `Hmsc()`.

Fix:
- Character columns in XData are converted to factors before Hmsc fitting.
- Character columns in TrData and studyDesign are also converted where appropriate.
- Generated `workflow_scripts/S1_define_models.R` includes the same conversion, so the exported S1-S7 scripts are reproducible.
- HMSC data check now reports categorical predictors that will be converted.

You can keep `substrate` in XFormula:
`~ pH + moisture + canopy + substrate + elevation`
