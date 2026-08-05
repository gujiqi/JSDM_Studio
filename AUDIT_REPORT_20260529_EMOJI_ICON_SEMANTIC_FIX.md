# JSDM Studio Icon Semantic Fix - 2026-05-29

## What was changed

- Replaced the automatic parameter-label icon system with colored semantic emoji icons in white rounded-square containers.
- Forced the key mappings requested by the user:
  - Y / species response matrix: 🌱
  - XData / environmental predictors: 🌡️
  - response family / distribution: 🧬
  - spatial / random effects: 🗺️
  - MCMC controls: ⏱️
  - diagnostics: 📈
  - model fit / validation: 🧪
  - output / ZIP: 📦
- Kept additional workflow semantics consistent:
  - traits: 🌿
  - phylogeny: 🌳
  - formulas: 🧮
  - prediction: 🔮
  - reports/scripts/tables: 📄
  - upload: 📤
  - preview: 👀
  - run: 🚀
  - save model: 💾
  - file list: 🗂️
- Changed CSS so parameter and command icons use a clean white background, soft border, subtle shadow and macOS-style rounded-square shape instead of transparent flat SVG line icons.

## Files changed

- `app.R`
  - Updated `.param-meaning-icon` and `.cmd-meaning-icon` styling.
  - Replaced the `param_icon_script` SVG filename registry with semantic colored emoji glyphs.
  - Improved label detection so `Y.csv`, `response matrix`, `species matrix` and `community matrix` map to the species-response icon, while `XData.csv` and environmental covariates map to the thermometer icon.

## Verification

- Parsed `app.R` successfully with R 4.5.3:
  - `C:\Program Files\R\R-4.5.3\bin\Rscript.exe`
- Started the Shiny app locally on `http://127.0.0.1:6883`; the server returned HTTP 200.
- Confirmed the served page includes the requested emoji icon glyphs and the white-background icon CSS.

## Follow-up cleanup

- Removed duplicate icons from the right-side `Value` label inside synchronized slider/numeric controls. Each synchronized parameter card now shows only one semantic icon.
- Removed the repeated HMSC input-map banner from the HMSC data-upload accordion because the HMSC workflow already has the same large HMSC hierarchy image above it.

## Remaining note

This update focuses on the current icon semantics and visual container issue. It does not replace the larger workflow-output and model-fitting audit already tracked in the earlier audit reports.
