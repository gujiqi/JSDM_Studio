# JSDM Studio Explanatory Interface Update

Date: 2026-08-11

This update strengthens JSDM Studio as a reviewer-oriented ecological modelling workbench. It draws design inspiration from clear single-engine teaching interfaces, while keeping JSDM Studio's original multi-engine architecture and avoiding copied wording.

## What Changed

- Updated the Home cover page to foreground review-grade explanations, separate engine contracts, the Universal Benchmark, executable scripts and auditable outputs.
- Added a compact parameter-literacy ribbon that explains the response matrix, predictors, random/spatial structure, MCMC diagnostics and reproducibility ZIP contract.
- Added an Engine Guide review layer that asks users to choose engines by scientific question, response scale, detection process, spatial structure, computational constraints and interpretation target.
- Added model-specific interpretation notes at the top of Hmsc, Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM, boral and Universal Benchmark workflows.
- Added a new Guides panel called "Reviewer-style parameter explanation layer" with original explanations for common controls and failure modes.
- Added an "Engine-specific interpretation boundaries" guide clarifying what each engine can and cannot support scientifically.

## Design Principles

JSDM Studio should not be only a GUI wrapper. For software review, each workflow must make three things explicit:

1. The data contract: what input shape, response scale and ID matching are required.
2. The ecological claim: what the fitted parameter or output can support.
3. The reproducibility contract: what script, diagnostic and ZIP output must exist after a run.

## Non-Copying Note

This update does not copy the reference HMSC Studio text or implementation. It adopts the broader strengths of explanatory interface design: clear method boundaries, parameter-level teaching, explicit limitations and reviewer-friendly language.

## Remaining Scope

The core fitting adapters were intentionally left unchanged in this update. Any future algorithmic changes should be tested separately with parse/source checks, strict Shiny input/output audits, HTTP smoke tests and workflow-specific synthetic cases.
