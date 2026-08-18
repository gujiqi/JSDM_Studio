# Shiny UI/Server Strict Audit

Date: 2026-08-11

## Scope

This audit checked the current `app.R` and Shiny runtime contract for JSDM Studio:

- every UI `inputId` referenced by server code
- every `input$...` reference backed by a UI control
- every UI `output` backed by a render or download handler
- every `downloadButton` backed by a `downloadHandler`
- ZIP/CSV download paths creating real files
- `app.R` parse/source startup
- local HTTP smoke startup

## Results

- `app.R` parse/source: passed.
- Strict Shiny binding audit: passed.
- UI input controls audited: 533.
- Server input references audited: 533.
- Missing UI inputs: 0.
- UI outputs audited: 87.
- Output handlers audited: 87.
- Outputs without handlers: 0.
- Download buttons audited: 10.
- Download handlers audited: 10.
- Download buttons without handlers: 0.
- Static bad download content blocks: 0.
- R backend files parsed: 8.
- Compare Models test: passed.
- HTTP smoke test: passed with HTTP 200 and expected page markers.

## Download Artifact Check

The following ZIP-style download paths were tested by forcing the same diagnostic fallback used when a user downloads before running a workflow:

- Hmsc
- Hmsc-HPC
- jSDM
- GJAM
- spOccupancy
- sjSDM
- boral
- Universal Benchmark

Each generated a real non-empty ZIP containing at least:

- `diagnostics/engine_status.json`
- `diagnostics/data_check_messages.csv`
- `standard/run_summary.csv`

The parameter dictionary CSV download path was also tested and wrote a real CSV with 536 parameter rows.

## Fix Status

No new Shiny UI/server binding defects were found in this pass. No `app.R` remediation was required.

## Remaining Risk

This audit verifies Shiny wiring and download behavior. It does not prove that every statistical engine will fit on every machine; fitted status still depends on installed engine dependencies such as Hmsc, pyhmsc/Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM/PyTorch and boral/JAGS.

