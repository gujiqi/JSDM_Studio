# JSDM Studio

JSDM Studio is a Shiny/WebView2 workbench for reproducible joint species distribution modelling workflows.

This GitHub source package is the lightweight developer version. It keeps the app source, R adapters, workflow scripts, documentation, icons, examples, installer scripts and the lightweight Hmsc-HPC Python source needed by the CPU pyhmsc workflow. It intentionally excludes large generated outputs, compiled installers, historical backup files and external reference PDFs/ZIPs.

## Included

- `app.R`
- `R/` engine adapters and helpers
- `workflow_scripts/`
- `workflow_templates/`
- `examples/`
- `docs/`
- `www/` and `assets/`
- `installer/` build scripts, excluding `installer/output/`
- `webview2_launcher/` C# source, excluding `bin/` and `obj/`
- `external_packages/hmsc-hpc-main/` lightweight source, excluding large example datasets

## Not Included

- `output/` benchmark and model-run outputs
- `installer/output/JSDMStudio_Setup.exe`
- compiled launcher binaries and WebView2 DLLs
- historical `app_before_*` backup files
- external package reference papers and archived source ZIPs
- temporary logs and smoke-test files

The complete local release ZIP is preserved separately as:

`JSDMStudio_FINAL_COMPLETE_MAX_20260603_PER_FILE_UPLOAD_OK.zip`

## Run Locally

On Windows:

1. Install R.
2. Run `install_packages.bat`.
3. Optional for Hmsc-HPC: install Python 3.10+ and run `install_hmschpc_python_packages.bat`.
4. Start the app with `Launch_in_default_browser.bat`.

You can also run directly from R:

```r
shiny::runApp(".")
```

## Build Installer

The GitHub source package does not include the compiled installer. To rebuild it locally:

```bat
installer\build_installer.bat
```

The generated installer will be written to `installer/output/`, which is ignored by Git.

## Output Policy

Model runs and Universal Benchmark runs can create many files. They should stay out of Git and be shared as release artifacts or separate archives when needed.
