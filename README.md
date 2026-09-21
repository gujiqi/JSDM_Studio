# JSDM Studio

JSDM Studio is a Windows desktop workbench for reproducible joint species distribution modelling. It provides guided workflows for Hmsc, Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM, boral, and a Universal Benchmark workflow.

> [!IMPORTANT]
> ## Windows Full Offline Installation (recommended)
>
> Download **one file** from the [latest GitHub Release](https://github.com/gujiqi/JSDM_Studio/releases/latest):
>
> `JSDMStudio_Full_Offline_Setup.exe`
>
> Then right-click it and choose **Run as administrator**. It is a self-contained offline installer; no companion `.bin` file is needed.
>
> This is the Windows x64 offline edition. It installs JSDM Studio, R 4.5.3, JAGS, Microsoft Visual C++ Redistributable, Python 3.10, WebView2, .NET Desktop Runtime 5, bundled R packages, the Hmsc-HPC CPU environment, and sjSDM's CPU PyTorch backend through reticulate. No package download is required during installation. Keep the installer open until it reports that dependency validation is complete. Allow at least 8 GB of free disk space and several minutes for the first installation.

### Verify the installation

Start JSDM Studio from the Start menu or optional desktop shortcut. The final dependency report is written under:

```text
C:\Program Files\JSDMStudio\diagnostics
```

An unsuccessful prerequisite check makes the setup fail clearly; it is not reported as a completed installation. For detailed troubleshooting, see [Windows Full Offline Installation](docs/WINDOWS_FULL_OFFLINE_INSTALL.md).

### Release integrity

The current Full Offline Edition has this SHA256 checksum:

| File | SHA256 |
| --- | --- |
| `JSDMStudio_Full_Offline_Setup.exe` | `9B5C1BE5EA3FC6EFC220A005A974B2BC6F2F6F72F6A974DFD17C745F0597BBEA` |

### Included workflows

* **Hmsc**: hierarchical JSDMs with traits, phylogeny, random effects, spatial structure, predictions, and diagnostic outputs.
* **Hmsc-HPC**: a CPU-oriented Python Hmsc-HPC workflow with reproducible Python and HDF5 outputs.
* **jSDM, GJAM, spOccupancy, sjSDM, and boral**: guided engine-specific inputs, diagnostics, standardized results, scripts, reports, and ZIP exports.
* **Universal Benchmark**: a shared latent ecological benchmark that generates engine-specific data and records which standardized results are comparable.

### macOS status

This release does **not** include a macOS `.dmg` or `.app`. The Windows installer cannot run on macOS. macOS users may inspect the source and use the Shiny app in a browser only after independently installing compatible R, Python, JAGS and engine dependencies; this route has not yet been packaged or validated as a supported one-click installation.

### Source and documentation

The repository contains source code, workflow scripts, examples, documentation, and installation material. Examples and generated outputs must be interpreted with their engine-specific assumptions and diagnostics; a completed workflow is not, by itself, evidence of statistical convergence or ecological adequacy.

For a local source checkout, the launch scripts are `Launch_JSDMStudio.bat` and `Start_WebView2_Window.bat`. The Full Offline Edition above is the recommended route for a new Windows computer.
