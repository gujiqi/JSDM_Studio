# JSDM Studio

JSDM Studio is a Windows desktop workbench for reproducible joint species distribution modelling. It provides guided workflows for Hmsc, Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM, boral, and a Universal Benchmark workflow.

> [!IMPORTANT]
> ## Windows Full Offline Installation (recommended)
>
> Download the **three files** from the [latest GitHub Release](https://github.com/gujiqi/JSDM_Studio/releases/latest) and keep them in the **same folder**:
>
> 1. `JSDMStudio_Full_Offline_Setup.exe`
> 2. `JSDMStudio_Full_Offline_Setup-1.bin`
> 3. `JSDMStudio_Full_Offline_Setup-2.bin`
>
> Then right-click **`JSDMStudio_Full_Offline_Setup.exe`** and choose **Run as administrator**. Do not open, rename, or move either `.bin` file: they are installer data volumes, not programs. The installer reads them automatically.
>
> This is the Windows x64 offline edition. It installs JSDM Studio, R 4.5.3, JAGS, Microsoft Visual C++ Redistributable, Python 3.10, WebView2, .NET Desktop Runtime 5, bundled R packages, CPU torch/libtorch, and the Hmsc-HPC CPU environment. No package download is required during installation. Allow at least 8 GB of free disk space and several minutes for the first installation.

### Verify the installation

Start JSDM Studio from the Start menu or optional desktop shortcut. The final dependency report is written under:

```text
%LOCALAPPDATA%\JSDMStudio\diagnostics
```

An unsuccessful prerequisite check makes the setup fail clearly; it is not reported as a completed installation. For detailed troubleshooting, see [Windows Full Offline Installation](docs/WINDOWS_FULL_OFFLINE_INSTALL.md).

### Release integrity

The current Full Offline Edition files have these SHA256 checksums:

| File | SHA256 |
| --- | --- |
| `JSDMStudio_Full_Offline_Setup.exe` | `83EF525D32AAA427FA23F032ADD46AB02DE91A62E5FA3D3A5F90B2568DE148A9` |
| `JSDMStudio_Full_Offline_Setup-1.bin` | `DDD79B4CCA7A71CD80D6A43CE4BC1774FDAD33D43BBF82C596274C398D355086` |
| `JSDMStudio_Full_Offline_Setup-2.bin` | `FD38B1C5FD4D3707102FCE5305EC9296ADD771B9F31349DF123B10630CCCA62B` |

### Included workflows

* **Hmsc**: hierarchical JSDMs with traits, phylogeny, random effects, spatial structure, predictions, and diagnostic outputs.
* **Hmsc-HPC**: a CPU-oriented Python Hmsc-HPC workflow with reproducible Python and HDF5 outputs.
* **jSDM, GJAM, spOccupancy, sjSDM, and boral**: guided engine-specific inputs, diagnostics, standardized results, scripts, reports, and ZIP exports.
* **Universal Benchmark**: a shared latent ecological benchmark that generates engine-specific data and records which standardized results are comparable.

### Source and documentation

The repository contains source code, workflow scripts, examples, documentation, and installation material. Examples and generated outputs must be interpreted with their engine-specific assumptions and diagnostics; a completed workflow is not, by itself, evidence of statistical convergence or ecological adequacy.

For a local source checkout, the launch scripts are `Launch_JSDMStudio.bat` and `Start_WebView2_Window.bat`. The Full Offline Edition above is the recommended route for a new Windows computer.
