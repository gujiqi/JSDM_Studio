# JSDM Studio Full Offline Installation for Windows

Use this edition when JSDM Studio must run on a Windows x64 computer without downloading R packages or Python packages from the internet.

## Install

1. Open the [latest release](https://github.com/gujiqi/JSDM_Studio/releases/latest).
2. Download `JSDMStudio_Full_Offline_Setup.exe` and every file named `JSDMStudio_Full_Offline_Setup-*.bin`.
3. Keep all of those files in the same folder. Do not rename the `.bin` files.
4. Right-click `JSDMStudio_Full_Offline_Setup.exe` and choose **Run as administrator**.
5. Follow the wizard. It can take several minutes while R packages, the CPU torch runtime and the Hmsc-HPC Python environment are installed.
6. Start JSDM Studio from the Start menu or the optional desktop shortcut.

The Full Offline Edition installs runtime requirements only: R 4.5.3, JAGS, Microsoft Visual C++ Redistributable, Python 3.10, WebView2, .NET Desktop Runtime 5, bundled R packages, CPU torch/libtorch and the Hmsc-HPC CPU Python environment. It does not install Rtools, the .NET SDK or Inno Setup.

Allow at least 8 GB of free disk space. The final validation report is written under:

```text
%LOCALAPPDATA%\JSDMStudio\diagnostics
```

A failed prerequisite step is reported as a failed setup rather than a completed installation.

## Integrity checks

The release page and its accompanying `JSDMStudio_Full_Offline_Setup_RELEASE.txt` provide SHA256 values for every required file. Verify them before installation when your institutional policy requires integrity checks.
