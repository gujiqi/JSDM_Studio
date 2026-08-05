# JSDM Studio

JSDM Studio is an Hmsc-powered graphical workbench for joint species distribution modelling.

It provides data upload, model specification, MCMC settings, output selection, progress display, automatic result export, diagnostics, and reports.

## Build the installer

1. Open the JSDMStudio folder.
2. Run `install_packages.bat`.
3. Run `Build_WebView2_Launcher.bat`.
4. Run `Start_WebView2_Window.bat` to test the desktop window.
5. Open `JSDMStudio\installer`.
6. Run `build_installer.bat`.
7. Find `JSDMStudio_Setup.exe` in `installer\output`.

## User requirements

Users need R for Windows. Microsoft Edge WebView2 Runtime is usually already included in Windows 10/11. Required R packages are checked and installed automatically when possible.
