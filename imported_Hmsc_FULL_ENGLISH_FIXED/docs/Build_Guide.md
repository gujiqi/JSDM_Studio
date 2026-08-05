# Build Guide: Generate JSDMStudio_Setup.exe from this ZIP

This guide is for the developer, not end users.

## Developer requirements

Install:

```text
R for Windows
Inno Setup
.NET SDK
```

## Build steps

After extracting the ZIP:

```text
1. Open the JSDMStudio folder.
2. Run install_packages.bat.
3. Run Build_WebView2_Launcher.bat.
4. Run Start_WebView2_Window.bat to test the desktop window.
5. Open JSDMStudio\installer.
6. Run build_installer.bat.
7. Find the final installer in JSDMStudio\installer\output.
```

Final file:

```text
JSDMStudio_Setup.exe
```

This is the file to distribute to users.

## Start menu entries

The installer creates:

```text
JSDM Studio - Window Mode
JSDM Studio - Browser Mode
Repair R packages
Open startup log
Troubleshoot R
Uninstall JSDM Studio
```
