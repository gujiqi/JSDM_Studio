# JSDM Studio User Installation Guide

## Do users only need the EXE?

Not completely.

This version is:

```text
Windows Setup.exe
+ WebView2 desktop window
+ local R/Shiny/Hmsc backend
```

Users need:

```text
1. R for Windows
2. JSDMStudio_Setup.exe
3. Microsoft Edge WebView2 Runtime, usually already included in Windows 10/11
```

R packages are checked and installed automatically when possible.

Required R packages:

```text
Hmsc
shiny
bslib
DT
yaml
ggplot2
zip
coda
ape
corrplot
htmltools
```

If automatic installation fails, run:

```text
Repair R packages
```

from the Start Menu.

## Recommended user workflow

```text
1. Install R for Windows.
2. Install JSDMStudio_Setup.exe.
3. Launch JSDM Studio from the desktop.
4. Wait for package checking on first launch.
5. Use example data for the first test.
```

If the WebView2 window fails, open:

```text
JSDM Studio Browser Mode
```
