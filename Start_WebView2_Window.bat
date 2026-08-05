@echo off
chcp 65001 > nul
cd /d "%~dp0"

if not exist JSDMStudioLauncher.exe (
  echo JSDMStudioLauncher.exe not found. Building it first...
  call Build_WebView2_Launcher.bat
)

if not exist JSDMStudioLauncher.exe (
  echo Launcher build failed.
  pause
  exit /b 1
)

if not exist Microsoft.Web.WebView2.WinForms.dll (
  echo WARNING: WebView2 dependency DLL not found in this folder.
  echo Rebuilding and copying all launcher dependencies...
  call Build_WebView2_Launcher.bat
)

echo Starting JSDM Studio WebView2 window...
echo If it closes immediately, open startup_log.txt in this folder.
echo.

start "" "%CD%\JSDMStudioLauncher.exe"

timeout /t 3 > nul
if exist startup_log.txt (
  echo.
  echo startup_log.txt exists. If the app did not open, check this file:
  echo %CD%\startup_log.txt
)

pause
