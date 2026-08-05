@echo off
chcp 65001 > nul
cd /d "%~dp0"

echo ==========================================
echo JSDM Studio WebView2 Debug Launch
echo ==========================================
echo.

if not exist JSDMStudioLauncher.exe (
  echo Launcher missing. Building...
  call Build_WebView2_Launcher.bat
)

echo Starting launcher in debug mode...
echo Do not close this window.
echo.

"%CD%\JSDMStudioLauncher.exe"

echo.
echo Launcher exited.
echo.

if exist startup_log.txt (
  echo Opening startup_log.txt...
  start notepad startup_log.txt
) else (
  echo startup_log.txt was not created.
)

pause
