@echo off
chcp 65001 > nul
setlocal EnableExtensions

set "ORIG_ROOT=%~dp0.."
for %%I in ("%ORIG_ROOT%") do set "ORIG_ROOT=%%~fI"

cd /d "%ORIG_ROOT%"

echo Building WebView2 launcher first...
call Build_WebView2_Launcher.bat
if not exist "JSDMStudioLauncher.exe" (
  echo ERROR: JSDMStudioLauncher.exe was not built.
  pause
  exit /b 1
)

REM WebView2 creates a runtime cache folder next to the launcher when the app
REM has been opened from the source tree. Do not package that cache; it can
REM contain very long paths and causes Inno Setup to fail on some Windows builds.
if exist "JSDMStudioLauncher.exe.WebView2" (
  echo Removing generated WebView2 runtime cache before packaging...
  rmdir /S /Q "JSDMStudioLauncher.exe.WebView2"
)

echo ==========================================
echo Building JSDMStudio_Setup.exe
echo ==========================================
echo.

set "ISCC="

where ISCC.exe >nul 2>nul
if %ERRORLEVEL%==0 (
  for /f "delims=" %%i in ('where ISCC.exe') do (
    set "ISCC=%%i"
    goto :found
  )
)

if exist "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" set "ISCC=C:\Program Files (x86)\Inno Setup 6\ISCC.exe"
if exist "C:\Program Files\Inno Setup 6\ISCC.exe" set "ISCC=C:\Program Files\Inno Setup 6\ISCC.exe"

:found
if "%ISCC%"=="" (
  echo ERROR: Inno Setup Compiler was not found.
  echo.
  echo Please install Inno Setup from:
  echo https://jrsoftware.org/isinfo.php
  echo.
  pause
  exit /b 1
)

echo Found Inno Setup Compiler:
echo "%ISCC%"
echo.

REM Build from a short staging path. The complete developer bundle contains
REM nested Hmsc-HPC reference files; compiling directly from a deep extracted
REM folder can hit Windows/Inno path handling limits.
set "STAGE_ROOT=%TEMP%\JSDMStudio_InnoStage"
if exist "%STAGE_ROOT%" rmdir /S /Q "%STAGE_ROOT%"
mkdir "%STAGE_ROOT%" >nul 2>nul

echo Staging source tree to short path:
echo "%STAGE_ROOT%"
robocopy "%ORIG_ROOT%" "%STAGE_ROOT%" /E ^
  /XD "JSDMStudioLauncher.exe.WebView2" "launcher_build" "output" "diagnostics" "input" "models" "tables" "plots" "predictions" "report" "results" "chains" "missing_data" "ordination" "sensitivity" "standard" "bin" "obj" "__pycache__" ^
  /XF "*.log" "*.err" "*.out" "*.tmp" "*.pyc" "*.pyo" "tmp_*" ^
  /NFL /NDL /NJH /NJS /NP >nul
set "ROBOCOPY_RC=%ERRORLEVEL%"
if %ROBOCOPY_RC% GTR 7 (
  echo ERROR: staging copy failed with robocopy code %ROBOCOPY_RC%.
  pause
  exit /b 1
)

cd /d "%STAGE_ROOT%\installer"
"%ISCC%" "JSDMStudio_Setup.iss"
set "BUILD_RC=%ERRORLEVEL%"

if %BUILD_RC% NEQ 0 (
  echo.
  echo Build failed.
  echo Staged source was left at:
  echo "%STAGE_ROOT%"
  pause
  exit /b 1
)

if not exist "%ORIG_ROOT%\installer\output" mkdir "%ORIG_ROOT%\installer\output"
copy /Y "%STAGE_ROOT%\installer\output\JSDMStudio_Setup.exe" "%ORIG_ROOT%\installer\output\" >nul

echo.
echo Build completed.
echo Installer:
echo "%ORIG_ROOT%\installer\output\JSDMStudio_Setup.exe"
echo.
pause
