@echo off
chcp 65001 > nul
setlocal enabledelayedexpansion
cd /d "%~dp0"

title JSDM Studio Launcher
echo ==========================================
echo JSDM Studio - Hmsc-powered Workbench
echo ==========================================
echo.

set "RSCRIPT="

REM 1. Try Rscript from PATH
where Rscript.exe >nul 2>nul
if %ERRORLEVEL%==0 (
  for /f "delims=" %%i in ('where Rscript.exe') do (
    set "RSCRIPT=%%i"
    goto :found_r
  )
)

REM 2. Try common system-wide R installations
for /d %%D in ("C:\Program Files\R\R-*") do (
  if exist "%%D\bin\Rscript.exe" (
    set "RSCRIPT=%%D\bin\Rscript.exe"
  )
)

REM 3. Try user-level R installations
for /d %%D in ("%LOCALAPPDATA%\Programs\R\R-*") do (
  if exist "%%D\bin\Rscript.exe" (
    set "RSCRIPT=%%D\bin\Rscript.exe"
  )
)

:found_r
if "%RSCRIPT%"=="" (
  echo ERROR: Rscript.exe was not found.
  echo.
  echo JSDM Studio needs R for Windows to run Hmsc.
  echo Please install R first, then launch JSDM Studio again.
  echo.
  echo Opening the official R download page...
  start "" "https://cran.r-project.org/bin/windows/base/"
  echo.
  pause
  exit /b 1
)

echo Found Rscript:
echo "%RSCRIPT%"
echo.

REM Check/install required R packages before launch.
echo Checking required R packages...
"%RSCRIPT%" check_packages.R
if %ERRORLEVEL% NEQ 0 (
  echo.
  echo Package check or installation failed.
  echo Please check your internet connection or run install_packages.bat manually.
  echo.
  pause
  exit /b 1
)

echo.
echo Starting JSDM Studio...
echo If the browser does not open automatically, look for the Shiny URL in this window.
echo.

"%RSCRIPT%" -e "shiny::runApp('.', launch.browser=TRUE)"

echo.
echo JSDM Studio has stopped.
pause
