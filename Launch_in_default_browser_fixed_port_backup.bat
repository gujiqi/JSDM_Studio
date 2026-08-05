@echo off
chcp 65001 > nul
cd /d "%~dp0"

title JSDM Studio Browser Mode
echo ==========================================
echo JSDM Studio - Browser Mode
echo ==========================================
echo.

set "RSCRIPT="
where Rscript.exe >nul 2>nul
if %ERRORLEVEL%==0 (
  for /f "delims=" %%i in ('where Rscript.exe') do (
    set "RSCRIPT=%%i"
    goto :found
  )
)
for /d %%D in ("C:\Program Files\R\R-*") do (
  if exist "%%D\bin\Rscript.exe" set "RSCRIPT=%%D\bin\Rscript.exe"
)
for /d %%D in ("%LOCALAPPDATA%\Programs\R\R-*") do (
  if exist "%%D\bin\Rscript.exe" set "RSCRIPT=%%D\bin\Rscript.exe"
)

:found
if "%RSCRIPT%"=="" (
  echo Rscript.exe not found. Please install R for Windows first.
  start "" "https://cran.r-project.org/bin/windows/base/"
  pause
  exit /b 1
)

echo Found Rscript:
echo "%RSCRIPT%"
echo.

echo Checking packages...
"%RSCRIPT%" check_packages.R
if %ERRORLEVEL% NEQ 0 (
  echo Package check failed. Try Repair_packages.bat.
  pause
  exit /b 1
)

echo Starting JSDM Studio in default browser...
"%RSCRIPT%" -e "shiny::runApp('.', host='127.0.0.1', port=3838, launch.browser=TRUE)"

pause
