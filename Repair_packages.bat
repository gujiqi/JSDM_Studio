@echo off
chcp 65001 > nul
cd /d "%~dp0"

title JSDM Studio Repair Packages
echo ==========================================
echo JSDM Studio - Repair / Install R Packages
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
echo Installing/checking packages. This may take several minutes...
"%RSCRIPT%" install_packages.R

echo.
echo Done. You can launch JSDM Studio again.
pause
