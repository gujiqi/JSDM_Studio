@echo off
set "RSCRIPT="

REM 1. Try PATH first
where Rscript.exe >nul 2>nul
if %ERRORLEVEL%==0 (
  for /f "delims=" %%i in ('where Rscript.exe') do (
    set "RSCRIPT=%%i"
    goto :found
  )
)

REM 2. Try common R installation folders
for /d %%D in ("C:\Program Files\R\R-*") do (
  if exist "%%D\bin\Rscript.exe" (
    set "RSCRIPT=%%D\bin\Rscript.exe"
  )
)

for /d %%D in ("%LOCALAPPDATA%\Programs\R\R-*") do (
  if exist "%%D\bin\Rscript.exe" (
    set "RSCRIPT=%%D\bin\Rscript.exe"
  )
)

:found
if "%RSCRIPT%"=="" (
  echo ERROR: Cannot find Rscript.exe.
  echo.
  echo Please check whether R is installed.
  echo Common path: C:\Program Files\R\R-4.3.3\bin\Rscript.exe
  echo.
  pause
  exit /b 1
)

echo Found Rscript:
echo "%RSCRIPT%"
echo.

cd /d "%~dp0"
echo Running JSDM Studio command-line workflow...
"%RSCRIPT%" main.R
pause
