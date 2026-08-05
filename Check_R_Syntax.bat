@echo off
chcp 65001 > nul
cd /d "%~dp0"
echo ==========================================
echo JSDMWorkbench - R Syntax Check
echo ==========================================
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
:found
if "%RSCRIPT%"=="" (
  echo Rscript.exe not found. Please install R for Windows.
  pause
  exit /b 1
)
"%RSCRIPT%" -e "parse('app.R'); cat('app.R syntax OK\n')"
if %ERRORLEVEL% NEQ 0 (
  echo app.R syntax check failed.
  pause
  exit /b 1
)
echo Syntax check passed.
pause
