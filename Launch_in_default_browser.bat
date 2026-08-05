@echo off
chcp 65001 > nul
cd /d "%~dp0"

echo ==========================================
echo JSDMWorkbench - Browser Mode Dynamic Port
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
:found
if "%RSCRIPT%"=="" (
  echo Rscript.exe not found. Please install R for Windows.
  pause
  exit /b 1
)

echo Found Rscript:
echo "%RSCRIPT%"
echo.

echo Starting JSDMWorkbench on a free local port...
"%RSCRIPT%" -e "port <- httpuv::randomPort(); cat('Selected port:', port, '\n'); shiny::runApp('.', host='127.0.0.1', port=port, launch.browser=TRUE)"

pause
