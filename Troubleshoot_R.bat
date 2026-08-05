@echo off
chcp 65001 > nul
cd /d "%~dp0"
echo JSDM Studio Troubleshooter
echo =========================
echo.
echo Checking Rscript...
where Rscript.exe
echo.
echo Common R folders:
dir "C:\Program Files\R" /b
echo.
echo Current folder:
echo %CD%
echo.
pause
