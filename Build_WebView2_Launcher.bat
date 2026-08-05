@echo off
chcp 65001 > nul
cd /d "%~dp0"

echo ==========================================
echo Building JSDM Studio WebView2 Launcher
echo ==========================================
echo.

where dotnet >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
  echo ERROR: dotnet was not found.
  echo.
  echo Please install .NET SDK first.
  echo.
  pause
  exit /b 1
)

cd webview2_launcher

echo Restoring and publishing launcher...
dotnet publish -c Release -r win-x64 --self-contained false -o ..\launcher_build

if %ERRORLEVEL% NEQ 0 (
  echo.
  echo Build failed.
  pause
  exit /b 1
)

cd ..

echo.
echo Copying launcher files and dependencies...
copy /Y launcher_build\*.* .\

if not exist JSDMStudioLauncher.exe (
  echo.
  echo ERROR: JSDMStudioLauncher.exe was not copied.
  pause
  exit /b 1
)

echo.
echo Build complete:
echo %CD%\JSDMStudioLauncher.exe
echo.
echo Important: WebView2 launcher also needs the DLL/deps/runtimeconfig files copied from launcher_build.
echo This script has copied them to the JSDMStudio folder.
echo.
pause
