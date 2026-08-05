@echo off
cd /d "%~dp0"
if exist startup_log.txt (
  start notepad startup_log.txt
) else (
  echo startup_log.txt does not exist yet.
  pause
)
