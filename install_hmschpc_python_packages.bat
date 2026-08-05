@echo off
setlocal
echo ==========================================
echo JSDM Studio - Hmsc-HPC Python packages
echo ==========================================
where python >nul 2>nul
if errorlevel 1 (
  echo Python was not found on PATH.
  echo Install Python 3.10+ and rerun this script.
  pause
  exit /b 1
)
python --version
python -m pip install --upgrade pip setuptools wheel
python -m pip install -e external_packages\hmsc-hpc-main
python -m pip install arviz biopython
python -c "import tensorflow, tensorflow_probability, pyhmsc, numpy, pandas; print('Hmsc-HPC Python stack OK')"
pause
