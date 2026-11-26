#!/usr/bin/env bash
set -euo pipefail

if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
  echo "The guided installer requires Ubuntu 24.04 LTS x86_64." >&2
  exit 1
fi
. /etc/os-release
if [[ "${ID:-}:${VERSION_ID:-}" != ubuntu:24.04 ]]; then
  echo "The guided installer targets Ubuntu 24.04 LTS. See INSTALL_LINUX.md for manual setup." >&2
  exit 1
fi
if ! grep -qw avx /proc/cpuinfo; then
  echo "The TensorFlow CPU backend requires AVX. Use manual setup for R-only engines." >&2
  exit 1
fi
if [[ -f app.R && -d installer/platforms ]]; then
  cd -- "$(pwd -P)"
else
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
fi
if [[ ! -f app.R || ! -f installer/platforms/install_linux_r_packages.R ]]; then
  echo "Extract the whole Linux ZIP and run the installer in the JSDMStudio folder." >&2
  exit 1
fi
system_packages=(
  ca-certificates fonts-noto-core fonts-noto-cjk
  r-base r-base-dev python3.12-venv python3.12-dev jags build-essential gfortran cmake pkg-config
  libcurl4-openssl-dev libssl-dev libxml2-dev libuv1-dev libgsl-dev libgdal-dev libgeos-dev
  libproj-dev libudunits2-dev libsqlite3-dev libfontconfig1-dev libfreetype-dev
  libharfbuzz-dev libfribidi-dev libpng-dev libtiff-dev libjpeg-dev libicu-dev
  libblas-dev liblapack-dev
)
system_missing=()
for package in "${system_packages[@]}"; do
  if ! status="$(dpkg-query -W -f='${Status}' "$package" 2>/dev/null)" ||
       [[ "$status" != "install ok installed" ]]; then
    system_missing+=("$package")
  fi
done
if [[ "${1:-}" == --check ]]; then
  echo "Installer platform check passed: Ubuntu 24.04 LTS x86_64 with AVX."
  if (( ${#system_missing[@]} )); then
    printf 'Missing system packages: %s\n' "${system_missing[*]}"
    exit 2
  fi
  echo "System package check passed. R/Python validation runs during installation."
  exit 0
fi
if [[ $# -gt 1 || ( $# -eq 1 && "$1" != --install-system ) ]]; then
  echo "Usage: bash Install_JSDM_Studio_Linux.sh [--check | --install-system]" >&2
  exit 1
fi
if (( EUID == 0 )); then
  echo "Run as your normal user. Only the apt step uses sudo." >&2
  exit 1
fi
pypi_index="${JSDM_STUDIO_PYPI_INDEX:-https://pypi.org/simple}"
if [[ "$pypi_index" != https://* ]]; then
  echo "JSDM_STUDIO_PYPI_INDEX must be an HTTPS package index." >&2
  exit 1
fi
support="${XDG_DATA_HOME:-$HOME/.local/share}/JSDMStudio"
mkdir -p "$support/runtimes"
log="$support/install-linux-$(date +%Y%m%d-%H%M%S).log"
exec > >(tee "$log") 2>&1
trap 'echo "Installation stopped. See $log for details." >&2' ERR
if (( ${#system_missing[@]} )); then
  if [[ "${1:-}" != --install-system ]]; then
    echo "Missing system packages: ${system_missing[*]}"
    echo "Re-run with --install-system to authorize their installation using apt."
    exit 2
  fi
  sudo apt-get update
  sudo apt-get install "${system_missing[@]}"
fi

release="$(mktemp -d "$support/runtimes/0.4.8-ubuntu2404.XXXXXXXX")"
unset R_HOME R_LIBS PYTHONHOME PYTHONPATH
export R_LIBS_USER="$release/r-library" R_LIBS_SITE=NULL
export RETICULATE_PYTHON="$release/python/bin/python"
export JSDM_HMSCHPC_PYTHON="$RETICULATE_PYTHON"
export PYTHONNOUSERSITE=1 PIP_CONFIG_FILE=/dev/null
mkdir -p "$R_LIBS_USER"
/usr/bin/python3.12 -I -m venv "$release/python"
pip_cpu() { "$RETICULATE_PYTHON" -I -m pip --isolated --disable-pip-version-check "$@"; }
# Fetch the CPU wheel from PyTorch, but resolve its general dependencies on PyPI.
pip_cpu install --no-deps --index-url https://download.pytorch.org/whl/cpu 'torch==2.5.1+cpu'
pip_cpu install --index-url "$pypi_index" 'torch==2.5.1' -r installer/platforms/requirements-linux-cpu.txt
pip_cpu check
"$RETICULATE_PYTHON" -I -c 'import sys; sys.path.insert(0, sys.argv[1]); import tensorflow, tensorflow_probability, pyhmsc, torch, torch_optimizer, pyro, madgrad; print("PYTHON_CPU_BACKENDS_OK")' \
  "$PWD/external_packages/hmsc-hpc-main"
/usr/bin/Rscript --vanilla installer/platforms/install_linux_r_packages.R
marker="$(mktemp "$support/.current-runtime-linux.XXXXXXXX")"
printf '%s\n' "$release" > "$marker"
mv -f "$marker" "$support/current-runtime-linux"
echo "Installation complete. Run: bash Launch_in_browser.sh"
