#!/bin/bash
set -euo pipefail

if [[ "$(/usr/bin/uname -s)" != Darwin ]]; then
  echo "This installation kit is for macOS." >&2
  exit 1
fi
if (( EUID == 0 )); then
  echo "Run this script as your normal user, not with sudo. It requests elevation only when needed." >&2
  exit 1
fi

kit="$(cd "$(dirname "$0")" && pwd -P)"
payload="$kit/JSDMStudio"
deps="$kit/dependencies"
source "$payload/installer/platforms/macos_runtime.sh"
arch="$(macos_native_arch)"
if [[ "$(/usr/bin/uname -m)" != "$arch" ]]; then
  exec /usr/bin/arch "-$arch" /bin/bash "$0" "$@"
fi
if [[ ! -f "$deps/installers/R-4.5.3-$arch.pkg" || ! -d "$deps/r-packages" ||
      ! -d "$deps/python-wheels" || ! -f "$deps/SHA256SUMS.txt" ]]; then
  echo "Use the complete macOS kit for $arch." >&2
  exit 1
fi
mac_version="$(/usr/bin/sw_vers -productVersion)"
if (( ${mac_version%%.*} < 14 )); then
  echo "This offline kit requires macOS 14 or later." >&2
  exit 1
fi

support="$HOME/Library/Application Support/JSDMStudio"
mkdir -p "$support/runtimes"
log="$support/install-macos-$(date +%Y%m%d-%H%M%S).log"
exec > >(tee "$log") 2>&1
trap 'echo "Installation stopped. See $log for details." >&2' ERR

echo "Installing JSDM Studio 0.4.8 for macOS $arch"
echo "Verifying the offline dependency payload"
(cd "$deps" && /usr/bin/shasum -a 256 -c SHA256SUMS.txt)
unset R_HOME R_LIBS PYTHONHOME PYTHONPATH
export R_LIBS_SITE=NULL PYTHONNOUSERSITE=1 PIP_CONFIG_FILE=/dev/null
install_pkg() {
  echo "Installing ${1##*/}. macOS may request an administrator password."
  /usr/bin/sudo /usr/sbin/installer -pkg "$1" -target /
}

if ! rscript="$(macos_rscript "$arch")"; then
  install_pkg "$deps/installers/R-4.5.3-$arch.pkg"
  rscript="$(macos_rscript "$arch")"
fi
if [[ ! -e /usr/local/lib/libjags.4.dylib ]]; then
  install_pkg "$deps/installers/JAGS-4.3.2.pkg"
fi
python310=/Library/Frameworks/Python.framework/Versions/3.10/bin/python3.10
if [[ ! -x "$python310" ]]; then
  install_pkg "$deps/installers/python-3.10.11-macos11.pkg"
fi
"$python310" -I -c 'import platform, ssl, sys, venv; assert platform.machine() == sys.argv[1] and sys.version_info[:2] == (3, 10), "Existing Python 3.10 is not native to this Mac"' "$arch"

# Keep the active runtime usable until the new application and dependencies pass.
release="$(mktemp -d "$support/runtimes/0.4.8-$arch.XXXXXXXX")"
export R_LIBS_USER="$release/r-library"
export RETICULATE_PYTHON="$release/python/bin/python"
export JSDM_HMSCHPC_PYTHON="$RETICULATE_PYTHON"
mkdir -p "$R_LIBS_USER"
"$python310" -I -m venv "$release/python"
echo "Installing the bundled Python packages"
pip_offline() { "$RETICULATE_PYTHON" -I -m pip --isolated --disable-pip-version-check "$@"; }
pip_offline install --no-index --find-links "$deps/python-wheels" --upgrade setuptools wheel
pip_offline install --no-index --find-links "$deps/python-wheels" \
  -r "$deps/requirements-hmschpc-cpu-macos.txt" -r "$deps/requirements-sjsdm-cpu-macos.txt"
pip_offline check
"$RETICULATE_PYTHON" -I -c 'import sys; sys.path.insert(0, sys.argv[1]); import tensorflow, tensorflow_probability, pyhmsc; print("HMSC_HPC_PYTHON_OK")' \
  "$payload/external_packages/hmsc-hpc-main"
"$RETICULATE_PYTHON" -I -c 'import torch, torch_optimizer, pyro, madgrad; print("SJSDM_PYTHON_OK")'

echo "Installing and checking the bundled R packages"
"$rscript" --vanilla "$payload/installer/platforms/install_macos_r_packages.R" "$deps/r-packages"

app_dir="$release/app"
shared="$support/app"
launch_dir="$HOME/Applications"
mkdir -p "$app_dir" "$shared" "$launch_dir"
/usr/bin/rsync -a --exclude=/config.yml --exclude=/output/ --exclude=/input/ --exclude=/diagnostics/ "$payload/" "$app_dir/"
if [[ ! -e "$shared/config.yml" ]]; then
  cp -n "$payload/config.yml" "$shared/config.yml"
fi
ln -s "$shared/config.yml" "$app_dir/config.yml"
for dir in input output diagnostics; do
  mkdir -p "$shared/$dir"
  ln -s "$shared/$dir" "$app_dir/$dir"
done
[[ -f "$app_dir/app.R" && -f "$app_dir/Launch_in_browser.sh" ]]
entry="$(mktemp "$launch_dir/.JSDMStudio-launch.XXXXXXXX")"
install -m 755 "$kit/Launch_JSDM_Studio.command" "$entry"
mv -f "$entry" "$launch_dir/JSDM Studio.command"
marker="$(mktemp "$support/.current-runtime-$arch.XXXXXXXX")"
printf '%s\n' "$release" > "$marker"
mv -f "$marker" "$support/current-runtime-$arch"
echo "Installation complete. Open JSDM Studio.command in your home Applications folder."
echo "The first launch opens your local browser."
