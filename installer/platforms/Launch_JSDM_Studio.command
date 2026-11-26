#!/bin/bash
set -euo pipefail

if [[ "$(/usr/bin/uname -s)" != Darwin ]]; then
  echo "This launcher requires macOS." >&2
  exit 1
fi
support="$HOME/Library/Application Support/JSDMStudio"
if [[ "$(/usr/sbin/sysctl -n hw.optional.arm64 2>/dev/null || true)" == 1 ]]; then
  arch=arm64
else
  arch=x86_64
fi
if [[ "$(/usr/bin/uname -m)" != "$arch" ]]; then
  exec /usr/bin/arch "-$arch" /bin/bash "$0" "$@"
fi
marker="$support/current-runtime-$arch"
if [[ -f "$marker" ]]; then
  IFS= read -r release < "$marker"
  if [[ "$release" != "$support/runtimes/"* || ! -f "$release/app/app.R" ||
        ! -d "$release/r-library" || ! -x "$release/python/bin/python" ]]; then
    echo "The verified runtime is missing. Run the macOS installer again." >&2
    exit 1
  fi
  app_dir="$release/app"
  export R_LIBS_USER="$release/r-library"
  export RETICULATE_PYTHON="$release/python/bin/python"
else
  app_dir="$support/app"
  export R_LIBS_USER="$support/r-library-4.5-$arch"
  export RETICULATE_PYTHON="$support/python-3.10-$arch/bin/python"
fi
if [[ ! -f "$app_dir/app.R" || ! -d "$R_LIBS_USER" || ! -x "$RETICULATE_PYTHON" ]]; then
  echo "Run Install_JSDM_Studio.command in the extracted kit first." >&2
  exit 1
fi
source "$app_dir/installer/platforms/macos_runtime.sh"
unset R_HOME R_LIBS PYTHONHOME PYTHONPATH
export R_LIBS_SITE=NULL PYTHONNOUSERSITE=1
if ! JSDM_RSCRIPT="$(macos_rscript "$arch")"; then
  echo "Matching R 4.5 is unavailable. Run Install_JSDM_Studio.command again." >&2
  exit 1
fi
export JSDM_RSCRIPT
export JSDM_HMSCHPC_PYTHON="$RETICULATE_PYTHON"
exec /bin/bash "$app_dir/Launch_in_browser.sh" "$@"
