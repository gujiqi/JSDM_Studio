#!/usr/bin/env bash
set -euo pipefail

cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
if [[ "$(uname -s)" == Linux && -z "${JSDM_RSCRIPT:-}" ]]; then
  support="${XDG_DATA_HOME:-$HOME/.local/share}/JSDMStudio"
  marker="$support/current-runtime-linux"
  if [[ -f "$marker" ]]; then
    IFS= read -r release < "$marker"
    if [[ "$release" != "$support/runtimes/"* || ! -d "$release/r-library" ||
          ! -x "$release/python/bin/python" ]]; then
      echo "The validated Linux runtime is missing. Run Install_JSDM_Studio_Linux.sh again." >&2
      exit 1
    fi
    unset R_HOME R_LIBS PYTHONHOME PYTHONPATH
    export JSDM_RSCRIPT=/usr/bin/Rscript
    export R_LIBS_USER="$release/r-library" R_LIBS_SITE=NULL
    export RETICULATE_PYTHON="$release/python/bin/python"
    export JSDM_HMSCHPC_PYTHON="$RETICULATE_PYTHON" PYTHONNOUSERSITE=1
  fi
fi
rscript="${JSDM_RSCRIPT:-$(command -v Rscript || true)}"
if [[ -z "$rscript" || ! -x "$rscript" ]]; then
  printf '%s\n' 'Rscript is required. Install R or set JSDM_RSCRIPT to its executable path.' >&2
  exit 1
fi
exec "$rscript" --vanilla run_shiny_backend.R --browser "$@"
