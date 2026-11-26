#!/bin/bash

macos_native_arch() {
  if [[ "$(/usr/sbin/sysctl -n hw.optional.arm64 2>/dev/null || true)" == 1 ]]; then
    printf '%s\n' arm64
  else
    /usr/bin/uname -m
  fi
}

macos_rscript() {
  local target_arch="$1" candidate
  local r_arch="$target_arch"
  if [[ "$target_arch" == arm64 ]]; then r_arch=aarch64; fi
  for candidate in \
    "${JSDM_RSCRIPT:-}" \
    "/Library/Frameworks/R.framework/Versions/4.5-$target_arch/Resources/bin/Rscript" \
    "/Library/Frameworks/R.framework/Resources/bin/Rscript" \
    "$(command -v Rscript || true)"
  do
    if [[ -x "$candidate" ]] && "$candidate" --vanilla -e \
      'a <- commandArgs(TRUE)[1]; quit(status = as.integer(!(getRversion() >= "4.5.0" && getRversion() < "4.6.0" && R.version$arch %in% c(a, if (a == "aarch64") "arm64" else a))))' \
      "$r_arch" >/dev/null 2>&1
    then
      printf '%s\n' "$candidate"
      return
    fi
  done
  return 1
}
