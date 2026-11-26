required <- c(
  "shiny", "httpuv", "bslib", "DT", "yaml", "htmltools", "jsonlite", "ggplot2", "zip", "digest",
  "Hmsc", "coda", "ape", "corrplot", "colorspace", "writexl", "vioplot",
  "jSDM", "gjam", "spOccupancy", "sjSDM", "boral", "rjags", "R2jags", "reticulate",
  "sf", "terra", "ragg", "svglite", "rmarkdown", "patchwork"
)
failed <- vapply(required, function(p) {
  tryCatch({
    suppressPackageStartupMessages(loadNamespace(p))
    ""
  }, error = conditionMessage)
}, character(1))
if (any(nzchar(failed))) {
  stop(paste(names(failed)[nzchar(failed)], failed[nzchar(failed)], sep = ": ", collapse = "\n"),
       call. = FALSE)
}
cat("JAGS:", as.character(rjags::jags.version()), "\n")
if (!isTRUE(sjSDM::is_torch_available())) stop("sjSDM could not load PyTorch.", call. = FALSE)
reticulate::py_run_string("import torch, torch_optimizer, pyro, madgrad")
cat("R_CORE_AND_ENGINE_PACKAGES_OK\n")
