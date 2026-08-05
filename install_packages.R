# install_packages.R
# Install or repair R packages used by JSDM Studio.

options(repos = c(CRAN = "https://cloud.r-project.org"))

core_packages <- c(
  "shiny", "httpuv", "bslib", "DT", "yaml", "htmltools", "jsonlite",
  "ggplot2", "zip", "Hmsc", "coda", "ape", "corrplot",
  "colorspace", "writexl", "vioplot"
)

optional_engine_packages <- c(
  "jSDM", "gjam", "spOccupancy", "sjSDM", "boral", "R2jags",
  "abind", "lme4", "foreach", "doParallel", "spAbundance",
  "reticulate", "mvtnorm", "Metrics", "mgcv", "checkmate", "qgam",
  "viridis", "scales", "beeswarm", "corpcor", "fishMod", "reshape2",
  "MASS", "RANN"
)

install_if_missing <- function(packages, required = TRUE) {
  failed <- character()
  for (pkg in unique(packages)) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      message("Package already installed: ", pkg)
      next
    }
    message("Installing package: ", pkg)
    ok <- tryCatch({
      install.packages(pkg, dependencies = TRUE)
      requireNamespace(pkg, quietly = TRUE)
    }, error = function(e) {
      message("Installation failed for ", pkg, ": ", conditionMessage(e))
      FALSE
    })
    if (!isTRUE(ok)) failed <- c(failed, pkg)
  }
  if (length(failed) > 0 && required) {
    stop("Required package installation failed: ", paste(failed, collapse = ", "), call. = FALSE)
  }
  failed
}

message("Installing/checking required packages for JSDM Studio...")
core_failed <- install_if_missing(core_packages, required = TRUE)

message("Installing/checking optional engine packages...")
optional_failed <- install_if_missing(optional_engine_packages, required = FALSE)

if (length(optional_failed) > 0) {
  message("Optional packages not installed: ", paste(optional_failed, collapse = ", "))
  message("The app can still start. Workflows that need those packages will report scaffold_only or fit_failed diagnostics.")
}

message("External software notes:")
message("- boral production fitting requires JAGS installed outside R.")
message("- sjSDM production fitting may require Python, reticulate, PyTorch, and optionally CUDA.")
message("- Hmsc-HPC workflow requires Python 3.10+ with numpy, pandas, patsy, PyYAML, h5py, scipy, tensorflow, tensorflow-probability, ujson; arviz and biopython are recommended for diagnostics/Newick trees.")
message("- Hmsc-HPC is bundled under external_packages/hmsc-hpc-main; the GUI runs it on CPU through PYTHONPATH.")
message("- Some packages may require Rtools on Windows.")
message("Package check complete.")
