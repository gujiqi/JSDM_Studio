library_dir <- Sys.getenv("R_LIBS_USER")
stopifnot(identical(Sys.info()[["sysname"]], "Linux"), dir.exists(library_dir))
.libPaths(library_dir, include.site = FALSE)
options(
  repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/2026-09-03"),
  HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(R.version$platform, R.version$arch)),
  timeout = 300
)
packages <- c(
  "shiny", "httpuv", "bslib", "DT", "yaml", "htmltools", "jsonlite", "ggplot2", "zip", "digest",
  "Hmsc", "coda", "ape", "corrplot", "colorspace", "writexl", "vioplot",
  "jSDM", "gjam", "spOccupancy", "sjSDM", "boral", "rjags", "R2jags",
  "reticulate", "sf", "terra", "ragg", "svglite", "rmarkdown", "patchwork"
)
available <- available.packages()
missing <- setdiff(packages, rownames(available))
if (length(missing)) stop("Packages not available for this R version: ", paste(missing, collapse = ", "))
if (available["sjSDM", "Version"] != "1.0.7") {
  stop("The repository does not offer the expected Python-backed sjSDM 1.0.7.", call. = FALSE)
}
install.packages(packages, lib = library_dir, available = available, dependencies = NA)
source("installer/platforms/verify_posix_runtime.R")
