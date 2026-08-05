# check_packages.R
# Called by launchers and troubleshooting scripts.

options(repos = c(CRAN = "https://cloud.r-project.org"))

required <- c(
  "shiny", "httpuv", "bslib", "DT", "yaml", "htmltools", "jsonlite",
  "ggplot2", "zip", "Hmsc", "coda", "ape", "corrplot"
)

optional <- c("jSDM", "gjam", "spOccupancy", "sjSDM", "boral", "R2jags")

missing_required <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_required) > 0) {
  message("Missing required packages: ", paste(missing_required, collapse = ", "))
  message("Running install_packages.R ...")
  source("install_packages.R", local = TRUE)
  missing_required <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
}

if (length(missing_required) > 0) {
  stop("Required packages are still missing: ", paste(missing_required, collapse = ", "), call. = FALSE)
}

missing_optional <- optional[!vapply(optional, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_optional) > 0) {
  message("Optional engine packages not available: ", paste(missing_optional, collapse = ", "))
  message("Those workflows will create diagnostics/scaffolds until their adapters and dependencies are installed.")
}

message("Required R packages are available.")
