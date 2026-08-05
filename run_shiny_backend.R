# run_shiny_backend.R
# Used by the WebView2 Windows launcher.
# It starts the Shiny app without opening a browser.
# It also checks and installs missing R packages on first run.

args <- commandArgs(trailingOnly = TRUE)
port <- 3838
if (length(args) >= 1 && grepl("^[0-9]+$", args[[1]])) {
  port <- as.integer(args[[1]])
}

log_file <- file.path(getwd(), "startup_log.txt")
log_msg <- function(...) {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", paste(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
}

options(shiny.host = "127.0.0.1")
options(shiny.port = port)
options(shiny.launch.browser = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

log_msg("Starting JSDM Studio Shiny backend")
log_msg("Working directory: ", getwd())
log_msg("Port: ", port)
log_msg("R version: ", R.version.string)

required <- c(
  "shiny", "httpuv", "bslib", "DT", "yaml", "htmltools", "jsonlite", "ggplot2", "zip",
  "Hmsc", "coda", "ape", "corrplot"
)

check_missing <- function() {
  required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
}

missing <- check_missing()

if (length(missing) > 0) {
  log_msg("Missing packages before auto-install: ", paste(missing, collapse = ", "))
  log_msg("Trying automatic installation via install_packages.R. This may take several minutes on first run.")

  if (file.exists("install_packages.R")) {
    tryCatch(
      {
        source("install_packages.R", local = TRUE)
      },
      error = function(e) {
        log_msg("Automatic package installation failed: ", conditionMessage(e))
      }
    )
  } else {
    log_msg("install_packages.R not found.")
  }

  missing <- check_missing()
}

if (length(missing) > 0) {
  log_msg("Packages still missing after installation attempt: ", paste(missing, collapse = ", "))
  stop("Missing required R packages after automatic installation attempt: ",
       paste(missing, collapse = ", "),
       ". Please run install_packages.bat manually and check internet connection.",
       call. = FALSE)
}

log_msg("All required packages are available.")
log_msg("Calling shiny::runApp")
shiny::runApp(".", host = "127.0.0.1", port = port, launch.browser = FALSE)
