# run_shiny_backend.R
# Starts the local Shiny backend for the desktop and browser launchers.

private_r <- FALSE
if (.Platform$OS.type == "windows") {
  Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
  Sys.setenv(PYTHONIOENCODING = "utf-8")
  private_r_home <- file.path(getwd(), "runtimes", "R-4.5.3")
  if (dir.exists(private_r_home)) {
    private_r <- TRUE
    if (tolower(normalizePath(R.home(), winslash = "/")) !=
        tolower(normalizePath(private_r_home, winslash = "/"))) {
      stop("Use the bundled Rscript via the Studio launcher.", call. = FALSE)
    }
    Sys.setenv(JAGS_HOME = file.path(getwd(), "runtimes", "JAGS-4.3.2"))
    .libPaths(R.home("library"), include.site = FALSE)
  }
  python_candidates <- c(
    file.path(getwd(), "runtimes", "python310", "python.exe"),
    file.path(getwd(), "runtimes", "hmschpc-python310", "Scripts", "python.exe")
  )
  shared_python <- python_candidates[file.exists(python_candidates)]
  if (length(shared_python)) {
    shared_python <- shared_python[[1L]]
    shared_python <- normalizePath(shared_python, winslash = "/", mustWork = TRUE)
    Sys.unsetenv(c("PYTHONHOME", "PYTHONPATH"))
    Sys.setenv(RETICULATE_PYTHON = shared_python, JSDM_HMSCHPC_PYTHON = shared_python,
               PYTHONNOUSERSITE = "1")
  }
}

args <- commandArgs(trailingOnly = TRUE)
browser_mode <- "--browser" %in% args
port_arg <- args[grepl("^[0-9]+$", args)]

log_dir <- file.path(getwd(), "diagnostics")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, "shiny_backend.log")
log_msg <- function(...) {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", paste(..., collapse = ""))
  cat(txt, "\n")
  tryCatch(cat(txt, "\n", file = log_file, append = TRUE),
           warning = function(w) message("Could not write backend log: ", conditionMessage(w)),
           error = function(e) message("Could not write backend log: ", conditionMessage(e)))
}

options(repos = c(CRAN = "https://cloud.r-project.org"))

log_msg("Starting JSDM Studio Shiny backend")
log_msg("Working directory: ", getwd())
log_msg("R version: ", R.version.string)

required <- c("shiny", "httpuv", "bslib", "DT", "yaml", "htmltools", "jsonlite", "ggplot2", "zip")
if (!browser_mode) required <- c(required, "Hmsc", "coda", "ape", "corrplot")

check_missing <- function() {
  required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
}

missing <- check_missing()

if (length(missing) > 0 && !browser_mode && !private_r) {
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
  log_msg("Missing R packages: ", paste(missing, collapse = ", "))
  stop("Missing R packages: ", paste(missing, collapse = ", "),
       if (private_r) ". Re-run the full offline installer." else
         ". Run the platform installer or Rscript install_packages.R.", call. = FALSE)
}

port <- if (length(port_arg)) as.integer(port_arg[[1]]) else if (browser_mode) httpuv::randomPort() else 3838L
if (is.na(port) || port < 1L || port > 65535L) stop("Invalid localhost port.", call. = FALSE)
log_msg("Listening at http://127.0.0.1:", port)
shiny::runApp(".", host = "127.0.0.1", port = port,
              launch.browser = browser_mode && !("--no-open" %in% args))
