audit_hmsc_convergence <- function(model_path, output_dir) {
  model_path <- normalizePath(model_path, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (!requireNamespace("Hmsc", quietly = TRUE)) stop("Hmsc is required", call. = FALSE)
  if (!requireNamespace("coda", quietly = TRUE)) stop("coda is required", call. = FALSE)

  model <- readRDS(model_path)
  posterior <- Hmsc::convertToCodaObject(model)

  audit_block <- function(block, block_name) {
  if (is.null(block) || !inherits(block, "mcmc.list")) {
    return(data.frame())
  }
  parameter_count <- NCOL(block[[1]])
  if (!length(parameter_count) || parameter_count < 1L) return(data.frame())

  psrf <- tryCatch(
    coda::gelman.diag(block, multivariate = FALSE, autoburnin = FALSE)$psrf,
    error = function(e) NULL
  )
  ess <- tryCatch(coda::effectiveSize(block), error = function(e) NULL)
  parameters <- colnames(block[[1]])
  if (is.null(parameters)) parameters <- paste0(block_name, "_", seq_len(parameter_count))

  data.frame(
    block = block_name,
    parameter = parameters,
    psrf_point = if (is.null(psrf)) NA_real_ else as.numeric(psrf[, 1]),
    psrf_upper = if (is.null(psrf)) NA_real_ else as.numeric(psrf[, 2]),
    effective_size = if (is.null(ess)) NA_real_ else as.numeric(ess),
    stringsAsFactors = FALSE
  )
  }

  blocks <- intersect(c("Beta", "Gamma", "Rho", "Omega", "Alpha"), names(posterior))
  details <- do.call(rbind, lapply(blocks, function(x) audit_block(posterior[[x]], x)))
  if (is.null(details)) details <- data.frame()

  write.csv(details, file.path(output_dir, "MCMC_parameter_level_diagnostics.csv"), row.names = FALSE)

  if (nrow(details)) {
    summary <- do.call(rbind, lapply(split(details, details$block), function(x) {
      data.frame(
        block = x$block[[1]],
        parameter_count = nrow(x),
        psrf_median = median(x$psrf_point, na.rm = TRUE),
        psrf_max = max(x$psrf_point, na.rm = TRUE),
        psrf_le_1_01 = mean(x$psrf_point <= 1.01, na.rm = TRUE),
        psrf_le_1_05 = mean(x$psrf_point <= 1.05, na.rm = TRUE),
        ess_median = median(x$effective_size, na.rm = TRUE),
        ess_min = min(x$effective_size, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    summary <- data.frame()
  }

  write.csv(summary, file.path(output_dir, "MCMC_block_diagnostics_summary.csv"), row.names = FALSE)
  invisible(list(details = details, summary = summary))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L) {
    stop("Usage: Rscript audit_hmsc_convergence.R MODEL_RDS OUTPUT_DIR", call. = FALSE)
  }
  audit_hmsc_convergence(args[[1]], args[[2]])
  cat("MODEL=", normalizePath(args[[1]], winslash = "/", mustWork = TRUE), "\n", sep = "")
  cat("DETAILS=", file.path(normalizePath(args[[2]], winslash = "/", mustWork = FALSE), "MCMC_parameter_level_diagnostics.csv"), "\n", sep = "")
  cat("SUMMARY=", file.path(normalizePath(args[[2]], winslash = "/", mustWork = FALSE), "MCMC_block_diagnostics_summary.csv"), "\n", sep = "")
}

if (!interactive() && identical(sys.nframe(), 0L)) main()
