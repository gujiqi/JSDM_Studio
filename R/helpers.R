# R/helpers.R
# Shared helper functions for JSDM Studio output workflows.
# These utilities are retained for backward compatibility with older Hmsc scripts.

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
}


safe_read_csv <- function(path, label = "file") {
  if (is.null(path) || is.na(path) || path == "" || !file.exists(path)) {
    stop(sprintf("Required file not found for %s: %s", label, path), call. = FALSE)
  }
  read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = TRUE)
}


safe_try <- function(expr, label = "step", log_fun = function(x) message(x)) {
  tryCatch(expr, error = function(e) {
    log_fun(paste0(label, " failed: ", e$message))
    NULL
  })
}

read_uploaded_csv <- function(file_input, label = "file") {
  if (is.null(file_input)) stop(sprintf("Please upload %s before running this workflow.", label), call. = FALSE)
  read.csv(file_input$datapath, row.names = 1, check.names = FALSE, stringsAsFactors = TRUE)
}

make_output_dir <- function(prefix = "hmsc") {
  out <- file.path("output", paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "models"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "results"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "inputs"), recursive = TRUE, showWarnings = FALSE)
  return(out)
}

write_lines <- function(x, file, append = FALSE) {
  cat(paste(x, collapse = "\n"), file = file, sep = "\n", append = append)
}

capture_to_file <- function(expr, file) {
  con <- file(file, open = "wt", encoding = "UTF-8")
  sink(con)
  on.exit({
    try(sink(), silent = TRUE)
    try(close(con), silent = TRUE)
  }, add = TRUE)
  force(expr)
}

validate_hmsc_inputs <- function(Y, XData, TrData = NULL, studyDesign = NULL, coordinates = NULL) {
  messages <- character()
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  suppressWarnings(storage.mode(Y) <- "numeric")
  if (!is.numeric(Y)) stop("Y.csv must be a numeric response matrix. Use 0/1 for presence-absence, non-negative integers for counts, or numeric continuous values for normal-response workflows.", call. = FALSE)
  if (nrow(Y) != nrow(XData)) stop("Y.csv and XData.csv must have the same number of rows; each row must represent the same sampling unit in the same order.", call. = FALSE)
  if (!is.null(rownames(Y)) && !is.null(rownames(XData)) && !all(rownames(Y) == rownames(XData))) {
    stop("Y.csv and XData.csv row names do not match. Align sampling-unit names or remove inconsistent row names before fitting.", call. = FALSE)
  }
  if (anyNA(Y)) stop("Y.csv contains missing values (NA). Handle missing responses before fitting this Hmsc workflow.", call. = FALSE)
  if (anyNA(XData)) messages <- c(messages, "XData.csv contains NA values. Hmsc fitting may fail unless missing predictors are removed or imputed.")

  if (!is.null(TrData)) {
    if (!all(colnames(Y) %in% rownames(TrData))) {
      missing_sp <- setdiff(colnames(Y), rownames(TrData))
      stop(paste0("traits/TrData.csv is missing species: ", paste(missing_sp, collapse = ", ")), call. = FALSE)
    }
    messages <- c(messages, sprintf("TrData check: %s species rows and %s trait columns.", nrow(TrData), ncol(TrData)))
  }
  if (!is.null(studyDesign)) {
    if (nrow(studyDesign) != nrow(Y)) messages <- c(messages, "studyDesign.csv row count does not match Y.csv; random-effect design must align with sampling units.")
  }
  if (!is.null(coordinates)) {
    if (nrow(coordinates) != nrow(Y)) messages <- c(messages, "coordinates.csv row count does not match Y.csv; spatial random effects require one coordinate row per sampling unit.")
  }

  messages <- c(messages, sprintf("Data check passed: %s sampling units, %s species/responses and %s environmental predictors.", nrow(Y), ncol(Y), ncol(XData)))
  return(messages)
}

save_input_summaries <- function(Y, XData, TrData, studyDesign, coordinates, outdir) {
  tables <- file.path(outdir, "tables")
  write.csv(data.frame(
    item = c("samples", "species", "environmental_variables",
             "traits_variables", "studyDesign_columns", "coordinate_columns"),
    value = c(nrow(Y), ncol(Y), ncol(XData),
              ifelse(is.null(TrData), 0, ncol(TrData)),
              ifelse(is.null(studyDesign), 0, ncol(studyDesign)),
              ifelse(is.null(coordinates), 0, ncol(coordinates)))
  ), file.path(tables, "data_dimensions.csv"), row.names = FALSE)

  prev <- colMeans(Y, na.rm = TRUE)
  occ <- colSums(Y > 0, na.rm = TRUE)
  write.csv(data.frame(species = colnames(Y), prevalence_or_mean = prev, occurrences_positive = occ),
            file.path(tables, "species_summary.csv"), row.names = FALSE)

  env_summary <- data.frame(variable = colnames(XData),
                            class = sapply(XData, function(z) paste(class(z), collapse = "/")),
                            missing = sapply(XData, function(z) sum(is.na(z))))
  write.csv(env_summary, file.path(tables, "environment_summary.csv"), row.names = FALSE)
}

hmsc_factor_preserve_order <- function(x) {
  x <- as.character(x)
  x[trimws(x) == ""] <- NA
  factor(x, levels = unique(x[!is.na(x)]))
}

hmsc_normalize_spatial_method <- function(spatial_method = "Full") {
  method <- toupper(trimws(as.character(spatial_method %||% "Full")))
  if (method %in% c("NNGP", "NEAREST_NEIGHBOR_GP", "NEAREST_NEIGHBOUR_GP")) return("NNGP")
  if (method %in% c("GPP", "GAUSSIAN_PREDICTIVE_PROCESS")) return("GPP")
  "Full"
}

hmsc_random_units <- function(n, preferred = NULL, prefix = "unit") {
  preferred <- as.character(preferred %||% character())
  if (length(preferred) == n && !anyDuplicated(preferred) && all(nzchar(preferred))) return(preferred)
  sprintf("%s_%03d", prefix, seq_len(n))
}

hmsc_factor_study_design <- function(studyDesign, n, units, ensure_col = "sample") {
  if (is.null(studyDesign)) {
    studyDesign <- data.frame(sample = units, stringsAsFactors = FALSE)
    rownames(studyDesign) <- units
  }
  studyDesign <- as.data.frame(studyDesign, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(studyDesign) != n) stop("studyDesign.csv rows must match Y rows.", call. = FALSE)
  if (nzchar(ensure_col %||% "") && !(ensure_col %in% names(studyDesign))) studyDesign[[ensure_col]] <- units
  for (nm in names(studyDesign)) studyDesign[[nm]] <- hmsc_factor_preserve_order(studyDesign[[nm]])
  studyDesign
}

hmsc_select_coordinate_columns <- function(coordinates, longitude_column = "x", latitude_column = "y") {
  nm <- names(coordinates)
  find_col <- function(x) {
    hit <- nm[tolower(nm) == tolower(trimws(as.character(x %||% "")))]
    if (length(hit) > 0) hit[1] else NA_character_
  }
  requested <- c(find_col(longitude_column), find_col(latitude_column))
  if (all(!is.na(requested))) return(requested)
  for (pair in list(c("longitude", "latitude"), c("lon", "lat"), c("x", "y"), c("easting", "northing"))) {
    hit <- c(find_col(pair[1]), find_col(pair[2]))
    if (all(!is.na(hit))) return(hit)
  }
  numeric_cols <- nm[vapply(coordinates, function(z) all(is.na(z) | !is.na(suppressWarnings(as.numeric(z)))), logical(1))]
  if (length(numeric_cols) >= 2) return(numeric_cols[1:2])
  stop("Coordinate columns were not found in coordinates.csv.", call. = FALSE)
}

build_random_levels <- function(mode, Y, studyDesign = NULL, coordinates = NULL,
                                random_effect_column = "sample",
                                longitude_column = "x", latitude_column = "y",
                                spatial_method = "Full", nNeighbours = 10) {
  mode <- tolower(trimws(as.character(mode %||% "none")))
  if (mode %in% c("spatial_full", "full")) { mode <- "spatial"; spatial_method <- "Full" }
  if (mode %in% c("spatial_nngp", "nngp")) { mode <- "spatial"; spatial_method <- "NNGP" }
  if (mode %in% c("spatial_gpp", "gpp")) { mode <- "spatial"; spatial_method <- "GPP" }
  spatial_method <- hmsc_normalize_spatial_method(spatial_method)

  n <- nrow(Y)
  y_units <- hmsc_random_units(n, rownames(Y), prefix = "sample")
  if (mode == "none") return(list(studyDesign = NULL, ranLevels = NULL))

  if (mode == "sample") {
    if (is.null(studyDesign)) stop("Sample random effect requires studyDesign.csv.", call. = FALSE)
    if (!nzchar(trimws(random_effect_column %||% ""))) random_effect_column <- "sample"
    studyDesign <- hmsc_factor_study_design(studyDesign, n, y_units, random_effect_column)
    if (!(random_effect_column %in% names(studyDesign))) stop(sprintf("Grouping column '%s' was not found in studyDesign.csv.", random_effect_column), call. = FALSE)
    rL <- Hmsc::HmscRandomLevel(units = levels(studyDesign[[random_effect_column]]))
    return(list(studyDesign = studyDesign, ranLevels = setNames(list(rL), random_effect_column)))
  }

  if (mode == "spatial") {
    if (is.null(coordinates)) stop("Spatial random effect requires coordinates.csv.", call. = FALSE)
    if (nrow(coordinates) != n) stop("coordinates.csv rows must match Y rows.", call. = FALSE)
    spatial_col <- "spatial_unit"
    spatial_units <- hmsc_random_units(n, rownames(coordinates), prefix = "spatial")
    studyDesign <- hmsc_factor_study_design(studyDesign, n, y_units, spatial_col)
    studyDesign[[spatial_col]] <- factor(spatial_units, levels = spatial_units)
    coord_cols <- hmsc_select_coordinate_columns(coordinates, longitude_column, latitude_column)
    xy <- as.matrix(coordinates[, coord_cols, drop = FALSE])
    suppressWarnings(storage.mode(xy) <- "numeric")
    if (any(!is.finite(xy))) stop("Spatial coordinate columns must be numeric and finite.", call. = FALSE)
    rownames(xy) <- spatial_units
    if (spatial_method == "NNGP") {
      if (nNeighbours < 1) stop("NNGP requires nNeighbours >= 1.", call. = FALSE)
      rL <- Hmsc::HmscRandomLevel(sData = xy, sMethod = "NNGP", nNeighbours = nNeighbours)
    } else if (spatial_method == "GPP") {
      knots <- Hmsc::constructKnots(as.data.frame(xy), nKnots = max(2L, min(10L, floor(sqrt(nrow(xy))))))
      rL <- Hmsc::HmscRandomLevel(sData = xy, sMethod = "GPP", sKnot = knots)
    } else {
      rL <- Hmsc::HmscRandomLevel(sData = xy, sMethod = "Full")
    }
    return(list(studyDesign = studyDesign, ranLevels = setNames(list(rL), spatial_col)))
  }

  stop("random_effects mode must be none, sample, spatial_full, spatial_nngp or spatial_gpp.", call. = FALSE)
}

model_structure_output <- function(m, outdir) {
  f <- file.path(outdir, "results", "model_structure.txt")
  capture_to_file({
    cat("Hmsc model structure\n====================\n\n")
    print(m)
    cat("\n\nKey slots\n---------\n")
    cat("ns species:", m$ns, "\n")
    cat("nc covariates:", m$nc, "\n")
    cat("nt traits:", m$nt, "\n")
    cat("nr random levels:", m$nr, "\n")
    cat("\nSpecies names:\n")
    print(m$spNames)
    cat("\nCovariate names:\n")
    print(m$covNames)
    cat("\nTrait names:\n")
    print(m$trNames)
  }, f)
}

save_model_fit_tables <- function(mf, prefix, outdir) {
  if (is.null(mf)) return(invisible(NULL))
  tbl_dir <- file.path(outdir, "tables")
  for (nm in names(mf)) {
    x <- mf[[nm]]
    if (is.null(x)) next
    df <- data.frame(species = m$spNames %||% seq_along(x), value = as.numeric(x))
    names(df)[2] <- nm
    write.csv(df, file.path(tbl_dir, paste0(prefix, "_", nm, ".csv")), row.names = FALSE)
  }
  invisible(NULL)
}

plot_model_fit_compare <- function(MF, MFCV, outdir, model_name = "model") {
  if (is.null(MF) || is.null(MFCV)) return(invisible(NULL))
  pdf(file.path(outdir, "plots", "model_fit_explanatory_vs_predictive.pdf"), width = 8, height = 7)
  on.exit(dev.off(), add = TRUE)
  metrics <- intersect(names(MF), names(MFCV))
  for (met in metrics) {
    x <- MF[[met]]
    y <- MFCV[[met]]
    if (is.null(x) || is.null(y) || !is.numeric(x) || !is.numeric(y)) next
    lim <- range(c(x, y), na.rm = TRUE)
    if (!is.finite(lim[1])) next
    if (met %in% c("AUC")) lim <- c(0,1)
    plot(x, y, xlab = paste("explanatory", met), ylab = paste("predictive", met),
         main = paste0(model_name, ": ", met, "\nmean explanatory=", round(mean(x, na.rm=TRUE),3),
                       ", mean predictive=", round(mean(y, na.rm=TRUE),3)),
         xlim = lim, ylim = lim, pch = 19, col = "#2C7A7B")
    abline(0, 1, lty = 2)
    if (met == "AUC") { abline(v = 0.5, lty = 3); abline(h = 0.5, lty = 3) }
    text(x, y, labels = seq_along(x), pos = 3, cex = 0.6)
  }
}

mcmc_diagnostics_rich <- function(m, outdir, cfg, log_fun) {
  log_fun("  MCMC Beta / Gamma / Omega / Rho / Alpha...")
  f <- file.path(outdir, "results", "MCMC_convergence.txt")
  mpost <- Hmsc::convertToCodaObject(m, spNamesNumbers = c(TRUE, FALSE), covNamesNumbers = c(TRUE, FALSE))
  capture_to_file({
    cat("MCMC Convergence statistics\n===========================\n\n")
    cat("PSRF/Gelman  1 \n\n")

    write_diag <- function(obj, name) {
      if (is.null(obj)) return(NULL)
      cat("\n", name, "\n", paste(rep("-", nchar(name)), collapse = ""), "\n", sep = "")
      ess <- safe_try(coda::effectiveSize(obj), paste(name, "effectiveSize"), function(x) cat(x, "\n"))
      psrf <- safe_try(coda::gelman.diag(obj, multivariate = FALSE)$psrf, paste(name, "gelman"), function(x) cat(x, "\n"))
      cat("\nEffective size summary:\n"); print(summary(as.numeric(ess)))
      cat("\nGelman PSRF summary:\n"); print(summary(as.numeric(psrf[,1])))
    }

    write_diag(mpost$Beta, "Beta")
    write_diag(mpost$Gamma, "Gamma")
    if (!is.null(mpost$Rho)) write_diag(mpost$Rho, "Rho")

    if (!is.null(mpost$Omega)) {
      cat("\nOmega\n-----\n")
      for (k in seq_along(mpost$Omega)) {
        cat("\nRandom level", k, "\n")
        obj <- mpost$Omega[[k]]
        if (length(obj) > 0) {
          psrf <- safe_try(coda::gelman.diag(obj, multivariate = FALSE)$psrf, "Omega gelman", function(x) cat(x, "\n"))
          print(summary(as.numeric(psrf[,1])))
        }
      }
    }

    if (!is.null(mpost$Alpha)) {
      cat("\nAlpha\n-----\n")
      for (k in seq_along(mpost$Alpha)) {
        obj <- mpost$Alpha[[k]]
        if (!is.null(obj)) {
          psrf <- safe_try(coda::gelman.diag(obj, multivariate = FALSE)$psrf, "Alpha gelman", function(x) cat(x, "\n"))
          print(psrf)
        }
      }
    }
  }, f)

  # plots
  pdf(file.path(outdir, "plots", "MCMC_traceplots.pdf"), width = 11, height = 8)
  on.exit(dev.off(), add = TRUE)
  if (!is.null(mpost$Beta)) plot(mpost$Beta, main = "Trace plots: Beta")
  if (!is.null(mpost$Gamma)) plot(mpost$Gamma, main = "Trace plots: Gamma")
  if (!is.null(mpost$Rho)) plot(mpost$Rho, main = "Trace plots: Rho")
}

parameter_estimates_rich <- function(m, outdir, cfg, log_fun) {
  log_fun(" Beta / Gamma /  / Omega ...")
  tbl_dir <- file.path(outdir, "tables")
  plot_dir <- file.path(outdir, "plots")
  result_txt <- file.path(outdir, "results", "parameter_estimates.txt")

  capture_to_file({
    cat("Parameter estimates\n===================\n\n")
    cat(" Pr(x>0)Pr(x<0) Beta/Gamma/Omega \n")
  }, result_txt)

  # Beta
  safe_try({
    if (m$nc > 0) {
      postBeta <- Hmsc::getPostEstimate(m, parName = "Beta")
      beta_mean <- as.data.frame(t(postBeta$mean))
      beta_support <- as.data.frame(t(postBeta$support))
      beta_supportNeg <- as.data.frame(t(postBeta$supportNeg))
      beta_mean <- cbind(Species = m$spNames, beta_mean)
      beta_support <- cbind(Species = m$spNames, beta_support)
      beta_supportNeg <- cbind(Species = m$spNames, beta_supportNeg)
      write.csv(beta_mean, file.path(tbl_dir, "parameter_estimates_Beta_mean.csv"), row.names = FALSE)
      write.csv(beta_support, file.path(tbl_dir, "parameter_estimates_Beta_Pr_positive.csv"), row.names = FALSE)
      write.csv(beta_supportNeg, file.path(tbl_dir, "parameter_estimates_Beta_Pr_negative.csv"), row.names = FALSE)

      pdf(file.path(plot_dir, "Beta_plot.pdf"), width = 11, height = 8)
      Hmsc::plotBeta(m, post = postBeta, supportLevel = cfg$outputs$support_level_beta %||% 0.95,
                     param = "Sign", covNamesNumbers = c(TRUE, FALSE),
                     spNamesNumbers = c(m$ns <= 30, FALSE), cex = c(0.6, 0.6, 0.8))
      title(main = "Beta: environmental responses", line = 2.5, cex.main = 0.9)
      dev.off()
    }
  }, "Beta ", log_fun)

  # Gamma
  safe_try({
    if (m$nt > 1 && m$nc > 1) {
      postGamma <- Hmsc::getPostEstimate(m, parName = "Gamma")
      gamma_mean <- as.data.frame(postGamma$mean)
      gamma_support <- as.data.frame(postGamma$support)
      gamma_supportNeg <- as.data.frame(postGamma$supportNeg)
      write.csv(gamma_mean, file.path(tbl_dir, "parameter_estimates_Gamma_mean.csv"))
      write.csv(gamma_support, file.path(tbl_dir, "parameter_estimates_Gamma_Pr_positive.csv"))
      write.csv(gamma_supportNeg, file.path(tbl_dir, "parameter_estimates_Gamma_Pr_negative.csv"))

      pdf(file.path(plot_dir, "Gamma_plot.pdf"), width = 11, height = 8)
      Hmsc::plotGamma(m, post = postGamma, supportLevel = cfg$outputs$support_level_gamma %||% 0.95,
                      param = "Sign", covNamesNumbers = c(TRUE, FALSE),
                      trNamesNumbers = c(m$nt < 21, FALSE), cex = c(0.6,0.6,0.8))
      title(main = "Gamma: trait effects on environmental responses", line = 2.5, cex.main = 0.9)
      dev.off()
    }
  }, "Gamma ", log_fun)

  # Variance partitioning
  safe_try({
    if (m$nr + m$nc > 1 && m$ns > 1) {
      predY <- Hmsc::computePredictedValues(m)
      VP <- Hmsc::computeVariancePartitioning(m)
      MF <- Hmsc::evaluateModelFit(hM = m, predY = predY)
      write.csv(VP$vals, file.path(tbl_dir, "variance_partitioning_explained.csv"))
      if (!is.null(VP$R2T$Beta)) write.csv(VP$R2T$Beta, file.path(tbl_dir, "variance_partitioning_R2T_Beta.csv"))
      if (!is.null(VP$R2T$Y)) write.csv(VP$R2T$Y, file.path(tbl_dir, "variance_partitioning_R2T_Y.csv"))

      pdf(file.path(plot_dir, "variance_partitioning.pdf"), width = 12, height = 7)
      Hmsc::plotVariancePartitioning(hM = m, VP = VP, main = "Proportion of explained variance",
                                     args.leg = list(bg = "white", cex = 0.7))
      dev.off()
    }
  }, "", log_fun)

  # Omega associations
  safe_try({
    if (m$nr > 0 && m$ns > 1) {
      assoc <- Hmsc::computeAssociations(m)
      saveRDS(assoc, file.path(outdir, "results", "species_associations.rds"))
      for (r in seq_along(assoc)) {
        rn <- names(m$ranLevels)[[r]] %||% paste0("randomLevel", r)
        write.csv(assoc[[r]]$mean, file.path(tbl_dir, paste0("Omega_", rn, "_mean.csv")))
        write.csv(assoc[[r]]$support, file.path(tbl_dir, paste0("Omega_", rn, "_Pr_positive.csv")))

        if (requireNamespace("corrplot", quietly = TRUE)) {
          toPlot <- ((assoc[[r]]$support > (cfg$outputs$support_level_omega %||% 0.9)) +
                       (assoc[[r]]$support < (1 - (cfg$outputs$support_level_omega %||% 0.9))) > 0) *
            sign(assoc[[r]]$mean)
          pdf(file.path(plot_dir, paste0("Omega_associations_", rn, ".pdf")), width = 10, height = 10)
          corrplot::corrplot(toPlot, method = "color",
                             col = grDevices::colorRampPalette(c("blue", "white", "red"))(3),
                             mar = c(0,0,1,0),
                             main = paste0("Species associations: ", rn))
          dev.off()
        }
      }
    }
  }, "Omega ", log_fun)
}

prediction_gradients_rich <- function(m, outdir, cfg, log_fun) {
  log_fun("  predictions.pdf...")
  covariates <- character()
  if (identical(as.character(m$XFormula), "~.")) {
    covariates <- colnames(m$XData)
  } else {
    covariates <- all.vars(m$XFormula)
  }
  covariates <- covariates[covariates %in% colnames(m$XData)]
  if (length(covariates) == 0) return(invisible(NULL))

  pdf(file.path(outdir, "plots", "predictions_environmental_gradients.pdf"), width = 10, height = 8)
  on.exit(dev.off(), add = TRUE)

  # choose example species: presence-absence chooses prevalence near 0.5; otherwise mean max
  ex.sp <- which.max(colMeans(m$Y, na.rm = TRUE))
  safe_try({
    if (!is.null(m$distr) && m$distr[1,1] == 2) ex.sp <- which.min(abs(colMeans(m$Y, na.rm=TRUE) - 0.5))
  }, "", log_fun)

  for (covariate in covariates) {
    safe_try({
      Gradient <- Hmsc::constructGradient(m, focalVariable = covariate)
      Gradient2 <- Hmsc::constructGradient(m, focalVariable = covariate, non.focalVariables = 1)
      predY <- Hmsc::predict.Hmsc(m, Gradient = Gradient, expected = TRUE)
      predY2 <- Hmsc::predict.Hmsc(m, Gradient = Gradient2, expected = TRUE)

      par(mfrow = c(2,1))
      pl <- Hmsc::plotGradient(m, Gradient, pred = predY, yshow = 0, measure = "S",
                               showData = TRUE, main = paste0(covariate, ": summed response, total effect"))
      if (inherits(pl, "ggplot")) print(pl + ggplot2::labs(title = paste0(covariate, ": summed response, total effect")))
      pl <- Hmsc::plotGradient(m, Gradient2, pred = predY2, yshow = 0, measure = "S",
                               showData = TRUE, main = paste0(covariate, ": summed response, marginal effect"))
      if (inherits(pl, "ggplot")) print(pl + ggplot2::labs(title = paste0(covariate, ": summed response, marginal effect")))

      par(mfrow = c(2,1))
      yshow <- if (!is.null(m$distr) && m$distr[1,1] == 2) c(-0.1, 1.1) else 0
      pl <- Hmsc::plotGradient(m, Gradient, pred = predY, yshow = yshow, measure = "Y", index = ex.sp,
                               showData = TRUE, main = paste0(covariate, ": example species total effect"))
      if (inherits(pl, "ggplot")) print(pl + ggplot2::labs(title = paste0(covariate, ": example species total effect")))
      pl <- Hmsc::plotGradient(m, Gradient2, pred = predY2, yshow = yshow, measure = "Y", index = ex.sp,
                               showData = TRUE, main = paste0(covariate, ": example species marginal effect"))
      if (inherits(pl, "ggplot")) print(pl + ggplot2::labs(title = paste0(covariate, ": example species marginal effect")))
    }, paste0(" ", covariate), log_fun)
  }
}

make_html_report <- function(outdir, cfg) {
  report <- file.path(outdir, "JSDM Studio_report.html")
  files <- list.files(outdir, recursive = TRUE)
  li <- paste0("<li>", htmltools::htmlEscape(files), "</li>", collapse = "\n")
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'><title>JSDM Studio Report</title>",
    "<style>body{font-family:Arial,'Microsoft YaHei',sans-serif;background:#f8fafc;color:#1f2937;margin:40px;}",
    ".card{background:white;border:1px solid #e5e7eb;border-radius:18px;padding:22px;margin:16px 0;box-shadow:0 10px 30px rgba(15,23,42,.08)}",
    "h1{color:#155e63}.pill{display:inline-block;background:#ccfbf1;color:#134e4a;padding:5px 10px;border-radius:999px;margin:3px}</style></head><body>",
    "<h1> JSDM Studio </h1>",
    "<div class='card'><h2></h2>",
    "<span class='pill'>Hmsc</span><span class='pill'>Shiny GUI</span><span class='pill'>Reproducible workflow</span>",
    "<p> JSDM Studio  used_config.yml</p></div>",
    "<div class='card'><h2></h2><ul>", li, "</ul></div>",
    "<div class='card'><h2></h2><ol>",
    "<li> MCMC_convergence.txt  MCMC_traceplots.pdf</li>",
    "<li> model_fit*.csv  model_fit_explanatory_vs_predictive.pdf</li>",
    "<li> Beta_plot.pdfGamma_plot.pdfvariance_partitioning.pdf  Omega </li>",
    "<li> samples/transient used_config.yml</li>",
    "</ol></div></body></html>"
  )
  writeLines(html, report, useBytes = TRUE)
}

fit_hmsc_workflow <- function(Y, XData, TrData = NULL, phyloTree = NULL, studyDesign = NULL,
                              coordinates = NULL, cfg = list(), outdir = NULL,
                              log_fun = function(...) message(...),
                              progress_fun = function(percent, stage, step = NULL, detail = "") {}) {
  if (is.null(outdir)) outdir <- make_output_dir()
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  set.seed(cfg$mcmc$seed %||% 123)
  progress_fun(3, "", 1, " YXDatatraits ")

  Y <- as.matrix(Y)
  suppressWarnings(storage.mode(Y) <- "numeric")
  validate_hmsc_inputs(Y, XData, TrData, studyDesign, coordinates)
  save_input_summaries(Y, XData, TrData, studyDesign, coordinates, outdir)
  progress_fun(8, "", 1, " tables/")

  log_fun(" ...")
  rstuff <- build_random_levels(
    mode = cfg$random_effects$mode %||% "sample",
    Y = Y, studyDesign = studyDesign, coordinates = coordinates,
    random_effect_column = cfg$random_effects$random_effect_column %||% "sample",
    longitude_column = cfg$random_effects$longitude_column %||% "x",
    latitude_column = cfg$random_effects$latitude_column %||% "y",
    spatial_method = cfg$random_effects$spatial_method %||% "Full",
    nNeighbours = cfg$random_effects$nNeighbours %||% 10
  )

  args <- list(
    Y = Y,
    XData = XData,
    XFormula = as.formula(cfg$model$XFormula %||% "~ ."),
    distr = cfg$model$distr %||% "normal"
  )
  if (!is.null(rstuff$ranLevels)) {
    args$studyDesign <- rstuff$studyDesign
    args$ranLevels <- rstuff$ranLevels
  }
  if (isTRUE(cfg$model$use_traits) && !is.null(TrData)) {
    args$TrData <- TrData[colnames(Y), , drop = FALSE]
    args$TrFormula <- as.formula(cfg$model$TrFormula %||% "~ .")
  }
  if (isTRUE(cfg$model$use_phylogeny) && !is.null(phyloTree)) {
    args$phyloTree <- phyloTree
  }

  progress_fun(12, " Hmsc ", 2, " XFormulatraits/phylogeny ")
  log_fun("  Hmsc ...")
  m <- do.call(Hmsc::Hmsc, args)
  model_structure_output(m, outdir)

  # Save unfitted model like the S1 template idea
  models <- list(main_model = m)
  save(models, file = file.path(outdir, "models", "unfitted_models.RData"))
  saveRDS(m, file.path(outdir, "models", "unfitted_model.rds"))

  progress_fun(20, "MCMC ", 3, " sampleMcmcHmsc ")
  log_fun("  MCMC...")
  samples <- cfg$mcmc$samples %||% 100
  thin <- cfg$mcmc$thin %||% 1
  nChains <- cfg$mcmc$nChains %||% 2
  transient <- cfg$mcmc$transient %||% 50
  nParallel <- cfg$mcmc$nParallel %||% 1

  total_iter <- transient + samples * thin
  log_fun(paste0(" MCMC  ", total_iter,
                 " ", nChains,
                 "thin=", thin,
                 "samples=", samples,
                 "transient=", transient,
                 ""))
  log_fun(" Hmsc  sampleMcmc GUI verbose ")

  m <- Hmsc::sampleMcmc(
    m, samples = samples, transient = transient, thin = thin,
    nChains = nChains, nParallel = nParallel,
    verbose = cfg$mcmc$verbose %||% 50
  )

  # Save fitted model like S2 naming convention
  models <- list(main_model = m)
  fitted_name <- paste0("models_thin_", thin, "_samples_", samples, "_chains_", nChains, ".RData")
  save(models, file = file.path(outdir, "models", fitted_name))
  saveRDS(m, file.path(outdir, "models", "hmsc_model.rds"))
  if (isTRUE(cfg$outputs$save_model)) saveRDS(m, file.path(outdir, "hmsc_model.rds"))
  progress_fun(45, "MCMC ", 3, "")

  results <- list(model = m)

  # Predictions and model fit
  progress_fun(50, "", 4, " predicted values ")
  if (isTRUE(cfg$outputs$compute_predicted_values) || isTRUE(cfg$outputs$evaluate_model_fit)) {
    log_fun(" computePredictedValues + evaluateModelFit...")
    predY <- Hmsc::computePredictedValues(m)
    saveRDS(predY, file.path(outdir, "results", "predicted_values.rds"))
    results$predY <- predY

    if (isTRUE(cfg$outputs$evaluate_model_fit)) {
      MF <- Hmsc::evaluateModelFit(hM = m, predY = predY)
      saveRDS(MF, file.path(outdir, "results", "model_fit_explanatory.rds"))
      capture.output(print(MF), file = file.path(outdir, "results", "model_fit_explanatory.txt"))
      for (nm in names(MF)) {
        if (is.numeric(MF[[nm]])) write.csv(data.frame(species = m$spNames, value = MF[[nm]]),
                                            file.path(outdir, "tables", paste0("model_fit_explanatory_", nm, ".csv")),
                                            row.names = FALSE)
      }
      results$MF <- MF
    }
  }

  progress_fun(62, "/WAIC", 5, "")
  if (isTRUE(cfg$outputs$compute_cv)) {
    log_fun(" createPartition + computePredictedValues(partition=...)")
    MFCV <- safe_try({
      partition <- Hmsc::createPartition(m, nfolds = cfg$outputs$nfolds %||% 2)
      cvpreds <- Hmsc::computePredictedValues(m, partition = partition, nParallel = nParallel)
      MFCV <- Hmsc::evaluateModelFit(hM = m, predY = cvpreds)
      saveRDS(cvpreds, file.path(outdir, "results", "predicted_values_cross_validation.rds"))
      saveRDS(MFCV, file.path(outdir, "results", "model_fit_cross_validation.rds"))
      capture.output(print(MFCV), file = file.path(outdir, "results", "model_fit_cross_validation.txt"))
      for (nm in names(MFCV)) {
        if (is.numeric(MFCV[[nm]])) write.csv(data.frame(species = m$spNames, value = MFCV[[nm]]),
                                              file.path(outdir, "tables", paste0("model_fit_cross_validation_", nm, ".csv")),
                                              row.names = FALSE)
      }
      if (!is.null(results$MF)) plot_model_fit_compare(results$MF, MFCV, outdir)
      MFCV
    }, "", log_fun)
    results$MFCV <- MFCV
  }

  if (isTRUE(cfg$outputs$compute_waic)) {
    log_fun("  WAIC...")
    WAIC <- safe_try(Hmsc::computeWAIC(m), "WAIC", log_fun)
    if (!is.null(WAIC)) {
      saveRDS(WAIC, file.path(outdir, "results", "WAIC.rds"))
      capture.output(print(WAIC), file = file.path(outdir, "results", "WAIC.txt"))
    }
    results$WAIC <- WAIC
  }

  progress_fun(72, "MCMC ", 6, " MCMC_convergence.txt  traceplots")
  if (isTRUE(cfg$outputs$compute_diagnostics)) {
    safe_try(mcmc_diagnostics_rich(m, outdir, cfg, log_fun), "MCMC ", log_fun)
  }

  progress_fun(80, "//Omega", 7, " BetaGamma")
  if (isTRUE(cfg$outputs$compute_parameters)) {
    safe_try(parameter_estimates_rich(m, outdir, cfg, log_fun), "", log_fun)
  }

  progress_fun(90, "", 8, " environmental gradient predictions")
  if (isTRUE(cfg$outputs$make_predictions)) {
    safe_try(prediction_gradients_rich(m, outdir, cfg, log_fun), "", log_fun)
  }

  # old single outputs retained for compatibility
  if (isTRUE(cfg$outputs$compute_variance_partitioning)) {
    safe_try({
      vp <- Hmsc::computeVariancePartitioning(m)
      saveRDS(vp, file.path(outdir, "results", "variance_partitioning.rds"))
      write.csv(vp$vals, file.path(outdir, "tables", "variance_partitioning_vals.csv"))
    }, " RDS/CSV", log_fun)
  }

  if (isTRUE(cfg$outputs$compute_associations)) {
    safe_try({
      assoc <- Hmsc::computeAssociations(m)
      saveRDS(assoc, file.path(outdir, "results", "species_associations.rds"))
    }, " RDS", log_fun)
  }

  progress_fun(96, "", 9, " HTML ")
  if (isTRUE(cfg$outputs$make_report)) {
    safe_try(make_html_report(outdir, cfg), "HTML ", log_fun)
  }

  writeLines(c(
    "JSDM Studio rich output run completed.",
    paste("Time:", Sys.time()),
    paste("Output directory:", outdir),
    "See used_config.yml, results/, tables/, plots/, models/."
  ), file.path(outdir, "RUN_COMPLETE.txt"))

  return(list(results = results, outdir = outdir))
}
