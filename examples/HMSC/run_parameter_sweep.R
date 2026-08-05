# HMSC parameter sweep for JSDM Studio.
# This is intentionally small-data / quick-MCMC. It checks that visible HMSC
# GUI parameters are wired into real Hmsc runs with small generated data.
# Run from the JSDMStudio project root:
# Rscript examples/Hmsc/run_parameter_sweep.R

app_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

e <- new.env(parent = globalenv())
sys.source(file.path(app_dir, "app.R"), envir = e)

write_csv_pair <- function(x, outdir, filename, row_names = TRUE) {
  if (is.null(x)) return(invisible(FALSE))
  write.csv(x, file.path(outdir, "data", filename), row.names = row_names)
  write.csv(x, file.path(outdir, "inputs", filename), row.names = row_names)
  TRUE
}

write_text_pair <- function(x, outdir, filename) {
  if (is.null(x)) return(invisible(FALSE))
  writeLines(x, file.path(outdir, "data", filename), useBytes = TRUE)
  writeLines(x, file.path(outdir, "inputs", filename), useBytes = TRUE)
  TRUE
}

base_cfg <- function(case_id) {
  list(
    project_name = paste0("JSDMStudio_HMSC_PARAMETER_SWEEP_", case_id),
    engine = "Hmsc",
    data = list(phylogeny = NULL),
    model = list(
      distr = "normal",
      XFormula = "~ temp + moisture + habitat",
      TrFormula = "~ height + life_form",
      use_traits = FALSE,
      use_phylogeny = FALSE,
      random_mode = "none",
      random_effect_column = "plot",
      spatial_method = "NNGP",
      nNeighbours = 2L,
      lon_col = "x",
      lat_col = "y",
      seed = 777L,
      XScale = TRUE,
      TrScale = TRUE,
      YScale = FALSE,
      truncateNumberOfFactors = TRUE,
      Loff_file = "",
      ranLevelsUsed = "",
      C_file = "",
      use_XRRR = FALSE,
      ncRRR = 1L,
      XRRRFormula = "~ rrr1",
      XRRRScale = TRUE,
      XRRR_file = "",
      random_level_type = "none",
      sMethod = "NNGP",
      random_N = 2L,
      longlat = FALSE,
      units_column = "plot",
      distMat_file = "",
      xData_file = "",
      sKnot_file = "",
      nfMin = 1L,
      nfMax = 3L,
      priors = list(setDefault = TRUE, a1 = NA, b1 = NA, a2 = NA, b2 = NA, alphapw = "")
    ),
    mcmc = list(
      samples = 4L,
      transient = 3L,
      thin = 1L,
      nChains = 2L,
      nParallel = 1L,
      verbose = 0L,
      preset = "Quick test",
      initPar = "fixed effects",
      alignPost = TRUE,
      updater = list(GammaEta = TRUE, Beta = TRUE, Gamma = TRUE, Omega = TRUE),
      sample_prior = FALSE,
      pool_chains = FALSE
    ),
    outputs = list(
      save_model = TRUE,
      predicted = TRUE,
      fit = TRUE,
      cv = TRUE,
      waic = FALSE,
      diagnostics = TRUE,
      parameters = TRUE,
      variance_partitioning = TRUE,
      omega = TRUE,
      gradients = TRUE,
      beta_support = TRUE,
      gamma_support = TRUE,
      omega_support = TRUE,
      convergence = list(showBeta = TRUE, showGamma = TRUE, showOmega = TRUE, maxOmega = 5L, showRho = TRUE, showAlpha = TRUE, effectiveSize = TRUE, gelmanPSRF = TRUE),
      plotting = list(var.part.order.explained = TRUE, var.part.order.raw = FALSE, show.sp.names.beta = TRUE, plotTree = FALSE, omega.order = "original", show.sp.names.omega = TRUE, plotBeta = TRUE, plotGamma = TRUE),
      predictions = list(species.list = "", trait.list = "", env.list = "", nfolds = 2L, partition_column = "plot", computeSAIR = FALSE)
    )
  )
}

make_normal_data <- function(n = 10L, p = 3L) {
  X <- data.frame(
    temp = seq(-1, 1, length.out = n),
    moisture = seq(0.2, 0.9, length.out = n),
    habitat = rep(c("dry", "wet"), length.out = n)
  )
  rownames(X) <- paste0("site_", seq_len(n))
  Y <- data.frame(
    sp_1 = 0.2 + 0.6 * X$temp - 0.2 * X$moisture,
    sp_2 = -0.1 + 0.3 * X$moisture + 0.2 * (X$habitat == "wet"),
    sp_3 = 0.4 - 0.5 * X$temp + 0.1 * X$moisture
  )
  Y <- Y[, seq_len(p), drop = FALSE]
  rownames(Y) <- rownames(X)
  study <- data.frame(plot = rep(paste0("plot_", 1:5), length.out = n), year = rep(c("y1", "y2"), length.out = n))
  rownames(study) <- rownames(X)
  coord <- data.frame(x = seq(0, 1, length.out = n), y = cos(seq_len(n) / 3))
  rownames(coord) <- rownames(X)
  traits <- data.frame(height = c(10, 20, 35)[seq_len(p)], life_form = c("grass", "forb", "shrub")[seq_len(p)])
  rownames(traits) <- names(Y)
  list(Y = Y, X = X, study = study, coord = coord, Tr = traits)
}

make_binary_data <- function(n = 10L, p = 3L) {
  d <- make_normal_data(n, p)
  z <- scale(d$X$temp + d$X$moisture)[, 1]
  Y <- sapply(seq_len(p), function(j) as.integer(z + sin(seq_len(n) / (j + 1)) > median(z)))
  Y[1, ] <- 0L
  Y[n, ] <- 1L
  Y <- as.data.frame(Y)
  names(Y) <- paste0("sp_", seq_len(p))
  rownames(Y) <- rownames(d$X)
  d$Y <- Y
  d
}

case_constructor_offsets <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("constructor_offsets")
  cfg$model$XScale <- FALSE
  cfg$model$YScale <- TRUE
  cfg$model$truncateNumberOfFactors <- FALSE
  cfg$model$Loff_file <- "Loff.csv"
  cfg$mcmc$thin <- 2L
  cfg$mcmc$transient <- 4L
  cfg$outputs$predictions$env.list <- "temp,moisture"
  Loff <- matrix(0.05, nrow = nrow(d$Y), ncol = ncol(d$Y), dimnames = dimnames(as.matrix(d$Y)))
  c(d, list(id = "constructor_offsets", cfg = cfg, extras = list(Loff.csv = Loff), parameters = c("XScale", "YScale", "truncateNumberOfFactors", "Loff_file", "thin", "transient", "env.list")))
}

case_xrrr <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("xrrr")
  cfg$model$use_XRRR <- TRUE
  cfg$model$XRRR_file <- "XRRRData.csv"
  cfg$model$XRRRFormula <- "~ rrr1"
  cfg$model$XRRRScale <- FALSE
  cfg$model$ncRRR <- 1L
  XRRRData <- data.frame(rrr1 = seq(-1, 1, length.out = nrow(d$Y)))
  rownames(XRRRData) <- rownames(d$Y)
  c(d, list(id = "xrrr", cfg = cfg, extras = list(XRRRData.csv = XRRRData), parameters = c("use_XRRR", "XRRR_file", "XRRRFormula", "XRRRScale", "ncRRR")))
}

case_sample_priors_partition <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("sample_priors_partition")
  cfg$model$random_mode <- "sample"
  cfg$model$random_effect_column <- "plot"
  cfg$model$ranLevelsUsed <- "plot"
  cfg$model$nfMin <- 1L
  cfg$model$nfMax <- 4L
  cfg$model$priors <- list(setDefault = FALSE, a1 = 20, b1 = 1, a2 = 20, b2 = 1, alphapw = "")
  cfg$outputs$predictions$partition_column <- "plot"
  cfg$mcmc$updater <- list(GammaEta = FALSE, Beta = FALSE, Gamma = FALSE, Omega = FALSE)
  c(d, list(id = "sample_priors_partition", cfg = cfg, parameters = c("random_mode=sample", "random_effect_column", "ranLevelsUsed", "nfMin", "nfMax", "a1", "b1", "a2", "b2", "partition_column", "updater toggles")))
}

case_spatial_full_alphapw <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("spatial_full_alphapw")
  cfg$model$random_mode <- "spatial_full"
  cfg$model$lon_col <- "x"
  cfg$model$lat_col <- "y"
  cfg$model$longlat <- FALSE
  cfg$model$priors$alphapw <- "0,1,1,0.5"
  c(d, list(id = "spatial_full_alphapw", cfg = cfg, parameters = c("random_mode=spatial_full", "lon_col", "lat_col", "longlat", "alphapw")))
}

case_spatial_nngp <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("spatial_nngp")
  cfg$model$random_mode <- "spatial_nngp"
  cfg$model$nNeighbours <- 2L
  c(d, list(id = "spatial_nngp", cfg = cfg, parameters = c("random_mode=spatial_nngp", "nNeighbours")))
}

case_spatial_gpp_knots <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("spatial_gpp_knots")
  cfg$model$random_mode <- "spatial_gpp"
  cfg$model$sKnot_file <- "knots.csv"
  knots <- data.frame(x = c(0.18, 0.52, 0.84), y = c(-0.72, 0.03, 0.68))
  rownames(knots) <- paste0("knot_", seq_len(nrow(knots)))
  c(d, list(id = "spatial_gpp_knots", cfg = cfg, extras = list(knots.csv = knots), parameters = c("random_mode=spatial_gpp", "sKnot_file")))
}

case_advanced_units <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("advanced_units")
  cfg$model$random_level_type <- "unstructured units"
  cfg$model$units_column <- "plot"
  c(d, list(id = "advanced_units", cfg = cfg, parameters = c("random_level_type=unstructured units", "units_column")))
}

case_advanced_distmat <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("advanced_distmat")
  cfg$model$random_level_type <- "distance matrix distMat"
  cfg$model$units_column <- "plot"
  cfg$model$distMat_file <- "distMat.csv"
  cfg$outputs$omega <- FALSE
  cfg$outputs$gradients <- FALSE
  lev <- unique(as.character(d$study$plot))
  distMat <- as.matrix(dist(seq_along(lev)))
  rownames(distMat) <- colnames(distMat) <- lev
  c(d, list(id = "advanced_distmat", cfg = cfg, extras = list(distMat.csv = distMat), parameters = c("random_level_type=distance matrix distMat", "distMat_file", "units_column", "omega=FALSE for distMat quick test", "gradients=FALSE")))
}

case_advanced_xdata <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("advanced_xdata")
  cfg$model$random_level_type <- "covariate-dependent xData"
  cfg$model$units_column <- "plot"
  cfg$model$xData_file <- "random_xData.csv"
  cfg$model$nfMin <- 2L
  cfg$model$nfMax <- 2L
  cfg$outputs$cv <- FALSE
  cfg$outputs$diagnostics <- FALSE
  cfg$outputs$gradients <- FALSE
  lev <- unique(as.character(d$study$plot))
  random_xData <- data.frame(plot_elevation = seq(100, 500, length.out = length(lev)))
  rownames(random_xData) <- lev
  c(d, list(id = "advanced_xdata", cfg = cfg, extras = list(random_xData.csv = random_xData), parameters = c("random_level_type=covariate-dependent xData", "xData_file", "units_column", "nfMin=nfMax=2", "cv=FALSE", "diagnostics=FALSE", "gradients=FALSE")))
}

case_advanced_N <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("advanced_N")
  cfg$model$random_level_type <- "N only"
  cfg$model$random_N <- 3L
  cfg$model$units_column <- "sample_block"
  c(d, list(id = "advanced_N", cfg = cfg, parameters = c("random_level_type=N only", "random_N", "units_column")))
}

case_phylogeny_C_traits <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("phylogeny_C_traits")
  cfg$model$use_traits <- TRUE
  cfg$model$use_phylogeny <- TRUE
  cfg$model$C_file <- "C_phylo.csv"
  cfg$outputs$plotting$plotTree <- TRUE
  C <- matrix(c(1, 0.7, 0.3, 0.7, 1, 0.4, 0.3, 0.4, 1), 3, 3)
  rownames(C) <- colnames(C) <- names(d$Y)
  c(d, list(id = "phylogeny_C_traits", cfg = cfg, extras = list(C_phylo.csv = C), parameters = c("use_traits", "TrFormula", "TrScale", "use_phylogeny", "C_file", "plotTree")))
}

case_phylogeny_newick_outputs <- function() {
  d <- make_binary_data()
  cfg <- base_cfg("phylogeny_newick_outputs")
  cfg$model$distr <- "probit"
  cfg$model$XFormula <- "~ temp + moisture + habitat"
  cfg$model$use_phylogeny <- TRUE
  cfg$data$phylogeny <- "tree.nwk"
  cfg$outputs$waic <- TRUE
  cfg$outputs$beta_support <- FALSE
  cfg$outputs$gamma_support <- FALSE
  cfg$outputs$omega_support <- FALSE
  cfg$outputs$plotting$omega.order <- "AOE"
  cfg$outputs$predictions$species.list <- "sp_1,sp_2"
  cfg$outputs$predictions$trait.list <- "height"
  cfg$outputs$predictions$env.list <- "temp"
  cfg$outputs$predictions$computeSAIR <- TRUE
  newick <- "((sp_1:0.2,sp_2:0.2):0.2,sp_3:0.4);"
  c(d, list(id = "phylogeny_newick_outputs", cfg = cfg, newick = newick, parameters = c("distr=probit", "phylogeny_file", "waic", "beta_support", "gamma_support", "omega_support", "omega.order", "species.list", "trait.list", "env.list", "computeSAIR")))
}

case_sample_prior_pool <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("sample_prior_pool")
  cfg$mcmc$sample_prior <- TRUE
  cfg$mcmc$pool_chains <- TRUE
  cfg$mcmc$alignPost <- FALSE
  cfg$mcmc$nChains <- 2L
  cfg$mcmc$nParallel <- 1L
  cfg$outputs$predicted <- FALSE
  cfg$outputs$fit <- FALSE
  cfg$outputs$cv <- FALSE
  cfg$outputs$diagnostics <- FALSE
  cfg$outputs$gradients <- FALSE
  c(d, list(id = "sample_prior_pool", cfg = cfg, parameters = c("sample_prior", "pool_chains", "alignPost", "nChains", "nParallel", "predicted=FALSE", "fit=FALSE", "cv=FALSE", "diagnostics=FALSE", "gradients=FALSE")))
}

case_output_and_plot_switches <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("output_and_plot_switches")
  cfg$outputs$save_model <- FALSE
  cfg$outputs$predicted <- FALSE
  cfg$outputs$fit <- FALSE
  cfg$outputs$cv <- FALSE
  cfg$outputs$waic <- TRUE
  cfg$outputs$diagnostics <- FALSE
  cfg$outputs$parameters <- FALSE
  cfg$outputs$variance_partitioning <- FALSE
  cfg$outputs$omega <- FALSE
  cfg$outputs$gradients <- FALSE
  cfg$outputs$beta_support <- FALSE
  cfg$outputs$gamma_support <- FALSE
  cfg$outputs$omega_support <- FALSE
  cfg$outputs$convergence <- list(showBeta = FALSE, showGamma = FALSE, showOmega = FALSE, maxOmega = 2L, showRho = FALSE, showAlpha = FALSE, effectiveSize = FALSE, gelmanPSRF = FALSE)
  cfg$outputs$plotting <- list(var.part.order.explained = FALSE, var.part.order.raw = TRUE, show.sp.names.beta = TRUE, plotTree = FALSE, omega.order = "FPC", show.sp.names.omega = FALSE, plotBeta = FALSE, plotGamma = FALSE)
  cfg$outputs$predictions$nfolds <- 3L
  c(d, list(id = "output_and_plot_switches", cfg = cfg, parameters = c("save_model=FALSE", "predicted=FALSE", "fit=FALSE", "cv=FALSE", "waic=TRUE", "diagnostics=FALSE", "parameters=FALSE", "variance_partitioning=FALSE", "omega=FALSE", "gradients=FALSE", "support toggles FALSE", "convergence toggles FALSE", "plotting toggles", "maxOmega", "nfolds")))
}

case_omega_plot_controls <- function() {
  d <- make_normal_data()
  cfg <- base_cfg("omega_plot_controls")
  cfg$model$random_mode <- "sample"
  cfg$model$random_effect_column <- "plot"
  cfg$model$units_column <- "plot"
  cfg$model$XFormula <- "~ temp + moisture + habitat"
  cfg$outputs$predicted <- TRUE
  cfg$outputs$fit <- TRUE
  cfg$outputs$cv <- FALSE
  cfg$outputs$parameters <- TRUE
  cfg$outputs$variance_partitioning <- TRUE
  cfg$outputs$omega <- TRUE
  cfg$outputs$convergence$maxOmega <- 2L
  cfg$outputs$plotting$var.part.order.explained <- FALSE
  cfg$outputs$plotting$var.part.order.raw <- TRUE
  cfg$outputs$plotting$show.sp.names.beta <- TRUE
  cfg$outputs$plotting$omega.order <- "alphabetical"
  cfg$outputs$plotting$show.sp.names.omega <- FALSE
  cfg$outputs$plotting$plotBeta <- TRUE
  cfg$outputs$plotting$plotGamma <- FALSE
  c(d, list(id = "omega_plot_controls", cfg = cfg, parameters = c("omega=TRUE", "maxOmega", "omega.order=alphabetical", "show.sp.names.omega=FALSE", "show.sp.names.beta=TRUE", "plotBeta=TRUE", "plotGamma=FALSE", "var.part.order.raw=TRUE")))
}

run_case <- function(case) {
  outdir <- e$make_engine_run_dir("Hmsc", paste0("JSDMStudio_HMSC_PARAMETER_SWEEP_", case$id))
  write_csv_pair(case$Y, outdir, "Y.csv")
  write_csv_pair(case$X, outdir, "XData.csv")
  write_csv_pair(case$Tr %||% NULL, outdir, "TrData.csv")
  write_csv_pair(case$study %||% NULL, outdir, "studyDesign.csv")
  write_csv_pair(case$coord %||% NULL, outdir, "coordinates.csv")
  if (!is.null(case$extras)) {
    for (nm in names(case$extras)) write_csv_pair(case$extras[[nm]], outdir, nm)
  }
  write_text_pair(case$newick %||% NULL, outdir, case$cfg$data$phylogeny %||% "tree.nwk")
  yaml::write_yaml(case$cfg, file.path(outdir, "used_config.yml"))
  check <- e$validate_hmsc(
    case$Y, case$X, case$Tr %||% NULL, case$study %||% NULL, case$coord %||% NULL,
    distr = case$cfg$model$distr,
    XFormula = case$cfg$model$XFormula,
    TrFormula = case$cfg$model$TrFormula,
    use_traits = isTRUE(case$cfg$model$use_traits),
    use_phylogeny = isTRUE(case$cfg$model$use_phylogeny),
    random_mode = case$cfg$model$random_mode,
    random_effect_column = case$cfg$model$random_effect_column,
    spatial_method = case$cfg$model$spatial_method,
    nNeighbours = case$cfg$model$nNeighbours,
    lon_col = case$cfg$model$lon_col,
    lat_col = case$cfg$model$lat_col,
    samples = case$cfg$mcmc$samples,
    transient = case$cfg$mcmc$transient,
    thin = case$cfg$mcmc$thin,
    nChains = case$cfg$mcmc$nChains,
    nParallel = case$cfg$mcmc$nParallel,
    nfMin = case$cfg$model$nfMin,
    nfMax = case$cfg$model$nfMax,
    nfolds = case$cfg$outputs$predictions$nfolds,
    partition_column = case$cfg$outputs$predictions$partition_column
  )
  write.csv(data.frame(message = check$messages), file.path(outdir, "diagnostics", "data_check_messages.csv"), row.names = FALSE)
  fit <- e$run_hmsc_s1s7_pipeline(outdir, case$cfg, case$Y, case$X, case$Tr %||% NULL, case$study %||% NULL, case$coord %||% NULL, log_fun = message)
  exported_status <- "not_run"
  exported_error <- ""
  if (isTRUE(case$run_exported %||% TRUE) && identical(fit$status, "fitted")) {
    script <- file.path(outdir, "reproducible_script", "run_this_HMSC_analysis.R")
    exported_log <- file.path(outdir, "diagnostics", "exported_script_parameter_sweep.log")
    exported_err <- file.path(outdir, "diagnostics", "exported_script_parameter_sweep_error.log")
    code <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script), stdout = exported_log, stderr = exported_err)
    exported_status <- if (identical(code, 0L)) "ok" else "failed"
    if (!identical(code, 0L) && file.exists(exported_err)) exported_error <- paste(readLines(exported_err, warn = FALSE), collapse = "; ")
  }
  required <- c("diagnostics/engine_status.json", "diagnostics/session_info.txt", "tables/HMSC_parameter_audit.csv", "workflow_scripts/S1_define_models.R", "reproducible_script/run_this_HMSC_analysis.R")
  if (identical(fit$status, "fitted")) {
    required <- c(required, "models/unfitted_models.RData", "standard/run_summary.csv", "standard/effects_long.csv", "standard/predictions_long.csv", "plots/plot_manifest.csv", "report/Hmsc_report.html")
    if (isTRUE(case$cfg$outputs$save_model %||% TRUE)) required <- c(required, "models/hmsc_model_main.rds")
  } else {
    required <- c(required, "diagnostics/HMSC_S1S7_error.txt")
  }
  present <- file.exists(file.path(outdir, required))
  zipfile <- tryCatch(e$make_zip(outdir), error = function(e) NA_character_)
  data.frame(
    case_id = case$id,
    expected_status = case$expected_status %||% "fitted",
    engine_status = fit$status,
    exported_script = exported_status,
    exported_error = exported_error,
    check_ok = isTRUE(check$ok),
    parameters_tested = paste(case$parameters, collapse = "; "),
    missing_required = paste(required[!present], collapse = "; "),
    errors = paste(as.character(fit$errors %||% character()), collapse = "; "),
    warnings = paste(as.character(fit$warnings %||% character()), collapse = "; "),
    outdir = outdir,
    zip = zipfile,
    zip_size = if (!is.na(zipfile) && file.exists(zipfile)) file.info(zipfile)$size else NA_real_,
    stringsAsFactors = FALSE
  )
}

cases <- list(
  case_constructor_offsets(),
  case_xrrr(),
  case_sample_priors_partition(),
  case_spatial_full_alphapw(),
  case_spatial_nngp(),
  case_spatial_gpp_knots(),
  case_advanced_units(),
  case_advanced_distmat(),
  case_advanced_xdata(),
  case_advanced_N(),
  case_phylogeny_C_traits(),
  case_phylogeny_newick_outputs(),
  case_sample_prior_pool(),
  case_output_and_plot_switches(),
  case_omega_plot_controls()
)

summary_rows <- do.call(rbind, lapply(cases, run_case))
summary_file <- file.path(app_dir, "output", paste0("Hmsc_parameter_sweep_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
write.csv(summary_rows, summary_file, row.names = FALSE)
print(summary_rows[, c("case_id", "expected_status", "engine_status", "exported_script", "missing_required")], row.names = FALSE)
cat("Parameter sweep summary:", summary_file, "\n")

unexpected <- subset(summary_rows, engine_status != expected_status | (expected_status == "fitted" & exported_script == "failed") | nzchar(missing_required))
if (nrow(unexpected) > 0) {
  writeLines(c("Unexpected HMSC parameter sweep results:", unexpected$case_id))
  stop("One or more HMSC parameter sweep cases failed unexpectedly.", call. = FALSE)
}
