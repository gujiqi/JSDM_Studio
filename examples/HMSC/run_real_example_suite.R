# Real Hmsc example suite for JSDM Studio.
# Run from the JSDMStudio project root:
# Rscript examples/Hmsc/run_real_example_suite.R

app_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

e <- new.env(parent = globalenv())
sys.source(file.path(app_dir, "app.R"), envir = e)

write_csv_pair <- function(x, outdir, filename) {
  if (is.null(x)) return(invisible(FALSE))
  write.csv(x, file.path(outdir, "data", filename), row.names = TRUE)
  write.csv(x, file.path(outdir, "inputs", filename), row.names = TRUE)
  TRUE
}

write_text_pair <- function(x, outdir, filename) {
  if (is.null(x)) return(invisible(FALSE))
  writeLines(x, file.path(outdir, "data", filename), useBytes = TRUE)
  writeLines(x, file.path(outdir, "inputs", filename), useBytes = TRUE)
  TRUE
}

write_table_pair <- function(x, outdir, filename) {
  if (is.null(x)) return(invisible(FALSE))
  write.csv(x, file.path(outdir, "data", filename), row.names = FALSE)
  write.csv(x, file.path(outdir, "inputs", filename), row.names = FALSE)
  TRUE
}

base_cfg <- function(case_id) {
  list(
    project_name = paste0("JSDMStudio_HMSC_SUITE_", case_id),
    engine = "Hmsc",
    data = list(phylogeny = NULL),
    model = list(
      distr = "probit",
      XFormula = "~ .",
      TrFormula = "~ .",
      use_traits = FALSE,
      use_phylogeny = FALSE,
      random_mode = "none",
      random_effect_column = "sample",
      spatial_method = "NNGP",
      nNeighbours = 3L,
      lon_col = "x",
      lat_col = "y",
      seed = 1234L,
      XScale = TRUE,
      TrScale = TRUE,
      YScale = FALSE,
      truncateNumberOfFactors = TRUE,
      Loff_file = "",
      ranLevelsUsed = "",
      C_file = "",
      use_XRRR = FALSE,
      ncRRR = 2L,
      XRRRFormula = "~ .",
      XRRRScale = TRUE,
      XRRR_file = "",
      random_level_type = "none",
      sMethod = "NNGP",
      random_N = 1L,
      longlat = FALSE,
      units_column = "sample",
      distMat_file = "",
      xData_file = "",
      sKnot_file = "",
      nfMin = 1L,
      nfMax = 4L,
      priors = list(setDefault = TRUE, a1 = NA, b1 = NA, a2 = NA, b2 = NA, alphapw = "")
    ),
    mcmc = list(
      samples = 5L,
      transient = 5L,
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
      convergence = list(showBeta = TRUE, showGamma = TRUE, showOmega = TRUE, maxOmega = 10L, showRho = TRUE, showAlpha = TRUE, effectiveSize = TRUE, gelmanPSRF = TRUE),
      plotting = list(var.part.order.explained = TRUE, var.part.order.raw = FALSE, show.sp.names.beta = FALSE, plotTree = FALSE, omega.order = "original", show.sp.names.omega = TRUE, plotBeta = TRUE, plotGamma = TRUE),
      predictions = list(species.list = "", trait.list = "", env.list = "", nfolds = 2L, partition_column = "sample", computeSAIR = FALSE)
    )
  )
}

make_binary_env <- function(n) {
  env <- data.frame(
    pH = seq(5.4, 7.2, length.out = n),
    moisture = seq(0.2, 0.9, length.out = n),
    substrate = rep(c("sand", "loam", "clay"), length.out = n)
  )
  rownames(env) <- paste0("site_", seq_len(n))
  env
}

make_hmsc_probit_case <- function() {
  n <- 12L
  env <- make_binary_env(n)
  score <- scale(env$pH)[, 1] + 0.8 * scale(env$moisture)[, 1]
  Y <- sapply(seq_len(4), function(j) as.integer(score + sin(seq_len(n) / (j + 1)) > median(score) + (j - 2.5) * 0.1))
  Y[1, ] <- 0L
  Y[n, ] <- 1L
  Y <- as.data.frame(Y)
  names(Y) <- paste0("sp_", seq_len(4))
  rownames(Y) <- rownames(env)
  list(id = "probit_linear_categorical", Y = Y, X = env, cfg = base_cfg("probit_linear_categorical"))
}

make_hmsc_sample_re_trait_case <- function() {
  n <- 12L
  env <- make_binary_env(n)
  score <- scale(env$pH)[, 1] - 0.4 * scale(env$moisture)[, 1]
  Y <- sapply(seq_len(4), function(j) as.integer(score + cos(seq_len(n) / (j + 1)) > median(score)))
  Y[1, ] <- 0L
  Y[n, ] <- 1L
  Y <- as.data.frame(Y)
  names(Y) <- paste0("sp_", seq_len(4))
  rownames(Y) <- rownames(env)
  study <- data.frame(plot = rep(paste0("plot_", 1:4), each = 3))
  rownames(study) <- rownames(env)
  traits <- data.frame(
    life_form = c("grass", "forb", "shrub", "grass"),
    height_mm = c(80, 120, 350, 95)
  )
  rownames(traits) <- names(Y)
  cfg <- base_cfg("sample_random_trait")
  cfg$model$random_mode <- "sample"
  cfg$model$random_effect_column <- "plot"
  cfg$model$use_traits <- TRUE
  cfg$model$TrFormula <- "~ life_form + height_mm"
  list(id = "sample_random_trait", Y = Y, X = env, Tr = traits, study = study, cfg = cfg)
}

make_hmsc_poisson_spatial_full_case <- function() {
  n <- 12L
  env <- data.frame(temp = seq(8, 18, length.out = n), moisture = seq(0.25, 0.85, length.out = n))
  rownames(env) <- paste0("plot_", seq_len(n))
  coord <- data.frame(x = rep(seq(0, 3), length.out = n), y = rep(seq(0, 2), each = 4, length.out = n))
  rownames(coord) <- rownames(env)
  lam <- exp(-0.7 + 0.05 * env$temp + 0.5 * env$moisture)
  Y <- sapply(seq_len(3), function(j) as.integer(round(lam + (seq_len(n) %% (j + 2)))))
  Y <- as.data.frame(Y)
  names(Y) <- paste0("sp_", seq_len(3))
  rownames(Y) <- rownames(env)
  cfg <- base_cfg("poisson_spatial_full")
  cfg$model$distr <- "poisson"
  cfg$model$random_mode <- "spatial_full"
  cfg$model$XFormula <- "~ temp + moisture"
  cfg$outputs$predictions$nfolds <- 2L
  list(id = "poisson_spatial_full", Y = Y, X = env, coord = coord, cfg = cfg)
}

make_hmsc_poisson_spatial_nngp_case <- function() {
  case <- make_hmsc_poisson_spatial_full_case()
  case$id <- "poisson_spatial_nngp"
  case$cfg <- base_cfg("poisson_spatial_nngp")
  case$cfg$model$distr <- "poisson"
  case$cfg$model$random_mode <- "spatial_nngp"
  case$cfg$model$spatial_method <- "NNGP"
  case$cfg$model$nNeighbours <- 3L
  case$cfg$model$XFormula <- "~ temp + moisture"
  case$cfg$outputs$predictions$nfolds <- 2L
  case
}

make_hmsc_normal_spatial_gpp_case <- function() {
  case <- make_hmsc_normal_case()
  coord <- data.frame(x = seq(0, 1, length.out = nrow(case$Y)), y = sin(seq_len(nrow(case$Y)) / 2))
  rownames(coord) <- rownames(case$Y)
  case$id <- "normal_spatial_gpp"
  case$coord <- coord
  case$cfg <- base_cfg("normal_spatial_gpp")
  case$cfg$model$distr <- "normal"
  case$cfg$model$random_mode <- "spatial_gpp"
  case$cfg$model$spatial_method <- "GPP"
  case$cfg$model$XFormula <- "~ temp + nutrient + habitat"
  case$cfg$outputs$predictions$nfolds <- 2L
  case
}

make_hmsc_normal_case <- function() {
  n <- 12L
  env <- data.frame(temp = seq(-1, 1, length.out = n), nutrient = sin(seq_len(n) / 3), habitat = rep(c("dry", "wet"), length.out = n))
  rownames(env) <- paste0("sample_", seq_len(n))
  Y <- data.frame(
    sp_1 = 0.5 + 0.8 * env$temp - 0.2 * env$nutrient,
    sp_2 = -0.2 + 0.4 * env$nutrient + 0.3 * (env$habitat == "wet"),
    sp_3 = 0.1 - 0.5 * env$temp + 0.2 * env$nutrient
  )
  rownames(Y) <- rownames(env)
  cfg <- base_cfg("normal_linear")
  cfg$model$distr <- "normal"
  cfg$model$XFormula <- "~ temp + nutrient + habitat"
  cfg$outputs$predictions$nfolds <- 2L
  list(id = "normal_linear", Y = Y, X = env, cfg = cfg)
}

make_hmsc_time_phylogeny_traits_case <- function() {
  n <- 16L
  env <- make_binary_env(n)
  coord <- data.frame(x = seq(0, 4, length.out = n), y = sin(seq_len(n) / 3))
  rownames(coord) <- rownames(env)
  idx <- seq_len(n)
  Y <- data.frame(
    sp_1 = 0.3 + 0.4 * scale(env$pH)[, 1] + 0.25 * sin(idx / 2),
    sp_2 = -0.2 + 0.5 * scale(env$moisture)[, 1] - 0.15 * cos(idx / 3),
    sp_3 = 0.1 - 0.2 * scale(env$pH)[, 1] + 0.35 * (env$substrate == "loam"),
    sp_4 = -0.1 + 0.3 * scale(env$moisture)[, 1] + 0.2 * (env$substrate == "clay")
  )
  rownames(Y) <- rownames(env)
  traits <- data.frame(
    life_form = c("grass", "forb", "shrub", "tree"),
    height_mm = c(70, 130, 420, 900),
    dispersal = c("wind", "animal", "animal", "wind")
  )
  rownames(traits) <- names(Y)
  cfg <- base_cfg("normal_time_phylogeny_traits")
  cfg$model$distr <- "normal"
  cfg$model$use_traits <- TRUE
  cfg$model$TrFormula <- "~ life_form + height_mm + dispersal"
  cfg$model$use_phylogeny <- TRUE
  cfg$model$phylogeny_type <- "time_calibrated_newick"
  cfg$model$random_mode <- "spatial_full"
  cfg$data$phylogeny <- "time_tree.nwk"
  cfg$model$C_file <- "time_tree_C.csv"
  cfg$outputs$plotting$plotTree <- TRUE
  newick <- "((sp_1:0.30,sp_2:0.30):0.20,(sp_3:0.25,sp_4:0.25):0.25);"
  Cmat <- matrix(
    c(1.0, 0.4, 0.0, 0.0,
      0.4, 1.0, 0.0, 0.0,
      0.0, 0.0, 1.0, 0.5,
      0.0, 0.0, 0.5, 1.0),
    nrow = 4L, byrow = TRUE,
    dimnames = list(names(Y), names(Y))
  )
  list(id = "normal_time_phylogeny_traits", Y = Y, X = env, Tr = traits, coord = coord, phylo_newick = newick, phylo_filename = "time_tree.nwk", C = Cmat, C_filename = "time_tree_C.csv", cfg = cfg, run_exported = TRUE)
}

make_hmsc_taxonomy_tree_spatial_traits_case <- function() {
  n <- 12L
  env <- data.frame(temp = seq(4, 16, length.out = n), moisture = seq(0.15, 0.95, length.out = n), substrate = rep(c("peat", "sand", "loam"), length.out = n))
  rownames(env) <- paste0("tax_site_", seq_len(n))
  Y <- data.frame(
    sp_1 = 0.2 + 0.3 * env$temp - 0.4 * env$moisture,
    sp_2 = 0.4 - 0.2 * env$temp + 0.7 * env$moisture,
    sp_3 = -0.1 + 0.1 * env$temp + 0.1 * (env$substrate == "peat"),
    sp_4 = 0.3 + 0.5 * env$moisture - 0.2 * (env$substrate == "sand")
  )
  rownames(Y) <- rownames(env)
  coord <- data.frame(x = seq(0, 3.3, length.out = n), y = cos(seq_len(n) / 3))
  rownames(coord) <- rownames(env)
  traits <- data.frame(
    photosynthetic_pathway = c("C3", "C3", "C4", "C3"),
    seed_mass_mg = c(1.2, 2.5, 0.8, 4.0)
  )
  rownames(traits) <- names(Y)
  taxonomy <- data.frame(
    species = names(Y),
    kingdom = rep("Plantae", 4),
    phylum = rep("Tracheophyta", 4),
    class = c("Magnoliopsida", "Magnoliopsida", "Liliopsida", "Magnoliopsida"),
    order = c("Asterales", "Asterales", "Poales", "Rosales"),
    family = c("Asteraceae", "Asteraceae", "Poaceae", "Rosaceae"),
    genus = c("Aster", "Solidago", "Poa", "Rosa"),
    species_rank = c("sp_1", "sp_2", "sp_3", "sp_4"),
    stringsAsFactors = FALSE
  )
  names(taxonomy)[names(taxonomy) == "species_rank"] <- "species_name"
  cfg <- base_cfg("normal_taxonomy_tree_spatial_traits")
  cfg$model$distr <- "normal"
  cfg$model$XFormula <- "~ temp + moisture + substrate"
  cfg$model$random_mode <- "spatial_full"
  cfg$model$use_traits <- TRUE
  cfg$model$TrFormula <- "~ photosynthetic_pathway + seed_mass_mg"
  cfg$model$use_phylogeny <- TRUE
  cfg$model$phylogeny_type <- "taxonomy_table"
  cfg$data$phylogeny <- "taxonomy_tree.csv"
  cfg$outputs$plotting$plotTree <- FALSE
  list(id = "normal_taxonomy_tree_spatial_traits", Y = Y, X = env, Tr = traits, coord = coord, taxonomy = taxonomy, taxonomy_filename = "taxonomy_tree.csv", cfg = cfg, run_exported = TRUE)
}

run_case <- function(case) {
  outdir <- e$make_engine_run_dir("Hmsc", paste0("JSDMStudio_HMSC_SUITE_", case$id))
  write_csv_pair(case$Y, outdir, "Y.csv")
  write_csv_pair(case$X, outdir, "XData.csv")
  write_csv_pair(case$Tr %||% NULL, outdir, "TrData.csv")
  write_csv_pair(case$study %||% NULL, outdir, "studyDesign.csv")
  write_csv_pair(case$coord %||% NULL, outdir, "coordinates.csv")
  write_text_pair(case$phylo_newick %||% NULL, outdir, case$phylo_filename %||% "time_tree.nwk")
  write_csv_pair(case$C %||% NULL, outdir, case$C_filename %||% "C_phylo.csv")
  write_table_pair(case$taxonomy %||% NULL, outdir, case$taxonomy_filename %||% "taxonomy_tree.csv")
  yaml::write_yaml(case$cfg, file.path(outdir, "used_config.yml"))
  check <- e$validate_hmsc(case$Y, case$X, case$Tr %||% NULL, case$study %||% NULL, case$coord %||% NULL,
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
                           partition_column = case$cfg$outputs$predictions$partition_column)
  write.csv(data.frame(message = check$messages), file.path(outdir, "diagnostics", "data_check_messages.csv"), row.names = FALSE)
  fit <- e$run_hmsc_s1s7_pipeline(outdir, case$cfg, case$Y, case$X, case$Tr %||% NULL, case$study %||% NULL, case$coord %||% NULL, log_fun = message)
  exported_status <- "not_run"
  exported_error <- ""
  if (isTRUE(case$run_exported %||% TRUE)) {
    script <- file.path(outdir, "reproducible_script", "run_this_HMSC_analysis.R")
    exported_log <- file.path(outdir, "diagnostics", "exported_script_run.log")
    exported_err <- file.path(outdir, "diagnostics", "exported_script_run_error.log")
    code <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script), stdout = exported_log, stderr = exported_err)
    exported_status <- if (identical(code, 0L)) "ok" else "failed"
    if (!identical(code, 0L) && file.exists(exported_err)) exported_error <- paste(readLines(exported_err, warn = FALSE), collapse = "; ")
  }
  zipfile <- tryCatch(e$make_zip(outdir), error = function(e) NA_character_)
  required <- c(
    "models/unfitted_models.RData",
    "models/hmsc_model_main.rds",
    "tables/S1_defined_models.csv",
    "tables/S1_phylogeny_summary.csv",
    "tables/S2_fit_models.csv",
    "tables/S4_model_fit_summary.csv",
    "results/MCMC_convergence.txt",
    "results/MCMC_convergence.pdf",
    "results/parameter_estimates.pdf",
    "predictions/predicted_values_main.rds",
    "diagnostics/engine_status.json",
    "diagnostics/session_info.txt",
    "plots/plot_manifest.csv",
    "report/Hmsc_report.html",
    "standard/run_summary.csv",
    "standard/fit_metrics.csv",
    "standard/effects_long.csv",
    "standard/predictions_long.csv",
    "standard/associations_long.csv",
    "standard/output_manifest.csv"
  )
  present <- file.exists(file.path(outdir, required))
  data.frame(
    case_id = case$id,
    distr = case$cfg$model$distr,
    random_mode = case$cfg$model$random_mode,
    use_traits = isTRUE(case$cfg$model$use_traits),
    check_ok = isTRUE(check$ok),
    engine_status = fit$status,
    exported_script = exported_status,
    exported_error = exported_error,
    errors = paste(as.character(fit$errors %||% character()), collapse = "; "),
    warnings = paste(as.character(fit$warnings %||% character()), collapse = "; "),
    missing_required = paste(required[!present], collapse = "; "),
    outdir = outdir,
    zip = zipfile,
    zip_size = if (!is.na(zipfile) && file.exists(zipfile)) file.info(zipfile)$size else NA_real_,
    stringsAsFactors = FALSE
  )
}

cases <- list(
  make_hmsc_probit_case(),
  make_hmsc_sample_re_trait_case(),
  make_hmsc_poisson_spatial_full_case(),
  make_hmsc_poisson_spatial_nngp_case(),
  make_hmsc_normal_spatial_gpp_case(),
  make_hmsc_time_phylogeny_traits_case(),
  make_hmsc_taxonomy_tree_spatial_traits_case(),
  make_hmsc_normal_case()
)

summary_rows <- do.call(rbind, lapply(cases, run_case))
summary_file <- file.path(app_dir, "output", paste0("Hmsc_real_example_suite_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
write.csv(summary_rows, summary_file, row.names = FALSE)
print(summary_rows[, c("case_id", "distr", "random_mode", "use_traits", "check_ok", "engine_status", "exported_script", "missing_required")], row.names = FALSE)
cat("Suite summary:", summary_file, "\n")

failed <- subset(summary_rows, !check_ok | engine_status != "fitted" | exported_script == "failed" | nzchar(missing_required))
if (nrow(failed) > 0) {
  writeLines(c("Failed Hmsc suite cases:", failed$case_id))
  stop("One or more Hmsc real example suite cases failed.", call. = FALSE)
}
