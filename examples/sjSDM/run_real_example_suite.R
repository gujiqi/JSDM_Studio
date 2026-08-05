# Real sjSDM example suite for JSDM Studio.
# Run from the JSDMStudio project root:
# Rscript examples/sjSDM/run_real_example_suite.R

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

base_cfg <- function(case_id) {
  list(
    project_name = paste0("JSDMStudio_SUITE_", case_id),
    engine = "sjSDM",
    model = list(
      family = "binomial_probit",
      env_model = "linear",
      env_formula = "~ .",
      spatial_model = "none",
      spatial_formula = "~ 0 + .",
      se = FALSE,
      iter = 5L,
      step_size = 4L,
      sampling = 100L,
      parallel = 0L,
      dtype = "float32",
      verbose = FALSE,
      seed = 1234L
    ),
    regularization_biotic = list(
      env_lambda = 0,
      env_alpha = 0.5,
      spatial_lambda = 0,
      spatial_alpha = 0.5,
      biotic_lambda = 0.01,
      biotic_alpha = 0.5,
      biotic_df = NA,
      on_diag = FALSE,
      reg_on_Cov = TRUE,
      inverse = FALSE,
      tune_regularization = FALSE,
      tune_steps = 0L,
      cv_k = 0L
    ),
    dnn_optimizer = list(
      hidden = "6,4",
      activation = "selu",
      dropout = 0,
      bias = TRUE,
      optimizer = "Adamax",
      learning_rate = 0.003,
      weight_decay = 0.001,
      device = "cpu"
    ),
    control = list(
      scheduler = 0L,
      lr_reduce_factor = 0.99,
      early_stopping_training = 0L,
      mixed = FALSE
    ),
    spatial_anova_metacommunity = list(
      generateSpatialEV = FALSE,
      spatial_ev_k = 4L,
      spatial_ev_threshold = 0,
      include_space_in_anova = FALSE,
      do_anova = TRUE,
      anova_samples = 100L,
      do_internal = TRUE,
      internal_fractions = "proportional",
      do_assembly = FALSE,
      assembly_predictor = ""
    ),
    outputs = list(
      real_fit = TRUE,
      predict = TRUE,
      Rsquared = TRUE,
      importance = TRUE,
      weights = TRUE,
      coef = TRUE,
      covariance_correlation = TRUE,
      residuals = TRUE,
      plots = TRUE,
      zip = TRUE,
      copy_inputs = TRUE,
      save_config = TRUE,
      report = TRUE
    )
  )
}

make_binomial_case <- function() {
  n <- 14L
  env <- data.frame(
    pH = seq(5.4, 7.2, length.out = n),
    moisture = seq(0.2, 0.85, length.out = n),
    canopy = round(seq(12, 83, length.out = n), 1),
    substrate = rep(c("sand", "loam", "clay"), length.out = n),
    elevation = seq(90, 220, length.out = n)
  )
  rownames(env) <- paste0("site_", seq_len(n))
  score <- scale(env$pH)[, 1] + 0.7 * scale(env$moisture)[, 1] - 0.02 * env$canopy
  Y <- sapply(seq_len(5), function(j) as.integer(score + sin(seq_len(n) / (j + 1)) + (j - 3) * 0.25 > quantile(score, 0.45)))
  Y[1, ] <- 0L
  Y[n, ] <- 1L
  Y <- as.data.frame(Y)
  names(Y) <- paste0("sp_", seq_len(5))
  rownames(Y) <- rownames(env)
  newdata <- env[c(2, 5, 9, 13), , drop = FALSE]
  cfg <- base_cfg("binomial_linear_categorical")
  cfg$spatial_anova_metacommunity$do_internal <- FALSE
  list(id = "binomial_linear_categorical", Y = Y, env = env, newdata = newdata, cfg = cfg)
}

make_poisson_spatial_case <- function() {
  n <- 16L
  env <- data.frame(
    temperature = seq(8, 19, length.out = n),
    forest = seq(0.15, 0.75, length.out = n),
    habitat = rep(c("wetland", "forest", "grassland", "edge"), length.out = n)
  )
  rownames(env) <- paste0("plot_", seq_len(n))
  spatial <- data.frame(x = rep(seq(0, 3), each = 4), y = rep(seq(0, 3), times = 4))
  rownames(spatial) <- rownames(env)
  lam <- exp(-0.6 + 0.06 * env$temperature + 0.5 * env$forest + 0.08 * spatial$x)
  Y <- sapply(seq_len(4), function(j) pmax(0L, as.integer(round(lam + (j - 1) * 0.4 + (seq_len(n) %% (j + 2))))))
  Y <- as.data.frame(Y)
  names(Y) <- paste0("sp_", seq_len(4))
  rownames(Y) <- rownames(env)
  cfg <- base_cfg("poisson_spatial_linear")
  cfg$model$family <- "poisson_log"
  cfg$model$spatial_model <- "linear"
  cfg$spatial_anova_metacommunity$include_space_in_anova <- TRUE
  cfg$spatial_anova_metacommunity$do_internal <- TRUE
  cfg$dnn_optimizer$optimizer <- "RMSprop"
  newdata <- env[c(3, 7, 11, 15), , drop = FALSE]
  new_spatial <- spatial[c(3, 7, 11, 15), , drop = FALSE]
  list(id = "poisson_spatial_linear", Y = Y, env = env, spatial = spatial, newdata = newdata, new_spatial = new_spatial, cfg = cfg)
}

make_gaussian_dnn_case <- function() {
  n <- 18L
  env <- data.frame(
    temp = seq(-1, 1, length.out = n),
    moisture = cos(seq_len(n) / 3),
    nutrient = sin(seq_len(n) / 4),
    disturbance = rep(c(0, 1), length.out = n)
  )
  rownames(env) <- paste0("sample_", seq_len(n))
  Y <- data.frame(
    sp_1 = 0.4 + 0.8 * env$temp - 0.2 * env$moisture,
    sp_2 = -0.1 + 0.5 * env$nutrient + 0.3 * env$disturbance,
    sp_3 = 0.2 - 0.6 * env$temp + 0.4 * env$moisture
  )
  rownames(Y) <- rownames(env)
  cfg <- base_cfg("gaussian_dnn_env")
  cfg$model$family <- "gaussian_identity"
  cfg$model$env_model <- "DNN"
  cfg$model$iter <- 6L
  cfg$model$step_size <- 6L
  cfg$dnn_optimizer$optimizer <- "SGD"
  cfg$dnn_optimizer$learning_rate <- 0.002
  cfg$outputs$importance <- FALSE
  cfg$spatial_anova_metacommunity$do_anova <- FALSE
  list(id = "gaussian_dnn_env", Y = Y, env = env, newdata = env[c(2, 6, 10, 14), , drop = FALSE], cfg = cfg)
}

make_binomial_logit_dnn_case <- function() {
  n <- 18L
  env <- data.frame(
    temp = seq(0.1, 1.8, length.out = n),
    moisture = seq(0.2, 0.9, length.out = n),
    shade = sin(seq_len(n) / 3),
    nutrient = cos(seq_len(n) / 4)
  )
  rownames(env) <- paste0("dnn_site_", seq_len(n))
  score <- scale(env$temp)[, 1] + 0.5 * scale(env$moisture)[, 1] - 0.25 * env$shade
  Y <- sapply(seq_len(4), function(j) as.integer(score + sin(seq_len(n) / (j + 1)) > median(score) + (j - 2.5) * 0.1))
  Y[1, ] <- 0L
  Y[n, ] <- 1L
  Y <- as.data.frame(Y)
  names(Y) <- paste0("sp_", seq_len(4))
  rownames(Y) <- rownames(env)
  cfg <- base_cfg("binomial_logit_dnn_scheduler")
  cfg$model$family <- "binomial_logit"
  cfg$model$env_model <- "DNN"
  cfg$model$iter <- 6L
  cfg$model$step_size <- 6L
  cfg$dnn_optimizer$activation <- "relu"
  cfg$dnn_optimizer$dropout <- 0.05
  cfg$dnn_optimizer$optimizer <- "AdaBound"
  cfg$control$scheduler <- 2L
  cfg$control$early_stopping_training <- 2L
  cfg$outputs$importance <- FALSE
  cfg$spatial_anova_metacommunity$do_anova <- FALSE
  list(id = "binomial_logit_dnn_scheduler", Y = Y, env = env, newdata = env[c(3, 8, 13, 17), , drop = FALSE], cfg = cfg)
}

make_spatial_dnn_case <- function() {
  n <- 16L
  env <- data.frame(
    pH = seq(5.5, 7.0, length.out = n),
    moisture = rep(seq(0.25, 0.85, length.out = 4), each = 4),
    canopy = rep(seq(20, 80, length.out = 4), times = 4)
  )
  rownames(env) <- paste0("sdnn_", seq_len(n))
  spatial <- data.frame(x = rep(seq(0, 3), each = 4), y = rep(seq(0, 3), times = 4))
  rownames(spatial) <- rownames(env)
  score <- scale(env$pH)[, 1] + 0.4 * spatial$x - 0.2 * spatial$y
  Y <- sapply(seq_len(4), function(j) as.integer(score + cos(seq_len(n) / (j + 1)) > median(score)))
  Y[1, ] <- 0L
  Y[n, ] <- 1L
  Y <- as.data.frame(Y)
  names(Y) <- paste0("sp_", seq_len(4))
  rownames(Y) <- rownames(env)
  cfg <- base_cfg("binomial_spatial_dnn_madgrad")
  cfg$model$family <- "binomial_probit"
  cfg$model$spatial_model <- "DNN"
  cfg$model$spatial_formula <- "~ 0 + ."
  cfg$model$iter <- 6L
  cfg$model$step_size <- 4L
  cfg$dnn_optimizer$hidden <- "5"
  cfg$dnn_optimizer$activation <- "tanh"
  cfg$dnn_optimizer$optimizer <- "madgrad"
  cfg$dnn_optimizer$dropout <- 0.02
  cfg$spatial_anova_metacommunity$include_space_in_anova <- TRUE
  cfg$spatial_anova_metacommunity$do_internal <- TRUE
  list(id = "binomial_spatial_dnn_madgrad", Y = Y, env = env, spatial = spatial, newdata = env[c(2, 6, 10, 14), , drop = FALSE], new_spatial = spatial[c(2, 6, 10, 14), , drop = FALSE], cfg = cfg)
}

make_nbinom_intercept_case <- function() {
  n <- 15L
  env <- data.frame(dummy = rep(1, n))
  rownames(env) <- paste0("quad_", seq_len(n))
  Y <- data.frame(
    sp_1 = c(0, 1, 1, 2, 0, 3, 1, 4, 2, 1, 5, 2, 3, 4, 1),
    sp_2 = c(2, 1, 0, 1, 3, 2, 4, 1, 5, 3, 2, 4, 6, 3, 2),
    sp_3 = c(1, 3, 2, 4, 1, 5, 3, 6, 2, 4, 5, 7, 3, 6, 4),
    sp_4 = c(0, 0, 1, 0, 2, 1, 3, 2, 1, 4, 2, 5, 3, 4, 2)
  )
  rownames(Y) <- rownames(env)
  cfg <- base_cfg("nbinom_intercept_only")
  cfg$model$family <- "nbinom"
  cfg$model$env_model <- "intercept-only"
  cfg$model$env_formula <- "~ 1"
  cfg$model$iter <- 5L
  cfg$model$step_size <- 5L
  cfg$outputs$importance <- FALSE
  cfg$spatial_anova_metacommunity$do_anova <- FALSE
  list(id = "nbinom_intercept_only", Y = Y, env = env, cfg = cfg)
}

make_spatial_ev_case <- function() {
  n <- 16L
  env <- data.frame(
    pH = seq(5.8, 7.1, length.out = n),
    moisture = seq(0.3, 0.7, length.out = n),
    habitat = rep(c("ridge", "valley"), length.out = n)
  )
  rownames(env) <- paste0("cell_", seq_len(n))
  spatial <- data.frame(x = rep(seq(0, 3), each = 4), y = rep(seq(0, 3), times = 4))
  rownames(spatial) <- rownames(env)
  lin <- env$pH - mean(env$pH) + 0.4 * spatial$x - 0.3 * spatial$y
  Y <- sapply(seq_len(4), function(j) as.integer(lin + cos(seq_len(n) / j) > median(lin)))
  Y[1, ] <- 0L
  Y[n, ] <- 1L
  Y <- as.data.frame(Y)
  names(Y) <- paste0("sp_", seq_len(4))
  rownames(Y) <- rownames(env)
  cfg <- base_cfg("binomial_spatial_ev")
  cfg$model$spatial_model <- "eigenvectors"
  cfg$spatial_anova_metacommunity$generateSpatialEV <- TRUE
  cfg$spatial_anova_metacommunity$spatial_ev_k <- 4L
  cfg$spatial_anova_metacommunity$include_space_in_anova <- TRUE
  cfg$outputs$predict <- TRUE
  list(id = "binomial_spatial_ev", Y = Y, env = env, spatial = spatial, cfg = cfg)
}

make_spatial_traits_assembly_se_case <- function() {
  n <- 14L
  env <- data.frame(
    pH = seq(5.6, 7.2, length.out = n),
    moisture = seq(0.25, 0.8, length.out = n),
    treatment = rep(c("control", "warming"), length.out = n)
  )
  rownames(env) <- paste0("trait_site_", seq_len(n))
  spatial <- data.frame(x = seq(0, 1, length.out = n), y = rep(c(0, 0.5), length.out = n))
  rownames(spatial) <- rownames(env)
  lin <- scale(env$pH)[, 1] + 0.3 * spatial$x - 0.2 * as.integer(env$treatment == "warming")
  Y <- sapply(seq_len(4), function(j) as.integer(lin + sin(seq_len(n) / (j + 1)) > median(lin) + (j - 2.5) * 0.1))
  Y[1, ] <- 0L
  Y[n, ] <- 1L
  Y <- as.data.frame(Y)
  names(Y) <- paste0("sp_", seq_len(4))
  rownames(Y) <- rownames(env)
  traits <- data.frame(
    body_size = c(1.1, 1.6, 2.0, 2.5),
    dispersal = c("low", "medium", "medium", "high"),
    trophic = c("herb", "pred", "herb", "pred")
  )
  rownames(traits) <- names(Y)
  groups <- data.frame(species = names(Y), group = c("forb", "grass", "forb", "shrub"))
  folds <- data.frame(site = rownames(Y), fold = rep(1:2, length.out = n))
  cfg <- base_cfg("binomial_spatial_traits_assembly_se")
  cfg$model$spatial_model <- "linear"
  cfg$model$se <- TRUE
  cfg$model$iter <- 5L
  cfg$model$step_size <- 4L
  cfg$spatial_anova_metacommunity$include_space_in_anova <- TRUE
  cfg$spatial_anova_metacommunity$do_internal <- TRUE
  cfg$spatial_anova_metacommunity$do_assembly <- TRUE
  cfg$spatial_anova_metacommunity$assembly_predictor <- "body_size"
  list(id = "binomial_spatial_traits_assembly_se", Y = Y, env = env, spatial = spatial,
       traits = traits, groups = groups, folds = folds, newdata = env[c(2, 5, 8, 12), , drop = FALSE],
       new_spatial = spatial[c(2, 5, 8, 12), , drop = FALSE], cfg = cfg)
}

make_cv_tuning_case <- function() {
  n <- 20L
  env <- data.frame(
    temperature = seq(7, 16, length.out = n),
    canopy = seq(0.1, 0.7, length.out = n)
  )
  rownames(env) <- paste0("cv_site_", seq_len(n))
  Y <- data.frame(
    sp_1 = rep(c(0L, 1L), length.out = n),
    sp_2 = rep(c(1L, 0L, 1L, 0L), length.out = n),
    sp_3 = rep(c(0L, 0L, 1L, 1L), length.out = n)
  )
  names(Y) <- paste0("sp_", seq_len(3))
  rownames(Y) <- rownames(env)
  cfg <- base_cfg("binomial_cv_tuning")
  cfg$model$iter <- 3L
  cfg$model$step_size <- 4L
  cfg$model$sampling <- 100L
  cfg$regularization_biotic$tune_regularization <- TRUE
  cfg$regularization_biotic$tune_steps <- 1L
  cfg$regularization_biotic$cv_k <- 2L
  cfg$spatial_anova_metacommunity$do_anova <- FALSE
  cfg$outputs$importance <- FALSE
  folds <- data.frame(site = rownames(Y), fold = rep(1:2, length.out = n))
  list(id = "binomial_cv_tuning", Y = Y, env = env, folds = folds, cfg = cfg)
}

run_case <- function(case) {
  outdir <- e$make_engine_run_dir("sjSDM", paste0("JSDMStudio_SUITE_", case$id))
  write_csv_pair(case$Y, outdir, "Y.csv")
  write_csv_pair(case$env, outdir, "env.csv")
  write_csv_pair(case$spatial %||% NULL, outdir, "spatial.csv")
  write_csv_pair(case$traits %||% NULL, outdir, "traits.csv")
  write_csv_pair(case$groups %||% NULL, outdir, "species_groups.csv")
  write_csv_pair(case$folds %||% NULL, outdir, "folds.csv")
  write_csv_pair(case$newdata %||% NULL, outdir, "newdata.csv")
  write_csv_pair(case$new_spatial %||% NULL, outdir, "new_spatial.csv")
  yaml::write_yaml(case$cfg, file.path(outdir, "used_config.yml"))
  e$write_sjsdm_reproducible_script(outdir)
  script_file <- file.path(outdir, "reproducible_script", "run_this_sjSDM_analysis.R")
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  log <- system2(rscript, shQuote(script_file), stdout = TRUE, stderr = TRUE)
  writeLines(as.character(log), file.path(outdir, "diagnostics", "suite_stdout_stderr.txt"))
  exit_status <- attr(log, "status") %||% 0
  status_path <- file.path(outdir, "diagnostics", "engine_status.json")
  status <- if (file.exists(status_path)) {
    tryCatch(jsonlite::fromJSON(status_path), error = function(e) list(status = "status_parse_failed", errors = conditionMessage(e)))
  } else {
    list(status = "missing_status", errors = "diagnostics/engine_status.json was not created")
  }
  zipfile <- tryCatch(e$make_zip(outdir), error = function(e) NA_character_)
  required <- c(
    "models/sjSDM_model.rds",
    "tables/coef_environment.csv",
    "tables/covariance_matrix.csv",
    "tables/correlation_matrix.csv",
    "predictions/predictions.csv",
    "tables/residuals.csv",
    "standard/run_summary.csv",
    "standard/effects_long.csv",
    "standard/predictions_long.csv",
    "standard/associations_long.csv"
  )
  if (isTRUE(case$cfg$model$se)) required <- c(required, "tables/standard_errors.csv", "tables/p_values.csv")
  if (!is.null(case$traits %||% NULL)) required <- c(required, "tables/traits_metadata.csv")
  if (!is.null(case$groups %||% NULL)) required <- c(required, "tables/species_groups.csv")
  if (!is.null(case$folds %||% NULL)) required <- c(required, "tables/folds_uploaded.csv")
  if (isTRUE(case$cfg$regularization_biotic$tune_regularization)) required <- c(required, "tables/sjSDM_cv_result_long.csv")
  if (isTRUE(case$cfg$spatial_anova_metacommunity$do_assembly)) required <- c(required, "internal_structure/assembly_effects.csv")
  present <- file.exists(file.path(outdir, required))
  empty_dirs <- GetEmpty <- vapply(list.dirs(outdir, recursive = TRUE, full.names = TRUE), function(d) {
    if (normalizePath(d, winslash = "/", mustWork = FALSE) == normalizePath(outdir, winslash = "/", mustWork = FALSE)) return(NA_character_)
    if (length(list.files(d, all.files = FALSE, no.. = TRUE)) == 0) substring(normalizePath(d, winslash = "/", mustWork = FALSE), nchar(normalizePath(outdir, winslash = "/", mustWork = FALSE)) + 2L) else NA_character_
  }, character(1))
  empty_dirs <- empty_dirs[!is.na(empty_dirs)]
  data.frame(
    case_id = case$id,
    family = case$cfg$model$family,
    env_model = case$cfg$model$env_model,
    spatial_model = case$cfg$model$spatial_model,
    generateSpatialEV = isTRUE(case$cfg$spatial_anova_metacommunity$generateSpatialEV),
    exit_status = as.integer(exit_status),
    engine_status = as.character(status$status %||% "unknown"),
    errors = paste(as.character(status$errors %||% character()), collapse = "; "),
    warnings = paste(as.character(status$warnings %||% character()), collapse = "; "),
    missing_required = paste(required[!present], collapse = "; "),
    empty_dirs = paste(empty_dirs, collapse = "; "),
    outdir = outdir,
    zip = zipfile,
    zip_size = if (!is.na(zipfile) && file.exists(zipfile)) file.info(zipfile)$size else NA_real_,
    stringsAsFactors = FALSE
  )
}

cases <- list(
  make_binomial_case(),
  make_poisson_spatial_case(),
  make_gaussian_dnn_case(),
  make_binomial_logit_dnn_case(),
  make_nbinom_intercept_case(),
  make_spatial_ev_case(),
  make_spatial_dnn_case(),
  make_spatial_traits_assembly_se_case(),
  make_cv_tuning_case()
)

summary_rows <- do.call(rbind, lapply(cases, run_case))
summary_file <- file.path(app_dir, "output", paste0("sjSDM_real_example_suite_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
write.csv(summary_rows, summary_file, row.names = FALSE)
write.csv(summary_rows, file.path(app_dir, "examples", "sjSDM", "last_real_example_suite_summary.csv"), row.names = FALSE)

print(summary_rows[, c("case_id", "family", "env_model", "spatial_model", "generateSpatialEV", "exit_status", "engine_status", "missing_required", "empty_dirs")], row.names = FALSE)
cat("Suite summary:", summary_file, "\n")

failed <- subset(summary_rows, exit_status != 0 | engine_status != "fitted" | nzchar(missing_required) | nzchar(empty_dirs))
if (nrow(failed) > 0) {
  writeLines(c("Failed sjSDM suite cases:", failed$case_id))
  stop("One or more sjSDM real example suite cases failed.", call. = FALSE)
}
