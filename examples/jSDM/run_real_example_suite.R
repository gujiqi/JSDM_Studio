#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x

args0 <- commandArgs(trailingOnly = FALSE)
script_arg <- args0[grepl("^--file=", args0)]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)) else getwd()
root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(root, "app.R"))) root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
setwd(root)

app_env <- new.env(parent = globalenv())
sys.source("app.R", app_env)
sys.source(file.path("R", "jsdm_adapter.R"), app_env)

rscript_bin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

make_site <- function(n, seed = 1) {
  set.seed(seed)
  dat <- data.frame(
    pH = round(rnorm(n, 6.5, 0.45), 3),
    moisture = round(runif(n, 0.15, 0.95), 3),
    canopy = round(runif(n, 0.05, 0.9), 3),
    habitat = sample(c("forest", "grassland", "wetland"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  rownames(dat) <- paste0("site_", seq_len(n))
  dat
}

make_traits <- function(s, seed = 2) {
  set.seed(seed)
  dat <- data.frame(
    life_form = rep(c("herb", "shrub", "tree"), length.out = s),
    height_mm = round(runif(s, 50, 900), 1),
    dispersal = rep(c("wind", "animal"), length.out = s),
    stringsAsFactors = FALSE
  )
  rownames(dat) <- paste0("sp_", seq_len(s))
  dat
}

linear_eta <- function(X, s, seed = 3) {
  set.seed(seed)
  Xnum <- scale(cbind(X$pH, X$moisture, X$canopy))
  B <- matrix(rnorm(ncol(Xnum) * s, 0, 0.75), ncol(Xnum), s)
  habitat_shift <- model.matrix(~ habitat, X)[, -1, drop = FALSE]
  H <- if (ncol(habitat_shift)) habitat_shift %*% matrix(rnorm(ncol(habitat_shift) * s, 0, 0.35), ncol(habitat_shift), s) else 0
  sweep(Xnum %*% B + H, 2, seq(-0.6, 0.6, length.out = s), "+")
}

make_binary <- function(X, s, seed = 4, link = "probit") {
  set.seed(seed)
  eta <- linear_eta(X, s, seed)
  p <- if (identical(link, "logit")) plogis(eta) else pnorm(eta)
  Y <- matrix(rbinom(length(p), 1, p), nrow = nrow(X), ncol = s)
  colnames(Y) <- paste0("sp_", seq_len(s))
  rownames(Y) <- rownames(X)
  for (j in seq_len(ncol(Y))) {
    if (all(Y[, j] == 0)) Y[sample(seq_len(nrow(Y)), 1), j] <- 1
    if (all(Y[, j] == 1)) Y[sample(seq_len(nrow(Y)), 1), j] <- 0
  }
  as.data.frame(Y, check.names = FALSE)
}

make_counts <- function(X, s, seed = 5) {
  set.seed(seed)
  eta <- linear_eta(X, s, seed) / 2
  Y <- matrix(rpois(length(eta), lambda = pmax(0.1, exp(eta))), nrow = nrow(X), ncol = s)
  colnames(Y) <- paste0("sp_", seq_len(s))
  rownames(Y) <- rownames(X)
  as.data.frame(Y, check.names = FALSE)
}

make_gaussian <- function(X, s, seed = 6) {
  set.seed(seed)
  eta <- linear_eta(X, s, seed)
  Y <- eta + matrix(rnorm(length(eta), 0, 0.35), nrow = nrow(X), ncol = s)
  colnames(Y) <- paste0("sp_", seq_len(s))
  rownames(Y) <- rownames(X)
  as.data.frame(round(Y, 3), check.names = FALSE)
}

base_config <- function(case) {
  list(
    project_name = paste0("jSDM_suite_", case$id),
    engine = "jSDM",
    data = list(Y = "Y.csv", XData = "XData.csv", trait_data = "trait_data.csv",
                long_format = "long_format.csv", trials = "trials.csv",
                newdata = "newdata.csv", prediction_ids = "prediction_ids.csv"),
    model = list(
      model_type = case$model_type,
      response_argument = case$response_argument %||% "",
      site_formula = case$site_formula %||% "~ pH + moisture + canopy + habitat",
      trait_formula = case$trait_formula %||% "~ life_form + height_mm + dispersal",
      n_latent = case$n_latent %||% 0,
      site_effect = case$site_effect %||% "none",
      trials = case$trials_scalar %||% 1,
      constrained_latent = identical(case$model_type, "binomial_probit_sp_constrained"),
      constrained_nchains = case$constrained_nchains %||% 2,
      long_site_col = "site",
      long_species_col = "species",
      long_response_col = "presence",
      scale_site_data = case$scale_site_data %||% TRUE,
      include_intercept = case$include_intercept %||% TRUE,
      allow_traits = !is.null(case$Tr)
    ),
    prediction = list(
      do_predict = case$do_predict %||% TRUE,
      predict_type = case$predict_type %||% "mean",
      predict_probs = case$predict_probs %||% "0.1,0.9",
      max_prediction_sites = case$max_prediction_sites %||% 6,
      Id_sites = case$Id_sites %||% "all",
      Id_species = case$Id_species %||% "all",
      prediction_histograms = TRUE
    ),
    diagnostics = list(
      cor_prob = 0.95,
      cor_type = case$cor_type %||% "mean",
      plot_residual_cor = case$plot_residual_cor %||% TRUE,
      plot_associations = case$plot_associations %||% TRUE,
      diag_beta = TRUE, diag_lambda = TRUE, diag_W = TRUE, diag_alpha = TRUE,
      diag_Valpha = TRUE, diag_V = TRUE, diag_deviance = TRUE, coda_summary = TRUE
    ),
    mcmc = list(burnin = 50, mcmc = 50, thin = 1, seed = case$seed %||% 1234,
                verbose = 0, preset = "Automated quick test", ropt = 0.44),
    starts = list(beta_start = 0, gamma_start = 0, lambda_start = 0, W_start = 0,
                  alpha_start = 0, V_alpha = 1, V_start = 1),
    priors = list(shape_Valpha = 0.5, rate_Valpha = 0.0005, shape_V = 0.5, rate_V = 0.0005,
                  mu_beta = 0, V_beta = 10, mu_gamma = 0, V_gamma = 10, mu_lambda = 0, V_lambda = 10),
    outputs = list(real_fit = TRUE, residual_cor = case$residual_cor %||% TRUE,
                   enviro_cor = case$enviro_cor %||% TRUE, traceplots = TRUE,
                   predictions = TRUE, save_model = case$save_model %||% TRUE,
                   report = TRUE, save_mcmc_rds = TRUE, csv_tables = TRUE,
                   figures = case$figures %||% TRUE, table_model_spec = TRUE,
                   table_beta = TRUE, table_lambda = TRUE, table_gamma = TRUE,
                   table_alpha = TRUE, table_predictions = TRUE,
                   copy_inputs = TRUE, save_config = TRUE, zip = TRUE)
  )
}

write_pair <- function(outdir, file, dat, row.names = TRUE) {
  if (is.null(dat)) return(invisible(NULL))
  write.csv(dat, file.path(outdir, "data", file), row.names = row.names)
  write.csv(dat, file.path(outdir, "inputs", file), row.names = row.names)
}

write_case_files <- function(outdir, case) {
  write_pair(outdir, "Y.csv", case$Y, row.names = TRUE)
  write_pair(outdir, "XData.csv", case$X, row.names = TRUE)
  write_pair(outdir, "trait_data.csv", case$Tr, row.names = TRUE)
  write_pair(outdir, "long_format.csv", case$long, row.names = FALSE)
  write_pair(outdir, "trials.csv", case$trials, row.names = FALSE)
  write_pair(outdir, "newdata.csv", case$newdata, row.names = TRUE)
  write_pair(outdir, "prediction_ids.csv", case$prediction_ids, row.names = FALSE)
}

run_case <- function(case) {
  outdir <- app_env$make_engine_run_dir("jSDM", paste0("JSDMStudio_JSDM_SUITE_", case$id))
  write_case_files(outdir, case)
  cfg <- base_config(case)
  yaml::write_yaml(cfg, file.path(outdir, "used_config.yml"))
  check <- app_env$validate_jsdm_full(
    case$Y, case$X, case$Tr, case$long, case$trials, case$newdata, case$prediction_ids,
    cfg$model$model_type, cfg$model$site_formula, cfg$model$trait_formula,
    cfg$model$n_latent, cfg$model$site_effect, cfg$mcmc$burnin, cfg$mcmc$mcmc, cfg$mcmc$thin,
    cfg$model$trials, cfg$model$allow_traits, cfg$prediction$do_predict,
    cfg$model$long_site_col, cfg$model$long_species_col, cfg$model$long_response_col,
    cfg$model$constrained_nchains
  )
  app_env$write_data_check_messages(outdir, check$messages)
  write.csv(app_env$data_summary(case$Y, case$X, case$Tr), file.path(outdir, "tables", "data_summary.csv"), row.names = FALSE)
  if (!isTRUE(check$ok)) {
    writeLines(check$messages, file.path(outdir, "diagnostics", "jSDM_check_failed.txt"))
    stop(sprintf("Case %s failed validation: %s", case$id, paste(check$messages, collapse = "; ")), call. = FALSE)
  }
  app_env$write_jsdm_reproducible_script(outdir)
  log_file <- file.path(outdir, "diagnostics", "suite_stdout_stderr.txt")
  t0 <- Sys.time()
  exit_code <- system2(rscript_bin, shQuote(file.path(outdir, "reproducible_script", "run_this_jSDM_analysis.R")),
                       stdout = log_file, stderr = log_file)
  runtime <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
  status_file <- file.path(outdir, "diagnostics", "engine_status.json")
  status <- if (file.exists(status_file) && requireNamespace("jsonlite", quietly = TRUE)) {
    tryCatch(as.character(jsonlite::fromJSON(status_file)$status), error = function(e) NA_character_)
  } else NA_character_
  zip_file <- tryCatch(app_env$make_zip(outdir), error = function(e) NA_character_)
  required <- c(
    "diagnostics/engine_status.json", "diagnostics/session_info.txt",
    "standard/run_summary.csv", "standard/effects_long.csv", "standard/predictions_long.csv",
    "standard/associations_long.csv", "standard/fit_metrics.csv",
    "results/README_jSDM_results.txt", "results/jSDM_run_summary.csv", "results/fit_metrics.csv",
    "tables/model_spec.csv", "tables/species_parameter_summary.csv", "tables/beta_summary.csv",
    "models/jsdm_model_primary.rds", "report/jSDM_report.html",
    "reproducible_script/run_this_jSDM_analysis.R", "workflow_scripts/run_jSDM_workflow.R"
  )
  if (isTRUE(cfg$outputs$save_model)) required <- c(required, "models/jsdm_model.rds")
  if (isTRUE(cfg$outputs$save_mcmc_rds)) required <- c(required, "mcmc/mcmc_sp.rds")
  missing <- required[!file.exists(file.path(outdir, required))]
  small <- required[file.exists(file.path(outdir, required)) & file.info(file.path(outdir, required))$size <= 0]
  data.frame(
    case_id = case$id,
    model_type = case$model_type,
    exit_code = as.integer(exit_code),
    status = status,
    runtime_seconds = runtime,
    output_dir = outdir,
    zip_file = zip_file,
    missing_required = paste(missing, collapse = ";"),
    empty_required = paste(small, collapse = ";"),
    pass = identical(as.integer(exit_code), 0L) && identical(status, "fitted") && !length(missing) && !length(small) && file.exists(zip_file),
    stringsAsFactors = FALSE
  )
}

X1 <- make_site(16, 10)
X2 <- make_site(14, 20)
X3 <- make_site(15, 30)
X4 <- make_site(14, 40)
X5 <- make_site(10, 50)
X6 <- make_site(15, 60)
X7 <- make_site(12, 70)

Y1 <- make_binary(X1, 4, 11)
Y5 <- make_binary(X5, 3, 55)
long5 <- do.call(rbind, lapply(seq_len(nrow(X5)), function(i) {
  data.frame(site = rownames(X5)[i], species = colnames(Y5),
             presence = as.numeric(Y5[i, ]), x1 = X5$pH[i], x2 = X5$moisture[i],
             stringsAsFactors = FALSE)
}))

cases <- list(
  list(id = "01_binomial_probit_traits_random_predict_ids", model_type = "binomial_probit",
       Y = Y1, X = X1, Tr = make_traits(4, 12),
       n_latent = 2, site_effect = "random", site_formula = "~ pH + moisture + canopy + habitat",
       trait_formula = "~ life_form + height_mm + dispersal",
       newdata = X1[1:5, ], prediction_ids = data.frame(site = rownames(X1)[1:5], species = c(colnames(Y1)[1:3], NA, NA))),
  list(id = "02_binomial_logit_trials_fixed", model_type = "binomial_logit",
       Y = {
         set.seed(22); trials <- rep(4, nrow(X2)); p <- plogis(linear_eta(X2, 3, 22) / 2)
         YY <- matrix(rbinom(length(p), size = rep(trials, 3), prob = p), nrow = nrow(X2), ncol = 3)
         colnames(YY) <- paste0("sp_", 1:3); rownames(YY) <- rownames(X2); as.data.frame(YY, check.names = FALSE)
       },
       X = X2, trials = data.frame(trials = rep(4, nrow(X2))), trials_scalar = 4,
       n_latent = 0, site_effect = "fixed", site_formula = "~ pH + moisture + habitat", predict_type = "mean"),
  list(id = "03_poisson_log_latent_random", model_type = "poisson_log",
       Y = make_counts(X3, 4, 33), X = X3, n_latent = 2, site_effect = "random",
       site_formula = "~ pH + moisture + canopy + habitat", predict_type = "mean"),
  list(id = "04_gaussian_traits_fixed", model_type = "gaussian",
       Y = make_gaussian(X4, 3, 44), X = X4, Tr = make_traits(3, 45), n_latent = 1,
       site_effect = "fixed", site_formula = "~ pH + moisture + canopy", trait_formula = "~ height_mm + life_form"),
  list(id = "05_binomial_probit_long_format_auto_species_formula", model_type = "binomial_probit_long_format",
       long = long5, X = NULL, Y = NULL, n_latent = 0, site_effect = "none",
       site_formula = "~ x1 + x2", do_predict = FALSE, residual_cor = FALSE, enviro_cor = FALSE),
  list(id = "06_binomial_probit_sp_constrained", model_type = "binomial_probit_sp_constrained",
       Y = make_binary(X6, 5, 66), X = X6, n_latent = 2, site_effect = "none",
       site_formula = "~ pH + moisture + canopy", constrained_nchains = 2, do_predict = FALSE),
  list(id = "07_output_switches_probit_none", model_type = "binomial_probit",
       Y = make_binary(X7, 3, 77), X = X7, n_latent = 0, site_effect = "none",
       site_formula = "~ pH + habitat", do_predict = FALSE, save_model = FALSE,
       residual_cor = FALSE, enviro_cor = FALSE, figures = FALSE)
)

results <- do.call(rbind, lapply(cases, run_case))
summary_file <- file.path("output", paste0("jSDM_real_example_suite_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
write.csv(results, summary_file, row.names = FALSE)
print(results[, c("case_id", "model_type", "exit_code", "status", "pass", "missing_required", "empty_required")])
cat("SUMMARY_FILE=", normalizePath(summary_file, winslash = "/", mustWork = FALSE), "\n", sep = "")

if (!all(results$pass)) {
  failed <- results[!results$pass, , drop = FALSE]
  stop(sprintf("jSDM suite failed for %d case(s). See summary: %s", nrow(failed), normalizePath(summary_file, winslash = "/", mustWork = FALSE)), call. = FALSE)
}
