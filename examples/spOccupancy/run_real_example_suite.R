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
sys.source(file.path("R", "spoccupancy_adapter.R"), app_env)
suppressPackageStartupMessages(library(spOccupancy))
rscript_bin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

det_frame <- function(mat, prefix = "obs") {
  out <- as.data.frame(mat, check.names = FALSE)
  names(out) <- paste0(prefix, "_rep", seq_len(ncol(out)))
  out
}

multi_y_wide <- function(arr) {
  sp <- dimnames(arr)[[1]] %||% paste0("sp", seq_len(dim(arr)[1]))
  R <- dim(arr)[3]
  out <- data.frame(row.names = dimnames(arr)[[2]] %||% paste0("site", seq_len(dim(arr)[2])))
  for (s in seq_along(sp)) for (r in seq_len(R)) out[[paste0(sp[s], "_rep", r)]] <- arr[s, , r]
  out
}

species_table <- function(arr) {
  data.frame(species = dimnames(arr)[[1]] %||% paste0("sp", seq_len(dim(arr)[1])), stringsAsFactors = FALSE)
}

base_config <- function(case) {
  list(
    project_name = paste0("spOccupancy_suite_", case$id),
    engine = "spOccupancy",
    data = list(
      y = "y.csv", occ.covs = "occ.covs.csv", det.covs = "det.covs.csv",
      coords = "coords.csv", species = "species.csv",
      integrated_sources = "integrated_sources.csv",
      newdata = "newdata.csv", newcoords = "newcoords.csv", folds = "folds.csv"
    ),
    model = list(
      model_type = case$model_type,
      occ.formula = case$occ_formula %||% "~ x1",
      det.formula = case$det_formula %||% "~ obs",
      formula = case$formula %||% "~ x1",
      data_structure = case$data_structure %||% "single-species replicated",
      range.ind = case$range_ind %||% FALSE
    ),
    spatial_latent_svc = list(
      cov.model = case$cov_model %||% "exponential",
      NNGP = case$NNGP %||% TRUE,
      n.neighbors = case$n_neighbors %||% 3,
      search.type = case$search_type %||% "cb",
      n.factors = case$n_factors %||% 2,
      svc.cols = case$svc_cols %||% "",
      ar1 = case$ar1 %||% FALSE,
      x.positive = case$x_positive %||% FALSE
    ),
    mcmc = list(
      n.batch = case$n_batch %||% 5,
      batch.length = case$batch_length %||% 5,
      n.burn = case$n_burn %||% 5,
      n.thin = case$n_thin %||% 1,
      n.chains = case$n_chains %||% 1,
      accept.rate = case$accept_rate %||% 0.43,
      n.report = case$n_report %||% 5,
      n.omp.threads = case$n_omp_threads %||% 1,
      verbose = case$verbose %||% FALSE,
      seed = case$seed %||% 1234,
      updateMCMC = case$updateMCMC %||% FALSE
    ),
    priors_inits_tuning = list(
      beta.normal = "mean=0,var=2.72",
      alpha.normal = "mean=0,var=2.72",
      community_priors = "beta.comm.normal mean=0,var=2.72; alpha.comm.normal mean=0,var=2.72",
      sigma.sq.ig = "a=2,b=1",
      phi.unif = "a=maxdist,b=maxdist/0.1",
      nu.unif = "0.5,2.5",
      inits = "auto",
      tuning = "phi=0.5",
      fix = TRUE
    ),
    validation_prediction_outputs = list(
      real_fit = TRUE,
      ppcOcc = case$ppcOcc %||% FALSE,
      waicOcc = case$waicOcc %||% TRUE,
      k.fold = case$k_fold %||% 0,
      k.fold.threads = 1,
      k.fold.seed = 100,
      k.fold.only = FALSE,
      predict = case$predict %||% TRUE,
      fitted = case$fitted %||% TRUE,
      save_model = TRUE,
      save_samples = TRUE,
      save_plots = case$save_plots %||% FALSE,
      zip = TRUE,
      copy_inputs = TRUE,
      save_config = TRUE,
      report = TRUE
    )
  )
}

write_pair <- function(outdir, file, dat, row.names = TRUE) {
  if (is.null(dat)) return(invisible(NULL))
  write.csv(dat, file.path(outdir, "data", file), row.names = row.names)
  write.csv(dat, file.path(outdir, "inputs", file), row.names = row.names)
}

write_case_files <- function(outdir, case) {
  write_pair(outdir, "y.csv", case$Y, TRUE)
  write_pair(outdir, "occ.covs.csv", case$occ, TRUE)
  write_pair(outdir, "det.covs.csv", case$det, TRUE)
  write_pair(outdir, "coords.csv", case$coords, TRUE)
  write_pair(outdir, "species.csv", case$species, FALSE)
  write_pair(outdir, "integrated_sources.csv", case$integrated, FALSE)
  write_pair(outdir, "newdata.csv", case$newdata, TRUE)
  write_pair(outdir, "newcoords.csv", case$newcoords, TRUE)
  write_pair(outdir, "folds.csv", case$folds, FALSE)
}

run_case <- function(case) {
  outdir <- app_env$make_engine_run_dir("spOccupancy", paste0("JSDMStudio_spOccupancy_SUITE_", case$id))
  write_case_files(outdir, case)
  cfg <- base_config(case)
  yaml::write_yaml(cfg, file.path(outdir, "used_config.yml"))
  check <- app_env$validate_spoccupancy_full(
    case$Y, case$occ, case$det, case$coords, case$species, case$integrated,
    case$newdata, case$newcoords, case$folds,
    cfg$model$model_type, cfg$mcmc$n.batch, cfg$mcmc$batch.length,
    cfg$mcmc$n.burn, cfg$mcmc$n.thin, cfg$mcmc$n.chains,
    cfg$spatial_latent_svc$n.factors, cfg$spatial_latent_svc$NNGP,
    cfg$spatial_latent_svc$n.neighbors, cfg$validation_prediction_outputs$k.fold,
    cfg$spatial_latent_svc$svc.cols, cfg$model$occ.formula,
    cfg$model$det.formula, cfg$model$formula, cfg$model$data_structure
  )
  app_env$write_data_check_messages(outdir, check$messages)
  if (!isTRUE(check$ok)) {
    writeLines(check$messages, file.path(outdir, "diagnostics", "spOccupancy_check_failed.txt"))
    stop(sprintf("Case %s failed validation: %s", case$id, paste(check$messages, collapse = "; ")), call. = FALSE)
  }
  app_env$write_spoccupancy_reproducible_script(outdir)
  log_file <- file.path(outdir, "diagnostics", "suite_stdout_stderr.txt")
  t0 <- Sys.time()
  exit_code <- system2(rscript_bin, shQuote(file.path(outdir, "reproducible_script", "run_this_spOccupancy_analysis.R")),
                       stdout = log_file, stderr = log_file)
  runtime <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
  status_file <- file.path(outdir, "diagnostics", "engine_status.json")
  status <- if (file.exists(status_file) && requireNamespace("jsonlite", quietly = TRUE)) {
    tryCatch(as.character(jsonlite::fromJSON(status_file)$status), error = function(e) NA_character_)
  } else NA_character_
  zip_file <- tryCatch(app_env$make_zip(outdir), error = function(e) NA_character_)
  required <- c(
    "diagnostics/engine_status.json", "diagnostics/session_info.txt",
    "standard/run_summary.csv", "standard/effects_long.csv",
    "standard/predictions_long.csv", "standard/associations_long.csv",
    "standard/fit_metrics.csv", "results/README_spOccupancy_results.txt",
    "results/run_summary.csv", "results/fit_metrics.csv",
    "tables/model_settings_used.csv", "tables/summary_beta.csv", "tables/summary_alpha.csv",
    "models/spOccupancy_model.rds", "models/spOccupancy_call_args.rds",
    "samples/posterior_sample_manifest.csv", "report/spOccupancy_report.html",
    "reproducible_script/run_this_spOccupancy_analysis.R",
    "workflow_scripts/run_spOccupancy_workflow.R"
  )
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
    pass = identical(as.integer(exit_code), 0L) && identical(status, "fitted") &&
      !length(missing) && !length(small) && file.exists(zip_file),
    stringsAsFactors = FALSE
  )
}

make_occ <- function(sim) {
  data.frame(x1 = sim$X[, 2], habitat = factor(ifelse(sim$X[, 2] > median(sim$X[, 2]), "high", "low")))
}

set.seed(11)
J <- 12
sim1 <- simOcc(4, 3, rep(3, J), 3, beta = c(0, 0.5), alpha = c(0, 0.2), sp = FALSE)
sim2 <- simOcc(4, 3, rep(3, J), 3, beta = c(0, 0.5), alpha = c(0, 0.2), sp = TRUE,
               cov.model = "exponential", sigma.sq = 1, phi = 3, nu = 0.5)
sim3 <- simMsOcc(3, 3, rep(3, 9), 3, N = 3,
                 beta = matrix(rnorm(3 * 2, 0, 0.4), 3, 2),
                 alpha = matrix(rnorm(3 * 2, 0, 0.4), 3, 2), sp = FALSE)
sim4 <- simMsOcc(3, 3, rep(3, 9), 3, N = 3,
                 beta = matrix(rnorm(3 * 2, 0, 0.4), 3, 2),
                 alpha = matrix(rnorm(3 * 2, 0, 0.4), 3, 2), sp = TRUE,
                 cov.model = "exponential", sigma.sq = rep(1, 3), phi = rep(3, 3), nu = rep(0.5, 3))
sim5 <- simMsOcc(3, 3, rep(3, 9), 3, N = 4,
                 beta = matrix(rnorm(4 * 2, 0, 0.4), 4, 2),
                 alpha = matrix(rnorm(4 * 2, 0, 0.4), 4, 2), sp = FALSE,
                 factor.model = TRUE, n.factors = 2)
sim6 <- simIntOcc(n.data = 2, J.x = 3, J.y = 3, J.obs = c(9, 9),
                  n.rep = list(rep(2, 9), rep(2, 9)), n.rep.max = c(2, 2),
                  beta = c(0, 0.5), alpha = list(c(0, 0.2), c(0, -0.1)), sp = FALSE)
sim7 <- simOcc(4, 3, rep(3, J), 3, beta = c(0, 0.5), alpha = c(0, 0.2), sp = TRUE,
               svc.cols = 1, cov.model = "exponential", sigma.sq = 1, phi = 3, nu = 0.5)

int_y <- data.frame(src1_rep1 = sim6$y[[1]][, 1], src1_rep2 = sim6$y[[1]][, 2],
                    src2_rep1 = sim6$y[[2]][, 1], src2_rep2 = sim6$y[[2]][, 2])
int_det <- data.frame(src1_obs_rep1 = sim6$X.p[[1]][, 1, 2], src1_obs_rep2 = sim6$X.p[[1]][, 2, 2],
                      src2_obs_rep1 = sim6$X.p[[2]][, 1, 2], src2_obs_rep2 = sim6$X.p[[2]][, 2, 2])

cases <- list(
  list(id = "01_PGOcc_single_predict_ppc_waic", model_type = "PGOcc",
       Y = as.data.frame(sim1$y), occ = make_occ(sim1), det = det_frame(sim1$X.p[, , 2]),
       newdata = make_occ(sim1)[1:4, , drop = FALSE], ppcOcc = TRUE),
  list(id = "02_spPGOcc_spatial_NNGP", model_type = "spPGOcc",
       Y = as.data.frame(sim2$y), occ = make_occ(sim2), det = det_frame(sim2$X.p[, , 2]),
       coords = as.data.frame(sim2$coords), newdata = make_occ(sim2)[1:4, , drop = FALSE],
       newcoords = as.data.frame(sim2$coords[1:4, , drop = FALSE])),
  list(id = "03_msPGOcc_multispecies", model_type = "msPGOcc",
       Y = multi_y_wide(sim3$y), occ = make_occ(sim3), det = det_frame(sim3$X.p[, , 2]),
       species = species_table(sim3$y), data_structure = "multi-species replicated"),
  list(id = "04_spMsPGOcc_spatial_multispecies", model_type = "spMsPGOcc",
       Y = multi_y_wide(sim4$y), occ = make_occ(sim4), det = det_frame(sim4$X.p[, , 2]),
       coords = as.data.frame(sim4$coords), species = species_table(sim4$y),
       data_structure = "multi-species replicated"),
  list(id = "05_lfMsPGOcc_latent_factor", model_type = "lfMsPGOcc",
       Y = multi_y_wide(sim5$y), occ = make_occ(sim5), det = det_frame(sim5$X.p[, , 2]),
       coords = as.data.frame(sim5$coords), species = species_table(sim5$y),
       newdata = make_occ(sim5)[1:3, , drop = FALSE],
       newcoords = as.data.frame(sim5$coords[1:3, , drop = FALSE] + 0.01),
       data_structure = "multi-species replicated", n_factors = 2),
  list(id = "06_intPGOcc_integrated_sources", model_type = "intPGOcc",
       Y = int_y, occ = data.frame(x1 = sim6$X.obs[, 2]), det = int_det,
       integrated = data.frame(source = c("source1", "source2"), prefix = c("src1", "src2")),
       data_structure = "integrated multiple sources"),
  list(id = "07_svcPGOcc_spatial_varying_coefficients", model_type = "svcPGOcc",
       Y = as.data.frame(sim7$y), occ = make_occ(sim7), det = det_frame(sim7$X.p[, , 2]),
       coords = as.data.frame(sim7$coords), svc_cols = "1",
       data_structure = "SVC / spatially varying coefficients")
)

results <- do.call(rbind, lapply(cases, run_case))
summary_file <- file.path("output", paste0("spOccupancy_real_example_suite_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
write.csv(results, summary_file, row.names = FALSE)
print(results[, c("case_id", "model_type", "exit_code", "status", "pass", "missing_required", "empty_required")])
cat("SUMMARY_FILE=", normalizePath(summary_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
if (!all(results$pass)) {
  failed <- results[!results$pass, , drop = FALSE]
  stop(sprintf("spOccupancy suite failed for %d case(s). See summary: %s",
               nrow(failed), normalizePath(summary_file, winslash = "/", mustWork = FALSE)), call. = FALSE)
}
