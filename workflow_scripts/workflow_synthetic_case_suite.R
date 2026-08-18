# Synthetic branch-coverage audit suite for JSDM Studio workflows.
# It reuses the Universal Benchmark data generator and the real engine adapters.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_arg <- args[startsWith(args, file_arg)]
script_path <- if (length(script_arg)) {
  normalizePath(sub(file_arg, "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
} else {
  normalizePath("workflow_scripts/workflow_synthetic_case_suite.R", winslash = "/", mustWork = FALSE)
}
app_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

load_runner_env <- function() {
  runner_path <- file.path(app_dir, "workflow_scripts", "universal_benchmark_runner.R")
  txt <- readLines(runner_path, warn = FALSE)
  txt <- txt[!grepl("^\\s*if \\(!interactive\\(\\)\\) main\\(\\)\\s*$", txt)]
  e <- new.env(parent = globalenv())
  eval(parse(text = paste(txt, collapse = "\n")), envir = e)
  e
}

runner_env <- load_runner_env()
app_env <- runner_env$source_app_env()

safe_yaml_write <- function(x, path) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("yaml package is required.", call. = FALSE)
  yaml::write_yaml(x, path)
}

safe_read_status <- function(outdir, engine) {
  tryCatch(runner_env$read_status(outdir, engine), error = function(e) {
    list(engine = engine, status = "missing_status", warnings = character(), errors = conditionMessage(e))
  })
}

write_pair <- function(x, outdir, filename, row.names = TRUE) {
  runner_env$write_csv_pair(x, outdir, filename, row.names = row.names)
}

write_text_pair <- function(text, outdir, filename) {
  runner_env$write_text_pair(text, outdir, filename)
}

contract_dirs <- c("inputs", "data", "models", "tables", "results", "plots", "predictions",
                   "diagnostics", "workflow_scripts", "reproducible_script", "standard", "report")
contract_files <- c("used_config.yml", "diagnostics/engine_status.json",
                    "diagnostics/data_check_messages.csv", "diagnostics/session_info.txt",
                    "standard/run_summary.csv", "standard/effects_long.csv",
                    "standard/predictions_long.csv", "standard/associations_long.csv",
                    "standard/fit_metrics.csv", "standard/diagnostics_long.csv",
                    "standard/output_manifest.csv",
                    "standard/effects_species_environment.csv",
                    "standard/predictions_site_species.csv",
                    "standard/associations_species_species.csv")

prepare_case <- function(root, engine, case_id, dat, Y = NULL) {
  outdir <- file.path(root, engine, case_id)
  if (dir.exists(outdir)) unlink(outdir, recursive = TRUE, force = TRUE)
  runner_env$ensure_dirs(outdir, engine)
  dat_case <- dat
  if (!is.null(Y)) dat_case$Y <- as.data.frame(Y, check.names = FALSE)
  runner_env$write_engine_inputs(outdir, dat_case, engine)
  list(outdir = outdir, dat = dat_case)
}

mark_case_status <- function(outdir, engine, status, warnings = character(), errors = character(),
                             Y = NULL, X = NULL, runtime_seconds = 0) {
  runner_env$write_status(outdir, engine, status = status, warnings = warnings,
                          errors = errors, runtime_seconds = runtime_seconds)
  app_env$write_standard_outputs(outdir, engine,
                                 list(status = status, warnings = warnings, errors = errors),
                                 Y = Y, X = X)
}

finish_case <- function(engine, case_id, outdir, dat_case, started, notes = "") {
  st <- safe_read_status(outdir, engine)
  if (identical(st$status %||% "missing_status", "missing_status")) {
    mark_case_status(outdir, engine, "fit_failed",
                     errors = "Engine finished without writing diagnostics/engine_status.json.",
                     Y = dat_case$Y, X = dat_case$X,
                     runtime_seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2))
    st <- safe_read_status(outdir, engine)
  }
  if (!file.exists(file.path(outdir, "standard", "run_summary.csv"))) {
    app_env$write_standard_outputs(outdir, engine, st, Y = dat_case$Y, X = dat_case$X)
  }
  runner_env$standardize_engine_synonyms(outdir, engine, dat_case)
  app_env$ensure_output_contract(outdir, engine, st, Y = dat_case$Y, X = dat_case$X)
  runner_env$fill_empty_dirs(outdir, engine)
  zipfile <- tryCatch(runner_env$make_zip_file(outdir), error = function(e) {
    writeLines(conditionMessage(e), file.path(outdir, "diagnostics", "zip_error.txt"))
    NA_character_
  })
  missing_dirs <- contract_dirs[!dir.exists(file.path(outdir, contract_dirs))]
  missing_files <- contract_files[!file.exists(file.path(outdir, contract_files))]
  zip_ok <- !is.na(zipfile) && file.exists(zipfile) && file.info(zipfile)$size > 0
  fail_reason <- paste(unlist(st$errors %||% character()), collapse = "; ")
  warn_text <- paste(unlist(st$warnings %||% character()), collapse = "; ")
  data.frame(
    engine = engine,
    case_id = case_id,
    status = st$status %||% "unknown",
    output_folder = normalizePath(outdir, winslash = "/", mustWork = FALSE),
    zip = if (zip_ok) normalizePath(zipfile, winslash = "/", mustWork = FALSE) else as.character(zipfile),
    zip_size_kb = if (zip_ok) round(file.info(zipfile)$size / 1024, 1) else NA_real_,
    file_count = length(list.files(outdir, recursive = TRUE, all.files = FALSE)),
    missing_count = length(c(missing_dirs, missing_files)),
    missing_items = paste(c(missing_dirs, missing_files), collapse = "; "),
    failure_reason = fail_reason,
    warnings = warn_text,
    notes = notes,
    runtime_seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2),
    stringsAsFactors = FALSE
  )
}

make_long_format <- function(Y, X) {
  Y <- as.data.frame(Y, check.names = FALSE)
  X <- as.data.frame(X, check.names = FALSE)
  sites <- rownames(Y)
  out <- expand.grid(site = sites, species = colnames(Y), KEEP.OUT.ATTRS = FALSE,
                     stringsAsFactors = FALSE)
  ym <- as.matrix(Y)
  out$presence <- as.integer(ym[cbind(match(out$site, sites), match(out$species, colnames(Y)))])
  xdf <- cbind(site = rownames(X), X, stringsAsFactors = FALSE)
  merge(out, xdf, by = "site", all.x = TRUE, sort = FALSE)
}

make_single_detection <- function(dat, species = colnames(dat$Y)[1]) {
  cols <- grep(paste0("^", species, "(_|\\.)rep"), names(dat$det), value = TRUE)
  if (!length(cols)) cols <- names(dat$det)[seq_len(min(3L, ncol(dat$det)))]
  y <- dat$det[, cols, drop = FALSE]
  names(y) <- paste0("rep", seq_len(ncol(y)))
  rownames(y) <- rownames(dat$det)
  y
}

make_gjam_mixed_y <- function(dat) {
  Y <- dat$Y
  Yn <- dat$Y_normal
  Yc <- dat$Y_count
  fc_raw <- Yn[, min(4, ncol(Yn))]
  fc <- (fc_raw - min(fc_raw)) / max(1e-6, diff(range(fc_raw)))
  fc1 <- pmin(0.95, pmax(0.05, 0.1 + 0.8 * fc))
  fc2 <- 1 - fc1
  oc_raw <- Yn[, min(6, ncol(Yn))]
  oc <- as.integer(cut(oc_raw, breaks = 4, labels = FALSE)) - 1L
  out <- data.frame(
    PA = Y[, 1],
    CON = Yn[, min(2, ncol(Yn))],
    CA = Yc[, min(3, ncol(Yc))],
    FC_1 = fc1,
    FC_2 = fc2,
    OC = oc,
    stringsAsFactors = FALSE
  )
  names(out) <- colnames(Y)[seq_len(ncol(out))]
  rownames(out) <- rownames(Y)
  out
}

run_hmsc_case <- function(root, case_id, dat, cfg, Y, notes = "") {
  started <- Sys.time()
  prep <- prepare_case(root, "Hmsc", case_id, dat, Y = Y)
  outdir <- prep$outdir
  dat_case <- prep$dat
  write_pair(dat$traits, outdir, "TrData.csv")
  write_pair(dat$study, outdir, "studyDesign.csv")
  write_pair(dat$coords, outdir, "coordinates.csv")
  write_pair(dat$phylo_cov, outdir, "phylo_cov.csv")
  write_text_pair(dat$newick, outdir, "phylogeny.nwk")
  safe_yaml_write(cfg, file.path(outdir, "used_config.yml"))
  check <- app_env$validate_hmsc(dat_case$Y, dat_case$X,
                                 if (isTRUE(cfg$model$use_traits)) dat$traits else NULL,
                                 dat$study, dat$coords,
                                 distr = cfg$model$distr,
                                 XFormula = cfg$model$XFormula,
                                 TrFormula = cfg$model$TrFormula,
                                 use_traits = isTRUE(cfg$model$use_traits),
                                 use_phylogeny = isTRUE(cfg$model$use_phylogeny),
                                 random_mode = cfg$model$random_mode,
                                 random_effect_column = cfg$model$random_effect_column,
                                 spatial_method = cfg$model$spatial_method,
                                 nNeighbours = cfg$model$nNeighbours,
                                 lon_col = cfg$model$lon_col,
                                 lat_col = cfg$model$lat_col,
                                 samples = cfg$mcmc$samples,
                                 transient = cfg$mcmc$transient,
                                 thin = cfg$mcmc$thin,
                                 nChains = cfg$mcmc$nChains,
                                 nParallel = cfg$mcmc$nParallel,
                                 nfMin = cfg$model$nfMin,
                                 nfMax = cfg$model$nfMax,
                                 nfolds = cfg$outputs$predictions$nfolds,
                                 partition_column = cfg$outputs$predictions$partition_column)
  app_env$write_data_check_messages(outdir, check$messages)
  if (!isTRUE(check$ok)) {
    mark_case_status(outdir, "Hmsc", "check_failed", warnings = check$messages,
                     errors = "Hmsc synthetic case check failed.", Y = dat_case$Y, X = dat_case$X)
  } else {
    tryCatch({
      app_env$run_hmsc_s1s7_pipeline(outdir, cfg, dat_case$Y, dat_case$X,
                                     if (isTRUE(cfg$model$use_traits)) dat$traits else NULL,
                                     dat$study, dat$coords, log_fun = message)
    }, error = function(e) {
      mark_case_status(outdir, "Hmsc", "fit_failed", warnings = check$messages,
                       errors = conditionMessage(e), Y = dat_case$Y, X = dat_case$X)
    })
  }
  finish_case("Hmsc", case_id, outdir, dat_case, started, notes)
}

run_hmschpc_case <- function(root, case_id, dat, cfg, Y, notes = "") {
  started <- Sys.time()
  prep <- prepare_case(root, "Hmsc-HPC", case_id, dat, Y = Y)
  outdir <- prep$outdir
  dat_case <- prep$dat
  safe_yaml_write(cfg, file.path(outdir, "used_config.yml"))
  check <- app_env$validate_hmschpc_full(
    Y = dat_case$Y, X = dat_case$X,
    Tr = if (isTRUE(cfg$model$use_traits)) dat$traits else NULL,
    study = dat$study, coord = dat$coords,
    phylo_cov = if (identical(cfg$model$phylogeny_mode, "covariance")) dat$phylo_cov else NULL,
    phylo_tree_file = if (identical(cfg$model$phylogeny_mode, "newick")) file.path(outdir, "data", "phylogeny.nwk") else NULL,
    newdata = dat$newdata,
    distribution = cfg$model$distribution,
    XFormula = cfg$model$XFormula,
    use_traits = isTRUE(cfg$model$use_traits),
    trait_formula = cfg$model$trait_formula,
    phylogeny_mode = cfg$model$phylogeny_mode,
    random_mode = cfg$random_effects$mode,
    random_column = cfg$random_effects$column,
    coord_x = cfg$random_effects$coord_x,
    coord_y = cfg$random_effects$coord_y,
    nf = cfg$random_effects$nf,
    nfMin = cfg$random_effects$nfMin,
    nfMax = cfg$random_effects$nfMax,
    samples = cfg$sampler$samples,
    transient = cfg$sampler$transient,
    thin = cfg$sampler$thin,
    chains = cfg$sampler$chains,
    verbose = cfg$sampler$verbose,
    run_sampler = isTRUE(cfg$sampler$run_sampler),
    random_slope_formula = cfg$random_effects$x_formula %||% "",
    random_name = cfg$random_effects$name,
    alpha = cfg$random_effects$alpha,
    chains_to_run = cfg$sampler$chains_to_run %||% integer(),
    seed = cfg$sampler$seed,
    precision = cfg$sampler$precision,
    truncated_normal_library = cfg$sampler$truncated_normal_library,
    hmcleapfrog = cfg$sampler$hmcleapfrog,
    hmcthin = cfg$sampler$hmcthin,
    python = cfg$runtime$python,
    python_source = cfg$runtime$python_source,
    predictions = cfg$outputs$predictions,
    diagnostics = cfg$outputs$diagnostics,
    plots = cfg$outputs$plots)
  app_env$write_data_check_messages(outdir, check$messages)
  if (!isTRUE(check$ok)) {
    mark_case_status(outdir, "Hmsc-HPC", "check_failed", warnings = check$messages,
                     errors = "Hmsc-HPC synthetic case check failed.", Y = dat_case$Y, X = dat_case$X)
  } else {
    tryCatch({
      app_env$run_hmschpc_workflow(outdir, cfg, dat_case$Y, dat_case$X,
                                   Tr = if (isTRUE(cfg$model$use_traits)) dat$traits else NULL,
                                   study = dat$study, coord = dat$coords,
                                   phylo_cov = if (identical(cfg$model$phylogeny_mode, "covariance")) dat$phylo_cov else NULL,
                                   phylo_tree_file = if (identical(cfg$model$phylogeny_mode, "newick")) file.path(outdir, "data", "phylogeny.nwk") else NULL,
                                   newdata = dat$newdata, log_fun = message)
    }, error = function(e) {
      mark_case_status(outdir, "Hmsc-HPC", "fit_failed", warnings = check$messages,
                       errors = conditionMessage(e), Y = dat_case$Y, X = dat_case$X)
    })
  }
  finish_case("Hmsc-HPC", case_id, outdir, dat_case, started, notes)
}

run_jsdm_case <- function(root, case_id, dat, cfg, Y = dat$Y, long = NULL, trials = dat$trial,
                          notes = "") {
  started <- Sys.time()
  prep <- prepare_case(root, "jSDM", case_id, dat, Y = Y)
  outdir <- prep$outdir
  dat_case <- prep$dat
  if (!is.null(long)) write_pair(long, outdir, "long_format.csv", row.names = FALSE)
  write_pair(trials, outdir, "trials.csv")
  safe_yaml_write(cfg, file.path(outdir, "used_config.yml"))
  check <- app_env$validate_jsdm_full(
    Y = if (identical(cfg$model$model_type, "binomial_probit_long_format")) NULL else dat_case$Y,
    X = dat_case$X,
    Tr = if (isTRUE(cfg$model$allow_traits)) dat$traits else NULL,
    long = if (identical(cfg$model$model_type, "binomial_probit_long_format")) long else NULL,
    trials = trials,
    newdata = dat$newdata,
    prediction_ids = NULL,
    model_type = cfg$model$model_type,
    site_formula = cfg$model$site_formula,
    trait_formula = cfg$model$trait_formula,
    n_latent = cfg$model$n_latent,
    site_effect = cfg$model$site_effect,
    burnin = cfg$mcmc$burnin,
    mcmc = cfg$mcmc$mcmc,
    thin = cfg$mcmc$thin,
    trials_scalar = cfg$model$trials,
    allow_traits = isTRUE(cfg$model$allow_traits),
    do_predict = isTRUE(cfg$prediction$do_predict),
    long_site_col = cfg$model$long_site_col,
    long_species_col = cfg$model$long_species_col,
    long_response_col = cfg$model$long_response_col,
    constrained_nchains = cfg$model$constrained_nchains)
  app_env$write_data_check_messages(outdir, check$messages)
  if (!isTRUE(check$ok)) {
    mark_case_status(outdir, "jSDM", "check_failed", warnings = check$messages,
                     errors = "jSDM synthetic case check failed.", Y = dat_case$Y, X = dat_case$X)
  } else {
    app_env$write_jsdm_reproducible_script(outdir)
    code <- runner_env$run_rscript(file.path(outdir, "reproducible_script", "run_this_jSDM_analysis.R"), outdir, "jSDM")
    if (!identical(as.integer(code), 0L) && !(safe_read_status(outdir, "jSDM")$status %in% c("fitted", "fit_failed", "check_failed"))) {
      mark_case_status(outdir, "jSDM", "fit_failed", warnings = check$messages,
                       errors = paste0("jSDM Rscript exited with status ", code, "."),
                       Y = dat_case$Y, X = dat_case$X)
    }
  }
  finish_case("jSDM", case_id, outdir, dat_case, started, notes)
}

run_gjam_case <- function(root, case_id, dat, cfg, Y, type_vec, notes = "") {
  started <- Sys.time()
  prep <- prepare_case(root, "GJAM", case_id, dat, Y = Y)
  outdir <- prep$outdir
  dat_case <- prep$dat
  type_tab <- data.frame(response = colnames(Y), typeName = type_vec, stringsAsFactors = FALSE)
  write_pair(type_tab, outdir, "typeNames.csv", row.names = FALSE)
  write_pair(dat$traits, outdir, "specByTrait.csv")
  write_pair(data.frame(trait = names(dat$traits), typeName = c("CAT", "CON", "CON"), stringsAsFactors = FALSE),
             outdir, "traitTypes.csv", row.names = FALSE)
  safe_yaml_write(cfg, file.path(outdir, "used_config.yml"))
  check <- app_env$validate_gjam_full(
    dat_case$Y, dat_case$X, type_table = type_tab,
    type_text = cfg$model$typeNames_text,
    single_type = cfg$model$type_single,
    fcgroups = cfg$response_types$FCgroups,
    ccgroups = cfg$response_types$CCgroups,
    ng = cfg$model$ng,
    burnin = cfg$model$burnin,
    holdoutN = cfg$model$holdoutN,
    random = cfg$model$random,
    notStandard = cfg$model$notStandard,
    formula_text = cfg$model$formula,
    newdata = dat$newdata,
    specByTrait = if (isTRUE(cfg$analysis$do_traits)) dat$traits else NULL,
    traitTypes = data.frame(trait = names(dat$traits), typeName = c("CAT", "CON", "CON")),
    do_traits = isTRUE(cfg$analysis$do_traits),
    trimY = isTRUE(cfg$response_types$trimY),
    trim_minObs = cfg$response_types$trim_minObs,
    REDUCT = isTRUE(cfg$model$REDUCT),
    reduct_N = cfg$model$reductList$N,
    reduct_r = cfg$model$reductList$r)
  app_env$write_data_check_messages(outdir, check$messages)
  if (!isTRUE(check$ok)) {
    mark_case_status(outdir, "GJAM", "check_failed", warnings = check$messages,
                     errors = "GJAM synthetic case check failed.", Y = dat_case$Y, X = dat_case$X)
  } else {
    app_env$write_gjam_reproducible_script(outdir)
    code <- runner_env$run_rscript(file.path(outdir, "reproducible_script", "run_this_GJAM_analysis.R"), outdir, "GJAM")
    if (!identical(as.integer(code), 0L) && !(safe_read_status(outdir, "GJAM")$status %in% c("fitted", "fit_failed", "check_failed"))) {
      mark_case_status(outdir, "GJAM", "fit_failed", warnings = check$messages,
                       errors = paste0("GJAM Rscript exited with status ", code, "."),
                       Y = dat_case$Y, X = dat_case$X)
    }
  }
  finish_case("GJAM", case_id, outdir, dat_case, started, notes)
}

run_spocc_case <- function(root, case_id, dat, cfg, y, species_table = NULL, notes = "") {
  started <- Sys.time()
  prep <- prepare_case(root, "spOccupancy", case_id, dat, Y = dat$Y)
  outdir <- prep$outdir
  dat_case <- prep$dat
  write_pair(y, outdir, "y.csv")
  write_pair(dat$occ_covs, outdir, "occ.covs.csv")
  write_pair(dat$det_cov, outdir, "det.covs.csv")
  write_pair(dat$coords, outdir, "coords.csv")
  if (is.null(species_table)) species_table <- data.frame(species = colnames(dat$Y), stringsAsFactors = FALSE)
  write_pair(species_table, outdir, "species.csv", row.names = FALSE)
  safe_yaml_write(cfg, file.path(outdir, "used_config.yml"))
  check <- app_env$validate_spoccupancy_full(
    y = y, occ = dat$occ_covs, det = dat$det_cov, coords = dat$coords,
    species = species_table, integrated = NULL, newdata = dat$newdata,
    newcoords = dat$newcoords, folds = dat$folds,
    model_type = cfg$model$model_type,
    n.batch = cfg$mcmc$n.batch,
    batch.length = cfg$mcmc$batch.length,
    n.burn = cfg$mcmc$n.burn,
    n.thin = cfg$mcmc$n.thin,
    n.chains = cfg$mcmc$n.chains,
    n.factors = cfg$spatial_latent_svc$n.factors,
    NNGP = isTRUE(cfg$spatial_latent_svc$NNGP),
    n.neighbors = cfg$spatial_latent_svc$n.neighbors,
    k.fold = cfg$validation_prediction_outputs$k.fold,
    svc.cols = cfg$spatial_latent_svc$svc.cols,
    occ.formula = cfg$model$occ.formula,
    det.formula = cfg$model$det.formula,
    formula = cfg$model$formula,
    data_structure = cfg$model$data_structure)
  app_env$write_data_check_messages(outdir, check$messages)
  if (!isTRUE(check$ok)) {
    mark_case_status(outdir, "spOccupancy", "check_failed", warnings = check$messages,
                     errors = "spOccupancy synthetic case check failed.", Y = dat_case$Y, X = dat_case$X)
  } else {
    app_env$write_spoccupancy_reproducible_script(outdir)
    code <- runner_env$run_rscript(file.path(outdir, "reproducible_script", "run_this_spOccupancy_analysis.R"), outdir, "spOccupancy")
    if (!identical(as.integer(code), 0L) && !(safe_read_status(outdir, "spOccupancy")$status %in% c("fitted", "fit_failed", "check_failed"))) {
      mark_case_status(outdir, "spOccupancy", "fit_failed", warnings = check$messages,
                       errors = paste0("spOccupancy Rscript exited with status ", code, "."),
                       Y = dat_case$Y, X = dat_case$X)
    }
  }
  finish_case("spOccupancy", case_id, outdir, dat_case, started, notes)
}

run_sjsdm_case <- function(root, case_id, dat, cfg, Y, notes = "") {
  started <- Sys.time()
  prep <- prepare_case(root, "sjSDM", case_id, dat, Y = Y)
  outdir <- prep$outdir
  dat_case <- prep$dat
  write_pair(dat_case$Y, outdir, "Y.csv")
  write_pair(dat$X, outdir, "env.csv")
  write_pair(dat$coords, outdir, "spatial.csv")
  write_pair(dat$traits, outdir, "traits.csv")
  write_pair(dat$newdata, outdir, "newdata.csv")
  write_pair(dat$coords[seq_len(nrow(dat$newdata)), , drop = FALSE], outdir, "new_spatial.csv")
  safe_yaml_write(cfg, file.path(outdir, "used_config.yml"))
  check <- app_env$validate_sjsdm(dat_case$Y, dat$X, dat$coords, dat$traits,
                                  family = cfg$model$family,
                                  env_formula = cfg$model$env_formula,
                                  spatial_formula = cfg$model$spatial_formula,
                                  env_model = cfg$model$env_model,
                                  spatial_model = cfg$model$spatial_model,
                                  newdata = dat$newdata,
                                  iter = cfg$model$iter,
                                  sampling = cfg$model$sampling,
                                  learning_rate = cfg$dnn_optimizer$learning_rate,
                                  step_size = cfg$model$step_size,
                                  parallel = cfg$model$parallel,
                                  cv_k = cfg$regularization_biotic$cv_k %||% 0L,
                                  dnn_hidden = cfg$dnn_optimizer$hidden,
                                  dropout = cfg$dnn_optimizer$dropout,
                                  device = cfg$dnn_optimizer$device,
                                  dtype = cfg$model$dtype,
                                  verbose = cfg$model$verbose,
                                  optimizer = cfg$dnn_optimizer$optimizer,
                                  weight_decay = cfg$dnn_optimizer$weight_decay,
                                  scheduler = cfg$control$scheduler,
                                  lr_reduce_factor = cfg$control$lr_reduce_factor,
                                  early_stopping_training = cfg$control$early_stopping_training,
                                  mixed = cfg$control$mixed,
                                  biotic_lambda = cfg$regularization_biotic$biotic_lambda,
                                  biotic_alpha = cfg$regularization_biotic$biotic_alpha,
                                  generate_spatial_ev = cfg$spatial_anova_metacommunity$generateSpatialEV,
                                  spatial_ev_threshold = cfg$spatial_anova_metacommunity$spatial_ev_threshold,
                                  anova_samples = cfg$spatial_anova_metacommunity$anova_samples,
                                  tune_steps = cfg$regularization_biotic$tune_steps %||% 0L)
  app_env$write_data_check_messages(outdir, check$messages)
  if (!isTRUE(check$ok)) {
    mark_case_status(outdir, "sjSDM", "check_failed", warnings = check$messages,
                     errors = "sjSDM synthetic case check failed.", Y = dat_case$Y, X = dat_case$X)
  } else {
    app_env$write_sjsdm_reproducible_script(outdir)
    code <- runner_env$run_rscript(file.path(outdir, "reproducible_script", "run_this_sjSDM_analysis.R"), outdir, "sjSDM")
    if (!identical(as.integer(code), 0L) && !(safe_read_status(outdir, "sjSDM")$status %in% c("fitted", "fit_failed", "check_failed"))) {
      mark_case_status(outdir, "sjSDM", "fit_failed", warnings = check$messages,
                       errors = paste0("sjSDM Rscript exited with status ", code, "."),
                       Y = dat_case$Y, X = dat_case$X)
    }
  }
  finish_case("sjSDM", case_id, outdir, dat_case, started, notes)
}

run_boral_case <- function(root, case_id, dat, cfg, Y, notes = "") {
  started <- Sys.time()
  prep <- prepare_case(root, "boral", case_id, dat, Y = Y)
  outdir <- prep$outdir
  dat_case <- prep$dat
  write_pair(dat_case$Y, outdir, "Y.csv")
  write_pair(dat$X, outdir, "XData.csv")
  write_pair(dat$traits, outdir, "traits.csv")
  write_pair(dat$trial, outdir, "trial.size.csv")
  safe_yaml_write(cfg, file.path(outdir, "used_config.yml"))
  check <- app_env$validate_boral(dat_case$Y, dat$X, traits = dat$traits,
                                  rowids = dat$rowids, ranefids = dat$ranefids,
                                  distmat = dat$distmat, offset = dat$offset,
                                  family = cfg$model$family,
                                  family_text = cfg$model$family_vector_override,
                                  num.lv = cfg$model$lv.control$num.lv,
                                  lv.type = cfg$model$lv.control$type,
                                  row.eff = cfg$model$row.eff,
                                  n.burnin = cfg$mcmc_prior$n.burnin,
                                  n.iteration = cfg$mcmc_prior$n.iteration,
                                  n.thin = cfg$mcmc_prior$n.thin,
                                  trial.size = cfg$model$trial.size,
                                  calc.ics = cfg$mcmc_prior$calc.ics,
                                  use_traits = cfg$traits_random_ssvs$use_traits,
                                  which.traits = cfg$traits_random_ssvs$which.traits,
                                  use_ssvs = cfg$traits_random_ssvs$use_ssvs,
                                  ssvs.index = cfg$traits_random_ssvs$ssvs.index,
                                  save.model = cfg$traits_random_ssvs$save.model,
                                  do.fit = cfg$model$do.fit,
                                  formula.X = cfg$model$formula.X)
  app_env$write_data_check_messages(outdir, check$messages)
  if (!isTRUE(check$ok)) {
    mark_case_status(outdir, "boral", "check_failed", warnings = check$messages,
                     errors = "boral synthetic case check failed.", Y = dat_case$Y, X = dat_case$X)
  } else {
    app_env$write_boral_reproducible_script(outdir)
    code <- runner_env$run_rscript(file.path(outdir, "reproducible_script", "run_this_boral_analysis.R"), outdir, "boral")
    if (!identical(as.integer(code), 0L) && !(safe_read_status(outdir, "boral")$status %in% c("fitted", "fit_failed", "check_failed", "model_defined"))) {
      mark_case_status(outdir, "boral", "fit_failed", warnings = check$messages,
                       errors = paste0("boral Rscript exited with status ", code, "."),
                       Y = dat_case$Y, X = dat_case$X)
    }
  }
  finish_case("boral", case_id, outdir, dat_case, started, notes)
}

run_suite <- function() {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  data_dir <- file.path(app_dir, "examples", "workflow_synthetic_case_suite")
  runner_env$generate_benchmark_data(data_dir, n_sites = 30L, n_species = 6L, n_visits = 3L, seed = 20260811L)
  dat <- runner_env$load_benchmark_data(data_dir, allow_generate = FALSE)
  root <- file.path(app_dir, "output", paste0("workflow_synthetic_case_audit_", timestamp))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  results <- list()
  add_result <- function(expr) {
    idx <- length(results) + 1L
    results[[idx]] <<- tryCatch(expr, error = function(e) {
      data.frame(engine = NA_character_, case_id = paste0("suite_error_", idx), status = "fit_failed",
                 output_folder = root, zip = NA_character_, zip_size_kb = NA_real_,
                 file_count = NA_integer_, missing_count = NA_integer_, missing_items = "",
                 failure_reason = conditionMessage(e), warnings = "", notes = "Unhandled suite-level error.",
                 runtime_seconds = NA_real_, stringsAsFactors = FALSE)
    })
  }

  hcfg <- runner_env$hmsc_cfg()
  hcfg$mcmc$samples <- 4L; hcfg$mcmc$transient <- 4L; hcfg$mcmc$nChains <- 1L
  hcfg$outputs$cv <- FALSE; hcfg$outputs$waic <- FALSE
  add_result(run_hmsc_case(root, "hmsc_probit_sample_traits_phylo", dat, hcfg, dat$Y,
                           "probit, sample random effect, traits, phylogenetic covariance."))
  hcfg2 <- hcfg
  hcfg2$model$distr <- "poisson"; hcfg2$model$use_traits <- FALSE; hcfg2$model$use_phylogeny <- FALSE
  hcfg2$model$random_mode <- "spatial_nngp"; hcfg2$model$spatial_method <- "NNGP"; hcfg2$model$nNeighbours <- 3L
  add_result(run_hmsc_case(root, "hmsc_poisson_spatial_nngp", dat, hcfg2, dat$Y_count,
                           "poisson counts, spatial NNGP random level, no traits/phylogeny."))
  hcfg_full <- hcfg
  hcfg_full$model$distr <- "probit"; hcfg_full$model$use_traits <- FALSE; hcfg_full$model$use_phylogeny <- FALSE
  hcfg_full$model$random_mode <- "spatial_full"; hcfg_full$model$spatial_method <- "Full"
  add_result(run_hmsc_case(root, "hmsc_probit_spatial_full", dat, hcfg_full, dat$Y,
                           "probit response, spatial Full Gaussian process random level."))
  hcfg3 <- hcfg
  hcfg3$model$distr <- "probit"; hcfg3$model$random_mode <- "spatial_gpp"; hcfg3$model$spatial_method <- "GPP"
  hcfg3$model$use_phylogeny <- FALSE
  add_result(run_hmsc_case(root, "hmsc_probit_spatial_gpp_traits", dat, hcfg3, dat$Y,
                           "probit response, spatial GPP, traits enabled."))
  hcfg4 <- hcfg
  hcfg4$model$distr <- "normal"; hcfg4$model$random_mode <- "none"
  hcfg4$model$use_traits <- FALSE; hcfg4$model$use_phylogeny <- FALSE
  add_result(run_hmsc_case(root, "hmsc_normal_no_random_guardrail", dat, hcfg4, dat$Y_normal,
                           "normal response, no random level; stable normal-family guardrail."))

  hpc <- runner_env$hmschpc_cfg(app_env)
  hpc$sampler$verbose <- 1L; hpc$sampler$samples <- 2L; hpc$sampler$transient <- 2L; hpc$sampler$chains <- 1L
  hpc$random_effects$nfMax <- 3L
  add_result(run_hmschpc_case(root, "hmschpc_probit_iid_traits_phylo", dat, hpc, dat$Y,
                              "CPU Hmsc-HPC probit, iid random level, traits, phylogenetic covariance."))
  hpc2 <- hpc
  hpc2$model$distribution <- "poisson"; hpc2$model$use_traits <- FALSE; hpc2$model$phylogeny_mode <- "none"
  hpc2$random_effects$mode <- "none"
  add_result(run_hmschpc_case(root, "hmschpc_poisson_no_random", dat, hpc2, dat$Y_count,
                              "CPU Hmsc-HPC poisson, no random level."))
  hpc3 <- hpc
  hpc3$model$distribution <- "normal"; hpc3$model$phylogeny_mode <- "none"
  hpc3$random_effects$mode <- "spatial_full"; hpc3$random_effects$name <- "spatial_site"
  add_result(run_hmschpc_case(root, "hmschpc_normal_spatial_full", dat, hpc3, dat$Y_normal,
                              "CPU Hmsc-HPC normal, spatial_full random level."))

  jcfg <- runner_env$jsdm_cfg()
  add_result(run_jsdm_case(root, "jsdm_binomial_probit_latent_random_traits", dat, jcfg, dat$Y,
                           notes = "jSDM probit matrix, latent variables, random site effect, traits."))
  jcfg2 <- jcfg
  jcfg2$model$model_type <- "binomial_logit"; jcfg2$model$n_latent <- 0L; jcfg2$model$site_effect <- "none"
  jcfg2$model$allow_traits <- FALSE
  add_result(run_jsdm_case(root, "jsdm_binomial_logit_trials", dat, jcfg2, dat$Y,
                           trials = dat$trial, notes = "jSDM binomial logit with trials."))
  jcfg3 <- jcfg
  jcfg3$model$model_type <- "poisson_log"; jcfg3$model$n_latent <- 0L; jcfg3$model$site_effect <- "none"
  jcfg3$model$allow_traits <- FALSE
  add_result(run_jsdm_case(root, "jsdm_poisson_log_counts", dat, jcfg3, dat$Y_count,
                           notes = "jSDM poisson log counts."))
  jcfg5 <- jcfg
  jcfg5$model$model_type <- "gaussian"; jcfg5$model$n_latent <- 0L; jcfg5$model$site_effect <- "none"
  jcfg5$model$allow_traits <- FALSE
  add_result(run_jsdm_case(root, "jsdm_gaussian_continuous", dat, jcfg5, dat$Y_normal,
                           notes = "jSDM gaussian continuous-response branch."))
  jcfg4 <- jcfg
  jcfg4$model$model_type <- "binomial_probit_long_format"; jcfg4$model$site_effect <- "none"
  jcfg4$model$n_latent <- 0L; jcfg4$model$allow_traits <- FALSE
  jcfg4$model$site_formula <- "~ species + species:pH + species:moisture + species:canopy + species:elevation + species:substrate"
  long <- make_long_format(dat$Y, dat$X)
  add_result(run_jsdm_case(root, "jsdm_binomial_probit_long_format", dat, jcfg4, dat$Y,
                           long = long, notes = "jSDM long-format probit branch."))

  gcfg <- runner_env$gjam_cfg(colnames(dat$Y))
  add_result(run_gjam_case(root, "gjam_pa_traits_prediction", dat, gcfg, dat$Y,
                           rep("PA", ncol(dat$Y)), "GJAM presence/absence with traits and prediction."))
  gcfg_stable <- runner_env$gjam_cfg(colnames(dat$Y))
  stable_mixed_y <- data.frame(PA = dat$Y[, 1], CON = dat$Y_normal[, min(2, ncol(dat$Y_normal))],
                               DA = dat$Y_count[, min(3, ncol(dat$Y_count))],
                               stringsAsFactors = FALSE, check.names = FALSE)
  names(stable_mixed_y) <- colnames(dat$Y)[seq_len(ncol(stable_mixed_y))]
  rownames(stable_mixed_y) <- rownames(dat$Y)
  stable_types <- c("PA", "CON", "DA")
  gcfg_stable$model$formula <- "~ pH + moisture"
  gcfg_stable$model$type_single <- "PA"
  gcfg_stable$model$typeNames_text <- paste(stable_types, collapse = ",")
  gcfg_stable$analysis$do_traits <- FALSE
  add_result(run_gjam_case(root, "gjam_mixed_pa_con_da_stable", dat, gcfg_stable, stable_mixed_y,
                           stable_types, "GJAM stable mixed response types: PA, CON and DA."))
  gcfg2 <- runner_env$gjam_cfg(colnames(dat$Y))
  ordinal <- as.integer(cut(dat$Y_normal[, min(4, ncol(dat$Y_normal))],
                            breaks = unique(quantile(dat$Y_normal[, min(4, ncol(dat$Y_normal))],
                                                     probs = seq(0, 1, length.out = 5), na.rm = TRUE)),
                            include.lowest = TRUE))
  ordinal[is.na(ordinal)] <- 1L
  ca_raw <- dat$Y_normal[, min(3, ncol(dat$Y_normal))]
  ca <- pmax(0, ca_raw - stats::quantile(ca_raw, probs = 0.35, na.rm = TRUE))
  oc_y <- data.frame(PA = dat$Y[, 1],
                     CON = dat$Y_normal[, min(2, ncol(dat$Y_normal))],
                     CA = ca,
                     OC = ordinal,
                     stringsAsFactors = FALSE, check.names = FALSE)
  names(oc_y) <- colnames(dat$Y)[seq_len(ncol(oc_y))]
  rownames(oc_y) <- rownames(dat$Y)
  oc_types <- c("PA", "CON", "CA", "OC")
  gcfg2$model$formula <- "~ pH + moisture"
  gcfg2$model$type_single <- "PA"
  gcfg2$model$typeNames_text <- paste(oc_types, collapse = ",")
  gcfg2$analysis$do_traits <- FALSE
  add_result(run_gjam_case(root, "gjam_mixed_pa_con_ca_oc_stable", dat, gcfg2, oc_y,
                           oc_types, "GJAM stable mixed response types: PA, CON, CA and OC."))
  fc_sim <- gjam::gjamSimData(n = nrow(dat$Y), S = min(4, ncol(dat$Y)), Q = 2, typeNames = "FC")
  fc_y <- as.data.frame(fc_sim$ydata, check.names = FALSE)
  names(fc_y) <- colnames(dat$Y)[seq_len(ncol(fc_y))]
  rownames(fc_y) <- rownames(dat$Y)
  dat_fc <- dat
  dat_fc$X <- as.data.frame(fc_sim$xdata, check.names = FALSE)
  rownames(dat_fc$X) <- rownames(dat$Y)
  dat_fc$newdata <- dat_fc$X[seq_len(min(8, nrow(dat_fc$X))), , drop = FALSE]
  rownames(dat_fc$newdata) <- sprintf("new_site_%02d", seq_len(nrow(dat_fc$newdata)))
  gcfg_fc <- runner_env$gjam_cfg(colnames(fc_y))
  gcfg_fc$model$formula <- paste(deparse(fc_sim$formula), collapse = "")
  gcfg_fc$model$type_single <- "FC"
  gcfg_fc$model$typeNames_text <- paste(as.character(fc_sim$typeNames), collapse = ",")
  gcfg_fc$response_types$FCgroups <- paste(attr(fc_sim$typeNames, "FCgroups") %||% rep(1, ncol(fc_y)), collapse = ",")
  gcfg_fc$analysis$do_traits <- FALSE
  add_result(run_gjam_case(root, "gjam_fractional_fc_single_response", dat_fc, gcfg_fc, fc_y,
                           rep("FC", ncol(fc_y)), "GJAM fractional-composition response branch with one stable FC group."))
  gcfg3 <- runner_env$gjam_cfg(colnames(dat$Y))
  gcfg3$model$type_single <- "DA"; gcfg3$model$typeNames_text <- paste(rep("DA", ncol(dat$Y_count)), collapse = ",")
  gcfg3$analysis$do_traits <- FALSE
  add_result(run_gjam_case(root, "gjam_discrete_abundance_counts", dat, gcfg3, dat$Y_count,
                           rep("DA", ncol(dat$Y_count)), "GJAM DA count-style response branch."))

  scfg <- runner_env$spocc_cfg()
  add_result(run_spocc_case(root, "spocc_msPGOcc_replicated_multi", dat, scfg, dat$det,
                            notes = "Multi-species replicated occupancy."))
  scfg2 <- scfg
  scfg2$model$model_type <- "PGOcc"; scfg2$model$data_structure <- "single-species replicated"
  scfg2$validation_prediction_outputs$waicOcc <- FALSE
  single_y <- make_single_detection(dat, colnames(dat$Y)[1])
  add_result(run_spocc_case(root, "spocc_PGOcc_single_species", dat, scfg2, single_y,
                            species_table = data.frame(species = colnames(dat$Y)[1], stringsAsFactors = FALSE),
                            notes = "Single-species non-spatial PGOcc."))
  scfg3 <- scfg2
  scfg3$model$model_type <- "spPGOcc"; scfg3$spatial_latent_svc$NNGP <- TRUE
  add_result(run_spocc_case(root, "spocc_spPGOcc_single_spatial", dat, scfg3, single_y,
                            species_table = data.frame(species = colnames(dat$Y)[1], stringsAsFactors = FALSE),
                            notes = "Single-species spatial PGOcc with coordinates."))
  scfg4 <- scfg
  scfg4$model$model_type <- "lfMsPGOcc"; scfg4$spatial_latent_svc$n.factors <- 2L
  add_result(run_spocc_case(root, "spocc_lfMsPGOcc_latent_factors", dat, scfg4, dat$det,
                            notes = "Multi-species latent-factor occupancy branch."))

  sj <- runner_env$sjsdm_cfg()
  add_result(run_sjsdm_case(root, "sjsdm_binomial_linear_spatial", dat, sj, dat$Y,
                            "sjSDM binomial-probit linear environment and spatial terms."))
  sj2 <- sj
  sj2$model$family <- "gaussian_identity"; sj2$model$spatial_model <- "none"; sj2$outputs$plots <- TRUE
  add_result(run_sjsdm_case(root, "sjsdm_gaussian_no_spatial", dat, sj2, dat$Y_normal,
                            "sjSDM gaussian response without spatial model."))
  sj3 <- sj
  sj3$model$env_model <- "DNN"; sj3$dnn_optimizer$hidden <- "4"; sj3$model$iter <- 3L; sj3$model$sampling <- 50L
  add_result(run_sjsdm_case(root, "sjsdm_binomial_dnn_environment", dat, sj3, dat$Y,
                            "sjSDM DNN environment branch on CPU."))

  bcfg <- runner_env$boral_cfg()
  add_result(run_boral_case(root, "boral_binomial_latent", dat, bcfg, dat$Y,
                            "boral binomial with one latent variable."))
  bcfg2 <- bcfg
  bcfg2$model$family <- "poisson"; bcfg2$model$row.eff <- "fixed"
  add_result(run_boral_case(root, "boral_poisson_row_effect", dat, bcfg2, dat$Y_count,
                            "boral poisson with fixed row effects."))
  bcfg3 <- bcfg
  bcfg3$model$family <- "normal"; bcfg3$model$use_offset <- TRUE; bcfg3$model$row.eff <- "none"
  add_result(run_boral_case(root, "boral_normal_offset", dat, bcfg3, dat$Y_normal,
                            "boral normal response with offset file present."))

  summary <- do.call(rbind, results)
  write.csv(summary, file.path(root, "workflow_synthetic_case_summary.csv"), row.names = FALSE)
  by_engine <- aggregate(case_id ~ engine + status, data = summary, FUN = length)
  names(by_engine)[names(by_engine) == "case_id"] <- "n_cases"
  write.csv(by_engine, file.path(root, "workflow_synthetic_case_status_by_engine.csv"), row.names = FALSE)
  report <- c(
    "# JSDM Studio Workflow Synthetic Case Audit",
    "",
    paste0("Generated: ", Sys.time()),
    paste0("Audit root: ", normalizePath(root, winslash = "/", mustWork = FALSE)),
    "",
    "## Status By Engine",
    "",
    paste(capture.output(print(by_engine, row.names = FALSE)), collapse = "\n"),
    "",
    "## Case Summary",
    "",
    paste(capture.output(print(summary[, c("engine", "case_id", "status", "file_count", "missing_count", "zip_size_kb", "failure_reason")],
                               row.names = FALSE)), collapse = "\n"),
    "",
    "Failures with status check_failed or fit_failed are acceptable only when diagnostics are written and ZIP/output contracts remain complete. They should not be displayed as Completed in the Shiny UI."
  )
  writeLines(report, file.path(root, "workflow_synthetic_case_audit_report.md"))
  cat("AUDIT_ROOT=", normalizePath(root, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("SUMMARY_CSV=", normalizePath(file.path(root, "workflow_synthetic_case_summary.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  invisible(summary)
}

if (!interactive()) run_suite()
