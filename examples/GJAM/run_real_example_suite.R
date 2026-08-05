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
sys.source(file.path("R", "gjam_adapter.R"), app_env)

suppressPackageStartupMessages(library(gjam))
rscript_bin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

make_sim <- function(types, n = 28, Q = 3, seed = 1, effort = NULL) {
  set.seed(seed)
  gjam::gjamSimData(n = n, S = length(types), Q = Q, typeNames = types, effort = effort)
}

formula_text <- function(f) paste(deparse(f), collapse = " ")

type_table <- function(Y, types) {
  data.frame(response = colnames(Y), typeName = as.character(types), stringsAsFactors = FALSE)
}

trait_inputs <- function(Y) {
  s <- ncol(Y)
  specByTrait <- data.frame(
    height = round(seq(0.5, 2.5, length.out = s), 3),
    guild = factor(rep(c("forb", "grass", "shrub"), length.out = s)),
    stringsAsFactors = FALSE
  )
  rownames(specByTrait) <- colnames(Y)
  traitTypes <- data.frame(trait = names(specByTrait), typeName = c("CON", "CAT"), stringsAsFactors = FALSE)
  list(specByTrait = specByTrait, traitTypes = traitTypes)
}

base_config <- function(case) {
  list(
    project_name = paste0("GJAM_suite_", case$id),
    engine = "GJAM",
    data = list(
      Y = "Y.csv", XData = "XData.csv", typeNames = "typeNames.csv",
      censor = "censor.csv", effort = "effort.csv", newdata = "newdata.csv",
      specByTrait = "specByTrait.csv", traitTypes = "traitTypes.csv",
      holdoutIndex = "holdoutIndex.csv"
    ),
    model = list(
      formula = case$formula %||% "~ .",
      type_single = case$type_single %||% "DA",
      typeNames_text = case$typeNames_text %||% "",
      notStandard = case$notStandard %||% "",
      ng = case$ng %||% 100,
      burnin = case$burnin %||% 20,
      holdoutN = case$holdoutN %||% 0,
      seed = case$seed %||% 1234,
      random = case$random %||% "",
      FULL = case$FULL %||% FALSE,
      PREDICTX = case$PREDICTX %||% TRUE,
      REDUCT = case$REDUCT %||% FALSE,
      reductList = list(N = case$reduct_N %||% 20, r = case$reduct_r %||% 3),
      ematAlpha = case$ematAlpha %||% 0.5,
      use_censor = !is.null(case$censor) || isTRUE(case$use_censor),
      use_effort = !is.null(case$effort)
    ),
    response_types = list(
      FCgroups = case$FCgroups %||% "",
      CCgroups = case$CCgroups %||% "",
      composition_reference = case$composition_reference %||% "adapter default",
      trimY = case$trimY %||% FALSE,
      trim_minObs = case$trim_minObs %||% 2
    ),
    priors_censoring = list(
      use_prior_template = case$use_prior_template %||% FALSE,
      prior_file = if (!is.null(case$prior)) "priorTemplate.csv" else NULL,
      prior_mode = case$prior_mode %||% "non-informative/default",
      censor_columns = case$censor_columns %||% "",
      censor_values = case$censor_values %||% "",
      censor_intervals = case$censor_intervals %||% ""
    ),
    analysis = list(
      do_predict = case$do_predict %||% TRUE,
      do_sensitivity = case$do_sensitivity %||% TRUE,
      do_ordination = case$do_ordination %||% FALSE,
      do_conditional = case$do_conditional %||% FALSE,
      do_iie = case$do_iie %||% FALSE,
      do_traits = !is.null(case$specByTrait),
      missingX = case$missingX %||% TRUE,
      missingY = case$missingY %||% TRUE,
      inverse_prediction = case$inverse_prediction %||% TRUE
    ),
    outputs = list(
      real_fit = TRUE, save_model = TRUE, save_chains = TRUE,
      save_parameters = TRUE, save_fit = TRUE, save_prediction = TRUE,
      save_missing = TRUE, save_plots = case$save_plots %||% TRUE,
      report = TRUE, zip = TRUE, copy_inputs = TRUE, save_config = TRUE,
      csv_tables = TRUE
    )
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
  write_pair(outdir, "typeNames.csv", type_table(case$Y, case$type_vec), row.names = FALSE)
  write_pair(outdir, "censor.csv", case$censor, row.names = FALSE)
  write_pair(outdir, "effort.csv", case$effort, row.names = FALSE)
  write_pair(outdir, "newdata.csv", case$newdata, row.names = TRUE)
  write_pair(outdir, "specByTrait.csv", case$specByTrait, row.names = TRUE)
  write_pair(outdir, "traitTypes.csv", case$traitTypes, row.names = FALSE)
  write_pair(outdir, "holdoutIndex.csv", case$holdoutIndex, row.names = FALSE)
  write_pair(outdir, "priorTemplate.csv", case$prior, row.names = FALSE)
}

run_case <- function(case) {
  outdir <- app_env$make_engine_run_dir("GJAM", paste0("JSDMStudio_GJAM_SUITE_", case$id))
  write_case_files(outdir, case)
  cfg <- base_config(case)
  yaml::write_yaml(cfg, file.path(outdir, "used_config.yml"))
  check <- app_env$validate_gjam_full(
    case$Y, case$X, type_table(case$Y, case$type_vec), cfg$model$typeNames_text,
    cfg$model$type_single, cfg$response_types$FCgroups, cfg$response_types$CCgroups,
    cfg$model$ng, cfg$model$burnin, cfg$model$holdoutN, cfg$model$random,
    cfg$model$notStandard, cfg$model$formula, case$censor, case$effort,
    case$newdata, case$specByTrait, case$traitTypes, case$holdoutIndex,
    case$prior, cfg$model$use_censor, cfg$model$use_effort,
    cfg$analysis$do_traits, cfg$response_types$trimY, cfg$response_types$trim_minObs,
    cfg$model$REDUCT, cfg$model$reductList$N, cfg$model$reductList$r
  )
  app_env$write_data_check_messages(outdir, check$messages)
  write.csv(app_env$data_summary(case$Y, case$X), file.path(outdir, "tables", "data_summary.csv"), row.names = FALSE)
  if (!isTRUE(check$ok)) {
    writeLines(check$messages, file.path(outdir, "diagnostics", "GJAM_check_failed.txt"))
    stop(sprintf("Case %s failed validation: %s", case$id, paste(check$messages, collapse = "; ")), call. = FALSE)
  }
  app_env$write_gjam_reproducible_script(outdir)
  log_file <- file.path(outdir, "diagnostics", "suite_stdout_stderr.txt")
  t0 <- Sys.time()
  exit_code <- system2(rscript_bin, shQuote(file.path(outdir, "reproducible_script", "run_this_GJAM_analysis.R")),
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
    "results/README_GJAM_results.txt", "results/run_summary.csv", "results/fit_metrics.csv",
    "tables/typeNames_used.csv", "tables/betaMu.csv", "tables/corMu.csv",
    "models/gjam_model.rds", "models/gjam_modelList.rds",
    "chains/chain_manifest.csv", "predictions/gjam_prediction_from_fit.rds",
    "report/GJAM_report.html",
    "reproducible_script/run_this_GJAM_analysis.R", "workflow_scripts/run_GJAM_workflow.R"
  )
  missing <- required[!file.exists(file.path(outdir, required))]
  small <- required[file.exists(file.path(outdir, required)) & file.info(file.path(outdir, required))$size <= 0]
  data.frame(
    case_id = case$id,
    typeNames = paste(case$type_vec, collapse = ","),
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

sim1 <- make_sim(c("DA","DA","OC","OC","CON","CA","PA"), n = 30, Q = 3, seed = 101)
sim2_eff <- list(columns = 1:4, values = round(runif(32, 0.6, 4), 2))
sim2 <- make_sim(rep("DA", 4), n = 32, Q = 3, seed = 102, effort = sim2_eff)
sim2$xdata$plot <- factor(rep(seq_len(8), each = 4))
sim3 <- make_sim(rep("CA", 4), n = 30, Q = 3, seed = 103)
sim4 <- make_sim(rep("FC", 3), n = 30, Q = 3, seed = 104)
sim5 <- make_sim(rep("CC", 3), n = 30, Q = 3, seed = 105)
sim6 <- make_sim(rep("CAT", 3), n = 30, Q = 3, seed = 106)
sim7 <- make_sim(rep("DA", 4), n = 28, Q = 3, seed = 107)
traits7 <- trait_inputs(sim7$ydata)
sim7$xdata[2, "x2"] <- NA
sim7$ydata[3, 1] <- NA

cases <- list(
  list(id = "01_mixed_scales_default_formula_predict_sensitivity",
       Y = sim1$ydata, X = sim1$xdata, type_vec = sim1$typeNames,
       formula = "~ .", newdata = sim1$xdata[1:6, , drop = FALSE],
       do_sensitivity = TRUE, do_predict = TRUE, save_plots = TRUE),
  list(id = "02_DA_effort_random_FULL_notStandard",
       Y = sim2$ydata, X = sim2$xdata, type_vec = sim2$typeNames,
       formula = "~ x2 + x3 + plot", effort = data.frame(effort = sim2_eff$values),
       random = "plot", FULL = TRUE, PREDICTX = FALSE, notStandard = "x2",
       do_predict = TRUE, save_plots = FALSE),
  list(id = "03_CA_censor_holdout_conditional",
       Y = sim3$ydata, X = sim3$xdata, type_vec = sim3$typeNames,
       formula = formula_text(sim3$formula), censor = data.frame(column = colnames(sim3$ydata)[1], value = 0.001, lower = -Inf, upper = 0.001),
       censor_columns = colnames(sim3$ydata)[1], censor_values = "0.001", censor_intervals = "-Inf,0.001",
       holdoutN = 3, use_prior_template = TRUE, prior_mode = "sign-constrained",
       do_conditional = TRUE, do_predict = TRUE),
  list(id = "04_FC_composition_groups_ordination",
       Y = sim4$ydata, X = sim4$xdata, type_vec = sim4$typeNames,
       formula = formula_text(sim4$formula), FCgroups = paste(attr(sim4$typeNames, "FCgroups"), collapse = ","),
       holdoutIndex = data.frame(index = c(1, 4, 7)), do_ordination = TRUE, do_predict = TRUE),
  list(id = "05_CC_composition_groups_iie",
       Y = sim5$ydata, X = sim5$xdata, type_vec = sim5$typeNames,
       formula = formula_text(sim5$formula), CCgroups = paste(attr(sim5$typeNames, "CCgroups"), collapse = ","),
       do_iie = TRUE, do_predict = TRUE),
  list(id = "06_CAT_categorical_safe_prior_skip",
       Y = sim6$ydata, X = sim6$xdata, type_vec = sim6$typeNames,
       formula = formula_text(sim6$formula), use_prior_template = TRUE,
       prior_mode = "sign-constrained", do_predict = TRUE, do_sensitivity = TRUE),
  list(id = "07_traits_missing_trimY",
       Y = sim7$ydata, X = sim7$xdata, type_vec = sim7$typeNames,
       formula = formula_text(sim7$formula), specByTrait = traits7$specByTrait,
       traitTypes = traits7$traitTypes, do_predict = TRUE, do_traits = TRUE,
       trimY = TRUE, trim_minObs = 2, missingX = TRUE, missingY = TRUE)
)

results <- do.call(rbind, lapply(cases, run_case))
summary_file <- file.path("output", paste0("GJAM_real_example_suite_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
write.csv(results, summary_file, row.names = FALSE)
print(results[, c("case_id", "typeNames", "exit_code", "status", "pass", "missing_required", "empty_required")])
cat("SUMMARY_FILE=", normalizePath(summary_file, winslash = "/", mustWork = FALSE), "\n", sep = "")

if (!all(results$pass)) {
  failed <- results[!results$pass, , drop = FALSE]
  stop(sprintf("GJAM suite failed for %d case(s). See summary: %s", nrow(failed), normalizePath(summary_file, winslash = "/", mustWork = FALSE)), call. = FALSE)
}
