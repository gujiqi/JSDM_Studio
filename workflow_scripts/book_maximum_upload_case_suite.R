# Maximum-upload branch audit for the converted book chapter cases.
# It uses the GUI upload-ready folders as the input source and the real
# workflow adapters as the fitting source.

options(stringsAsFactors = FALSE)

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_arg <- args0[startsWith(args0, file_arg)]
script_path <- if (length(script_arg)) {
  normalizePath(sub(file_arg, "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
} else {
  normalizePath("workflow_scripts/book_maximum_upload_case_suite.R", winslash = "/", mustWork = FALSE)
}
app_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) stop("Cannot locate JSDM Studio app.R.", call. = FALSE)
setwd(app_dir)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(cases = "ch06_plant_traits_whittaker,ch07_deadwood_fungi,ch11_finnish_birds",
              smoke = FALSE)
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "true"
    out[[key]] <- val
  }
  out$cases <- trimws(unlist(strsplit(out$cases, ",")))
  out$cases <- out$cases[nzchar(out$cases)]
  out$smoke <- tolower(as.character(out$smoke)) %in% c("true", "t", "1", "yes", "y")
  out
}

load_suite_env <- function() {
  suite_path <- file.path(app_dir, "workflow_scripts", "workflow_synthetic_case_suite.R")
  txt <- readLines(suite_path, warn = FALSE)
  txt <- txt[!grepl("^\\s*if \\(!interactive\\(\\)\\) run_suite\\(\\)\\s*$", txt)]
  e <- new.env(parent = globalenv())
  eval(parse(text = paste(txt, collapse = "\n")), envir = e)
  e
}

senv <- load_suite_env()
runner_env <- senv$runner_env
app_env <- senv$app_env

upload_root <- file.path(app_dir, "examples", "00_UPLOAD_READY_book_chapter_cases")
raw_root <- file.path(app_dir, "examples", "book_chapter_benchmarks")

read_csv_auto <- function(path, row.names = 1) {
  if (!file.exists(path)) stop("Missing upload file: ", path, call. = FALSE)
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(x) && !ncol(x)) return(x)
  if (!is.null(row.names) && ncol(x) >= row.names) {
    rn <- x[[row.names]]
    if (!anyDuplicated(rn) && all(nzchar(as.character(rn)))) {
      rownames(x) <- as.character(rn)
      x <- x[-row.names]
    }
  }
  x
}

read_text <- function(path) {
  if (!file.exists(path)) stop("Missing upload file: ", path, call. = FALSE)
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

load_engine_upload_dat <- function(case_id, engine_folder) {
  base <- runner_env$load_benchmark_data(file.path(raw_root, case_id), allow_generate = FALSE)
  case_dir <- file.path(upload_root, case_id)
  eng_dir <- file.path(case_dir, engine_folder)
  if (!dir.exists(eng_dir)) stop("Missing engine upload folder: ", eng_dir, call. = FALSE)
  dat <- base
  if (engine_folder == "01_Hmsc") {
    dat$Y <- read_csv_auto(file.path(eng_dir, "Y.csv"))
    dat$X <- read_csv_auto(file.path(eng_dir, "XData.csv"))
    dat$traits <- read_csv_auto(file.path(eng_dir, "TrData.csv"))
    dat$study <- read_csv_auto(file.path(eng_dir, "studyDesign.csv"))
    dat$coords <- read_csv_auto(file.path(eng_dir, "coordinates.csv"))
    dat$phylo_cov <- read_csv_auto(file.path(eng_dir, "phylo_cov.csv"))
    dat$newick <- read_text(file.path(eng_dir, "phylogeny.nwk"))
  } else if (engine_folder == "02_Hmsc-HPC") {
    dat$Y <- read_csv_auto(file.path(eng_dir, "Y.csv"))
    dat$X <- read_csv_auto(file.path(eng_dir, "XData.csv"))
    dat$traits <- read_csv_auto(file.path(eng_dir, "traits.csv"))
    dat$study <- read_csv_auto(file.path(eng_dir, "studyDesign.csv"))
    dat$coords <- read_csv_auto(file.path(eng_dir, "coordinates.csv"))
    dat$phylo_cov <- read_csv_auto(file.path(eng_dir, "phylo_cov.csv"))
    dat$newdata <- read_csv_auto(file.path(eng_dir, "newdata.csv"))
  } else if (engine_folder == "03_jSDM") {
    dat$Y <- read_csv_auto(file.path(eng_dir, "Y.csv"))
    dat$X <- read_csv_auto(file.path(eng_dir, "XData.csv"))
    dat$traits <- read_csv_auto(file.path(eng_dir, "trait_data.csv"))
    dat$trial <- read_csv_auto(file.path(eng_dir, "trials.csv"))
    dat$newdata <- read_csv_auto(file.path(eng_dir, "newdata.csv"))
    dat$long_format <- read_csv_auto(file.path(eng_dir, "long_format.csv"), row.names = NULL)
    dat$prediction_ids <- read_csv_auto(file.path(eng_dir, "prediction_ids.csv"), row.names = NULL)
  } else if (engine_folder == "04_GJAM") {
    dat$Y <- read_csv_auto(file.path(eng_dir, "Y.csv"))
    dat$X <- read_csv_auto(file.path(eng_dir, "XData.csv"))
    dat$traits <- read_csv_auto(file.path(eng_dir, "specByTrait.csv"))
    dat$typeNames <- read_csv_auto(file.path(eng_dir, "typeNames.csv"), row.names = NULL)
    dat$traitTypes <- read_csv_auto(file.path(eng_dir, "traitTypes.csv"), row.names = NULL)
    dat$effort <- read_csv_auto(file.path(eng_dir, "effort.csv"))
    dat$newdata <- read_csv_auto(file.path(eng_dir, "newdata.csv"))
    dat$holdoutIndex <- read_csv_auto(file.path(eng_dir, "holdoutIndex.csv"), row.names = NULL)
    dat$censor <- read_csv_auto(file.path(eng_dir, "censor.csv"), row.names = NULL)
  } else if (engine_folder == "05_spOccupancy") {
    dat$det <- read_csv_auto(file.path(eng_dir, "y.csv"))
    dat$occ_covs <- read_csv_auto(file.path(eng_dir, "occ.covs.csv"))
    dat$det_cov <- read_csv_auto(file.path(eng_dir, "det.covs.csv"))
    dat$coords <- read_csv_auto(file.path(eng_dir, "coords.csv"))
    dat$species_table <- read_csv_auto(file.path(eng_dir, "species.csv"), row.names = NULL)
    dat$integrated_sources <- read_csv_auto(file.path(eng_dir, "integrated_sources.csv"), row.names = NULL)
    dat$newdata <- read_csv_auto(file.path(eng_dir, "newdata.csv"))
    dat$newcoords <- read_csv_auto(file.path(eng_dir, "newcoords.csv"))
    dat$folds <- read_csv_auto(file.path(eng_dir, "folds.csv"), row.names = NULL)
  } else if (engine_folder == "06_sjSDM") {
    dat$Y <- read_csv_auto(file.path(eng_dir, "Y.csv"))
    dat$X <- read_csv_auto(file.path(eng_dir, "env.csv"))
    dat$coords <- read_csv_auto(file.path(eng_dir, "spatial.csv"))
    dat$traits <- read_csv_auto(file.path(eng_dir, "traits.csv"))
    dat$newdata <- read_csv_auto(file.path(eng_dir, "newdata.csv"))
  } else if (engine_folder == "07_boral") {
    dat$Y <- read_csv_auto(file.path(eng_dir, "Y.csv"))
    dat$X <- read_csv_auto(file.path(eng_dir, "XData.csv"))
    dat$traits <- read_csv_auto(file.path(eng_dir, "traits.csv"))
    dat$rowids <- read_csv_auto(file.path(eng_dir, "row.ids.csv"))
    dat$ranefids <- read_csv_auto(file.path(eng_dir, "ranef.ids.csv"))
    dat$distmat <- read_csv_auto(file.path(eng_dir, "distmat.csv"))
    dat$offset <- read_csv_auto(file.path(eng_dir, "offset.csv"))
    dat$newdata <- read_csv_auto(file.path(eng_dir, "newdata.csv"))
    dat$trial <- read_csv_auto(file.path(eng_dir, "trial.size.csv"))
  }
  dat$upload_engine_folder <- eng_dir
  dat
}

copy_max_upload_files <- function(outdir, dat) {
  src <- dat$upload_engine_folder
  files <- list.files(src, full.names = TRUE, recursive = FALSE)
  rows <- lapply(files, function(f) {
    rel <- basename(f)
    for (root_name in c("inputs", "data")) {
      dest <- file.path(outdir, root_name, "maximum_upload_files", rel)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      file.copy(f, dest, overwrite = TRUE)
    }
    data.frame(file = rel, size_bytes = file.info(f)$size, copied = TRUE, stringsAsFactors = FALSE)
  })
  tab <- if (length(rows)) do.call(rbind, rows) else data.frame(file = character(), size_bytes = numeric(), copied = logical())
  write.csv(tab, file.path(outdir, "tables", "maximum_upload_file_check.csv"), row.names = FALSE)
  writeLines(c("Maximum-upload files copied from:", normalizePath(src, winslash = "/", mustWork = FALSE),
               "These are the GUI upload-ready files used/inspected for this branch test."),
             file.path(outdir, "inputs", "README_maximum_upload.txt"))
  invisible(nrow(tab))
}

refresh_result_row <- function(res, dat) {
  outdir <- res$output_folder[1]
  copied <- copy_max_upload_files(outdir, dat)
  runner_env$fill_empty_dirs(outdir, res$engine[1])
  zipfile <- tryCatch(runner_env$make_zip_file(outdir), error = function(e) {
    writeLines(conditionMessage(e), file.path(outdir, "diagnostics", "zip_refresh_error.txt"))
    NA_character_
  })
  missing_dirs <- senv$contract_dirs[!dir.exists(file.path(outdir, senv$contract_dirs))]
  missing_files <- senv$contract_files[!file.exists(file.path(outdir, senv$contract_files))]
  res$zip <- if (!is.na(zipfile)) normalizePath(zipfile, winslash = "/", mustWork = FALSE) else NA_character_
  res$zip_size_kb <- if (!is.na(zipfile) && file.exists(zipfile)) round(file.info(zipfile)$size / 1024, 1) else NA_real_
  res$file_count <- length(list.files(outdir, recursive = TRUE, all.files = FALSE))
  res$missing_count <- length(c(missing_dirs, missing_files))
  res$missing_items <- paste(c(missing_dirs, missing_files), collapse = "; ")
  res$maximum_upload_files_copied <- copied
  res
}

add_result <- function(results, expr, dat) {
  idx <- length(results) + 1L
  results[[idx]] <- tryCatch(refresh_result_row(force(expr), dat), error = function(e) {
    data.frame(engine = NA_character_, case_id = paste0("suite_error_", idx), status = "fit_failed",
               output_folder = NA_character_, zip = NA_character_, zip_size_kb = NA_real_,
               file_count = NA_integer_, missing_count = NA_integer_, missing_items = "",
               failure_reason = conditionMessage(e), warnings = "", notes = "Unhandled suite-level error.",
               runtime_seconds = NA_real_, maximum_upload_files_copied = NA_integer_,
               stringsAsFactors = FALSE)
  })
  results
}

run_one_case <- function(root, case_id, smoke = FALSE) {
  message("Running maximum-upload case: ", case_id)
  results <- list()

  hdat <- load_engine_upload_dat(case_id, "01_Hmsc")
  hcfg <- runner_env$hmsc_cfg()
  hcfg$mcmc$samples <- 4L; hcfg$mcmc$transient <- 4L; hcfg$mcmc$nChains <- 1L
  hcfg$outputs$cv <- FALSE; hcfg$outputs$waic <- FALSE
  results <- add_result(results, senv$run_hmsc_case(root, paste(case_id, "hmsc_sample_traits_phylo", sep = "__"), hdat, hcfg, hdat$Y,
                                                    "maximum upload: sample random effect, traits, phylogenetic covariance."), hdat)
  if (!smoke) {
    hcfg_full <- hcfg; hcfg_full$model$use_traits <- FALSE; hcfg_full$model$use_phylogeny <- FALSE
    hcfg_full$model$random_mode <- "spatial_full"; hcfg_full$model$spatial_method <- "Full"
    results <- add_result(results, senv$run_hmsc_case(root, paste(case_id, "hmsc_spatial_full", sep = "__"), hdat, hcfg_full, hdat$Y,
                                                      "maximum upload: spatial Full branch."), hdat)
    hcfg_nngp <- hcfg_full; hcfg_nngp$model$distr <- "poisson"; hcfg_nngp$model$random_mode <- "spatial_nngp"
    hcfg_nngp$model$spatial_method <- "NNGP"; hcfg_nngp$model$nNeighbours <- 3L
    results <- add_result(results, senv$run_hmsc_case(root, paste(case_id, "hmsc_poisson_spatial_nngp", sep = "__"), hdat, hcfg_nngp, hdat$Y_count,
                                                      "maximum upload: poisson Y from common_inputs, spatial NNGP."), hdat)
    hcfg_gpp <- hcfg; hcfg_gpp$model$use_phylogeny <- FALSE; hcfg_gpp$model$random_mode <- "spatial_gpp"; hcfg_gpp$model$spatial_method <- "GPP"
    results <- add_result(results, senv$run_hmsc_case(root, paste(case_id, "hmsc_spatial_gpp_traits", sep = "__"), hdat, hcfg_gpp, hdat$Y,
                                                      "maximum upload: spatial GPP with traits."), hdat)
    hcfg_norm <- hcfg; hcfg_norm$model$distr <- "normal"; hcfg_norm$model$random_mode <- "none"
    hcfg_norm$model$use_traits <- FALSE; hcfg_norm$model$use_phylogeny <- FALSE
    results <- add_result(results, senv$run_hmsc_case(root, paste(case_id, "hmsc_normal_no_random", sep = "__"), hdat, hcfg_norm, hdat$Y_normal,
                                                      "maximum upload: normal Y from common_inputs, no random level."), hdat)
  }

  hpcdat <- load_engine_upload_dat(case_id, "02_Hmsc-HPC")
  hpc <- runner_env$hmschpc_cfg(app_env)
  hpc$sampler$verbose <- 1L; hpc$sampler$samples <- 2L; hpc$sampler$transient <- 2L; hpc$sampler$chains <- 1L
  hpc$random_effects$nfMax <- 3L
  results <- add_result(results, senv$run_hmschpc_case(root, paste(case_id, "hmschpc_iid_traits_phylo", sep = "__"), hpcdat, hpc, hpcdat$Y,
                                                       "maximum upload: CPU iid random level, traits, phylogenetic covariance."), hpcdat)
  if (!smoke) {
    hpc2 <- hpc; hpc2$model$distribution <- "poisson"; hpc2$model$use_traits <- FALSE
    hpc2$model$phylogeny_mode <- "none"; hpc2$random_effects$mode <- "none"
    results <- add_result(results, senv$run_hmschpc_case(root, paste(case_id, "hmschpc_poisson_no_random", sep = "__"), hpcdat, hpc2, hpcdat$Y_count,
                                                         "maximum upload: poisson Y from common_inputs, no random level."), hpcdat)
    hpc3 <- hpc; hpc3$model$distribution <- "normal"; hpc3$model$phylogeny_mode <- "none"
    hpc3$random_effects$mode <- "spatial_full"; hpc3$random_effects$name <- "spatial_site"
    results <- add_result(results, senv$run_hmschpc_case(root, paste(case_id, "hmschpc_normal_spatial_full", sep = "__"), hpcdat, hpc3, hpcdat$Y_normal,
                                                         "maximum upload: normal Y from common_inputs, spatial_full."), hpcdat)
  }

  jdat <- load_engine_upload_dat(case_id, "03_jSDM")
  jcfg <- runner_env$jsdm_cfg()
  results <- add_result(results, senv$run_jsdm_case(root, paste(case_id, "jsdm_binomial_probit_traits", sep = "__"), jdat, jcfg, jdat$Y,
                                                    notes = "maximum upload: probit matrix, latent variables, traits."), jdat)
  if (!smoke) {
    jcfg2 <- jcfg; jcfg2$model$model_type <- "binomial_logit"; jcfg2$model$n_latent <- 0L; jcfg2$model$site_effect <- "none"; jcfg2$model$allow_traits <- FALSE
    results <- add_result(results, senv$run_jsdm_case(root, paste(case_id, "jsdm_binomial_logit_trials", sep = "__"), jdat, jcfg2, jdat$Y,
                                                      trials = jdat$trial, notes = "maximum upload: binomial logit with trials."), jdat)
    jcfg3 <- jcfg; jcfg3$model$model_type <- "poisson_log"; jcfg3$model$n_latent <- 0L; jcfg3$model$site_effect <- "none"; jcfg3$model$allow_traits <- FALSE
    results <- add_result(results, senv$run_jsdm_case(root, paste(case_id, "jsdm_poisson_log", sep = "__"), jdat, jcfg3, jdat$Y_count,
                                                      notes = "maximum upload: poisson Y from common_inputs."), jdat)
    jcfg4 <- jcfg; jcfg4$model$model_type <- "gaussian"; jcfg4$model$n_latent <- 0L; jcfg4$model$site_effect <- "none"; jcfg4$model$allow_traits <- FALSE
    results <- add_result(results, senv$run_jsdm_case(root, paste(case_id, "jsdm_gaussian", sep = "__"), jdat, jcfg4, jdat$Y_normal,
                                                      notes = "maximum upload: gaussian Y from common_inputs."), jdat)
    jcfg5 <- jcfg; jcfg5$model$model_type <- "binomial_probit_long_format"; jcfg5$model$site_effect <- "none"; jcfg5$model$n_latent <- 0L; jcfg5$model$allow_traits <- FALSE
    jcfg5$model$site_formula <- "~ species + species:pH + species:moisture + species:canopy + species:elevation + species:substrate"
    results <- add_result(results, senv$run_jsdm_case(root, paste(case_id, "jsdm_long_format", sep = "__"), jdat, jcfg5, jdat$Y,
                                                      long = jdat$long_format, notes = "maximum upload: long-format file branch."), jdat)
  }

  gdat <- load_engine_upload_dat(case_id, "04_GJAM")
  gcfg <- runner_env$gjam_cfg(colnames(gdat$Y))
  results <- add_result(results, senv$run_gjam_case(root, paste(case_id, "gjam_pa_traits", sep = "__"), gdat, gcfg, gdat$Y,
                                                    rep("PA", ncol(gdat$Y)), "maximum upload: PA with traits and prediction."), gdat)
  if (!smoke) {
    gcfg2 <- runner_env$gjam_cfg(colnames(gdat$Y))
    stable_mixed_y <- data.frame(PA = gdat$Y[, 1], CON = gdat$Y_normal[, min(2, ncol(gdat$Y_normal))],
                                 DA = gdat$Y_count[, min(3, ncol(gdat$Y_count))],
                                 stringsAsFactors = FALSE, check.names = FALSE)
    names(stable_mixed_y) <- colnames(gdat$Y)[seq_len(ncol(stable_mixed_y))]
    rownames(stable_mixed_y) <- rownames(gdat$Y)
    stable_types <- c("PA", "CON", "DA")
    gcfg2$model$formula <- "~ pH + moisture"; gcfg2$model$typeNames_text <- paste(stable_types, collapse = ","); gcfg2$analysis$do_traits <- FALSE
    results <- add_result(results, senv$run_gjam_case(root, paste(case_id, "gjam_mixed_pa_con_da", sep = "__"), gdat, gcfg2, stable_mixed_y,
                                                      stable_types, "maximum upload: stable mixed PA/CON/DA."), gdat)
    gcfg3 <- runner_env$gjam_cfg(colnames(gdat$Y_count)); gcfg3$model$type_single <- "DA"; gcfg3$model$typeNames_text <- paste(rep("DA", ncol(gdat$Y_count)), collapse = ","); gcfg3$analysis$do_traits <- FALSE
    results <- add_result(results, senv$run_gjam_case(root, paste(case_id, "gjam_da_counts", sep = "__"), gdat, gcfg3, gdat$Y_count,
                                                      rep("DA", ncol(gdat$Y_count)), "maximum upload: DA count branch."), gdat)
  }

  sodat <- load_engine_upload_dat(case_id, "05_spOccupancy")
  scfg <- runner_env$spocc_cfg()
  results <- add_result(results, senv$run_spocc_case(root, paste(case_id, "spocc_msPGOcc", sep = "__"), sodat, scfg, sodat$det,
                                                     species_table = sodat$species_table, notes = "maximum upload: multi-species replicated occupancy."), sodat)
  if (!smoke) {
    single_species <- colnames(sodat$Y)[1]
    single_y <- senv$make_single_detection(sodat, single_species)
    sp_tab <- data.frame(species = single_species, stringsAsFactors = FALSE)
    scfg2 <- scfg; scfg2$model$model_type <- "PGOcc"; scfg2$model$data_structure <- "single-species replicated"; scfg2$validation_prediction_outputs$waicOcc <- FALSE
    results <- add_result(results, senv$run_spocc_case(root, paste(case_id, "spocc_PGOcc_single", sep = "__"), sodat, scfg2, single_y,
                                                       species_table = sp_tab, notes = "maximum upload: single-species PGOcc."), sodat)
    scfg3 <- scfg2; scfg3$model$model_type <- "spPGOcc"; scfg3$spatial_latent_svc$NNGP <- TRUE
    results <- add_result(results, senv$run_spocc_case(root, paste(case_id, "spocc_spPGOcc_single", sep = "__"), sodat, scfg3, single_y,
                                                       species_table = sp_tab, notes = "maximum upload: single-species spatial PGOcc."), sodat)
    scfg4 <- scfg; scfg4$model$model_type <- "lfMsPGOcc"; scfg4$spatial_latent_svc$n.factors <- 2L
    results <- add_result(results, senv$run_spocc_case(root, paste(case_id, "spocc_lfMsPGOcc", sep = "__"), sodat, scfg4, sodat$det,
                                                       species_table = sodat$species_table, notes = "maximum upload: latent-factor multi-species branch."), sodat)
  }

  sjdat <- load_engine_upload_dat(case_id, "06_sjSDM")
  sj <- runner_env$sjsdm_cfg()
  results <- add_result(results, senv$run_sjsdm_case(root, paste(case_id, "sjsdm_binomial_linear_spatial", sep = "__"), sjdat, sj, sjdat$Y,
                                                     "maximum upload: binomial linear env + spatial."), sjdat)
  if (!smoke) {
    sj2 <- sj; sj2$model$family <- "gaussian_identity"; sj2$model$spatial_model <- "none"
    results <- add_result(results, senv$run_sjsdm_case(root, paste(case_id, "sjsdm_gaussian_no_spatial", sep = "__"), sjdat, sj2, sjdat$Y_normal,
                                                       "maximum upload: gaussian Y from common_inputs, no spatial."), sjdat)
    sj3 <- sj; sj3$model$env_model <- "DNN"; sj3$dnn_optimizer$hidden <- "4"; sj3$model$iter <- 3L; sj3$model$sampling <- 50L
    results <- add_result(results, senv$run_sjsdm_case(root, paste(case_id, "sjsdm_dnn_environment", sep = "__"), sjdat, sj3, sjdat$Y,
                                                       "maximum upload: DNN environment branch on CPU."), sjdat)
  }

  bdat <- load_engine_upload_dat(case_id, "07_boral")
  bcfg <- runner_env$boral_cfg()
  results <- add_result(results, senv$run_boral_case(root, paste(case_id, "boral_binomial", sep = "__"), bdat, bcfg, bdat$Y,
                                                     "maximum upload: binomial latent variable branch."), bdat)
  if (!smoke) {
    bcfg2 <- bcfg; bcfg2$model$family <- "poisson"; bcfg2$model$row.eff <- "fixed"
    results <- add_result(results, senv$run_boral_case(root, paste(case_id, "boral_poisson_row_effect", sep = "__"), bdat, bcfg2, bdat$Y_count,
                                                       "maximum upload: poisson Y from common_inputs, fixed row effects."), bdat)
    bcfg3 <- bcfg; bcfg3$model$family <- "normal"; bcfg3$model$use_offset <- TRUE; bcfg3$model$row.eff <- "none"
    results <- add_result(results, senv$run_boral_case(root, paste(case_id, "boral_normal_offset", sep = "__"), bdat, bcfg3, bdat$Y_normal,
                                                       "maximum upload: normal Y from common_inputs with offset file present."), bdat)
  }
  do.call(rbind, results)
}

write_report <- function(root, summary) {
  by_engine <- aggregate(case_id ~ engine + status, data = summary, FUN = length)
  names(by_engine)[names(by_engine) == "case_id"] <- "n_cases"
  write.csv(by_engine, file.path(root, "maximum_upload_status_by_engine.csv"), row.names = FALSE)
  write.csv(summary, file.path(root, "maximum_upload_case_summary.csv"), row.names = FALSE)
  bad <- summary[!(summary$status %in% c("fitted", "model_defined")) | summary$missing_count > 0, , drop = FALSE]
  lines <- c(
    "# JSDM Studio Maximum Upload Book Case Audit",
    "",
    paste0("Generated: ", Sys.time()),
    paste0("Audit root: ", normalizePath(root, winslash = "/", mustWork = FALSE)),
    "",
    "## Scope",
    "",
    "The suite uses examples/00_UPLOAD_READY_book_chapter_cases and reads each engine-specific upload folder.",
    "All available files for each engine are copied into each output ZIP under inputs/maximum_upload_files and data/maximum_upload_files.",
    "Mutually exclusive model options are tested as separate branches.",
    "",
    "## Status By Engine",
    "",
    paste(capture.output(print(by_engine, row.names = FALSE)), collapse = "\n"),
    "",
    "## Non-Fitted Or Incomplete Cases",
    "",
    if (nrow(bad)) paste(capture.output(print(bad[, c("engine", "case_id", "status", "missing_count", "failure_reason")], row.names = FALSE)), collapse = "\n") else "None.",
    "",
    "## Full Case Summary",
    "",
    paste(capture.output(print(summary[, c("engine", "case_id", "status", "file_count", "missing_count", "zip_size_kb", "maximum_upload_files_copied", "failure_reason")],
                               row.names = FALSE)), collapse = "\n")
  )
  writeLines(lines, file.path(root, "maximum_upload_audit_report.md"), useBytes = TRUE)
}

run_suite <- function() {
  cli <- parse_cli()
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  root <- file.path(app_dir, "output", paste0("book_maximum_upload_audit_", timestamp))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  all_results <- list()
  for (case_id in cli$cases) {
    if (!dir.exists(file.path(upload_root, case_id))) stop("Missing upload-ready case: ", case_id, call. = FALSE)
    all_results[[case_id]] <- run_one_case(root, case_id, smoke = cli$smoke)
  }
  summary <- do.call(rbind, all_results)
  rownames(summary) <- NULL
  write_report(root, summary)
  zipfile <- runner_env$make_zip_file(root, file.path(app_dir, "output", paste0("JSDMStudio_book_maximum_upload_audit_", timestamp, ".zip")))
  sha <- if (requireNamespace("digest", quietly = TRUE)) digest::digest(file = zipfile, algo = "sha256") else NA_character_
  writeLines(c(paste0("AUDIT_ROOT=", normalizePath(root, winslash = "/", mustWork = FALSE)),
               paste0("SUMMARY_CSV=", normalizePath(file.path(root, "maximum_upload_case_summary.csv"), winslash = "/", mustWork = FALSE)),
               paste0("AUDIT_ZIP=", normalizePath(zipfile, winslash = "/", mustWork = FALSE)),
               paste0("SHA256=", sha)),
             file.path(root, "maximum_upload_audit_paths.txt"), useBytes = TRUE)
  cat("AUDIT_ROOT=", normalizePath(root, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("SUMMARY_CSV=", normalizePath(file.path(root, "maximum_upload_case_summary.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("AUDIT_ZIP=", normalizePath(zipfile, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("SHA256=", sha, "\n", sep = "")
  invisible(summary)
}

if (!interactive()) run_suite()
