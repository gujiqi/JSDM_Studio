# Real Hmsc-HPC / pyhmsc workflow suite for JSDM Studio.
# It creates small synthetic ecological datasets and runs the same adapter used by the Shiny workflow.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (length(x) == 1 && is.atomic(x) && is.na(x)) return(y)
  x
}

script_path <- {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- args[startsWith(args, file_arg)]
  if (length(hit) > 0) normalizePath(sub(file_arg, "", hit[[1]]), winslash = "/", mustWork = FALSE) else
    normalizePath("examples/HmscHPC/run_real_example_suite.R", winslash = "/", mustWork = FALSE)
}
example_dir <- dirname(script_path)
app_dir <- normalizePath(file.path(example_dir, "..", ".."), winslash = "/", mustWork = TRUE)
setwd(app_dir)

source(file.path(app_dir, "R", "hmschpc_adapter.R"))

timestamp_id <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

engine_output_dirs <- function(key) {
  c("inputs", "data", "models", "tables", "results", "plots", "predictions",
    "diagnostics", "workflow_scripts", "reproducible_script", "standard",
    "report", "samples", "compile", "hdf5", "python_logs")
}

make_dirs <- function(outdir) {
  for (d in engine_output_dirs("hmschpc")) {
    dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
  }
}

write_engine_status <- function(outdir, status) {
  dir.create(file.path(outdir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
  status$timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(status, file.path(outdir, "diagnostics", "engine_status.json"),
                         auto_unbox = TRUE, pretty = TRUE)
  }
  invisible(status)
}

make_zip <- function(outdir) {
  zipfile <- paste0(normalizePath(outdir, winslash = "/", mustWork = TRUE), ".zip")
  if (file.exists(zipfile)) unlink(zipfile)
  files <- list.files(outdir, recursive = TRUE, all.files = FALSE, full.names = FALSE)
  if (!length(files)) stop("ZIP creation failed because output folder has no files: ", outdir, call. = FALSE)
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(zipfile, files = files, root = outdir, recurse = FALSE, include_directories = FALSE, mode = "mirror")
  } else {
    old <- setwd(outdir)
    on.exit(setwd(old), add = TRUE)
    utils::zip(zipfile, files = files, flags = "-r9X")
  }
  normalizePath(zipfile, winslash = "/", mustWork = FALSE)
}

read_status <- function(outdir) {
  path <- file.path(outdir, "diagnostics", "engine_status.json")
  if (file.exists(path) && requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::fromJSON(path, simplifyVector = FALSE)
  } else {
    list(status = "missing_status")
  }
}

empty_dirs <- function(outdir) {
  dirs <- list.dirs(outdir, recursive = TRUE, full.names = TRUE)
  rel <- character()
  for (d in dirs) {
    if (identical(normalizePath(d, winslash = "/", mustWork = FALSE),
                  normalizePath(outdir, winslash = "/", mustWork = FALSE))) next
    if (length(list.files(d, all.files = FALSE, no.. = TRUE)) == 0) {
      rel <- c(rel, gsub("\\\\", "/", sub(paste0("^", normalizePath(outdir, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(d, winslash = "/", mustWork = FALSE))))
    }
  }
  rel
}

required_files_for_status <- function(status, predictions = TRUE) {
  base <- c(
    "used_config.yml",
    "workflow_scripts/hmschpc_model.yaml",
    "workflow_scripts/run_HmscHPC_workflow.R",
    "reproducible_script/run_this_HmscHPC_analysis.R",
    "reproducible_script/run_this_HmscHPC_analysis.py",
    "diagnostics/engine_status.json",
    "standard/run_summary.csv",
    "standard/effects_long.csv",
    "standard/predictions_long.csv",
    "standard/fit_metrics.csv",
    "standard/associations_long.csv",
    "standard/diagnostics_long.csv",
    "standard/output_manifest.csv",
    "tables/Hmsc-HPC_result_workflow_map.csv",
    "tables/HmscHPC_S1S7_result_index.csv",
    "workflow_scripts/S1_define_models.py",
    "workflow_scripts/S2_fit_models.py",
    "workflow_scripts/S3_evaluate_convergence.py",
    "workflow_scripts/S4_compute_model_fit.py",
    "workflow_scripts/S5_show_model_fit.py",
    "workflow_scripts/S6_show_parameter_estimates.py",
    "workflow_scripts/S7_make_predictions.py",
    "results/S1_model_definition.csv",
    "results/S2_fit_models.csv",
    "results/S3_convergence_summary.csv",
    "results/S4_model_fit_summary.csv",
    "results/S5_model_fit_prediction_summary.csv",
    "results/S6_parameter_estimates_Beta.csv",
    "results/S7_predictions_training.csv",
    "diagnostics/data_check_messages.csv",
    "diagnostics/session_info.txt",
    "report/Hmsc-HPC_report.html",
    "results/README_HmscHPC_results.txt"
  )
  if (identical(status, "fitted")) {
    out <- c(base,
             "models/compiled_model/init.json",
             "samples/posterior.h5",
             "tables/Beta_summary.csv")
    if (isTRUE(predictions)) out <- c(out, "predictions/predicted_mean.csv")
    out
  } else if (identical(status, "model_defined")) {
    c(base, "models/compiled_model/init.json")
  } else {
    base
  }
}

make_synthetic <- function(seed = 1, distribution = "poisson", n = 9, s = 3) {
  set.seed(seed)
  site_id <- paste0("site_", seq_len(n))
  sp_id <- paste0("sp_", seq_len(s))
  substrate <- rep(c("sand", "clay", "peat"), length.out = n)
  forest_cover <- round(runif(n, 0.15, 0.9), 3)
  elevation <- round(as.numeric(scale(seq_len(n) + rnorm(n, sd = 0.4))), 3)
  moisture <- round(runif(n, 0.1, 1), 3)
  X <- data.frame(forest_cover = forest_cover, elevation = elevation,
                  moisture = moisture, substrate = substrate,
                  stringsAsFactors = FALSE)
  rownames(X) <- site_id
  sp_offset <- seq(-0.35, 0.35, length.out = s)
  eta <- outer(0.4 + 0.9 * forest_cover - 0.45 * elevation + 0.25 * moisture +
                 ifelse(substrate == "peat", 0.35, ifelse(substrate == "clay", -0.15, 0)),
               rep(1, s)) + matrix(rep(sp_offset, each = n), nrow = n)
  if (distribution == "poisson") {
    Y <- matrix(rpois(n * s, lambda = pmax(exp(eta), 0.05)), nrow = n)
  } else if (distribution == "probit") {
    Y <- matrix(rbinom(n * s, size = 1, prob = pmin(pmax(plogis(eta), 0.02), 0.98)), nrow = n)
  } else {
    Y <- eta + matrix(rnorm(n * s, sd = 0.35), nrow = n)
  }
  Y <- as.data.frame(Y, check.names = FALSE)
  colnames(Y) <- sp_id
  rownames(Y) <- site_id
  traits <- data.frame(
    body_size = round(seq(0.5, 1.5, length.out = s), 3),
    forest_specialist = rep(c(0, 1), length.out = s),
    life_form = rep(c("herb", "shrub", "tree"), length.out = s),
    stringsAsFactors = FALSE
  )
  rownames(traits) <- sp_id
  phylo_cov <- matrix(0.25, nrow = s, ncol = s, dimnames = list(sp_id, sp_id))
  diag(phylo_cov) <- 1
  study <- data.frame(
    plot = rep(paste0("plot_", seq_len(ceiling(n / 3))), each = 3, length.out = n),
    block = rep(c("north", "south", "east"), length.out = n),
    xcoord = round(seq(0, 1, length.out = n) + rnorm(n, sd = 0.02), 3),
    ycoord = round(runif(n, 0, 1), 3),
    stringsAsFactors = FALSE
  )
  rownames(study) <- site_id
  coord <- study[, c("xcoord", "ycoord")]
  newdata <- X[seq_len(min(3, n)), , drop = FALSE]
  rownames(newdata) <- paste0("new_site_", seq_len(nrow(newdata)))
  list(Y = Y, X = X, traits = traits, phylo_cov = phylo_cov, study = study,
       coord = coord, newdata = newdata,
       tree = "((sp_1:0.4,sp_2:0.4):0.6,sp_3:1.0);")
}

base_config <- function(distribution = "poisson", XFormula = "~ forest_cover + elevation + moisture",
                        use_traits = FALSE, trait_formula = "~ body_size + forest_specialist",
                        phylogeny_mode = "none", random_mode = "none",
                        random_slope_formula = "", run_sampler = TRUE,
                        random_name = "plot", nf = 1, nfMin = 1, nfMax = 4, alpha = 1,
                        chains = 1, chains_to_run = NULL,
                        samples = 3, transient = 3, thin = 1, verbose = 3,
                        seed = 1234, precision = 64, tnlib = "tf",
                        hmcleapfrog = 5, hmcthin = 0,
                        update_beta_eta = FALSE, save_eta = TRUE,
                        eager = FALSE, profile = FALSE,
                        predictions = TRUE, diagnostics = TRUE, plots = TRUE) {
  list(
    engine = "Hmsc-HPC",
    inputs = list(Y = "synthetic_Y.csv", XData = "synthetic_XData.csv",
                  traits = if (use_traits) "synthetic_traits.csv" else NULL,
                  studyDesign = if (!identical(random_mode, "none")) "synthetic_studyDesign.csv" else NULL,
                  coordinates = if (identical(random_mode, "spatial_full")) "synthetic_coordinates.csv" else NULL,
                  phylo_cov = if (identical(phylogeny_mode, "covariance")) "synthetic_phylo_cov.csv" else NULL,
                  phylo_tree = if (identical(phylogeny_mode, "newick")) "synthetic_tree.nwk" else NULL,
                  newdata = "synthetic_newdata.csv"),
    model = list(distribution = distribution, XFormula = XFormula,
                 use_traits = use_traits, trait_formula = trait_formula,
                 phylogeny_mode = phylogeny_mode),
    random_effects = list(mode = random_mode, name = random_name, column = "plot",
                          coord_x = "xcoord", coord_y = "ycoord",
                          x_formula = random_slope_formula,
                          nf = nf, nfMin = nfMin, nfMax = nfMax, alpha = alpha),
    sampler = list(run_sampler = run_sampler, samples = samples, transient = transient,
                   thin = thin, chains = chains, chains_to_run = chains_to_run,
                   verbose = verbose, seed = seed, precision = precision,
                   truncated_normal_library = tnlib, hmcleapfrog = hmcleapfrog,
                   hmcthin = hmcthin, update_beta_eta = update_beta_eta,
                   save_eta = save_eta, eager = eager, profile = profile),
    outputs = list(predictions = predictions, diagnostics = diagnostics, plots = plots, zip = TRUE),
    runtime = list(python = hmschpc_python_bin(list()),
                   python_source = hmschpc_python_source_dir())
  )
}

case_plan <- list(
  list(id = "01_fixed_poisson_categorical", distribution = "poisson",
       cfg = base_config(distribution = "poisson",
                         XFormula = "~ forest_cover + elevation + moisture + C(substrate)",
                         chains = 2, chains_to_run = c(0, 1), samples = 3, transient = 3,
                         seed = 101, tnlib = "tf")),
  list(id = "02_fixed_probit_binary", distribution = "probit",
       cfg = base_config(distribution = "probit", XFormula = "~ forest_cover + elevation",
                         samples = 3, transient = 3, thin = 1, seed = 102,
                         precision = 32, tnlib = "tf")),
  list(id = "03_fixed_normal_hmc_toggle", distribution = "normal",
       cfg = base_config(distribution = "normal", XFormula = "~ forest_cover + elevation + moisture",
                         samples = 3, transient = 4, thin = 2, seed = 103,
                         hmcthin = 1, diagnostics = FALSE, plots = FALSE)),
  list(id = "04_traits_phylo_cov", distribution = "poisson",
       cfg = base_config(distribution = "poisson", use_traits = TRUE,
                         trait_formula = "~ body_size + forest_specialist + C(life_form)",
                         phylogeny_mode = "covariance", samples = 3, transient = 3,
                         seed = 104, tnlib = "scipy")),
  list(id = "05_iid_random_intercept", distribution = "poisson",
       cfg = base_config(distribution = "poisson", random_mode = "iid",
                         random_name = "plot_iid",
                         samples = 3, transient = 3, seed = 105,
                         update_beta_eta = TRUE, save_eta = TRUE)),
  list(id = "06_spatial_full_random_intercept", distribution = "poisson",
       cfg = base_config(distribution = "poisson", random_mode = "spatial_full",
                         random_name = "spatial_plot", alpha = 0.5,
                         samples = 3, transient = 3, seed = 106,
                         hmcleapfrog = 6, save_eta = TRUE)),
  list(id = "07_phylo_newick_fixed", distribution = "poisson",
       cfg = base_config(distribution = "poisson", phylogeny_mode = "newick",
                         samples = 3, transient = 3, seed = 107, tnlib = "tf")),
  list(id = "08_random_slope_compile_only", distribution = "poisson",
       cfg = base_config(distribution = "poisson", random_mode = "random_slope_iid",
                         random_slope_formula = "~ xcoord + ycoord",
                         run_sampler = FALSE, chains = 1, samples = 3, transient = 3,
                         seed = 108)),
  list(id = "09_gaussian_alias_tfd_predictions_off", distribution = "normal",
       cfg = base_config(distribution = "gaussian",
                         XFormula = "~ forest_cover + elevation + moisture + C(substrate)",
                         samples = 1, transient = 1, thin = 1, seed = 109,
                         precision = 64, tnlib = "tfd", hmcleapfrog = 3,
                         hmcthin = 0, save_eta = FALSE, eager = TRUE,
                         profile = TRUE, predictions = FALSE))
)

run_case <- function(case) {
  dat <- make_synthetic(seed = as.integer(gsub("\\D", "", case$id)) + 100,
                        distribution = case$distribution, n = 9, s = 3)
  outdir <- file.path(app_dir, "output", paste0("HmscHPC_real_", case$id, "_", timestamp_id()))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  make_dirs(outdir)
  tree_file <- NULL
  if (identical(case$cfg$model$phylogeny_mode, "newick")) {
    tree_file <- file.path(outdir, "synthetic_tree.nwk")
    writeLines(dat$tree, tree_file)
  }
  check <- validate_hmschpc_full(
    dat$Y, dat$X, dat$traits, dat$study, dat$coord,
    dat$phylo_cov, tree_file, dat$newdata,
    distribution = case$cfg$model$distribution,
    XFormula = case$cfg$model$XFormula,
    use_traits = isTRUE(case$cfg$model$use_traits),
    trait_formula = case$cfg$model$trait_formula,
    phylogeny_mode = case$cfg$model$phylogeny_mode,
    random_mode = case$cfg$random_effects$mode,
    random_name = case$cfg$random_effects$name,
    random_column = case$cfg$random_effects$column,
    coord_x = case$cfg$random_effects$coord_x,
    coord_y = case$cfg$random_effects$coord_y,
    nf = case$cfg$random_effects$nf,
    nfMin = case$cfg$random_effects$nfMin,
    nfMax = case$cfg$random_effects$nfMax,
    alpha = case$cfg$random_effects$alpha,
    samples = case$cfg$sampler$samples,
    transient = case$cfg$sampler$transient,
    thin = case$cfg$sampler$thin,
    chains = case$cfg$sampler$chains,
    chains_to_run = case$cfg$sampler$chains_to_run %||% integer(),
    verbose = case$cfg$sampler$verbose,
    seed = case$cfg$sampler$seed,
    precision = case$cfg$sampler$precision,
    truncated_normal_library = case$cfg$sampler$truncated_normal_library,
    hmcleapfrog = case$cfg$sampler$hmcleapfrog,
    hmcthin = case$cfg$sampler$hmcthin,
    run_sampler = case$cfg$sampler$run_sampler,
    random_slope_formula = case$cfg$random_effects$x_formula,
    python = case$cfg$runtime$python,
    python_source = case$cfg$runtime$python_source,
    predictions = case$cfg$outputs$predictions,
    diagnostics = case$cfg$outputs$diagnostics,
    plots = case$cfg$outputs$plots
  )
  if (!isTRUE(check$ok)) {
    writeLines(check$messages, file.path(outdir, "diagnostics", "data_check_messages.txt"))
    write_engine_status(outdir, list(engine = "Hmsc-HPC", status = "check_failed",
                                     errors = check$messages, warnings = character()))
    zipfile <- make_zip(outdir)
    return(data.frame(case = case$id, status = "check_failed", output_folder = outdir,
                      zip = zipfile, zip_size_kb = round(file.info(zipfile)$size / 1024, 1),
                      file_count = length(list.files(outdir, recursive = TRUE)),
                      missing_required = paste(required_files_for_status("check_failed")[!file.exists(file.path(outdir, required_files_for_status("check_failed")))], collapse = ";"),
                      empty_dirs = paste(empty_dirs(outdir), collapse = ";"),
                      stringsAsFactors = FALSE))
  }
  status <- run_hmschpc_workflow(
    outdir, case$cfg, dat$Y, dat$X,
    Tr = if (isTRUE(case$cfg$model$use_traits)) dat$traits else NULL,
    study = if (!identical(case$cfg$random_effects$mode, "none")) dat$study else NULL,
    coord = if (identical(case$cfg$random_effects$mode, "spatial_full")) dat$coord else NULL,
    phylo_cov = if (identical(case$cfg$model$phylogeny_mode, "covariance")) dat$phylo_cov else NULL,
    phylo_tree_file = if (identical(case$cfg$model$phylogeny_mode, "newick")) tree_file else NULL,
    newdata = dat$newdata,
    log_fun = function(txt) message("[", case$id, "] ", txt)
  )
  actual <- read_status(outdir)$status %||% status$status %||% "unknown"
  zipfile <- make_zip(outdir)
  req <- required_files_for_status(actual, predictions = isTRUE(case$cfg$outputs$predictions))
  miss <- req[!file.exists(file.path(outdir, req))]
  data.frame(case = case$id, status = actual, output_folder = normalizePath(outdir, winslash = "/", mustWork = FALSE),
             zip = zipfile, zip_size_kb = round(file.info(zipfile)$size / 1024, 1),
             file_count = length(list.files(outdir, recursive = TRUE)),
             missing_required = paste(miss, collapse = ";"),
             empty_dirs = paste(empty_dirs(outdir), collapse = ";"),
             stringsAsFactors = FALSE)
}

results <- do.call(rbind, lapply(case_plan, run_case))
summary_path <- file.path(app_dir, "output", paste0("HmscHPC_real_example_suite_summary_", timestamp_id(), ".csv"))
write.csv(results, summary_path, row.names = FALSE)
write.csv(results, file.path(example_dir, "last_real_example_suite_summary.csv"), row.names = FALSE)

audit_path <- file.path(example_dir, "HmscHPC_WORKFLOW_AUDIT.md")
lines <- c(
  "# Hmsc-HPC workflow audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "The suite runs synthetic Hmsc-HPC cases through the JSDM Studio R adapter and Python runner.",
  "Expected successful statuses are `fitted` for sampler cases and `model_defined` for the guarded random-slope compile-only case.",
  "",
  "## Case summary",
  "",
  paste(capture.output(print(results[, c("case", "status", "zip_size_kb", "file_count", "missing_required", "empty_dirs")], row.names = FALSE)), collapse = "\n"),
  "",
  "## Diagnostics",
  "",
  "For any failed case, open `diagnostics/engine_status.json` first, then `diagnostics/HmscHPC_error.txt`, `diagnostics/HmscHPC_compile.log`, `diagnostics/HmscHPC_validate_init.log` and `diagnostics/HmscHPC_sample.log` in that output folder.",
  "",
  paste0("Machine summary CSV: ", normalizePath(summary_path, winslash = "/", mustWork = FALSE))
)
writeLines(lines, audit_path)

print(results[, c("case", "status", "zip_size_kb", "file_count", "missing_required", "empty_dirs")], row.names = FALSE)
if (any(!results$status %in% c("fitted", "model_defined"))) {
  stop("One or more Hmsc-HPC cases failed. See ", summary_path, call. = FALSE)
}
if (any(nzchar(results$missing_required))) {
  stop("One or more Hmsc-HPC cases missed required output files. See ", summary_path, call. = FALSE)
}
message("Hmsc-HPC real example suite completed: ", summary_path)
