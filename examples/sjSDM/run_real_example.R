# Real sjSDM smoke test for JSDM Studio.
# Run from the JSDMStudio project root:
# Rscript examples/sjSDM/run_real_example.R

app_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
e <- new.env(parent = globalenv())
sys.source(file.path(app_dir, "app.R"), envir = e)

outdir <- e$make_engine_run_dir("sjSDM", "JSDMStudio_REAL_EXAMPLE")
for (fn in c("Y.csv", "env.csv", "newdata.csv")) {
  src <- file.path(app_dir, "examples", "sjSDM", fn)
  if (file.exists(src)) {
    file.copy(src, file.path(outdir, "data", fn), overwrite = TRUE)
    file.copy(src, file.path(outdir, "inputs", fn), overwrite = TRUE)
  }
}

cfg <- list(
  project_name = "JSDMStudio_REAL_EXAMPLE",
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
    hidden = "10,10,10",
    activation = "selu",
    dropout = 0,
    bias = TRUE,
    optimizer = "Adamax",
    learning_rate = 0.003,
    weight_decay = 0.002,
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
    spatial_ev_k = 20L,
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

yaml::write_yaml(cfg, file.path(outdir, "used_config.yml"))
e$write_sjsdm_reproducible_script(outdir)

script_file <- file.path(outdir, "reproducible_script", "run_this_sjSDM_analysis.R")
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
log <- system2(rscript, shQuote(script_file), stdout = TRUE, stderr = TRUE)
writeLines(as.character(log), file.path(outdir, "diagnostics", "example_stdout_stderr.txt"))
exit_status <- attr(log, "status") %||% 0

zipfile <- e$make_zip(outdir)
cat("Output folder:", outdir, "\n")
cat("ZIP:", zipfile, "\n")
cat("Exit status:", exit_status, "\n")
if (!identical(as.integer(exit_status), 0L)) {
  stop("Real sjSDM example failed. See diagnostics/example_stdout_stderr.txt in the output folder.", call. = FALSE)
}
