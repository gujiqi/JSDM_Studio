neutralize <- function() {
  if (.Platform$OS.type == "windows") {
    Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
  }
}
neutralize()

app_env <- new.env(parent = globalenv())
sys.source("app.R", envir = app_env)

Y <- app_env$read_csv_safe("examples/HMSC/Y.csv")
X <- app_env$read_csv_safe("examples/HMSC/XData.csv")
newdata <- app_env$read_csv_safe("examples/HMSC/newdata.csv")

cfg <- list(
  project_name = "PredictionSmoke",
  engine = "Hmsc",
  data = list(Y = "Y.csv", XData = "XData.csv", newdata = "newdata.csv"),
  model = list(
    distr = "probit", XFormula = "~ .", TrFormula = "~ .",
    use_traits = FALSE, use_phylogeny = FALSE,
    random_mode = "none", random_effect_column = "sample",
    spatial_method = "NNGP", nNeighbours = 5,
    lon_col = "longitude", lat_col = "latitude",
    seed = 1234, XScale = TRUE, TrScale = TRUE, YScale = FALSE,
    truncateNumberOfFactors = TRUE,
    Loff_file = "", ranLevelsUsed = "", C_file = "",
    use_XRRR = FALSE, ncRRR = 2, XRRRFormula = "~ .", XRRRScale = TRUE, XRRR_file = "",
    random_level_type = "none", sMethod = "NNGP", random_N = 1, longlat = FALSE,
    units_column = "sample", distMat_file = "", xData_file = "", sKnot_file = "",
    nfMin = 1, nfMax = 10,
    priors = list(setDefault = TRUE, a1 = NA, b1 = NA, a2 = NA, b2 = NA, alphapw = "")
  ),
  mcmc = list(
    samples = 3, transient = 1, thin = 1, nChains = 1, nParallel = 1,
    verbose = 0, preset = "Smoke", initPar = "fixed effects", alignPost = FALSE,
    updater = list(GammaEta = TRUE, Beta = TRUE, Gamma = TRUE, Omega = TRUE),
    sample_prior = FALSE, pool_chains = FALSE
  ),
  outputs = list(
    save_model = TRUE, predicted = TRUE, fit = FALSE, cv = FALSE, waic = FALSE,
    diagnostics = FALSE, parameters = FALSE, variance_partitioning = FALSE, omega = FALSE,
    gradients = FALSE, beta_support = FALSE, gamma_support = FALSE, omega_support = FALSE,
    convergence = list(showBeta = FALSE, showGamma = FALSE, showOmega = FALSE, maxOmega = 10,
                       showRho = FALSE, showAlpha = FALSE, effectiveSize = FALSE, gelmanPSRF = FALSE),
    plotting = list(var.part.order.explained = TRUE, var.part.order.raw = FALSE,
                    show.sp.names.beta = FALSE, plotTree = FALSE, omega.order = "original",
                    show.sp.names.omega = TRUE, plotBeta = FALSE, plotGamma = FALSE),
    predictions = list(species.list = "", trait.list = "", env.list = "", nfolds = 2,
                       partition_column = "sample", computeSAIR = FALSE)
  )
)

outdir <- app_env$make_engine_run_dir("Hmsc_prediction_smoke", "PredictionSmoke")
res <- app_env$run_hmsc_s1s7_pipeline(
  outdir = outdir, cfg = cfg, Y = Y, XData = X, newdata = newdata,
  log_fun = function(txt) message(txt)
)

std <- file.path(outdir, "standard", "newdata_predictions_long.csv")
pred <- file.path(outdir, "predictions", "newdata_predicted_mean_main.csv")
manifest <- file.path(outdir, "predictions", "newdata_prediction_manifest.csv")

stopifnot(identical(res$status, "fitted"))
stopifnot(file.exists(std), file.exists(pred), file.exists(manifest))
tab <- read.csv(std, check.names = FALSE)
stopifnot(nrow(tab) == nrow(newdata) * ncol(Y))
stopifnot(all(tab$prediction_context == "newdata"))
zipfile <- app_env$make_zip(outdir)
stopifnot(file.exists(zipfile), file.info(zipfile)$size > 0)

cat("HMSC_NEWDATA_PREDICTION_SMOKE_OK\n")
cat("OUTDIR=", normalizePath(outdir, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("NEWDATA_STANDARD_ROWS=", nrow(tab), "\n", sep = "")
cat("NEWDATA_PRED_MEAN=", normalizePath(pred, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("ZIP=", normalizePath(zipfile, winslash = "/", mustWork = TRUE), "\n", sep = "")
