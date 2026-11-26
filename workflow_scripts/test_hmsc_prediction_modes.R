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
study <- app_env$read_csv_safe("examples/HMSC/studyDesign.csv")
coords <- app_env$read_csv_safe("examples/HMSC/coordinates.csv")
newdata <- app_env$read_csv_safe("examples/HMSC/newdata.csv")
newstudy <- app_env$read_csv_safe("examples/HMSC/prediction_studyDesign.csv")
newcoords <- app_env$read_csv_safe("examples/HMSC/newcoords.csv")

base_cfg <- function(random_mode) {
  list(
    project_name = paste0("PredictionMode_", random_mode),
    engine = "Hmsc",
    data = list(Y = "Y.csv", XData = "XData.csv", newdata = "newdata.csv",
                prediction_studyDesign = "prediction_studyDesign.csv", newcoords = "newcoords.csv"),
    model = list(
      distr = "probit", XFormula = "~ .", TrFormula = "~ .",
      use_traits = FALSE, use_phylogeny = FALSE,
      random_mode = random_mode, random_effect_column = "plot",
      spatial_method = "Full", nNeighbours = 5,
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
      samples = 5, transient = 10, thin = 1, nChains = 1, nParallel = 1,
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
}

cases <- list(
  no_random = list(mode = "none", study = NULL, coords = NULL, newstudy = NULL, newcoords = NULL),
  sample_grouping = list(mode = "sample", study = study, coords = NULL, newstudy = newstudy, newcoords = NULL),
  spatial_full = list(mode = "spatial_full", study = study, coords = coords, newstudy = NULL, newcoords = newcoords)
)

summary_rows <- list()
for (nm in names(cases)) {
  case <- cases[[nm]]
  cfg <- base_cfg(case$mode)
  outdir <- app_env$make_engine_run_dir("Hmsc_prediction_mode", paste0("PredictionMode_", nm))
  res <- app_env$run_hmsc_s1s7_pipeline(
    outdir = outdir, cfg = cfg, Y = Y, XData = X,
    studyDesign = case$study, coord = case$coords,
    newdata = newdata, newStudyDesign = case$newstudy, newcoord = case$newcoords,
    log_fun = function(txt) message("[", nm, "] ", txt)
  )
  std <- file.path(outdir, "standard", "newdata_predictions_long.csv")
  zipfile <- app_env$make_zip(outdir)
  rows <- if (file.exists(std)) nrow(read.csv(std, check.names = FALSE)) else 0L
  stopifnot(identical(res$status, "fitted"))
  stopifnot(rows == nrow(newdata) * ncol(Y))
  stopifnot(file.exists(zipfile), file.info(zipfile)$size > 0)
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    case = nm, random_mode = case$mode, status = res$status,
    newdata_prediction_rows = rows,
    outdir = normalizePath(outdir, winslash = "/", mustWork = TRUE),
    zip = normalizePath(zipfile, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, summary_rows)
write.csv(summary, file.path("output", "latest_hmsc_prediction_modes_summary.csv"), row.names = FALSE)
cat("HMSC_PREDICTION_MODES_OK\n")
print(summary)
