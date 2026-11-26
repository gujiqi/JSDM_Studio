neutralize <- function() {
  if (.Platform$OS.type == "windows") {
    Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
  }
}
neutralize()

app_env <- new.env(parent = globalenv())
sys.source("app.R", envir = app_env)

case_dir <- file.path("examples", "00_UPLOAD_READY_book_chapter_cases", "ch11_finnish_birds", "00_common_inputs")
Y <- app_env$read_csv_safe(file.path(case_dir, "Y_occurrence.csv"))
X <- app_env$read_csv_safe(file.path(case_dir, "XData.csv"))
study <- app_env$read_csv_safe(file.path(case_dir, "studyDesign.csv"))
coords <- app_env$read_csv_safe(file.path(case_dir, "coordinates.csv"))
newdata <- app_env$read_csv_safe(file.path(case_dir, "newdata.csv"))
newcoords <- app_env$read_csv_safe(file.path(case_dir, "newcoords.csv"))

cfg <- list(
  project_name = "BookCh11_HMSC_spatial_newcoords_prediction",
  engine = "Hmsc",
  data = list(
    Y = "Y_occurrence.csv",
    XData = "XData.csv",
    studyDesign = "studyDesign.csv",
    coordinates = "coordinates.csv",
    newdata = "newdata.csv",
    newcoords = "newcoords.csv"
  ),
  model = list(
    distr = "probit",
    XFormula = "~ pH + moisture + canopy + elevation + substrate",
    TrFormula = "~ .",
    use_traits = FALSE,
    use_phylogeny = FALSE,
    random_mode = "spatial_full",
    random_effect_column = "sample",
    spatial_method = "Full",
    nNeighbours = 5,
    lon_col = "x",
    lat_col = "y",
    longlat = FALSE,
    seed = 1234,
    XScale = TRUE,
    TrScale = TRUE,
    YScale = FALSE,
    truncateNumberOfFactors = TRUE,
    Loff_file = "",
    ranLevelsUsed = "",
    C_file = "",
    use_XRRR = FALSE,
    ncRRR = 2,
    XRRRFormula = "~ .",
    XRRRScale = TRUE,
    XRRR_file = "",
    random_level_type = "none",
    sMethod = "NNGP",
    random_N = 1,
    units_column = "sample",
    distMat_file = "",
    xData_file = "",
    sKnot_file = "",
    nfMin = 1,
    nfMax = 10,
    priors = list(setDefault = TRUE, a1 = NA, b1 = NA, a2 = NA, b2 = NA, alphapw = "")
  ),
  mcmc = list(
    samples = 5,
    transient = 10,
    thin = 1,
    nChains = 1,
    nParallel = 1,
    verbose = 0,
    preset = "Book chapter 11 smoke",
    initPar = "fixed effects",
    alignPost = FALSE,
    updater = list(GammaEta = TRUE, Beta = TRUE, Gamma = TRUE, Omega = TRUE),
    sample_prior = FALSE,
    pool_chains = FALSE
  ),
  outputs = list(
    save_model = TRUE,
    predicted = TRUE,
    fit = FALSE,
    cv = FALSE,
    waic = FALSE,
    diagnostics = FALSE,
    parameters = FALSE,
    variance_partitioning = FALSE,
    omega = FALSE,
    gradients = FALSE,
    beta_support = FALSE,
    gamma_support = FALSE,
    omega_support = FALSE,
    convergence = list(showBeta = FALSE, showGamma = FALSE, showOmega = FALSE, maxOmega = 10,
                       showRho = FALSE, showAlpha = FALSE, effectiveSize = FALSE, gelmanPSRF = FALSE),
    plotting = list(var.part.order.explained = TRUE, var.part.order.raw = FALSE,
                    show.sp.names.beta = FALSE, plotTree = FALSE, omega.order = "original",
                    show.sp.names.omega = TRUE, plotBeta = FALSE, plotGamma = FALSE),
    predictions = list(species.list = "", trait.list = "", env.list = "", nfolds = 2,
                       partition_column = "sample", computeSAIR = FALSE)
  )
)

expected_rows <- nrow(newdata) * ncol(Y)
outdir <- app_env$make_engine_run_dir("Hmsc_book_ch11_prediction", "BookCh11SpatialPrediction")
res <- app_env$run_hmsc_s1s7_pipeline(
  outdir = outdir,
  cfg = cfg,
  Y = Y,
  XData = X,
  studyDesign = study,
  coord = coords,
  newdata = newdata,
  newcoord = newcoords,
  log_fun = function(txt) message("[book_ch11_spatial_prediction] ", txt)
)

std_file <- file.path(outdir, "standard", "newdata_predictions_long.csv")
mean_file <- file.path(outdir, "predictions", "newdata_predicted_mean_main.csv")
manifest_file <- file.path(outdir, "predictions", "newdata_prediction_manifest.csv")
zipfile <- app_env$make_zip(outdir)

stopifnot(identical(res$status, "fitted"))
stopifnot(file.exists(std_file))
stopifnot(file.exists(mean_file))
stopifnot(file.exists(manifest_file))
stopifnot(file.exists(zipfile), file.info(zipfile)$size > 0)

std <- read.csv(std_file, check.names = FALSE)
pred_mean <- read.csv(mean_file, check.names = FALSE)
manifest <- read.csv(manifest_file, check.names = FALSE)

stopifnot(nrow(std) == expected_rows)
stopifnot(all(std$prediction_context == "newdata"))
stopifnot(nrow(pred_mean) == nrow(newdata))

summary <- data.frame(
  case = "ch11_finnish_birds",
  engine = "Hmsc",
  random_mode = "spatial_full",
  status = res$status,
  training_sites = nrow(Y),
  species = ncol(Y),
  prediction_sites = nrow(newdata),
  expected_prediction_rows = expected_rows,
  actual_prediction_rows = nrow(std),
  manifest_status = paste(unique(manifest$status), collapse = ";"),
  outdir = normalizePath(outdir, winslash = "/", mustWork = TRUE),
  zip = normalizePath(zipfile, winslash = "/", mustWork = TRUE),
  stringsAsFactors = FALSE
)
write.csv(summary, file.path("output", "latest_hmsc_ch11_spatial_prediction_summary.csv"), row.names = FALSE)

cat("HMSC_CH11_SPATIAL_NEWDATA_PREDICTION_OK\n")
print(summary)
cat("\nFirst prediction rows:\n")
print(utils::head(std, 10))
