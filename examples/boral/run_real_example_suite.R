# Seven-case boral workflow smoke suite for JSDM Studio.
# It exercises data upload branches and writes real output folders/ZIPs.
# If system JAGS/boral is not available, cases are expected to end as fit_failed
# with complete diagnostics rather than empty or fake output.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_arg <- args[startsWith(args, file_arg)]
this_file <- if (length(script_arg) > 0) normalizePath(sub(file_arg, "", script_arg[[1]]), winslash = "/", mustWork = TRUE) else normalizePath("examples/boral/run_real_example_suite.R", winslash = "/", mustWork = FALSE)
example_dir <- dirname(this_file)
app_dir <- normalizePath(file.path(example_dir, "..", ".."), winslash = "/", mustWork = TRUE)
setwd(app_dir)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (length(x) == 1 && is.atomic(x) && is.na(x)) return(y)
  x
}
source(file.path(app_dir, "R", "boral_adapter.R"), local = FALSE)

write_csv <- function(x, path, row_names = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = row_names)
}

make_dirs <- function(outdir) {
  dirs <- c("", "inputs", "data", "models", "results", "tables", "plots", "diagnostics",
            "report", "predictions", "reproducible_script", "standard", "workflow_scripts",
            "jags", "mcmc", "ordination", "residuals", "random_effects", "variable_selection")
  for (d in dirs) dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
}

make_zip <- function(outdir) {
  zipfile <- paste0(outdir, ".zip")
  if (file.exists(zipfile)) unlink(zipfile)
  rel_files <- list.files(outdir, recursive = TRUE, all.files = FALSE, no.. = TRUE)
  if (!length(rel_files)) stop("ZIP creation failed because output folder has no files: ", outdir, call. = FALSE)
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(zipfile = zipfile, files = rel_files, recurse = FALSE, root = outdir, mode = "mirror")
  } else {
    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(outdir)
    utils::zip(zipfile, files = rel_files, flags = "-r9Xq")
  }
  if (!file.exists(zipfile) || file.info(zipfile)$size <= 0) stop("ZIP creation failed for ", outdir, call. = FALSE)
  zipfile
}

write_check_messages <- function(outdir, messages) {
  write.csv(data.frame(message = messages, stringsAsFactors = FALSE),
            file.path(outdir, "diagnostics", "data_check_messages.csv"), row.names = FALSE)
}

base_config <- function(name) {
  list(
    project_name = paste0("boral_suite_", name),
    engine = "boral",
    data = list(Y = "Y.csv", XData = "XData.csv", traits = "traits.csv",
                row.ids = "row.ids.csv", ranef.ids = "ranef.ids.csv",
                distmat = "distmat.csv", offset = "offset.csv",
                newdata = "newdata.csv", trial.size.file = "trial.size.csv"),
    model = list(model_mode = "correlated response GLMs", family = "poisson",
                 family_vector_override = "", lv.control = list(num.lv = 1, type = "independent"),
                 model.name = "jagsboralmodel.txt", formula.X = "~ .", X.ind = "",
                 trial.size = 1, row.eff = "none", use_offset = FALSE, do.fit = TRUE),
    traits_random_ssvs = list(use_traits = FALSE, which.traits = "", traits_no_intercept = TRUE,
                              use_ranef = FALSE, use_ssvs = FALSE, ssvs.index = "",
                              ssvs.traitsindex = "", ssvs.g = 1e-6, save.model = FALSE),
    mcmc_prior = list(n.burnin = 50, n.iteration = 120, n.thin = 5, seed = 1234,
                      prior.type = "normal,normal,normal,uniform", hypparams = "10,10,10,30",
                      calc.ics = FALSE, save_mcmc_samples = TRUE, save_hpd = TRUE, save_dic = TRUE),
    diagnostics_outputs = list(summary = TRUE, residual_plot = TRUE, lvsplot = TRUE,
                               ind.spp = "default", ranefsplot = FALSE, coefsplot = TRUE,
                               enviro_cor = TRUE, residual_cor = TRUE, varpart = FALSE,
                               predict = TRUE, fitted = TRUE, tidyboral = TRUE),
    outputs = list(save_model = TRUE, save_jags = TRUE, save_tables = TRUE,
                   save_plots = TRUE, report = TRUE, zip = TRUE,
                   copy_inputs = TRUE, save_config = TRUE)
  )
}

set.seed(42)
n <- 10
S <- 4
site_names <- paste0("site_", seq_len(n))
sp_names <- paste0("sp_", seq_len(S))
X <- data.frame(pH = round(runif(n, 5.5, 7.5), 2),
                moisture = round(runif(n), 3),
                substrate = rep(c("sand", "clay"), length.out = n),
                canopy = round(runif(n, 0.2, 0.9), 3),
                row.names = site_names)
newdata <- X[1:3, , drop = FALSE]
rowids <- data.frame(site_group = rep(1:3, length.out = n), year = rep(c("y1", "y2"), length.out = n), row.names = site_names)
ranefids <- data.frame(observer = rep(c("a", "b", "c"), length.out = n), row.names = site_names)
traits <- data.frame(body_size = round(runif(S, 0.2, 1.4), 2),
                     life_form = rep(c("annual", "perennial"), length.out = S),
                     dispersal = round(runif(S), 2),
                     row.names = sp_names)
coord <- cbind(x = seq_len(n), y = sin(seq_len(n)))
distmat <- as.matrix(dist(coord))
offset <- matrix(log(runif(n * S, 0.8, 1.2)), nrow = n, ncol = S, dimnames = list(site_names, sp_names))

cases <- list()

Y1 <- matrix(rpois(n * S, lambda = 2), nrow = n, dimnames = list(site_names, sp_names))
cfg1 <- base_config("poisson_lv_covariates")
cfg1$model$family <- "poisson"
cfg1$model$lv.control$num.lv <- 2
cases[[length(cases) + 1]] <- list(name = "01_poisson_lv_covariates", cfg = cfg1, Y = Y1, X = X, newdata = newdata)

Y2 <- matrix(rbinom(n * S, size = 5, prob = 0.45), nrow = n, dimnames = list(site_names, sp_names))
cfg2 <- base_config("binomial_trials_xind")
cfg2$model$family <- "binomial"
cfg2$model$lv.control$num.lv <- 0
cfg2$model$trial.size <- 5
cfg2$model$X.ind <- "1,1,1,1"
trial2 <- data.frame(t(rep(5, S)))
names(trial2) <- sp_names
cases[[length(cases) + 1]] <- list(name = "02_binomial_trials_xind", cfg = cfg2, Y = Y2, X = X, trials = trial2, newdata = newdata)

Y3 <- matrix(rnorm(n * S), nrow = n, dimnames = list(site_names, sp_names))
cfg3 <- base_config("normal_pure_ordination")
cfg3$model$family <- "normal"
cfg3$model$formula.X <- ""
cfg3$model$lv.control$num.lv <- 2
cfg3$diagnostics_outputs$predict <- FALSE
cases[[length(cases) + 1]] <- list(name = "03_normal_pure_ordination", cfg = cfg3, Y = Y3)

Y4 <- matrix(rnbinom(n * S, mu = 3, size = 2), nrow = n, dimnames = list(site_names, sp_names))
cfg4 <- base_config("negative_binomial_row_offset")
cfg4$model$family <- "negative.binomial"
cfg4$model$row.eff <- "random"
cfg4$model$use_offset <- TRUE
cfg4$diagnostics_outputs$ranefsplot <- FALSE
cases[[length(cases) + 1]] <- list(name = "04_negative_binomial_row_offset", cfg = cfg4, Y = Y4, X = X, rowids = rowids, offset = offset, newdata = newdata)

Y5 <- matrix(rpois(n * S, lambda = 2), nrow = n, dimnames = list(site_names, sp_names))
cfg5 <- base_config("spatial_exponential_distmat")
cfg5$model$family <- "poisson"
cfg5$model$lv.control$num.lv <- 1
cfg5$model$lv.control$type <- "exponential"
cases[[length(cases) + 1]] <- list(name = "05_spatial_exponential_distmat", cfg = cfg5, Y = Y5, X = X, distmat = distmat, newdata = newdata)

Y6 <- matrix(rpois(n * S, lambda = 2), nrow = n, dimnames = list(site_names, sp_names))
cfg6 <- base_config("traits_fourth_corner_ssvs")
cfg6$model$family <- "poisson"
cfg6$model$model_mode <- "trait/fourth-corner model"
cfg6$traits_random_ssvs$use_traits <- TRUE
cfg6$traits_random_ssvs$use_ssvs <- TRUE
cfg6$diagnostics_outputs$varpart <- TRUE
cases[[length(cases) + 1]] <- list(name = "06_traits_fourth_corner_ssvs", cfg = cfg6, Y = Y6, X = X, traits = traits, newdata = newdata)

Y7 <- cbind(
  sp_poisson = rpois(n, 2),
  sp_binom = rbinom(n, 5, 0.5),
  sp_normal = rnorm(n),
  sp_beta = pmin(pmax(rbeta(n, 2, 5), 0.01), 0.99),
  sp_ordinal = sample(1:3, n, replace = TRUE),
  sp_lnormal = rlnorm(n, 0, 0.3)
)
rownames(Y7) <- site_names
cfg7 <- base_config("mixed_family_full_upload_define")
cfg7$model$family_vector_override <- "poisson,binomial,normal,beta,ordinal,lnormal"
cfg7$model$trial.size <- 5
cfg7$model$lv.control$num.lv <- 1
cfg7$model$lv.control$type <- "squared.exponential"
cfg7$model$use_offset <- TRUE
cfg7$model$do.fit <- FALSE
cfg7$traits_random_ssvs$use_ranef <- TRUE
cfg7$diagnostics_outputs$predict <- FALSE
offset7 <- matrix(0, nrow = n, ncol = ncol(Y7), dimnames = dimnames(Y7))
cases[[length(cases) + 1]] <- list(name = "07_mixed_family_full_upload_define", cfg = cfg7, Y = Y7, X = X, rowids = rowids, ranefids = ranefids, distmat = distmat, offset = offset7, traits = traits[seq_len(ncol(Y7) %% S + 1), , drop = FALSE], trials = data.frame(t(c(0, 5, 0, 0, 0, 0))), newdata = newdata)

suite_root <- file.path(app_dir, "output", paste0("boral_example_suite_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(suite_root, recursive = TRUE, showWarnings = FALSE)

run_one <- function(case) {
  outdir <- file.path(suite_root, case$name)
  make_dirs(outdir)
  data_dir <- file.path(outdir, "data")
  input_dir <- file.path(outdir, "inputs")
  write_csv(case$Y, file.path(data_dir, "Y.csv"))
  write_csv(case$Y, file.path(input_dir, "Y.csv"))
  if (!is.null(case$X)) { write_csv(case$X, file.path(data_dir, "XData.csv")); write_csv(case$X, file.path(input_dir, "XData.csv")) }
  if (!is.null(case$traits)) { write_csv(case$traits, file.path(data_dir, "traits.csv")); write_csv(case$traits, file.path(input_dir, "traits.csv")) }
  if (!is.null(case$rowids)) { write_csv(case$rowids, file.path(data_dir, "row.ids.csv")); write_csv(case$rowids, file.path(input_dir, "row.ids.csv")) }
  if (!is.null(case$ranefids)) { write_csv(case$ranefids, file.path(data_dir, "ranef.ids.csv")); write_csv(case$ranefids, file.path(input_dir, "ranef.ids.csv")) }
  if (!is.null(case$distmat)) { write_csv(case$distmat, file.path(data_dir, "distmat.csv")); write_csv(case$distmat, file.path(input_dir, "distmat.csv")) }
  if (!is.null(case$offset)) { write_csv(case$offset, file.path(data_dir, "offset.csv")); write_csv(case$offset, file.path(input_dir, "offset.csv")) }
  if (!is.null(case$newdata)) { write_csv(case$newdata, file.path(data_dir, "newdata.csv")); write_csv(case$newdata, file.path(input_dir, "newdata.csv")) }
  if (!is.null(case$trials)) { write_csv(case$trials, file.path(data_dir, "trial.size.csv"), row_names = FALSE); write_csv(case$trials, file.path(input_dir, "trial.size.csv"), row_names = FALSE) }
  yaml::write_yaml(case$cfg, file.path(outdir, "used_config.yml"))
  write_check_messages(outdir, c("Synthetic boral suite case generated.", paste0("Case: ", case$name)))
  script <- write_boral_reproducible_script(outdir)
  stdout <- file.path(outdir, "diagnostics", "suite_stdout.log")
  stderr <- file.path(outdir, "diagnostics", "suite_stderr.log")
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  exit <- system2(rscript, shQuote(script), stdout = stdout, stderr = stderr)
  status_file <- file.path(outdir, "diagnostics", "engine_status.json")
  st <- if (file.exists(status_file)) tryCatch(jsonlite::fromJSON(status_file), error = function(e) list(status = "status_unreadable")) else list(status = "missing_status")
  zip <- make_zip(outdir)
  files <- list.files(outdir, recursive = TRUE)
  required <- c("used_config.yml", "diagnostics/engine_status.json", "diagnostics/data_check_messages.csv",
                "diagnostics/session_info.txt", "standard/run_summary.csv", "standard/output_manifest.csv",
                "tables/boral_result_workflow_map.csv", "reproducible_script/run_this_boral_analysis.R",
                "report/boral_report.html")
  missing <- required[!file.exists(file.path(outdir, required))]
  empty_dirs <- names(which(vapply(c("models", "tables", "plots", "diagnostics", "report", "predictions", "standard", "jags", "mcmc", "ordination", "residuals", "random_effects", "variable_selection"),
                                   function(d) length(list.files(file.path(outdir, d), all.files = FALSE, no.. = TRUE)) == 0, logical(1))))
  data.frame(case = case$name, status = st$status %||% "unknown", exit_code = exit,
             outdir = outdir, zip = zip, file_count = length(files),
             missing_required = paste(missing, collapse = ";"),
             empty_dirs = paste(empty_dirs, collapse = ";"),
             has_dependency_error = file.exists(file.path(outdir, "diagnostics", "boral_dependency_error.txt")),
             stringsAsFactors = FALSE)
}

results <- do.call(rbind, lapply(cases, run_one))
write.csv(results, file.path(example_dir, "last_real_example_suite_summary.csv"), row.names = FALSE)
write.csv(results, file.path(suite_root, "boral_suite_summary.csv"), row.names = FALSE)
print(results)
cat("boral example suite root:", suite_root, "\n")
