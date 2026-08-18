# Universal Benchmark runner for JSDM Studio.
# This script is intentionally command-line runnable so the GUI benchmark can be
# reproduced without Shiny.

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (length(x) == 1 && is.atomic(x) && is.na(x)) return(y)
  x
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return(default)
  if (is.logical(x)) return(isTRUE(x))
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

timestamp_id <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(action = "all", n_sites = 36L, n_species = 6L, n_visits = 3L,
              seed = 20260601L, engines = "Hmsc,Hmsc-HPC,jSDM,GJAM,spOccupancy,sjSDM,boral",
              master_dir = "", input_dir = "", target_dir = "")
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "true"
    out[[key]] <- val
  }
  out$n_sites <- as.integer(out$n_sites)
  out$n_species <- as.integer(out$n_species)
  out$n_visits <- as.integer(out$n_visits)
  out$seed <- as.integer(out$seed)
  out$engines <- trimws(unlist(strsplit(as.character(out$engines), ",")))
  out$engines <- out$engines[nzchar(out$engines)]
  out
}

script_path <- {
  args0 <- commandArgs(trailingOnly = FALSE)
  hit <- args0[startsWith(args0, "--file=")]
  if (length(hit)) normalizePath(sub("^--file=", "", hit[[1]]), winslash = "/", mustWork = FALSE)
  else normalizePath("workflow_scripts/universal_benchmark_runner.R", winslash = "/", mustWork = FALSE)
}
script_dir <- dirname(script_path)
app_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) {
  app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
if (!file.exists(file.path(app_dir, "app.R"))) {
  stop("Cannot locate JSDM Studio app.R from universal benchmark runner.", call. = FALSE)
}
setwd(app_dir)

safe_require <- function(pkg) requireNamespace(pkg, quietly = TRUE)

ensure_dirs <- function(outdir, engine = NULL) {
  dirs <- c("", "inputs", "data", "models", "tables", "results", "plots", "predictions",
            "diagnostics", "workflow_scripts", "reproducible_script", "standard", "report")
  key <- tolower(gsub("[^a-z0-9]+", "", engine %||% ""))
  if (grepl("hmschpc", key)) dirs <- c(dirs, "samples", "compile", "hdf5", "python_logs")
  if (grepl("jsdm", key) && !grepl("sjsdm", key)) dirs <- c(dirs, "mcmc")
  if (grepl("gjam", key)) dirs <- c(dirs, "chains", "sensitivity", "ordination", "missing_data")
  if (grepl("spoccupancy|spocc", key)) dirs <- c(dirs, "samples", "spatial", "model_assessment")
  if (grepl("sjsdm", key)) dirs <- c(dirs, "anova", "internal_structure", "importance", "weights", "residuals")
  if (grepl("boral", key)) dirs <- c(dirs, "jags", "mcmc", "ordination", "residuals", "random_effects", "variable_selection")
  if (grepl("hmsc", key) && !grepl("hmschpc", key)) dirs <- c(dirs, "samples")
  for (d in unique(dirs)) dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
  invisible(outdir)
}

write_csv_pair <- function(x, outdir, filename, row.names = TRUE) {
  if (is.null(x)) return(invisible(FALSE))
  write.csv(x, file.path(outdir, "data", filename), row.names = row.names)
  write.csv(x, file.path(outdir, "inputs", filename), row.names = row.names)
  invisible(TRUE)
}

write_text_pair <- function(x, outdir, filename) {
  if (is.null(x)) return(invisible(FALSE))
  writeLines(as.character(x), file.path(outdir, "data", filename), useBytes = TRUE)
  writeLines(as.character(x), file.path(outdir, "inputs", filename), useBytes = TRUE)
  invisible(TRUE)
}

make_zip_file <- function(outdir, zipfile = paste0(outdir, ".zip")) {
  if (!dir.exists(outdir)) stop("Cannot zip missing output folder: ", outdir, call. = FALSE)
  if (file.exists(zipfile)) unlink(zipfile)
  dir.create(dirname(zipfile), recursive = TRUE, showWarnings = FALSE)
  manifest <- file.path(outdir, "standard", "output_manifest.csv")
  dir.create(dirname(manifest), recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(file = list.files(outdir, recursive = TRUE), stringsAsFactors = FALSE),
            manifest, row.names = FALSE)
  rel_files <- list.files(outdir, recursive = TRUE, all.files = FALSE, no.. = TRUE)
  if (!length(rel_files)) stop("Cannot create ZIP because output folder is empty: ", outdir, call. = FALSE)
  if (safe_require("zip")) {
    zip::zipr(zipfile = zipfile, files = rel_files, root = outdir, recurse = FALSE, mode = "mirror")
  } else {
    old <- getwd()
    on.exit(setwd(old), add = TRUE)
    setwd(outdir)
    utils::zip(zipfile, files = rel_files, flags = "-r9Xq")
  }
  if (!file.exists(zipfile) || file.info(zipfile)$size <= 0) {
    stop("ZIP creation failed or produced an empty file: ", zipfile, call. = FALSE)
  }
  normalizePath(zipfile, winslash = "/", mustWork = FALSE)
}

write_status <- function(outdir, engine, status = "fit_failed", warnings = character(), errors = character(),
                         runtime_seconds = 0) {
  dir.create(file.path(outdir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  payload <- list(engine = engine, status = status, warnings = as.character(warnings),
                  errors = as.character(errors), runtime_seconds = runtime_seconds,
                  timestamp = as.character(Sys.time()))
  if (safe_require("jsonlite")) {
    writeLines(jsonlite::toJSON(payload, pretty = TRUE, auto_unbox = TRUE),
               file.path(outdir, "diagnostics", "engine_status.json"))
  } else {
    capture.output(str(payload), file = file.path(outdir, "diagnostics", "engine_status.json"))
  }
  write.csv(data.frame(engine = engine, status = status,
                       message = paste(c(warnings, errors), collapse = "; "),
                       stringsAsFactors = FALSE),
            file.path(outdir, "tables", "engine_status.csv"), row.names = FALSE)
  if (length(errors)) {
    writeLines(errors, file.path(outdir, "diagnostics", paste0(gsub("[^A-Za-z0-9]+", "_", engine), "_error.txt")))
  }
  try(writeLines(capture.output(utils::sessionInfo()), file.path(outdir, "diagnostics", "session_info.txt")), silent = TRUE)
  invisible(payload)
}

read_status <- function(outdir, engine = "unknown") {
  path <- file.path(outdir, "diagnostics", "engine_status.json")
  if (file.exists(path) && safe_require("jsonlite")) {
    x <- tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(x$status)) return(x)
  }
  run <- file.path(outdir, "standard", "run_summary.csv")
  if (file.exists(run)) {
    r <- tryCatch(read.csv(run, check.names = FALSE), error = function(e) NULL)
    if (!is.null(r) && nrow(r) && "status" %in% names(r)) return(list(engine = engine, status = as.character(r$status[1])))
  }
  list(engine = engine, status = "missing_status", warnings = character(), errors = character())
}

fill_empty_dirs <- function(outdir, engine) {
  dirs <- list.dirs(outdir, recursive = TRUE, full.names = TRUE)
  for (d in dirs) {
    if (normalizePath(d, winslash = "/", mustWork = FALSE) == normalizePath(outdir, winslash = "/", mustWork = FALSE)) next
    if (length(list.files(d, all.files = FALSE, no.. = TRUE)) == 0) {
      writeLines(c(paste0(engine, " benchmark output folder"),
                   "No engine-specific file was produced in this folder for the selected settings.",
                   "See diagnostics/engine_status.json and standard/*.csv."),
                 file.path(d, "README.txt"))
    }
  }
}

copy_dir_contents <- function(src, dest) {
  if (!dir.exists(src)) stop("Cannot copy missing directory: ", src, call. = FALSE)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(src, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  for (f in files) {
    rel <- substring(normalizePath(f, winslash = "/", mustWork = FALSE),
                     nchar(normalizePath(src, winslash = "/", mustWork = FALSE)) + 2L)
    target <- file.path(dest, rel)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    file.copy(f, target, overwrite = TRUE, copy.date = TRUE)
  }
  invisible(length(files))
}

make_tree_newick <- function(species) {
  parts <- species
  while (length(parts) > 1) {
    paired <- character()
    i <- 1L
    while (i <= length(parts)) {
      if (i == length(parts)) paired <- c(paired, parts[i])
      else paired <- c(paired, paste0("(", parts[i], ":0.2,", parts[i + 1L], ":0.2):0.2"))
      i <- i + 2L
    }
    parts <- paired
  }
  paste0(parts, ";")
}

generate_benchmark_data <- function(target_dir = file.path(app_dir, "examples", "universal_benchmark"),
                                    n_sites = 36L, n_species = 6L, n_visits = 3L, seed = 20260601L) {
  n_sites <- max(30L, min(40L, as.integer(n_sites)))
  n_species <- max(5L, min(8L, as.integer(n_species)))
  n_visits <- max(2L, min(5L, as.integer(n_visits)))
  set.seed(seed)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(target_dir, "truth"), recursive = TRUE, showWarnings = FALSE)

  site_id <- sprintf("site_%02d", seq_len(n_sites))
  species <- sprintf("sp_%02d", seq_len(n_species))
  substrate <- rep(c("sand", "loam", "clay", "rock"), length.out = n_sites)
  X <- data.frame(
    pH = round(seq(5.4, 7.5, length.out = n_sites) + rnorm(n_sites, 0, 0.05), 3),
    moisture = round(seq(0.15, 0.9, length.out = n_sites) + 0.05 * sin(seq_len(n_sites)), 3),
    canopy = round(seq(10, 85, length.out = n_sites) + rnorm(n_sites, 0, 2), 3),
    elevation = round(seq(80, 460, length.out = n_sites) + rnorm(n_sites, 0, 5), 3),
    substrate = substrate,
    stringsAsFactors = FALSE
  )
  rownames(X) <- site_id
  X_scaled <- data.frame(
    pH = as.numeric(scale(X$pH)),
    moisture = as.numeric(scale(X$moisture)),
    canopy = as.numeric(scale(X$canopy)),
    elevation = as.numeric(scale(X$elevation)),
    substrate_clay = as.numeric(X$substrate == "clay"),
    substrate_rock = as.numeric(X$substrate == "rock")
  )
  predictors <- names(X_scaled)

  beta <- matrix(0, nrow = n_species, ncol = length(predictors), dimnames = list(species, predictors))
  for (j in seq_len(n_species)) {
    beta[j, ] <- c(
      0.65 - 0.08 * j,
      -0.35 + 0.10 * j,
      ifelse(j %% 2 == 0, 0.38, -0.25),
      0.22 * sin(j),
      ifelse(j %% 3 == 0, 0.45, -0.10),
      ifelse(j %% 2 == 1, -0.32, 0.20)
    )
  }
  intercept <- seq(-0.85, 0.65, length.out = n_species)
  assoc <- matrix(0, n_species, n_species, dimnames = list(species, species))
  for (i in seq_len(n_species)) for (j in seq_len(n_species)) assoc[i, j] <- 0.45 ^ abs(i - j)
  diag(assoc) <- 1
  latent <- matrix(rnorm(n_sites * n_species, 0, 0.20), n_sites, n_species)
  eta <- as.matrix(X_scaled) %*% t(beta) + matrix(rep(intercept, each = n_sites), n_sites, n_species) + latent
  colnames(eta) <- species
  rownames(eta) <- site_id
  prob <- pmin(pmax(plogis(eta), 0.05), 0.95)
  Y_occ <- matrix(rbinom(n_sites * n_species, 1, as.vector(prob)), n_sites, n_species,
                  dimnames = list(site_id, species))
  for (j in seq_len(n_species)) {
    m <- mean(Y_occ[, j])
    if (m < 0.12 || m > 0.88) Y_occ[, j] <- as.integer(prob[, j] > stats::quantile(prob[, j], 0.5))
  }
  lambda <- pmax(exp(0.4 + eta / 2), 0.05)
  Y_count <- matrix(rpois(n_sites * n_species, lambda = as.vector(lambda)), n_sites, n_species,
                    dimnames = list(site_id, species))
  Y_normal <- eta + matrix(rnorm(n_sites * n_species, 0, 0.25), n_sites, n_species)
  dimnames(Y_normal) <- list(site_id, species)

  traits <- data.frame(
    life_form = rep(c("grass", "forb", "shrub", "tree"), length.out = n_species),
    height_mm = round(seq(40, 420, length.out = n_species), 1),
    dispersal = round(seq(0.15, 0.95, length.out = n_species), 3),
    stringsAsFactors = FALSE
  )
  rownames(traits) <- species
  study <- data.frame(
    sample = site_id,
    plot = rep(sprintf("plot_%02d", seq_len(ceiling(n_sites / 3))), each = 3, length.out = n_sites),
    year = rep(c("year_1", "year_2", "year_3"), length.out = n_sites),
    stringsAsFactors = FALSE
  )
  rownames(study) <- site_id
  coords <- data.frame(
    x = round(rep(seq(0, 5, length.out = 6), length.out = n_sites) + rnorm(n_sites, 0, 0.02), 4),
    y = round(rep(seq(0, 5, length.out = ceiling(n_sites / 6)), each = 6, length.out = n_sites) + rnorm(n_sites, 0, 0.02), 4)
  )
  rownames(coords) <- site_id
  newdata <- X[seq_len(min(8L, n_sites)), , drop = FALSE]
  rownames(newdata) <- sprintf("new_site_%02d", seq_len(nrow(newdata)))
  newcoords <- coords[seq_len(nrow(newdata)), , drop = FALSE]
  rownames(newcoords) <- rownames(newdata)
  folds <- data.frame(site_id = site_id, fold = rep(seq_len(3), length.out = n_sites), stringsAsFactors = FALSE)
  trial <- as.data.frame(matrix(1, nrow = n_sites, ncol = n_species, dimnames = list(site_id, species)))
  offset <- as.data.frame(matrix(0, nrow = n_sites, ncol = n_species, dimnames = list(site_id, species)))
  rowids <- data.frame(row_id = site_id, block = study$plot, stringsAsFactors = FALSE)
  rownames(rowids) <- site_id
  ranefids <- data.frame(observer = rep(c("obs_a", "obs_b", "obs_c"), length.out = n_sites), stringsAsFactors = FALSE)
  rownames(ranefids) <- site_id
  distmat <- as.matrix(stats::dist(coords))
  rownames(distmat) <- colnames(distmat) <- site_id
  phylo_cov <- assoc
  newick <- make_tree_newick(species)

  det_cov <- data.frame(row.names = site_id)
  for (r in seq_len(n_visits)) det_cov[[paste0("obs_rep", r)]] <- round(runif(n_sites, 0.2, 1.0), 3)
  det <- data.frame(row.names = site_id)
  for (sp in species) {
    for (r in seq_len(n_visits)) {
      pdet <- plogis(-0.2 + 0.7 * det_cov[[paste0("obs_rep", r)]] + 0.25 * X_scaled$moisture)
      det[[paste0(sp, "_rep", r)]] <- as.integer(Y_occ[, sp] * rbinom(n_sites, 1, pdet))
    }
  }
  occ_covs <- X

  write.csv(Y_occ, file.path(target_dir, "Y_occurrence.csv"))
  write.csv(Y_count, file.path(target_dir, "Y_count.csv"))
  write.csv(round(Y_normal, 4), file.path(target_dir, "Y_normal.csv"))
  write.csv(X, file.path(target_dir, "XData.csv"))
  write.csv(traits, file.path(target_dir, "traits.csv"))
  write.csv(study, file.path(target_dir, "studyDesign.csv"))
  write.csv(coords, file.path(target_dir, "coordinates.csv"))
  writeLines(newick, file.path(target_dir, "phylogeny.nwk"), useBytes = TRUE)
  write.csv(phylo_cov, file.path(target_dir, "phylo_cov.csv"))
  write.csv(newdata, file.path(target_dir, "newdata.csv"))
  write.csv(newcoords, file.path(target_dir, "newcoords.csv"))
  write.csv(folds, file.path(target_dir, "folds.csv"), row.names = FALSE)
  write.csv(trial, file.path(target_dir, "trial.size.csv"))
  write.csv(offset, file.path(target_dir, "offset.csv"))
  write.csv(rowids, file.path(target_dir, "row.ids.csv"))
  write.csv(ranefids, file.path(target_dir, "ranef.ids.csv"))
  write.csv(distmat, file.path(target_dir, "distmat.csv"))
  write.csv(det, file.path(target_dir, "spOccupancy_y_detection.csv"))
  write.csv(occ_covs, file.path(target_dir, "occ.covs.csv"))
  write.csv(det_cov, file.path(target_dir, "det.covs.csv"))

  true_eff <- as.data.frame(as.table(beta), stringsAsFactors = FALSE)
  names(true_eff) <- c("species", "predictor", "true_effect")
  true_assoc <- as.data.frame(as.table(assoc), stringsAsFactors = FALSE)
  names(true_assoc) <- c("species_i", "species_j", "true_association")
  true_assoc <- true_assoc[as.integer(factor(true_assoc$species_i, levels = species)) <
                             as.integer(factor(true_assoc$species_j, levels = species)), , drop = FALSE]
  pred_truth <- as.data.frame(as.table(prob), stringsAsFactors = FALSE)
  names(pred_truth) <- c("site_id", "species", "truth_probability")
  write.csv(round(eta, 4), file.path(target_dir, "truth", "latent_occurrence.csv"))
  write.csv(round(prob, 6), file.path(target_dir, "truth", "occurrence_probability.csv"))
  write.csv(true_eff, file.path(target_dir, "truth", "true_environment_effects.csv"), row.names = FALSE)
  write.csv(true_assoc, file.path(target_dir, "truth", "true_species_associations.csv"), row.names = FALSE)
  write.csv(pred_truth, file.path(target_dir, "truth", "true_predictions.csv"), row.names = FALSE)
  cfg <- list(seed = seed, n_sites = n_sites, n_species = n_species, n_visits = n_visits,
              site_id = site_id, species = species, predictors = predictors,
              note = "All engine-specific inputs are derived from the same latent ecological truth.")
  if (safe_require("yaml")) yaml::write_yaml(cfg, file.path(target_dir, "truth", "data_generation_config.yml"))
  writeLines(c(
    "# Rebuild the JSDM Studio universal benchmark data.",
    "runner <- file.path('..', '..', 'workflow_scripts', 'universal_benchmark_runner.R')",
    "system2(file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript'),",
    "        c(runner, '--action=generate'))"
  ), file.path(target_dir, "truth", "data_generation_script.R"))
  write.csv(data.frame(file = list.files(target_dir, recursive = TRUE), stringsAsFactors = FALSE),
            file.path(target_dir, "benchmark_input_manifest.csv"), row.names = FALSE)
  list(dir = normalizePath(target_dir, winslash = "/", mustWork = FALSE), Y = Y_occ, X = X,
       traits = traits, study = study, coords = coords, phylo_cov = phylo_cov, newick = newick,
       newdata = newdata, newcoords = newcoords, det = det, occ_covs = occ_covs, det_cov = det_cov,
       rowids = rowids, ranefids = ranefids, distmat = distmat, offset = offset, trial = trial,
       truth_effects = true_eff, truth_predictions = pred_truth, truth_associations = true_assoc)
}

benchmark_required_files <- function() {
  c("Y_occurrence.csv", "Y_count.csv", "Y_normal.csv", "XData.csv", "traits.csv",
    "studyDesign.csv", "coordinates.csv", "phylogeny.nwk", "phylo_cov.csv",
    "newdata.csv", "newcoords.csv", "folds.csv", "trial.size.csv", "offset.csv",
    "row.ids.csv", "ranef.ids.csv", "distmat.csv", "spOccupancy_y_detection.csv",
    "occ.covs.csv", "det.covs.csv", file.path("truth", "latent_occurrence.csv"),
    file.path("truth", "occurrence_probability.csv"),
    file.path("truth", "true_environment_effects.csv"),
    file.path("truth", "true_species_associations.csv"),
    file.path("truth", "true_predictions.csv"),
    file.path("truth", "data_generation_config.yml"),
    file.path("truth", "data_generation_script.R"))
}

validate_benchmark_data_dir <- function(dir, allow_generate = FALSE) {
  if (!dir.exists(dir)) {
    if (allow_generate) return(invisible(FALSE))
    stop("Benchmark input directory does not exist: ", dir, call. = FALSE)
  }
  missing <- benchmark_required_files()[!file.exists(file.path(dir, benchmark_required_files()))]
  if (length(missing)) {
    if (allow_generate) return(invisible(FALSE))
    stop("Benchmark input directory is missing required files: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  TRUE
}

load_benchmark_data <- function(dir = file.path(app_dir, "examples", "universal_benchmark"),
                                allow_generate = TRUE) {
  if (!file.exists(file.path(dir, "Y_occurrence.csv"))) {
    if (allow_generate) generate_benchmark_data(dir)
    else stop("Benchmark input directory is missing Y_occurrence.csv: ", dir, call. = FALSE)
  }
  validate_benchmark_data_dir(dir, allow_generate = FALSE)
  read_csv <- function(path, row.names = 1) read.csv(file.path(dir, path), row.names = row.names, check.names = FALSE, stringsAsFactors = FALSE)
  list(
    dir = normalizePath(dir, winslash = "/", mustWork = FALSE),
    Y = read_csv("Y_occurrence.csv"),
    Y_count = read_csv("Y_count.csv"),
    Y_normal = read_csv("Y_normal.csv"),
    X = read_csv("XData.csv"),
    traits = read_csv("traits.csv"),
    study = read_csv("studyDesign.csv"),
    coords = read_csv("coordinates.csv"),
    phylo_cov = read_csv("phylo_cov.csv"),
    newick = paste(readLines(file.path(dir, "phylogeny.nwk"), warn = FALSE), collapse = ""),
    newdata = read_csv("newdata.csv"),
    newcoords = read_csv("newcoords.csv"),
    folds = read.csv(file.path(dir, "folds.csv"), check.names = FALSE, stringsAsFactors = FALSE),
    trial = read_csv("trial.size.csv"),
    offset = read_csv("offset.csv"),
    rowids = read_csv("row.ids.csv"),
    ranefids = read_csv("ranef.ids.csv"),
    distmat = read_csv("distmat.csv"),
    det = read_csv("spOccupancy_y_detection.csv"),
    occ_covs = read_csv("occ.covs.csv"),
    det_cov = read_csv("det.covs.csv"),
    truth_effects = read.csv(file.path(dir, "truth", "true_environment_effects.csv"), check.names = FALSE, stringsAsFactors = FALSE),
    truth_predictions = read.csv(file.path(dir, "truth", "true_predictions.csv"), check.names = FALSE, stringsAsFactors = FALSE),
    truth_associations = read.csv(file.path(dir, "truth", "true_species_associations.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  )
}

source_app_env <- function() {
  e <- new.env(parent = globalenv())
  sys.source(file.path(app_dir, "app.R"), envir = e)
  for (f in c("jsdm_adapter.R", "gjam_adapter.R", "spoccupancy_adapter.R", "boral_adapter.R", "hmschpc_adapter.R")) {
    path <- file.path(app_dir, "R", f)
    if (file.exists(path)) sys.source(path, envir = e)
  }
  e
}

dependency_preflight <- function(master_dir = NULL) {
  pkgs <- c("shiny", "DT", "yaml", "jsonlite", "zip", "Hmsc", "coda", "ape", "jSDM",
            "gjam", "spOccupancy", "sjSDM", "reticulate", "torch", "boral", "rjags", "R2jags")
  rows <- lapply(pkgs, function(pkg) {
    ok <- safe_require(pkg)
    data.frame(component = pkg, kind = "R package", available = ok,
               version = if (ok) tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "") else "",
               note = if (ok) "available" else "not installed or not loadable",
               stringsAsFactors = FALSE)
  })
  py <- Sys.which("python")
  rows[[length(rows) + 1L]] <- data.frame(component = "python", kind = "runtime", available = nzchar(py),
                                          version = py, note = if (nzchar(py)) "found on PATH" else "not found on PATH",
                                          stringsAsFactors = FALSE)
  python_module_row <- function(module, import_name = module) {
    if (!nzchar(py)) {
      return(data.frame(component = paste0("Hmsc-HPC Python:", module),
                        kind = "Python package", available = FALSE, version = "",
                        note = "Python runtime not found", stringsAsFactors = FALSE))
    }
    source_dir <- normalizePath(file.path(app_dir, "external_packages", "hmsc-hpc-main"),
                                winslash = "/", mustWork = FALSE)
    code <- sprintf(
      "import os, sys, importlib; os.environ['TF_CPP_MIN_LOG_LEVEL']='2'; sys.path.insert(0, %s); m=importlib.import_module(%s); print(getattr(m, '__version__', 'available'))",
      shQuote(source_dir), shQuote(import_name)
    )
    out <- tryCatch(system2(py, c("-c", shQuote(code)), stdout = TRUE, stderr = FALSE),
                    error = function(e) structure(conditionMessage(e), status = 1L))
    ok <- is.null(attr(out, "status")) || identical(attr(out, "status"), 0L)
    data.frame(component = paste0("Hmsc-HPC Python:", module),
               kind = "Python package", available = ok,
               version = if (ok && length(out)) trimws(out[[1]]) else "",
               note = if (ok) "available" else "install with install_hmschpc_python_packages.bat or selected Python",
               stringsAsFactors = FALSE)
  }
  for (item in list(c("numpy", "numpy"), c("pandas", "pandas"), c("patsy", "patsy"),
                    c("PyYAML", "yaml"), c("h5py", "h5py"), c("scipy", "scipy"),
                    c("tensorflow", "tensorflow"), c("tensorflow_probability", "tensorflow_probability"),
                    c("ujson", "ujson"), c("pyhmsc", "pyhmsc"))) {
    rows[[length(rows) + 1L]] <- python_module_row(item[[1]], item[[2]])
  }
  rows[[length(rows) + 1L]] <- data.frame(component = "JAGS", kind = "system dependency",
                                          available = safe_require("rjags"),
                                          version = if (safe_require("rjags")) "rjags loadable" else "",
                                          note = "boral requires system JAGS plus rjags/R2jags",
                                          stringsAsFactors = FALSE)
  out <- do.call(rbind, rows)
  if (!is.null(master_dir)) {
    dir.create(file.path(master_dir, "07_reports"), recursive = TRUE, showWarnings = FALSE)
    write.csv(out, file.path(master_dir, "07_reports", "dependency_preflight_report.csv"), row.names = FALSE)
  }
  out
}

hmsc_cfg <- function() {
  list(
    project_name = "Universal_Benchmark_Hmsc",
    engine = "Hmsc",
    data = list(phylogeny = "phylogeny.nwk"),
    model = list(
      distr = "probit", XFormula = "~ pH + moisture + canopy + elevation + substrate",
      TrFormula = "~ life_form + height_mm + dispersal", use_traits = TRUE,
      use_phylogeny = TRUE, random_mode = "sample", random_effect_column = "plot",
      spatial_method = "Full", nNeighbours = 10L, lon_col = "x", lat_col = "y",
      seed = 20260601L, XScale = TRUE, TrScale = TRUE, YScale = FALSE,
      truncateNumberOfFactors = TRUE, Loff_file = "", ranLevelsUsed = "",
      C_file = "phylo_cov.csv", use_XRRR = FALSE, ncRRR = 2L,
      XRRRFormula = "~ .", XRRRScale = TRUE, XRRR_file = "",
      random_level_type = "none", sMethod = "Full", random_N = 1L,
      longlat = FALSE, units_column = "sample", distMat_file = "",
      xData_file = "", sKnot_file = "", nfMin = 1L, nfMax = 4L,
      priors = list(setDefault = TRUE, a1 = NA, b1 = NA, a2 = NA, b2 = NA, alphapw = "")
    ),
    mcmc = list(samples = 5L, transient = 5L, thin = 1L, nChains = 2L, nParallel = 1L,
                verbose = 0L, preset = "Universal quick benchmark", initPar = "fixed effects",
                alignPost = TRUE, updater = list(GammaEta = TRUE, Beta = TRUE, Gamma = TRUE, Omega = TRUE),
                sample_prior = FALSE, pool_chains = FALSE),
    outputs = list(save_model = TRUE, predicted = TRUE, fit = TRUE, cv = FALSE, waic = FALSE,
                   diagnostics = TRUE, parameters = TRUE, variance_partitioning = TRUE,
                   omega = TRUE, gradients = FALSE, beta_support = TRUE, gamma_support = TRUE,
                   omega_support = TRUE,
                   convergence = list(showBeta = TRUE, showGamma = TRUE, showOmega = TRUE,
                                      maxOmega = 10L, showRho = TRUE, showAlpha = TRUE,
                                      effectiveSize = TRUE, gelmanPSRF = TRUE),
                   plotting = list(var.part.order.explained = TRUE, var.part.order.raw = FALSE,
                                   show.sp.names.beta = FALSE, plotTree = FALSE,
                                   omega.order = "original", show.sp.names.omega = TRUE,
                                   plotBeta = TRUE, plotGamma = TRUE),
                   predictions = list(species.list = "", trait.list = "", env.list = "",
                                      nfolds = 2L, partition_column = "plot", computeSAIR = FALSE))
  )
}

hmschpc_cfg <- function(app_env) {
  list(
    engine = "Hmsc-HPC",
    inputs = list(Y = "Y.csv", XData = "XData.csv", traits = "traits.csv",
                  studyDesign = "studyDesign.csv", coordinates = NULL,
                  phylo_cov = "phylo_cov.csv", phylo_tree = NULL, newdata = "newdata.csv"),
    model = list(distribution = "probit", XFormula = "~ pH + moisture + canopy + elevation + C(substrate)",
                 use_traits = TRUE, trait_formula = "~ height_mm + dispersal + C(life_form)",
                 phylogeny_mode = "covariance"),
    random_effects = list(mode = "iid", name = "plot_iid", column = "plot",
                          coord_x = "x", coord_y = "y", x_formula = "",
                          nf = 1L, nfMin = 1L, nfMax = 4L, alpha = 1),
    sampler = list(run_sampler = TRUE, samples = 3L, transient = 3L, thin = 1L, chains = 1L,
                   chains_to_run = NULL, verbose = 0L, seed = 20260601L, precision = 32L,
                   truncated_normal_library = "tf", hmcleapfrog = 5L, hmcthin = 0L,
                   update_beta_eta = FALSE, save_eta = TRUE, eager = FALSE, profile = FALSE),
    outputs = list(predictions = TRUE, diagnostics = TRUE, plots = TRUE, zip = TRUE),
    runtime = list(python = app_env$hmschpc_python_bin(list()),
                   python_source = app_env$hmschpc_python_source_dir())
  )
}

jsdm_cfg <- function() {
  list(
    project_name = "Universal_Benchmark_jSDM", engine = "jSDM",
    data = list(Y = "Y.csv", XData = "XData.csv", trait_data = "trait_data.csv",
                long_format = "long_format.csv", trials = "trials.csv",
                newdata = "newdata.csv", prediction_ids = "prediction_ids.csv"),
    model = list(model_type = "binomial_probit", response_argument = "",
                 site_formula = "~ pH + moisture + canopy + elevation + substrate",
                 trait_formula = "~ life_form + height_mm + dispersal",
                 n_latent = 2L, site_effect = "random", trials = 1L,
                 constrained_latent = FALSE, constrained_nchains = 2L,
                 long_site_col = "site", long_species_col = "species",
                 long_response_col = "presence", scale_site_data = TRUE,
                 include_intercept = TRUE, allow_traits = TRUE),
    prediction = list(do_predict = TRUE, predict_type = "mean", predict_probs = "0.1,0.9",
                      max_prediction_sites = 8L, Id_sites = "all", Id_species = "all",
                      prediction_histograms = TRUE),
    diagnostics = list(cor_prob = 0.95, cor_type = "mean", plot_residual_cor = TRUE,
                       plot_associations = TRUE, diag_beta = TRUE, diag_lambda = TRUE,
                       diag_W = TRUE, diag_alpha = TRUE, diag_Valpha = TRUE,
                       diag_V = TRUE, diag_deviance = TRUE, coda_summary = TRUE),
    mcmc = list(burnin = 50L, mcmc = 50L, thin = 1L, seed = 20260601L,
                verbose = 0L, preset = "Universal quick benchmark", ropt = 0.44),
    starts = list(beta_start = 0, gamma_start = 0, lambda_start = 0,
                  W_start = 0, alpha_start = 0, V_alpha = 1, V_start = 1),
    priors = list(shape_Valpha = 0.5, rate_Valpha = 0.0005, shape_V = 0.5,
                  rate_V = 0.0005, mu_beta = 0, V_beta = 10, mu_gamma = 0,
                  V_gamma = 10, mu_lambda = 0, V_lambda = 10),
    outputs = list(real_fit = TRUE, residual_cor = TRUE, enviro_cor = TRUE,
                   traceplots = TRUE, predictions = TRUE, save_model = TRUE,
                   report = TRUE, save_mcmc_rds = TRUE, csv_tables = TRUE,
                   figures = TRUE, table_model_spec = TRUE, table_beta = TRUE,
                   table_lambda = TRUE, table_gamma = TRUE, table_alpha = TRUE,
                   table_predictions = TRUE, copy_inputs = TRUE, save_config = TRUE,
                   zip = TRUE)
  )
}

gjam_cfg <- function(species) {
  list(
    project_name = "Universal_Benchmark_GJAM", engine = "GJAM",
    data = list(Y = "Y.csv", XData = "XData.csv", typeNames = "typeNames.csv",
                censor = "censor.csv", effort = "effort.csv", newdata = "newdata.csv",
                specByTrait = "specByTrait.csv", traitTypes = "traitTypes.csv",
                holdoutIndex = "holdoutIndex.csv"),
    model = list(formula = "~ pH + moisture + canopy + elevation + substrate",
                 type_single = "PA", typeNames_text = paste(rep("PA", length(species)), collapse = ","),
                 notStandard = "", ng = 100L, burnin = 20L, holdoutN = 0L,
                 seed = 20260601L, random = "", FULL = FALSE, PREDICTX = TRUE,
                 REDUCT = FALSE, reductList = list(N = 20L, r = 3L), ematAlpha = 0.5,
                 use_censor = FALSE, use_effort = FALSE),
    response_types = list(FCgroups = "", CCgroups = "", composition_reference = "not used",
                          trimY = FALSE, trim_minObs = 2L),
    priors_censoring = list(use_prior_template = FALSE, prior_file = NULL,
                            prior_mode = "default", censor_columns = "",
                            censor_values = "", censor_intervals = ""),
    analysis = list(do_predict = TRUE, do_sensitivity = TRUE, do_ordination = FALSE,
                    do_conditional = FALSE, do_iie = FALSE, do_traits = TRUE,
                    missingX = FALSE, missingY = FALSE, inverse_prediction = FALSE),
    outputs = list(real_fit = TRUE, save_model = TRUE, save_chains = TRUE,
                   save_parameters = TRUE, save_fit = TRUE, save_prediction = TRUE,
                   save_missing = TRUE, save_plots = TRUE, report = TRUE, zip = TRUE,
                   copy_inputs = TRUE, save_config = TRUE, csv_tables = TRUE)
  )
}

spocc_cfg <- function() {
  list(
    project_name = "Universal_Benchmark_spOccupancy", engine = "spOccupancy",
    data = list(y = "y.csv", occ.covs = "occ.covs.csv", det.covs = "det.covs.csv",
                coords = "coords.csv", species = "species.csv",
                integrated_sources = "integrated_sources.csv",
                newdata = "newdata.csv", newcoords = "newcoords.csv", folds = "folds.csv"),
    model = list(model_type = "msPGOcc", occ.formula = "~ pH + moisture + substrate",
                 det.formula = "~ obs", formula = "~ pH + moisture",
                 data_structure = "multi-species replicated", range.ind = FALSE),
    spatial_latent_svc = list(cov.model = "exponential", NNGP = TRUE, n.neighbors = 3L,
                              search.type = "cb", n.factors = 2L, svc.cols = "",
                              ar1 = FALSE, x.positive = FALSE),
    mcmc = list(n.batch = 5L, batch.length = 5L, n.burn = 5L, n.thin = 1L,
                n.chains = 1L, accept.rate = 0.43, n.report = 5L,
                n.omp.threads = 1L, verbose = FALSE, seed = 20260601L,
                updateMCMC = FALSE),
    priors_inits_tuning = list(beta.normal = "mean=0,var=2.72",
                               alpha.normal = "mean=0,var=2.72",
                               community_priors = "default", sigma.sq.ig = "a=2,b=1",
                               phi.unif = "maxdist", nu.unif = "0.5,2.5",
                               inits = "auto", tuning = "phi=0.5", fix = TRUE),
    validation_prediction_outputs = list(real_fit = TRUE, ppcOcc = FALSE,
                                         waicOcc = TRUE, k.fold = 0L,
                                         k.fold.threads = 1L, k.fold.seed = 100L,
                                         k.fold.only = FALSE, predict = FALSE,
                                         fitted = TRUE, save_model = TRUE,
                                         save_samples = TRUE, save_plots = FALSE,
                                         zip = TRUE, copy_inputs = TRUE,
                                         save_config = TRUE, report = TRUE)
  )
}

sjsdm_cfg <- function() {
  list(
    project_name = "Universal_Benchmark_sjSDM", engine = "sjSDM",
    model = list(family = "binomial_probit", env_model = "linear",
                 env_formula = "~ pH + moisture + canopy + elevation + substrate",
                 spatial_model = "linear", spatial_formula = "~ 0 + .",
                 se = FALSE, iter = 5L, step_size = 8L, sampling = 100L,
                 parallel = 0L, dtype = "float32", verbose = FALSE, seed = 20260601L),
    regularization_biotic = list(env_lambda = 0, env_alpha = 0.5,
                                 spatial_lambda = 0, spatial_alpha = 0.5,
                                 biotic_lambda = 0.01, biotic_alpha = 0.5,
                                 biotic_df = NA, on_diag = FALSE,
                                 reg_on_Cov = TRUE, inverse = FALSE,
                                 tune_regularization = FALSE, tune_steps = 0L, cv_k = 0L),
    dnn_optimizer = list(hidden = "6,4", activation = "selu", dropout = 0,
                         bias = TRUE, optimizer = "Adamax", learning_rate = 0.003,
                         weight_decay = 0.001, device = "cpu"),
    control = list(scheduler = 0L, lr_reduce_factor = 0.99,
                   early_stopping_training = 0L, mixed = FALSE),
    spatial_anova_metacommunity = list(generateSpatialEV = FALSE, spatial_ev_k = 4L,
                                       spatial_ev_threshold = 0, include_space_in_anova = TRUE,
                                       do_anova = TRUE, anova_samples = 100L,
                                       do_internal = TRUE, internal_fractions = "proportional",
                                       do_assembly = FALSE, assembly_predictor = ""),
    outputs = list(real_fit = TRUE, predict = TRUE, Rsquared = TRUE, importance = TRUE,
                   weights = TRUE, coef = TRUE, covariance_correlation = TRUE,
                   residuals = TRUE, plots = TRUE, zip = TRUE, copy_inputs = TRUE,
                   save_config = TRUE, report = TRUE)
  )
}

boral_cfg <- function() {
  list(
    project_name = "Universal_Benchmark_boral", engine = "boral",
    data = list(Y = "Y.csv", XData = "XData.csv", traits = "traits.csv",
                row.ids = "row.ids.csv", ranef.ids = "ranef.ids.csv",
                distmat = "distmat.csv", offset = "offset.csv",
                newdata = "newdata.csv", trial.size.file = "trial.size.csv"),
    model = list(model_mode = "correlated response GLMs", family = "binomial",
                 family_vector_override = "", lv.control = list(num.lv = 1L, type = "independent"),
                 model.name = "jagsboralmodel.txt", formula.X = "~ pH + moisture + canopy + elevation + substrate",
                 X.ind = "", trial.size = 1L, row.eff = "none",
                 use_offset = FALSE, do.fit = TRUE),
    traits_random_ssvs = list(use_traits = FALSE, which.traits = "",
                              traits_no_intercept = TRUE, use_ranef = FALSE,
                              use_ssvs = FALSE, ssvs.index = "",
                              ssvs.traitsindex = "", ssvs.g = 1e-6,
                              save.model = FALSE),
    mcmc_prior = list(n.burnin = 50L, n.iteration = 120L, n.thin = 5L,
                      seed = 20260601L, prior.type = "normal,normal,normal,uniform",
                      hypparams = "10,10,10,30", calc.ics = FALSE,
                      save_mcmc_samples = TRUE, save_hpd = TRUE, save_dic = TRUE),
    diagnostics_outputs = list(summary = TRUE, residual_plot = TRUE, lvsplot = TRUE,
                               ind.spp = "default", ranefsplot = FALSE,
                               coefsplot = TRUE, enviro_cor = TRUE,
                               residual_cor = TRUE, varpart = FALSE,
                               predict = TRUE, fitted = TRUE, tidyboral = TRUE),
    outputs = list(save_model = TRUE, save_jags = TRUE, save_tables = TRUE,
                   save_plots = TRUE, report = TRUE, zip = TRUE,
                   copy_inputs = TRUE, save_config = TRUE)
  )
}

run_rscript <- function(script, outdir, engine) {
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  stdout <- file.path(outdir, "diagnostics", paste0(gsub("[^A-Za-z0-9]+", "_", engine), "_stdout.log"))
  stderr <- file.path(outdir, "diagnostics", paste0(gsub("[^A-Za-z0-9]+", "_", engine), "_stderr.log"))
  system2(rscript, shQuote(script), stdout = stdout, stderr = stderr)
}

write_engine_inputs <- function(outdir, dat, engine) {
  write_csv_pair(dat$Y, outdir, "Y.csv")
  write_csv_pair(dat$X, outdir, "XData.csv")
  write_csv_pair(dat$traits, outdir, if (engine == "jSDM") "trait_data.csv" else "traits.csv")
  write_csv_pair(dat$study, outdir, "studyDesign.csv")
  write_csv_pair(dat$coords, outdir, if (engine == "spOccupancy") "coords.csv" else "coordinates.csv")
  write_csv_pair(dat$phylo_cov, outdir, "phylo_cov.csv")
  write_text_pair(dat$newick, outdir, "phylogeny.nwk")
  write_csv_pair(dat$newdata, outdir, "newdata.csv")
  write_csv_pair(dat$newcoords, outdir, "newcoords.csv")
  write_csv_pair(dat$folds, outdir, "folds.csv", row.names = FALSE)
  write_csv_pair(dat$trial, outdir, "trial.size.csv")
  write_csv_pair(dat$offset, outdir, "offset.csv")
  write_csv_pair(dat$rowids, outdir, "row.ids.csv")
  write_csv_pair(dat$ranefids, outdir, "ranef.ids.csv")
  write_csv_pair(dat$distmat, outdir, "distmat.csv")
  invisible(TRUE)
}

standardize_engine_synonyms <- function(outdir, engine, dat) {
  dir.create(file.path(outdir, "standard"), recursive = TRUE, showWarnings = FALSE)
  st <- read_status(outdir, engine)
  status <- st$status %||% "unknown"
  species <- colnames(dat$Y)
  predictors <- unique(dat$truth_effects$predictor)
  read_std <- function(name) {
    p <- file.path(outdir, "standard", name)
    if (file.exists(p)) tryCatch(read.csv(p, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  }
  effects_long <- read_std("effects_long.csv")
  effects <- NULL
  if (!is.null(effects_long) && nrow(effects_long) && "estimate" %in% names(effects_long)) {
    sp_col <- intersect(c("species", "response_id", "response", "column_id"), names(effects_long))[1]
    pr_col <- intersect(c("predictor", "parameter", "row_id"), names(effects_long))[1]
    if (!is.na(sp_col) && !is.na(pr_col)) {
      effects <- data.frame(
        engine = engine,
        species = as.character(effects_long[[sp_col]]),
        predictor = as.character(effects_long[[pr_col]]),
        estimate = suppressWarnings(as.numeric(effects_long$estimate)),
        lower = if ("lower" %in% names(effects_long)) suppressWarnings(as.numeric(effects_long$lower)) else NA_real_,
        upper = if ("upper" %in% names(effects_long)) suppressWarnings(as.numeric(effects_long$upper)) else NA_real_,
        statistic = NA_real_,
        p_or_support = NA_real_,
        effect_type = if (identical(engine, "spOccupancy")) "occurrence" else "environment",
        scale = switch(engine, GJAM = "observation_or_latent_GJAM", spOccupancy = "occupancy_link", "engine_link"),
        comparable = identical(status, "fitted"),
        note = paste0("Mapped from standard/effects_long.csv; status=", status),
        stringsAsFactors = FALSE
      )
      effects <- effects[nzchar(effects$species) & nzchar(effects$predictor), , drop = FALSE]
    }
  }
  if (is.null(effects) || !nrow(effects)) {
    base <- expand.grid(species = species, predictor = predictors, stringsAsFactors = FALSE)
    effects <- data.frame(engine = engine, base, estimate = NA_real_, lower = NA_real_,
                          upper = NA_real_, statistic = NA_real_, p_or_support = NA_real_,
                          effect_type = if (identical(engine, "spOccupancy")) "occurrence" else "environment",
                          scale = "not_fitted_or_not_exported", comparable = FALSE,
                          note = paste0("No fitted effect table exported; status=", status),
                          stringsAsFactors = FALSE)
  }
  write.csv(effects, file.path(outdir, "standard", "effects_species_environment.csv"), row.names = FALSE)

  pred_long <- read_std("predictions_long.csv")
  preds <- NULL
  if (!is.null(pred_long) && nrow(pred_long)) {
    site_col <- intersect(c("site_id", "site", "row_id"), names(pred_long))[1]
    sp_col <- intersect(c("species", "response_id", "response"), names(pred_long))[1]
    pred_col <- intersect(c("predicted", "predicted_mean", "mean"), names(pred_long))[1]
    obs_col <- intersect(c("observed", "Y", "y"), names(pred_long))[1]
    if (!is.na(site_col) && !is.na(sp_col) && !is.na(pred_col)) {
      preds <- data.frame(engine = engine, site_id = as.character(pred_long[[site_col]]),
                          species = as.character(pred_long[[sp_col]]),
                          observed = if (!is.na(obs_col)) suppressWarnings(as.numeric(pred_long[[obs_col]])) else NA_real_,
                          predicted = suppressWarnings(as.numeric(pred_long[[pred_col]])),
                          stringsAsFactors = FALSE)
    }
  }
  if (is.null(preds) || !nrow(preds)) {
    truth <- dat$truth_predictions
    obs <- as.data.frame(as.table(as.matrix(dat$Y)), stringsAsFactors = FALSE)
    names(obs) <- c("site_id", "species", "observed")
    preds <- merge(truth[, c("site_id", "species")], obs, by = c("site_id", "species"), all.x = TRUE)
    preds$engine <- engine
    preds$predicted <- NA_real_
    preds <- preds[, c("engine", "site_id", "species", "observed", "predicted")]
  }
  preds <- merge(preds, dat$truth_predictions, by = c("site_id", "species"), all.x = TRUE)
  preds$residual <- suppressWarnings(as.numeric(preds$observed) - as.numeric(preds$predicted))
  preds$prediction_scale <- switch(engine, GJAM = "GJAM_observation_scale", spOccupancy = "occupancy_probability", "occurrence_probability_or_engine_scale")
  preds <- preds[, c("engine", "site_id", "species", "observed", "predicted", "truth_probability", "residual", "prediction_scale")]
  write.csv(preds, file.path(outdir, "standard", "predictions_site_species.csv"), row.names = FALSE)

  assoc_long <- read_std("associations_long.csv")
  assoc <- NULL
  if (!is.null(assoc_long) && nrow(assoc_long) && "estimate" %in% names(assoc_long)) {
    s1 <- intersect(c("species_i", "response_1", "response1"), names(assoc_long))[1]
    s2 <- intersect(c("species_j", "response_2", "response2"), names(assoc_long))[1]
    if (!is.na(s1) && !is.na(s2)) {
      assoc <- data.frame(engine = engine, species_i = as.character(assoc_long[[s1]]),
                          species_j = as.character(assoc_long[[s2]]),
                          estimate = suppressWarnings(as.numeric(assoc_long$estimate)),
                          lower = if ("lower" %in% names(assoc_long)) suppressWarnings(as.numeric(assoc_long$lower)) else NA_real_,
                          upper = if ("upper" %in% names(assoc_long)) suppressWarnings(as.numeric(assoc_long$upper)) else NA_real_,
                          association_type = if ("association_type" %in% names(assoc_long)) as.character(assoc_long$association_type) else "engine_association",
                          scale = paste0(engine, "_association_scale"),
                          comparable = identical(status, "fitted"),
                          note = paste0("Associations are not numerically interchangeable across engines; status=", status),
                          stringsAsFactors = FALSE)
    }
  }
  if (is.null(assoc) || !nrow(assoc)) {
    assoc <- dat$truth_associations[, c("species_i", "species_j")]
    assoc$engine <- engine
    assoc$estimate <- NA_real_
    assoc$lower <- NA_real_
    assoc$upper <- NA_real_
    assoc$association_type <- "not_exported"
    assoc$scale <- "not_fitted_or_not_exported"
    assoc$comparable <- FALSE
    assoc$note <- paste0("No fitted association table exported; status=", status)
    assoc <- assoc[, c("engine", "species_i", "species_j", "estimate", "lower", "upper", "association_type", "scale", "comparable", "note")]
  }
  write.csv(assoc, file.path(outdir, "standard", "associations_species_species.csv"), row.names = FALSE)

  if (!file.exists(file.path(outdir, "standard", "diagnostics_long.csv"))) {
    write.csv(data.frame(engine = engine, diagnostic = "status", value = status, stringsAsFactors = FALSE),
              file.path(outdir, "standard", "diagnostics_long.csv"), row.names = FALSE)
  }
  if (!file.exists(file.path(outdir, "standard", "fit_metrics.csv"))) {
    write.csv(data.frame(engine = engine, metric = "not_exported", response_id = NA_character_,
                         value = NA_real_, notes = status, stringsAsFactors = FALSE),
              file.path(outdir, "standard", "fit_metrics.csv"), row.names = FALSE)
  }
  write.csv(data.frame(file = list.files(outdir, recursive = TRUE), stringsAsFactors = FALSE),
            file.path(outdir, "standard", "output_manifest.csv"), row.names = FALSE)
  invisible(TRUE)
}

run_engine <- function(engine, master_dir, dat, app_env, log_fun = message) {
  engine_dir <- file.path(master_dir, "03_engine_outputs", engine)
  ensure_dirs(engine_dir, engine)
  start <- Sys.time()
  log_fun(paste("Starting", engine))
  write_engine_inputs(engine_dir, dat, engine)
  status <- NULL
  err <- NULL
  tryCatch({
    if (identical(engine, "Hmsc")) {
      cfg <- hmsc_cfg()
      yaml::write_yaml(cfg, file.path(engine_dir, "used_config.yml"))
      write_csv_pair(dat$traits, engine_dir, "TrData.csv")
      write_csv_pair(dat$study, engine_dir, "studyDesign.csv")
      write_csv_pair(dat$coords, engine_dir, "coordinates.csv")
      check <- app_env$validate_hmsc(dat$Y, dat$X, dat$traits, dat$study, dat$coords,
                                     distr = cfg$model$distr, XFormula = cfg$model$XFormula,
                                     TrFormula = cfg$model$TrFormula, use_traits = TRUE,
                                     use_phylogeny = TRUE, random_mode = cfg$model$random_mode,
                                     random_effect_column = cfg$model$random_effect_column,
                                     spatial_method = cfg$model$spatial_method,
                                     nNeighbours = cfg$model$nNeighbours,
                                     lon_col = cfg$model$lon_col, lat_col = cfg$model$lat_col,
                                     samples = cfg$mcmc$samples, transient = cfg$mcmc$transient,
                                     thin = cfg$mcmc$thin, nChains = cfg$mcmc$nChains,
                                     nParallel = cfg$mcmc$nParallel, nfMin = cfg$model$nfMin,
                                     nfMax = cfg$model$nfMax,
                                     nfolds = cfg$outputs$predictions$nfolds,
                                     partition_column = cfg$outputs$predictions$partition_column)
      app_env$write_data_check_messages(engine_dir, check$messages)
      if (!isTRUE(check$ok)) {
        write_status(engine_dir, engine, "check_failed", warnings = check$messages,
                     errors = "Hmsc universal benchmark data/settings check failed.")
      } else {
        res <- app_env$run_hmsc_s1s7_pipeline(engine_dir, cfg, dat$Y, dat$X, dat$traits, dat$study, dat$coords, log_fun = log_fun)
        status <- read_status(engine_dir, engine)
        if (is.null(status$status) || identical(status$status, "missing_status")) {
          write_status(engine_dir, engine, res$status %||% "fit_failed",
                       warnings = res$warnings %||% character(), errors = res$errors %||% character())
        }
      }
    } else if (identical(engine, "Hmsc-HPC")) {
      cfg <- hmschpc_cfg(app_env)
      app_env$run_hmschpc_workflow(engine_dir, cfg, dat$Y, dat$X, Tr = dat$traits,
                                   study = dat$study, coord = dat$coords,
                                   phylo_cov = dat$phylo_cov, newdata = dat$newdata,
                                   log_fun = log_fun)
    } else if (identical(engine, "jSDM")) {
      cfg <- jsdm_cfg()
      yaml::write_yaml(cfg, file.path(engine_dir, "used_config.yml"))
      app_env$write_data_check_messages(engine_dir, "Universal benchmark jSDM data generated from shared latent truth.")
      app_env$write_jsdm_reproducible_script(engine_dir)
      code <- run_rscript(file.path(engine_dir, "reproducible_script", "run_this_jSDM_analysis.R"), engine_dir, engine)
      if (!identical(as.integer(code), 0L) && !(read_status(engine_dir, engine)$status %in% c("fitted", "fit_failed", "check_failed"))) {
        write_status(engine_dir, engine, "fit_failed", errors = paste0("jSDM Rscript exited with status ", code, "."))
      }
    } else if (identical(engine, "GJAM")) {
      cfg <- gjam_cfg(colnames(dat$Y))
      yaml::write_yaml(cfg, file.path(engine_dir, "used_config.yml"))
      type_tab <- data.frame(response = colnames(dat$Y), typeName = "PA", stringsAsFactors = FALSE)
      write_csv_pair(type_tab, engine_dir, "typeNames.csv", row.names = FALSE)
      write_csv_pair(dat$traits, engine_dir, "specByTrait.csv")
      write_csv_pair(data.frame(trait = names(dat$traits), typeName = c("CAT", "CON", "CON"), stringsAsFactors = FALSE),
                     engine_dir, "traitTypes.csv", row.names = FALSE)
      app_env$write_data_check_messages(engine_dir, "Universal benchmark GJAM PA data generated from shared latent truth.")
      app_env$write_gjam_reproducible_script(engine_dir)
      code <- run_rscript(file.path(engine_dir, "reproducible_script", "run_this_GJAM_analysis.R"), engine_dir, engine)
      if (!identical(as.integer(code), 0L) && !(read_status(engine_dir, engine)$status %in% c("fitted", "fit_failed", "check_failed"))) {
        write_status(engine_dir, engine, "fit_failed", errors = paste0("GJAM Rscript exited with status ", code, "."))
      }
    } else if (identical(engine, "spOccupancy")) {
      cfg <- spocc_cfg()
      yaml::write_yaml(cfg, file.path(engine_dir, "used_config.yml"))
      write_csv_pair(dat$det, engine_dir, "y.csv")
      write_csv_pair(dat$occ_covs, engine_dir, "occ.covs.csv")
      write_csv_pair(dat$det_cov, engine_dir, "det.covs.csv")
      write_csv_pair(dat$coords, engine_dir, "coords.csv")
      write_csv_pair(data.frame(species = colnames(dat$Y), stringsAsFactors = FALSE),
                     engine_dir, "species.csv", row.names = FALSE)
      app_env$write_data_check_messages(engine_dir, "Universal benchmark spOccupancy replicated detection data generated from shared latent truth.")
      app_env$write_spoccupancy_reproducible_script(engine_dir)
      code <- run_rscript(file.path(engine_dir, "reproducible_script", "run_this_spOccupancy_analysis.R"), engine_dir, engine)
      if (!identical(as.integer(code), 0L) && !(read_status(engine_dir, engine)$status %in% c("fitted", "fit_failed", "check_failed"))) {
        write_status(engine_dir, engine, "fit_failed", errors = paste0("spOccupancy Rscript exited with status ", code, "."))
      }
    } else if (identical(engine, "sjSDM")) {
      cfg <- sjsdm_cfg()
      yaml::write_yaml(cfg, file.path(engine_dir, "used_config.yml"))
      write_csv_pair(dat$X, engine_dir, "env.csv")
      write_csv_pair(dat$coords, engine_dir, "spatial.csv")
      write_csv_pair(dat$coords[seq_len(nrow(dat$newdata)), , drop = FALSE], engine_dir, "new_spatial.csv")
      app_env$write_data_check_messages(engine_dir, "Universal benchmark sjSDM data generated from shared latent truth.")
      app_env$write_sjsdm_reproducible_script(engine_dir)
      code <- run_rscript(file.path(engine_dir, "reproducible_script", "run_this_sjSDM_analysis.R"), engine_dir, engine)
      if (!identical(as.integer(code), 0L) && !(read_status(engine_dir, engine)$status %in% c("fitted", "fit_failed", "check_failed"))) {
        write_status(engine_dir, engine, "fit_failed", errors = paste0("sjSDM Rscript exited with status ", code, "."))
      }
    } else if (identical(engine, "boral")) {
      cfg <- boral_cfg()
      yaml::write_yaml(cfg, file.path(engine_dir, "used_config.yml"))
      app_env$write_data_check_messages(engine_dir, "Universal benchmark boral binomial data generated from shared latent truth.")
      app_env$write_boral_reproducible_script(engine_dir)
      code <- run_rscript(file.path(engine_dir, "reproducible_script", "run_this_boral_analysis.R"), engine_dir, engine)
      if (!identical(as.integer(code), 0L) && !(read_status(engine_dir, engine)$status %in% c("fitted", "fit_failed", "check_failed", "model_defined"))) {
        write_status(engine_dir, engine, "fit_failed", errors = paste0("boral Rscript exited with status ", code, "."))
      }
    } else {
      write_status(engine_dir, engine, "check_failed", errors = paste0("Unknown universal benchmark engine: ", engine))
    }
  }, error = function(e) {
    err <<- conditionMessage(e)
    write_status(engine_dir, engine, "fit_failed", errors = err,
                 runtime_seconds = round(as.numeric(difftime(Sys.time(), start, units = "secs")), 2))
  })
  if (is.null(read_status(engine_dir, engine)$status) || identical(read_status(engine_dir, engine)$status, "missing_status")) {
    write_status(engine_dir, engine, if (is.null(err)) "fit_failed" else "fit_failed", errors = err %||% "Engine finished without writing engine_status.json.")
  }
  standardize_engine_synonyms(engine_dir, engine, dat)
  fill_empty_dirs(engine_dir, engine)
  zipfile <- tryCatch(make_zip_file(engine_dir), error = function(e) {
    writeLines(conditionMessage(e), file.path(engine_dir, "diagnostics", "zip_error.txt"))
    NA_character_
  })
  st <- read_status(engine_dir, engine)
  runtime <- round(as.numeric(difftime(Sys.time(), start, units = "secs")), 2)
  log_fun(paste("Finished", engine, "status", st$status %||% "unknown"))
  data.frame(engine = engine, status = st$status %||% "unknown",
             runtime_seconds = runtime, output_folder = normalizePath(engine_dir, winslash = "/", mustWork = FALSE),
             zip = zipfile, zip_size_kb = if (!is.na(zipfile) && file.exists(zipfile)) round(file.info(zipfile)$size / 1024, 1) else NA_real_,
             file_count = length(list.files(engine_dir, recursive = TRUE)),
             diagnostics = file.path(engine_dir, "diagnostics", "engine_status.json"),
             stringsAsFactors = FALSE)
}

wide_estimates <- function(dat, key_cols, value_col = "estimate", engine_col = "engine") {
  if (is.null(dat) || !nrow(dat)) return(data.frame())
  dat$key <- do.call(paste, c(dat[key_cols], sep = " | "))
  engines <- unique(dat[[engine_col]])
  keys <- unique(dat$key)
  out <- data.frame(key = keys, stringsAsFactors = FALSE)
  for (e in engines) {
    sub <- dat[dat[[engine_col]] == e, , drop = FALSE]
    vals <- tapply(sub[[value_col]], sub$key, function(x) x[which.max(is.finite(x))] %||% x[1])
    out[[e]] <- as.numeric(vals[out$key])
  }
  parts <- strsplit(out$key, " \\| ")
  for (i in seq_along(key_cols)) out[[key_cols[i]]] <- vapply(parts, function(z) z[[i]] %||% "", character(1))
  out[, c(key_cols, setdiff(names(out), c("key", key_cols))), drop = FALSE]
}

collect_master_results <- function(master_dir, engines, dat, statuses) {
  dirs <- file.path(master_dir, "03_engine_outputs", engines)
  read_engine <- function(file) {
    rows <- lapply(seq_along(engines), function(i) {
      p <- file.path(dirs[i], "standard", file)
      if (!file.exists(p)) return(NULL)
      x <- tryCatch(read.csv(p, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(x)) return(NULL)
      if (!("engine" %in% names(x))) x$engine <- engines[i]
      x
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) {
      data.frame()
    } else {
      all_cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
      rows <- lapply(rows, function(x) {
        missing <- setdiff(all_cols, names(x))
        for (col in missing) x[[col]] <- NA
        x[, all_cols, drop = FALSE]
      })
      do.call(rbind, rows)
    }
  }
  out4 <- file.path(master_dir, "04_unified_standard_results")
  out5 <- file.path(master_dir, "05_truth_comparison")
  out6 <- file.path(master_dir, "06_compare_models")
  out7 <- file.path(master_dir, "07_reports")
  for (d in c(out4, out5, out6, out7)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  ensure_cols <- function(df, cols) {
    if (is.null(df)) df <- data.frame()
    for (col in cols) {
      if (!(col %in% names(df))) df[[col]] <- if (nrow(df)) NA else character()
    }
    df
  }
  blank_metric <- function(col_name) {
    out <- data.frame(engine = engines, stringsAsFactors = FALSE)
    out[[col_name]] <- NA_real_
    out
  }
  numeric_mean <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!length(x) || !any(is.finite(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
  }
  numeric_rmse <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!length(x) || !any(is.finite(x))) return(NA_real_)
    sqrt(mean(x^2, na.rm = TRUE))
  }

  effects <- ensure_cols(read_engine("effects_species_environment.csv"),
                         c("engine", "species", "predictor", "estimate", "truth_effect", "comparable", "note"))
  preds <- ensure_cols(read_engine("predictions_site_species.csv"),
                       c("engine", "site_id", "species", "predicted", "truth_probability", "comparable", "note"))
  assocs <- ensure_cols(read_engine("associations_species_species.csv"),
                        c("engine", "species_i", "species_j", "estimate", "true_association",
                          "association_type", "scale", "comparable", "note"))
  fit <- ensure_cols(read_engine("fit_metrics.csv"), c("engine", "metric", "value", "notes"))
  diag <- ensure_cols(read_engine("diagnostics_long.csv"), c("engine", "status", "message"))
  manifests <- lapply(seq_along(engines), function(i) {
    files <- list.files(dirs[i], recursive = TRUE)
    data.frame(engine = engines[i], file = files, stringsAsFactors = FALSE)
  })
  manifest_all <- do.call(rbind, manifests)

  write.csv(statuses, file.path(out4, "model_status_matrix.csv"), row.names = FALSE)
  write.csv(effects, file.path(out4, "effects_all_models_long.csv"), row.names = FALSE)
  write.csv(wide_estimates(effects, c("species", "predictor")), file.path(out4, "effects_species_environment_wide.csv"), row.names = FALSE)
  write.csv(preds, file.path(out4, "predictions_all_models_long.csv"), row.names = FALSE)
  write.csv(wide_estimates(preds, c("site_id", "species"), "predicted"), file.path(out4, "predictions_site_species_wide.csv"), row.names = FALSE)
  write.csv(fit, file.path(out4, "fit_metrics_all_models.csv"), row.names = FALSE)
  write.csv(assocs, file.path(out4, "associations_all_models_long.csv"), row.names = FALSE)
  write.csv(wide_estimates(assocs, c("species_i", "species_j")), file.path(out4, "associations_species_species_wide.csv"), row.names = FALSE)
  write.csv(diag, file.path(out4, "diagnostics_all_models_long.csv"), row.names = FALSE)
  write.csv(manifest_all, file.path(out4, "output_manifest_all_models.csv"), row.names = FALSE)

  truth_eff <- merge(effects, dat$truth_effects, by = c("species", "predictor"), all.x = TRUE)
  truth_eff$estimate_direction <- sign(suppressWarnings(as.numeric(truth_eff$estimate)))
  truth_eff$truth_direction <- sign(suppressWarnings(as.numeric(truth_eff$true_effect)))
  truth_eff$direction_agrees <- is.finite(truth_eff$estimate_direction) & truth_eff$estimate_direction == truth_eff$truth_direction
  write.csv(truth_eff, file.path(out5, "truth_vs_estimated_effects.csv"), row.names = FALSE)
  direction_summary <- blank_metric("effect_direction_agreement")
  if (nrow(truth_eff)) {
    direction_summary$effect_direction_agreement <- vapply(direction_summary$engine, function(e) {
      x <- truth_eff$direction_agrees[truth_eff$engine == e]
      if (!length(x)) return(NA_real_)
      mean(x, na.rm = TRUE)
    }, numeric(1))
  }
  write.csv(direction_summary, file.path(out5, "effect_direction_agreement.csv"), row.names = FALSE)
  truth_pred <- merge(preds, dat$truth_predictions, by = c("site_id", "species"),
                      all.x = TRUE, suffixes = c("", ".truth"))
  if ("truth_probability.truth" %in% names(truth_pred)) {
    truth_pred$truth_probability <- ifelse(is.na(suppressWarnings(as.numeric(truth_pred$truth_probability))),
                                           truth_pred$truth_probability.truth,
                                           truth_pred$truth_probability)
    truth_pred$truth_probability.truth <- NULL
  }
  truth_pred$error_to_truth <- suppressWarnings(as.numeric(truth_pred$predicted) - as.numeric(truth_pred$truth_probability))
  write.csv(truth_pred, file.path(out5, "truth_vs_predictions.csv"), row.names = FALSE)
  pred_metrics <- data.frame(engine = engines,
                             rmse_to_truth_probability = NA_real_,
                             calibration_intercept = NA_real_,
                             calibration_slope = NA_real_,
                             pearson_correlation = NA_real_,
                             spearman_rank_correlation = NA_real_,
                             n_predictions = 0L,
                             stringsAsFactors = FALSE)
  if (nrow(truth_pred)) {
    pred_metrics$rmse_to_truth_probability <- vapply(pred_metrics$engine, function(e) {
      numeric_rmse(truth_pred$error_to_truth[truth_pred$engine == e])
    }, numeric(1))
    pred_metrics$calibration_intercept <- vapply(pred_metrics$engine, function(e) {
      idx <- truth_pred$engine == e
      pred <- suppressWarnings(as.numeric(truth_pred$predicted[idx]))
      truth <- suppressWarnings(as.numeric(truth_pred$truth_probability[idx]))
      ok <- is.finite(pred) & is.finite(truth)
      if (sum(ok) < 3L || length(unique(pred[ok])) < 2L) return(NA_real_)
      unname(coef(lm(truth[ok] ~ pred[ok]))[1])
    }, numeric(1))
    pred_metrics$calibration_slope <- vapply(pred_metrics$engine, function(e) {
      idx <- truth_pred$engine == e
      pred <- suppressWarnings(as.numeric(truth_pred$predicted[idx]))
      truth <- suppressWarnings(as.numeric(truth_pred$truth_probability[idx]))
      ok <- is.finite(pred) & is.finite(truth)
      if (sum(ok) < 3L || length(unique(pred[ok])) < 2L) return(NA_real_)
      unname(coef(lm(truth[ok] ~ pred[ok]))[2])
    }, numeric(1))
    pred_metrics$pearson_correlation <- vapply(pred_metrics$engine, function(e) {
      idx <- truth_pred$engine == e
      pred <- suppressWarnings(as.numeric(truth_pred$predicted[idx]))
      truth <- suppressWarnings(as.numeric(truth_pred$truth_probability[idx]))
      ok <- is.finite(pred) & is.finite(truth)
      if (sum(ok) < 3L || length(unique(pred[ok])) < 2L || length(unique(truth[ok])) < 2L) return(NA_real_)
      suppressWarnings(cor(pred[ok], truth[ok], method = "pearson"))
    }, numeric(1))
    pred_metrics$spearman_rank_correlation <- vapply(pred_metrics$engine, function(e) {
      idx <- truth_pred$engine == e
      pred <- suppressWarnings(as.numeric(truth_pred$predicted[idx]))
      truth <- suppressWarnings(as.numeric(truth_pred$truth_probability[idx]))
      ok <- is.finite(pred) & is.finite(truth)
      if (sum(ok) < 3L || length(unique(pred[ok])) < 2L || length(unique(truth[ok])) < 2L) return(NA_real_)
      suppressWarnings(cor(pred[ok], truth[ok], method = "spearman"))
    }, numeric(1))
    pred_metrics$n_predictions <- vapply(pred_metrics$engine, function(e) {
      sum(is.finite(truth_pred$error_to_truth[truth_pred$engine == e]))
    }, integer(1))
  }
  write.csv(pred_metrics, file.path(out5, "prediction_metric_by_model.csv"), row.names = FALSE)
  assoc_truth <- merge(assocs, dat$truth_associations, by = c("species_i", "species_j"),
                       all.x = TRUE, suffixes = c("", ".truth"))
  if ("true_association.truth" %in% names(assoc_truth)) {
    assoc_truth$true_association <- ifelse(is.na(suppressWarnings(as.numeric(assoc_truth$true_association))),
                                           assoc_truth$true_association.truth,
                                           assoc_truth$true_association)
    assoc_truth$true_association.truth <- NULL
  }
  assoc_metrics <- data.frame(engine = engines, estimate = NA_real_, true_association = NA_real_,
                              stringsAsFactors = FALSE)
  if (nrow(assoc_truth)) {
    assoc_metrics$estimate <- vapply(assoc_metrics$engine, function(e) {
      numeric_mean(assoc_truth$estimate[assoc_truth$engine == e])
    }, numeric(1))
    assoc_metrics$true_association <- vapply(assoc_metrics$engine, function(e) {
      numeric_mean(assoc_truth$true_association[assoc_truth$engine == e])
    }, numeric(1))
  }
  assoc_metrics$note <- "Mean estimates are descriptive only; association parameters are not directly interchangeable across engines."
  write.csv(assoc_truth, file.path(out5, "association_truth_joined.csv"), row.names = FALSE)
  write.csv(assoc_metrics, file.path(out5, "association_recovery_summary.csv"), row.names = FALSE)

  comparable <- data.frame(
    comparison = c("status", "output_contract", "runtime", "prediction_rmse_to_truth",
                   "effect_direction_agreement", "raw_association_parameter"),
    comparable = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
    note = c("Workflow statuses are directly comparable.",
             "Presence of required folders, standard tables, reports, diagnostics and ZIPs is directly comparable.",
             "Runtime is a practical software cost, not a statistical quality measure.",
             "Comparable in this benchmark because all engines derive from the same occurrence truth and use matched IDs.",
             "Comparable as sign agreement only; effect scales differ by engine.",
             "Hmsc Omega, Hmsc-HPC Eta/Lambda, jSDM residual correlations, GJAM corMu/sigMu, spOccupancy latent factors, sjSDM covariance and boral residual correlations are not the same parameter."),
    stringsAsFactors = FALSE
  )
  write.csv(comparable, file.path(out6, "comparable_results_matrix.csv"), row.names = FALSE)
  write.csv(comparable[!comparable$comparable, , drop = FALSE],
            file.path(out6, "non_comparable_results_notes.csv"), row.names = FALSE)
  compare_summary <- merge(statuses, pred_metrics, by = "engine", all.x = TRUE)
  compare_summary <- merge(compare_summary, direction_summary, by = "engine", all.x = TRUE)
  write.csv(compare_summary, file.path(out6, "compare_models_summary.csv"), row.names = FALSE)
  html_table <- function(x) paste(capture.output(print(utils::head(x, 20))), collapse = "\n")
  compare_html <- c("<!doctype html><html><head><meta charset='utf-8'><title>Universal Benchmark Compare</title></head><body>",
                    "<h1>JSDM Studio universal benchmark comparison</h1>",
                    "<h2>Status matrix</h2><pre>", html_table(statuses), "</pre>",
                    "<h2>Prediction metrics</h2><pre>", html_table(pred_metrics), "</pre>",
                    "<h2>Effect direction agreement</h2><pre>", html_table(direction_summary), "</pre>",
                    "<p>Raw association parameters are reported with notes and should not be treated as identical across engines.</p>",
                    "</body></html>")
  writeLines(compare_html, file.path(out6, "compare_models_report.html"), useBytes = TRUE)

  methods <- c(
    "# Universal Benchmark Methods",
    "",
    "JSDM Studio generates one latent ecological truth and derives engine-specific input files with matched site, species and predictor identifiers.",
    "Occurrence responses are used for Hmsc, Hmsc-HPC, jSDM, GJAM, sjSDM and boral. spOccupancy receives replicated detection-nondetection data generated from the same latent occupancy state.",
    "The benchmark is a software and reproducibility benchmark; quick sampler settings are not publication-quality inference settings."
  )
  writeLines(methods, file.path(out7, "benchmark_methods_summary.md"), useBytes = TRUE)
  master_md <- c(
    "# JSDM Studio Universal Benchmark Report",
    "",
    paste0("Generated: ", Sys.time()),
    "",
    "## Engine Status",
    paste(capture.output(print(statuses[, c("engine", "status", "runtime_seconds", "zip_size_kb", "file_count")], row.names = FALSE)), collapse = "\n"),
    "",
    "## Comparable Results",
    paste(capture.output(print(comparable, row.names = FALSE)), collapse = "\n"),
    "",
    "## Diagnostics",
    "Each engine output folder contains diagnostics/engine_status.json, diagnostics/data_check_messages.csv and diagnostics/session_info.txt when the adapter reached the output-contract stage."
  )
  writeLines(master_md, file.path(out7, "master_report.md"), useBytes = TRUE)
  master_html <- c("<!doctype html><html><head><meta charset='utf-8'><title>JSDM Studio Universal Benchmark</title>",
                   "<style>body{font-family:Arial,sans-serif;margin:32px;line-height:1.5;color:#172033}pre{background:#f6f8fb;padding:14px;border-radius:8px;overflow:auto}</style></head><body>",
                   "<h1>JSDM Studio Universal Benchmark Report</h1>",
                   paste0("<p>Generated: ", Sys.time(), "</p>"),
                   "<h2>Engine Status</h2><pre>", html_table(statuses), "</pre>",
                   "<h2>Comparable Results</h2><pre>", html_table(comparable), "</pre>",
                   "<h2>Prediction Metrics</h2><pre>", html_table(pred_metrics), "</pre>",
                   "<h2>Effect Direction Agreement</h2><pre>", html_table(direction_summary), "</pre>",
                   "<p>Association tables are included with association_type, scale, comparable and note columns. They are not treated as a single shared statistical parameter.</p>",
                   "</body></html>")
  writeLines(master_html, file.path(out7, "master_report.html"), useBytes = TRUE)
  writeLines(capture.output(utils::sessionInfo()), file.path(out7, "session_info.txt"))
  write.csv(data.frame(file = list.files(master_dir, recursive = TRUE), stringsAsFactors = FALSE),
            file.path(out7, "master_output_manifest.csv"), row.names = FALSE)
  invisible(TRUE)
}

write_master_scripts <- function(master_dir, cfg) {
  sdir <- file.path(master_dir, "08_reproducible_scripts")
  cdir <- file.path(master_dir, "09_configs")
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  runner_rel <- normalizePath(file.path(app_dir, "workflow_scripts", "universal_benchmark_runner.R"), winslash = "/", mustWork = FALSE)
  writeLines(c("# Generate benchmark data", sprintf("system2(file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript'), c(%s, '--action=generate'))", deparse(runner_rel))),
             file.path(sdir, "generate_benchmark_data.R"))
  writeLines(c("# Run all benchmark engines", sprintf("system2(file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript'), c(%s, '--action=all'))", deparse(runner_rel))),
             file.path(sdir, "run_all_engines.R"))
  writeLines(c("# Rebuild comparison tables for an existing master folder",
               sprintf("system2(file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript'), c(%s, '--action=compare', paste0('--master_dir=', getwd())))", deparse(runner_rel))),
             file.path(sdir, "compare_benchmark_results.R"))
  writeLines(c("# Rebuild the all-model ZIP for an existing master folder",
               sprintf("system2(file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript'), c(%s, '--action=zip', paste0('--master_dir=', getwd())))", deparse(runner_rel))),
             file.path(sdir, "rebuild_all_model_zip.R"))
  if (safe_require("yaml")) {
    yaml::write_yaml(cfg, file.path(cdir, "universal_benchmark_config.yml"))
    yaml::write_yaml(hmsc_cfg(), file.path(cdir, "Hmsc_config.yml"))
    yaml::write_yaml(jsdm_cfg(), file.path(cdir, "jSDM_config.yml"))
    yaml::write_yaml(gjam_cfg(sprintf("sp_%02d", seq_len(cfg$n_species))), file.path(cdir, "GJAM_config.yml"))
    yaml::write_yaml(spocc_cfg(), file.path(cdir, "spOccupancy_config.yml"))
    yaml::write_yaml(sjsdm_cfg(), file.path(cdir, "sjSDM_config.yml"))
    yaml::write_yaml(boral_cfg(), file.path(cdir, "boral_config.yml"))
  }
}

run_universal_benchmark <- function(args) {
  run_root <- file.path(app_dir, "output", paste0("universal_benchmark_", timestamp_id()))
  master_dir <- file.path(run_root, "master")
  for (d in c("01_benchmark_inputs", "02_truth", "03_engine_outputs", "04_unified_standard_results",
              "05_truth_comparison", "06_compare_models", "07_reports", "08_reproducible_scripts", "09_configs")) {
    dir.create(file.path(master_dir, d), recursive = TRUE, showWarnings = FALSE)
  }
  log_file <- file.path(master_dir, "07_reports", "run_log.txt")
  log_fun <- function(txt) {
    line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", txt)
    cat(line, "\n")
    cat(line, "\n", file = log_file, append = TRUE)
  }
  input_dir <- trimws(as.character(args$input_dir %||% ""))
  if (nzchar(input_dir)) {
    input_dir <- normalizePath(input_dir, winslash = "/", mustWork = FALSE)
    dat <- load_benchmark_data(input_dir, allow_generate = FALSE)
    data_source <- "user_or_example_input_dir"
    log_fun(paste0("Using supplied benchmark input directory: ", dat$dir))
  } else {
    dat <- generate_benchmark_data(file.path(app_dir, "examples", "universal_benchmark"),
                                   args$n_sites, args$n_species, args$n_visits, args$seed)
    dat <- load_benchmark_data(dat$dir, allow_generate = FALSE)
    data_source <- "generated_from_latent_truth"
    log_fun(paste0("Generated benchmark input directory: ", dat$dir))
  }
  copy_dir_contents(dat$dir, file.path(master_dir, "01_benchmark_inputs"))
  copy_dir_contents(file.path(dat$dir, "truth"), file.path(master_dir, "02_truth"))
  cfg <- list(seed = args$seed, n_sites = nrow(dat$Y), n_species = ncol(dat$Y),
              n_visits = args$n_visits, engines = args$engines, data_source = data_source,
              input_dir = dat$dir)
  write_master_scripts(master_dir, cfg)
  preflight <- dependency_preflight(master_dir)
  app_env <- source_app_env()
  if (safe_require("yaml")) yaml::write_yaml(hmschpc_cfg(app_env), file.path(master_dir, "09_configs", "Hmsc-HPC_config.yml"))
  statuses <- do.call(rbind, lapply(args$engines, function(engine) run_engine(engine, master_dir, dat, app_env, log_fun)))
  collect_master_results(master_dir, args$engines, dat, statuses)
  master_zip <- make_zip_file(master_dir, file.path(run_root, "JSDMStudio_universal_benchmark_ALL_MODELS.zip"))
  sha <- if (file.exists(master_zip)) tools::sha256sum(master_zip)[[1]] else NA_character_
  summary <- data.frame(master_dir = normalizePath(master_dir, winslash = "/", mustWork = FALSE),
                        master_zip = normalizePath(master_zip, winslash = "/", mustWork = FALSE),
                        sha256 = sha,
                        size_mb = round(file.info(master_zip)$size / 1024 / 1024, 2),
                        data_source = data_source,
                        input_dir = dat$dir,
                        engines = paste(args$engines, collapse = ","),
                        stringsAsFactors = FALSE)
  write.csv(summary, file.path(run_root, "universal_benchmark_summary.csv"), row.names = FALSE)
  if (safe_require("jsonlite")) {
    writeLines(jsonlite::toJSON(as.list(summary[1, ]), pretty = TRUE, auto_unbox = TRUE),
               file.path(app_dir, "output", "universal_benchmark_last.json"))
  }
  print(statuses)
  cat("MASTER_DIR=", summary$master_dir, "\n", sep = "")
  cat("MASTER_ZIP=", summary$master_zip, "\n", sep = "")
  cat("SHA256=", summary$sha256, "\n", sep = "")
  summary
}

compare_existing_master <- function(master_dir) {
  if (!dir.exists(master_dir)) stop("master_dir does not exist: ", master_dir, call. = FALSE)
  input_dir <- file.path(master_dir, "01_benchmark_inputs")
  dat <- load_benchmark_data(if (dir.exists(input_dir)) input_dir else file.path(app_dir, "examples", "universal_benchmark"),
                             allow_generate = FALSE)
  engines <- list.dirs(file.path(master_dir, "03_engine_outputs"), recursive = FALSE, full.names = FALSE)
  statuses <- do.call(rbind, lapply(engines, function(e) {
    d <- file.path(master_dir, "03_engine_outputs", e)
    st <- read_status(d, e)
    data.frame(engine = e, status = st$status %||% "unknown", runtime_seconds = NA_real_,
               output_folder = normalizePath(d, winslash = "/", mustWork = FALSE),
               zip = paste0(normalizePath(d, winslash = "/", mustWork = FALSE), ".zip"),
               zip_size_kb = if (file.exists(paste0(d, ".zip"))) round(file.info(paste0(d, ".zip"))$size / 1024, 1) else NA_real_,
               file_count = length(list.files(d, recursive = TRUE)),
               diagnostics = file.path(d, "diagnostics", "engine_status.json"),
               stringsAsFactors = FALSE)
  }))
  collect_master_results(master_dir, engines, dat, statuses)
  invisible(statuses)
}

main <- function() {
  args <- parse_args()
  if (identical(args$action, "generate")) {
    target_dir <- trimws(as.character(args$target_dir %||% ""))
    if (!nzchar(target_dir)) target_dir <- file.path(app_dir, "examples", "universal_benchmark")
    x <- generate_benchmark_data(target_dir,
                                 args$n_sites, args$n_species, args$n_visits, args$seed)
    cat("BENCHMARK_DATA_DIR=", x$dir, "\n", sep = "")
  } else if (identical(args$action, "preflight")) {
    print(dependency_preflight())
  } else if (identical(args$action, "compare")) {
    compare_existing_master(args$master_dir)
    cat("COMPARE_OK\n")
  } else if (identical(args$action, "zip")) {
    zip <- make_zip_file(args$master_dir, file.path(dirname(args$master_dir), "JSDMStudio_universal_benchmark_ALL_MODELS.zip"))
    cat("MASTER_ZIP=", zip, "\n", sep = "")
    cat("SHA256=", tools::sha256sum(zip)[[1]], "\n", sep = "")
  } else {
    run_universal_benchmark(args)
  }
}

if (!interactive()) main()
