# Universal Benchmark runner for JSDM Studio.
# This script is intentionally command-line runnable so the GUI benchmark can be
# reproduced without Shiny.

options(stringsAsFactors = FALSE)

rscript_clean_env <- function() {
  if (.Platform$OS.type == "windows") {
    Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
    Sys.setenv(PYTHONIOENCODING = "utf-8")
  }
  invisible(TRUE)
}
rscript_clean_env()

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
  out <- list(action = "all", n_sites = 120L, n_species = 12L, n_visits = 3L,
              seed = 20260601L, engines = "Hmsc,Hmsc-HPC,jSDM,GJAM,spOccupancy,sjSDM,boral",
              master_dir = "", input_dir = "", target_dir = "",
              allow_large = "true",
              generating_structure = "probit_latent",
              truth_profile = "recovery", effect_scale = "1.60", association_strength = "1.65",
              residual_sd = "1.05", observation_noise_sd = "0.15",
              hmsc_samples = "500", hmsc_transient = "500", hmsc_thin = "1", hmsc_nChains = "2",
              hmschpc_samples = "500", hmschpc_transient = "500", hmschpc_thin = "1", hmschpc_chains = "1",
              jsdm_burnin = "200", jsdm_mcmc = "300", jsdm_thin = "5",
              gjam_ng = "700", gjam_burnin = "150",
              spocc_n_batch = "30", spocc_batch_length = "10", spocc_n_burn = "100", spocc_n_thin = "5", spocc_n_chains = "1",
              sjsdm_iter = "100", sjsdm_sampling = "1000", sjsdm_step_size = "32",
              sjsdm_learning_rate = "0.010", sjsdm_weight_decay = "0", sjsdm_biotic_lambda = "0.001",
              boral_n_burnin = "250", boral_n_iteration = "750", boral_n_thin = "10")
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
  out$allow_large <- as_bool(out$allow_large, default = TRUE)
  out$engines <- trimws(unlist(strsplit(as.character(out$engines), ",")))
  out$engines <- out$engines[nzchar(out$engines)]
  if (!out$generating_structure %in% c("probit_latent", "logit_occupancy_detection")) {
    stop("Unknown generating_structure: ", out$generating_structure, call. = FALSE)
  }
  out
}

arg_int <- function(args, key, default, min_value = NULL, max_value = NULL) {
  val <- args[[key]] %||% ""
  if (!nzchar(trimws(as.character(val)))) return(as.integer(default))
  out <- suppressWarnings(as.integer(val))
  if (!is.finite(out)) return(as.integer(default))
  if (!is.null(min_value)) out <- max(as.integer(min_value), out)
  if (!is.null(max_value)) out <- min(as.integer(max_value), out)
  as.integer(out)
}

arg_num <- function(args, key, default, min_value = NULL, max_value = NULL) {
  val <- args[[key]] %||% ""
  if (!nzchar(trimws(as.character(val)))) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(val))
  if (!is.finite(out)) return(as.numeric(default))
  if (!is.null(min_value)) out <- max(as.numeric(min_value), out)
  if (!is.null(max_value)) out <- min(as.numeric(max_value), out)
  as.numeric(out)
}

script_path <- {
  args0 <- commandArgs(trailingOnly = FALSE)
  hit <- if (sys.nframe() == 0L) args0[startsWith(args0, "--file=")] else character()
  path <- if (length(hit)) sub("^--file=", "", hit[[1]]) else "workflow_scripts/universal_benchmark_runner.R"
  if (!file.exists(path)) path <- gsub("~+~", " ", path, fixed = TRUE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
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
source(file.path(app_dir, "R", "empirical_benchmark.R"), local = TRUE)

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
  scientific_caveats <- c(
    "A fitted software-test run means the workflow and output contract executed; publication inference still requires adequate sampler settings and convergence diagnostics.",
    "Association parameters are engine-specific and should be compared only where standard tables explicitly mark them as comparable."
  )
  payload <- list(engine = engine, status = status, warnings = as.character(warnings),
                  errors = as.character(errors), runtime_seconds = runtime_seconds,
                  scientific_caveats = scientific_caveats,
                  operational_risks = character(),
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
  write.csv(data.frame(
    engine = engine,
    caveat = c("software_test_vs_publication_inference", "engine_specific_association_scales"),
    control = c("Increase sampler settings and inspect diagnostics before ecological inference.",
                "Use comparable flags and notes in standard association tables."),
    note = scientific_caveats,
    stringsAsFactors = FALSE
  ), file.path(outdir, "diagnostics", "scientific_caveats.csv"), row.names = FALSE)
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

truth_profile_defaults <- function(profile = "recovery") {
  profile <- tolower(trimws(as.character(profile %||% "recovery")))
  if (profile %in% c("quick", "software", "software_check")) {
    return(list(effect_scale = 1.05, association_strength = 0.80,
                residual_sd = 0.55, observation_noise_sd = 0.25))
  }
  if (profile %in% c("association", "association_stress", "assoc")) {
    return(list(effect_scale = 1.15, association_strength = 1.45,
                residual_sd = 1.05, observation_noise_sd = 0.18))
  }
  list(effect_scale = 1.60, association_strength = 1.65,
       residual_sd = 1.05, observation_noise_sd = 0.15)
}

stable_scale <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

make_benchmark_predictors <- function(n_sites) {
  trend <- stable_scale(seq(-1, 1, length.out = n_sites))
  wave <- stable_scale(sin(seq(0, 4 * pi, length.out = n_sites)))
  z1 <- stable_scale(0.20 * trend + rnorm(n_sites))
  z2 <- stable_scale(0.20 * wave + rnorm(n_sites))
  z3 <- stable_scale(0.15 * trend - 0.10 * z1 + rnorm(n_sites))
  z4 <- stable_scale(-0.15 * wave + 0.15 * z2 + rnorm(n_sites))
  substrate <- sample(rep(c("clay", "loam", "rock", "sand"), length.out = n_sites), n_sites)
  data.frame(
    pH = round(6.45 + 0.45 * z1, 3),
    moisture = round(pmin(pmax(stats::plogis(0.15 + 0.85 * z2), 0.05), 0.95), 3),
    canopy = round(pmin(pmax(52 + 18 * z3, 5), 95), 3),
    elevation = round(260 + 85 * z4, 3),
    substrate = substrate,
    stringsAsFactors = FALSE
  )
}

make_effect_truth <- function(n_species, species, predictors, effect_scale = 1.35) {
  theta <- seq(0, 2 * pi, length.out = n_species + 1L)[seq_len(n_species)]
  continuous <- cbind(
    pH = 0.95 * sin(theta),
    moisture = 0.90 * cos(theta),
    canopy = rep(c(-0.85, 0.85), length.out = n_species),
    elevation = 0.75 * sin(2 * theta + 0.35)
  )
  substrate_level <- cbind(
    substrate_clay = rep(c(0.65, -0.35, 0.20), length.out = n_species),
    substrate_loam = rep(c(0.15, 0.40, -0.25, -0.10), length.out = n_species),
    substrate_rock = rep(c(-0.55, 0.30, 0.60, -0.20), length.out = n_species),
    substrate_sand = rep(c(0.00, -0.20, 0.25, 0.45), length.out = n_species)
  )
  substrate_contrast <- cbind(
    substrate_loam = substrate_level[, "substrate_loam"] - substrate_level[, "substrate_clay"],
    substrate_rock = substrate_level[, "substrate_rock"] - substrate_level[, "substrate_clay"],
    substrate_sand = substrate_level[, "substrate_sand"] - substrate_level[, "substrate_clay"]
  )
  beta <- cbind(continuous, substrate_contrast)
  substrate_level_effects <- effect_scale * substrate_level
  dimnames(substrate_level_effects) <- list(species, colnames(substrate_level))
  beta <- effect_scale * beta[, predictors, drop = FALSE]
  dimnames(beta) <- list(species, predictors)
  attr(beta, "substrate_level_effects") <- substrate_level_effects
  beta
}

make_residual_association_truth <- function(n_sites, n_species, site_id, species,
                                            association_strength = 1.15,
                                            residual_sd = 0.85) {
  theta <- seq(0, 2 * pi, length.out = n_species + 1L)[seq_len(n_species)]
  L <- cbind(
    cos(theta),
    sin(theta),
    rep(c(-1, 0, 1), length.out = n_species)
  )
  L <- apply(L, 2, stable_scale)
  L <- association_strength * L / sqrt(ncol(L))
  unique_sd <- 0.45
  corr <- stats::cov2cor(tcrossprod(L) + diag(unique_sd^2, n_species))
  corr <- (corr + t(corr)) / 2
  diag(corr) <- 1
  eig <- eigen(corr, symmetric = TRUE)
  eig$values[eig$values < 1e-6] <- 1e-6
  corr <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  corr <- stats::cov2cor(corr)
  dimnames(corr) <- list(species, species)
  Sigma <- residual_sd^2 * corr
  latent <- matrix(rnorm(n_sites * n_species), n_sites, n_species) %*% chol(Sigma)
  dimnames(latent) <- list(site_id, species)
  list(correlation = corr, covariance = Sigma, latent = latent)
}

make_taxonomy_correlation_reference <- function(taxonomy, species) {
  rank_cols <- intersect(c("kingdom", "phylum", "class", "order", "family", "genus"), names(taxonomy))
  n <- length(species)
  C <- matrix(0, n, n, dimnames = list(species, species))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      C[i, j] <- if (i == j) {
        1
      } else {
        mean(as.character(taxonomy[i, rank_cols, drop = TRUE]) ==
               as.character(taxonomy[j, rank_cols, drop = TRUE]), na.rm = TRUE)
      }
    }
  }
  C <- 0.02 + 0.96 * C
  diag(C) <- 1
  C <- (C + t(C)) / 2
  eig <- eigen(C, symmetric = TRUE)
  eig$values[eig$values < 1e-6] <- 1e-6
  C <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  C <- stats::cov2cor(C)
  dimnames(C) <- list(species, species)
  C
}

calibrate_binary_intercept <- function(linear_predictor, target_prevalence,
                                       link = c("probit", "logit")) {
  link <- match.arg(link)
  target_prevalence <- min(max(as.numeric(target_prevalence), 0.02), 0.98)
  inv_link <- if (identical(link, "probit")) stats::pnorm else stats::plogis
  f <- function(b) mean(inv_link(linear_predictor + b), na.rm = TRUE) - target_prevalence
  stats::uniroot(f, interval = c(-12, 12), extendInt = "yes")$root
}

generate_benchmark_data <- function(target_dir = file.path(app_dir, "examples", "universal_benchmark"),
                                    n_sites = 120L, n_species = 12L, n_visits = 3L, seed = 20260601L,
                                    allow_large = TRUE, truth_profile = "recovery",
                                    generating_structure = "probit_latent",
                                    effect_scale = NA_real_, association_strength = NA_real_,
                                    residual_sd = NA_real_, observation_noise_sd = NA_real_) {
  generating_structure <- match.arg(generating_structure,
                                    c("probit_latent", "logit_occupancy_detection"))
  if (isTRUE(allow_large)) {
    n_sites <- max(30L, min(1000L, as.integer(n_sites)))
    n_species <- max(5L, min(100L, as.integer(n_species)))
  } else {
    n_sites <- max(30L, min(40L, as.integer(n_sites)))
    n_species <- max(5L, min(8L, as.integer(n_species)))
  }
  n_visits <- max(2L, min(5L, as.integer(n_visits)))
  set.seed(seed)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(target_dir, "truth"), recursive = TRUE, showWarnings = FALSE)
  unlink(file.path(target_dir, "phylo_cov.csv"), force = TRUE)
  defaults <- truth_profile_defaults(truth_profile)
  if (!is.finite(effect_scale)) effect_scale <- defaults$effect_scale
  if (!is.finite(association_strength)) association_strength <- defaults$association_strength
  if (!is.finite(residual_sd)) residual_sd <- defaults$residual_sd
  if (!is.finite(observation_noise_sd)) observation_noise_sd <- defaults$observation_noise_sd

  site_id <- sprintf("site_%02d", seq_len(n_sites))
  species <- sprintf("sp_%02d", seq_len(n_species))
  X <- make_benchmark_predictors(n_sites)
  rownames(X) <- site_id
  X_scaled <- data.frame(
    pH = stable_scale(X$pH),
    moisture = stable_scale(X$moisture),
    canopy = stable_scale(X$canopy),
    elevation = stable_scale(X$elevation),
    substrate_loam = as.numeric(X$substrate == "loam"),
    substrate_rock = as.numeric(X$substrate == "rock"),
    substrate_sand = as.numeric(X$substrate == "sand")
  )
  predictors <- names(X_scaled)

  beta <- make_effect_truth(n_species, species, predictors, effect_scale = effect_scale)
  resid_truth <- make_residual_association_truth(n_sites, n_species, site_id, species,
                                                 association_strength = association_strength,
                                                 residual_sd = residual_sd)
  assoc <- resid_truth$correlation
  fixed_eta <- as.matrix(X_scaled) %*% t(beta) + resid_truth$latent
  target_prevalence <- seq(0.30, 0.70, length.out = n_species)
  occurrence_link <- if (identical(generating_structure, "probit_latent")) "probit" else "logit"
  intercept <- vapply(seq_len(n_species), function(j) {
    calibrate_binary_intercept(fixed_eta[, j], target_prevalence[j], link = occurrence_link)
  }, numeric(1))
  eta <- sweep(fixed_eta, 2, intercept, "+")
  colnames(eta) <- species
  rownames(eta) <- site_id
  inverse_link <- if (identical(occurrence_link, "probit")) stats::pnorm else stats::plogis
  prob <- pmin(pmax(inverse_link(eta), 0.02), 0.98)
  Y_occ <- matrix(rbinom(n_sites * n_species, 1, as.vector(prob)), n_sites, n_species,
                  dimnames = list(site_id, species))
  lambda <- pmax(exp(0.4 + eta / 2), 0.05)
  Y_count <- matrix(rpois(n_sites * n_species, lambda = as.vector(lambda)), n_sites, n_species,
                    dimnames = list(site_id, species))
  Y_normal <- eta + matrix(rnorm(n_sites * n_species, 0, observation_noise_sd), n_sites, n_species)
  dimnames(Y_normal) <- list(site_id, species)

  traits <- data.frame(
    life_form = rep(c("grass", "forb", "shrub", "tree"), length.out = n_species),
    height_mm = round(seq(40, 420, length.out = n_species), 1),
    dispersal = round(seq(0.15, 0.95, length.out = n_species), 3),
    stringsAsFactors = FALSE
  )
  rownames(traits) <- species
  taxonomy <- data.frame(
    kingdom = rep("Synthetic", n_species),
    phylum = rep("Benchmarkophyta", n_species),
    class = rep(c("Class_A", "Class_B"), length.out = n_species),
    order = rep(c("Order_1", "Order_2", "Order_3"), length.out = n_species),
    family = rep(sprintf("Family_%02d", seq_len(ceiling(n_species / 2))), each = 2, length.out = n_species),
    genus = paste0("Genus_", seq_len(n_species)),
    species = species,
    stringsAsFactors = FALSE
  )
  rownames(taxonomy) <- species
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
  newick <- make_tree_newick(species)
  phylo_cov <- make_taxonomy_correlation_reference(taxonomy, species)

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

  write_indexed_csv <- function(x, path, id_name) {
    df <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
    ids <- rownames(df)
    if (is.null(ids) || any(!nzchar(ids))) ids <- seq_len(nrow(df))
    out <- cbind(stats::setNames(data.frame(ids, check.names = FALSE), id_name), df)
    write.csv(out, path, row.names = FALSE)
  }

  write_indexed_csv(Y_occ, file.path(target_dir, "Y_occurrence.csv"), "site_id")
  write_indexed_csv(Y_count, file.path(target_dir, "Y_count.csv"), "site_id")
  write_indexed_csv(round(Y_normal, 4), file.path(target_dir, "Y_normal.csv"), "site_id")
  write_indexed_csv(X, file.path(target_dir, "XData.csv"), "site_id")
  write_indexed_csv(traits, file.path(target_dir, "traits.csv"), "species")
  write_indexed_csv(taxonomy, file.path(target_dir, "taxonomy.csv"), "species")
  write_indexed_csv(study, file.path(target_dir, "studyDesign.csv"), "site_id")
  write_indexed_csv(coords, file.path(target_dir, "coordinates.csv"), "site_id")
  writeLines(newick, file.path(target_dir, "phylogeny.nwk"), useBytes = TRUE)
  write_indexed_csv(newdata, file.path(target_dir, "newdata.csv"), "site_id")
  write_indexed_csv(newcoords, file.path(target_dir, "newcoords.csv"), "site_id")
  write.csv(folds, file.path(target_dir, "folds.csv"), row.names = FALSE)
  write_indexed_csv(trial, file.path(target_dir, "trial.size.csv"), "site_id")
  write_indexed_csv(offset, file.path(target_dir, "offset.csv"), "site_id")
  write.csv(rowids, file.path(target_dir, "row.ids.csv"), row.names = FALSE)
  write_indexed_csv(ranefids, file.path(target_dir, "ranef.ids.csv"), "site_id")
  write_indexed_csv(distmat, file.path(target_dir, "distmat.csv"), "site_id")
  write_indexed_csv(det, file.path(target_dir, "spOccupancy_y_detection.csv"), "site_id")
  write_indexed_csv(occ_covs, file.path(target_dir, "occ.covs.csv"), "site_id")
  write_indexed_csv(det_cov, file.path(target_dir, "det.covs.csv"), "site_id")

  true_eff <- as.data.frame(as.table(beta), stringsAsFactors = FALSE)
  names(true_eff) <- c("species", "predictor", "true_effect")
  substrate_level_effects <- attr(beta, "substrate_level_effects")
  true_level_eff <- as.data.frame(as.table(substrate_level_effects), stringsAsFactors = FALSE)
  names(true_level_eff) <- c("species", "predictor", "true_level_effect")
  true_level_eff$note <- "Raw substrate level effects used to derive treatment contrasts; true_environment_effects.csv stores model-identifiable contrasts with clay as the baseline."
  true_assoc <- as.data.frame(as.table(assoc), stringsAsFactors = FALSE)
  names(true_assoc) <- c("species_i", "species_j", "true_association")
  true_assoc <- true_assoc[as.integer(factor(true_assoc$species_i, levels = species)) <
                             as.integer(factor(true_assoc$species_j, levels = species)), , drop = FALSE]
  pred_truth <- as.data.frame(as.table(prob), stringsAsFactors = FALSE)
  names(pred_truth) <- c("site_id", "species", "truth_probability")
  write_indexed_csv(round(eta, 4), file.path(target_dir, "truth", "latent_occurrence.csv"), "site_id")
  write_indexed_csv(round(prob, 6), file.path(target_dir, "truth", "occurrence_probability.csv"), "site_id")
  write.csv(true_eff, file.path(target_dir, "truth", "true_environment_effects.csv"), row.names = FALSE)
  write.csv(true_level_eff, file.path(target_dir, "truth", "true_environment_level_effects.csv"), row.names = FALSE)
  write.csv(true_assoc, file.path(target_dir, "truth", "true_species_associations.csv"), row.names = FALSE)
  write.csv(pred_truth, file.path(target_dir, "truth", "true_predictions.csv"), row.names = FALSE)
  write_indexed_csv(round(stats::cor(as.matrix(X_scaled)), 4), file.path(target_dir, "truth", "predictor_correlation.csv"), "predictor")
  write_indexed_csv(round(resid_truth$covariance, 6), file.path(target_dir, "truth", "residual_association_covariance.csv"), "species")
  write_indexed_csv(phylo_cov, file.path(target_dir, "truth", "derived_phylogeny_correlation_reference.csv"), "species")
  cfg <- list(seed = seed, n_sites = n_sites, n_species = n_species, n_visits = n_visits,
              allow_large = isTRUE(allow_large),
              generating_structure = generating_structure,
              generator_implementation = "JSDM Studio structural simulator; not a native package simulate method",
              truth_profile = truth_profile,
              effect_scale = effect_scale,
              association_strength = association_strength,
              residual_sd = residual_sd,
              observation_noise_sd = observation_noise_sd,
              occurrence_link = occurrence_link,
              site_id = site_id, species = species, predictors = predictors,
              categorical_baseline = "substrate clay",
              note = "All engine-specific inputs are derived from one explicitly named generating structure. Probit latent structure approximates the binary link used in the benchmark Hmsc, Hmsc-HPC, jSDM and sjSDM configurations; logit occupancy-detection structure approximates the benchmark spOccupancy observation process. Neither structure is an exact native simulator for any fitted package. GJAM and boral are cross-structure fits. Residual association truth is generated in the latent residual process; phylogeny/taxonomy relatedness is recorded separately. Categorical substrate truth is exported as model-identifiable treatment contrasts in true_environment_effects.csv and raw level effects in true_environment_level_effects.csv.")
  if (safe_require("yaml")) yaml::write_yaml(cfg, file.path(target_dir, "truth", "data_generation_config.yml"))
  writeLines(c(
    "# Rebuild the JSDM Studio universal benchmark data.",
    "if (.Platform$OS.type == 'windows') Sys.unsetenv(c('LC_ALL','LC_COLLATE','LC_CTYPE','LC_MONETARY','LC_TIME','LANG','LANGUAGE'))",
    "runner <- file.path('..', '..', 'workflow_scripts', 'universal_benchmark_runner.R')",
    "system2(file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript'),",
    paste0("        c(runner, '--action=generate', '--n_sites=", n_sites,
           "', '--n_species=", n_species,
           "', '--n_visits=", n_visits,
           "', '--seed=", seed,
           "', '--allow_large=", tolower(as.character(isTRUE(allow_large))),
           "', '--generating_structure=", generating_structure,
           "', '--truth_profile=", truth_profile,
           "', '--effect_scale=", effect_scale,
           "', '--association_strength=", association_strength,
           "', '--residual_sd=", residual_sd,
           "', '--observation_noise_sd=", observation_noise_sd, "'))")
  ), file.path(target_dir, "truth", "data_generation_script.R"))
  write.csv(data.frame(file = list.files(target_dir, recursive = TRUE), stringsAsFactors = FALSE),
            file.path(target_dir, "benchmark_input_manifest.csv"), row.names = FALSE)
  list(dir = normalizePath(target_dir, winslash = "/", mustWork = FALSE), Y = Y_occ, X = X,
       traits = traits, taxonomy = taxonomy, study = study, coords = coords, phylo_cov = phylo_cov, newick = newick,
       newdata = newdata, newcoords = newcoords, det = det, occ_covs = occ_covs, det_cov = det_cov,
       rowids = rowids, ranefids = ranefids, distmat = distmat, offset = offset, trial = trial,
       truth_effects = true_eff, truth_predictions = pred_truth, truth_associations = true_assoc)
}

benchmark_required_files <- function() {
  c("Y_occurrence.csv", "Y_count.csv", "Y_normal.csv", "XData.csv", "traits.csv",
    "studyDesign.csv", "coordinates.csv", "phylogeny.nwk", "taxonomy.csv",
    "newdata.csv", "newcoords.csv", "folds.csv", "trial.size.csv", "offset.csv",
    "row.ids.csv", "ranef.ids.csv", "distmat.csv", "spOccupancy_y_detection.csv",
    "occ.covs.csv", "det.covs.csv", file.path("truth", "latent_occurrence.csv"),
    file.path("truth", "occurrence_probability.csv"),
    file.path("truth", "true_environment_effects.csv"),
    file.path("truth", "true_species_associations.csv"),
    file.path("truth", "true_predictions.csv"),
    file.path("truth", "derived_phylogeny_correlation_reference.csv"),
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
    taxonomy = read_csv("taxonomy.csv"),
    study = read_csv("studyDesign.csv"),
    coords = read_csv("coordinates.csv"),
    phylo_cov = read_csv(file.path("truth", "derived_phylogeny_correlation_reference.csv")),
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

hmsc_cfg <- function(args = NULL) {
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
      C_file = "", use_XRRR = FALSE, ncRRR = 2L,
      XRRRFormula = "~ .", XRRRScale = TRUE, XRRR_file = "",
      random_level_type = "none", sMethod = "Full", random_N = 1L,
      longlat = FALSE, units_column = "sample", distMat_file = "",
      xData_file = "", sKnot_file = "", nfMin = 1L, nfMax = 4L,
      priors = list(setDefault = TRUE, a1 = NA, b1 = NA, a2 = NA, b2 = NA, alphapw = "")
    ),
    mcmc = list(samples = arg_int(args, "hmsc_samples", 500L, min_value = 1L),
                transient = arg_int(args, "hmsc_transient", 500L, min_value = 0L),
                thin = arg_int(args, "hmsc_thin", 1L, min_value = 1L),
                nChains = arg_int(args, "hmsc_nChains", 2L, min_value = 1L),
                nParallel = 1L,
                verbose = 0L, preset = "Universal reviewer benchmark", initPar = "fixed effects",
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

hmschpc_cfg <- function(app_env, args = NULL) {
  list(
    engine = "Hmsc-HPC",
    inputs = list(Y = "Y.csv", XData = "XData.csv", traits = "traits.csv",
                  studyDesign = "studyDesign.csv", coordinates = NULL,
                  taxonomy = NULL, phylo_tree = "phylo_tree.nwk", newdata = "newdata.csv"),
    model = list(distribution = "probit", XFormula = "~ pH + moisture + canopy + elevation + C(substrate)",
                 use_traits = TRUE, trait_formula = "~ height_mm + dispersal + C(life_form)",
                 phylogeny_mode = "newick"),
    random_effects = list(mode = "iid", name = "plot_iid", column = "plot",
                          coord_x = "x", coord_y = "y", x_formula = "",
                          nf = 1L, nfMin = 1L, nfMax = 4L, alpha = 1),
    sampler = list(run_sampler = TRUE,
                   samples = arg_int(args, "hmschpc_samples", 500L, min_value = 1L),
                   transient = arg_int(args, "hmschpc_transient", 500L, min_value = 0L),
                   thin = arg_int(args, "hmschpc_thin", 1L, min_value = 1L),
                   chains = arg_int(args, "hmschpc_chains", 1L, min_value = 1L),
                   chains_to_run = NULL, verbose = 5L, seed = 20260601L, precision = 32L,
                   truncated_normal_library = "tf", hmcleapfrog = 5L, hmcthin = 0L,
                   update_beta_eta = FALSE, save_eta = TRUE, eager = FALSE, profile = FALSE),
    outputs = list(predictions = TRUE, diagnostics = TRUE, plots = TRUE, zip = TRUE),
    runtime = list(python = app_env$hmschpc_python_bin(list()),
                   python_source = app_env$hmschpc_python_source_dir())
  )
}

jsdm_cfg <- function(args = NULL) {
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
    mcmc = list(burnin = arg_int(args, "jsdm_burnin", 200L, min_value = 0L),
                mcmc = arg_int(args, "jsdm_mcmc", 300L, min_value = 1L),
                thin = arg_int(args, "jsdm_thin", 5L, min_value = 1L),
                seed = 20260601L,
                verbose = 0L, preset = "Universal reviewer benchmark", ropt = 0.44),
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

gjam_cfg <- function(species, args = NULL) {
  list(
    project_name = "Universal_Benchmark_GJAM", engine = "GJAM",
    data = list(Y = "Y.csv", XData = "XData.csv", typeNames = "typeNames.csv",
                censor = "censor.csv", effort = "effort.csv", newdata = "newdata.csv",
                specByTrait = "specByTrait.csv", traitTypes = "traitTypes.csv",
                holdoutIndex = "holdoutIndex.csv"),
    model = list(formula = "~ pH + moisture + canopy + elevation + substrate",
                 type_single = "PA", typeNames_text = paste(rep("PA", length(species)), collapse = ","),
                 notStandard = "",
                 ng = arg_int(args, "gjam_ng", 700L, min_value = 10L),
                 burnin = arg_int(args, "gjam_burnin", 150L, min_value = 0L),
                 holdoutN = 0L,
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

spocc_cfg <- function(args = NULL) {
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
    mcmc = list(n.batch = arg_int(args, "spocc_n_batch", 30L, min_value = 1L),
                batch.length = arg_int(args, "spocc_batch_length", 10L, min_value = 1L),
                n.burn = arg_int(args, "spocc_n_burn", 100L, min_value = 0L),
                n.thin = arg_int(args, "spocc_n_thin", 5L, min_value = 1L),
                n.chains = arg_int(args, "spocc_n_chains", 1L, min_value = 1L),
                accept.rate = 0.43, n.report = 5L,
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

sjsdm_cfg <- function(args = NULL) {
  list(
    project_name = "Universal_Benchmark_sjSDM", engine = "sjSDM",
    model = list(family = "binomial_probit", env_model = "linear",
                 env_formula = "~ pH + moisture + canopy + elevation + substrate",
                 spatial_model = "none", spatial_formula = "~ 0 + .",
                 se = FALSE,
                 iter = arg_int(args, "sjsdm_iter", 100L, min_value = 1L),
                 step_size = arg_int(args, "sjsdm_step_size", 32L, min_value = 1L),
                 sampling = arg_int(args, "sjsdm_sampling", 1000L, min_value = 1L),
                 parallel = 0L, dtype = "float32", verbose = FALSE, seed = 20260601L),
    regularization_biotic = list(env_lambda = 0, env_alpha = 0.5,
                                 spatial_lambda = 0, spatial_alpha = 0.5,
                                 biotic_lambda = arg_num(args, "sjsdm_biotic_lambda", 0.001, min_value = 0),
                                 biotic_alpha = 0.5,
                                 biotic_df = NA, on_diag = FALSE,
                                 reg_on_Cov = TRUE, inverse = FALSE,
                                 tune_regularization = FALSE, tune_steps = 0L, cv_k = 0L),
    dnn_optimizer = list(hidden = "6,4", activation = "selu", dropout = 0,
                         bias = TRUE, optimizer = "Adamax",
                         learning_rate = arg_num(args, "sjsdm_learning_rate", 0.010, min_value = 1e-6),
                         weight_decay = arg_num(args, "sjsdm_weight_decay", 0, min_value = 0),
                         device = "cpu"),
    control = list(scheduler = 0L, lr_reduce_factor = 0.99,
                   early_stopping_training = 0L, mixed = FALSE),
    spatial_anova_metacommunity = list(generateSpatialEV = FALSE, spatial_ev_k = 4L,
                                       spatial_ev_threshold = 0, include_space_in_anova = FALSE,
                                       do_anova = TRUE, anova_samples = 100L,
                                       do_internal = FALSE, internal_fractions = "proportional",
                                       do_assembly = FALSE, assembly_predictor = ""),
    outputs = list(real_fit = TRUE, predict = TRUE, Rsquared = TRUE, importance = TRUE,
                   weights = TRUE, coef = TRUE, covariance_correlation = TRUE,
                   residuals = TRUE, plots = TRUE, zip = TRUE, copy_inputs = TRUE,
                   save_config = TRUE, report = TRUE)
  )
}

boral_cfg <- function(args = NULL) {
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
    mcmc_prior = list(n.burnin = arg_int(args, "boral_n_burnin", 250L, min_value = 0L),
                      n.iteration = arg_int(args, "boral_n_iteration", 750L, min_value = 1L),
                      n.thin = arg_int(args, "boral_n_thin", 10L, min_value = 1L),
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
  rscript_clean_env()
  system2(rscript, shQuote(script), stdout = stdout, stderr = stderr)
}

write_engine_inputs <- function(outdir, dat, engine) {
  write_csv_pair(dat$Y, outdir, "Y.csv")
  write_csv_pair(dat$X, outdir, "XData.csv")
  write_csv_pair(dat$traits, outdir, if (engine == "jSDM") "trait_data.csv" else "traits.csv")
  write_csv_pair(dat$study, outdir, "studyDesign.csv")
  write_csv_pair(dat$coords, outdir, if (engine == "spOccupancy") "coords.csv" else "coordinates.csv")
  write_csv_pair(dat$taxonomy, outdir, "taxonomy.csv")
  write_text_pair(dat$newick, outdir, if (engine == "Hmsc-HPC") "phylo_tree.nwk" else "phylogeny.nwk")
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

scale_numeric_frame <- function(train, new = NULL) {
  train <- as.data.frame(train, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(new)) new <- as.data.frame(new, check.names = FALSE, stringsAsFactors = FALSE)
  numeric_cols <- names(train)[vapply(train, is.numeric, logical(1))]
  meta <- data.frame(column = character(), center = numeric(), scale = numeric(),
                     transformed = logical(), note = character(), stringsAsFactors = FALSE)
  for (nm in numeric_cols) {
    center <- mean(train[[nm]], na.rm = TRUE)
    scale <- stats::sd(train[[nm]], na.rm = TRUE)
    if (!is.finite(center)) center <- 0
    if (!is.finite(scale) || scale <= 0) scale <- 1
    train[[nm]] <- (train[[nm]] - center) / scale
    if (!is.null(new) && nm %in% names(new)) {
      new[[nm]] <- (suppressWarnings(as.numeric(new[[nm]])) - center) / scale
    }
    meta <- rbind(meta, data.frame(column = nm, center = center, scale = scale,
                                   transformed = TRUE,
                                   note = "Numeric column standardized for optimizer stability and benchmark truth-scale matching.",
                                   stringsAsFactors = FALSE))
  }
  if (!length(numeric_cols)) {
    meta <- data.frame(column = NA_character_, center = NA_real_, scale = NA_real_,
                       transformed = FALSE, note = "No numeric columns were standardized.",
                       stringsAsFactors = FALSE)
  }
  list(train = train, new = new, metadata = meta)
}

standardize_engine_synonyms <- function(outdir, engine, dat) {
  dir.create(file.path(outdir, "standard"), recursive = TRUE, showWarnings = FALSE)
  st <- read_status(outdir, engine)
  status <- st$status %||% "unknown"
  species <- colnames(dat$Y)
  canonical_species <- function(values) {
    values <- as.character(values)
    key <- function(x) gsub("[^a-z0-9]", "", tolower(x))
    source_keys <- key(species)
    if (anyDuplicated(source_keys)) return(values)
    mapped <- species[match(key(values), source_keys)]
    ifelse(is.na(mapped), values, mapped)
  }
  empirical <- isTRUE(dat$empirical)
  predictors <- if (empirical) names(dat$X) else unique(dat$truth_effects$predictor)
  read_std <- function(name) {
    p <- file.path(outdir, "standard", name)
    if (file.exists(p)) tryCatch(read.csv(p, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  }
  if (identical(engine, "spOccupancy")) {
    split_species_parameter <- function(x) {
      x <- as.character(x)
      matched <- vapply(x, function(value) {
        hits <- species[endsWith(value, paste0("-", species))]
        if (length(hits) == 1L) hits else NA_character_
      }, character(1))
      term <- ifelse(is.na(matched), NA_character_,
                     substr(x, 1L, nchar(x) - nchar(matched) - 1L))
      list(species = matched, term = term)
    }
    canonical_predictor <- function(x) {
      key <- function(v) gsub("[^a-z0-9]", "", tolower(v))
      keys <- key(predictors)
      if (anyDuplicated(keys)) return(x)
      mapped <- predictors[match(key(x), keys)]
      ifelse(is.na(mapped), x, mapped)
    }
    occupancy <- read_std("effects_long.csv")
    if (!is.null(occupancy) && nrow(occupancy) &&
        all(c("response_id", "predictor") %in% names(occupancy))) {
      parts <- split_species_parameter(occupancy$predictor)
      keep <- occupancy$response_id == "occupancy" & !is.na(parts$species) &
        parts$term != "(Intercept)"
      occupancy <- occupancy[keep, , drop = FALSE]
      occupancy$response_id <- parts$species[keep]
      occupancy$predictor <- canonical_predictor(parts$term[keep])
      write.csv(occupancy, file.path(outdir, "standard", "effects_long.csv"), row.names = FALSE)
    }
    alpha_path <- file.path(outdir, "tables", "summary_alpha.csv")
    if (file.exists(alpha_path)) {
      alpha <- tryCatch(read.csv(alpha_path, check.names = FALSE, stringsAsFactors = FALSE),
                        error = function(e) NULL)
      if (!is.null(alpha) && all(c("parameter", "mean", "q025", "q975") %in% names(alpha))) {
        parts <- split_species_parameter(alpha$parameter)
        keep <- !is.na(parts$species)
        detection <- data.frame(engine = engine, species = parts$species[keep],
                                predictor = parts$term[keep],
                                estimate = alpha$mean[keep], lower = alpha$q025[keep],
                                upper = alpha$q975[keep], effect_type = "detection",
                                scale = "detection_logit_link", comparable = FALSE,
                                note = "Detection-process coefficient; not an occurrence environmental effect.",
                                stringsAsFactors = FALSE)
        write.csv(detection, file.path(outdir, "standard", "detection_effects.csv"), row.names = FALSE)
      }
    }
  }
  if (identical(engine, "boral")) {
    coef_path <- file.path(outdir, "tables", "X.coefs.mean.rds")
    if (file.exists(coef_path)) {
      coef_mat <- tryCatch(readRDS(coef_path), error = function(e) NULL)
      if (!is.null(coef_mat) && length(dim(coef_mat)) == 2L) {
        coef_mat <- as.matrix(coef_mat)
        if (is.null(colnames(coef_mat)) && ncol(coef_mat) == length(predictors)) {
          colnames(coef_mat) <- predictors
        }
        eff <- as.data.frame(as.table(coef_mat), stringsAsFactors = FALSE)
        names(eff) <- c("response_id", "predictor", "estimate")
        eff$engine <- "boral"
        eff$direction <- ifelse(is.na(eff$estimate), NA_character_,
                                ifelse(eff$estimate > 0, "positive",
                                       ifelse(eff$estimate < 0, "negative", "zero")))
        eff$lower <- NA_real_
        eff$upper <- NA_real_
        eff$notes <- "Mapped from tables/X.coefs.mean.rds after fitted boral model."
        eff <- eff[, c("engine", "response_id", "predictor", "direction", "estimate", "lower", "upper", "notes")]
        write.csv(eff, file.path(outdir, "standard", "effects_long.csv"), row.names = FALSE)
      }
    }
    cor_path <- file.path(outdir, "tables", "residual_cor.csv")
    if (file.exists(cor_path)) {
      cor_tab <- tryCatch(read.csv(cor_path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(cor_tab) && all(c("Var1", "Var2", "correlation") %in% names(cor_tab))) {
        upper <- !is.na(cor_tab$Var1) & !is.na(cor_tab$Var2) &
          as.character(cor_tab$Var1) < as.character(cor_tab$Var2)
        if ("component" %in% names(cor_tab)) upper <- upper & cor_tab$component == "cor"
        assoc <- cor_tab[upper, , drop = FALSE]
        if (nrow(assoc)) {
          assoc_long <- data.frame(
            engine = "boral",
            species_i = as.character(assoc$Var1),
            species_j = as.character(assoc$Var2),
            association_type = "boral_residual_correlation",
            estimate = suppressWarnings(as.numeric(assoc$correlation)),
            lower = NA_real_,
            upper = NA_real_,
            notes = "Mapped from tables/residual_cor.csv; not directly equivalent to Hmsc Omega or other association parameters.",
            stringsAsFactors = FALSE
          )
          write.csv(assoc_long, file.path(outdir, "standard", "associations_long.csv"), row.names = FALSE)
        }
      }
    }
  }
  effects_long <- read_std("effects_long.csv")
  effects <- NULL
  if (!is.null(effects_long) && nrow(effects_long) && "estimate" %in% names(effects_long)) {
    sp_col <- intersect(c("species", "response_id", "response", "column_id"), names(effects_long))[1]
    pr_col <- intersect(c("predictor", "parameter", "row_id"), names(effects_long))[1]
    if (!is.na(sp_col) && !is.na(pr_col)) {
      effects <- data.frame(
        engine = engine,
        species = canonical_species(effects_long[[sp_col]]),
        predictor = as.character(effects_long[[pr_col]]),
        estimate = suppressWarnings(as.numeric(effects_long$estimate)),
        lower = if ("lower" %in% names(effects_long)) suppressWarnings(as.numeric(effects_long$lower)) else NA_real_,
        upper = if ("upper" %in% names(effects_long)) suppressWarnings(as.numeric(effects_long$upper)) else NA_real_,
        statistic = NA_real_,
        p_or_support = NA_real_,
        effect_type = if (identical(engine, "spOccupancy")) "occurrence" else "environment",
        scale = switch(engine, GJAM = "observation_or_latent_GJAM", spOccupancy = "occupancy_link", "engine_link"),
        comparable = identical(status, "fitted") && !empirical,
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
                          species = canonical_species(pred_long[[sp_col]]),
                          observed = if (!is.na(obs_col)) suppressWarnings(as.numeric(pred_long[[obs_col]])) else NA_real_,
                          predicted = suppressWarnings(as.numeric(pred_long[[pred_col]])),
                          source_prediction_scale = if ("prediction_scale" %in% names(pred_long))
                            as.character(pred_long$prediction_scale) else NA_character_,
                          stringsAsFactors = FALSE)
    }
  }
  if (is.null(preds) || !nrow(preds)) {
    obs <- as.data.frame(as.table(as.matrix(dat$Y)), stringsAsFactors = FALSE)
    names(obs) <- c("site_id", "species", "observed")
    preds <- if (empirical) obs else merge(dat$truth_predictions[, c("site_id", "species")],
                                          obs, by = c("site_id", "species"), all.x = TRUE)
    preds$engine <- engine
    preds$predicted <- NA_real_
    preds <- preds[, c("engine", "site_id", "species", "observed", "predicted")]
  }
  if (!empirical) preds <- merge(preds, dat$truth_predictions, by = c("site_id", "species"), all.x = TRUE)
  preds$residual <- suppressWarnings(as.numeric(preds$observed) - as.numeric(preds$predicted))
  fallback_scale <- switch(engine,
                           GJAM = "PA_observation_probability",
                           spOccupancy = "occupancy_probability",
                           `Hmsc-HPC` = "unverified_engine_scale",
                           "occurrence_probability")
  preds$prediction_scale <- if ("source_prediction_scale" %in% names(preds)) {
    ifelse(!is.na(preds$source_prediction_scale) & nzchar(preds$source_prediction_scale),
           preds$source_prediction_scale, fallback_scale)
  } else fallback_scale
  preds <- preds[, c("engine", "site_id", "species", "observed", "predicted",
                     if (!empirical) "truth_probability", "residual", "prediction_scale")]
  write.csv(preds, file.path(outdir, "standard", "predictions_site_species.csv"), row.names = FALSE)

  assoc_long <- read_std("associations_long.csv")
  assoc <- NULL
  if (!is.null(assoc_long) && nrow(assoc_long) && "estimate" %in% names(assoc_long)) {
    s1 <- intersect(c("species_i", "response_1", "response1"), names(assoc_long))[1]
    s2 <- intersect(c("species_j", "response_2", "response2"), names(assoc_long))[1]
    if (!is.na(s1) && !is.na(s2)) {
      assoc_type <- if ("association_type" %in% names(assoc_long)) as.character(assoc_long$association_type) else "engine_association"
      assoc_scale <- paste0(engine, "_association_scale")
      assoc_note <- paste0("Associations are not numerically interchangeable across engines; status=", status)
      if (identical(engine, "Hmsc-HPC")) {
        assoc_type <- sub("^lambda_correlation", "derived_association_from_lambda", assoc_type)
        assoc_scale <- "Hmsc-HPC_derived_species_association_scale"
        assoc_note <- paste0(
          "Derived species-by-species association summary from Hmsc-HPC random-level loadings, not raw Lambda; ",
          "associations are not numerically interchangeable across engines; status=", status
        )
      }
      assoc <- data.frame(engine = engine, species_i = canonical_species(assoc_long[[s1]]),
                          species_j = canonical_species(assoc_long[[s2]]),
                          estimate = suppressWarnings(as.numeric(assoc_long$estimate)),
                          lower = if ("lower" %in% names(assoc_long)) suppressWarnings(as.numeric(assoc_long$lower)) else NA_real_,
                          upper = if ("upper" %in% names(assoc_long)) suppressWarnings(as.numeric(assoc_long$upper)) else NA_real_,
                          association_type = assoc_type,
                          scale = assoc_scale,
                          comparable = identical(status, "fitted") && !empirical,
                          note = assoc_note,
                          stringsAsFactors = FALSE)
    }
  }
  if (is.null(assoc) || !nrow(assoc)) {
    assoc <- if (empirical) data.frame(species_i = NA_character_, species_j = NA_character_)
      else dat$truth_associations[, c("species_i", "species_j")]
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
  if (!is.null(dat$truth_predictions) && "reference_method" %in% names(dat$truth_predictions)) {
    reference_methods <- sort(unique(as.character(dat$truth_predictions$reference_method)))
    reference_methods <- reference_methods[nzchar(reference_methods) & !is.na(reference_methods)]
    write.csv(data.frame(
      field = "truth_probability",
      reference_method = reference_methods,
      interpretation = ifelse(
        grepl("species_prevalence", reference_methods),
        "Observed training-set species prevalence; an audit reference, not simulated ecological truth.",
        ifelse(
          grepl("species_mean_response", reference_methods),
          "Observed training-set mean response on the fitted response scale; an audit reference, not an occurrence probability or simulated truth.",
          "Reference supplied by the benchmark data generator; consult the benchmark configuration for its statistical scale."
        )
      ),
      newdata_rule = "NA for prediction-only rows unless external or simulated truth is available.",
      stringsAsFactors = FALSE
    ), file.path(outdir, "standard", "reference_values_metadata.csv"), row.names = FALSE)
  }
  write.csv(data.frame(file = list.files(outdir, recursive = TRUE), stringsAsFactors = FALSE),
            file.path(outdir, "standard", "output_manifest.csv"), row.names = FALSE)
  invisible(TRUE)
}

run_engine <- function(engine, master_dir, dat, app_env, log_fun = message, args = NULL) {
  engine_dir <- file.path(master_dir, "03_engine_outputs", engine)
  ensure_dirs(engine_dir, engine)
  start <- Sys.time()
  log_fun(paste("Starting", engine))
  write_engine_inputs(engine_dir, dat, engine)
  status <- NULL
  err <- NULL
  captured_warnings <- character()
  empirical <- isTRUE(dat$empirical)
  adapt_cfg <- function(cfg) if (empirical) empirical_engine_config(cfg, engine, dat) else cfg
  source_note <- if (empirical) "Observed binary incidence; see the empirical benchmark input and design records."
    else "Synthetic benchmark input from the simulation contract."
  tryCatch(withCallingHandlers({
    if (empirical && identical(engine, "spOccupancy") && is.null(dat$det)) {
      write_status(engine_dir, engine, "not_applicable", errors = "No observed visits.csv was supplied. Detection histories were not generated.")
    } else if (identical(engine, "Hmsc")) {
      cfg <- adapt_cfg(hmsc_cfg(args))
      yaml::write_yaml(cfg, file.path(engine_dir, "used_config.yml"))
      write_csv_pair(dat$traits, engine_dir, "TrData.csv")
      write_csv_pair(dat$study, engine_dir, "studyDesign.csv")
      write_csv_pair(dat$coords, engine_dir, "coordinates.csv")
      write_csv_pair(dat$taxonomy, engine_dir, "taxonomy.csv")
      check <- app_env$validate_hmsc(dat$Y, dat$X, dat$traits, dat$study, dat$coords,
                                     distr = cfg$model$distr, XFormula = cfg$model$XFormula,
                                     newdata = if (empirical) dat$newdata, newstudy = if (empirical) dat$newstudy,
                                     TrFormula = cfg$model$TrFormula, use_traits = cfg$model$use_traits,
                                     use_phylogeny = cfg$model$use_phylogeny, random_mode = cfg$model$random_mode,
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
        res <- app_env$run_hmsc_s1s7_pipeline(engine_dir, cfg, dat$Y, dat$X, dat$traits, dat$study, dat$coords,
                   newdata = if (empirical) dat$newdata, newStudyDesign = if (empirical) dat$newstudy,
                   log_fun = log_fun)
        status <- read_status(engine_dir, engine)
        if (is.null(status$status) || identical(status$status, "missing_status")) {
          write_status(engine_dir, engine, res$status %||% "fit_failed",
                       warnings = res$warnings %||% character(), errors = res$errors %||% character())
        }
      }
    } else if (identical(engine, "Hmsc-HPC")) {
      cfg <- adapt_cfg(hmschpc_cfg(app_env, args))
      hpc_tree_file <- if (!empirical) file.path(engine_dir, "data", "phylo_tree.nwk")
      app_env$run_hmschpc_workflow(engine_dir, cfg, dat$Y, dat$X, Tr = dat$traits,
                                   study = dat$study, coord = dat$coords,
                                   phylo_tree_file = hpc_tree_file,
                                   taxonomy = dat$taxonomy, newdata = dat$newdata,
                                   log_fun = log_fun)
    } else if (identical(engine, "jSDM")) {
      cfg <- adapt_cfg(jsdm_cfg(args))
      yaml::write_yaml(cfg, file.path(engine_dir, "used_config.yml"))
      app_env$write_data_check_messages(engine_dir, source_note)
      app_env$write_jsdm_reproducible_script(engine_dir)
      code <- run_rscript(file.path(engine_dir, "reproducible_script", "run_this_jSDM_analysis.R"), engine_dir, engine)
      if (!identical(as.integer(code), 0L) && !(read_status(engine_dir, engine)$status %in% c("fitted", "fit_failed", "check_failed"))) {
        write_status(engine_dir, engine, "fit_failed", errors = paste0("jSDM Rscript exited with status ", code, "."))
      }
    } else if (identical(engine, "GJAM")) {
      cfg <- adapt_cfg(gjam_cfg(colnames(dat$Y), args))
      yaml::write_yaml(cfg, file.path(engine_dir, "used_config.yml"))
      type_tab <- data.frame(response = colnames(dat$Y), typeName = "PA", stringsAsFactors = FALSE)
      write_csv_pair(type_tab, engine_dir, "typeNames.csv", row.names = FALSE)
      write_csv_pair(dat$traits, engine_dir, "specByTrait.csv")
      if (!empirical) write_csv_pair(data.frame(trait = names(dat$traits), typeName = c("CAT", "CON", "CON"), stringsAsFactors = FALSE),
                                    engine_dir, "traitTypes.csv", row.names = FALSE)
      app_env$write_data_check_messages(engine_dir, source_note)
      app_env$write_gjam_reproducible_script(engine_dir)
      code <- run_rscript(file.path(engine_dir, "reproducible_script", "run_this_GJAM_analysis.R"), engine_dir, engine)
      if (!identical(as.integer(code), 0L) && !(read_status(engine_dir, engine)$status %in% c("fitted", "fit_failed", "check_failed"))) {
        write_status(engine_dir, engine, "fit_failed", errors = paste0("GJAM Rscript exited with status ", code, "."))
      }
    } else if (identical(engine, "spOccupancy")) {
      cfg <- adapt_cfg(spocc_cfg(args))
      yaml::write_yaml(cfg, file.path(engine_dir, "used_config.yml"))
      write_csv_pair(dat$det, engine_dir, "y.csv")
      write_csv_pair(dat$occ_covs, engine_dir, "occ.covs.csv")
      write_csv_pair(dat$det_cov, engine_dir, "det.covs.csv")
      write_csv_pair(dat$coords, engine_dir, "coords.csv")
      write_csv_pair(data.frame(species = colnames(dat$Y), stringsAsFactors = FALSE),
                     engine_dir, "species.csv", row.names = FALSE)
      app_env$write_data_check_messages(engine_dir, if (empirical) "Observed repeat-visit detection histories." else source_note)
      app_env$write_spoccupancy_reproducible_script(engine_dir)
      code <- run_rscript(file.path(engine_dir, "reproducible_script", "run_this_spOccupancy_analysis.R"), engine_dir, engine)
      if (!identical(as.integer(code), 0L) && !(read_status(engine_dir, engine)$status %in% c("fitted", "fit_failed", "check_failed"))) {
        write_status(engine_dir, engine, "fit_failed", errors = paste0("spOccupancy Rscript exited with status ", code, "."))
      }
    } else if (identical(engine, "sjSDM")) {
      cfg <- adapt_cfg(sjsdm_cfg(args))
      yaml::write_yaml(cfg, file.path(engine_dir, "used_config.yml"))
      env_scaled <- if (empirical) list(train = dat$X, new = dat$newdata, metadata = dat$predictor_map)
        else scale_numeric_frame(dat$X, dat$newdata)
      spatial_scaled <- if (empirical) list(train = NULL, new = NULL, metadata = data.frame(note = "No spatial term in the empirical baseline."))
        else scale_numeric_frame(dat$coords, dat$newcoords)
      write_csv_pair(env_scaled$train, engine_dir, "env.csv")
      write_csv_pair(env_scaled$new, engine_dir, "newdata.csv")
      write_csv_pair(spatial_scaled$train, engine_dir, "spatial.csv")
      write_csv_pair(spatial_scaled$new, engine_dir, "new_spatial.csv")
      write.csv(env_scaled$metadata, file.path(engine_dir, "tables", "sjSDM_env_scaling.csv"), row.names = FALSE)
      write.csv(spatial_scaled$metadata, file.path(engine_dir, "tables", "sjSDM_spatial_scaling.csv"), row.names = FALSE)
      app_env$write_data_check_messages(engine_dir, source_note)
      app_env$write_sjsdm_reproducible_script(engine_dir)
      code <- run_rscript(file.path(engine_dir, "reproducible_script", "run_this_sjSDM_analysis.R"), engine_dir, engine)
      if (!identical(as.integer(code), 0L) && !(read_status(engine_dir, engine)$status %in% c("fitted", "fit_failed", "check_failed"))) {
        write_status(engine_dir, engine, "fit_failed", errors = paste0("sjSDM Rscript exited with status ", code, "."))
      }
    } else if (identical(engine, "boral")) {
      cfg <- adapt_cfg(boral_cfg(args))
      yaml::write_yaml(cfg, file.path(engine_dir, "used_config.yml"))
      app_env$write_data_check_messages(engine_dir, source_note)
      app_env$write_boral_reproducible_script(engine_dir)
      code <- run_rscript(file.path(engine_dir, "reproducible_script", "run_this_boral_analysis.R"), engine_dir, engine)
      if (!identical(as.integer(code), 0L) && !(read_status(engine_dir, engine)$status %in% c("fitted", "fit_failed", "check_failed", "model_defined"))) {
        write_status(engine_dir, engine, "fit_failed", errors = paste0("boral Rscript exited with status ", code, "."))
      }
    } else {
      write_status(engine_dir, engine, "check_failed", errors = paste0("Unknown universal benchmark engine: ", engine))
    }
  }, warning = function(w) {
    captured_warnings <<- unique(c(captured_warnings, conditionMessage(w)))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    err <<- conditionMessage(e)
    write_status(engine_dir, engine, "fit_failed", errors = err,
                 runtime_seconds = round(as.numeric(difftime(Sys.time(), start, units = "secs")), 2))
  })
  if (length(captured_warnings) > 0) {
    dir.create(file.path(engine_dir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
    write.csv(data.frame(engine = engine, warning = captured_warnings, stringsAsFactors = FALSE),
              file.path(engine_dir, "diagnostics", "captured_package_warnings.csv"), row.names = FALSE)
    caveat_note <- paste0("Package warning captured during benchmark execution: ", captured_warnings)
    caveat_file <- file.path(engine_dir, "diagnostics", "scientific_caveats.csv")
    old_caveats <- if (file.exists(caveat_file)) {
      tryCatch(read.csv(caveat_file, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) data.frame())
    } else data.frame()
    add_caveats <- data.frame(engine = engine,
                              caveat = paste0("captured_package_warning_", seq_along(captured_warnings)),
                              control = "Review captured package warnings before ecological interpretation; tiny benchmark data can trigger benign numerical warnings.",
                              note = caveat_note,
                              stringsAsFactors = FALSE)
    write.csv(rbind(old_caveats, add_caveats), caveat_file, row.names = FALSE)
    status_file <- file.path(engine_dir, "diagnostics", "engine_status.json")
    if (file.exists(status_file) && safe_require("jsonlite")) {
      status_payload <- tryCatch(jsonlite::fromJSON(status_file, simplifyVector = FALSE), error = function(e) NULL)
      if (!is.null(status_payload)) {
        status_payload$scientific_caveats <- unique(c(status_payload$scientific_caveats %||% character(), caveat_note))
        status_payload$operational_risks <- status_payload$operational_risks %||% character()
        writeLines(jsonlite::toJSON(status_payload, pretty = TRUE, auto_unbox = TRUE), status_file)
      }
    }
  }
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
  generation_cfg <- if (safe_require("yaml")) {
    tryCatch(yaml::read_yaml(file.path(dat$dir, "truth", "data_generation_config.yml")),
             error = function(e) list())
  } else list()
  generating_structure <- generation_cfg$generating_structure %||% "unknown"
  alignment <- function(engine) {
    if (identical(generating_structure, "probit_latent") &&
        engine %in% c("Hmsc", "Hmsc-HPC", "jSDM", "sjSDM")) return("link_aligned_template")
    if (identical(generating_structure, "logit_occupancy_detection") &&
        identical(engine, "spOccupancy")) return("observation_process_aligned_template")
    return("cross_structure_or_unknown")
  }
  statuses$generating_structure <- generating_structure
  statuses$alignment <- vapply(statuses$engine, alignment, character(1))
  statuses$native_package_simulation <- FALSE
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
  detection_effects <- read_engine("detection_effects.csv")
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
  write.csv(statuses[, c("engine", "generating_structure", "alignment", "native_package_simulation"), drop = FALSE],
            file.path(out4, "generating_structure_alignment.csv"), row.names = FALSE)
  write.csv(effects, file.path(out4, "effects_all_models_long.csv"), row.names = FALSE)
  write.csv(detection_effects, file.path(out4, "detection_effects_all_models.csv"), row.names = FALSE)
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
  probability_scale <- truth_pred$prediction_scale %in%
    c("occurrence_probability", "occupancy_probability", "PA_observation_probability")
  pred_numeric <- suppressWarnings(as.numeric(truth_pred$predicted))
  truth_pred$error_to_truth[!probability_scale | !is.finite(pred_numeric) |
                              pred_numeric < 0 | pred_numeric > 1] <- NA_real_
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
      idx <- truth_pred$engine == e & is.finite(truth_pred$error_to_truth)
      pred <- suppressWarnings(as.numeric(truth_pred$predicted[idx]))
      truth <- suppressWarnings(as.numeric(truth_pred$truth_probability[idx]))
      ok <- is.finite(pred) & is.finite(truth)
      if (sum(ok) < 3L || length(unique(pred[ok])) < 2L) return(NA_real_)
      unname(coef(lm(truth[ok] ~ pred[ok]))[1])
    }, numeric(1))
    pred_metrics$calibration_slope <- vapply(pred_metrics$engine, function(e) {
      idx <- truth_pred$engine == e & is.finite(truth_pred$error_to_truth)
      pred <- suppressWarnings(as.numeric(truth_pred$predicted[idx]))
      truth <- suppressWarnings(as.numeric(truth_pred$truth_probability[idx]))
      ok <- is.finite(pred) & is.finite(truth)
      if (sum(ok) < 3L || length(unique(pred[ok])) < 2L) return(NA_real_)
      unname(coef(lm(truth[ok] ~ pred[ok]))[2])
    }, numeric(1))
    pred_metrics$pearson_correlation <- vapply(pred_metrics$engine, function(e) {
      idx <- truth_pred$engine == e & is.finite(truth_pred$error_to_truth)
      pred <- suppressWarnings(as.numeric(truth_pred$predicted[idx]))
      truth <- suppressWarnings(as.numeric(truth_pred$truth_probability[idx]))
      ok <- is.finite(pred) & is.finite(truth)
      if (sum(ok) < 3L || length(unique(pred[ok])) < 2L || length(unique(truth[ok])) < 2L) return(NA_real_)
      suppressWarnings(cor(pred[ok], truth[ok], method = "pearson"))
    }, numeric(1))
    pred_metrics$spearman_rank_correlation <- vapply(pred_metrics$engine, function(e) {
      idx <- truth_pred$engine == e & is.finite(truth_pred$error_to_truth)
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
             "Probability-scale prediction errors currently use the same in-sample simulated truth and matched IDs; they are not an independent hold-out validation or an engine ranking.",
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
                    paste0("<p>Generating structure: ", generating_structure,
                           ". These are Studio structural simulations, not native draws from each fitted package.</p>"),
                    "<h2>Status matrix</h2><pre>", html_table(statuses), "</pre>",
                    "<h2>Prediction metrics</h2><pre>", html_table(pred_metrics), "</pre>",
                    "<h2>Effect direction agreement</h2><pre>", html_table(direction_summary), "</pre>",
                    "<p>Raw association parameters are reported with notes and should not be treated as identical across engines.</p>",
                    "</body></html>")
  writeLines(compare_html, file.path(out6, "compare_models_report.html"), useBytes = TRUE)

  methods <- c(
    "# Universal Benchmark Methods",
    "",
    paste0("Generating structure: ", generating_structure, ". This is a JSDM Studio structural simulator, not the native simulate function of any included package."),
    "The probit_latent option uses a correlated Gaussian latent term and Bernoulli observations with a probit inverse link. The logit_occupancy_detection option uses the same explicit linear predictor with a logistic occupancy inverse link and a separate logistic replicated-detection process.",
    "The current truth_vs_predictions.csv table evaluates training-site probability recovery. Newdata prediction files are exported separately; they are not an independently simulated test set in this design.",
    "The benchmark distinguishes link/observation-process-aligned templates from cross-structure fits in 04_unified_standard_results/generating_structure_alignment.csv. Alignment is partial, not a claim that the full package generative model is recovered.",
    "JSDM Studio derives engine-specific input files with matched site, species and predictor identifiers. For a fairer methodological study, vary generating structures, replicate seeds, inspect convergence, and report within-structure and cross-structure results separately.",
    "The default recovery profile uses lower-collinearity environmental predictors, stronger non-zero species-by-predictor effects, and a residual multivariate normal latent process whose correlation matrix is written to truth/true_species_associations.csv.",
    "The phylogeny/taxonomy relatedness reference is generated separately from taxonomy ranks and is not reused as the residual association truth.",
    "Occurrence responses are used for Hmsc, Hmsc-HPC, jSDM, GJAM, sjSDM and boral. spOccupancy receives replicated detection-nondetection data generated from the same latent occupancy state.",
    "The benchmark is a software and reproducibility benchmark. Sampler settings are recorded in 09_configs/universal_benchmark_config.yml and should be increased, replicated and checked with convergence diagnostics before manuscript-scale statistical claims."
  )
  writeLines(methods, file.path(out7, "benchmark_methods_summary.md"), useBytes = TRUE)
  master_md <- c(
    "# JSDM Studio Universal Benchmark Report",
    "",
    paste0("Generated: ", Sys.time()),
    paste0("Generating structure: ", generating_structure, " (Studio structural simulator; not native package simulation)"),
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
  quote_r_vec <- function(x) {
    quote_one <- function(s) paste0("\"", gsub("\"", "\\\\\"", as.character(s), fixed = TRUE), "\"")
    paste0("c(", paste(vapply(x, quote_one, character(1)), collapse = ", "), ")")
  }
  runner_args <- cfg$runner_args %||% list()
  common_flags <- c(
    paste0("--n_sites=", cfg$n_sites),
    paste0("--n_species=", cfg$n_species),
    paste0("--n_visits=", cfg$n_visits),
    paste0("--seed=", cfg$seed),
    paste0("--generating_structure=", cfg$generating_structure),
    paste0("--engines=", paste(cfg$engines, collapse = ",")),
    if (nzchar(runner_args$truth_profile %||% "")) paste0("--truth_profile=", runner_args$truth_profile) else character(),
    if (nzchar(runner_args$effect_scale %||% "")) paste0("--effect_scale=", runner_args$effect_scale) else character(),
    if (nzchar(runner_args$association_strength %||% "")) paste0("--association_strength=", runner_args$association_strength) else character(),
    if (nzchar(runner_args$residual_sd %||% "")) paste0("--residual_sd=", runner_args$residual_sd) else character(),
    if (nzchar(runner_args$observation_noise_sd %||% "")) paste0("--observation_noise_sd=", runner_args$observation_noise_sd) else character(),
    if (nzchar(runner_args$hmsc_samples %||% "")) paste0("--hmsc_samples=", runner_args$hmsc_samples) else character(),
    if (nzchar(runner_args$hmsc_transient %||% "")) paste0("--hmsc_transient=", runner_args$hmsc_transient) else character(),
    if (nzchar(runner_args$hmsc_thin %||% "")) paste0("--hmsc_thin=", runner_args$hmsc_thin) else character(),
    if (nzchar(runner_args$hmsc_nChains %||% "")) paste0("--hmsc_nChains=", runner_args$hmsc_nChains) else character(),
    if (nzchar(runner_args$hmschpc_samples %||% "")) paste0("--hmschpc_samples=", runner_args$hmschpc_samples) else character(),
    if (nzchar(runner_args$hmschpc_transient %||% "")) paste0("--hmschpc_transient=", runner_args$hmschpc_transient) else character(),
    if (nzchar(runner_args$hmschpc_thin %||% "")) paste0("--hmschpc_thin=", runner_args$hmschpc_thin) else character(),
    if (nzchar(runner_args$hmschpc_chains %||% "")) paste0("--hmschpc_chains=", runner_args$hmschpc_chains) else character(),
    if (nzchar(runner_args$jsdm_burnin %||% "")) paste0("--jsdm_burnin=", runner_args$jsdm_burnin) else character(),
    if (nzchar(runner_args$jsdm_mcmc %||% "")) paste0("--jsdm_mcmc=", runner_args$jsdm_mcmc) else character(),
    if (nzchar(runner_args$jsdm_thin %||% "")) paste0("--jsdm_thin=", runner_args$jsdm_thin) else character(),
    if (nzchar(runner_args$gjam_ng %||% "")) paste0("--gjam_ng=", runner_args$gjam_ng) else character(),
    if (nzchar(runner_args$gjam_burnin %||% "")) paste0("--gjam_burnin=", runner_args$gjam_burnin) else character(),
    if (nzchar(runner_args$spocc_n_batch %||% "")) paste0("--spocc_n_batch=", runner_args$spocc_n_batch) else character(),
    if (nzchar(runner_args$spocc_batch_length %||% "")) paste0("--spocc_batch_length=", runner_args$spocc_batch_length) else character(),
    if (nzchar(runner_args$spocc_n_burn %||% "")) paste0("--spocc_n_burn=", runner_args$spocc_n_burn) else character(),
    if (nzchar(runner_args$spocc_n_thin %||% "")) paste0("--spocc_n_thin=", runner_args$spocc_n_thin) else character(),
    if (nzchar(runner_args$spocc_n_chains %||% "")) paste0("--spocc_n_chains=", runner_args$spocc_n_chains) else character(),
    if (nzchar(runner_args$sjsdm_iter %||% "")) paste0("--sjsdm_iter=", runner_args$sjsdm_iter) else character(),
    if (nzchar(runner_args$sjsdm_sampling %||% "")) paste0("--sjsdm_sampling=", runner_args$sjsdm_sampling) else character(),
    if (nzchar(runner_args$sjsdm_step_size %||% "")) paste0("--sjsdm_step_size=", runner_args$sjsdm_step_size) else character(),
    if (nzchar(runner_args$sjsdm_learning_rate %||% "")) paste0("--sjsdm_learning_rate=", runner_args$sjsdm_learning_rate) else character(),
    if (nzchar(runner_args$sjsdm_weight_decay %||% "")) paste0("--sjsdm_weight_decay=", runner_args$sjsdm_weight_decay) else character(),
    if (nzchar(runner_args$sjsdm_biotic_lambda %||% "")) paste0("--sjsdm_biotic_lambda=", runner_args$sjsdm_biotic_lambda) else character(),
    if (nzchar(runner_args$boral_n_burnin %||% "")) paste0("--boral_n_burnin=", runner_args$boral_n_burnin) else character(),
    if (nzchar(runner_args$boral_n_iteration %||% "")) paste0("--boral_n_iteration=", runner_args$boral_n_iteration) else character(),
    if (nzchar(runner_args$boral_n_thin %||% "")) paste0("--boral_n_thin=", runner_args$boral_n_thin) else character(),
    if (isTRUE(runner_args$allow_large %||% FALSE)) "--allow_large=true" else character()
  )
  writeLines(c("# Generate benchmark data",
               "if (.Platform$OS.type == 'windows') Sys.unsetenv(c('LC_ALL','LC_COLLATE','LC_CTYPE','LC_MONETARY','LC_TIME','LANG','LANGUAGE'))",
               sprintf("system2(file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript'), shQuote(%s))",
                       quote_r_vec(c(runner_rel, "--action=generate", common_flags)))),
             file.path(sdir, "generate_benchmark_data.R"))
  writeLines(c("# Run all benchmark engines",
               "if (.Platform$OS.type == 'windows') Sys.unsetenv(c('LC_ALL','LC_COLLATE','LC_CTYPE','LC_MONETARY','LC_TIME','LANG','LANGUAGE'))",
               sprintf("system2(file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript'), shQuote(%s))",
                       quote_r_vec(c(runner_rel, "--action=all", common_flags)))),
             file.path(sdir, "run_all_engines.R"))
  writeLines(c("# Rebuild comparison tables for an existing master folder",
               "if (.Platform$OS.type == 'windows') Sys.unsetenv(c('LC_ALL','LC_COLLATE','LC_CTYPE','LC_MONETARY','LC_TIME','LANG','LANGUAGE'))",
               sprintf("system2(file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript'), shQuote(c(%s, '--action=compare', paste0('--master_dir=', getwd()))))", deparse(runner_rel))),
             file.path(sdir, "compare_benchmark_results.R"))
  writeLines(c("# Rebuild the all-model ZIP for an existing master folder",
               "if (.Platform$OS.type == 'windows') Sys.unsetenv(c('LC_ALL','LC_COLLATE','LC_CTYPE','LC_MONETARY','LC_TIME','LANG','LANGUAGE'))",
               sprintf("system2(file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript'), shQuote(c(%s, '--action=zip', paste0('--master_dir=', getwd()))))", deparse(runner_rel))),
             file.path(sdir, "rebuild_all_model_zip.R"))
  if (safe_require("yaml")) {
    yaml::write_yaml(cfg, file.path(cdir, "universal_benchmark_config.yml"))
    yaml::write_yaml(hmsc_cfg(runner_args), file.path(cdir, "Hmsc_config.yml"))
    yaml::write_yaml(jsdm_cfg(runner_args), file.path(cdir, "jSDM_config.yml"))
    yaml::write_yaml(gjam_cfg(sprintf("sp_%02d", seq_len(cfg$n_species)), runner_args), file.path(cdir, "GJAM_config.yml"))
    yaml::write_yaml(spocc_cfg(runner_args), file.path(cdir, "spOccupancy_config.yml"))
    yaml::write_yaml(sjsdm_cfg(runner_args), file.path(cdir, "sjSDM_config.yml"))
    yaml::write_yaml(boral_cfg(runner_args), file.path(cdir, "boral_config.yml"))
  }
}

run_universal_benchmark <- function(args) {
  run_root <- file.path(app_dir, "output", paste0("universal_benchmark_", timestamp_id(), "_", args$generating_structure))
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
  truth_defaults <- truth_profile_defaults(args$truth_profile)
  truth_params <- list(
    effect_scale = arg_num(args, "effect_scale", truth_defaults$effect_scale, min_value = 0.1, max_value = 5),
    association_strength = arg_num(args, "association_strength", truth_defaults$association_strength, min_value = 0, max_value = 5),
    residual_sd = arg_num(args, "residual_sd", truth_defaults$residual_sd, min_value = 0.01, max_value = 5),
    observation_noise_sd = arg_num(args, "observation_noise_sd", truth_defaults$observation_noise_sd, min_value = 0.001, max_value = 5)
  )
  if (nzchar(input_dir)) {
    input_dir <- normalizePath(input_dir, winslash = "/", mustWork = FALSE)
    dat <- load_benchmark_data(input_dir, allow_generate = FALSE)
    data_source <- "user_or_example_input_dir"
    log_fun(paste0("Using supplied benchmark input directory: ", dat$dir))
  } else {
    target_dir <- trimws(as.character(args$target_dir %||% ""))
    if (!nzchar(target_dir)) target_dir <- file.path(app_dir, "examples", "universal_benchmark")
    dat <- generate_benchmark_data(target_dir,
                                   args$n_sites, args$n_species, args$n_visits, args$seed,
                                   allow_large = args$allow_large,
                                   truth_profile = args$truth_profile,
                                   generating_structure = args$generating_structure,
                                   effect_scale = truth_params$effect_scale,
                                   association_strength = truth_params$association_strength,
                                   residual_sd = truth_params$residual_sd,
                                   observation_noise_sd = truth_params$observation_noise_sd)
    dat <- load_benchmark_data(dat$dir, allow_generate = FALSE)
    data_source <- "generated_from_latent_truth"
    log_fun(paste0("Generated benchmark input directory: ", dat$dir))
  }
  copy_dir_contents(dat$dir, file.path(master_dir, "01_benchmark_inputs"))
  copy_dir_contents(file.path(dat$dir, "truth"), file.path(master_dir, "02_truth"))
  cfg <- list(seed = args$seed, n_sites = nrow(dat$Y), n_species = ncol(dat$Y),
              n_visits = args$n_visits, engines = args$engines, data_source = data_source,
              generating_structure = if (safe_require("yaml")) {
                yaml::read_yaml(file.path(dat$dir, "truth", "data_generation_config.yml"))$generating_structure %||% "unknown"
              } else "unknown",
              allow_large = isTRUE(args$allow_large),
              truth_profile = args$truth_profile,
              truth_parameters = truth_params,
              input_dir = dat$dir,
              runner_args = args,
              sampler_overrides = list(
                Hmsc = list(samples = arg_int(args, "hmsc_samples", 500L, min_value = 1L),
                            transient = arg_int(args, "hmsc_transient", 500L, min_value = 0L),
                            thin = arg_int(args, "hmsc_thin", 1L, min_value = 1L),
                            nChains = arg_int(args, "hmsc_nChains", 2L, min_value = 1L)),
                Hmsc_HPC = list(samples = arg_int(args, "hmschpc_samples", 500L, min_value = 1L),
                                transient = arg_int(args, "hmschpc_transient", 500L, min_value = 0L),
                                thin = arg_int(args, "hmschpc_thin", 1L, min_value = 1L),
                                chains = arg_int(args, "hmschpc_chains", 1L, min_value = 1L)),
                jSDM = list(burnin = arg_int(args, "jsdm_burnin", 200L, min_value = 0L),
                            mcmc = arg_int(args, "jsdm_mcmc", 300L, min_value = 1L),
                            thin = arg_int(args, "jsdm_thin", 5L, min_value = 1L)),
                GJAM = list(ng = arg_int(args, "gjam_ng", 700L, min_value = 10L),
                            burnin = arg_int(args, "gjam_burnin", 150L, min_value = 0L)),
                spOccupancy = list(n.batch = arg_int(args, "spocc_n_batch", 30L, min_value = 1L),
                                   batch.length = arg_int(args, "spocc_batch_length", 10L, min_value = 1L),
                                   n.burn = arg_int(args, "spocc_n_burn", 100L, min_value = 0L),
                                   n.thin = arg_int(args, "spocc_n_thin", 5L, min_value = 1L),
                                   n.chains = arg_int(args, "spocc_n_chains", 1L, min_value = 1L)),
                sjSDM = list(iter = arg_int(args, "sjsdm_iter", 100L, min_value = 1L),
                             sampling = arg_int(args, "sjsdm_sampling", 1000L, min_value = 1L),
                             step_size = arg_int(args, "sjsdm_step_size", 32L, min_value = 1L),
                             learning_rate = arg_num(args, "sjsdm_learning_rate", 0.010, min_value = 1e-6),
                             weight_decay = arg_num(args, "sjsdm_weight_decay", 0, min_value = 0),
                             biotic_lambda = arg_num(args, "sjsdm_biotic_lambda", 0.001, min_value = 0)),
                boral = list(n.burnin = arg_int(args, "boral_n_burnin", 250L, min_value = 0L),
                             n.iteration = arg_int(args, "boral_n_iteration", 750L, min_value = 1L),
                             n.thin = arg_int(args, "boral_n_thin", 10L, min_value = 1L))
              ))
  write_master_scripts(master_dir, cfg)
  preflight <- dependency_preflight(master_dir)
  app_env <- source_app_env()
  if (safe_require("yaml")) yaml::write_yaml(hmschpc_cfg(app_env, args), file.path(master_dir, "09_configs", "Hmsc-HPC_config.yml"))
  statuses <- do.call(rbind, lapply(args$engines, function(engine) run_engine(engine, master_dir, dat, app_env, log_fun, args = args)))
  collect_master_results(master_dir, args$engines, dat, statuses)
  master_zip <- make_zip_file(master_dir, file.path(run_root, "JSDMStudio_universal_benchmark_ALL_MODELS.zip"))
  sha <- if (file.exists(master_zip)) digest::digest(master_zip, algo = "sha256", file = TRUE) else NA_character_
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

write_design_matrix_recovery <- function(matrix_dir) {
  p <- file.path(matrix_dir, "design_matrix_status.csv")
  if (!file.exists(p)) stop("Missing design_matrix_status.csv", call. = FALSE)
  status <- read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
  read_one <- function(path, engine, column) {
    if (!file.exists(path)) return(NA_real_)
    x <- tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
                  error = function(e) NULL)
    if (is.null(x) || !all(c("engine", column) %in% names(x))) return(NA_real_)
    value <- suppressWarnings(as.numeric(x[[column]][x$engine == engine]))
    if (!length(value)) NA_real_ else value[[1]]
  }
  status$effect_direction_agreement <- NA_real_
  status$in_sample_probability_rmse <- NA_real_
  status$n_training_site_predictions <- NA_real_
  for (i in seq_len(nrow(status))) {
    master <- file.path(dirname(status$master_zip[[i]]), "master")
    root <- file.path(master, "05_truth_comparison")
    status$effect_direction_agreement[[i]] <- read_one(
      file.path(root, "effect_direction_agreement.csv"), status$engine[[i]],
      "effect_direction_agreement")
    status$in_sample_probability_rmse[[i]] <- read_one(
      file.path(root, "prediction_metric_by_model.csv"), status$engine[[i]],
      "rmse_to_truth_probability")
    status$n_training_site_predictions[[i]] <- read_one(
      file.path(root, "prediction_metric_by_model.csv"), status$engine[[i]],
      "n_predictions")
  }
  status$evaluation_scope <- "one-seed in-sample technical check; not a model ranking"
  write.csv(status[, c("generating_structure", "engine", "alignment", "status",
                       "effect_direction_agreement", "in_sample_probability_rmse",
                       "n_training_site_predictions", "runtime_seconds", "evaluation_scope")],
            file.path(matrix_dir, "design_matrix_recovery.csv"), row.names = FALSE)
  invisible(status)
}

run_design_matrix <- function(args) {
  if (nzchar(trimws(args$input_dir %||% ""))) {
    stop("A generating-structure matrix requires generated data; omit input_dir.", call. = FALSE)
  }
  matrix_dir <- file.path(app_dir, "output", paste0("universal_design_matrix_", timestamp_id()))
  dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
  structures <- c("probit_latent", "logit_occupancy_detection")
  rows <- lapply(structures, function(structure) {
    run_args <- args
    run_args$generating_structure <- structure
    run_args$target_dir <- file.path(matrix_dir, "generated_inputs", structure)
    result <- run_universal_benchmark(run_args)
    status_path <- file.path(result$master_dir, "04_unified_standard_results", "model_status_matrix.csv")
    status <- read.csv(status_path, check.names = FALSE, stringsAsFactors = FALSE)
    status$master_zip <- result$master_zip
    status$sha256 <- result$sha256
    status
  })
  matrix_status <- do.call(rbind, rows)
  write.csv(matrix_status, file.path(matrix_dir, "design_matrix_status.csv"), row.names = FALSE)
  write_design_matrix_recovery(matrix_dir)
  writeLines(c(
    "# Generating-structure design matrix",
    "",
    "The same random seed, site/species counts, signal settings and engine settings were applied to each generating structure.",
    "Rows marked link_aligned_template or observation_process_aligned_template share only stated parts of the generating mechanism with the fitted engine. They are not draws from that package's native simulator.",
    "Cross-structure results are sensitivity checks. Do not rank engines from a single seed, short chains, or non-equivalent association estimands.",
    "For a publication study, repeat seeds, check convergence and add true native-package simulation where supported."
  ), file.path(matrix_dir, "DESIGN_AND_LIMITATIONS.md"), useBytes = TRUE)
  cat("DESIGN_MATRIX_DIR=", normalizePath(matrix_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
  invisible(matrix_status)
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
                                 args$n_sites, args$n_species, args$n_visits, args$seed,
                                 allow_large = args$allow_large,
                                 truth_profile = args$truth_profile,
                                 generating_structure = args$generating_structure,
                                 effect_scale = arg_num(args, "effect_scale", NA_real_, min_value = 0.1, max_value = 5),
                                 association_strength = arg_num(args, "association_strength", NA_real_, min_value = 0, max_value = 5),
                                 residual_sd = arg_num(args, "residual_sd", NA_real_, min_value = 0.01, max_value = 5),
                                 observation_noise_sd = arg_num(args, "observation_noise_sd", NA_real_, min_value = 0.001, max_value = 5))
    cat("BENCHMARK_DATA_DIR=", x$dir, "\n", sep = "")
  } else if (identical(args$action, "empirical")) {
    run_empirical_benchmark(args)
  } else if (identical(args$action, "preflight")) {
    print(dependency_preflight())
  } else if (identical(args$action, "matrix")) {
    run_design_matrix(args)
  } else if (identical(args$action, "matrix_compare")) {
    write_design_matrix_recovery(args$master_dir)
    cat("DESIGN_MATRIX_COMPARE_OK\n")
  } else if (identical(args$action, "compare")) {
    compare_existing_master(args$master_dir)
    cat("COMPARE_OK\n")
  } else if (identical(args$action, "zip")) {
    zip <- make_zip_file(args$master_dir, file.path(dirname(args$master_dir), "JSDMStudio_universal_benchmark_ALL_MODELS.zip"))
    cat("MASTER_ZIP=", zip, "\n", sep = "")
    cat("SHA256=", digest::digest(zip, algo = "sha256", file = TRUE), "\n", sep = "")
  } else if (identical(args$action, "all")) {
    run_universal_benchmark(args)
  } else {
    stop("Unknown benchmark action: ", args$action, call. = FALSE)
  }
}

if (sys.nframe() == 0L) main()
