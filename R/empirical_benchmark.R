empirical_input_names <- c(
  y = "Y.csv", x = "XData.csv", visits = "visits.csv",
  new_x = "newdata.csv", test_y = "test_y.csv"
)

empirical_read_csv <- function(path, key = "site_id", optional = FALSE, unique_key = TRUE) {
  if (!file.exists(path) && optional) return(NULL)
  label <- basename(path)
  if (!file.exists(path)) stop(label, " is required.", call. = FALSE)
  dat <- read.csv(path, row.names = NULL, check.names = FALSE, colClasses = "character",
                  fileEncoding = "UTF-8-BOM", na.strings = c("", "NA"))
  if (!nrow(dat) || any(!nzchar(names(dat))) || anyDuplicated(names(dat)) ||
      ncol(dat) < 2L || !identical(names(dat)[1], key)) {
    stop(label, " needs a first column named ", key, " and distinct data columns.", call. = FALSE)
  }
  ids <- dat[[key]]
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || any(ids != trimws(ids)) ||
      (unique_key && anyDuplicated(ids))) {
    stop(label, " has missing, duplicate or whitespace-padded IDs.", call. = FALSE)
  }
  if (unique_key) rownames(dat) <- ids
  dat
}

empirical_align <- function(dat, ids, columns, label) {
  if (!setequal(rownames(dat), ids) || !setequal(names(dat), columns)) {
    stop(label, " site IDs or column names do not match its companion table.", call. = FALSE)
  }
  dat[ids, columns, drop = FALSE]
}

empirical_binary <- function(dat, label) {
  mat <- as.matrix(dat)
  if (anyNA(mat) || any(!mat %in% c("0", "1"))) {
    stop(label, " must contain observed 0/1 records only (no missing values).", call. = FALSE)
  }
  storage.mode(mat) <- "integer"
  mat
}

empirical_design <- function(x, new_x) {
  a <- list()
  b <- list()
  map <- list()
  for (nm in names(x)) {
    v <- x[[nm]]
    vn <- if (is.null(new_x)) NULL else new_x[[nm]]
    if (anyNA(v) || anyNA(vn) || any(!nzchar(trimws(c(v, vn))))) {
      stop("Missing covariate values in '", nm, "'.", call. = FALSE)
    }
    numeric_v <- suppressWarnings(as.numeric(v))
    if (all(!is.na(numeric_v))) {
      numeric_new <- suppressWarnings(as.numeric(vn))
      center <- mean(numeric_v)
      scale <- stats::sd(numeric_v)
      if (!is.finite(center) || !is.finite(scale) || scale <= 0 || any(!is.finite(numeric_new))) {
        stop("Covariate '", nm, "' needs finite numeric values and nonzero training variance.", call. = FALSE)
      }
      a <- c(a, list((numeric_v - center) / scale))
      b <- c(b, list((numeric_new - center) / scale))
      map[[length(map) + 1L]] <- data.frame(predictor = nm, input_column = nm, level = "",
        reference = "", center = center, scale = scale)
    } else {
      levels <- sort(unique(v), method = "radix")
      if (length(levels) < 2L || any(!vn %in% levels)) {
        stop("Factor '", nm, "' needs at least two training levels; newdata cannot introduce levels.", call. = FALSE)
      }
      for (level in levels[-1]) {
        a <- c(a, list(as.integer(v == level)))
        b <- c(b, list(as.integer(vn == level)))
        map[[length(map) + 1L]] <- data.frame(predictor = paste(nm, level, sep = " = "),
          input_column = nm, level = level, reference = levels[1], center = 0, scale = 1)
      }
    }
  }
  map <- do.call(rbind, map)
  map$engine_column <- sprintf("x%03d", seq_len(nrow(map)))
  X <- as.data.frame(setNames(a, map$engine_column))
  if (nrow(X) <= ncol(X) + 1L || qr(cbind(1, X))$rank < ncol(X) + 1L) {
    stop("Predictors are collinear or too numerous for the training sites. Remove redundant columns.", call. = FALSE)
  }
  rownames(X) <- rownames(x)
  Xnew <- if (!is.null(new_x)) as.data.frame(setNames(b, map$engine_column), row.names = rownames(new_x))
  list(train = X, newdata = Xnew, map = map)
}

load_empirical_benchmark <- function(dir) {
  files <- file.path(dir, empirical_input_names)
  names(files) <- names(empirical_input_names)
  y <- empirical_read_csv(files[["y"]])
  x <- empirical_read_csv(files[["x"]])
  ids <- rownames(y)
  species <- names(y)[-1]
  y <- empirical_binary(y[-1], "Y.csv")
  if (nrow(y) < 10L || ncol(y) < 2L || any(colSums(y) == 0 | colSums(y) == nrow(y))) {
    stop("Use at least 10 training sites and 2 species, with both outcomes observed for every species.", call. = FALSE)
  }
  x <- empirical_align(x[-1], ids, names(x)[-1], "XData.csv")
  new_x <- empirical_read_csv(files[["new_x"]], optional = TRUE)
  test_y <- empirical_read_csv(files[["test_y"]], optional = TRUE)
  new_ids <- rownames(new_x)
  if (!is.null(new_x)) {
    if (any(new_ids %in% ids)) stop("newdata.csv must contain site IDs absent from the training data.", call. = FALSE)
    new_x <- empirical_align(new_x[-1], new_ids, names(x), "newdata.csv")
  }
  if (!is.null(test_y)) {
    if (is.null(new_x)) stop("test_y.csv requires newdata.csv.", call. = FALSE)
    test_y <- empirical_binary(empirical_align(test_y[-1], new_ids, species, "test_y.csv"), "test_y.csv")
  }
  design <- empirical_design(x, new_x)
  site_key <- sprintf("site%05d", seq_along(ids))
  new_key <- if (length(new_ids)) sprintf("project%05d", seq_along(new_ids)) else character()
  species_key <- sprintf("sp%05d", seq_along(species))
  det <- NULL
  visit_count <- NULL
  visits <- empirical_read_csv(files[["visits"]], optional = TRUE, unique_key = FALSE)
  if (!is.null(visits)) {
    if (ncol(visits) < 4L || !identical(names(visits)[2], "visit_id") ||
        !setequal(names(visits)[-(1:2)], species) || !setequal(visits$site_id, ids) ||
        anyNA(visits$visit_id) || any(!nzchar(trimws(visits$visit_id))) ||
        anyDuplicated(visits[c("site_id", "visit_id")])) {
      stop("visits.csv requires site_id, visit_id, and the same species as Y.csv; use unique site/visit pairs.", call. = FALSE)
    }
    visit_y <- empirical_binary(visits[species], "visits.csv")
    groups <- split(seq_len(nrow(visits)), factor(visits$site_id, levels = ids))
    visit_count <- lengths(groups)
    if (min(visit_count) < 2L) stop("visits.csv must contain at least two actual surveys per training site.", call. = FALSE)
    pooled <- t(vapply(groups, function(i) as.integer(colSums(visit_y[i, , drop = FALSE]) > 0L),
      integer(length(species))))
    if (any(pooled != y)) stop("Y.csv must equal whether each species was detected in any recorded visit.", call. = FALSE)
    reps <- seq_len(max(visit_count))
    det_names <- as.vector(outer(species_key, paste0("_rep", reps), paste0))
    det <- matrix(NA_integer_, length(ids), length(det_names), dimnames = list(site_key, det_names))
    for (i in seq_along(groups)) {
      rows <- groups[[i]]
      det[i, seq_len(length(rows) * length(species))] <- as.vector(t(visit_y[rows, , drop = FALSE]))
    }
  }
  dimnames(y) <- list(site_key, species_key)
  rownames(design$train) <- site_key
  if (!is.null(design$newdata)) rownames(design$newdata) <- new_key
  if (!is.null(test_y)) dimnames(test_y) <- list(new_key, species_key)
  list(
    empirical = TRUE, Y = as.data.frame(y), X = design$train, newdata = design$newdata, test_y = test_y,
    study = data.frame(sample = site_key, row.names = site_key),
    newstudy = if (length(new_key)) data.frame(sample = new_key, row.names = new_key),
    trial = matrix(1L, nrow(y), ncol(y), dimnames = dimnames(y)),
    det = det, occ_covs = design$train, visit_count = visit_count, paths = files[file.exists(files)],
    species_map = data.frame(engine_id = species_key, species = species),
    site_map = rbind(data.frame(engine_id = site_key, site_id = ids, scope = "training"),
                     data.frame(engine_id = new_key, site_id = as.character(new_ids), scope = rep("newdata", length(new_key)))),
    predictor_map = design$map
  )
}

empirical_profile <- function(profile = "standard") {
  profiles <- c("quick", "standard")
  if (!profile %in% profiles) stop("Unknown empirical run preset.", call. = FALSE)
  if (profile == "standard") return(list(hmschpc_chains = 2L, spocc_n_chains = 2L,
    jsdm_mcmc = 1000L, jsdm_burnin = 500L, spocc_n_batch = 100L, spocc_n_burn = 300L))
  list(hmsc_samples = 20L, hmsc_transient = 30L, hmsc_nChains = 2L,
       hmschpc_samples = 20L, hmschpc_transient = 30L, jsdm_burnin = 50L,
       jsdm_mcmc = 100L, jsdm_thin = 2L, gjam_ng = 150L, gjam_burnin = 50L,
       spocc_n_batch = 20L, spocc_n_burn = 50L, spocc_n_thin = 2L,
       sjsdm_iter = 40L, sjsdm_sampling = 200L,
       boral_n_burnin = 100L, boral_n_iteration = 300L, boral_n_thin = 5L)
}

empirical_engine_config <- function(cfg, engine, dat) {
  rhs <- paste("~", paste(names(dat$X), collapse = " + "))
  switch(engine,
    Hmsc = {
      cfg$model$XFormula <- rhs
      cfg$model$XScale <- FALSE
      cfg$model$use_traits <- cfg$model$use_phylogeny <- FALSE
      cfg$data$phylogeny <- NULL
      cfg$model$random_effect_column <- "sample"
      cfg$outputs$predictions$partition_column <- "sample"
    },
    "Hmsc-HPC" = {
      cfg$model$XFormula <- rhs
      cfg$model$use_traits <- FALSE
      cfg$model$phylogeny_mode <- "none"
      cfg$inputs$phylo_tree <- cfg$inputs$traits <- NULL
      cfg$random_effects$column <- "sample"
      cfg$random_effects$name <- "sample"
    },
    jSDM = {
      cfg$model$site_formula <- rhs
      cfg$model$allow_traits <- cfg$model$scale_site_data <- FALSE
      cfg$model$site_effect <- "none"
      cfg$model$n_latent <- min(2L, ncol(dat$Y) - 1L)
      cfg$prediction$max_prediction_sites <- max(nrow(dat$Y), NROW(dat$newdata))
      cfg$prediction$do_predict <- is.null(dat$newdata)
      cfg$prediction$prediction_histograms <- FALSE
    },
    GJAM = {
      cfg$model$formula <- rhs
      cfg$model$notStandard <- paste(names(dat$X), collapse = ",")
      cfg$analysis$do_traits <- cfg$model$PREDICTX <- FALSE
    },
    spOccupancy = {
      cfg$model$occ.formula <- rhs
      cfg$model$det.formula <- "~ 1"
      cfg$validation_prediction_outputs$predict <- !is.null(dat$newdata)
    },
    sjSDM = cfg$model$env_formula <- rhs,
    boral = {
      cfg$model$formula.X <- rhs
      cfg$model$lv.control$num.lv <- min(2L, ncol(dat$Y) - 1L)
    }
  )
  cfg$empirical_benchmark <- list(response = "observed_incidence", validation = if (is.null(dat$test_y)) "not_supplied" else "independent_test_sites")
  cfg
}

empirical_jsdm_predictions <- function(fit, x, species) {
  # Integrate probit probabilities over a new site's independent standard-normal scores.
  design <- model.matrix(fit$model_spec$site_formula, x)
  np <- ncol(design)
  k <- as.integer(fit$model_spec$n_latent)
  if (!identical(fit$model_spec$link, "probit") || length(fit$mcmc.sp) != length(species) || k < 1L) {
    stop("Unexpected jSDM model specification for new-site prediction.", call. = FALSE)
  }
  result <- vapply(seq_along(species), function(j) {
    m <- as.matrix(fit$mcmc.sp[[j]])
    if (ncol(m) < np + k) stop("jSDM coefficient draws are missing.", call. = FALSE)
    scale <- sqrt(1 + rowSums(m[, np + seq_len(k), drop = FALSE]^2))
    rowMeans(stats::pnorm(sweep(design %*% t(m[, seq_len(np), drop = FALSE]), 2L, scale, "/")))
  }, numeric(nrow(design)))
  matrix(result, nrow(design), length(species), dimnames = list(rownames(x), species))
}

empirical_new_predictions <- function(outdir, engine, dat, exported) {
  ids <- rownames(dat$newdata)
  species <- names(dat$Y)
  long <- function(p) {
    p <- as.matrix(p)
    if (!identical(dim(p), c(length(ids), length(species)))) {
      stop("New-site predictions have unexpected dimensions.", call. = FALSE)
    }
    if (!is.null(colnames(p)) && setequal(colnames(p), species)) p <- p[, species, drop = FALSE]
    data.frame(site_id = rep(ids, times = length(species)), species = rep(species, each = length(ids)),
               predicted = as.vector(p))
  }
  read_prediction <- function(file) readRDS(file.path(outdir, "predictions", file))
  selected <- exported[exported$site_id %in% ids, c("site_id", "species", "predicted"), drop = FALSE]
  pred <- switch(engine,
    jSDM = long(empirical_jsdm_predictions(readRDS(file.path(outdir, "models", "jsdm_model_primary.rds")), dat$newdata, species)),
    GJAM = {
      manifest <- read.csv(file.path(outdir, "predictions", "gjamPredict_newdata_manifest.csv"))
      if (!identical(manifest$prediction_source, "newdata") || !identical(manifest$status, "saved")) {
        stop("GJAM did not export new-site predictions: ", manifest$reason, call. = FALSE)
      }
      obj <- read_prediction("gjamPredict.rds")
      long(obj$ypredMu %||% obj$sdList$yMu)
    },
    spOccupancy = {
      psi <- read_prediction("spOccupancy_predict.rds")$psi.0.samples
      if (length(dim(psi)) != 3L || !identical(dim(psi)[-1], c(length(species), length(ids)))) {
        stop("spOccupancy did not return the expected draw-by-species-by-site occupancy array.", call. = FALSE)
      }
      long(t(apply(psi, c(2, 3), mean)))
    },
    sjSDM = {
      m <- read.csv(file.path(outdir, "predictions", "newdata_predictions.csv"), row.names = 1, check.names = FALSE)
      long(m)
    },
    boral = {
      fit <- readRDS(file.path(outdir, "models", "boral_model.rds"))
      pred <- withr::with_seed(20260601L, boral::predict.boral(fit, newX = as.matrix(dat$newdata),
                              predict.type = "marginal", scale = "response"))
      long(if (is.list(pred)) pred$linpred else pred)
    },
    selected
  )
  if (nrow(pred) != length(ids) * length(species) || anyDuplicated(pred[c("site_id", "species")]) ||
      any(!pred$site_id %in% ids) || any(!pred$species %in% species) ||
      !all(is.finite(pred$predicted) & pred$predicted >= 0 & pred$predicted <= 1)) {
    stop("Complete, matched new-site probability predictions were not available.", call. = FALSE)
  }
  pred
}

empirical_collect <- function(master, dat, statuses) {
  out <- file.path(master, "04_unified_standard_results")
  compare <- file.path(master, "06_compare_models")
  read_engine <- function(dir, file) {
    path <- file.path(dir, "standard", file)
    if (!file.exists(path)) return(NULL)
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  decode <- function(x, keys, names) {
    key <- match(x, keys)
    ifelse(is.na(key), as.character(x), names[key])
  }
  species_names <- function(x) decode(x, dat$species_map$engine_id, dat$species_map$species)
  predictor_names <- function(x) decode(x, dat$predictor_map$engine_column, dat$predictor_map$predictor)
  effects <- assoc <- predictions <- detection <- metrics <- availability <- list()
  for (i in seq_len(nrow(statuses))) {
    engine <- statuses$engine[i]
    dir <- statuses$output_folder[i]
    state <- statuses$status[i]
    native <- read_status(dir, engine)
    reason <- if (state == "fitted") "" else paste(unlist(native$errors), collapse = "; ")
    fx <- read_engine(dir, "effects_species_environment.csv")
    fx <- fx[fx$species %in% names(dat$Y) & is.finite(fx$estimate), , drop = FALSE]
    omega <- read_engine(dir, "associations_species_species.csv")
    omega <- omega[omega$species_i %in% names(dat$Y) & omega$species_j %in% names(dat$Y) &
      is.finite(omega$estimate), , drop = FALSE]
    p <- read_engine(dir, "predictions_site_species.csv")
    train <- p[p$site_id %in% rownames(dat$Y) & p$species %in% names(dat$Y) &
      is.finite(p$predicted), c("site_id", "species", "predicted"), drop = FALSE]
    held <- NULL
    pred_note <- if (is.null(dat$newdata)) "newdata.csv not supplied" else "No fitted model"
    if (state == "fitted" && !is.null(dat$newdata)) {
      held <- tryCatch(empirical_new_predictions(dir, engine, dat, p), error = identity)
      if (inherits(held, "error")) {
        pred_note <- conditionMessage(held)
        held <- NULL
      } else pred_note <- "Complete"
    }
    if (state == "fitted") {
      fx$engine_predictor <- fx$predictor
      fx$species <- species_names(fx$species)
      fx$predictor <- predictor_names(fx$predictor)
      fx$scale <- rep(switch(engine, Hmsc = "probit_link", "Hmsc-HPC" = "probit_link",
          jSDM = "probit_link", sjSDM = "probit_model", spOccupancy = "occupancy_logit_link",
          boral = "binomial_logit_link", GJAM = "latent_GJAM"), nrow(fx))
      effects[[i]] <- fx
      omega$species_i <- species_names(omega$species_i)
      omega$species_j <- species_names(omega$species_j)
      assoc[[i]] <- omega
      if (engine == "spOccupancy") {
        alpha <- read_engine(dir, "detection_effects.csv")
        if (!is.null(alpha) && nrow(alpha)) {
          alpha$species <- species_names(alpha$species)
          detection[[length(detection) + 1L]] <- alpha
        }
      }
    }
    combine_predictions <- function(tab, scope) {
      if (is.null(tab) || !nrow(tab)) return(NULL)
      d <- if (scope == "training") as.matrix(dat$Y) else dat$test_y
      observed <- if (is.null(d)) rep(NA_integer_, nrow(tab)) else d[cbind(match(tab$site_id, rownames(d)),
                                                                            match(tab$species, colnames(d)))]
      data.frame(engine = engine,
        site_id = decode(tab$site_id, dat$site_map$engine_id, dat$site_map$site_id),
        species = species_names(tab$species), observed_detection = observed,
        predicted = tab$predicted, prediction_context = scope,
        prediction_scale = if (engine == "spOccupancy") "latent_occupancy_probability" else "observed_binary_response_probability")
    }
    scope <- if (is.null(dat$test_y)) "projection" else "held_out"
    test <- combine_predictions(held, scope)
    predictions[[i]] <- rbind(combine_predictions(train, "training"), test)
    metric_reason <- if (engine == "spOccupancy") "Occupancy probability is not the observed detection target."
      else if (is.null(dat$test_y)) "No independent test responses supplied."
      else if (is.null(test)) pred_note else ""
    brier <- auc <- NA_real_
    n_auc <- 0L
    if (!nzchar(metric_reason)) {
      brier <- mean((test$predicted - test$observed_detection)^2)
      aucs <- vapply(split(test, test$species), function(s) {
        n1 <- sum(s$observed_detection == 1L)
        n0 <- nrow(s) - n1
        if (!n1 || !n0) return(NA_real_)
        (sum(rank(s$predicted)[s$observed_detection == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
      }, numeric(1))
      n_auc <- sum(is.finite(aucs))
      if (n_auc) auc <- mean(aucs, na.rm = TRUE)
    }
    metrics[[i]] <- data.frame(engine = engine, held_out_Brier = brier, mean_species_AUC = auc,
      species_with_AUC = n_auc, evaluated_site_species = if (is.finite(brier)) nrow(test) else 0L,
      note = metric_reason)
    availability[[i]] <- data.frame(engine = engine, status = state, effects = nrow(fx),
      association_entries = nrow(omega), training_predictions = nrow(train), new_site_predictions = NROW(test),
      runtime_seconds = statuses$runtime_seconds[i], reason = reason, prediction_note = pred_note)
  }
  write_table <- function(parts, path, empty) {
    nonempty <- Filter(function(p) NROW(p) > 0L, parts)
    write.csv(if (length(nonempty)) do.call(rbind, nonempty) else empty, path, row.names = FALSE)
  }
  write_table(effects, file.path(out, "effects_all_models_long.csv"), data.frame(engine = character()))
  write_table(assoc, file.path(out, "associations_all_models_long.csv"), data.frame(engine = character()))
  write_table(predictions, file.path(out, "predictions_all_models_long.csv"), data.frame(engine = character()))
  write_table(detection, file.path(out, "detection_effects.csv"), data.frame(engine = character()))
  write.csv(do.call(rbind, metrics), file.path(compare, "universal_benchmark_comparison.csv"), row.names = FALSE)
  write.csv(do.call(rbind, availability), file.path(out, "engine_status_matrix.csv"), row.names = FALSE)
  rules <- data.frame(
    output = c("Environmental effects", "Predictions", "Occupancy", "Associations", "Runtime"),
    comparison_rule = c(
      "Report direction with model and link. Do not equate magnitudes or uncertainties across engines.",
      "Brier and mean species AUC use only complete predictions on the same held-out observations, when supplied. AUC excludes single-class species.",
      "spOccupancy estimates latent occupancy from actual visits. Its probabilities and detection effects have separate interpretations.",
      "Species-association definitions differ across engines. No common residual-association ranking or recovery metric.",
      "Descriptive wall time for the recorded engine configurations and machine."))
  write.csv(rules, file.path(compare, "comparability_rules.csv"), row.names = FALSE)
}

run_empirical_benchmark <- function(args) {
  profile <- args$profile %||% "standard"
  engines <- c("Hmsc", "Hmsc-HPC", "jSDM", "GJAM", "spOccupancy", "sjSDM", "boral")
  if (any(!args$engines %in% engines)) stop("Unrecognised empirical engine.", call. = FALSE)
  engines <- intersect(engines, args$engines)
  if (!length(engines)) stop("Select an engine.", call. = FALSE)
  dat <- load_empirical_benchmark(args$input_dir)
  output <- file.path(app_dir, "output")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  root <- tempfile("empirical_benchmark_", tmpdir = output)
  master <- file.path(root, "master")
  dirs <- c("01_benchmark_inputs", "03_engine_outputs", "04_unified_standard_results", "06_compare_models",
            "07_reports", "08_reproducible_scripts", "09_configs")
  for (dir in dirs) dir.create(file.path(master, dir), recursive = TRUE)
  file.copy(dat$paths, file.path(master, "01_benchmark_inputs"))
  config <- file.path(master, "09_configs")
  fits <- utils::modifyList(parse_args(character()), empirical_profile(profile))
  for (map in c("site_map", "species_map", "predictor_map")) {
    write.csv(dat[[map]], file.path(config, paste0(map, ".csv")), row.names = FALSE)
  }
  yaml::write_yaml(list(data_source = "empirical", profile = profile, engines = engines, sampler_settings = fits,
    training_sites = nrow(dat$Y), species = ncol(dat$Y), prediction_sites = NROW(dat$newdata),
    independent_test_observations = !is.null(dat$test_y), observed_visits = dat$visit_count,
    scope = "Nonspatial environmental baseline. No simulated truth or synthetic visits."),
    file.path(config, "universal_benchmark_config.yml"))
  writeLines(c(
    "Empirical binary-incidence comparison", "",
    "Original observations: 01_benchmark_inputs. The maps and fit settings are in 09_configs.",
    "Complete per-engine artifacts and runnable engine scripts: 03_engine_outputs.",
    "These runs do not estimate truth recovery. Interpret native diagnostics before ecological inference.",
    "spOccupancy is not applicable without repeated observed visits; occupancy is not observed incidence.",
    "Baseline: common site predictors, no phylogenetic, trait or spatial terms. Use the individual workflows for those terms.",
    "Projection notes preserve engine-specific prediction scope. Empirical Brier/AUC require disjoint observed test sites."
  ), file.path(master, "07_reports", "README.txt"))
  writeLines(c(
    "# Run from the root of the extracted master ZIP with R and the engine dependencies installed.",
    "studio <- Sys.getenv('JSDM_STUDIO_DIR')",
    "stopifnot(nzchar(studio), file.exists(file.path(studio, 'app.R')), dir.exists('01_benchmark_inputs'))",
    "runner <- file.path(studio, 'workflow_scripts', 'universal_benchmark_runner.R')",
    paste0("flags <- c('--action=empirical', '--profile=", profile, "', '--engines=", paste(engines, collapse = ","), "')"),
    "status <- system2(file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript'),",
    "                  c(shQuote(runner), flags, paste0('--input_dir=', shQuote(normalizePath('01_benchmark_inputs', winslash='/')))))",
    "quit(status = status)"
  ), file.path(master, "08_reproducible_scripts", "replay_empirical_benchmark.R"))
  log <- function(message) {
    cat(message, "\n")
    cat(message, "\n", file = file.path(master, "07_reports", "run_log.txt"), append = TRUE)
  }
  app_env <- source_app_env()
  statuses <- do.call(rbind, lapply(engines, run_engine, master_dir = master, dat = dat,
                                  app_env = app_env, log_fun = log, args = fits))
  write.csv(statuses, file.path(master, "07_reports", "engine_workflow_status.csv"), row.names = FALSE)
  empirical_collect(master, dat, statuses)
  archive <- make_zip_file(master, file.path(root, "JSDMStudio_Empirical_Benchmark.zip"))
  cat("MASTER_DIR=", normalizePath(master, winslash = "/"), "\n", sep = "")
  cat("MASTER_ZIP=", archive, "\n", sep = "")
  invisible(master)
}
