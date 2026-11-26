# Final article case-study runner for JSDM Studio.
# Produces stable, labelled CASE_STUDY outputs for manuscript figures/tables.

options(stringsAsFactors = FALSE)

if (.Platform$OS.type == "windows") {
  Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
}

args0 <- commandArgs(trailingOnly = FALSE)
script_arg <- args0[startsWith(args0, "--file=")]
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
} else {
  normalizePath("workflow_scripts/final_case_study_runner.R", winslash = "/", mustWork = FALSE)
}
app_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) stop("Cannot locate JSDM Studio app.R.", call. = FALSE)
setwd(app_dir)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.atomic(x) && is.na(x))) y else x
}

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(
    hmsc_samples = "500", hmsc_transient = "500", hmsc_thin = "2", hmsc_chains = "2",
    hmsc_large_samples = "300", hmsc_large_transient = "300",
    hmsc_bryophyte_samples = "20", hmsc_bryophyte_transient = "20",
    hpc_samples = "400", hpc_transient = "400",
    jsdm_burnin = "250", jsdm_mcmc = "500", jsdm_thin = "5",
    gjam_ng = "1000", gjam_burnin = "300",
    spocc_n_batch = "50", spocc_batch_length = "10", spocc_n_burn = "150", spocc_n_thin = "5",
    sjsdm_iter = "200", sjsdm_sampling = "1000",
    boral_n_burnin = "300", boral_n_iteration = "1200", boral_n_thin = "10",
    max_species_non_hmsc = "40",
    run_recommended = "true"
  )
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "true"
    out[[key]] <- val
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x)) return(default)
  x <- tolower(trimws(as.character(x)[1]))
  if (!nzchar(x)) return(default)
  x %in% c("true", "t", "1", "yes", "y")
}

arg_int <- function(args, key, default, min_value = NULL, max_value = NULL) {
  val <- suppressWarnings(as.integer(args[[key]] %||% default))
  if (!is.finite(val)) val <- as.integer(default)
  if (!is.null(min_value)) val <- max(as.integer(min_value), val)
  if (!is.null(max_value)) val <- min(as.integer(max_value), val)
  as.integer(val)
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

actual_root <- file.path(app_dir, "examples", "00_ACTUAL_UPLOAD_DATA")
book_root <- file.path(app_dir, "examples", "00_UPLOAD_READY_book_chapter_cases")

read_csv_auto <- function(path, row.names = 1) {
  if (!file.exists(path)) stop("Missing input file: ", path, call. = FALSE)
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(x) && !ncol(x)) return(x)
  if (!is.null(row.names) && ncol(x) >= row.names) {
    rn <- x[[row.names]]
    rn_chr <- as.character(rn)
    if (!anyDuplicated(rn_chr) && all(nzchar(rn_chr))) {
      rownames(x) <- rn_chr
      x <- x[-row.names]
    }
  }
  x
}

read_csv_optional <- function(path, row.names = 1) {
  if (!file.exists(path)) return(NULL)
  read_csv_auto(path, row.names = row.names)
}

read_text_optional <- function(path) {
  if (!file.exists(path)) return("")
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

numeric_df <- function(x) {
  if (is.null(x)) return(NULL)
  rn <- rownames(x)
  cn <- colnames(x)
  y <- as.data.frame(lapply(x, function(z) suppressWarnings(as.numeric(z))),
                     check.names = FALSE, stringsAsFactors = FALSE)
  rownames(y) <- rn
  colnames(y) <- cn
  y
}

numeric_predictors <- function(X) {
  if (is.null(X) || !ncol(X)) return(character())
  names(X)[vapply(X, function(z) {
    zz <- suppressWarnings(as.numeric(z))
    all(!is.na(zz)) && is.finite(stats::var(zz, na.rm = TRUE)) && stats::var(zz, na.rm = TRUE) > 0
  }, logical(1))]
}

valid_python_names <- function(x) {
  x[grepl("^[A-Za-z_][A-Za-z0-9_]*$", x)]
}

make_r_formula <- function(X, max_terms = 5L) {
  nms <- names(X)
  nms <- nms[nzchar(nms)]
  if (!length(nms)) return("~ 1")
  num <- numeric_predictors(X)
  keep <- unique(c(head(num, max_terms), head(setdiff(nms, num), max(0L, max_terms - length(head(num, max_terms))))))
  keep <- head(keep, max_terms)
  if (!length(keep)) return("~ 1")
  paste("~", paste(sprintf("`%s`", keep), collapse = " + "))
}

make_patsy_formula <- function(X, max_terms = 5L) {
  nms <- valid_python_names(names(X))
  nms <- nms[nzchar(nms)]
  if (!length(nms)) return("~ 1")
  num <- valid_python_names(numeric_predictors(X))
  keep <- unique(c(head(num, max_terms), head(setdiff(nms, num), max(0L, max_terms - length(head(num, max_terms))))))
  keep <- head(keep, max_terms)
  if (!length(keep)) return("~ 1")
  paste("~", paste(keep, collapse = " + "))
}

choose_coord_cols <- function(coord) {
  if (is.null(coord)) return(c("x", "y"))
  nm <- names(coord)
  low <- tolower(nm)
  pairs <- list(c("longitude", "latitude"), c("lon", "lat"), c("x", "y"), c("easting", "northing"))
  for (p in pairs) {
    hit <- nm[match(p, low)]
    if (length(hit) == 2 && all(!is.na(hit))) return(hit)
  }
  num <- nm[vapply(coord, function(z) all(is.na(z) | !is.na(suppressWarnings(as.numeric(z)))), logical(1))]
  if (length(num) >= 2) return(num[1:2])
  c(nm[1] %||% "x", nm[2] %||% "y")
}

choose_random_col <- function(study) {
  if (is.null(study) || !ncol(study)) return("sample")
  pref <- c("route", "plot", "spatial_block", "site_group", "site", "sample", names(study)[1])
  hit <- pref[pref %in% names(study)]
  hit[1] %||% names(study)[1]
}

make_reference_truth <- function(dat, label = "case_study_observed_reference") {
  if (is.null(dat$Y) || is.null(dat$X)) return(dat)
  Y <- as.data.frame(dat$Y, check.names = FALSE)
  X <- as.data.frame(dat$X, check.names = FALSE)
  species <- colnames(Y)
  predictors <- numeric_predictors(X)
  if (!length(predictors)) predictors <- head(names(X), min(3L, ncol(X)))
  eff <- expand.grid(species = species, predictor = predictors, stringsAsFactors = FALSE)
  eff$true_effect <- NA_real_
  eff$reference_method <- NA_character_
  response_kind <- function(y) {
    y <- y[is.finite(y)]
    if (length(y) && all(y %in% c(0, 1))) return("binary")
    if (length(y) && all(y >= 0) && all(abs(y - round(y)) < 1e-8)) return("count")
    "continuous"
  }
  species_kind <- setNames(vapply(Y, function(z) response_kind(suppressWarnings(as.numeric(z))), character(1)), species)
  for (i in seq_len(nrow(eff))) {
    sp <- eff$species[[i]]
    pr <- eff$predictor[[i]]
    yy <- suppressWarnings(as.numeric(Y[[sp]]))
    xx <- suppressWarnings(as.numeric(X[[pr]]))
    kind <- species_kind[[sp]]
    eff$reference_method[[i]] <- paste0(label, "_", kind, "_univariate_reference")
    ok <- is.finite(yy) & is.finite(xx)
    if (sum(ok) >= 3 && length(unique(yy[ok])) >= 2 && stats::var(xx[ok], na.rm = TRUE) > 0) {
      fit <- tryCatch(suppressWarnings(switch(
        kind,
        binary = stats::glm(yy[ok] ~ xx[ok], family = stats::binomial()),
        count = stats::glm(yy[ok] ~ xx[ok], family = stats::poisson()),
        continuous = stats::lm(yy[ok] ~ xx[ok])
      )), error = function(e) NULL)
      if (!is.null(fit) && length(stats::coef(fit)) >= 2) eff$true_effect[[i]] <- unname(stats::coef(fit)[[2]])
      if (!is.finite(eff$true_effect[[i]])) eff$true_effect[[i]] <- suppressWarnings(stats::cor(yy[ok], xx[ok], use = "pairwise.complete.obs"))
    }
  }
  pred <- as.data.frame(as.table(as.matrix(Y)), stringsAsFactors = FALSE)
  names(pred) <- c("site_id", "species", "observed")
  species_mean <- colMeans(as.matrix(Y), na.rm = TRUE)
  pred$truth_probability <- as.numeric(species_mean[as.character(pred$species)])
  pred$reference_method <- ifelse(species_kind[as.character(pred$species)] == "binary",
                                  paste0(label, "_species_prevalence_reference"),
                                  paste0(label, "_species_mean_response_reference"))
  pairs <- utils::combn(species, 2, simplify = FALSE)
  assoc <- do.call(rbind, lapply(pairs, function(p) {
    data.frame(species_i = p[[1]], species_j = p[[2]],
               true_association = suppressWarnings(stats::cor(as.numeric(Y[[p[[1]]]]),
                                                              as.numeric(Y[[p[[2]]]]),
                                                              use = "pairwise.complete.obs")),
               reference_method = paste0(label, "_observed_response_correlation"),
               stringsAsFactors = FALSE)
  }))
  dat$truth_effects <- eff
  dat$truth_predictions <- pred[, c("site_id", "species", "truth_probability", "reference_method")]
  dat$truth_associations <- assoc
  dat
}

derive_spocc_Y <- function(det, species_table = NULL) {
  det <- as.data.frame(det, check.names = FALSE)
  species <- if (!is.null(species_table) && nrow(species_table)) as.character(species_table[[1]]) else {
    unique(sub("(_rep[0-9]+|\\.rep[0-9]+)$", "", colnames(det)))
  }
  out <- matrix(NA_real_, nrow(det), length(species), dimnames = list(rownames(det), species))
  for (sp in species) {
    cols <- grep(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", sp), "(_|\\.)rep"), colnames(det), value = TRUE)
    if (!length(cols)) cols <- grep(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", sp)), colnames(det), value = TRUE)
    if (length(cols)) out[, sp] <- apply(as.matrix(det[, cols, drop = FALSE]), 1, max, na.rm = TRUE)
  }
  as.data.frame(out, check.names = FALSE)
}

load_actual_engine_dat <- function(case_folder, engine) {
  case_dir <- file.path(actual_root, case_folder)
  eng_dir <- file.path(case_dir, engine)
  if (!dir.exists(eng_dir)) stop("Missing engine folder: ", eng_dir, call. = FALSE)
  dat <- list(upload_engine_folder = eng_dir, case_input_root = case_dir)
  if (identical(engine, "Hmsc")) {
    dat$Y <- numeric_df(read_csv_auto(file.path(eng_dir, "Y.csv")))
    dat$X <- read_csv_auto(file.path(eng_dir, "XData.csv"))
    dat$traits <- read_csv_optional(file.path(eng_dir, "TrData.csv"))
    dat$study <- read_csv_optional(file.path(eng_dir, "studyDesign.csv"))
    dat$coords <- read_csv_optional(file.path(eng_dir, "coordinates.csv"))
    dat$taxonomy <- read_csv_optional(file.path(eng_dir, "taxonomy.csv"), row.names = NULL)
    dat$newdata <- read_csv_optional(file.path(eng_dir, "newdata.csv"))
    dat$newcoords <- read_csv_optional(file.path(eng_dir, "newcoords.csv"))
    dat$newstudy <- read_csv_optional(file.path(eng_dir, "prediction_studyDesign.csv"))
    dat$newick <- read_text_optional(file.path(eng_dir, "phylogeny.nwk"))
  } else if (identical(engine, "Hmsc-HPC")) {
    dat$Y <- numeric_df(read_csv_auto(file.path(eng_dir, "Y.csv")))
    dat$X <- read_csv_auto(file.path(eng_dir, "XData.csv"))
    dat$traits <- read_csv_optional(file.path(eng_dir, "traits.csv"))
    dat$study <- read_csv_optional(file.path(eng_dir, "studyDesign.csv"))
    dat$coords <- read_csv_optional(file.path(eng_dir, "coordinates.csv"))
    dat$taxonomy <- read_csv_optional(file.path(eng_dir, "taxonomy.csv"), row.names = NULL)
    dat$newdata <- read_csv_optional(file.path(eng_dir, "newdata.csv"))
    dat$newick <- read_text_optional(file.path(eng_dir, "phylo_tree.nwk"))
    if (!nzchar(dat$newick)) dat$newick <- read_text_optional(file.path(eng_dir, "phylogeny.nwk"))
  } else if (identical(engine, "jSDM")) {
    dat$Y <- numeric_df(read_csv_auto(file.path(eng_dir, "Y.csv")))
    dat$X <- read_csv_auto(file.path(eng_dir, "XData.csv"))
    dat$traits <- read_csv_optional(file.path(eng_dir, "trait_data.csv"))
    dat$trial <- numeric_df(read_csv_optional(file.path(eng_dir, "trials.csv")))
    dat$newdata <- read_csv_optional(file.path(eng_dir, "newdata.csv"))
    dat$long_format <- read_csv_optional(file.path(eng_dir, "long_format.csv"), row.names = NULL)
    dat$prediction_ids <- read_csv_optional(file.path(eng_dir, "prediction_ids.csv"), row.names = NULL)
  } else if (identical(engine, "GJAM")) {
    dat$Y <- numeric_df(read_csv_auto(file.path(eng_dir, "Y.csv")))
    dat$X <- read_csv_auto(file.path(eng_dir, "XData.csv"))
    dat$traits <- read_csv_optional(file.path(eng_dir, "specByTrait.csv"))
    dat$typeNames <- read_csv_optional(file.path(eng_dir, "typeNames.csv"), row.names = NULL)
    dat$traitTypes <- read_csv_optional(file.path(eng_dir, "traitTypes.csv"), row.names = NULL)
    dat$newdata <- read_csv_optional(file.path(eng_dir, "newdata.csv"))
    dat$effort <- read_csv_optional(file.path(eng_dir, "effort.csv"))
    dat$censor <- read_csv_optional(file.path(eng_dir, "censor.csv"), row.names = NULL)
    dat$holdoutIndex <- read_csv_optional(file.path(eng_dir, "holdoutIndex.csv"), row.names = NULL)
  } else if (identical(engine, "spOccupancy")) {
    dat$det <- numeric_df(read_csv_auto(file.path(eng_dir, "y.csv")))
    dat$occ_covs <- read_csv_auto(file.path(eng_dir, "occ.covs.csv"))
    dat$det_cov <- read_csv_auto(file.path(eng_dir, "det.covs.csv"))
    dat$coords <- read_csv_optional(file.path(eng_dir, "coords.csv"))
    dat$species_table <- read_csv_optional(file.path(eng_dir, "species.csv"), row.names = NULL)
    dat$newdata <- read_csv_optional(file.path(eng_dir, "newdata.csv"))
    dat$newcoords <- read_csv_optional(file.path(eng_dir, "newcoords.csv"))
    dat$folds <- read_csv_optional(file.path(eng_dir, "folds.csv"), row.names = NULL)
    dat$Y <- derive_spocc_Y(dat$det, dat$species_table)
    dat$X <- dat$occ_covs
  } else if (identical(engine, "sjSDM")) {
    dat$Y <- numeric_df(read_csv_auto(file.path(eng_dir, "Y.csv")))
    dat$X <- read_csv_auto(file.path(eng_dir, "env.csv"))
    dat$coords <- read_csv_optional(file.path(eng_dir, "spatial.csv"))
    dat$traits <- read_csv_optional(file.path(eng_dir, "traits.csv"))
    dat$newdata <- read_csv_optional(file.path(eng_dir, "newdata.csv"))
    dat$new_spatial <- read_csv_optional(file.path(eng_dir, "new_spatial.csv"))
  } else if (identical(engine, "boral")) {
    dat$Y <- numeric_df(read_csv_auto(file.path(eng_dir, "Y.csv")))
    dat$X <- read_csv_auto(file.path(eng_dir, "XData.csv"))
    dat$traits <- read_csv_optional(file.path(eng_dir, "traits.csv"))
    dat$rowids <- read_csv_optional(file.path(eng_dir, "row.ids.csv"))
    dat$ranefids <- read_csv_optional(file.path(eng_dir, "ranef.ids.csv"))
    dat$distmat <- numeric_df(read_csv_optional(file.path(eng_dir, "distmat.csv")))
    dat$offset <- numeric_df(read_csv_optional(file.path(eng_dir, "offset.csv")))
    dat$newdata <- read_csv_optional(file.path(eng_dir, "newdata.csv"))
    dat$trial <- numeric_df(read_csv_optional(file.path(eng_dir, "trial.size.csv")))
  } else {
    stop("Unsupported engine: ", engine, call. = FALSE)
  }
  if (is.null(dat$trial) && !is.null(dat$Y)) {
    dat$trial <- as.data.frame(matrix(1, nrow(dat$Y), ncol(dat$Y), dimnames = dimnames(as.matrix(dat$Y))),
                               check.names = FALSE)
  }
  if (is.null(dat$offset) && !is.null(dat$Y)) {
    dat$offset <- as.data.frame(matrix(0, nrow(dat$Y), ncol(dat$Y), dimnames = dimnames(as.matrix(dat$Y))),
                                check.names = FALSE)
  }
  dat <- make_reference_truth(dat, "actual_upload_case_reference")
  dat
}

load_book_ch11_full_hmsc <- function() {
  eng_dir <- file.path(book_root, "ch11_finnish_birds", "01_Hmsc_full_finland_prediction")
  dat <- list(upload_engine_folder = eng_dir, case_input_root = dirname(eng_dir))
  dat$Y <- numeric_df(read_csv_auto(file.path(eng_dir, "Y.csv")))
  dat$X <- read_csv_auto(file.path(eng_dir, "XData.csv"))
  dat$traits <- read_csv_auto(file.path(eng_dir, "TrData.csv"))
  dat$study <- read_csv_auto(file.path(eng_dir, "studyDesign.csv"))
  dat$coords <- read_csv_auto(file.path(eng_dir, "coordinates.csv"))
  dat$newdata <- read_csv_auto(file.path(eng_dir, "newdata.csv"))
  dat$newcoords <- read_csv_auto(file.path(eng_dir, "newcoords.csv"))
  dat$newstudy <- read_csv_optional(file.path(eng_dir, "prediction_studyDesign.csv"))
  dat$newick <- read_text_optional(file.path(eng_dir, "phylogeny.nwk"))
  dat$taxonomy <- NULL
  dat$trial <- as.data.frame(matrix(1, nrow(dat$Y), ncol(dat$Y), dimnames = dimnames(as.matrix(dat$Y))),
                             check.names = FALSE)
  dat$offset <- as.data.frame(matrix(0, nrow(dat$Y), ncol(dat$Y), dimnames = dimnames(as.matrix(dat$Y))),
                              check.names = FALSE)
  dat <- make_reference_truth(dat, "book_ch11_observed_reference")
  dat
}

subset_dat_species <- function(dat, max_species = 40L) {
  if (is.null(dat$Y) || ncol(dat$Y) <= max_species) return(dat)
  yy <- as.matrix(dat$Y)
  prev <- colMeans(yy, na.rm = TRUE)
  ok <- is.finite(prev) & prev > 0 & prev < 1
  candidates <- colnames(dat$Y)[ok]
  if (length(candidates) < min(max_species, ncol(dat$Y))) candidates <- colnames(dat$Y)[is.finite(prev)]
  score <- abs(prev[candidates] - 0.5)
  keep <- candidates[order(score, candidates)][seq_len(min(max_species, length(candidates)))]
  subset_cols <- function(x) {
    if (is.null(x)) return(NULL)
    hit <- intersect(keep, colnames(x))
    if (!length(hit)) return(x)
    x[, hit, drop = FALSE]
  }
  dat$Y <- subset_cols(dat$Y)
  dat$trial <- subset_cols(dat$trial)
  dat$offset <- subset_cols(dat$offset)
  if (!is.null(dat$traits) && nrow(dat$traits)) {
    hit <- intersect(keep, rownames(dat$traits))
    if (length(hit)) dat$traits <- dat$traits[hit, , drop = FALSE]
  }
  if (!is.null(dat$det) && ncol(dat$det)) {
    keep_prefix <- paste0(keep, "_rep")
    det_keep <- vapply(colnames(dat$det), function(z) any(startsWith(z, keep_prefix)), logical(1))
    if (any(det_keep)) dat$det <- dat$det[, det_keep, drop = FALSE]
  }
  dat$fit_species_subset_info <- data.frame(
    uploaded_species_n = ncol(yy), fitted_species_n = length(keep),
    subset_applied = TRUE,
    selection_rule = paste0("Final article non-Hmsc engine subset: up to ", max_species,
                            " species closest to prevalence 0.5; full upload files copied unchanged."),
    stringsAsFactors = FALSE
  )
  dat$fit_species_subset_species <- data.frame(species = keep, prevalence = as.numeric(prev[keep]),
                                               stringsAsFactors = FALSE)
  dat <- make_reference_truth(dat, "subset_observed_reference")
  dat
}

copy_upload_files <- function(outdir, dat) {
  src <- dat$upload_engine_folder
  if (is.null(src) || !dir.exists(src)) return(0L)
  files <- list.files(src, full.names = TRUE, recursive = TRUE, all.files = FALSE)
  rows <- lapply(files, function(f) {
    rel <- substring(normalizePath(f, winslash = "/", mustWork = FALSE),
                     nchar(normalizePath(src, winslash = "/", mustWork = FALSE)) + 2L)
    for (root_name in c("inputs", "data")) {
      dest <- file.path(outdir, root_name, "maximum_upload_files", rel)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      file.copy(f, dest, overwrite = TRUE)
    }
    data.frame(file = rel, size_bytes = file.info(f)$size, copied = TRUE, stringsAsFactors = FALSE)
  })
  tab <- if (length(rows)) do.call(rbind, rows) else data.frame(file = character(), size_bytes = numeric(), copied = logical())
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  write.csv(tab, file.path(outdir, "tables", "CASE_STUDY_uploaded_file_manifest.csv"), row.names = FALSE)
  writeLines(c("Final CASE_STUDY input source:",
               normalizePath(src, winslash = "/", mustWork = FALSE),
               "All upload-ready files are copied unchanged for auditability."),
             file.path(outdir, "inputs", "README_CASE_STUDY_INPUTS.txt"), useBytes = TRUE)
  nrow(tab)
}

refresh_result <- function(res, dat, case_group, scientific_question, recommended_role) {
  outdir <- res$output_folder[1]
  copied <- copy_upload_files(outdir, dat)
  if (!is.null(dat$fit_species_subset_info)) {
    write.csv(dat$fit_species_subset_info, file.path(outdir, "tables", "modeling_species_subset_summary.csv"), row.names = FALSE)
    write.csv(dat$fit_species_subset_species, file.path(outdir, "tables", "modeling_species_subset_species.csv"), row.names = FALSE)
  }
  write.csv(data.frame(
    case_group = case_group,
    scientific_question = scientific_question,
    recommended_role = recommended_role,
    engine = res$engine[1],
    status = res$status[1],
    stringsAsFactors = FALSE
  ), file.path(outdir, "tables", "CASE_STUDY_context.csv"), row.names = FALSE)
  runner_env$fill_empty_dirs(outdir, res$engine[1])
  st <- senv$safe_read_status(outdir, res$engine[1])
  app_env$ensure_output_contract(outdir, res$engine[1], st, Y = dat$Y, X = dat$X %||% dat$occ_covs)
  runner_env$standardize_engine_synonyms(outdir, res$engine[1], dat)
  zipfile <- tryCatch(runner_env$make_zip_file(outdir), error = function(e) {
    writeLines(conditionMessage(e), file.path(outdir, "diagnostics", "zip_refresh_error.txt"), useBytes = TRUE)
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
  res$case_group <- case_group
  res$scientific_question <- scientific_question
  res$recommended_role <- recommended_role
  res
}

bind_rows_fill <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  rows <- rows[vapply(rows, nrow, integer(1)) > 0]
  if (!length(rows)) return(data.frame())
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(cols, names(x))
    for (m in missing) x[[m]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

collect_case_study_standard <- function(root, summary = NULL) {
  if (is.null(summary)) {
    summary_file <- file.path(root, "CASE_STUDY_RUN_SUMMARY.csv")
    if (!file.exists(summary_file)) stop("Missing CASE_STUDY_RUN_SUMMARY.csv.", call. = FALSE)
    summary <- read.csv(summary_file, check.names = FALSE, stringsAsFactors = FALSE)
  }
  std_root <- file.path(root, "standard")
  dir.create(std_root, recursive = TRUE, showWarnings = FALSE)
  standard_files <- c(
    effects_long = "effects_long.csv",
    effects_species_environment = "effects_species_environment.csv",
    predictions_long = "predictions_long.csv",
    predictions_site_species = "predictions_site_species.csv",
    associations_long = "associations_long.csv",
    associations_species_species = "associations_species_species.csv",
    fit_metrics = "fit_metrics.csv",
    diagnostics_long = "diagnostics_long.csv",
    run_summary = "run_summary.csv",
    output_manifest = "output_manifest.csv"
  )
  for (nm in names(standard_files)) {
    rows <- list()
    for (i in seq_len(nrow(summary))) {
      outdir <- summary$output_folder[i]
      if (!nzchar(outdir %||% "") || !dir.exists(outdir)) next
      f <- file.path(outdir, "standard", standard_files[[nm]])
      if (!file.exists(f)) next
      x <- tryCatch(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(x) || !nrow(x)) next
      x$case_group <- summary$case_group[i]
      x$case_id <- summary$case_id[i]
      x$case_engine <- summary$engine[i]
      x$case_status <- summary$status[i]
      rows[[length(rows) + 1L]] <- x
    }
    write.csv(bind_rows_fill(rows), file.path(std_root, paste0("CASE_STUDY_", nm, "_all_engines.csv")), row.names = FALSE)
  }
  artifact_rows <- list()
  for (i in seq_len(nrow(summary))) {
    outdir <- summary$output_folder[i]
    if (!nzchar(outdir %||% "") || !dir.exists(outdir)) next
    files <- list.files(outdir, recursive = TRUE, full.names = FALSE)
    keep <- files[grepl("\\.(csv|png|pdf|html|rds|Rdata|txt|yml|json|zip)$", files, ignore.case = TRUE)]
    if (!length(keep)) next
    artifact_rows[[length(artifact_rows) + 1L]] <- data.frame(
      case_group = summary$case_group[i],
      engine = summary$engine[i],
      case_id = summary$case_id[i],
      status = summary$status[i],
      file = keep,
      full_path = file.path(outdir, keep),
      stringsAsFactors = FALSE
    )
  }
  artifact_index <- bind_rows_fill(artifact_rows)
  write.csv(artifact_index, file.path(std_root, "CASE_STUDY_ARTICLE_OUTPUT_INDEX.csv"), row.names = FALSE)
  write.csv(summary, file.path(std_root, "CASE_STUDY_RUN_SUMMARY_COPY.csv"), row.names = FALSE)
  md <- c(
    "# CASE STUDY Article Output Index",
    "",
    paste0("Generated: ", Sys.time()),
    paste0("Output root: ", normalizePath(root, winslash = "/", mustWork = FALSE)),
    "",
    "Use the combined `standard/CASE_STUDY_*_all_engines.csv` files for manuscript tables and figures.",
    "Use `standard/CASE_STUDY_ARTICLE_OUTPUT_INDEX.csv` to locate individual engine plots, prediction tables, reports and model objects.",
    "",
    "## Fitted Case Summary",
    "",
    paste(capture.output(print(summary[, c("case_group", "engine", "status", "file_count", "missing_count", "zip_size_kb", "runtime_seconds")], row.names = FALSE)), collapse = "\n")
  )
  writeLines(md, file.path(root, "CASE_STUDY_ARTICLE_OUTPUT_INDEX.md"), useBytes = TRUE)
  invisible(file.path(std_root, "CASE_STUDY_ARTICLE_OUTPUT_INDEX.csv"))
}

safe_add <- function(results, expr, dat, case_group, scientific_question, recommended_role, engine, case_id) {
  idx <- length(results) + 1L
  results[[idx]] <- tryCatch(
    refresh_result(force(expr), dat, case_group, scientific_question, recommended_role),
    error = function(e) {
      data.frame(engine = engine, case_id = case_id, status = "fit_failed",
                 output_folder = NA_character_, zip = NA_character_, zip_size_kb = NA_real_,
                 file_count = NA_integer_, missing_count = NA_integer_, missing_items = "",
                 failure_reason = conditionMessage(e), warnings = "",
                 notes = "Unhandled final case-study runner error.",
                 runtime_seconds = NA_real_, maximum_upload_files_copied = NA_integer_,
                 case_group = case_group, scientific_question = scientific_question,
                 recommended_role = recommended_role, stringsAsFactors = FALSE)
    }
  )
  results
}

configure_hmsc <- function(dat, args, large = FALSE, spatial = FALSE, full_bird = FALSE) {
  cfg <- runner_env$hmsc_cfg()
  xy <- choose_coord_cols(dat$coords)
  group_col <- choose_random_col(dat$study)
  cfg$model$distr <- "probit"
  cfg$model$XFormula <- if (full_bird) "~ hab + poly(clim, degree = 2, raw = TRUE)" else make_r_formula(dat$X, max_terms = 5L)
  cfg$model$TrFormula <- if (full_bird) "~ Migration + LogMass" else make_r_formula(dat$traits %||% data.frame(), max_terms = 5L)
  cfg$model$use_traits <- !is.null(dat$traits) && nrow(dat$traits) == ncol(dat$Y) && ncol(dat$traits) > 0
  cfg$model$use_phylogeny <- nzchar(dat$newick %||% "")
  cfg$data$phylogeny <- if (cfg$model$use_phylogeny) "phylogeny.nwk" else ""
  cfg$model$random_effect_column <- group_col
  cfg$model$lon_col <- xy[1]
  cfg$model$lat_col <- xy[2]
  cfg$model$random_mode <- if (spatial) "spatial_full" else if (!is.null(dat$study) && group_col %in% names(dat$study)) "sample" else "none"
  cfg$model$spatial_method <- if (spatial) "Full" else "Full"
  cfg$model$nfMin <- 1L
  cfg$model$nfMax <- min(if (large) 20L else 10L, max(1L, ncol(dat$Y) - 1L))
  cfg$mcmc$samples <- if (large) arg_int(args, "hmsc_large_samples", 300L, 1L) else arg_int(args, "hmsc_samples", 500L, 1L)
  cfg$mcmc$transient <- if (large) arg_int(args, "hmsc_large_transient", 300L, 0L) else arg_int(args, "hmsc_transient", 500L, 0L)
  cfg$mcmc$thin <- arg_int(args, "hmsc_thin", 2L, 1L)
  cfg$mcmc$nChains <- arg_int(args, "hmsc_chains", 2L, 1L)
  cfg$mcmc$nParallel <- 1L
  cfg$mcmc$verbose <- 0L
  cfg$mcmc$preset <- "Final article case-study stable run"
  cfg$outputs$predicted <- TRUE
  cfg$outputs$fit <- TRUE
  cfg$outputs$cv <- FALSE
  cfg$outputs$waic <- FALSE
  cfg$outputs$parameters <- TRUE
  cfg$outputs$variance_partitioning <- TRUE
  cfg$outputs$omega <- TRUE
  cfg$outputs$gradients <- TRUE
  cfg$outputs$plotting$plotBeta <- TRUE
  cfg$outputs$plotting$plotGamma <- TRUE
  cfg$outputs$plotting$plotTree <- cfg$model$use_phylogeny
  cfg$outputs$predictions$nfolds <- 2L
  cfg$outputs$predictions$partition_column <- group_col
  cfg$outputs$predictions$species.list <- if (full_bird && "Corvus_monedula" %in% colnames(dat$Y)) "Corvus_monedula" else colnames(dat$Y)[1]
  cfg$outputs$predictions$trait.list <- if (full_bird && "LogMass" %in% names(dat$traits)) "LogMass" else ""
  cfg$outputs$predictions$env.list <- if (full_bird) "clim" else paste(head(names(dat$X), 2), collapse = ",")
  cfg
}

run_hmsc_with_prediction <- function(root, case_id, dat, cfg, notes = "") {
  started <- Sys.time()
  prep <- senv$prepare_case(root, "Hmsc", case_id, dat, Y = dat$Y)
  outdir <- normalizePath(prep$outdir, winslash = "/", mustWork = FALSE)
  dat_case <- prep$dat
  new_study_use <- dat$newstudy
  if (tolower(as.character(cfg$model$random_mode %||% "")) %in% c("spatial", "spatial_full", "spatial_nngp", "spatial_gpp")) {
    new_study_use <- NULL
  }
  if (!is.null(dat$traits)) senv$write_pair(dat$traits, outdir, "TrData.csv")
  if (!is.null(dat$study)) senv$write_pair(dat$study, outdir, "studyDesign.csv")
  if (!is.null(dat$coords)) senv$write_pair(dat$coords, outdir, "coordinates.csv")
  if (!is.null(dat$taxonomy)) senv$write_pair(dat$taxonomy, outdir, "taxonomy.csv", row.names = FALSE)
  if (nzchar(dat$newick %||% "")) senv$write_text_pair(dat$newick, outdir, "phylogeny.nwk")
  if (!is.null(dat$newdata)) senv$write_pair(dat$newdata, outdir, "newdata.csv")
  if (!is.null(dat$newcoords)) senv$write_pair(dat$newcoords, outdir, "newcoords.csv")
  if (!is.null(new_study_use)) senv$write_pair(new_study_use, outdir, "prediction_studyDesign.csv")
  yaml::write_yaml(cfg, file.path(outdir, "used_config.yml"))
  check <- app_env$validate_hmsc(dat_case$Y, dat_case$X,
                                 if (isTRUE(cfg$model$use_traits)) dat$traits else NULL,
                                 dat$study, dat$coords,
                                 distr = cfg$model$distr,
                                 XFormula = cfg$model$XFormula,
                                 TrFormula = cfg$model$TrFormula,
                                 use_traits = isTRUE(cfg$model$use_traits),
                                 use_phylogeny = isTRUE(cfg$model$use_phylogeny),
                                 random_mode = cfg$model$random_mode,
                                 random_effect_column = cfg$model$random_effect_column,
                                 spatial_method = cfg$model$spatial_method,
                                 nNeighbours = cfg$model$nNeighbours,
                                 lon_col = cfg$model$lon_col,
                                 lat_col = cfg$model$lat_col,
                                 samples = cfg$mcmc$samples,
                                 transient = cfg$mcmc$transient,
                                 thin = cfg$mcmc$thin,
                                 nChains = cfg$mcmc$nChains,
                                 nParallel = cfg$mcmc$nParallel,
                                 nfMin = cfg$model$nfMin,
                                 nfMax = cfg$model$nfMax,
                                 nfolds = cfg$outputs$predictions$nfolds,
                                 partition_column = cfg$outputs$predictions$partition_column)
  app_env$write_data_check_messages(outdir, check$messages)
  if (!isTRUE(check$ok)) {
    senv$mark_case_status(outdir, "Hmsc", "check_failed", warnings = check$messages,
                          errors = "Hmsc final case-study check failed.",
                          Y = dat_case$Y, X = dat_case$X)
  } else {
    tryCatch({
      app_env$run_hmsc_s1s7_pipeline(outdir, cfg, dat_case$Y, dat_case$X,
                                     TrData = if (isTRUE(cfg$model$use_traits)) dat$traits else NULL,
                                     studyDesign = dat$study,
                                     coord = dat$coords,
                                     newdata = dat$newdata,
                                     newStudyDesign = new_study_use,
                                     newcoord = dat$newcoords,
                                     log_fun = message)
    }, error = function(e) {
      senv$mark_case_status(outdir, "Hmsc", "fit_failed", warnings = check$messages,
                            errors = conditionMessage(e), Y = dat_case$Y, X = dat_case$X)
    })
  }
  senv$finish_case("Hmsc", case_id, outdir, dat_case, started, notes)
}

configure_hpc <- function(dat, args, spatial = FALSE) {
  cfg <- runner_env$hmschpc_cfg(app_env)
  xy <- choose_coord_cols(dat$coords)
  group_col <- choose_random_col(dat$study)
  cfg$model$distribution <- "probit"
  cfg$model$XFormula <- make_patsy_formula(dat$X, max_terms = 5L)
  cfg$model$trait_formula <- make_patsy_formula(dat$traits %||% data.frame(), max_terms = 5L)
  cfg$model$use_traits <- !is.null(dat$traits) && nrow(dat$traits) == ncol(dat$Y) && ncol(dat$traits) > 0
  cfg$model$phylogeny_mode <- if (nzchar(dat$newick %||% "")) "newick" else if (!is.null(dat$taxonomy)) "taxonomy" else "none"
  cfg$random_effects$mode <- if (spatial) "spatial_full" else if (!is.null(dat$study) && group_col %in% names(dat$study)) "iid" else "none"
  cfg$random_effects$column <- group_col
  cfg$random_effects$name <- if (spatial) "spatial_site" else paste0(group_col, "_iid")
  cfg$random_effects$coord_x <- xy[1]
  cfg$random_effects$coord_y <- xy[2]
  cfg$random_effects$nf <- 1L
  cfg$random_effects$nfMin <- 1L
  cfg$random_effects$nfMax <- min(10L, max(1L, ncol(dat$Y) - 1L))
  cfg$sampler$samples <- arg_int(args, "hpc_samples", 400L, 1L)
  cfg$sampler$transient <- arg_int(args, "hpc_transient", 400L, 0L)
  cfg$sampler$thin <- 1L
  cfg$sampler$chains <- 1L
  cfg$sampler$verbose <- 1L
  cfg$runtime$python <- app_env$hmschpc_python_bin(list())
  cfg$runtime$python_source <- app_env$hmschpc_python_source_dir()
  cfg
}

configure_jsdm <- function(dat, args) {
  cfg <- runner_env$jsdm_cfg()
  cfg$model$model_type <- "binomial_probit"
  cfg$model$site_formula <- "~ ."
  cfg$model$trait_formula <- "~ ."
  cfg$model$allow_traits <- FALSE
  cfg$model$site_effect <- "none"
  cfg$model$n_latent <- min(2L, max(0L, ncol(dat$Y) - 1L))
  cfg$mcmc$burnin <- arg_int(args, "jsdm_burnin", 250L, 0L)
  cfg$mcmc$mcmc <- arg_int(args, "jsdm_mcmc", 500L, 1L)
  cfg$mcmc$thin <- arg_int(args, "jsdm_thin", 5L, 1L)
  cfg$prediction$do_predict <- TRUE
  cfg
}

configure_gjam <- function(dat, args) {
  cfg <- runner_env$gjam_cfg(colnames(dat$Y))
  cfg$model$formula <- make_r_formula(dat$X, max_terms = 4L)
  type_col <- intersect(c("type", "typeName"), names(dat$typeNames))[1]
  type_vec <- if (!is.na(type_col)) as.character(dat$typeNames[[type_col]]) else rep("PA", ncol(dat$Y))
  if (length(type_vec) != ncol(dat$Y)) type_vec <- rep(type_vec[1] %||% "PA", ncol(dat$Y))
  cfg$model$type_single <- type_vec[1]
  cfg$model$typeNames_text <- paste(type_vec, collapse = ",")
  cfg$model$ng <- arg_int(args, "gjam_ng", 1000L, 10L)
  cfg$model$burnin <- arg_int(args, "gjam_burnin", 300L, 0L)
  cfg$analysis$do_traits <- !is.null(dat$traits) && ncol(dat$traits) > 0
  cfg$analysis$do_sensitivity <- TRUE
  cfg$analysis$do_predict <- TRUE
  cfg
}

trim_gjam_traits <- function(dat) {
  if (is.null(dat$traits) || !ncol(dat$traits)) return(dat)
  keep <- seq_len(min(3L, ncol(dat$traits)))
  dat$traits <- dat$traits[, keep, drop = FALSE]
  dat
}

make_gjam_safe_names <- function(nms) {
  x <- gsub("[^A-Za-z0-9]+", "", as.character(nms))
  x[!nzchar(x)] <- paste0("V", which(!nzchar(x)))
  x <- ifelse(grepl("^[A-Za-z]", x), x, paste0("V", x))
  make.unique(x, sep = "")
}

sanitize_gjam_predictors <- function(dat) {
  if (is.null(dat$X) || !ncol(dat$X)) return(dat)
  old <- names(dat$X)
  safe <- make_gjam_safe_names(old)
  names(dat$X) <- safe
  if (!is.null(dat$newdata) && ncol(dat$newdata)) {
    hit <- match(old, names(dat$newdata))
    if (all(!is.na(hit))) names(dat$newdata)[hit] <- safe
  }
  dat$predictor_name_map <- data.frame(original_predictor = old, gjam_formula_name = safe,
                                       stringsAsFactors = FALSE)
  dat
}

configure_spocc <- function(dat, args) {
  cfg <- runner_env$spocc_cfg()
  cfg$model$model_type <- "msPGOcc"
  cfg$model$data_structure <- "multi-species replicated"
  cfg$model$occ.formula <- make_r_formula(dat$occ_covs, max_terms = 4L)
  det_names <- names(dat$det_cov)
  if (length(det_names) && all(grepl("(__rep|_rep|\\.rep)[0-9]+$", det_names))) {
    cfg$model$det.formula <- "~ obs"
  } else {
    cfg$model$det.formula <- make_r_formula(dat$det_cov, max_terms = 2L)
  }
  cfg$model$formula <- cfg$model$occ.formula
  cfg$mcmc$n.batch <- arg_int(args, "spocc_n_batch", 50L, 1L)
  cfg$mcmc$batch.length <- arg_int(args, "spocc_batch_length", 10L, 1L)
  cfg$mcmc$n.burn <- arg_int(args, "spocc_n_burn", 150L, 0L)
  cfg$mcmc$n.thin <- arg_int(args, "spocc_n_thin", 5L, 1L)
  cfg$mcmc$n.chains <- 1L
  cfg$mcmc$n.report <- 10L
  cfg$spatial_latent_svc$n.factors <- min(2L, max(1L, ncol(dat$Y) - 1L))
  cfg$spatial_latent_svc$NNGP <- TRUE
  cfg$spatial_latent_svc$n.neighbors <- 3L
  cfg$validation_prediction_outputs$waicOcc <- TRUE
  cfg$validation_prediction_outputs$predict <- !is.null(dat$newdata)
  cfg$validation_prediction_outputs$save_plots <- TRUE
  cfg
}

configure_sjsdm <- function(dat, args) {
  cfg <- runner_env$sjsdm_cfg()
  cfg$model$family <- "binomial_probit"
  cfg$model$env_model <- "linear"
  cfg$model$env_formula <- "~ ."
  cfg$model$spatial_model <- if (!is.null(dat$coords)) "linear" else "none"
  cfg$model$spatial_formula <- "~ 0 + ."
  cfg$model$iter <- arg_int(args, "sjsdm_iter", 200L, 1L)
  cfg$model$sampling <- arg_int(args, "sjsdm_sampling", 1000L, 1L)
  cfg$model$step_size <- min(32L, max(4L, nrow(dat$Y)))
  cfg$model$parallel <- 0L
  cfg$dnn_optimizer$device <- "cpu"
  cfg$outputs$plots <- TRUE
  cfg
}

configure_boral <- function(dat, args) {
  cfg <- runner_env$boral_cfg()
  cfg$model$family <- "binomial"
  cfg$model$formula.X <- make_r_formula(dat$X, max_terms = 4L)
  cfg$model$trial.size <- 1L
  cfg$model$lv.control$num.lv <- min(2L, max(1L, ncol(dat$Y) - 1L))
  cfg$model$row.eff <- "none"
  cfg$traits_random_ssvs$use_traits <- FALSE
  cfg$mcmc_prior$n.burnin <- arg_int(args, "boral_n_burnin", 300L, 0L)
  cfg$mcmc_prior$n.iteration <- arg_int(args, "boral_n_iteration", 1200L, 1L)
  cfg$mcmc_prior$n.thin <- arg_int(args, "boral_n_thin", 10L, 1L)
  cfg$mcmc_prior$calc.ics <- FALSE
  cfg
}

write_dependency_preflight <- function(root) {
  pkgs <- c("Hmsc", "gjam", "spOccupancy", "sjSDM", "boral", "jSDM", "rjags", "R2jags",
            "ape", "ggplot2", "yaml", "jsonlite", "zip")
  tab <- data.frame(
    package = pkgs,
    available = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE),
    version = vapply(pkgs, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "", character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(tab, file.path(root, "CASE_STUDY_dependency_preflight.csv"), row.names = FALSE)
  tab
}

run_final_cases <- function() {
  args <- parse_cli()
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  root <- file.path(app_dir, "output", paste0("CASE_STUDIES_FINAL_", timestamp))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  preflight <- write_dependency_preflight(root)

  case_plan <- data.frame(
    case_group = c(
      "CASE_STUDY_01_VascularPlants",
      "CASE_STUDY_02_DeadwoodFungi",
      "CASE_STUDY_03_Birds_FinlandPrediction",
      "CASE_STUDY_04_Bryophytes_TibetSpatial",
      "CASE_STUDY_05_FreshwaterMacroinvertebrates"
    ),
    taxonomic_group = c("Vascular plants", "Deadwood fungi", "Birds", "Bryophytes", "Freshwater macroinvertebrates"),
    input_case = c(
      "01_vascular_plants__ch06_plant_traits_whittaker",
      "02_deadwood_fungi__ch07_deadwood_fungi",
      "ch11_finnish_birds/full_finland_prediction",
      "04_bryophytes__tibet_bryophyte_30site_spatial_example",
      "05_freshwater_macroinvertebrates__cestes_freshwater_macroinvertebrates_doledec1996"
    ),
    scientific_question = c(
      "Trait filtering and environmental selection in plant communities",
      "Substrate-driven fungal community assembly on deadwood",
      "Spatial prediction of Finnish bird occurrence and richness",
      "Alpine bryophyte spatial community assembly across Tibetan environmental gradients",
      "Freshwater metacommunity structure with occupancy and detection-aware comparison"
    ),
    hmsc_role = c(
      "Primary trait/phylogeny baseline",
      "Primary hierarchical baseline",
      "Primary full prediction workflow",
      "Primary high-diversity trait workflow; phylogeny retained as input",
      "Primary hierarchical baseline"
    ),
    recommended_engines = c("jSDM; boral", "GJAM; boral", "Hmsc-HPC", "Hmsc-HPC; sjSDM", "spOccupancy; GJAM"),
    stringsAsFactors = FALSE
  )
  write.csv(case_plan, file.path(root, "CASE_STUDY_PLAN.csv"), row.names = FALSE)

  all_results <- list()

  run_hmsc_actual <- function(plan_row, actual_folder, spatial = FALSE, large = FALSE,
                              disable_phylogeny = FALSE, random_init = FALSE,
                              samples_override = NULL, transient_override = NULL,
                              light_convergence = FALSE) {
    dat <- load_actual_engine_dat(actual_folder, "Hmsc")
    cfg <- configure_hmsc(dat, args, large = large, spatial = spatial, full_bird = FALSE)
    if (isTRUE(disable_phylogeny)) {
      cfg$model$use_phylogeny <- FALSE
      cfg$data$phylogeny <- ""
      cfg$outputs$plotting$plotTree <- FALSE
    }
    if (isTRUE(random_init)) cfg$mcmc$initPar <- ""
    if (!is.null(samples_override)) cfg$mcmc$samples <- as.integer(samples_override)
    if (!is.null(transient_override)) cfg$mcmc$transient <- as.integer(transient_override)
    if (isTRUE(light_convergence)) {
      cfg$outputs$convergence$showBeta <- FALSE
      cfg$outputs$convergence$showGamma <- FALSE
      cfg$outputs$convergence$showOmega <- FALSE
      cfg$outputs$convergence$showRho <- FALSE
      cfg$outputs$convergence$showAlpha <- FALSE
    }
    case_id <- paste0(plan_row$case_group, "__Hmsc_final")
    all_results <<- safe_add(all_results,
      run_hmsc_with_prediction(root, case_id, dat, cfg,
                               paste0(plan_row$hmsc_role, "; final case-study Hmsc run with prediction outputs where newdata is available.",
                                      if (isTRUE(disable_phylogeny)) " Hmsc phylogeny was disabled by this runner configuration; the raw tree or taxonomy input is still archived with inputs." else "",
                                      if (isTRUE(random_init)) " Random initial values were used because fixed-effect SSDM initialisation is unstable for this sparse high-diversity matrix." else "",
                                      if (isTRUE(light_convergence)) " Detailed per-parameter convergence plots were disabled for this high-dimensional demonstration to keep the run reproducible on desktop memory." else "")),
      dat, plan_row$case_group, plan_row$scientific_question, plan_row$hmsc_role, "Hmsc", case_id)
  }

  # Hmsc is run for every final article case.
  run_hmsc_actual(case_plan[1, ], case_plan$input_case[1], spatial = FALSE, large = FALSE)
  run_hmsc_actual(case_plan[2, ], case_plan$input_case[2], spatial = FALSE, large = FALSE)

  bird_dat <- load_book_ch11_full_hmsc()
  bird_cfg <- configure_hmsc(bird_dat, args, large = TRUE, spatial = TRUE, full_bird = TRUE)
  bird_id <- "CASE_STUDY_03_Birds_FinlandPrediction__Hmsc_full_book_prediction"
  all_results <- safe_add(all_results,
    run_hmsc_with_prediction(root, bird_id, bird_dat, bird_cfg,
                             "Full book Chapter 11 Hmsc run: 137 routes, 50 bird species and 9454 Finland prediction grid rows."),
    bird_dat, case_plan$case_group[3], case_plan$scientific_question[3],
    case_plan$hmsc_role[3], "Hmsc", bird_id)

  run_hmsc_actual(case_plan[4, ], case_plan$input_case[4], spatial = FALSE, large = TRUE,
                  disable_phylogeny = FALSE, random_init = TRUE,
                  samples_override = arg_int(args, "hmsc_bryophyte_samples", 20L, 1L),
                  transient_override = arg_int(args, "hmsc_bryophyte_transient", 20L, 0L),
                  light_convergence = TRUE)
  run_hmsc_actual(case_plan[5, ], case_plan$input_case[5], spatial = TRUE, large = FALSE)

  if (as_bool(args$run_recommended, TRUE)) {
    max_sp <- arg_int(args, "max_species_non_hmsc", 40L, 1L)

    # Vascular plants: jSDM and boral as lighter comparison/ordination engines.
    dat <- subset_dat_species(load_actual_engine_dat(case_plan$input_case[1], "jSDM"), max_sp)
    cfg <- configure_jsdm(dat, args)
    cid <- "CASE_STUDY_01_VascularPlants__jSDM_final"
    all_results <- safe_add(all_results,
      senv$run_jsdm_case(root, cid, dat, cfg, dat$Y, trials = dat$trial,
                         notes = "Recommended comparison engine for plant occurrence effects."),
      dat, case_plan$case_group[1], case_plan$scientific_question[1], "Recommended comparison engine", "jSDM", cid)

    dat <- subset_dat_species(load_actual_engine_dat(case_plan$input_case[1], "boral"), max_sp)
    cfg <- configure_boral(dat, args)
    cid <- "CASE_STUDY_01_VascularPlants__boral_final"
    all_results <- safe_add(all_results,
      senv$run_boral_case(root, cid, dat, cfg, dat$Y,
                          notes = "Recommended ordination/residual-association engine for plant occurrence data."),
      dat, case_plan$case_group[1], case_plan$scientific_question[1], "Recommended ordination engine", "boral", cid)

    # Deadwood fungi: GJAM and boral.
    dat <- subset_dat_species(sanitize_gjam_predictors(trim_gjam_traits(load_actual_engine_dat(case_plan$input_case[2], "GJAM"))), max_sp)
    cfg <- configure_gjam(dat, args)
    type_col <- intersect(c("type", "typeName"), names(dat$typeNames))[1]
    type_vec <- if (!is.na(type_col)) as.character(dat$typeNames[[type_col]]) else rep("PA", ncol(dat$Y))
    if (length(type_vec) != ncol(dat$Y)) type_vec <- rep(type_vec[1] %||% "PA", ncol(dat$Y))
    cid <- "CASE_STUDY_02_DeadwoodFungi__GJAM_final"
    all_results <- safe_add(all_results,
      senv$run_gjam_case(root, cid, dat, cfg, dat$Y, type_vec,
                         notes = "Recommended mixed/observation-scale engine for fungal community response."),
      dat, case_plan$case_group[2], case_plan$scientific_question[2], "Recommended observation-scale engine", "GJAM", cid)

    dat <- subset_dat_species(load_actual_engine_dat(case_plan$input_case[2], "boral"), max_sp)
    cfg <- configure_boral(dat, args)
    cid <- "CASE_STUDY_02_DeadwoodFungi__boral_final"
    all_results <- safe_add(all_results,
      senv$run_boral_case(root, cid, dat, cfg, dat$Y,
                          notes = "Recommended ordination/residual-association engine for fungal occurrence data."),
      dat, case_plan$case_group[2], case_plan$scientific_question[2], "Recommended ordination engine", "boral", cid)

    # Birds: Hmsc-HPC CPU subset on upload-ready reduced input.
    dat <- subset_dat_species(load_actual_engine_dat("03_birds__ch11_finnish_birds", "Hmsc-HPC"), max_sp)
    cfg <- configure_hpc(dat, args, spatial = FALSE)
    cid <- "CASE_STUDY_03_Birds__HmscHPC_cpu_reduced_final"
    all_results <- safe_add(all_results,
      senv$run_hmschpc_case(root, cid, dat, cfg, dat$Y,
                            notes = "Recommended CPU pyhmsc/HDF5 comparison branch for the reduced bird upload-ready case."),
      dat, case_plan$case_group[3], case_plan$scientific_question[3], "Recommended CPU Hmsc-HPC comparison branch", "Hmsc-HPC", cid)

    # Bryophytes: Hmsc-HPC and sjSDM, with species subset for non-Hmsc stability.
    dat <- subset_dat_species(load_actual_engine_dat(case_plan$input_case[4], "Hmsc-HPC"), max_sp)
    cfg <- configure_hpc(dat, args, spatial = FALSE)
    cid <- "CASE_STUDY_04_Bryophytes__HmscHPC_cpu_subset_final"
    all_results <- safe_add(all_results,
      senv$run_hmschpc_case(root, cid, dat, cfg, dat$Y,
                            notes = "Recommended CPU Hmsc-HPC branch with Newick phylogeny retained; full 190-species upload copied, fitted subset used for stable final run and the tree is pruned internally to the fitted species."),
      dat, case_plan$case_group[4], case_plan$scientific_question[4], "Recommended CPU Hmsc-HPC branch", "Hmsc-HPC", cid)

    dat <- subset_dat_species(load_actual_engine_dat(case_plan$input_case[4], "sjSDM"), max_sp)
    cfg <- configure_sjsdm(dat, args)
    cid <- "CASE_STUDY_04_Bryophytes__sjSDM_cpu_subset_final"
    all_results <- safe_add(all_results,
      senv$run_sjsdm_case(root, cid, dat, cfg, dat$Y,
                          notes = "Recommended fast neural/biotic-structure branch; full 190-species upload copied, fitted subset used for stable final run."),
      dat, case_plan$case_group[4], case_plan$scientific_question[4], "Recommended fast high-diversity comparison branch", "sjSDM", cid)

    # Freshwater macroinvertebrates: detection-aware spOccupancy and GJAM comparison.
    dat <- subset_dat_species(load_actual_engine_dat(case_plan$input_case[5], "spOccupancy"), max_sp)
    cfg <- configure_spocc(dat, args)
    sp_tab <- data.frame(species = colnames(dat$Y), stringsAsFactors = FALSE)
    cid <- "CASE_STUDY_05_FreshwaterMacroinvertebrates__spOccupancy_final"
    all_results <- safe_add(all_results,
      senv$run_spocc_case(root, cid, dat, cfg, dat$det,
                          species_table = sp_tab,
                          notes = "Recommended detection-aware occupancy branch for replicated survey input."),
      dat, case_plan$case_group[5], case_plan$scientific_question[5], "Recommended detection-aware occupancy branch", "spOccupancy", cid)

    dat <- subset_dat_species(sanitize_gjam_predictors(trim_gjam_traits(load_actual_engine_dat(case_plan$input_case[5], "GJAM"))), max_sp)
    cfg <- configure_gjam(dat, args)
    type_col <- intersect(c("type", "typeName"), names(dat$typeNames))[1]
    type_vec <- if (!is.na(type_col)) as.character(dat$typeNames[[type_col]]) else rep("PA", ncol(dat$Y))
    if (length(type_vec) != ncol(dat$Y)) type_vec <- rep(type_vec[1] %||% "PA", ncol(dat$Y))
    cid <- "CASE_STUDY_05_FreshwaterMacroinvertebrates__GJAM_final"
    all_results <- safe_add(all_results,
      senv$run_gjam_case(root, cid, dat, cfg, dat$Y, type_vec,
                         notes = "Recommended observation-scale comparison branch for freshwater community response."),
      dat, case_plan$case_group[5], case_plan$scientific_question[5], "Recommended observation-scale comparison branch", "GJAM", cid)
  }

  summary <- do.call(rbind, all_results)
  rownames(summary) <- NULL
  write.csv(summary, file.path(root, "CASE_STUDY_RUN_SUMMARY.csv"), row.names = FALSE)

  by_case <- aggregate(engine ~ case_group + status, summary, length)
  names(by_case)[names(by_case) == "engine"] <- "n_engine_runs"
  by_engine <- aggregate(case_id ~ engine + status, summary, length)
  names(by_engine)[names(by_engine) == "case_id"] <- "n_case_runs"
  write.csv(by_case, file.path(root, "CASE_STUDY_STATUS_BY_CASE.csv"), row.names = FALSE)
  write.csv(by_engine, file.path(root, "CASE_STUDY_STATUS_BY_ENGINE.csv"), row.names = FALSE)
  collect_case_study_standard(root, summary)

  bad <- summary[!(summary$status %in% c("fitted", "model_defined")) | summary$missing_count > 0, , drop = FALSE]
  report <- c(
    "# JSDM Studio Final Article Case-Study Outputs",
    "",
    paste0("Generated: ", Sys.time()),
    paste0("Output root: ", normalizePath(root, winslash = "/", mustWork = FALSE)),
    "",
    "## Scope",
    "",
    "This run produces labelled CASE_STUDY output folders for manuscript case-study figures and tables.",
    "Hmsc is run for every selected case. Recommended additional engines are run where they are scientifically appropriate and stable.",
    "Full upload-ready input folders are copied into each engine ZIP under inputs/maximum_upload_files and data/maximum_upload_files.",
    "",
    "## Dependency Preflight",
    "",
    paste(capture.output(print(preflight, row.names = FALSE)), collapse = "\n"),
    "",
    "## Case Plan",
    "",
    paste(capture.output(print(case_plan, row.names = FALSE)), collapse = "\n"),
    "",
    "## Status By Case",
    "",
    paste(capture.output(print(by_case, row.names = FALSE)), collapse = "\n"),
    "",
    "## Status By Engine",
    "",
    paste(capture.output(print(by_engine, row.names = FALSE)), collapse = "\n"),
    "",
    "## Non-Fitted Or Incomplete Outputs",
    "",
    if (nrow(bad)) paste(capture.output(print(bad[, c("engine", "case_id", "status", "missing_count", "failure_reason")], row.names = FALSE)), collapse = "\n") else "None.",
    "",
    "## Full Output Summary",
    "",
    paste(capture.output(print(summary[, c("case_group", "engine", "case_id", "status", "file_count", "missing_count", "zip_size_kb", "output_folder", "zip")], row.names = FALSE)), collapse = "\n")
  )
  writeLines(report, file.path(root, "CASE_STUDY_FINAL_REPORT.md"), useBytes = TRUE)

  master_zip <- tryCatch(runner_env$make_zip_file(root, file.path(dirname(root), paste0(basename(root), "_ALL_OUTPUTS.zip"))),
                         error = function(e) {
                           writeLines(conditionMessage(e), file.path(root, "CASE_STUDY_master_zip_error.txt"), useBytes = TRUE)
                           NA_character_
                         })
  cat("CASE_STUDY_ROOT=", normalizePath(root, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("CASE_STUDY_MASTER_ZIP=", if (!is.na(master_zip)) normalizePath(master_zip, winslash = "/", mustWork = FALSE) else NA_character_, "\n", sep = "")
  cat("CASE_STUDY_SUMMARY=", normalizePath(file.path(root, "CASE_STUDY_RUN_SUMMARY.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  invisible(summary)
}

if (!interactive() && identical(sys.nframe(), 0L)) run_final_cases()
