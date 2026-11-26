# Maximum-upload quick audit for public integrated JSDM Studio example cases.
# The script gives every case a readable name, maps the common upload files to
# each engine's native input names, runs one stable branch per engine, and writes
# honest status/diagnostics/ZIP outputs.

options(stringsAsFactors = FALSE)

if (.Platform$OS.type == "windows") {
  Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
}

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_arg <- args0[startsWith(args0, file_arg)]
script_path <- if (length(script_arg)) {
  normalizePath(sub(file_arg, "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
} else {
  normalizePath("workflow_scripts/public_maximum_upload_case_suite.R", winslash = "/", mustWork = FALSE)
}
app_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) stop("Cannot locate JSDM Studio app.R.", call. = FALSE)
setwd(app_dir)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(cases = "all",
              engines = "Hmsc,Hmsc-HPC,jSDM,GJAM,spOccupancy,sjSDM,boral",
              max_cases = NA_character_)
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "true"
    out[[key]] <- val
  }
  out$engines <- trimws(unlist(strsplit(out$engines, ",")))
  out$engines <- out$engines[nzchar(out$engines)]
  out$max_cases <- suppressWarnings(as.integer(out$max_cases))
  out
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

upload_root <- file.path(app_dir, "examples", "00_UPLOAD_READY_public_integrated_cases")
if (!dir.exists(upload_root)) stop("Missing public integrated cases folder: ", upload_root, call. = FALSE)

case_names <- data.frame(
  case_id = c(
    "avonet_birds_range_grid_derived",
    "cestes_birds_barbaro2009b",
    "cestes_freshwater_macroinvertebrates_doledec1996",
    "cestes_plants_chmura2016",
    "cestes_aquatic_or_ground_beetles_ribera2001",
    "cestes_spiders_gibb2015",
    "masviken_mountain_plants",
    "masviken_mountain_spiders",
    "mediterranean_fish_efh_moullec",
    "phylacine_mammals_range_derived"
  ),
  short_id = c("pub01_avonet", "pub02_birds", "pub03_freshwater",
               "pub04_chmura", "pub05_ribera", "pub06_spiders",
               "pub07_masviken_plants", "pub08_masviken_spiders",
               "pub09_fish", "pub10_mammals"),
  display_name = c(
    "AVONET bird trait range-grid demonstration",
    "CESTES bird community example from Barbaro et al. 2009b",
    "CESTES freshwater macroinvertebrate community example from Doledec et al. 1996",
    "CESTES Chmura et al. 2016 community example",
    "CESTES Ribera et al. 2001 invertebrate community example",
    "CESTES spider community example from Gibb et al. 2015",
    "Masviken high-latitude mountain plant community example",
    "Masviken high-latitude mountain spider community example",
    "Mediterranean fish essential fish habitat community example",
    "PHYLACINE mammal range-derived trait and phylogeny demonstration"
  ),
  taxon_group = c("birds", "birds", "freshwater macroinvertebrates",
                  "CESTES community", "invertebrates", "spiders",
                  "plants", "spiders", "marine fish", "mammals"),
  source_kind = c("range-derived trait demonstration", "observed community",
                  "observed community", "observed community", "observed community",
                  "observed community", "observed community", "observed community",
                  "observed community", "range-derived demonstration"),
  recommended_use = c(
    "Software upload/trait/phylogeny demonstration; not a local field survey benchmark.",
    "Maximum-upload smoke test for multi-engine occurrence workflows.",
    "Maximum-upload smoke test for freshwater occurrence workflows.",
    "Larger CESTES occurrence workflow demonstration with all species retained after site selection.",
    "Larger CESTES invertebrate occurrence workflow demonstration with all species retained after site selection.",
    "Maximum-upload smoke test for arthropod occurrence workflows.",
    "Mountain plant occurrence workflow demonstration with derived upload covariates.",
    "Mountain spider occurrence workflow demonstration with derived upload covariates.",
    "Marine fish occurrence/count/trait workflow demonstration.",
    "Software upload/trait/phylogeny demonstration; not a local field survey benchmark."
  ),
  stringsAsFactors = FALSE
)

ensure_case_names <- function() {
  write.csv(case_names, file.path(upload_root, "PUBLIC_CASE_NAMES.csv"), row.names = FALSE)
  for (i in seq_len(nrow(case_names))) {
    case_id <- case_names$case_id[[i]]
    case_dir <- file.path(upload_root, case_id)
    if (!dir.exists(case_dir)) next
    for (subdir in c("00_common_inputs", sprintf("%02d_%s", 1:7,
                                                 c("Hmsc", "Hmsc-HPC", "jSDM", "GJAM",
                                                   "spOccupancy", "sjSDM", "boral")))) {
      d <- file.path(case_dir, subdir)
      if (!dir.exists(d)) next
      manifest <- file.path(d, "dataset_manifest.csv")
      if (file.exists(manifest)) {
        m <- tryCatch(read.csv(manifest, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(m) && all(c("field", "value") %in% names(m))) {
          m <- m[!m$field %in% c("display_name", "taxon_group", "source_kind", "recommended_use"), , drop = FALSE]
          m <- rbind(
            data.frame(field = "display_name", value = case_names$display_name[[i]], stringsAsFactors = FALSE),
            data.frame(field = "taxon_group", value = case_names$taxon_group[[i]], stringsAsFactors = FALSE),
            data.frame(field = "source_kind", value = case_names$source_kind[[i]], stringsAsFactors = FALSE),
            data.frame(field = "recommended_use", value = case_names$recommended_use[[i]], stringsAsFactors = FALSE),
            m
          )
          write.csv(m, manifest, row.names = FALSE)
        }
      }
      readme <- file.path(d, "README_UPLOAD_FIRST.md")
      if (file.exists(readme)) {
        body <- readLines(readme, warn = FALSE)
        body <- body[!grepl("^Display name:|^Taxon group:|^Source kind:|^Recommended use:", body)]
        header <- c(
          paste0("Display name: ", case_names$display_name[[i]]),
          paste0("Taxon group: ", case_names$taxon_group[[i]]),
          paste0("Source kind: ", case_names$source_kind[[i]]),
          paste0("Recommended use: ", case_names$recommended_use[[i]]),
          ""
        )
        writeLines(c(header, body), readme, useBytes = TRUE)
      }
    }
  }
  invisible(case_names)
}

read_csv_auto <- function(path, row.names = 1) {
  if (!file.exists(path)) stop("Missing upload file: ", path, call. = FALSE)
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(x) && !ncol(x)) return(x)
  if (!is.null(row.names) && ncol(x) >= row.names) {
    rn <- x[[row.names]]
    if (!anyDuplicated(rn) && all(nzchar(as.character(rn)))) {
      rownames(x) <- as.character(rn)
      x <- x[-row.names]
    }
  }
  x
}

read_text <- function(path) {
  if (!file.exists(path)) stop("Missing upload file: ", path, call. = FALSE)
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

numeric_df <- function(x) {
  rn <- rownames(x)
  cn <- colnames(x)
  y <- as.data.frame(lapply(x, function(z) suppressWarnings(as.numeric(z))),
                     check.names = FALSE, stringsAsFactors = FALSE)
  rownames(y) <- rn
  colnames(y) <- cn
  y
}

numeric_predictors <- function(X) {
  names(X)[vapply(X, function(z) {
    zz <- suppressWarnings(as.numeric(z))
    all(!is.na(zz)) && stats::var(zz, na.rm = TRUE) > 0
  }, logical(1))]
}

make_formula_dot <- function() "~ ."

make_stable_formula <- function(X, max_terms = 3L) {
  num <- numeric_predictors(X)
  if (!length(num)) return("~ 1")
  keep <- head(num, max_terms)
  paste("~", paste(sprintf("`%s`", keep), collapse = " + "))
}

make_patsy_stable_formula <- function(X, max_terms = 3L) {
  num <- numeric_predictors(X)
  num <- num[grepl("^[A-Za-z_][A-Za-z0-9_]*$", num)]
  if (!length(num)) return("~ 1")
  keep <- head(num, max_terms)
  paste("~", paste(keep, collapse = " + "))
}

make_reference_truth <- function(dat) {
  Y <- as.data.frame(dat$Y, check.names = FALSE)
  X <- as.data.frame(dat$X, check.names = FALSE)
  sites <- rownames(Y)
  species <- colnames(Y)
  num <- numeric_predictors(X)
  if (!length(num)) num <- names(X)[seq_len(min(2L, ncol(X)))]
  eff <- expand.grid(species = species, predictor = num, stringsAsFactors = FALSE)
  eff$true_effect <- NA_real_
  eff$reference_method <- NA_character_
  response_kind <- function(y) {
    y <- y[is.finite(y)]
    if (length(y) && all(y %in% c(0, 1))) return("binary")
    if (length(y) && all(y >= 0) && all(abs(y - round(y)) < 1e-8)) return("count")
    "continuous"
  }
  species_kind <- setNames(vapply(Y, function(z) response_kind(suppressWarnings(as.numeric(z))), character(1)), species)
  for (k in seq_len(nrow(eff))) {
    sp <- eff$species[[k]]
    pr <- eff$predictor[[k]]
    yy <- suppressWarnings(as.numeric(Y[[sp]]))
    xx <- suppressWarnings(as.numeric(X[[pr]]))
    kind <- species_kind[[sp]]
    eff$reference_method[[k]] <- paste0("public_case_", kind, "_univariate_reference")
    ok <- is.finite(yy) & is.finite(xx)
    if (sum(ok) >= 3 && length(unique(yy[ok])) >= 2 && stats::var(xx[ok], na.rm = TRUE) > 0) {
      fit <- tryCatch(suppressWarnings(switch(
        kind,
        binary = stats::glm(yy[ok] ~ xx[ok], family = stats::binomial()),
        count = stats::glm(yy[ok] ~ xx[ok], family = stats::poisson()),
        continuous = stats::lm(yy[ok] ~ xx[ok])
      )), error = function(e) NULL)
      if (!is.null(fit) && length(stats::coef(fit)) >= 2) eff$true_effect[[k]] <- unname(stats::coef(fit)[[2]])
      if (!is.finite(eff$true_effect[[k]])) {
        eff$true_effect[[k]] <- suppressWarnings(stats::cor(yy[ok], xx[ok], use = "pairwise.complete.obs"))
      }
    }
  }
  pred <- as.data.frame(as.table(as.matrix(Y)), stringsAsFactors = FALSE)
  names(pred) <- c("site_id", "species", "observed")
  species_mean <- colMeans(as.matrix(Y), na.rm = TRUE)
  pred$truth_probability <- as.numeric(species_mean[as.character(pred$species)])
  pred$reference_method <- ifelse(species_kind[as.character(pred$species)] == "binary",
                                  "species_prevalence_reference_for_public_case",
                                  "species_mean_response_reference_for_public_case")
  assoc <- utils::combn(species, 2, simplify = FALSE)
  assoc <- do.call(rbind, lapply(assoc, function(p) {
    data.frame(species_i = p[[1]], species_j = p[[2]],
               true_association = suppressWarnings(stats::cor(as.numeric(Y[[p[[1]]]]),
                                                              as.numeric(Y[[p[[2]]]]),
                                                              use = "pairwise.complete.obs")),
               reference_method = "observed_response_correlation_reference",
               stringsAsFactors = FALSE)
  }))
  dat$truth_effects <- eff
  dat$truth_predictions <- pred[, c("site_id", "species", "truth_probability", "reference_method")]
  dat$truth_associations <- assoc
  dat
}

load_public_dat <- function(case_id, engine_folder = "00_common_inputs") {
  case_dir <- file.path(upload_root, case_id)
  src <- file.path(case_dir, engine_folder)
  if (!dir.exists(src)) src <- file.path(case_dir, "00_common_inputs")
  if (!dir.exists(src)) stop("Missing case input folder: ", src, call. = FALSE)
  dat <- list()
  dat$Y <- numeric_df(read_csv_auto(file.path(src, "Y_occurrence.csv")))
  dat$Y_count <- numeric_df(read_csv_auto(file.path(src, "Y_count.csv")))
  dat$Y_normal <- numeric_df(read_csv_auto(file.path(src, "Y_normal.csv")))
  dat$X <- read_csv_auto(file.path(src, "XData.csv"))
  dat$traits <- read_csv_auto(file.path(src, "traits.csv"))
  dat$study <- read_csv_auto(file.path(src, "studyDesign.csv"))
  dat$coords <- read_csv_auto(file.path(src, "coordinates.csv"))
  taxonomy_path <- file.path(src, "taxonomy.csv")
  dat$taxonomy <- if (file.exists(taxonomy_path)) {
    read_csv_auto(taxonomy_path, row.names = NULL)
  } else {
    data.frame(taxon_id = colnames(dat$Y),
               kingdom = "Unknown", phylum = "Unknown", class = "Unknown",
               order = paste0("Order_", rep(seq_len(3), length.out = ncol(dat$Y))),
               family = paste0("Family_", rep(seq_len(4), length.out = ncol(dat$Y))),
               genus = paste0("Genus_", seq_len(ncol(dat$Y))),
               species = colnames(dat$Y), stringsAsFactors = FALSE)
  }
  dat$newick <- read_text(file.path(src, "phylogeny.nwk"))
  dat$newdata <- read_csv_auto(file.path(src, "newdata.csv"))
  dat$newcoords <- read_csv_auto(file.path(src, "newcoords.csv"))
  dat$folds <- read_csv_auto(file.path(src, "folds.csv"), row.names = NULL)
  dat$trial <- numeric_df(read_csv_auto(file.path(src, "trial.size.csv")))
  dat$offset <- numeric_df(read_csv_auto(file.path(src, "offset.csv")))
  dat$rowids <- read_csv_auto(file.path(src, "row.ids.csv"))
  dat$ranefids <- read_csv_auto(file.path(src, "ranef.ids.csv"))
  dat$distmat <- numeric_df(read_csv_auto(file.path(src, "distmat.csv")))
  dat$det <- numeric_df(read_csv_auto(file.path(src, "spOccupancy_y_detection.csv")))
  dat$occ_covs <- read_csv_auto(file.path(src, "occ.covs.csv"))
  dat$det_cov <- read_csv_auto(file.path(src, "det.covs.csv"))
  dat$upload_engine_folder <- src
  dat
}

select_fit_species <- function(Y, max_species = 30L) {
  species <- colnames(Y)
  if (length(species) <= max_species) return(species)
  yy <- as.matrix(Y)
  prev <- colMeans(yy, na.rm = TRUE)
  ok <- is.finite(prev) & prev > 0 & prev < 1
  candidates <- species[ok]
  if (length(candidates) < min(max_species, length(species))) {
    candidates <- species[is.finite(prev) & prev > 0]
  }
  if (length(candidates) < min(max_species, length(species))) {
    candidates <- species
  }
  score <- abs(prev[candidates] - 0.5)
  candidates[order(score, candidates)][seq_len(min(max_species, length(candidates)))]
}

subset_public_dat_for_fit <- function(dat, engine, max_species = 30L) {
  original_species <- colnames(dat$Y)
  keep <- select_fit_species(dat$Y, max_species = max_species)
  subset_needed <- length(original_species) > length(keep)
  subset_cols <- function(x) {
    hit <- intersect(keep, colnames(x))
    x[, hit, drop = FALSE]
  }
  dat$Y <- subset_cols(dat$Y)
  dat$Y_count <- subset_cols(dat$Y_count)
  dat$Y_normal <- subset_cols(dat$Y_normal)
  dat$trial <- subset_cols(dat$trial)
  dat$offset <- subset_cols(dat$offset)
  if (!is.null(dat$traits) && nrow(dat$traits)) {
    dat$traits <- dat$traits[intersect(keep, rownames(dat$traits)), , drop = FALSE]
  }
  if (!is.null(dat$taxonomy) && nrow(dat$taxonomy)) {
    nm_low <- tolower(names(dat$taxonomy))
    id_col <- names(dat$taxonomy)[nm_low %in% c("species", "species_name", "species_id", "taxon", "taxon_id", "taxon_name", "response", "response_id")]
    if (length(id_col)) dat$taxonomy <- dat$taxonomy[as.character(dat$taxonomy[[id_col[1]]]) %in% keep, , drop = FALSE]
  }
  if (!is.null(dat$det) && ncol(dat$det)) {
    keep_prefix <- paste0(keep, "_rep")
    det_keep <- vapply(colnames(dat$det), function(z) any(startsWith(z, keep_prefix)), logical(1))
    dat$det <- dat$det[, det_keep, drop = FALSE]
  }
  dat$fit_species_subset_info <- data.frame(
    engine = engine,
    uploaded_species_n = length(original_species),
    fitted_species_n = length(keep),
    subset_applied = subset_needed,
    selection_rule = if (subset_needed) {
      paste0("All upload files are copied unchanged. Quick engine fitting uses up to ",
             max_species, " species closest to prevalence 0.5 to keep public maximum-upload tests stable.")
    } else {
      "All uploaded species were used for fitting."
    },
    stringsAsFactors = FALSE
  )
  dat$fit_species_subset_species <- data.frame(
    engine = engine,
    species = keep,
    prevalence = as.numeric(colMeans(as.matrix(dat$Y), na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  dat <- make_reference_truth(dat)
  dat
}

write_public_truth_files <- function(root, case_id, dat) {
  truth_dir <- file.path(root, "00_case_truth_reference", case_id)
  dir.create(truth_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(dat$truth_effects, file.path(truth_dir, "reference_effects_species_environment.csv"), row.names = FALSE)
  write.csv(dat$truth_predictions, file.path(truth_dir, "reference_predictions_site_species.csv"), row.names = FALSE)
  write.csv(dat$truth_associations, file.path(truth_dir, "reference_associations_species_species.csv"), row.names = FALSE)
  writeLines(c(
    "These are public-case reference tables, not known latent ecological truth.",
    "Effect references are simple GLM/correlation summaries from the uploaded occurrence matrix and predictors.",
    "Prediction references are species prevalence values for smoke-test comparison."
  ), file.path(truth_dir, "README_reference_truth.txt"), useBytes = TRUE)
}

copy_max_upload_files <- function(outdir, dat) {
  src <- dat$upload_engine_folder
  files <- list.files(src, full.names = TRUE, recursive = FALSE)
  files <- files[!tolower(basename(files)) %in% c("phylo_cov.csv", "c_phylo.csv")]
  rows <- lapply(files, function(f) {
    rel <- basename(f)
    for (root_name in c("inputs", "data")) {
      dest <- file.path(outdir, root_name, "maximum_upload_files", rel)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      file.copy(f, dest, overwrite = TRUE)
    }
    data.frame(file = rel, size_bytes = file.info(f)$size, copied = TRUE, stringsAsFactors = FALSE)
  })
  tab <- if (length(rows)) do.call(rbind, rows) else data.frame(file = character(), size_bytes = numeric(), copied = logical())
  write.csv(tab, file.path(outdir, "tables", "maximum_upload_file_check.csv"), row.names = FALSE)
  writeLines(c("Maximum-upload files copied from:", normalizePath(src, winslash = "/", mustWork = FALSE),
               "These files are copied into the ZIP for auditability. The engine fit uses the compatible subset required by that package."),
             file.path(outdir, "inputs", "README_maximum_upload.txt"), useBytes = TRUE)
  if (!is.null(dat$fit_species_subset_info)) {
    write.csv(dat$fit_species_subset_info, file.path(outdir, "tables", "modeling_species_subset_summary.csv"), row.names = FALSE)
    write.csv(dat$fit_species_subset_species, file.path(outdir, "tables", "modeling_species_subset_species.csv"), row.names = FALSE)
    writeLines(c(
      "Maximum-upload public case handling:",
      paste0("Uploaded species: ", dat$fit_species_subset_info$uploaded_species_n[[1]]),
      paste0("Fitted species in this quick branch: ", dat$fit_species_subset_info$fitted_species_n[[1]]),
      paste0("Subset applied: ", dat$fit_species_subset_info$subset_applied[[1]]),
      dat$fit_species_subset_info$selection_rule[[1]]
    ), file.path(outdir, "diagnostics", "modeling_species_subset_note.txt"), useBytes = TRUE)
  }
  if (!is.null(dat$gjam_name_mapping)) {
    write.csv(dat$gjam_name_mapping, file.path(outdir, "tables", "GJAM_name_mapping.csv"), row.names = FALSE)
    write.csv(dat$gjam_name_mapping, file.path(outdir, "inputs", "GJAM_name_mapping.csv"), row.names = FALSE)
  }
  invisible(nrow(tab))
}

refresh_result_row <- function(res, dat, case_display_name = "", source_case_id = "") {
  outdir <- res$output_folder[1]
  copied <- copy_max_upload_files(outdir, dat)
  runner_env$fill_empty_dirs(outdir, res$engine[1])
  zipfile <- tryCatch(runner_env$make_zip_file(outdir), error = function(e) {
    writeLines(conditionMessage(e), file.path(outdir, "diagnostics", "zip_refresh_error.txt"))
    NA_character_
  })
  missing_dirs <- senv$contract_dirs[!dir.exists(file.path(outdir, senv$contract_dirs))]
  missing_files <- senv$contract_files[!file.exists(file.path(outdir, senv$contract_files))]
  res$case_display_name <- case_display_name
  res$source_case_id <- source_case_id
  res$zip <- if (!is.na(zipfile)) normalizePath(zipfile, winslash = "/", mustWork = FALSE) else NA_character_
  res$zip_size_kb <- if (!is.na(zipfile) && file.exists(zipfile)) round(file.info(zipfile)$size / 1024, 1) else NA_real_
  res$file_count <- length(list.files(outdir, recursive = TRUE, all.files = FALSE))
  res$missing_count <- length(c(missing_dirs, missing_files))
  res$missing_items <- paste(c(missing_dirs, missing_files), collapse = "; ")
  res$maximum_upload_files_copied <- copied
  res
}

case_meta <- function(case_id) {
  hit <- case_names[case_names$case_id == case_id, , drop = FALSE]
  if (nrow(hit)) hit[1, , drop = FALSE] else data.frame(case_id = case_id, display_name = case_id,
                                                        taxon_group = "", source_kind = "",
                                                        recommended_use = "", stringsAsFactors = FALSE)
}

short_case_slug <- function(case_id) {
  hit <- case_names$short_id[match(case_id, case_names$case_id)]
  ifelse(is.na(hit) | !nzchar(hit), gsub("[^A-Za-z0-9]+", "_", case_id), hit)
}

safe_add_result <- function(results, expr, dat, case_id, engine) {
  idx <- length(results) + 1L
  meta <- case_meta(case_id)
  results[[idx]] <- tryCatch(refresh_result_row(force(expr), dat, meta$display_name[[1]], case_id), error = function(e) {
    data.frame(engine = engine, case_id = case_id, status = "fit_failed",
               output_folder = NA_character_, zip = NA_character_, zip_size_kb = NA_real_,
               file_count = NA_integer_, missing_count = NA_integer_, missing_items = "",
               failure_reason = conditionMessage(e), warnings = "",
               notes = "Unhandled public maximum-upload suite-level error.",
               runtime_seconds = NA_real_, source_case_id = case_id,
               case_display_name = meta$display_name[[1]],
               maximum_upload_files_copied = NA_integer_, stringsAsFactors = FALSE)
  })
  results
}

adapt_hmsc_cfg <- function(dat) {
  cfg <- runner_env$hmsc_cfg()
  cfg$model$XFormula <- make_stable_formula(dat$X, max_terms = 1L)
  cfg$model$TrFormula <- make_formula_dot()
  cfg$model$random_mode <- "none"
  cfg$model$random_effect_column <- "sample"
  cfg$model$use_traits <- FALSE
  cfg$model$use_phylogeny <- FALSE
  cfg$model$distr <- "probit"
  cfg$model$C_file <- ""
  cfg$model$nfMax <- min(4L, max(1L, ncol(dat$Y) - 1L))
  cfg$mcmc$samples <- 3L
  cfg$mcmc$transient <- 3L
  cfg$mcmc$thin <- 1L
  cfg$mcmc$nChains <- 1L
  cfg$mcmc$nParallel <- 1L
  cfg$outputs$cv <- FALSE
  cfg$outputs$waic <- FALSE
  cfg$outputs$gradients <- FALSE
  cfg$outputs$predictions$nfolds <- 2L
  cfg$outputs$predictions$partition_column <- "sample"
  cfg
}

adapt_hmschpc_cfg <- function(dat) {
  cfg <- runner_env$hmschpc_cfg(app_env)
  cfg$model$XFormula <- make_patsy_stable_formula(dat$X, max_terms = 1L)
  cfg$model$trait_formula <- make_formula_dot()
  cfg$model$distribution <- "probit"
  cfg$model$use_traits <- FALSE
  cfg$model$phylogeny_mode <- "none"
  cfg$random_effects$mode <- "none"
  cfg$random_effects$column <- "sample"
  cfg$random_effects$name <- "sample_iid"
  cfg$random_effects$nf <- 1L
  cfg$random_effects$nfMin <- 1L
  cfg$random_effects$nfMax <- min(4L, max(1L, ncol(dat$Y) - 1L))
  cfg$sampler$samples <- 2L
  cfg$sampler$transient <- 2L
  cfg$sampler$chains <- 1L
  cfg$sampler$verbose <- 1L
  cfg$runtime$python <- app_env$hmschpc_python_bin(list())
  cfg$runtime$python_source <- app_env$hmschpc_python_source_dir()
  cfg
}

adapt_jsdm_cfg <- function(dat) {
  cfg <- runner_env$jsdm_cfg()
  cfg$model$site_formula <- make_formula_dot()
  cfg$model$trait_formula <- make_formula_dot()
  cfg$model$model_type <- "binomial_probit"
  cfg$model$n_latent <- min(2L, max(0L, ncol(dat$Y) - 1L))
  cfg$model$site_effect <- "none"
  cfg$model$allow_traits <- FALSE
  cfg$mcmc$burnin <- 50L
  cfg$mcmc$mcmc <- 50L
  cfg$mcmc$thin <- 1L
  cfg$prediction$do_predict <- TRUE
  cfg
}

adapt_gjam_cfg <- function(dat) {
  cfg <- runner_env$gjam_cfg(colnames(dat$Y))
  first_numeric <- numeric_predictors(dat$X)[1] %||% names(dat$X)[1]
  cfg$model$formula <- if (!is.na(first_numeric) && nzchar(first_numeric)) paste("~", first_numeric) else "~ 1"
  cfg$model$type_single <- "PA"
  cfg$model$typeNames_text <- paste(rep("PA", ncol(dat$Y)), collapse = ",")
  cfg$model$ng <- 50L
  cfg$model$burnin <- 10L
  cfg$analysis$do_traits <- FALSE
  cfg$analysis$do_sensitivity <- FALSE
  cfg$analysis$do_predict <- TRUE
  cfg
}

gjam_trait_subset <- function(dat) {
  tr <- as.data.frame(dat$traits, check.names = FALSE)
  keep <- seq_len(min(3L, ncol(tr)))
  tr <- tr[, keep, drop = FALSE]
  while (ncol(tr) < 3L) {
    tr[[paste0("placeholder_trait_", ncol(tr) + 1L)]] <- seq_len(nrow(tr))
  }
  rownames(tr) <- rownames(dat$traits)
  dat$traits <- tr
  dat
}

gjam_safe_names <- function(dat) {
  make_safe <- function(prefix, n) paste0(prefix, seq_len(n))
  old_species <- colnames(dat$Y)
  safe_species <- make_safe("gsp", length(old_species))
  old_predictors <- colnames(dat$X)
  safe_predictors <- make_safe("gx", length(old_predictors))
  colnames(dat$Y) <- safe_species
  colnames(dat$Y_count) <- safe_species
  colnames(dat$Y_normal) <- safe_species
  colnames(dat$trial) <- safe_species
  colnames(dat$offset) <- safe_species
  if (!is.null(dat$taxonomy) && nrow(dat$taxonomy)) {
    nm_low <- tolower(names(dat$taxonomy))
    id_col <- names(dat$taxonomy)[nm_low %in% c("species", "species_name", "species_id", "taxon", "taxon_id", "taxon_name", "response", "response_id")]
    if (length(id_col)) dat$taxonomy[[id_col[1]]] <- safe_species[match(dat$taxonomy[[id_col[1]]], old_species)]
    if ("species" %in% names(dat$taxonomy)) dat$taxonomy$species <- safe_species[match(dat$taxonomy$species, old_species)]
  }
  rownames(dat$traits) <- safe_species
  colnames(dat$X) <- safe_predictors
  colnames(dat$newdata) <- safe_predictors
  colnames(dat$occ_covs) <- safe_predictors
  dat$truth_effects$species <- safe_species[match(dat$truth_effects$species, old_species)]
  dat$truth_effects$predictor <- safe_predictors[match(dat$truth_effects$predictor, old_predictors)]
  dat$truth_predictions$species <- safe_species[match(dat$truth_predictions$species, old_species)]
  dat$truth_associations$species_i <- safe_species[match(dat$truth_associations$species_i, old_species)]
  dat$truth_associations$species_j <- safe_species[match(dat$truth_associations$species_j, old_species)]
  dat$gjam_name_mapping <- rbind(
    data.frame(kind = "species", original = old_species, safe = safe_species, stringsAsFactors = FALSE),
    data.frame(kind = "predictor", original = old_predictors, safe = safe_predictors, stringsAsFactors = FALSE)
  )
  dat
}

adapt_spocc_cfg <- function(dat) {
  cfg <- runner_env$spocc_cfg()
  cfg$model$model_type <- "msPGOcc"
  cfg$model$data_structure <- "multi-species replicated"
  cfg$model$occ.formula <- make_formula_dot()
  cfg$model$det.formula <- "~ obs"
  cfg$model$formula <- make_formula_dot()
  cfg$mcmc$n.batch <- 3L
  cfg$mcmc$batch.length <- 3L
  cfg$mcmc$n.burn <- 3L
  cfg$mcmc$n.thin <- 1L
  cfg$mcmc$n.chains <- 1L
  cfg$mcmc$n.report <- 3L
  cfg$spatial_latent_svc$n.factors <- min(2L, max(1L, ncol(dat$Y) - 1L))
  cfg$spatial_latent_svc$NNGP <- TRUE
  cfg$spatial_latent_svc$n.neighbors <- 3L
  cfg$validation_prediction_outputs$waicOcc <- FALSE
  cfg$validation_prediction_outputs$save_plots <- FALSE
  cfg$validation_prediction_outputs$predict <- FALSE
  cfg
}

adapt_sjsdm_cfg <- function(dat) {
  cfg <- runner_env$sjsdm_cfg()
  cfg$model$family <- "binomial_probit"
  cfg$model$env_model <- "linear"
  cfg$model$env_formula <- make_formula_dot()
  cfg$model$spatial_model <- "linear"
  cfg$model$spatial_formula <- "~ 0 + ."
  cfg$model$iter <- 3L
  cfg$model$sampling <- 50L
  cfg$model$step_size <- min(8L, nrow(dat$Y))
  cfg$model$parallel <- 0L
  cfg$dnn_optimizer$device <- "cpu"
  cfg$outputs$plots <- FALSE
  cfg
}

adapt_boral_cfg <- function(dat) {
  cfg <- runner_env$boral_cfg()
  cfg$model$family <- "binomial"
  cfg$model$formula.X <- make_formula_dot()
  cfg$model$trial.size <- 1L
  cfg$model$lv.control$num.lv <- min(1L, max(0L, ncol(dat$Y) - 1L))
  cfg$model$row.eff <- "none"
  cfg$traits_random_ssvs$use_traits <- FALSE
  cfg$mcmc_prior$n.burnin <- 10L
  cfg$mcmc_prior$n.iteration <- 30L
  cfg$mcmc_prior$n.thin <- 1L
  cfg$mcmc_prior$calc.ics <- FALSE
  cfg
}

run_one_case <- function(root, case_id, engines) {
  meta <- case_meta(case_id)
  slug <- short_case_slug(case_id)
  results <- list()
  case_dir <- file.path(upload_root, case_id)
  if (!dir.exists(case_dir)) stop("Missing public case: ", case_id, call. = FALSE)

  if ("Hmsc" %in% engines) {
    dat <- load_public_dat(case_id, "01_Hmsc")
    dat <- subset_public_dat_for_fit(dat, "Hmsc")
    write_public_truth_files(root, case_id, dat)
    cfg <- adapt_hmsc_cfg(dat)
    results <- safe_add_result(results,
      senv$run_hmsc_case(root, paste0(slug, "__hmsc_fixed"),
                         dat, cfg, dat$Y,
                         paste0(meta$display_name[[1]], "; maximum-upload quick branch: stable fixed-effect probit with one numeric environmental predictor. Full XData, random-effect, traits and phylogeny files are uploaded/copied but not all enabled in this all-case smoke branch.")),
      dat, case_id, "Hmsc")
  }

  if ("Hmsc-HPC" %in% engines) {
    dat <- load_public_dat(case_id, "02_Hmsc-HPC")
    dat <- subset_public_dat_for_fit(dat, "Hmsc-HPC")
    cfg <- adapt_hmschpc_cfg(dat)
    results <- safe_add_result(results,
      senv$run_hmschpc_case(root, paste0(slug, "__hpc_fixed"),
                            dat, cfg, dat$Y,
                            paste0(meta$display_name[[1]], "; maximum-upload quick branch: CPU probit fixed-effect fit. Traits, phylogeny and random-effect files are uploaded/copied but not enabled in this stable smoke branch.")),
      dat, case_id, "Hmsc-HPC")
  }

  if ("jSDM" %in% engines) {
    dat <- load_public_dat(case_id, "03_jSDM")
    dat <- subset_public_dat_for_fit(dat, "jSDM")
    cfg <- adapt_jsdm_cfg(dat)
    results <- safe_add_result(results,
      senv$run_jsdm_case(root, paste0(slug, "__jsdm_probit"),
                         dat, cfg, dat$Y,
                         trials = dat$trial,
                         notes = paste0(meta$display_name[[1]], "; maximum-upload quick branch: binomial probit matrix input.")),
      dat, case_id, "jSDM")
  }

  if ("GJAM" %in% engines) {
    dat <- load_public_dat(case_id, "04_GJAM")
    dat <- subset_public_dat_for_fit(dat, "GJAM")
    dat <- gjam_safe_names(dat)
    dat <- gjam_trait_subset(dat)
    cfg <- adapt_gjam_cfg(dat)
    results <- safe_add_result(results,
      senv$run_gjam_case(root, paste0(slug, "__gjam_pa"),
                         dat, cfg, dat$Y, rep("PA", ncol(dat$Y)),
                         paste0(meta$display_name[[1]], "; maximum-upload quick branch: PA typeNames.")),
      dat, case_id, "GJAM")
  }

  if ("spOccupancy" %in% engines) {
    dat <- load_public_dat(case_id, "05_spOccupancy")
    dat <- subset_public_dat_for_fit(dat, "spOccupancy")
    cfg <- adapt_spocc_cfg(dat)
    sp_tab <- data.frame(species = colnames(dat$Y), stringsAsFactors = FALSE)
    results <- safe_add_result(results,
      senv$run_spocc_case(root, paste0(slug, "__spocc_mspgocc"),
                          dat, cfg, dat$det, species_table = sp_tab,
                          notes = paste0(meta$display_name[[1]], "; maximum-upload quick branch: replicated multi-species detection/occupancy.")),
      dat, case_id, "spOccupancy")
  }

  if ("sjSDM" %in% engines) {
    dat <- load_public_dat(case_id, "06_sjSDM")
    dat <- subset_public_dat_for_fit(dat, "sjSDM")
    cfg <- adapt_sjsdm_cfg(dat)
    results <- safe_add_result(results,
      senv$run_sjsdm_case(root, paste0(slug, "__sjsdm_linear"),
                          dat, cfg, dat$Y,
                          paste0(meta$display_name[[1]], "; maximum-upload quick branch: CPU binomial linear env/spatial model.")),
      dat, case_id, "sjSDM")
  }

  if ("boral" %in% engines) {
    dat <- load_public_dat(case_id, "07_boral")
    dat <- subset_public_dat_for_fit(dat, "boral")
    cfg <- adapt_boral_cfg(dat)
    results <- safe_add_result(results,
      senv$run_boral_case(root, paste0(slug, "__boral_binomial"),
                          dat, cfg, dat$Y,
                          paste0(meta$display_name[[1]], "; maximum-upload quick branch: binomial/trial.size=1 latent model.")),
      dat, case_id, "boral")
  }

  do.call(rbind, results)
}

write_report <- function(root, summary) {
  by_engine <- aggregate(case_id ~ engine + status, data = summary, FUN = length)
  names(by_engine)[names(by_engine) == "case_id"] <- "n_cases"
  by_case <- aggregate(engine ~ case_id + case_display_name + status, data = summary, FUN = length)
  names(by_case)[names(by_case) == "engine"] <- "n_engines"
  write.csv(by_engine, file.path(root, "public_maximum_upload_status_by_engine.csv"), row.names = FALSE)
  write.csv(by_case, file.path(root, "public_maximum_upload_status_by_case.csv"), row.names = FALSE)

  lines <- c(
    "# Public Integrated Case Maximum-Upload Audit",
    "",
    paste0("Generated: ", Sys.time()),
    paste0("Audit root: ", normalizePath(root, winslash = "/", mustWork = FALSE)),
    "",
    "## What Was Tested",
    "",
    "Each listed public case was given a display name and run through one stable quick branch per selected engine.",
    "All files in that engine's upload-ready folder were copied into inputs/maximum_upload_files and data/maximum_upload_files.",
    "The fitted branch used the engine-compatible subset of those files. This avoids forcing incompatible CSVs into engines that cannot use them.",
    "",
    "## Important Interpretation",
    "",
    "For public cases, the reference truth tables are derived from observed occurrence/covariate summaries, not known latent ecological truth.",
    "A status of fitted means the adapter wrote fitted outputs. check_failed and fit_failed are honest failures and should be inspected through diagnostics/.",
    "",
    "## Status By Engine",
    "",
    paste(capture.output(print(by_engine, row.names = FALSE)), collapse = "\n"),
    "",
    "## Case Summary",
    "",
    paste(capture.output(print(summary[, c("case_id", "case_display_name", "engine", "status", "file_count", "missing_count", "zip_size_kb", "maximum_upload_files_copied", "failure_reason")],
                               row.names = FALSE)), collapse = "\n")
  )
  writeLines(lines, file.path(root, "public_maximum_upload_audit_report.md"), useBytes = TRUE)
}

run_suite <- function() {
  cli <- parse_cli()
  ensure_case_names()
  all_cases <- list.dirs(upload_root, full.names = FALSE, recursive = FALSE)
  all_cases <- all_cases[all_cases %in% case_names$case_id]
  cases <- if (identical(tolower(cli$cases), "all")) all_cases else trimws(unlist(strsplit(cli$cases, ",")))
  cases <- cases[nzchar(cases) & cases %in% all_cases]
  if (length(cli$max_cases) == 1L && !is.na(cli$max_cases) && cli$max_cases > 0L) {
    cases <- head(cases, cli$max_cases)
  }
  if (!length(cases)) stop("No public cases selected.", call. = FALSE)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  root <- file.path(app_dir, "output", paste0("public_maximum_upload_audit_", timestamp))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  write.csv(case_names, file.path(root, "PUBLIC_CASE_NAMES.csv"), row.names = FALSE)
  results <- list()
  for (case_id in cases) {
    message("Running public maximum-upload case: ", case_id)
    results[[case_id]] <- run_one_case(root, case_id, cli$engines)
  }
  summary <- do.call(rbind, results)
  rownames(summary) <- NULL
  write.csv(summary, file.path(root, "public_maximum_upload_case_summary.csv"), row.names = FALSE)
  write_report(root, summary)
  zipfile <- runner_env$make_zip_file(root, file.path(app_dir, "output", paste0("JSDMStudio_public_maximum_upload_audit_", timestamp, ".zip")))
  sha <- if (file.exists(zipfile)) digest::digest(zipfile, algo = "sha256", file = TRUE) else NA_character_
  writeLines(c(
    paste0("AUDIT_ROOT=", normalizePath(root, winslash = "/", mustWork = FALSE)),
    paste0("SUMMARY_CSV=", normalizePath(file.path(root, "public_maximum_upload_case_summary.csv"), winslash = "/", mustWork = FALSE)),
    paste0("AUDIT_ZIP=", normalizePath(zipfile, winslash = "/", mustWork = FALSE)),
    paste0("AUDIT_ZIP_SHA256=", sha)
  ), file.path(root, "public_maximum_upload_audit_paths.txt"), useBytes = TRUE)
  cat("AUDIT_ROOT=", normalizePath(root, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("SUMMARY_CSV=", normalizePath(file.path(root, "public_maximum_upload_case_summary.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("AUDIT_ZIP=", normalizePath(zipfile, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("AUDIT_ZIP_SHA256=", sha, "\n", sep = "")
  attr(summary, "audit_root") <- normalizePath(root, winslash = "/", mustWork = FALSE)
  invisible(summary)
}

if (!interactive()) {
  script_warnings <- character()
  suite_summary <- withCallingHandlers(
    run_suite(),
    warning = function(w) {
      script_warnings <<- c(script_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  audit_root <- attr(suite_summary, "audit_root") %||% file.path(app_dir, "output")
  if (length(script_warnings)) {
    writeLines(unique(script_warnings), file.path(audit_root, "public_maximum_upload_script_warnings.txt"), useBytes = TRUE)
  }
}
