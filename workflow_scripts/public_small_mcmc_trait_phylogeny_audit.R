#!/usr/bin/env Rscript

# Small but real MCMC audit for public traits + phylogeny branches.
# This is a sanity MCMC run, not a publication-level convergence run.

if (.Platform$OS.type == "windows") {
  Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else getwd()
app_dir <- normalizePath(file.path(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)), ".."),
                         winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
setwd(app_dir)

suite_file <- file.path(app_dir, "workflow_scripts", "public_maximum_upload_case_suite.R")
suite_lines <- readLines(suite_file, warn = FALSE)
suite_lines <- suite_lines[!grepl("^\\s*if\\s*\\(!interactive\\(\\)\\)\\s*run_suite\\(\\)\\s*$", suite_lines)]
penv <- new.env(parent = globalenv())
eval(parse(text = paste(suite_lines, collapse = "\n")), envir = penv)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
root <- file.path(app_dir, "output", paste0("public_small_mcmc_trait_phylogeny_audit_", timestamp))
dir.create(root, recursive = TRUE, showWarnings = FALSE)

stable_trait_formula <- function(Tr, max_terms = 2L) {
  if (is.null(Tr) || ncol(Tr) == 0L) return("~ 1")
  tr <- as.data.frame(Tr, check.names = FALSE)
  num <- names(tr)[vapply(tr, is.numeric, logical(1))]
  if (length(num) == 0L) return("~ 1")
  paste("~", paste(utils::head(num, max_terms), collapse = " + "))
}

repair_phylo_cov <- function(C, species) {
  info <- list(repaired = FALSE, min_eigen_before = NA_real_, min_eigen_after = NA_real_, note = "not_checked")
  C <- as.matrix(C)
  storage.mode(C) <- "numeric"
  if (!all(species %in% rownames(C)) || !all(species %in% colnames(C))) {
    stop("phylo_cov names do not match fitted species.")
  }
  C <- C[species, species, drop = FALSE]
  C[!is.finite(C)] <- 0
  C <- (C + t(C)) / 2
  diag(C) <- 1
  eig <- tryCatch(eigen(C, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
  info$min_eigen_before <- suppressWarnings(min(eig, na.rm = TRUE))
  if (!is.finite(info$min_eigen_before) || info$min_eigen_before <= 1e-8) {
    jitter <- if (is.finite(info$min_eigen_before)) abs(info$min_eigen_before) + 1e-6 else 1e-6
    C <- C + diag(jitter, nrow(C))
    d <- sqrt(diag(C))
    C <- C / (d %o% d)
    C <- (C + t(C)) / 2
    diag(C) <- 1
    info$repaired <- TRUE
    info$note <- paste0("positive-definite jitter applied: ", signif(jitter, 4))
  } else {
    info$note <- "positive_definite"
  }
  eig2 <- tryCatch(eigen(C, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
  info$min_eigen_after <- suppressWarnings(min(eig2, na.rm = TRUE))
  rownames(C) <- colnames(C) <- species
  list(C = C, info = info)
}

finalize <- function(res, dat, case_id, engine, trait_formula, phylo_info, mcmc_label) {
  outdir <- res$output_folder[[1]]
  res <- penv$refresh_result_row(res, dat, penv$case_meta(case_id)$display_name[[1]], case_id)
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
  tab <- data.frame(
    engine = engine,
    source_case_id = case_id,
    mcmc_label = mcmc_label,
    fitted_species_n = ncol(dat$Y),
    trait_rows = nrow(dat$traits),
    trait_cols = ncol(dat$traits),
    trait_formula = trait_formula,
    phylogeny_input = "taxonomy.csv",
    derived_matrix_expected = TRUE,
    phylo_min_eigen_before = NA_real_,
    phylo_min_eigen_after = NA_real_,
    phylo_note = phylo_info$note,
    stringsAsFactors = FALSE
  )
  write.csv(tab, file.path(outdir, "tables", "small_mcmc_trait_phylogeny_use_summary.csv"), row.names = FALSE)
  files <- list.files(outdir, recursive = TRUE)
  res$has_engine_status <- file.exists(file.path(outdir, "diagnostics", "engine_status.json"))
  res$has_gamma_output <- any(grepl("Gamma|trait", files, ignore.case = TRUE))
  res$has_phylogeny_summary <- file.exists(file.path(outdir, "tables", "S1_phylogeny_summary.csv")) ||
    any(grepl("phylo", files, ignore.case = TRUE))
  res$has_standard_effects <- file.exists(file.path(outdir, "standard", "effects_species_environment.csv"))
  res$has_standard_predictions <- file.exists(file.path(outdir, "standard", "predictions_site_species.csv"))
  zipfile <- tryCatch(penv$runner_env$make_zip_file(outdir), error = function(e) NA_character_)
  if (!is.na(zipfile)) {
    res$zip <- normalizePath(zipfile, winslash = "/", mustWork = FALSE)
    res$zip_size_kb <- round(file.info(zipfile)$size / 1024, 1)
  }
  res$file_count <- length(list.files(outdir, recursive = TRUE, all.files = FALSE))
  res
}

run_hmsc <- function(case_id) {
  dat <- penv$load_public_dat(case_id, "01_Hmsc")
  dat <- penv$subset_public_dat_for_fit(dat, "Hmsc_small_mcmc_traits_phylogeny", max_species = 12L)
  phylo_info <- list(repaired = FALSE, min_eigen_before = NA_real_, min_eigen_after = NA_real_,
                     note = "taxonomy rank table used; derived matrix written by JSDM Studio")
  cfg <- penv$adapt_hmsc_cfg(dat)
  cfg$model$use_traits <- TRUE
  cfg$model$use_phylogeny <- TRUE
  cfg$model$TrFormula <- stable_trait_formula(dat$traits, max_terms = 2L)
  cfg$data$phylogeny <- "taxonomy.csv"
  cfg$model$C_file <- ""
  cfg$model$random_mode <- "none"
  cfg$model$nfMax <- min(4L, max(1L, ncol(dat$Y) - 1L))
  cfg$mcmc$samples <- 50L
  cfg$mcmc$transient <- 50L
  cfg$mcmc$thin <- 1L
  cfg$mcmc$nChains <- 2L
  cfg$mcmc$nParallel <- 1L
  cfg$mcmc$verbose <- 10L
  cfg$outputs$cv <- FALSE
  cfg$outputs$waic <- FALSE
  cfg$outputs$gradients <- FALSE
  slug <- penv$short_case_slug(case_id)
  res <- penv$senv$run_hmsc_case(root, paste0(slug, "__hmsc_small_mcmc_traits_phylo"),
                                 dat, cfg, dat$Y,
                                 "Small MCMC Hmsc sanity run: 50 samples, 50 transient, 2 chains, traits and taxonomy-derived phylogeny enabled.")
  finalize(res, dat, case_id, "Hmsc", cfg$model$TrFormula, phylo_info, "samples=50; transient=50; thin=1; chains=2")
}

run_hpc <- function(case_id) {
  dat <- penv$load_public_dat(case_id, "02_Hmsc-HPC")
  dat <- penv$subset_public_dat_for_fit(dat, "Hmsc-HPC_small_mcmc_traits_phylogeny", max_species = 12L)
  phylo_info <- list(repaired = FALSE, min_eigen_before = NA_real_, min_eigen_after = NA_real_,
                     note = "taxonomy rank table used; derived matrix written by JSDM Studio")
  cfg <- penv$adapt_hmschpc_cfg(dat)
  cfg$model$use_traits <- TRUE
  cfg$model$trait_formula <- stable_trait_formula(dat$traits, max_terms = 2L)
  cfg$model$phylogeny_mode <- "taxonomy"
  cfg$random_effects$mode <- "none"
  cfg$random_effects$name <- "no_random"
  cfg$sampler$samples <- 30L
  cfg$sampler$transient <- 30L
  cfg$sampler$thin <- 1L
  cfg$sampler$chains <- 2L
  cfg$sampler$verbose <- 10L
  slug <- penv$short_case_slug(case_id)
  res <- penv$senv$run_hmschpc_case(root, paste0(slug, "__hpc_small_mcmc_traits_phylo"),
                                    dat, cfg, dat$Y,
                                    "Small MCMC Hmsc-HPC CPU sanity run: 30 samples, 30 transient, 2 chains, traits and taxonomy-derived phylogeny enabled.")
  finalize(res, dat, case_id, "Hmsc-HPC", cfg$model$trait_formula, phylo_info, "samples=30; transient=30; thin=1; chains=2")
}

case_ids <- c("cestes_birds_barbaro2009b", "masviken_mountain_plants", "phylacine_mammals_range_derived")
results <- list()
for (case_id in case_ids) {
  message("Running small MCMC Hmsc case: ", case_id)
  results[[length(results) + 1L]] <- tryCatch(run_hmsc(case_id), error = function(e) {
    data.frame(engine = "Hmsc", case_id = case_id, status = "fit_failed", output_folder = NA_character_,
               zip = NA_character_, zip_size_kb = NA_real_, file_count = NA_integer_, missing_count = NA_integer_,
               missing_items = "", failure_reason = conditionMessage(e), warnings = "",
               notes = "Unhandled small MCMC Hmsc audit error.", runtime_seconds = NA_real_, stringsAsFactors = FALSE)
  })
  message("Running small MCMC Hmsc-HPC case: ", case_id)
  results[[length(results) + 1L]] <- tryCatch(run_hpc(case_id), error = function(e) {
    data.frame(engine = "Hmsc-HPC", case_id = case_id, status = "fit_failed", output_folder = NA_character_,
               zip = NA_character_, zip_size_kb = NA_real_, file_count = NA_integer_, missing_count = NA_integer_,
               missing_items = "", failure_reason = conditionMessage(e), warnings = "",
               notes = "Unhandled small MCMC Hmsc-HPC audit error.", runtime_seconds = NA_real_, stringsAsFactors = FALSE)
  })
}

summary <- do.call(rbind, results)
write.csv(summary, file.path(root, "public_small_mcmc_trait_phylogeny_summary.csv"), row.names = FALSE)
writeLines(c(
  "# Public small MCMC traits/phylogeny audit",
  "",
  "This is a sanity MCMC run, not a publication-level convergence analysis.",
  "Hmsc settings: samples=50, transient=50, thin=1, chains=2.",
  "Hmsc-HPC settings: samples=30, transient=30, thin=1, chains=2.",
  "All tested branches enable traits and taxonomy-derived phylogeny.",
  "",
  "See public_small_mcmc_trait_phylogeny_summary.csv for status."
), file.path(root, "README_public_small_mcmc_trait_phylogeny_audit.md"), useBytes = TRUE)

zipfile <- penv$runner_env$make_zip_file(root, file.path(app_dir, "output", paste0("JSDMStudio_public_small_mcmc_trait_phylogeny_audit_", timestamp, ".zip")))
cat("PUBLIC_SMALL_MCMC_AUDIT_ROOT=", normalizePath(root, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("PUBLIC_SMALL_MCMC_AUDIT_ZIP=", normalizePath(zipfile, winslash = "/", mustWork = FALSE), "\n", sep = "")
print(table(summary$engine, summary$status, useNA = "ifany"))
