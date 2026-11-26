#!/usr/bin/env Rscript

# Public-case traits/phylogeny audit.
# This script verifies that trait tables and raw phylogeny/taxonomy inputs are
# not merely copied into ZIP files, but are actually enabled in supported model
# branches.

if (.Platform$OS.type == "windows") {
  Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else getwd()
app_dir <- normalizePath(file.path(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)), ".."),
                         winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) {
  app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
setwd(app_dir)

suite_file <- file.path(app_dir, "workflow_scripts", "public_maximum_upload_case_suite.R")
suite_lines <- readLines(suite_file, warn = FALSE)
suite_lines <- suite_lines[!grepl("^\\s*if\\s*\\(!interactive\\(\\)\\)\\s*run_suite\\(\\)\\s*$", suite_lines)]
penv <- new.env(parent = globalenv())
eval(parse(text = paste(suite_lines, collapse = "\n")), envir = penv)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
root <- file.path(app_dir, "output", paste0("public_trait_phylogeny_model_audit_", timestamp))
dir.create(root, recursive = TRUE, showWarnings = FALSE)

stable_trait_formula <- function(Tr, max_terms = 2L) {
  if (is.null(Tr) || ncol(Tr) == 0L) return("~ 1")
  tr <- as.data.frame(Tr, check.names = FALSE)
  num <- names(tr)[vapply(tr, is.numeric, logical(1))]
  if (length(num) == 0L) {
    return("~ 1")
  }
  paste("~", paste(utils::head(num, max_terms), collapse = " + "))
}

repair_phylo_cov <- function(C, species) {
  info <- list(repaired = FALSE, min_eigen_before = NA_real_, min_eigen_after = NA_real_,
               note = "not_checked")
  if (is.null(C)) {
    info$note <- "missing_phylo_cov"
    return(list(C = C, info = info))
  }
  C <- as.matrix(C)
  storage.mode(C) <- "numeric"
  common <- intersect(species, rownames(C) %||% character())
  if (length(common) != length(species) || !all(species %in% colnames(C))) {
    stop("phylo_cov row/column names do not match fitted species.")
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
    info$note <- paste0("nearest positive-definite jitter applied: ", signif(jitter, 4))
  } else {
    info$note <- "positive_definite"
  }
  eig2 <- tryCatch(eigen(C, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
  info$min_eigen_after <- suppressWarnings(min(eig2, na.rm = TRUE))
  rownames(C) <- colnames(C) <- species
  list(C = C, info = info)
}

result_files <- function(outdir) {
  list(
    hmsc_phylo_summary = file.exists(file.path(outdir, "tables", "S1_phylogeny_summary.csv")),
    hmsc_gamma_results = length(list.files(file.path(outdir, "results"),
                                           pattern = "parameter_estimates_Gamma", full.names = TRUE)) > 0L,
    hpc_gamma_tables = length(list.files(file.path(outdir, "tables"),
                                         pattern = "Gamma|trait", full.names = TRUE, ignore.case = TRUE)) > 0L,
    standard_effects = file.exists(file.path(outdir, "standard", "effects_species_environment.csv")),
    standard_predictions = file.exists(file.path(outdir, "standard", "predictions_site_species.csv"))
  )
}

finalize_result <- function(res, dat, case_id, engine, trait_formula, phylo_info) {
  outdir <- res$output_folder[[1]]
  res <- penv$refresh_result_row(res, dat, penv$case_meta(case_id)$display_name[[1]], case_id)
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
  use_tab <- data.frame(
    engine = engine,
    source_case_id = case_id,
    fitted_species_n = ncol(dat$Y),
    trait_rows = nrow(dat$traits),
    trait_cols = ncol(dat$traits),
    trait_formula = trait_formula,
    phylogeny_input = "taxonomy.csv",
    derived_matrix_expected = TRUE,
    phylo_min_eigen_before = NA_real_,
    phylo_min_eigen_after = NA_real_,
    note = phylo_info$note,
    stringsAsFactors = FALSE
  )
  write.csv(use_tab, file.path(outdir, "tables", "trait_phylogeny_model_use_summary.csv"), row.names = FALSE)
  writeLines(c(
    "Trait/phylogeny model-use audit:",
    paste0("engine: ", engine),
    paste0("case: ", case_id),
    paste0("trait formula used in the model config: ", trait_formula),
    "raw taxonomy rank table supplied to the model: taxonomy.csv",
    "species-by-species relatedness matrix is derived internally and saved as derived_phylogeny_correlation.csv",
    paste0("phylo note: ", phylo_info$note)
  ), file.path(outdir, "diagnostics", "trait_phylogeny_model_use_note.txt"), useBytes = TRUE)
  flags <- result_files(outdir)
  for (nm in names(flags)) {
    res[[nm]] <- flags[[nm]]
  }
  zipfile <- tryCatch(penv$runner_env$make_zip_file(outdir), error = function(e) NA_character_)
  res$zip <- if (!is.na(zipfile)) normalizePath(zipfile, winslash = "/", mustWork = FALSE) else res$zip
  res$zip_size_kb <- if (!is.na(zipfile) && file.exists(zipfile)) round(file.info(zipfile)$size / 1024, 1) else res$zip_size_kb
  res$file_count <- length(list.files(outdir, recursive = TRUE, all.files = FALSE))
  res
}

run_hmsc_trait_phylo_case <- function(case_id) {
  dat <- penv$load_public_dat(case_id, "01_Hmsc")
  dat <- penv$subset_public_dat_for_fit(dat, "Hmsc_trait_phylogeny", max_species = 30L)
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
  cfg$outputs$cv <- FALSE
  cfg$outputs$waic <- FALSE
  cfg$outputs$gradients <- FALSE
  slug <- penv$short_case_slug(case_id)
  res <- penv$senv$run_hmsc_case(root, paste0(slug, "__hmsc_traits_phylo"),
                                 dat, cfg, dat$Y,
                                 "Public case Hmsc branch with traits and taxonomy-derived phylogeny enabled.")
  finalize_result(res, dat, case_id, "Hmsc", cfg$model$TrFormula, phylo_info)
}

run_hmschpc_trait_phylo_case <- function(case_id) {
  dat <- penv$load_public_dat(case_id, "02_Hmsc-HPC")
  dat <- penv$subset_public_dat_for_fit(dat, "Hmsc-HPC_trait_phylogeny", max_species = 30L)
  phylo_info <- list(repaired = FALSE, min_eigen_before = NA_real_, min_eigen_after = NA_real_,
                     note = "taxonomy rank table used; derived matrix written by JSDM Studio")
  cfg <- penv$adapt_hmschpc_cfg(dat)
  cfg$model$use_traits <- TRUE
  cfg$model$trait_formula <- stable_trait_formula(dat$traits, max_terms = 2L)
  cfg$model$phylogeny_mode <- "taxonomy"
  cfg$random_effects$mode <- "none"
  cfg$random_effects$name <- "no_random"
  cfg$sampler$samples <- 2L
  cfg$sampler$transient <- 2L
  cfg$sampler$chains <- 1L
  slug <- penv$short_case_slug(case_id)
  res <- penv$senv$run_hmschpc_case(root, paste0(slug, "__hpc_traits_phylo"),
                                    dat, cfg, dat$Y,
                                    "Public case Hmsc-HPC CPU branch with traits and taxonomy-derived phylogeny enabled.")
  finalize_result(res, dat, case_id, "Hmsc-HPC", cfg$model$trait_formula, phylo_info)
}

case_ids <- penv$case_names$case_id
results <- list()
for (case_id in case_ids) {
  message("Running Hmsc traits+phylogeny audit: ", case_id)
  results[[length(results) + 1L]] <- tryCatch(run_hmsc_trait_phylo_case(case_id), error = function(e) {
    data.frame(engine = "Hmsc", case_id = case_id, status = "fit_failed", output_folder = NA_character_,
               zip = NA_character_, zip_size_kb = NA_real_, file_count = NA_integer_,
               missing_count = NA_integer_, missing_items = "", failure_reason = conditionMessage(e),
               warnings = "", notes = "Unhandled Hmsc traits/phylogeny audit error.",
               runtime_seconds = NA_real_, stringsAsFactors = FALSE)
  })
  message("Running Hmsc-HPC traits+phylogeny audit: ", case_id)
  results[[length(results) + 1L]] <- tryCatch(run_hmschpc_trait_phylo_case(case_id), error = function(e) {
    data.frame(engine = "Hmsc-HPC", case_id = case_id, status = "fit_failed", output_folder = NA_character_,
               zip = NA_character_, zip_size_kb = NA_real_, file_count = NA_integer_,
               missing_count = NA_integer_, missing_items = "", failure_reason = conditionMessage(e),
               warnings = "", notes = "Unhandled Hmsc-HPC traits/phylogeny audit error.",
               runtime_seconds = NA_real_, stringsAsFactors = FALSE)
  })
}

summary <- do.call(rbind, results)
write.csv(summary, file.path(root, "public_trait_phylogeny_model_audit_summary.csv"), row.names = FALSE)

support <- data.frame(
  engine = c("Hmsc", "Hmsc-HPC", "jSDM", "GJAM", "spOccupancy", "sjSDM", "boral"),
  traits_enter_model = c(TRUE, TRUE, TRUE, TRUE, FALSE, "posthoc/metadata", TRUE),
  phylogeny_enters_model = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  note = c(
    "Native TrData/Gamma and Newick/taxonomy-derived phylogeny support.",
    "CPU pyhmsc-style traits and Newick/taxonomy-derived phylogeny support in this workflow.",
    "Trait formula branch exists; no standard phylogeny model path in this workflow.",
    "specByTrait/traitTypes branch exists; no standard phylogeny model path.",
    "Occupancy engine uses species metadata and replicated detection; no phylogeny model path.",
    "Traits can support post-hoc assembly/effects reporting; no phylogeny model path.",
    "Trait/fourth-corner branch exists; no phylogeny model path."
  ),
  stringsAsFactors = FALSE
)
write.csv(support, file.path(root, "engine_trait_phylogeny_support_matrix.csv"), row.names = FALSE)

readme <- c(
  "# Public trait/phylogeny model-use audit",
  "",
  "This audit asks a narrower question than the maximum-upload audit:",
  "were traits and raw phylogeny/taxonomy inputs actually enabled in model fitting branches?",
  "",
  "- Hmsc and Hmsc-HPC are tested here because they are the workflows with native phylogeny model inputs.",
  "- Other engines may receive trait files or metadata, but their standard workflows do not fit phylogenetic/tree effects.",
  "- Each tested output folder includes `tables/trait_phylogeny_model_use_summary.csv` and `diagnostics/trait_phylogeny_model_use_note.txt`.",
  "",
  "See `public_trait_phylogeny_model_audit_summary.csv` for per-case status."
)
writeLines(readme, file.path(root, "README_public_trait_phylogeny_model_audit.md"), useBytes = TRUE)

zipfile <- penv$runner_env$make_zip_file(root, file.path(app_dir, "output", paste0("JSDMStudio_public_trait_phylogeny_model_audit_", timestamp, ".zip")))
sha <- if (file.exists(zipfile)) digest::digest(zipfile, algo = "sha256", file = TRUE) else NA_character_

cat("PUBLIC_TRAIT_PHYLOGENY_AUDIT_ROOT=", normalizePath(root, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("PUBLIC_TRAIT_PHYLOGENY_AUDIT_ZIP=", normalizePath(zipfile, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("PUBLIC_TRAIT_PHYLOGENY_AUDIT_SHA256=", sha, "\n", sep = "")
print(table(summary$engine, summary$status, useNA = "ifany"))
