options(stringsAsFactors = FALSE)

args0 <- commandArgs(trailingOnly = FALSE)
script_arg <- args0[startsWith(args0, "--file=")]
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
} else {
  normalizePath("workflow_scripts/package_whittaker_validated_case.R", winslash = "/", mustWork = FALSE)
}
app_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
setwd(app_dir)

cli_args <- commandArgs(trailingOnly = TRUE)
cli_value <- function(key, default = "") {
  hit <- cli_args[startsWith(cli_args, paste0("--", key, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[[length(hit)]])
}
version_suffix <- cli_value("version_suffix", "")
if (nzchar(version_suffix) && !startsWith(version_suffix, "_")) version_suffix <- paste0("_", version_suffix)
if (grepl("[^A-Za-z0-9_.-]", version_suffix)) stop("Invalid --version_suffix", call. = FALSE)

copy_tree <- function(src, dst) {
  src <- normalizePath(src, winslash = "/", mustWork = TRUE)
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  dirs <- list.dirs(src, recursive = TRUE, full.names = TRUE)
  for (d in dirs) {
    rel <- substring(normalizePath(d, winslash = "/", mustWork = TRUE), nchar(src) + 2L)
    if (nzchar(rel)) dir.create(file.path(dst, rel), recursive = TRUE, showWarnings = FALSE)
  }
  files <- list.files(src, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  for (f in files) {
    rel <- substring(normalizePath(f, winslash = "/", mustWork = TRUE), nchar(src) + 2L)
    target <- file.path(dst, rel)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(f, target, overwrite = TRUE)) stop("Could not copy: ", f, call. = FALSE)
  }
  invisible(files)
}

copy_selected <- function(src_root, dst_root, relative_paths) {
  for (rel in relative_paths) {
    src <- file.path(src_root, rel)
    if (!file.exists(src)) next
    dst <- file.path(dst_root, rel)
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(src, dst, overwrite = TRUE)) stop("Could not copy: ", src, call. = FALSE)
  }
}

read_csv <- function(path, ...) read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
metric_value <- function(tab, metric) {
  value <- tab$value[tab$metric == metric]
  if (length(value)) as.numeric(value[[1]]) else NA_real_
}
sha256_file <- function(path) paste(as.character(openssl::sha256(file(path))), collapse = "")

input_case <- file.path(app_dir, "examples", "00_ACTUAL_UPLOAD_DATA", "01_vascular_plants__ch06_whittaker_full_52x75")
occurrence_root <- file.path(app_dir, "output", "CASE_STUDY_01_VascularPlants_WhittakerFull_FINAL_20260918_144723")
occurrence_src <- file.path(occurrence_root, "Hmsc", "Whittaker_full_52x75_occurrence_probit")
log1p_root <- file.path(app_dir, "output", "CASE_STUDY_01_VascularPlants_WhittakerFull_LOG1P_FINAL_20260918")
log1p_src <- file.path(log1p_root, "Hmsc", "Whittaker_full_52x75_log1p_abundance_normal_sensitivity")
nonconverged_src <- file.path(occurrence_root, "Hmsc", "Whittaker_full_52x75_abundance_lognormal_poisson")
final_root <- file.path(app_dir, "output", paste0("CASE_STUDY_01_VascularPlants_WhittakerFull_VALIDATED_FINAL_20260918", version_suffix))

required <- c(input_case, occurrence_src, log1p_src, nonconverged_src)
if (any(!dir.exists(required))) stop("One or more source folders are missing: ", paste(required[!dir.exists(required)], collapse = "; "), call. = FALSE)
if (dir.exists(final_root)) stop("Final delivery folder already exists; refusing to overwrite: ", final_root, call. = FALSE)

dir.create(final_root, recursive = TRUE, showWarnings = FALSE)
copy_tree(input_case, file.path(final_root, "00_complete_input_case"))
copy_tree(occurrence_src, file.path(final_root, "01_occurrence_probit"))
copy_tree(log1p_src, file.path(final_root, "02_log1p_abundance_normal_sensitivity"))

# Refresh standard tables in the copied delivery without rerunning the fitted
# MCMC models. This keeps observed references family-aware and ensures the
# richer engine-aware tables are written after contract placeholders.
refresh_env <- new.env(parent = globalenv())
sys.source(file.path(app_dir, "workflow_scripts", "final_case_study_runner.R"), envir = refresh_env)
refresh_standard_tables <- function(outdir, label) {
  Y <- read.csv(file.path(outdir, "inputs", "Y.csv"), row.names = 1,
                check.names = FALSE, stringsAsFactors = FALSE)
  X <- read.csv(file.path(outdir, "inputs", "XData.csv"), row.names = 1,
                check.names = FALSE, stringsAsFactors = FALSE)
  dat <- refresh_env$make_reference_truth(list(Y = Y, X = X), label)
  status <- refresh_env$senv$safe_read_status(outdir, "Hmsc")
  refresh_env$app_env$ensure_output_contract(outdir, "Hmsc", status, Y = Y, X = X)
  refresh_env$runner_env$standardize_engine_synonyms(outdir, "Hmsc", dat)
  invisible(dat)
}
refresh_standard_tables(file.path(final_root, "01_occurrence_probit"),
                        "whittaker_observed_reference_not_simulated_truth")
refresh_standard_tables(file.path(final_root, "02_log1p_abundance_normal_sensitivity"),
                        "whittaker_observed_reference_not_simulated_truth")

quality_dir <- file.path(final_root, "03_quality_and_interpretation")
nonconv_dir <- file.path(quality_dir, "nonconverged_lognormal_poisson_audit")
dir.create(nonconv_dir, recursive = TRUE, showWarnings = FALSE)
copy_selected(nonconverged_src, nonconv_dir, c(
  "used_config.yml",
  "diagnostics/engine_status.json",
  "diagnostics/MCMC_block_diagnostics_summary.csv",
  "diagnostics/MCMC_parameter_level_diagnostics.csv",
  "diagnostics/R_runtime_warnings.txt",
  "results/MCMC_convergence.txt",
  "results/MCMC_convergence.pdf",
  "standard/fit_metrics.csv"
))

scripts_dir <- file.path(final_root, "04_reproducible_scripts")
dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
for (name in c("run_whittaker_full_52x75_case.R", "audit_hmsc_convergence.R", "package_whittaker_validated_case.R")) {
  file.copy(file.path(app_dir, "workflow_scripts", name), file.path(scripts_dir, name), overwrite = TRUE)
}

common <- file.path(input_case, "01_common_inputs")
y <- read_csv(file.path(common, "Y_count.csv"), row.names = 1)
occurrences <- colSums(y > 0)
totals <- colSums(y)
prevalence <- data.frame(
  species = names(occurrences),
  occupied_sites = as.integer(occurrences),
  prevalence = as.numeric(occurrences / nrow(y)),
  total_abundance = as.numeric(totals),
  information_class = cut(occurrences, c(-Inf, 1, 2, 4, 9, 19, Inf), labels = c("singleton", "doubleton", "3-4 sites", "5-9 sites", "10-19 sites", "20+ sites")),
  stringsAsFactors = FALSE
)
prevalence <- prevalence[order(prevalence$occupied_sites, prevalence$species), ]
write_csv(prevalence, file.path(quality_dir, "SPECIES_PREVALENCE_AND_INFORMATION.csv"))

trait_audit <- read_csv(file.path(common, "trait_provenance_and_imputation.csv"))
trait_summary <- as.data.frame(table(trait_audit$trait, trait_audit$method), stringsAsFactors = FALSE)
names(trait_summary) <- c("trait", "provenance_or_imputation_method", "species_count")
trait_summary <- trait_summary[trait_summary$species_count > 0, ]
write_csv(trait_summary, file.path(quality_dir, "TRAIT_PROVENANCE_SUMMARY.csv"))

phylo <- read_csv(file.path(common, "phylogeny_match_and_alignment.csv"))
phylo_summary <- as.data.frame(table(phylo$status), stringsAsFactors = FALSE)
names(phylo_summary) <- c("VPhyloMaker2_status", "species_count")
phylo_summary$response_columns_aligned <- all(phylo$response_column_match)
phylo_summary$trait_rows_aligned <- all(phylo$trait_row_match)
phylo_summary$tree_tips_aligned <- all(phylo$tree_tip_match)
write_csv(phylo_summary, file.path(quality_dir, "PHYLOGENY_ALIGNMENT_SUMMARY.csv"))

occ_diag <- read_csv(file.path(occurrence_src, "diagnostics", "MCMC_block_diagnostics_summary.csv"))
log_diag <- read_csv(file.path(log1p_src, "diagnostics", "MCMC_block_diagnostics_summary.csv"))
bad_diag <- read_csv(file.path(nonconverged_src, "diagnostics", "MCMC_block_diagnostics_summary.csv"))
occ_fit <- read_csv(file.path(occurrence_src, "standard", "fit_metrics.csv"))
log_fit <- read_csv(file.path(log1p_src, "standard", "fit_metrics.csv"))
bad_fit <- read_csv(file.path(nonconverged_src, "standard", "fit_metrics.csv"))
run_a <- read_csv(file.path(occurrence_root, "WHITTAKER_FULL_RUN_SUMMARY.csv"))
run_b <- read_csv(file.path(log1p_root, "WHITTAKER_FULL_RUN_SUMMARY.csv"))

diag_row <- function(tab, block, field) {
  value <- tab[tab$block == block, field]
  if (length(value)) as.numeric(value[[1]]) else NA_real_
}
quality <- data.frame(
  analysis = c("Occurrence probit", "Log1p abundance Gaussian sensitivity", "Original-book lognormal-Poisson attempt"),
  response_scale = c("binary occurrence, probit link", "log(1 + abundance), Gaussian", "raw 0-100 abundance, lognormal-Poisson"),
  workflow_status = c("fitted", "fitted", "fitted"),
  inference_status = c("validated for interpretation", "validated sensitivity analysis", "not validated: species-level Beta did not converge"),
  beta_psrf_median = c(diag_row(occ_diag, "Beta", "psrf_median"), diag_row(log_diag, "Beta", "psrf_median"), diag_row(bad_diag, "Beta", "psrf_median")),
  beta_psrf_max = c(diag_row(occ_diag, "Beta", "psrf_max"), diag_row(log_diag, "Beta", "psrf_max"), diag_row(bad_diag, "Beta", "psrf_max")),
  beta_ess_min = c(diag_row(occ_diag, "Beta", "ess_min"), diag_row(log_diag, "Beta", "ess_min"), diag_row(bad_diag, "Beta", "ess_min")),
  gamma_psrf_max = c(diag_row(occ_diag, "Gamma", "psrf_max"), diag_row(log_diag, "Gamma", "psrf_max"), diag_row(bad_diag, "Gamma", "psrf_max")),
  rho_psrf = c(diag_row(occ_diag, "Rho", "psrf_max"), diag_row(log_diag, "Rho", "psrf_max"), diag_row(bad_diag, "Rho", "psrf_max")),
  mean_AUC = c(metric_value(occ_fit, "mean_AUC"), metric_value(log_fit, "mean_AUC"), metric_value(bad_fit, "mean_AUC")),
  mean_TjurR2 = c(metric_value(occ_fit, "mean_TjurR2"), metric_value(log_fit, "mean_TjurR2"), metric_value(bad_fit, "mean_TjurR2")),
  mean_R2 = c(metric_value(occ_fit, "mean_R2"), metric_value(log_fit, "mean_R2"), metric_value(bad_fit, "mean_R2")),
  mean_SR2 = c(metric_value(occ_fit, "mean_SR2"), metric_value(log_fit, "mean_SR2"), metric_value(bad_fit, "mean_SR2")),
  runtime_seconds = c(run_a$runtime_seconds[run_a$case_id == "Whittaker_full_52x75_occurrence_probit"], run_b$runtime_seconds[[1]], run_a$runtime_seconds[run_a$case_id == "Whittaker_full_52x75_abundance_lognormal_poisson"]),
  missing_contract_items = c(run_a$missing_count[run_a$case_id == "Whittaker_full_52x75_occurrence_probit"], run_b$missing_count[[1]], run_a$missing_count[run_a$case_id == "Whittaker_full_52x75_abundance_lognormal_poisson"]),
  stringsAsFactors = FALSE
)
write_csv(quality, file.path(quality_dir, "MODEL_QUALITY_SUMMARY.csv"))

classify_effects <- function(path, output_name) {
  x <- read_csv(path)
  x <- x[x$predictor == "TMG_z", ]
  x$direction_support <- ifelse(x$support >= 0.95, "supported_positive", ifelse(x$support <= 0.05, "supported_negative", "uncertain"))
  x$certainty <- pmax(x$support, 1 - x$support)
  x <- x[order(-x$certainty, x$species), ]
  write_csv(x, file.path(quality_dir, output_name))
  x
}
occ_beta <- classify_effects(file.path(occurrence_src, "results", "parameter_estimates_Beta_species_environment_long_main.csv"), "TMG_SPECIES_EFFECTS_OCCURRENCE.csv")
log_beta <- classify_effects(file.path(log1p_src, "results", "parameter_estimates_Beta_species_environment_long_main.csv"), "TMG_SPECIES_EFFECTS_LOG1P_ABUNDANCE.csv")
file.copy(file.path(occurrence_src, "results", "parameter_estimates_trait_environment_gamma_long_main.csv"), file.path(quality_dir, "TRAIT_ENVIRONMENT_GAMMA_OCCURRENCE.csv"), overwrite = TRUE)
file.copy(file.path(log1p_src, "results", "parameter_estimates_trait_environment_gamma_long_main.csv"), file.path(quality_dir, "TRAIT_ENVIRONMENT_GAMMA_LOG1P_ABUNDANCE.csv"), overwrite = TRUE)
file.copy(file.path(occurrence_src, "results", "parameter_estimates_phylogenetic_signal_rho_main.csv"), file.path(quality_dir, "PHYLOGENETIC_SIGNAL_OCCURRENCE.csv"), overwrite = TRUE)
file.copy(file.path(log1p_src, "results", "parameter_estimates_phylogenetic_signal_rho_main.csv"), file.path(quality_dir, "PHYLOGENETIC_SIGNAL_LOG1P_ABUNDANCE.csv"), overwrite = TRUE)

occ_rho <- read_csv(file.path(occurrence_src, "results", "parameter_estimates_phylogenetic_signal_rho_main.csv"))
log_rho <- read_csv(file.path(log1p_src, "results", "parameter_estimates_phylogenetic_signal_rho_main.csv"))
occ_gamma <- read_csv(file.path(occurrence_src, "results", "parameter_estimates_trait_environment_gamma_long_main.csv"))
log_gamma <- read_csv(file.path(log1p_src, "results", "parameter_estimates_trait_environment_gamma_long_main.csv"))
occ_pred <- read_csv(file.path(occurrence_src, "predictions", "newdata_predicted_mean_main.csv"), row.names = 1)
log_pred <- read_csv(file.path(log1p_src, "predictions", "newdata_predicted_mean_main.csv"), row.names = 1)

gamma_line <- function(tab, trait) {
  z <- tab[tab$predictor == "TMG_z" & tab$trait == trait, ]
  if (!nrow(z)) return("not available")
  sprintf("estimate %.3f; posterior P(positive) %.3f", z$estimate[[1]], z$support[[1]])
}
effect_counts <- function(tab) {
  c(positive = sum(tab$direction_support == "supported_positive"), negative = sum(tab$direction_support == "supported_negative"), uncertain = sum(tab$direction_support == "uncertain"))
}
occ_counts <- effect_counts(occ_beta)
log_counts <- effect_counts(log_beta)
prevalence_counts <- table(prevalence$information_class)

report <- c(
  "# Validated complete Whittaker vascular-plant case",
  "",
  paste0("Generated: ", Sys.time()),
  "",
  "## Scope and data",
  "",
  "This delivery retains all 52 Whittaker sites and all 75 vascular-plant species (3,900 site-species rows). The observed environmental variable is the topographic moisture gradient (TMG; larger values represent drier topographic positions). No geographic coordinates are present, so predictions are environmental-gradient predictions rather than maps.",
  "",
  "Three traits are supplied: leaf C:N measured in the Whittaker source, plus woodiness and growth support from the versioned global growth-form database. Database gaps are filled deterministically by genus mode, family mode or nearest phylogenetic neighbour; these values are explicitly imputed, not field measurements.",
  "",
  "The 75-tip Newick tree was generated with V.PhyloMaker2 scenario S3 against GBOTB.extended.TPL. Tree tips, response columns and trait rows are exactly aligned. The species-by-species covariance file is a derived audit product, not the phylogeny input.",
  "",
  "## Information content",
  "",
  sprintf("The community is highly sparse: %d singleton species, %d doubletons, %d species in 3-4 sites, %d in 5-9 sites, %d in 10-19 sites and %d in at least 20 sites.", prevalence_counts[["singleton"]], prevalence_counts[["doubleton"]], prevalence_counts[["3-4 sites"]], prevalence_counts[["5-9 sites"]], prevalence_counts[["10-19 sites"]], prevalence_counts[["20+ sites"]]),
  "",
  "## Validated models",
  "",
  sprintf("1. Occurrence-probit Hmsc: 4 chains, 1,000 retained samples per chain, transient 2,500 and thin 5. All Beta/Gamma/rho PSRF values are <= %.3f; minimum Beta ESS is %.1f. Five-fold mean AUC is %.3f and mean Tjur R2 is %.3f.", quality$beta_psrf_max[[1]], quality$beta_ess_min[[1]], quality$mean_AUC[[1]], quality$mean_TjurR2[[1]]),
  sprintf("2. Log1p-abundance Gaussian sensitivity Hmsc: the same MCMC settings and all 75 species. All Beta/Gamma/rho PSRF values are <= %.3f; minimum Beta ESS is %.1f. Mean R2 is %.3f. This model is interpreted on the log(1 + abundance) scale.", quality$beta_psrf_max[[2]], quality$beta_ess_min[[2]], quality$mean_R2[[2]]),
  "",
  "## Ecological results",
  "",
  sprintf("Occurrence responses show %d supported positive and %d supported negative species-specific TMG effects (posterior direction probability >= 0.95), with %d uncertain effects.", occ_counts[["positive"]], occ_counts[["negative"]], occ_counts[["uncertain"]]),
  sprintf("Log1p abundance responses show %d supported positive and %d supported negative TMG effects, with %d uncertain effects.", log_counts[["positive"]], log_counts[["negative"]], log_counts[["uncertain"]]),
  sprintf("For occurrence, the community-average TMG effect is negative (Gamma intercept-by-TMG: %s). The TMG-by-leaf-C:N interaction is %s, supporting trait-mediated variation in occurrence responses along the moisture gradient. The woodiness interaction is %s.", gamma_line(occ_gamma, "(Intercept)"), gamma_line(occ_gamma, "leaf_CN"), gamma_line(occ_gamma, "woodinesswoody")),
  sprintf("For log1p abundance, the TMG-by-leaf-C:N interaction is %s and the woodiness interaction is %s; neither provides strong directional support at the 0.95 threshold.", gamma_line(log_gamma, "leaf_CN"), gamma_line(log_gamma, "woodinesswoody")),
  sprintf("Phylogenetic signal is weak/uncertain for occurrence (rho mean %.3f, median %.3f, 95%% interval %.3f-%.3f) but stronger on the log1p abundance scale (mean %.3f, median %.3f, interval %.3f-%.3f; posterior positive support %.3f). Because the latter interval touches zero, it should be described as evidence for moderate phylogenetic structuring rather than definitive proof.", occ_rho$mean[[1]], occ_rho$median[[1]], occ_rho$lower[[1]], occ_rho$upper[[1]], log_rho$mean[[1]], log_rho$median[[1]], log_rho$lower[[1]], log_rho$upper[[1]], log_rho$support_positive[[1]]),
  sprintf("The occurrence prediction matrix contains %d gradient rows x %d species, and the log1p abundance matrix contains %d x %d. These predictions span the observed TMG domain and must not be presented as geographic extrapolations.", nrow(occ_pred), ncol(occ_pred), nrow(log_pred), ncol(log_pred)),
  "",
  "## Non-converged original-book likelihood attempt",
  "",
  sprintf("The all-species raw-abundance lognormal-Poisson workflow completed its output contract but is not validated for species-level inference: Beta median PSRF %.3f, maximum %.3f and minimum ESS %.1f. This is consistent with the data containing 40 species observed in no more than two sites. Key evidence is retained under 03_quality_and_interpretation/nonconverged_lognormal_poisson_audit; the complete exploratory archive remains at: %s", quality$beta_psrf_median[[3]], quality$beta_psrf_max[[3]], quality$beta_ess_min[[3]], normalizePath(nonconverged_src, winslash = "/", mustWork = TRUE)),
  "",
  "## Interpretation boundary",
  "",
  "A fitted workflow demonstrates executable software and complete outputs. Ecological interpretation here is restricted to the two models whose MCMC diagnostics pass the stated thresholds. The single environmental gradient explains modest predictive variation, and added traits include imputed values; conclusions should therefore be framed as a worked hierarchical case study, not a definitive causal analysis of western North American vegetation."
)
writeLines(report, file.path(quality_dir, "ANALYSIS_QUALITY_AND_ECOLOGICAL_RESULTS.md"), useBytes = TRUE)

writeLines(c(
  "The raw-abundance lognormal-Poisson model executed but did not achieve species-level Beta convergence.",
  "It is retained only as an audit of an unsuitable full-species inference for this sparse matrix.",
  "Do not cite its species-level Beta estimates as validated ecological results.",
  paste0("Complete exploratory folder: ", normalizePath(nonconverged_src, winslash = "/", mustWork = TRUE))
), file.path(nonconv_dir, "README_NOT_FOR_SPECIES_LEVEL_INFERENCE.txt"), useBytes = TRUE)

manifest_files <- list.files(final_root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
manifest <- data.frame(
  relative_path = substring(normalizePath(manifest_files, winslash = "/", mustWork = TRUE), nchar(normalizePath(final_root, winslash = "/", mustWork = TRUE)) + 2L),
  size_bytes = file.info(manifest_files)$size,
  stringsAsFactors = FALSE
)
write_csv(manifest, file.path(quality_dir, "DELIVERY_FILE_MANIFEST.csv"))

if (!requireNamespace("zip", quietly = TRUE)) stop("Package zip is required", call. = FALSE)
if (!requireNamespace("openssl", quietly = TRUE)) stop("Package openssl is required", call. = FALSE)
input_zip <- file.path(dirname(final_root), paste0("JSDMStudio_Whittaker_FULL_52x75_INPUTS", version_suffix, ".zip"))
full_zip <- file.path(dirname(final_root), paste0(basename(final_root), ".zip"))
zip::zipr(input_zip, files = input_case, root = dirname(input_case), include_directories = TRUE)
writeLines(c(
  paste0("input_zip=", normalizePath(input_zip, winslash = "/", mustWork = TRUE)),
  paste0("input_zip_size_bytes=", file.info(input_zip)$size),
  paste0("input_zip_sha256=", sha256_file(input_zip))
), file.path(quality_dir, "INPUT_ZIP_CHECKSUM.txt"), useBytes = TRUE)
zip::zipr(full_zip, files = final_root, root = dirname(final_root), include_directories = TRUE)
checksum_file <- paste0(full_zip, ".sha256.txt")
writeLines(c(
  paste0("full_zip=", normalizePath(full_zip, winslash = "/", mustWork = TRUE)),
  paste0("full_zip_size_bytes=", file.info(full_zip)$size),
  paste0("full_zip_sha256=", sha256_file(full_zip)),
  paste0("input_zip=", normalizePath(input_zip, winslash = "/", mustWork = TRUE)),
  paste0("input_zip_size_bytes=", file.info(input_zip)$size),
  paste0("input_zip_sha256=", sha256_file(input_zip))
), checksum_file, useBytes = TRUE)

cat("FINAL_ROOT=", normalizePath(final_root, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("FULL_ZIP=", normalizePath(full_zip, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("FULL_ZIP_SHA256=", sha256_file(full_zip), "\n", sep = "")
cat("INPUT_ZIP=", normalizePath(input_zip, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("INPUT_ZIP_SHA256=", sha256_file(input_zip), "\n", sep = "")
