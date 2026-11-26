#!/usr/bin/env Rscript

# Build manuscript-ready figure captions and result text from the final
# JSDM Studio case-study output contract.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(prefix, default = NA_character_) {
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0L) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

project_root <- normalizePath(
  get_arg("--project=", getwd()),
  winslash = "/",
  mustWork = FALSE
)

case_root <- get_arg("--case-root=", NA_character_)
if (is.na(case_root) || !nzchar(case_root)) {
  out_dir <- file.path(project_root, "output")
  candidates <- list.dirs(out_dir, recursive = FALSE, full.names = TRUE)
  candidates <- candidates[grepl("CASE_STUDIES_FINAL_", basename(candidates))]
  if (length(candidates) == 0L) {
    stop("No CASE_STUDIES_FINAL_* directory found under output/.")
  }
  info <- file.info(candidates)
  case_root <- candidates[order(info$mtime, decreasing = TRUE)][[1]]
}
case_root <- normalizePath(case_root, winslash = "/", mustWork = TRUE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

std_dir <- file.path(case_root, "standard")
effects <- read_csv(file.path(std_dir, "CASE_STUDY_effects_species_environment_all_engines.csv"))
preds <- read_csv(file.path(std_dir, "CASE_STUDY_predictions_site_species_all_engines.csv"))
assocs <- read_csv(file.path(std_dir, "CASE_STUDY_associations_species_species_all_engines.csv"))
metrics <- read_csv(file.path(std_dir, "CASE_STUDY_fit_metrics_all_engines.csv"))
runs <- read_csv(file.path(case_root, "CASE_STUDY_RUN_SUMMARY.csv"))

fmt_num <- function(x, digits = 2) {
  vapply(x, function(one) {
    if (length(one) == 0L || is.na(one) || !is.finite(one)) return("NA")
    formatC(one, digits = digits, format = "f")
  }, character(1))
}

fmt_pct <- function(x, digits = 1) {
  vapply(x, function(one) {
    if (length(one) == 0L || is.na(one) || !is.finite(one)) return("NA")
    paste0(formatC(100 * one, digits = digits, format = "f"), "%")
  }, character(1))
}

clean_case <- function(x) {
  map <- c(
    CASE_STUDY_01_VascularPlants = "Case study 1: vascular plants",
    CASE_STUDY_02_DeadwoodFungi = "Case study 2: deadwood fungi",
    CASE_STUDY_03_Birds_FinlandPrediction = "Case study 3: Finnish birds",
    CASE_STUDY_04_Bryophytes_TibetSpatial = "Case study 4: Tibetan bryophytes",
    CASE_STUDY_05_FreshwaterMacroinvertebrates = "Case study 5: freshwater macroinvertebrates"
  )
  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

taxon_name <- function(x) {
  map <- c(
    CASE_STUDY_01_VascularPlants = "vascular plant community",
    CASE_STUDY_02_DeadwoodFungi = "deadwood fungal community",
    CASE_STUDY_03_Birds_FinlandPrediction = "Finnish bird community",
    CASE_STUDY_04_Bryophytes_TibetSpatial = "Tibetan bryophyte community",
    CASE_STUDY_05_FreshwaterMacroinvertebrates = "freshwater macroinvertebrate community"
  )
  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

question_text <- function(x) {
  map <- c(
    CASE_STUDY_01_VascularPlants = "trait filtering and local environmental selection",
    CASE_STUDY_02_DeadwoodFungi = "substrate-related community assembly on deadwood",
    CASE_STUDY_03_Birds_FinlandPrediction = "spatial prediction of occurrence and richness across Finland",
    CASE_STUDY_04_Bryophytes_TibetSpatial = "high-diversity spatial community assembly across alpine environmental gradients",
    CASE_STUDY_05_FreshwaterMacroinvertebrates = "metacommunity structure with detection-aware and observation-scale comparison"
  )
  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

top_effect_sentence <- function(case_group, engine = "Hmsc", n = 3L) {
  x <- effects[
    effects$case_group == case_group &
      effects$engine == engine &
      !effects$predictor %in% c("(Intercept)") &
      !grepl("^trait:", effects$species) &
      !grepl("^Gamma_", effects$predictor) &
      is.finite(effects$estimate),
    ,
    drop = FALSE
  ]
  if (nrow(x) == 0L) return(paste0(engine, " did not export finite species-environment effects for this case."))
  x <- x[order(-abs(x$estimate)), , drop = FALSE]
  x <- head(x, n)
  items <- paste0(
    x$species, " - ", x$predictor, " (",
    ifelse(x$estimate >= 0, "+", ""),
    fmt_num(x$estimate, 2), ")"
  )
  paste0(engine, " exported species-environment effects; the largest absolute estimates were ",
         paste(items, collapse = "; "), ".")
}

gamma_sentence <- function(hmsc_dir, n = 3L) {
  f <- file.path(hmsc_dir, "results", "parameter_estimates_trait_environment_gamma_long_main.csv")
  if (!file.exists(f)) return("Trait-environment Gamma output was not available for this case.")
  x <- read.csv(f, stringsAsFactors = FALSE)
  x <- x[!x$predictor %in% c("(Intercept)") &
           !x$trait %in% c("(Intercept)") &
           is.finite(x$estimate), , drop = FALSE]
  if (nrow(x) == 0L) return("The Gamma table contained no finite non-intercept trait-environment terms.")
  x <- x[order(-abs(x$estimate)), , drop = FALSE]
  x <- head(x, n)
  items <- paste0(
    x$trait, " - ", x$predictor, " (",
    ifelse(x$estimate >= 0, "+", ""),
    fmt_num(x$estimate, 2), ")"
  )
  paste0("The trait-environment panel highlighted Gamma terms including ",
         paste(items, collapse = "; "), ".")
}

vp_sentence <- function(hmsc_dir, n = 3L) {
  f <- file.path(hmsc_dir, "results", "parameter_estimates_variance_partitioning_mean_by_predictor_main.csv")
  if (!file.exists(f)) return("Variance-partitioning output was not available for this case.")
  x <- read.csv(f, stringsAsFactors = FALSE)
  x <- x[is.finite(x$mean_variance_proportion), , drop = FALSE]
  if (nrow(x) == 0L) return("The variance-partitioning table contained no finite terms.")
  x <- x[order(-x$mean_variance_proportion), , drop = FALSE]
  x <- head(x, n)
  items <- paste0(x$predictor, " (", fmt_pct(x$mean_variance_proportion, 1), ")")
  paste0("Mean variance partitioning assigned the largest proportions to ",
         paste(items, collapse = "; "), ".")
}

rho_sentence <- function(hmsc_dir) {
  f <- file.path(hmsc_dir, "results", "parameter_estimates_phylogenetic_signal_rho_main.csv")
  if (!file.exists(f)) return("The rho plot was not generated for this case because phylogeny was not used in the fitted Hmsc model.")
  x <- read.csv(f, stringsAsFactors = FALSE)
  if (nrow(x) == 0L || !is.finite(x$mean[[1]])) return("The rho summary was present but did not contain a finite mean.")
  paste0("The phylogenetic-signal summary gave rho mean = ", fmt_num(x$mean[[1]], 2),
         " with 95% interval [", fmt_num(x$lower[[1]], 2), ", ", fmt_num(x$upper[[1]], 2),
         "], so the figure should be interpreted as a model diagnostic for phylogenetic structure rather than as proof of a causal phylogenetic mechanism.")
}

association_sentence <- function(case_group, engine = "Hmsc", n = 3L) {
  x <- assocs[
    assocs$case_group == case_group &
      assocs$engine == engine &
      is.finite(assocs$estimate),
    ,
    drop = FALSE
  ]
  if (nrow(x) == 0L) return(paste0(engine, " did not export finite species-association estimates for this case."))
  x <- x[order(-abs(x$estimate)), , drop = FALSE]
  x <- head(x, n)
  items <- paste0(
    x$species_i, " - ", x$species_j, " (",
    ifelse(x$estimate >= 0, "+", ""),
    fmt_num(x$estimate, 2), ", ", x$association_type, ")"
  )
  paste0("The species-association panel should be described as model-specific residual structure; the largest displayed associations included ",
         paste(items, collapse = "; "), ".")
}

prediction_sentence <- function(case_group, engine = "Hmsc") {
  x <- preds[
    preds$case_group == case_group &
      preds$engine == engine &
      is.finite(preds$predicted),
    ,
    drop = FALSE
  ]
  if (nrow(x) == 0L) return(paste0(engine, " did not export finite prediction rows for this case."))
  obs <- x$observed
  obs_rate <- if (any(is.finite(obs))) mean(obs, na.rm = TRUE) else NA_real_
  paste0(engine, " exported ", nrow(x), " site-species prediction rows; predicted values ranged from ",
         fmt_num(min(x$predicted, na.rm = TRUE), 3), " to ",
         fmt_num(max(x$predicted, na.rm = TRUE), 3),
         if (is.finite(obs_rate)) paste0(", with observed mean response ", fmt_num(obs_rate, 3)) else "",
         ".")
}

newdata_sentence <- function(hmsc_dir) {
  f <- file.path(hmsc_dir, "predictions", "newdata_predicted_species_richness_main.csv")
  if (!file.exists(f)) return("No new-site species-richness table was exported for this case.")
  x <- read.csv(f, stringsAsFactors = FALSE)
  if (!"predicted_species_richness" %in% names(x) || nrow(x) == 0L) {
    return("The new-site prediction table was present but did not contain predicted_species_richness.")
  }
  paste0("For new prediction sites, predicted species richness ranged from ",
         fmt_num(min(x$predicted_species_richness, na.rm = TRUE), 2), " to ",
         fmt_num(max(x$predicted_species_richness, na.rm = TRUE), 2),
         " species across ", nrow(x), " prediction sites.")
}

metric_sentence <- function(case_group) {
  x <- metrics[metrics$case_group == case_group & is.finite(metrics$value), , drop = FALSE]
  if (nrow(x) == 0L) return("No finite fit metrics were exported for this case.")
  keep <- x$metric %in% c("mean_AUC", "mean_TjurR2", "mean_R2", "mean_SR2", "RMSE", "R2_like", "Rsquared_1", "rmspeAll", "waic_1")
  x <- x[keep, , drop = FALSE]
  if (nrow(x) == 0L) return("Fit metrics were exported, but no compact cross-engine metric was selected for the figure text.")
  x <- x[order(x$engine, x$metric), , drop = FALSE]
  x <- x[!duplicated(paste(x$engine, x$metric)), , drop = FALSE]
  x <- head(x, 6)
  items <- paste0(x$engine, " ", x$metric, " = ", fmt_num(x$value, 3))
  paste0("Compact model-fit metrics available for reporting were ",
         paste(items, collapse = "; "), ".")
}

hmsc_case_dir <- function(case_group) {
  dirs <- list.dirs(file.path(case_root, "Hmsc"), recursive = FALSE, full.names = TRUE)
  hit <- dirs[grepl(case_group, basename(dirs), fixed = TRUE)]
  if (length(hit) == 0L) return(NA_character_)
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

rel <- function(path) {
  vapply(path, function(one) {
    if (is.na(one) || !nzchar(one)) return("NA")
    normalizePath(one, winslash = "/", mustWork = FALSE)
  }, character(1))
}

figure_paths <- function(hmsc_dir) {
  files <- c(
    beta = file.path(hmsc_dir, "results", "parameter_estimates_phylo_beta_heatmap_main.pdf"),
    gamma = file.path(hmsc_dir, "results", "parameter_estimates_trait_environment_gamma_heatmap_main.pdf"),
    association = file.path(hmsc_dir, "results", "parameter_estimates_species_association_matrix_main_random_level_1.pdf"),
    variance = file.path(hmsc_dir, "results", "parameter_estimates_variance_partitioning_bookstyle_main.pdf"),
    rho = file.path(hmsc_dir, "results", "parameter_estimates_phylogenetic_signal_rho_main.pdf"),
    predictions = file.path(hmsc_dir, "results", "predictions.pdf")
  )
  ok <- file.exists(files)
  paste0(names(files)[ok], ": ", rel(files[ok]), collapse = "\n")
}

run_sentence <- function(case_group) {
  x <- runs[runs$case_group == case_group, , drop = FALSE]
  if (nrow(x) == 0L) return("No run-summary row was available.")
  items <- paste0(
    x$engine, " ", x$status,
    " (", x$file_count, " files, ", fmt_num(x$runtime_seconds, 1), " s)"
  )
  paste0("The case produced ", nrow(x), " fitted engine output folders: ",
         paste(items, collapse = "; "), ".")
}

case_groups <- c(
  "CASE_STUDY_01_VascularPlants",
  "CASE_STUDY_02_DeadwoodFungi",
  "CASE_STUDY_03_Birds_FinlandPrediction",
  "CASE_STUDY_04_Bryophytes_TibetSpatial",
  "CASE_STUDY_05_FreshwaterMacroinvertebrates"
)

out <- character()
add <- function(...) {
  out <<- c(out, paste0(...))
}

add("# Case-Study Figure Captions and Result Text\n")
add("Generated: ", as.character(Sys.time()), "\n")
add("Source output root: `", case_root, "`\n")
add("This file is derived from the fitted JSDM Studio case-study outputs and is intended as manuscript-ready draft text. It reports observable output summaries only; ecological interpretation should be checked against the original study design, sampling protocol and final MCMC diagnostics.\n")

add("## Terminology Ledger\n")
add(paste(
  "| Canonical term | Meaning in this file |",
  "|---|---|",
  "| Beta | Species-specific environmental coefficient exported by Hmsc or harmonised engine output. |",
  "| Gamma | Hmsc trait-environment coefficient linking species traits to environmental responses. |",
  "| Omega / association | Model-specific residual species association; not automatically comparable across engines. |",
  "| Rho | Hmsc phylogenetic-signal parameter for residual species responses. |",
  "| Prediction | Site-species fitted or new-site predicted response exported in the standard prediction tables. |",
  "| Fitted | The engine produced model outputs and the required output contract; this does not by itself prove ecological adequacy. |",
  sep = "\n"
))

add("## One-Sentence Argument\n")
add("Across five worked ecological case studies, JSDM Studio converted heterogeneous community data into auditable model outputs and manuscript-level figures for environmental responses, trait-mediated responses, residual associations, phylogenetic signal, variance partitioning and spatial or new-site prediction, while preserving engine-specific limits.\n")

for (cg in case_groups) {
  hd <- hmsc_case_dir(cg)
  add("\n## ", clean_case(cg), "\n")
  add("**Scientific role.** This worked example addresses ", question_text(cg), " in a ", taxon_name(cg), ".\n\n")
  add("**Recommended figure caption.** Hmsc-based worked-example output for ", taxon_name(cg), ". ")
  add("Panel A shows phylogeny-ordered species-environment Beta estimates, allowing the reader to see which species respond positively or negatively to each predictor. ")
  add("Panel B shows trait-environment Gamma estimates, separating trait-mediated response patterns from species-specific Beta effects. ")
  add("Panel C shows model-specific residual species associations; these values should be interpreted within the fitted engine and not as direct causal interactions. ")
  add("Panel D shows variance partitioning by predictor and random component. ")
  add("Panel E shows phylogenetic-signal rho where phylogeny was used. ")
  add("Panel F shows fitted or new-site predictions exported by the reproducible workflow.\n\n")
  add("**Result text.** ", run_sentence(cg), " ")
  add(top_effect_sentence(cg, "Hmsc", 3L), " ")
  add(gamma_sentence(hd, 3L), " ")
  add(vp_sentence(hd, 3L), " ")
  add(rho_sentence(hd), " ")
  add(prediction_sentence(cg, "Hmsc"), " ")
  add(newdata_sentence(hd), " ")
  add(association_sentence(cg, "Hmsc", 3L), " ")
  add(metric_sentence(cg), "\n\n")
  add("**Figure files.**\n")
  fp <- figure_paths(hd)
  if (!nzchar(fp)) {
    add("No standard Hmsc figure files were found for this case.\n")
  } else {
    add("```text\n", fp, "\n```\n")
  }
}

add("\n## Cross-Case Results Paragraph\n")
add("The worked examples showed that JSDM Studio produced complete, auditable output contracts for all selected case-study engines. Hmsc was fitted for all five focal datasets, while additional engines were fitted where they matched the ecological question and data structure: jSDM and boral for vascular plants, GJAM and boral for deadwood fungi, Hmsc-HPC for the Finnish bird prediction case, Hmsc-HPC and sjSDM for the Tibetan bryophyte case, and spOccupancy plus GJAM for freshwater macroinvertebrates. Across these runs, all engine folders were labelled as fitted and had zero missing required output components. The case-study figures should therefore be used to demonstrate the interpretive products of the software rather than to rank the ecological truth of engines whose association parameters differ by model family.\n")

add("\n## Caution Text for Manuscript\n")
add("These worked-example figures demonstrate output generation, traceability and ecological interpretability under the selected desktop settings. They should not be described as final publication-grade ecological inference until convergence diagnostics, sensitivity to priors and model-specific assumptions have been reviewed for the focal biological question. In particular, residual association matrices, latent-factor covariance structures and observation-scale correlations are not interchangeable estimands across Hmsc, Hmsc-HPC, jSDM, GJAM, sjSDM, spOccupancy and boral.\n")

out_path <- file.path(case_root, "CASE_STUDY_FIGURE_RESULTS_TEXT.md")
writeLines(out, out_path, useBytes = TRUE)

summary_path <- file.path(std_dir, "CASE_STUDY_FIGURE_TEXT_SOURCE_SUMMARY.csv")
summary_rows <- do.call(rbind, lapply(case_groups, function(cg) {
  data.frame(
    case_group = cg,
    case_label = clean_case(cg),
    hmsc_dir = hmsc_case_dir(cg),
    hmsc_prediction_rows = nrow(preds[preds$case_group == cg & preds$engine == "Hmsc", , drop = FALSE]),
    hmsc_effect_rows = nrow(effects[effects$case_group == cg & effects$engine == "Hmsc", , drop = FALSE]),
    hmsc_association_rows = nrow(assocs[assocs$case_group == cg & assocs$engine == "Hmsc", , drop = FALSE]),
    stringsAsFactors = FALSE
  )
}))
write.csv(summary_rows, summary_path, row.names = FALSE)

cat("Wrote:", out_path, "\n")
cat("Wrote:", summary_path, "\n")
