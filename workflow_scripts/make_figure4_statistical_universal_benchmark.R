#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
})

engine_order <- c("Hmsc", "Hmsc-HPC", "jSDM", "GJAM", "spOccupancy", "sjSDM", "boral")
engine_palette <- c(
  "Hmsc" = "#1f78b4",
  "Hmsc-HPC" = "#4f9a8b",
  "jSDM" = "#8b63ad",
  "GJAM" = "#d55e00",
  "spOccupancy" = "#5f7fa8",
  "sjSDM" = "#2ca25f",
  "boral" = "#7f7f7f"
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (length(x) == 1 && is.atomic(x) && is.na(x)) return(y)
  x
}

args <- commandArgs(trailingOnly = TRUE)
named_args <- list()
pos_args <- character()
for (a in args) {
  if (grepl("^--", a)) {
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "true"
    named_args[[key]] <- val
  } else {
    pos_args <- c(pos_args, a)
  }
}
cmd <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", cmd, value = TRUE))
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg[1], winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
project_root <- if (nzchar(named_args$project_root %||% "")) {
  normalizePath(named_args$project_root, winslash = "/", mustWork = TRUE)
} else if (length(pos_args) >= 1) {
  normalizePath(pos_args[1], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
}

output_root <- file.path(project_root, "docs", "paper_figures", "figure4_statistical_universal_benchmark")
source_data_dir <- file.path(output_root, "source_data")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)

safe_read_csv <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

to_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

standardize_vec <- function(x) {
  x <- to_num(x)
  if (sum(is.finite(x)) < 2) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

normalize_species <- function(x) {
  v <- trimws(as.character(x))
  v <- gsub("\\.", "_", v)
  is_sp <- grepl("^sp_?\\d+$", v, ignore.case = TRUE)
  if (any(is_sp, na.rm = TRUE)) {
    nums <- sub("^sp_?", "", v[is_sp], ignore.case = TRUE)
    v[is_sp] <- sprintf("sp_%02d", suppressWarnings(as.integer(nums)))
  }
  v
}

normalize_site_id <- function(x) {
  v <- trimws(as.character(x))
  hit_numeric <- grepl("^\\d+$", v)
  if (any(hit_numeric, na.rm = TRUE)) {
    v[hit_numeric] <- sprintf("site_%02d", suppressWarnings(as.integer(v[hit_numeric])))
  }
  hit_site <- grepl("^site_?\\d+$", v, ignore.case = TRUE)
  if (any(hit_site, na.rm = TRUE)) {
    nums <- sub("^site_?", "", v[hit_site], ignore.case = TRUE)
    v[hit_site] <- sprintf("site_%02d", suppressWarnings(as.integer(nums)))
  }
  v
}

extract_spocc_species <- function(predictor) {
  v <- trimws(as.character(predictor))
  hit <- grepl("-sp_?\\d+$", v, ignore.case = TRUE)
  out <- rep(NA_character_, length(v))
  out[hit] <- sub("^.*-(sp_?\\d+)$", "\\1", v[hit], ignore.case = TRUE)
  normalize_species(out)
}

normalize_predictor <- function(x) {
  v <- trimws(as.character(x))
  v <- sub("-sp_?\\d+$", "", v, ignore.case = TRUE)
  v <- gsub("^C\\(substrate\\)\\[T\\.(.+)\\]$", "substrate_\\1", v)
  v <- gsub("^substrate(clay|loam|rock|sand)$", "substrate_\\1", v, ignore.case = TRUE)
  v <- ifelse(tolower(v) %in% c("intercept", "(intercept)"), "(Intercept)", v)
  v
}

normalize_probability <- function(pred, link_hint = "auto") {
  pred_num <- to_num(pred)
  finite_pred <- pred_num[is.finite(pred_num)]
  if (!length(finite_pred)) {
    return(list(value = pred_num, transform = "not_available"))
  }
  if (min(finite_pred) < 0 || max(finite_pred) > 1) {
    if (identical(link_hint, "probit")) {
      return(list(value = stats::pnorm(pred_num), transform = "pnorm_from_probit_link_scale"))
    }
    return(list(value = stats::plogis(pred_num), transform = "plogis_from_logit_link_scale"))
  }
  list(value = pred_num, transform = "identity_probability_scale")
}

manual_auc <- function(obs, score) {
  ok <- is.finite(obs) & is.finite(score)
  obs <- obs[ok]
  score <- score[ok]
  if (length(obs) < 3 || length(unique(obs)) < 2) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  n_pos <- sum(obs == 1)
  n_neg <- sum(obs == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  (sum(ranks[obs == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

latest_universal_master <- function(project_root) {
  dirs <- list.dirs(file.path(project_root, "output"), recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^universal_benchmark_\\d", basename(dirs))]
  required <- file.path(dirs, "master", "04_unified_standard_results", "model_status_matrix.csv")
  dirs <- dirs[file.exists(required)]
  if (!length(dirs)) stop("No universal_benchmark_* output with unified standard files was found.")
  dirs[which.max(file.info(dirs)$mtime)]
}

read_occurrence_observed <- function(path) {
  y <- read.csv(path, row.names = 1, check.names = FALSE)
  y$site_id <- rownames(y)
  as_tibble(y) %>%
    pivot_longer(cols = -site_id, names_to = "species", values_to = "observed_input") %>%
    mutate(
      site_id_norm = normalize_site_id(site_id),
      species_norm = normalize_species(species),
      observed_input = to_num(observed_input)
    ) %>%
    select(site_id_norm, species_norm, observed_input)
}

if (nzchar(named_args$master_dir %||% "")) {
  master_dir <- normalizePath(named_args$master_dir, winslash = "/", mustWork = TRUE)
  benchmark_dir <- dirname(master_dir)
} else {
  benchmark_dir <- latest_universal_master(project_root)
  master_dir <- file.path(benchmark_dir, "master")
}
benchmark_run <- basename(benchmark_dir)

paths <- list(
  effects = file.path(master_dir, "04_unified_standard_results", "effects_all_models_long.csv"),
  predictions = file.path(master_dir, "04_unified_standard_results", "predictions_all_models_long.csv"),
  associations = file.path(master_dir, "04_unified_standard_results", "associations_all_models_long.csv"),
  status = file.path(master_dir, "04_unified_standard_results", "model_status_matrix.csv"),
  truth_effects = file.path(master_dir, "02_truth", "true_environment_effects.csv"),
  truth_predictions = file.path(master_dir, "02_truth", "true_predictions.csv"),
  truth_associations = file.path(master_dir, "02_truth", "true_species_associations.csv"),
  y_occurrence = file.path(master_dir, "01_benchmark_inputs", "Y_occurrence.csv")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing)) stop("Missing required Figure 4 source files: ", paste(missing, collapse = ", "))

engine_frame <- tibble(engine = factor(engine_order, levels = engine_order))

truth_effects <- safe_read_csv(paths$truth_effects) %>%
  mutate(
    species_norm = normalize_species(species),
    predictor_norm = normalize_predictor(predictor),
    true_effect = to_num(true_effect)
  ) %>%
  select(species_norm, predictor_norm, true_effect)

main_effect_predictors <- c("pH", "moisture", "canopy", "elevation")

effects_raw <- safe_read_csv(paths$effects) %>%
  mutate(
    estimate = to_num(estimate),
    spocc_species = extract_spocc_species(predictor),
    species_norm = ifelse(engine == "spOccupancy" & !is.na(spocc_species), spocc_species, normalize_species(species)),
    predictor_norm = normalize_predictor(predictor),
    engine = factor(engine, levels = engine_order)
  )

effects_joined <- effects_raw %>%
  filter(effect_type %in% c("environment", "occurrence"), predictor_norm != "(Intercept)", is.finite(estimate)) %>%
  group_by(engine, species_norm, predictor_norm, scale) %>%
  summarise(estimate = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
  inner_join(truth_effects, by = c("species_norm", "predictor_norm")) %>%
  filter(predictor_norm %in% main_effect_predictors) %>%
  filter(is.finite(true_effect), is.finite(estimate)) %>%
  group_by(engine) %>%
  mutate(
    estimate_std = standardize_vec(estimate),
    truth_std = standardize_vec(true_effect),
    sign_agrees = sign(estimate) == sign(true_effect)
  ) %>%
  ungroup()

effect_metrics <- effects_joined %>%
  group_by(engine) %>%
  summarise(
    n_effects = n(),
    effect_correlation = ifelse(n() >= 3, suppressWarnings(cor(true_effect, estimate)), NA_real_),
    standardized_rmse = sqrt(mean((estimate_std - truth_std)^2, na.rm = TRUE)),
    sign_accuracy = mean(sign_agrees, na.rm = TRUE),
    note = "Computed after species and continuous predictor ID matching; categorical factor contrasts remain in the standard tables but are excluded from the main cross-engine recovery score to avoid baseline-coding artefacts.",
    .groups = "drop"
  ) %>%
  right_join(engine_frame, by = "engine") %>%
  mutate(
    n_effects = replace_na(n_effects, 0L),
    note = ifelse(n_effects == 0, "No matched finite species-by-predictor effects were available.", note)
  ) %>%
  arrange(engine)

effect_metrics_long <- effect_metrics %>%
  select(engine, n_effects, effect_correlation, standardized_rmse, sign_accuracy, note) %>%
  pivot_longer(
    cols = c(effect_correlation, standardized_rmse, sign_accuracy),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_label = recode(
      metric,
      effect_correlation = "Effect r",
      standardized_rmse = "Std. RMSE",
      sign_accuracy = "Sign accuracy"
    ),
    metric_label = factor(metric_label, levels = c("Effect r", "Std. RMSE", "Sign accuracy")),
    engine = factor(engine, levels = engine_order)
  )

truth_predictions <- safe_read_csv(paths$truth_predictions) %>%
  mutate(
    site_id_norm = normalize_site_id(site_id),
    species_norm = normalize_species(species),
    truth_probability_truth = to_num(truth_probability)
  ) %>%
  select(site_id_norm, species_norm, truth_probability_truth)

observed_occurrence <- read_occurrence_observed(paths$y_occurrence)

predictions_raw <- safe_read_csv(paths$predictions) %>%
  mutate(
    engine = factor(engine, levels = engine_order),
    site_id_norm = normalize_site_id(site_id),
    species_norm = normalize_species(species),
    observed = to_num(observed),
    truth_probability_existing = to_num(truth_probability),
    predicted_raw = to_num(predicted)
  ) %>%
  left_join(truth_predictions, by = c("site_id_norm", "species_norm")) %>%
  left_join(observed_occurrence, by = c("site_id_norm", "species_norm")) %>%
  mutate(
    observed_used = coalesce(observed, observed_input),
    truth_probability_used = coalesce(truth_probability_existing, truth_probability_truth)
  )

prediction_enriched <- predictions_raw %>%
  group_split(engine) %>%
  lapply(function(df) {
    engine_name <- as.character(df$engine[1] %||% "")
    link_hint <- if (engine_name %in% c("Hmsc-HPC", "sjSDM")) "probit" else "logit"
    transformed <- normalize_probability(df$predicted_raw, link_hint = link_hint)
    df$predicted_probability_used <- transformed$value
    df$prediction_transform <- transformed$transform
    df
  }) %>%
  bind_rows()

prediction_availability <- prediction_enriched %>%
  group_by(engine) %>%
  summarise(
    n_rows = n(),
    n_finite_prediction = sum(is.finite(predicted_probability_used)),
    n_binary_observed = sum(is.finite(observed_used) & observed_used %in% c(0, 1)),
    n_truth_probability = sum(is.finite(truth_probability_used)),
    transforms = paste(sort(unique(prediction_transform)), collapse = "; "),
    .groups = "drop"
  ) %>%
  right_join(engine_frame, by = "engine") %>%
  mutate(
    across(c(n_rows, n_finite_prediction, n_binary_observed, n_truth_probability), ~ replace_na(., 0L)),
    transforms = replace_na(transforms, "not_available"),
    occurrence_metric_status = case_when(
      n_finite_prediction >= 3 & n_binary_observed >= 3 ~ "available",
      n_rows == 0 ~ "not_available_no_standard_prediction_rows",
      n_finite_prediction == 0 ~ "not_available_no_finite_predicted_probability",
      n_binary_observed < 3 ~ "not_available_no_binary_observed_values",
      TRUE ~ "not_available"
    ),
    note = case_when(
      occurrence_metric_status == "available" ~ "Occurrence prediction metrics computed after site/species ID matching.",
      engine == "spOccupancy" ~ "Current standard prediction table did not contain finite species-level predicted occupancy probabilities.",
      TRUE ~ "Prediction output was present but could not be matched to finite observed occurrence values."
    )
  ) %>%
  arrange(engine)

prediction_metrics <- prediction_enriched %>%
  group_by(engine) %>%
  summarise(
    n_occurrence_eval = sum(is.finite(predicted_probability_used) & is.finite(observed_used) & observed_used %in% c(0, 1)),
    auc_occurrence = manual_auc(observed_used, predicted_probability_used),
    brier_occurrence = ifelse(n_occurrence_eval >= 3, mean((predicted_probability_used - observed_used)^2, na.rm = TRUE), NA_real_),
    probability_rmse_to_truth = ifelse(
      sum(is.finite(predicted_probability_used) & is.finite(truth_probability_used)) >= 3,
      sqrt(mean((predicted_probability_used - truth_probability_used)^2, na.rm = TRUE)),
      NA_real_
    ),
    truth_probability_correlation = ifelse(
      sum(is.finite(predicted_probability_used) & is.finite(truth_probability_used)) >= 3,
      suppressWarnings(cor(predicted_probability_used, truth_probability_used, use = "complete.obs")),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  right_join(engine_frame, by = "engine") %>%
  left_join(prediction_availability %>% select(engine, occurrence_metric_status, transforms, note), by = "engine") %>%
  mutate(
    n_occurrence_eval = replace_na(n_occurrence_eval, 0L),
    occurrence_metric_status = replace_na(occurrence_metric_status, "not_available"),
    transforms = replace_na(transforms, "not_available"),
    note = replace_na(note, "No occurrence prediction metrics could be computed.")
  ) %>%
  arrange(engine)

prediction_metrics_long <- prediction_metrics %>%
  select(engine, n_occurrence_eval, occurrence_metric_status, auc_occurrence, brier_occurrence, probability_rmse_to_truth) %>%
  pivot_longer(
    cols = c(auc_occurrence, brier_occurrence, probability_rmse_to_truth),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    response_family = "Occurrence",
    metric_label = recode(
      metric,
      auc_occurrence = "AUC",
      brier_occurrence = "Brier",
      probability_rmse_to_truth = "Prob. RMSE"
    ),
    metric_label = factor(metric_label, levels = c("AUC", "Brier", "Prob. RMSE")),
    engine = factor(engine, levels = engine_order)
  )

prediction_family_coverage <- tibble(
  response_family = c("Occurrence", "Count", "Continuous"),
  figure4_metric = c("AUC, Brier, probability RMSE", "RMSE, MAE pending", "RMSE, R2 pending"),
  status_in_current_run = c("plotted", "not_run_in_current_universal_benchmark", "not_run_in_current_universal_benchmark"),
  note = c(
    "Computed from occurrence-compatible unified prediction tables.",
    "Universal Benchmark has count input files, but this master run did not include count-family prediction recovery tables.",
    "Universal Benchmark has continuous input files, but this master run did not include continuous-family prediction recovery tables."
  )
)

truth_associations <- safe_read_csv(paths$truth_associations) %>%
  mutate(
    species_i_norm = normalize_species(species_i),
    species_j_norm = normalize_species(species_j),
    pair_i = pmin(species_i_norm, species_j_norm),
    pair_j = pmax(species_i_norm, species_j_norm),
    true_association = to_num(true_association)
  ) %>%
  select(pair_i, pair_j, true_association) %>%
  distinct()

associations_raw <- safe_read_csv(paths$associations) %>%
  mutate(
    engine = factor(engine, levels = engine_order),
    association_type = ifelse(
      as.character(engine) == "Hmsc-HPC" & grepl("^lambda_correlation", association_type),
      sub("^lambda_correlation", "derived_association_from_lambda", association_type),
      association_type
    ),
    scale = ifelse(
      as.character(engine) == "Hmsc-HPC" & grepl("^derived_association_from_lambda", association_type),
      "Hmsc-HPC_derived_species_association_scale",
      scale
    ),
    species_i_norm = normalize_species(species_i),
    species_j_norm = normalize_species(species_j),
    pair_i = pmin(species_i_norm, species_j_norm),
    pair_j = pmax(species_i_norm, species_j_norm),
    estimate = to_num(estimate)
  )

association_joined <- associations_raw %>%
  filter(association_type != "not_available", is.finite(estimate), species_i_norm != species_j_norm) %>%
  group_by(engine, association_type, scale, pair_i, pair_j) %>%
  summarise(estimate = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
  inner_join(truth_associations, by = c("pair_i", "pair_j")) %>%
  filter(is.finite(true_association), is.finite(estimate)) %>%
  mutate(sign_agrees = sign(estimate) == sign(true_association))

association_label_map <- c(
  "Omega" = "Hmsc Omega",
  "lambda_correlation:plot_iid" = "HPC derived assoc.",
  "derived_association_from_lambda:plot_iid" = "HPC derived assoc.",
  "residual_cor" = "jSDM resid.",
  "corMu" = "GJAM corMu",
  "covariance" = "sjSDM cov.",
  "boral_residual_correlation" = "boral resid."
)

association_metrics <- association_joined %>%
  group_by(engine, association_type, scale) %>%
  summarise(
    n_pairs = n(),
    association_structure_correlation = ifelse(n() >= 3, suppressWarnings(cor(estimate, true_association)), NA_real_),
    sign_accuracy = mean(sign_agrees, na.rm = TRUE),
    association_mae = mean(abs(estimate - true_association), na.rm = TRUE),
    comparable = "structure_only",
    note = "Lower-triangle species-pair structure only; raw association parameters are not interchangeable across engines.",
    .groups = "drop"
  ) %>%
  mutate(
    engine = factor(engine, levels = engine_order),
    association_label = recode(association_type, !!!association_label_map, .default = association_type)
  )

association_availability <- association_metrics %>%
  group_by(engine) %>%
  summarise(
    n_association_sets = n(),
    n_pairs_total = sum(n_pairs, na.rm = TRUE),
    association_types = paste(association_label, collapse = "; "),
    association_metric_status = ifelse(any(is.finite(association_structure_correlation)), "available", "not_available"),
    note = paste(unique(note), collapse = " "),
    .groups = "drop"
  ) %>%
  right_join(engine_frame, by = "engine") %>%
  mutate(
    n_association_sets = replace_na(n_association_sets, 0L),
    n_pairs_total = replace_na(n_pairs_total, 0L),
    association_types = replace_na(association_types, "none"),
    association_metric_status = replace_na(association_metric_status, "not_available"),
    note = replace_na(note, "No finite matched species-by-species association summary was available for this engine.")
  ) %>%
  arrange(engine)

association_metrics_long <- association_metrics %>%
  select(engine, association_label, association_type, scale, n_pairs, association_structure_correlation, sign_accuracy, comparable, note) %>%
  pivot_longer(
    cols = c(association_structure_correlation, sign_accuracy),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_label = recode(
      metric,
      association_structure_correlation = "Structure r",
      sign_accuracy = "Sign accuracy"
    ),
    metric_label = factor(metric_label, levels = c("Structure r", "Sign accuracy")),
    engine = factor(engine, levels = engine_order)
  ) %>%
  right_join(expand_grid(
    engine = factor(engine_order, levels = engine_order),
    metric_label = factor(c("Structure r", "Sign accuracy"), levels = c("Structure r", "Sign accuracy"))
  ), by = c("engine", "metric_label")) %>%
  left_join(association_availability %>% select(engine, association_metric_status, availability_note = note), by = "engine")

runtime_metrics <- safe_read_csv(paths$status) %>%
  mutate(
    engine = factor(engine, levels = engine_order),
    runtime_seconds = to_num(runtime_seconds),
    file_count = to_num(file_count),
    zip_size_kb = to_num(zip_size_kb),
    runtime_context = "single_small_sequential_desktop_run",
    note = "Runtime is shown only for the current small sequential desktop benchmark; it is not a scaling claim."
  ) %>%
  right_join(engine_frame, by = "engine") %>%
  arrange(engine)

scaling_path <- file.path(source_data_dir, "figure4_runtime_scaling_metrics.csv")
runtime_scaling <- if (file.exists(scaling_path)) {
  safe_read_csv(scaling_path) %>%
    mutate(
      engine = factor(engine, levels = engine_order),
      runtime_seconds = to_num(runtime_seconds),
      n_sites = to_num(n_sites),
      n_species = to_num(n_species),
      size_order = to_num(size_order),
      status = as.character(status),
      size_label = as.character(size_label)
    ) %>%
    filter(is.finite(runtime_seconds), runtime_seconds > 0, !is.na(engine), nzchar(size_label))
} else {
  tibble()
}

runtime_scaling_summary <- if (nrow(runtime_scaling)) {
  runtime_scaling %>%
    group_by(engine, size_label, size_order, n_sites, n_species) %>%
    summarise(
      median_runtime = median(runtime_seconds, na.rm = TRUE),
      runtime_min = min(runtime_seconds, na.rm = TRUE),
      runtime_max = max(runtime_seconds, na.rm = TRUE),
      n_runs = n(),
      fitted_runs = sum(status == "fitted", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(size_order, engine) %>%
    mutate(size_label = factor(size_label, levels = unique(size_label[order(size_order)])))
} else {
  tibble()
}

readr::write_csv(effects_joined, file.path(source_data_dir, "figure4_effect_recovery_joined_points.csv"))
readr::write_csv(effects_joined, file.path(output_root, "figure4_environment_effect_recovery_points.csv"))
readr::write_csv(effect_metrics, file.path(output_root, "figure4_environment_effect_recovery_metrics.csv"))
readr::write_csv(effect_metrics, file.path(source_data_dir, "figure4_environment_effect_recovery_metrics.csv"))
readr::write_csv(effect_metrics_long, file.path(output_root, "figure4_environment_effect_recovery_metrics_long.csv"))
readr::write_csv(effect_metrics_long, file.path(source_data_dir, "figure4_environment_effect_recovery_metrics_long.csv"))
readr::write_csv(prediction_enriched, file.path(source_data_dir, "figure4_prediction_recovery_joined_points.csv"))
readr::write_csv(prediction_enriched, file.path(output_root, "figure4_prediction_recovery_points.csv"))
readr::write_csv(prediction_metrics, file.path(output_root, "figure4_prediction_recovery_metrics.csv"))
readr::write_csv(prediction_metrics, file.path(source_data_dir, "figure4_prediction_recovery_metrics.csv"))
readr::write_csv(prediction_metrics_long, file.path(output_root, "figure4_prediction_recovery_metrics_long.csv"))
readr::write_csv(prediction_metrics_long, file.path(source_data_dir, "figure4_prediction_recovery_metrics_long.csv"))
readr::write_csv(prediction_availability, file.path(output_root, "figure4_prediction_availability.csv"))
readr::write_csv(prediction_availability, file.path(source_data_dir, "figure4_prediction_availability.csv"))
readr::write_csv(prediction_family_coverage, file.path(output_root, "figure4_prediction_family_coverage.csv"))
readr::write_csv(prediction_family_coverage, file.path(source_data_dir, "figure4_prediction_family_coverage.csv"))
readr::write_csv(association_joined, file.path(source_data_dir, "figure4_association_recovery_joined_pairs.csv"))
readr::write_csv(association_joined, file.path(output_root, "figure4_association_recovery_points.csv"))
readr::write_csv(association_metrics, file.path(output_root, "figure4_association_recovery_metrics.csv"))
readr::write_csv(association_metrics, file.path(source_data_dir, "figure4_association_recovery_metrics.csv"))
readr::write_csv(association_metrics_long, file.path(output_root, "figure4_association_recovery_metrics_long.csv"))
readr::write_csv(association_metrics_long, file.path(source_data_dir, "figure4_association_recovery_metrics_long.csv"))
readr::write_csv(association_availability, file.path(output_root, "figure4_association_availability.csv"))
readr::write_csv(association_availability, file.path(source_data_dir, "figure4_association_availability.csv"))
readr::write_csv(runtime_metrics, file.path(output_root, "figure4_runtime_metrics.csv"))
readr::write_csv(runtime_metrics, file.path(source_data_dir, "figure4_runtime_metrics.csv"))
readr::write_csv(runtime_scaling, file.path(output_root, "figure4_runtime_scaling_metrics.csv"))
readr::write_csv(runtime_scaling_summary, file.path(output_root, "figure4_runtime_scaling_summary.csv"))
readr::write_csv(runtime_scaling_summary, file.path(source_data_dir, "figure4_runtime_scaling_summary.csv"))

read_optional_csv <- function(path) {
  if (!file.exists(path)) return(tibble())
  safe_read_csv(path)
}

stat_effect_scaling <- read_optional_csv(file.path(source_data_dir, "figure4_statistical_scaling_effect_metrics.csv"))
stat_effect_scaling_long <- if (nrow(stat_effect_scaling)) {
  stat_effect_scaling %>%
    mutate(
      engine = factor(engine, levels = engine_order),
      size_order = to_num(size_order),
      size_label = as.character(size_label),
      effect_correlation = to_num(effect_correlation),
      standardized_rmse = to_num(standardized_rmse),
      sign_accuracy = to_num(sign_accuracy)
    ) %>%
    select(engine, size_order, size_label, n_sites, n_species, replicate,
           effect_correlation, standardized_rmse, sign_accuracy) %>%
    pivot_longer(cols = c(effect_correlation, standardized_rmse, sign_accuracy),
                 names_to = "metric", values_to = "value") %>%
    group_by(engine, size_order, size_label, metric) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      value_min = min(value, na.rm = TRUE),
      value_max = max(value, na.rm = TRUE),
      n_runs = sum(is.finite(value)),
      .groups = "drop"
    ) %>%
    mutate(
      metric_label = recode(metric,
                            effect_correlation = "Effect r",
                            standardized_rmse = "Std. RMSE",
                            sign_accuracy = "Sign accuracy"),
      metric_label = factor(metric_label, levels = c("Effect r", "Std. RMSE", "Sign accuracy")),
      size_label = factor(size_label, levels = unique(size_label[order(size_order)]))
    )
} else {
  tibble()
}

stat_prediction_scaling <- read_optional_csv(file.path(source_data_dir, "figure4_statistical_scaling_prediction_metrics.csv"))
stat_prediction_scaling_long <- if (nrow(stat_prediction_scaling)) {
  stat_prediction_scaling %>%
    mutate(
      engine = factor(engine, levels = engine_order),
      size_order = to_num(size_order),
      size_label = as.character(size_label),
      auc_occurrence = to_num(auc_occurrence),
      brier_occurrence = to_num(brier_occurrence),
      probability_rmse_to_truth = to_num(probability_rmse_to_truth)
    ) %>%
    select(engine, size_order, size_label, n_sites, n_species, replicate,
           auc_occurrence, brier_occurrence, probability_rmse_to_truth) %>%
    pivot_longer(cols = c(auc_occurrence, brier_occurrence, probability_rmse_to_truth),
                 names_to = "metric", values_to = "value") %>%
    group_by(engine, size_order, size_label, metric) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      value_min = min(value, na.rm = TRUE),
      value_max = max(value, na.rm = TRUE),
      n_runs = sum(is.finite(value)),
      .groups = "drop"
    ) %>%
    mutate(
      metric_label = recode(metric,
                            auc_occurrence = "AUC",
                            brier_occurrence = "Brier",
                            probability_rmse_to_truth = "Prob. RMSE"),
      metric_label = factor(metric_label, levels = c("AUC", "Brier", "Prob. RMSE")),
      size_label = factor(size_label, levels = unique(size_label[order(size_order)]))
    )
} else {
  tibble()
}

stat_association_scaling <- read_optional_csv(file.path(source_data_dir, "figure4_statistical_scaling_association_metrics.csv"))
if (nrow(stat_association_scaling) && !("association_label" %in% names(stat_association_scaling))) {
  stat_association_scaling$association_label <- if ("association_type" %in% names(stat_association_scaling)) {
    stat_association_scaling$association_type
  } else {
    NA_character_
  }
}
stat_association_scaling_long <- if (nrow(stat_association_scaling)) {
  stat_association_scaling %>%
    mutate(
      engine = factor(engine, levels = engine_order),
      size_order = to_num(size_order),
      size_label = as.character(size_label),
      association_structure_correlation = to_num(association_structure_correlation),
      sign_accuracy = to_num(sign_accuracy),
      association_label = as.character(association_label %||% association_type)
    ) %>%
    select(engine, association_label, association_type, size_order, size_label,
           n_sites, n_species, replicate, association_structure_correlation, sign_accuracy) %>%
    pivot_longer(cols = c(association_structure_correlation, sign_accuracy),
                 names_to = "metric", values_to = "value") %>%
    group_by(engine, association_label, association_type, size_order, size_label, metric) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      value_min = min(value, na.rm = TRUE),
      value_max = max(value, na.rm = TRUE),
      n_runs = sum(is.finite(value)),
      .groups = "drop"
    ) %>%
    mutate(
      metric_label = recode(metric,
                            association_structure_correlation = "Structure r",
                            sign_accuracy = "Sign accuracy"),
      metric_label = factor(metric_label, levels = c("Structure r", "Sign accuracy")),
      size_label = factor(size_label, levels = unique(size_label[order(size_order)]))
    ) %>%
    filter(is.finite(median_value))
} else {
  tibble()
}

theme_publication <- function(base_size = 17) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size * 1.08, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size * 0.78, color = "grey30", margin = margin(b = 8)),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = base_size * 0.72),
      axis.text.y = element_text(size = base_size * 0.78),
      strip.background = element_rect(fill = "grey94", color = "grey78", linewidth = 0.4),
      strip.text = element_text(face = "bold", size = base_size * 0.78),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
      legend.position = "none",
      plot.margin = margin(8, 12, 8, 8)
    )
}

ideal_effect <- tibble(
  metric_label = factor(c("Effect r", "Std. RMSE", "Sign accuracy"), levels = c("Effect r", "Std. RMSE", "Sign accuracy")),
  ideal = c(1, 0, 1)
)

p_a <- if (nrow(stat_effect_scaling_long)) {
  ggplot(stat_effect_scaling_long,
         aes(x = size_label, y = median_value, color = engine, group = engine)) +
    geom_hline(data = ideal_effect, aes(yintercept = ideal), inherit.aes = FALSE,
               linetype = "dashed", color = "grey55", linewidth = 0.55) +
    geom_line(linewidth = 0.75, alpha = 0.9, na.rm = TRUE) +
    geom_errorbar(aes(ymin = value_min, ymax = value_max), width = 0.08,
                  linewidth = 0.45, alpha = 0.6, na.rm = TRUE) +
    geom_point(size = 3.0, na.rm = TRUE) +
    facet_wrap(~ metric_label, scales = "free_y", nrow = 1) +
    scale_color_manual(values = engine_palette, drop = FALSE) +
    labs(
      title = "Environmental-effect recovery",
      subtitle = "Continuous predictors; categorical contrasts are exported separately",
      x = NULL,
      y = "Metric value"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 22, hjust = 1))
} else {
  ggplot(effect_metrics_long, aes(x = engine, y = value, color = engine)) +
    geom_hline(data = ideal_effect, aes(yintercept = ideal), inherit.aes = FALSE, linetype = "dashed", color = "grey55", linewidth = 0.55) +
    geom_point(size = 3.8, na.rm = TRUE) +
    geom_text(
      data = effect_metrics_long %>% filter(!is.finite(value)),
      aes(y = 0, label = "NA"),
      color = "grey45", size = 4.2, fontface = "bold"
    ) +
    facet_wrap(~ metric_label, scales = "free_y", nrow = 1) +
    scale_color_manual(values = engine_palette, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = "Environmental-effect recovery",
      subtitle = "Continuous predictors; categorical contrasts are exported separately",
      x = NULL,
      y = "Metric value"
    ) +
    theme_publication()
}

ideal_prediction <- tibble(
  metric_label = factor(c("AUC", "Brier", "Prob. RMSE"), levels = c("AUC", "Brier", "Prob. RMSE")),
  ideal = c(1, 0, 0)
)

p_b <- if (nrow(stat_prediction_scaling_long)) {
  ggplot(stat_prediction_scaling_long,
         aes(x = size_label, y = median_value, color = engine, group = engine)) +
    geom_hline(data = ideal_prediction, aes(yintercept = ideal), inherit.aes = FALSE,
               linetype = "dashed", color = "grey55", linewidth = 0.55) +
    geom_line(linewidth = 0.75, alpha = 0.9, na.rm = TRUE) +
    geom_errorbar(aes(ymin = value_min, ymax = value_max), width = 0.08,
                  linewidth = 0.45, alpha = 0.6, na.rm = TRUE) +
    geom_point(size = 3.0, na.rm = TRUE) +
    facet_wrap(~ metric_label, scales = "free_y", nrow = 1) +
    scale_color_manual(values = engine_palette, drop = FALSE) +
    labs(
      title = "Prediction recovery",
      subtitle = "Occurrence prediction across benchmark sizes",
      x = NULL,
      y = "Metric value"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 22, hjust = 1))
} else {
  ggplot(prediction_metrics_long, aes(x = engine, y = value, color = engine)) +
    geom_hline(data = ideal_prediction, aes(yintercept = ideal), inherit.aes = FALSE, linetype = "dashed", color = "grey55", linewidth = 0.55) +
    geom_point(size = 3.8, na.rm = TRUE) +
    geom_text(
      data = prediction_metrics_long %>% filter(!is.finite(value)),
      aes(y = 0, label = "NA"),
      color = "grey45", size = 4.2, fontface = "bold"
    ) +
    facet_wrap(~ metric_label, scales = "free_y", nrow = 1) +
    scale_color_manual(values = engine_palette, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = "Prediction recovery",
      subtitle = "Occurrence metrics from current run;\ncount and continuous checks are recorded separately",
      x = NULL,
      y = "Metric value"
    ) +
    theme_publication()
}

ideal_association <- tibble(
  metric_label = factor(c("Structure r", "Sign accuracy"), levels = c("Structure r", "Sign accuracy")),
  ideal = c(1, 1)
)

p_c <- if (nrow(stat_association_scaling_long)) {
  ggplot(
    stat_association_scaling_long,
    aes(
      x = size_label,
      y = median_value,
      color = engine,
      shape = association_label,
      group = interaction(engine, association_label)
    )
  ) +
    geom_hline(data = ideal_association, aes(yintercept = ideal), inherit.aes = FALSE,
               linetype = "dashed", color = "grey55", linewidth = 0.55) +
    geom_line(linewidth = 0.75, alpha = 0.85, na.rm = TRUE) +
    geom_errorbar(aes(ymin = value_min, ymax = value_max), width = 0.08,
                  linewidth = 0.45, alpha = 0.6, na.rm = TRUE) +
    geom_point(size = 3.1, stroke = 1.0, na.rm = TRUE) +
    facet_wrap(~ metric_label, scales = "free_y", nrow = 1) +
    scale_color_manual(values = engine_palette, drop = FALSE) +
    scale_shape_manual(
      values = c(
        "Hmsc Omega" = 15,
        "HPC derived assoc." = 17,
        "jSDM resid." = 3,
        "GJAM corMu" = 16,
        "sjSDM cov." = 7,
        "boral resid." = 18
      ),
      drop = FALSE,
      na.translate = FALSE
    ) +
    labs(
      title = "Association-structure recovery",
      subtitle = "Lower-triangle structure across benchmark sizes;\nestimands remain engine-specific",
      x = NULL,
      y = "Metric value"
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 22, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      legend.key.width = unit(1.2, "lines"),
      legend.margin = margin(t = 2, r = 4, b = 2, l = 4),
      legend.box.margin = margin(2, 2, 2, 2)
    ) +
    guides(
      color = "none",
      shape = guide_legend(title = NULL, nrow = 1, byrow = TRUE)
    )
} else {
  ggplot(association_metrics_long, aes(x = engine, y = value, color = engine, shape = association_label)) +
    geom_hline(data = ideal_association, aes(yintercept = ideal), inherit.aes = FALSE, linetype = "dashed", color = "grey55", linewidth = 0.55) +
    geom_point(size = 4.1, stroke = 1.1, na.rm = TRUE) +
    geom_text(
      data = association_metrics_long %>% filter(!is.finite(value)),
      aes(y = 0, label = "NA"),
      color = "grey45", size = 4.2, fontface = "bold"
    ) +
    facet_wrap(~ metric_label, scales = "free_y", nrow = 1) +
    scale_color_manual(values = engine_palette, drop = FALSE) +
    scale_shape_manual(
      values = c(
        "Hmsc Omega" = 15,
        "HPC derived assoc." = 17,
        "jSDM resid." = 3,
        "GJAM corMu" = 16,
        "sjSDM cov." = 7,
        "boral resid." = 18
      ),
      drop = FALSE,
      na.translate = FALSE
    ) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = "Association-structure recovery",
      subtitle = "Lower-triangle structure only;\nraw association estimands are not interchangeable",
      x = NULL,
      y = "Metric value"
    ) +
    theme_publication() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      legend.key.width = unit(1.2, "lines"),
      legend.margin = margin(t = 2, r = 4, b = 2, l = 4),
      legend.box.margin = margin(2, 2, 2, 2)
    ) +
    guides(
      color = "none",
      shape = guide_legend(title = NULL, nrow = 1, byrow = TRUE)
    )
}

p_d <- if (nrow(runtime_scaling_summary)) {
  ggplot(runtime_scaling_summary,
         aes(x = size_label, y = median_runtime, color = engine, group = engine)) +
    geom_line(linewidth = 0.85, alpha = 0.9, na.rm = TRUE) +
    geom_errorbar(aes(ymin = runtime_min, ymax = runtime_max), width = 0.10,
                  linewidth = 0.55, alpha = 0.75, na.rm = TRUE) +
    geom_point(size = 3.4, na.rm = TRUE) +
    scale_color_manual(values = engine_palette, drop = FALSE) +
    scale_y_log10(labels = label_number(accuracy = 0.1)) +
    labs(
      title = "Computational scaling",
      subtitle = "Sequential benchmark sizes;\npoints are medians and bars are run ranges",
      x = "Benchmark size",
      y = "Seconds, log scale"
    ) +
    theme_publication() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(angle = 18, hjust = 1)
    )
} else {
  ggplot(runtime_metrics, aes(x = engine, y = runtime_seconds, color = engine)) +
    geom_point(size = 4.2, na.rm = TRUE) +
    scale_color_manual(values = engine_palette, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_y_log10(breaks = c(0.5, 1, 3, 10, 30), labels = label_number(accuracy = 0.1)) +
    labs(
      title = "Runtime",
      subtitle = "Single small sequential desktop run;\nuse scaling runs for speed claims",
      x = NULL,
      y = "Seconds, log scale"
    ) +
    theme_publication()
}

figure4 <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    title = "Statistical Universal Benchmark",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 26, hjust = 0, margin = margin(b = 8)),
      plot.tag = element_text(face = "bold", size = 22),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

ggsave(file.path(output_root, "Figure4_Statistical_Universal_Benchmark.pdf"), figure4, width = 16.2, height = 10.8, units = "in", device = cairo_pdf)
ggsave(file.path(output_root, "Figure4_Statistical_Universal_Benchmark.png"), figure4, width = 16.2, height = 10.8, units = "in", dpi = 360)
ggsave(file.path(output_root, "Figure4_Statistical_Universal_Benchmark.tiff"), figure4, width = 16.2, height = 10.8, units = "in", dpi = 600, compression = "lzw")
ggsave(file.path(output_root, "Figure4_Statistical_Universal_Benchmark.svg"), figure4, width = 16.2, height = 10.8, units = "in")

runtime_caption <- if (nrow(runtime_scaling_summary)) {
  "D, Computational scaling summarizes sequential small, medium, large and ultra-large benchmark runs. Points show median elapsed runtime and bars show the observed range across repeated runs. These timings are software-environment specific and should be interpreted together with package versions, sampler settings and hardware."
} else {
  "D, Runtime shows elapsed seconds for a single small sequential desktop benchmark. It is useful as a reproducibility record, but not as a general scaling or speed ranking. Scaling claims require small, medium, large and ultra-large benchmark sizes with repeated runs."
}

effect_caption <- if (nrow(stat_effect_scaling_long)) {
  "A, Environmental-effect recovery is summarized across small, medium, large and ultra-large latent-truth benchmarks by effect correlation, standardized RMSE and sign accuracy after matching engine-specific estimates to the corresponding simulated species-by-continuous-predictor effects. Categorical substrate contrasts are exported in the standard tables but excluded from the main cross-engine recovery score because factor baselines differ across engines. Trends with increasing sample size are emphasized over a single small-run score."
} else {
  "A, Environmental-effect recovery is summarized by effect correlation, standardized RMSE and sign accuracy after matching engine-specific estimates to the corresponding simulated species-by-continuous-predictor effects. Categorical substrate contrasts are exported in the standard tables but excluded from the main cross-engine recovery score because factor baselines differ across engines. Correlation and sign accuracy are the safest cross-engine summaries; standardized RMSE should be interpreted only after compatible scaling."
}

prediction_caption <- if (nrow(stat_prediction_scaling_long)) {
  "B, Prediction recovery is summarized across benchmark sizes using occurrence metrics: AUC, Brier score and RMSE to the simulated occurrence probability. Lower Brier and probability RMSE are better; higher AUC is better. Count and continuous response-family checks are recorded separately when exported by the engine-specific benchmark runs."
} else {
  "B, Prediction recovery uses occurrence metrics from the current Universal Benchmark run: AUC, Brier score and RMSE to the simulated occurrence probability. Lower Brier and probability RMSE are better; higher AUC is better. Count and continuous response-family checks are listed in `figure4_prediction_family_coverage.csv` because this master run did not produce family-specific count or continuous prediction recovery tables."
}

association_caption <- if (nrow(stat_association_scaling_long)) {
  "C, Association recovery is restricted to species-by-species structure summaries across benchmark sizes. Hmsc-HPC is labelled as a derived association summary, not raw Lambda. These values compare lower-triangle structure with the simulated association matrix; they are not evidence that Hmsc Omega, GJAM corMu, sjSDM covariance, jSDM residual correlation and boral residual correlation are the same estimand."
} else {
  "C, Association recovery is restricted to species-by-species structure summaries. Hmsc-HPC is labelled as a derived association summary, not raw Lambda. These values compare lower-triangle structure with the simulated association matrix; they are not evidence that Hmsc Omega, GJAM corMu, sjSDM covariance, jSDM residual correlation and boral residual correlation are the same estimand. spOccupancy is retained as not available when no finite comparable residual association matrix is exported in this run."
}

caption <- c(
  "# Figure 4. Statistical Universal Benchmark",
  "",
  effect_caption,
  "",
  prediction_caption,
  "",
  association_caption,
  "",
  runtime_caption,
  "",
  "Status note: fitted means the engine produced fitted outputs and unified tables. It does not prove convergence, identifiability or ecological adequacy."
)
writeLines(caption, file.path(output_root, "Figure4_caption_draft.md"), useBytes = TRUE)

source_note <- c(
  "Figure 4 source note",
  paste0("Universal benchmark run: ", benchmark_run),
  "Source tables are read from the selected universal_benchmark_* master directory under output/.",
  "All plotted source tables are exported beside this figure and under source_data/.",
  "The benchmark binary occurrence truth is generated on a probit latent scale, then engine-specific predictions are mapped back to occurrence probabilities where possible.",
  "Hmsc-HPC association is treated as a derived species-by-species association summary; raw Lambda is not plotted as an association estimand.",
  "sjSDM standard predictions are exported with predict(type='link'), which is the probability scale in sjSDM 1.0.7; Hmsc-HPC probit-scale predictions are converted with pnorm when the exported values are outside [0, 1], and other link-scale values use the logit fallback unless the standard table records a probability scale.",
  "spOccupancy prediction or association values are shown as NA when the standard output lacks finite comparable species-level quantities.",
  if (nrow(stat_effect_scaling_long)) "Panel A uses statistical scaling source data from source_data/figure4_statistical_scaling_effect_metrics.csv." else "Panel A uses the selected benchmark run only.",
  "Panel A excludes categorical substrate contrasts from the plotted recovery metrics; complete engine-specific categorical effects remain in standard/effects_species_environment.csv and the truth contrast table.",
  if (nrow(stat_prediction_scaling_long)) "Panel B uses statistical scaling source data from source_data/figure4_statistical_scaling_prediction_metrics.csv." else "Panel B uses the selected benchmark run only.",
  if (nrow(stat_association_scaling_long)) "Panel C uses statistical scaling source data from source_data/figure4_statistical_scaling_association_metrics.csv." else "Panel C uses the selected benchmark run only.",
  "Count and continuous family metrics require separate family-specific benchmark outputs unless exported by the statistical scaling run.",
  if (nrow(runtime_scaling_summary)) paste0("Runtime scaling source table: ", scaling_path) else "No runtime scaling source table was found; panel D uses the selected benchmark run only."
)
writeLines(source_note, file.path(output_root, "Figure4_source_note.txt"), useBytes = TRUE)
writeLines(c("# Figure 4 source note", "", source_note[-1]), file.path(output_root, "Figure4_Statistical_Universal_Benchmark_source_note.md"), useBytes = TRUE)

manifest <- tibble(
  file = c(
    "Figure4_Statistical_Universal_Benchmark.pdf",
    "Figure4_Statistical_Universal_Benchmark.png",
    "Figure4_Statistical_Universal_Benchmark.tiff",
    "Figure4_Statistical_Universal_Benchmark.svg",
    "Figure4_caption_draft.md",
    "Figure4_source_note.txt",
    "figure4_environment_effect_recovery_metrics.csv",
    "figure4_prediction_recovery_metrics.csv",
    "figure4_prediction_availability.csv",
    "figure4_prediction_family_coverage.csv",
    "figure4_association_recovery_metrics.csv",
    "figure4_association_availability.csv",
    "figure4_runtime_metrics.csv",
    "figure4_runtime_scaling_metrics.csv",
    "figure4_runtime_scaling_summary.csv",
    "source_data/figure4_statistical_scaling_effect_metrics.csv",
    "source_data/figure4_statistical_scaling_prediction_metrics.csv",
    "source_data/figure4_statistical_scaling_association_metrics.csv",
    "source_data/figure4_statistical_scaling_status.csv",
    "source_data/figure4_statistical_scaling_run_log.csv"
  ),
  role = c(
    "vector_figure",
    "raster_preview",
    "print_raster",
    "editable_vector",
    "caption",
    "source_note",
    "panel_a_metrics",
    "panel_b_metrics",
    "panel_b_availability",
    "panel_b_family_coverage",
    "panel_c_metrics",
    "panel_c_availability",
    "panel_d_metrics",
    "panel_d_scaling_raw",
    "panel_d_scaling_summary",
    "panel_a_statistical_scaling",
    "panel_b_statistical_scaling",
    "panel_c_statistical_scaling",
    "statistical_scaling_status",
    "statistical_scaling_run_log"
  )
)
readr::write_csv(manifest, file.path(output_root, "Figure4_manifest.csv"))

message("Figure 4 regenerated in: ", output_root)
message("Benchmark run: ", benchmark_run)
