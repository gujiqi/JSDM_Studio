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
main_effect_predictors <- c("pH", "moisture", "canopy", "elevation")

engine_palette <- c(
  "Hmsc" = "#1f78b4",
  "Hmsc-HPC" = "#4f9a8b",
  "jSDM" = "#8b63ad",
  "GJAM" = "#d55e00",
  "spOccupancy" = "#5f7fa8",
  "sjSDM" = "#2ca25f",
  "boral" = "#7f7f7f"
)

size_palette <- c(
  "Built-in" = "#222222",
  "Small" = "#8c8c8c",
  "Medium" = "#5f7fa8",
  "Large" = "#2ca25f",
  "Ultra-large" = "#d55e00"
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (length(x) == 1 && is.atomic(x) && is.na(x)) return(y)
  x
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(project_root = "", scaling_source_dir = "", run_log = "", builtin_master_dir = "")
  pos <- character()
  for (a in args) {
    if (grepl("^--", a)) {
      kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
      key <- kv[[1]]
      val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "true"
      out[[key]] <- val
    } else {
      pos <- c(pos, a)
    }
  }
  if (!nzchar(out$project_root) && length(pos)) out$project_root <- pos[[1]]
  out
}

cmd <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", cmd, value = TRUE))
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg[1], winslash = "/", mustWork = TRUE))
} else {
  getwd()
}

args <- parse_args()
project_root <- if (nzchar(args$project_root)) {
  normalizePath(args$project_root, winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
}

output_root <- file.path(project_root, "docs", "paper_figures", "figure5_universal_benchmark_cross_model_comparison")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

figure4_source_data <- if (nzchar(args$scaling_source_dir)) {
  normalizePath(args$scaling_source_dir, winslash = "/", mustWork = TRUE)
} else {
  file.path(project_root, "docs", "paper_figures", "figure4_statistical_universal_benchmark", "source_data")
}
run_log_path <- if (nzchar(args$run_log)) {
  normalizePath(args$run_log, winslash = "/", mustWork = TRUE)
} else {
  file.path(figure4_source_data, "figure4_statistical_scaling_run_log.csv")
}

missing <- c(run_log_path)
missing <- missing[!file.exists(missing)]
if (length(missing)) {
  stop("Missing required Figure 5 scaling source files: ", paste(missing, collapse = ", "))
}

safe_read_csv <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

to_num <- function(x) suppressWarnings(as.numeric(x))

standardize_vec <- function(x) {
  x <- to_num(x)
  if (sum(is.finite(x)) < 2) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
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

read_occurrence_observed <- function(path) {
  if (!file.exists(path)) return(tibble(site_id_norm = character(), species_norm = character(), observed_input = numeric()))
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

format_size_label <- function(size_label, n_sites, n_species) {
  paste0(size_label, "\n", n_sites, "x", n_species)
}

read_effect_points_one <- function(row) {
  master_dir <- normalizePath(row$master_dir, winslash = "/", mustWork = TRUE)
  effects_path <- file.path(master_dir, "04_unified_standard_results", "effects_all_models_long.csv")
  truth_path <- file.path(master_dir, "02_truth", "true_environment_effects.csv")
  if (!file.exists(effects_path) || !file.exists(truth_path)) {
    warning("Skipping Figure 5 effect points for missing master files: ", master_dir)
    return(tibble())
  }

  truth_effects <- safe_read_csv(truth_path) %>%
    mutate(
      species_norm = normalize_species(species),
      predictor_norm = normalize_predictor(predictor),
      true_effect = to_num(true_effect)
    ) %>%
    select(species_norm, predictor_norm, true_effect)

  safe_read_csv(effects_path) %>%
    mutate(
      estimate_raw = to_num(estimate),
      spocc_species = extract_spocc_species(predictor),
      species_norm = ifelse(engine == "spOccupancy" & !is.na(spocc_species),
                            spocc_species, normalize_species(species)),
      predictor_norm = normalize_predictor(predictor),
      engine = factor(engine, levels = engine_order)
    ) %>%
    filter(effect_type %in% c("environment", "occurrence"),
           predictor_norm != "(Intercept)",
           predictor_norm %in% main_effect_predictors,
           is.finite(estimate_raw)) %>%
    group_by(engine, species_norm, predictor_norm, scale) %>%
    summarise(estimate_raw = mean(estimate_raw, na.rm = TRUE), .groups = "drop") %>%
    inner_join(truth_effects, by = c("species_norm", "predictor_norm")) %>%
    filter(is.finite(true_effect), is.finite(estimate_raw)) %>%
    group_by(engine) %>%
    mutate(
      true_effect_standardised = standardize_vec(true_effect),
      estimated_effect_standardised = standardize_vec(estimate_raw)
    ) %>%
    ungroup() %>%
    filter(is.finite(true_effect_standardised), is.finite(estimated_effect_standardised)) %>%
    mutate(
      size_order = as.integer(row$size_order),
      size_label = as.character(row$size_label),
      n_sites = as.integer(row$n_sites),
      n_species = as.integer(row$n_species),
      replicate = as.integer(row$replicate),
      benchmark_label = format_size_label(size_label, n_sites, n_species),
      benchmark_run = basename(dirname(master_dir)),
      master_dir = master_dir,
      .before = 1
    )
}

engine_frame <- tibble(engine = factor(engine_order, levels = engine_order))

read_config_scalar <- function(path, key) {
  if (!file.exists(path)) return(NA_character_)
  txt <- readLines(path, warn = FALSE)
  hit <- grep(paste0("^", key, ":\\s*"), txt, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub(paste0("^", key, ":\\s*"), "", hit[[1]]))
}

count_fitted_engines <- function(master_dir) {
  status_file <- file.path(master_dir, "04_unified_standard_results", "model_status_matrix.csv")
  if (!file.exists(status_file)) return(0L)
  status_tbl <- tryCatch(safe_read_csv(status_file), error = function(e) tibble())
  if (!nrow(status_tbl) || !"status" %in% names(status_tbl)) return(0L)
  sum(status_tbl$status == "fitted", na.rm = TRUE)
}

find_latest_builtin_master <- function(project_root) {
  explicit_default <- file.path(project_root, "output", "universal_benchmark_figure5_builtin_20260831", "master")
  if (dir.exists(explicit_default) && count_fitted_engines(explicit_default) == length(engine_order)) {
    return(normalizePath(explicit_default, winslash = "/", mustWork = TRUE))
  }

  dirs <- list.dirs(file.path(project_root, "output"), recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^universal_benchmark_", basename(dirs))]
  masters <- file.path(dirs, "master")
  cfgs <- file.path(masters, "09_configs", "universal_benchmark_config.yml")
  keep <- vapply(seq_along(masters), function(i) {
    dir.exists(masters[[i]]) &&
      identical(read_config_scalar(cfgs[[i]], "n_sites"), "120") &&
      identical(read_config_scalar(cfgs[[i]], "n_species"), "12") &&
      count_fitted_engines(masters[[i]]) == length(engine_order)
  }, logical(1))
  masters <- masters[keep]
  if (!length(masters)) return(NA_character_)
  normalizePath(masters[which.max(file.info(dirname(masters))$mtime)], winslash = "/", mustWork = TRUE)
}

read_prediction_metrics_one <- function(row) {
  master_dir <- normalizePath(row$master_dir, winslash = "/", mustWork = TRUE)
  pred_path <- file.path(master_dir, "04_unified_standard_results", "predictions_all_models_long.csv")
  truth_path <- file.path(master_dir, "02_truth", "true_predictions.csv")
  y_path <- file.path(master_dir, "01_benchmark_inputs", "Y_occurrence.csv")
  if (!file.exists(pred_path) || !file.exists(truth_path)) {
    return(engine_frame %>%
             mutate(
               n_occurrence_eval = 0L,
               auc_occurrence = NA_real_,
               brier_occurrence = NA_real_,
               probability_rmse_to_truth = NA_real_,
               truth_probability_correlation = NA_real_,
               occurrence_metric_status = "missing_prediction_source",
               transforms = "not_available",
               note = "Missing standard prediction or truth prediction table."
             ))
  }

  truth_predictions <- safe_read_csv(truth_path) %>%
    mutate(
      site_id_norm = normalize_site_id(site_id),
      species_norm = normalize_species(species),
      truth_probability_truth = to_num(truth_probability)
    ) %>%
    select(site_id_norm, species_norm, truth_probability_truth)

  observed_occurrence <- read_occurrence_observed(y_path)

  predictions_raw <- safe_read_csv(pred_path) %>%
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
        TRUE ~ "Prediction output was present but could not be matched to finite observed occurrence values."
      )
    )

  prediction_enriched %>%
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
      size_order = as.integer(row$size_order),
      size_label = as.character(row$size_label),
      n_sites = as.integer(row$n_sites),
      n_species = as.integer(row$n_species),
      replicate = as.integer(row$replicate),
      benchmark_label = format_size_label(size_label, n_sites, n_species),
      benchmark_run = basename(dirname(master_dir)),
      master_dir = master_dir,
      .before = 1
    )
}

read_status_one <- function(row) {
  master_dir <- normalizePath(row$master_dir, winslash = "/", mustWork = TRUE)
  status_file <- file.path(master_dir, "04_unified_standard_results", "model_status_matrix.csv")
  if (!file.exists(status_file)) {
    return(engine_frame %>%
             mutate(status = "missing_status_file", file_count = NA_integer_, missing_count = NA_integer_))
  }
  safe_read_csv(status_file) %>%
    mutate(
      engine = factor(engine, levels = engine_order),
      size_order = as.integer(row$size_order),
      size_label = as.character(row$size_label),
      n_sites = as.integer(row$n_sites),
      n_species = as.integer(row$n_species),
      replicate = as.integer(row$replicate),
      benchmark_label = format_size_label(size_label, n_sites, n_species),
      benchmark_run = basename(dirname(master_dir)),
      master_dir = master_dir,
      .before = 1
    )
}

run_log <- safe_read_csv(run_log_path) %>%
  arrange(size_order, replicate) %>%
  distinct(size_order, size_label, n_sites, n_species, replicate, master_dir, .keep_all = TRUE)

builtin_master <- if (nzchar(args$builtin_master_dir)) {
  normalizePath(args$builtin_master_dir, winslash = "/", mustWork = TRUE)
} else {
  find_latest_builtin_master(project_root)
}
if (!is.na(builtin_master) && dir.exists(builtin_master) &&
    !normalizePath(builtin_master, winslash = "/", mustWork = FALSE) %in% normalizePath(run_log$master_dir, winslash = "/", mustWork = FALSE)) {
  run_log <- bind_rows(
    tibble(
      size_order = 0L,
      size_label = "Built-in",
      n_sites = 120L,
      n_species = 12L,
      replicate = 1L,
      seed = NA_integer_,
      runner_status = 0L,
      master_dir = builtin_master,
      runner_log = NA_character_,
      figure_status = NA_integer_,
      figure_log = NA_character_
    ),
    run_log
  )
}

run_log <- run_log %>%
  arrange(size_order, replicate) %>%
  mutate(size_order = as.integer(size_order))

if (!nrow(run_log)) stop("Figure 5 run log did not contain any master directories.")
if (any(!dir.exists(run_log$master_dir))) {
  stop("One or more Figure 5 master directories from the run log do not exist.")
}

declared_size_specs <- tibble(
  size_label = c("Built-in", "Small", "Medium", "Large", "Ultra-large"),
  dataset = c("Built-in example", "Figure5 Small", "Figure5 Medium", "Figure5 Large", "Figure5 Ultra-large"),
  manuscript_input_files = c(32L, 32L, 32L, 32L, 32L),
  manuscript_raw_input_size_mb = c(0.39, 0.06, 0.37, 1.52, 4.02)
)

size_table <- run_log %>%
  group_by(size_order, size_label, n_sites, n_species) %>%
  summarise(master_dir = first(master_dir), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    visits = 3L,
    prediction_sites = 8L,
    measured_master_input_files = {
      files <- c(
        list.files(file.path(master_dir, "01_benchmark_inputs"), recursive = TRUE, full.names = TRUE),
        list.files(file.path(master_dir, "02_truth"), recursive = TRUE, full.names = TRUE)
      )
      sum(file.exists(files) & !file.info(files)$isdir)
    },
    measured_master_input_size_mb = {
      files <- c(
        list.files(file.path(master_dir, "01_benchmark_inputs"), recursive = TRUE, full.names = TRUE),
        list.files(file.path(master_dir, "02_truth"), recursive = TRUE, full.names = TRUE)
      )
      files <- files[file.exists(files) & !file.info(files)$isdir]
      round(sum(file.info(files)$size, na.rm = TRUE) / 1024^2, 2)
    },
    benchmark_label = format_size_label(size_label, n_sites, n_species)
  ) %>%
  ungroup() %>%
  left_join(declared_size_specs, by = "size_label") %>%
  mutate(
    raw_input_files = coalesce(manuscript_input_files, measured_master_input_files),
    raw_input_size_mb = coalesce(manuscript_raw_input_size_mb, measured_master_input_size_mb)
  )

size_levels <- as.character(size_table$size_label[order(size_table$size_order)])
benchmark_levels <- as.character(size_table$benchmark_label[order(size_table$size_order)])
size_table <- size_table %>%
  mutate(
    benchmark_label = factor(benchmark_label, levels = benchmark_levels),
    size_label = factor(size_label, levels = size_levels)
  )

effect_points <- bind_rows(lapply(seq_len(nrow(run_log)), function(i) read_effect_points_one(run_log[i, , drop = FALSE]))) %>%
  mutate(
    engine = factor(engine, levels = engine_order),
    size_label = factor(size_label, levels = as.character(size_table$size_label)),
    benchmark_label = factor(benchmark_label, levels = as.character(size_table$benchmark_label))
  )

if (!nrow(effect_points)) stop("No finite Figure 5 effect points could be assembled from the scaling masters.")

effect_metrics <- effect_points %>%
  group_by(size_order, size_label, n_sites, n_species, benchmark_label, engine) %>%
  summarise(
    n_effects = n(),
    raw_effect_correlation = ifelse(n() >= 3, suppressWarnings(cor(true_effect, estimate_raw)), NA_real_),
    standardised_rmse = sqrt(mean((estimated_effect_standardised - true_effect_standardised)^2)),
    direction_agreement = mean(sign(estimate_raw) == sign(true_effect)),
    .groups = "drop"
  ) %>%
  arrange(size_order, engine)

prediction_metrics <- bind_rows(lapply(seq_len(nrow(run_log)), function(i) read_prediction_metrics_one(run_log[i, , drop = FALSE]))) %>%
  mutate(
    engine = factor(engine, levels = engine_order),
    size_label = factor(size_label, levels = size_levels),
    benchmark_label = factor(benchmark_label, levels = benchmark_levels),
    auc_occurrence = to_num(auc_occurrence),
    brier_occurrence = to_num(brier_occurrence),
    probability_rmse_to_truth = to_num(probability_rmse_to_truth),
    truth_probability_correlation = to_num(truth_probability_correlation)
  ) %>%
  arrange(size_order, engine)

prediction_metrics_long <- prediction_metrics %>%
  select(size_order, size_label, n_sites, n_species, benchmark_label, replicate, engine,
         n_occurrence_eval, auc_occurrence, brier_occurrence,
         probability_rmse_to_truth, truth_probability_correlation) %>%
  pivot_longer(
    cols = c(auc_occurrence, brier_occurrence, probability_rmse_to_truth),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_label = recode(
      metric,
      auc_occurrence = "AUC",
      brier_occurrence = "Brier",
      probability_rmse_to_truth = "Prob. RMSE"
    ),
    metric_label = factor(metric_label, levels = c("AUC", "Brier", "Prob. RMSE")),
    metric_direction = ifelse(metric_label == "AUC", "higher_is_better", "lower_is_better")
  ) %>%
  filter(is.finite(value))

status_by_size <- bind_rows(lapply(seq_len(nrow(run_log)), function(i) read_status_one(run_log[i, , drop = FALSE]))) %>%
  mutate(
    engine = factor(engine, levels = engine_order),
    size_label = factor(size_label, levels = size_levels),
    benchmark_label = factor(benchmark_label, levels = benchmark_levels)
  )

stale_legacy_outputs <- c(
  "figure5_panelC_output_contract_detail.csv",
  "figure5_panelC_output_completeness_summary.csv",
  "figure5_panelB_prediction_metric_availability.csv",
  "figure5_panelB_prediction_recovery_points.csv"
)
unlink(file.path(output_root, stale_legacy_outputs), force = TRUE)

readr::write_csv(size_table, file.path(output_root, "figure5_benchmark_size_table.csv"))
readr::write_csv(status_by_size, file.path(output_root, "figure5_model_status_by_size.csv"))
readr::write_csv(effect_points, file.path(output_root, "figure5_panelA_effect_recovery_points.csv"))
readr::write_csv(effect_metrics, file.path(output_root, "figure5_panelA_effect_recovery_metrics.csv"))
readr::write_csv(prediction_metrics, file.path(output_root, "figure5_panelB_prediction_recovery_metrics.csv"))
readr::write_csv(prediction_metrics_long, file.path(output_root, "figure5_panelB_prediction_recovery_metrics_long.csv"))

theme_paper <- function(base_size = 9.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 32, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_size * 0.82),
      strip.background = element_rect(fill = "grey94", color = "grey80", linewidth = 0.3),
      strip.text = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey91", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 2.0),
      plot.subtitle = element_text(size = base_size, color = "grey25"),
      plot.margin = margin(5, 8, 5, 5)
    )
}

axis_limit <- max(abs(c(effect_points$true_effect_standardised,
                        effect_points$estimated_effect_standardised)), na.rm = TRUE)
axis_limit <- max(3.2, min(5.0, ceiling(axis_limit * 10) / 10))

pA <- ggplot(effect_points, aes(true_effect_standardised, estimated_effect_standardised, color = size_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.40, color = "grey35") +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey74") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey74") +
  geom_point(alpha = 0.42, size = 0.92, stroke = 0) +
  facet_wrap(~ engine, nrow = 2) +
  scale_color_manual(values = size_palette, drop = FALSE) +
  coord_equal(xlim = c(-axis_limit, axis_limit), ylim = c(-axis_limit, axis_limit)) +
  labs(
    title = "Environmental-effect recovery",
    subtitle = "Species x continuous-predictor effects; colours show benchmark dataset",
    x = "True effect (standardised)",
    y = "Estimated effect (standardised)",
    color = "Benchmark"
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 0.85, size = 2.1))) +
  theme_paper() +
  theme(
    axis.text.x = element_text(angle = 0),
    axis.title.y = element_text(margin = margin(r = 6)),
    legend.position = "right",
    legend.box.margin = margin(t = 0, b = 0),
    legend.margin = margin(t = 0, b = 0)
  )

ideal_prediction <- tibble(
  metric_label = factor(c("AUC", "Brier", "Prob. RMSE"), levels = c("AUC", "Brier", "Prob. RMSE")),
  ideal = c(1, 0, 0)
)

pB <- ggplot(prediction_metrics_long,
             aes(benchmark_label, value, color = engine, group = engine)) +
  geom_hline(data = ideal_prediction, aes(yintercept = ideal), inherit.aes = FALSE,
             linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_line(linewidth = 0.70, alpha = 0.88, na.rm = TRUE) +
  geom_point(size = 2.15, na.rm = TRUE) +
  facet_wrap(~ metric_label, scales = "free_y", nrow = 1) +
  scale_color_manual(values = engine_palette, drop = FALSE) +
  labs(
    title = "Prediction recovery",
    subtitle = "Occurrence probability benchmark; AUC higher, Brier and probability RMSE lower",
    x = "Benchmark size",
    y = "Metric value",
    color = "Engine"
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_paper() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "bottom",
    legend.text = element_text(size = 7.8)
  )

figure5 <- pA / pB +
  plot_layout(heights = c(1.55, 1.05), guides = "keep") +
  plot_annotation(
    title = "Universal Benchmark Cross-model Result Comparison",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 17, family = "Arial", hjust = 0),
      plot.tag = element_text(face = "bold", size = 15, family = "Arial"),
      plot.margin = margin(7, 7, 7, 7)
    )
  )

width_in <- 9.2
height_in <- 7.35
ggsave(file.path(output_root, "Figure5_Universal_Benchmark_Cross_Model_Comparison.pdf"),
       figure5, width = width_in, height = height_in, units = "in", device = cairo_pdf)
ggsave(file.path(output_root, "Figure5_Universal_Benchmark_Cross_Model_Comparison.png"),
       figure5, width = width_in, height = height_in, units = "in", dpi = 450, bg = "white")
ggsave(file.path(output_root, "Figure5_Universal_Benchmark_Cross_Model_Comparison.tiff"),
       figure5, width = width_in, height = height_in, units = "in", dpi = 600,
       compression = "lzw", bg = "white")
ggsave(file.path(output_root, "Figure5_Universal_Benchmark_Cross_Model_Comparison.svg"),
       figure5, width = width_in, height = height_in, units = "in", bg = "white")

caption <- c(
  "**Figure 5. Universal Benchmark cross-model result comparison.** JSDM Studio maps engine-specific outputs back to shared site, species and predictor identifiers across five latent-truth benchmark datasets: Built-in example, 120 sites x 12 species; Figure5 Small, 30 x 10; Figure5 Medium, 100 x 20; Figure5 Large, 200 x 50; and Figure5 Ultra-large, 300 x 100, each with three detection visits and eight prediction sites. **A**, true environmental effects plotted against estimated effects for continuous predictors after within-engine and within-size standardisation. Each point is one species-by-predictor effect; the dashed line indicates the one-to-one relationship expected under perfect recovery. **B**, occurrence prediction recovery across the same benchmark datasets, summarized by AUC, Brier score and probability RMSE against the simulated occurrence probability. AUC is higher-is-better, whereas Brier score and probability RMSE are lower-is-better. Coefficient standardisation is used because engines report effects on different link or observation scales. The previous output-completeness panel was removed; output completeness is reported in the reproducibility validation figure and audit tables rather than treated as a statistical recovery metric."
)
writeLines(caption, file.path(output_root, "Figure5_caption_draft.md"), useBytes = TRUE)

note <- c(
  "# Figure 5 source note",
  "",
  paste0("Project root: ", project_root),
  paste0("Figure 4 scaling source data: ", figure4_source_data),
  paste0("Scaling run log: ", run_log_path),
  "",
  "This Figure 5 revision uses one 120 x 12 Built-in benchmark run plus the four true model runs used for the revised Figure 4 scaling benchmark.",
  "Panel A reads each size-specific Universal Benchmark master directory, joins engine-standard effects to true species-by-predictor effects, and plots only continuous predictors: pH, moisture, canopy and elevation.",
  "Panel B recomputes occurrence-compatible AUC, Brier score and probability RMSE directly from each master directory's unified prediction output and simulated truth files.",
  "The previous output-completeness Panel C was removed. Output completeness belongs in the reproducibility validation figure, not in this cross-model statistical-result comparison.",
  "All source data for the plotted panels are exported in this directory."
)
writeLines(note, file.path(output_root, "Figure5_Universal_Benchmark_Cross_Model_Comparison_source_note.md"), useBytes = TRUE)

manifest <- tibble(
  file = c(
    "Figure5_Universal_Benchmark_Cross_Model_Comparison.pdf",
    "Figure5_Universal_Benchmark_Cross_Model_Comparison.png",
    "Figure5_Universal_Benchmark_Cross_Model_Comparison.tiff",
    "Figure5_Universal_Benchmark_Cross_Model_Comparison.svg",
    "Figure5_caption_draft.md",
    "Figure5_Universal_Benchmark_Cross_Model_Comparison_source_note.md",
    "figure5_benchmark_size_table.csv",
    "figure5_model_status_by_size.csv",
    "figure5_panelA_effect_recovery_points.csv",
    "figure5_panelA_effect_recovery_metrics.csv",
    "figure5_panelB_prediction_recovery_metrics.csv",
    "figure5_panelB_prediction_recovery_metrics_long.csv"
  ),
  role = c(
    "vector_figure",
    "raster_preview",
    "print_raster",
    "editable_vector",
    "caption",
    "source_note",
    "benchmark_size_table",
    "model_status_by_size",
    "panel_a_points",
    "panel_a_metrics",
    "panel_b_metrics",
    "panel_b_metrics_long"
  )
)
readr::write_csv(manifest, file.path(output_root, "Figure5_manifest.csv"))

message("Figure 5 written to: ", output_root)
message("Scaling source: ", figure4_source_data)
message("Benchmark sizes: ", paste(sprintf("%s=%sx%s", size_table$size_label, size_table$n_sites, size_table$n_species), collapse = "; "))
