#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
cmd <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", cmd, value = TRUE))
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1], winslash = "/", mustWork = TRUE)) else getwd()
project_root <- if (length(args) >= 1) normalizePath(args[1], winslash = "/", mustWork = TRUE) else normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

output_root <- file.path(project_root, "docs", "paper_figures", "figure6_hmsc_worked_example_ecological_results")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

safe_read_csv <- function(path) readr::read_csv(path, show_col_types = FALSE, progress = FALSE)

latest_hmsc_ch11 <- function(project_root) {
  dirs <- list.dirs(file.path(project_root, "output"), recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^BookCh11_Hmsc_full_Finland_prediction_\\d", basename(dirs))]
  required <- file.path(dirs, "standard", "newdata_predictions_long.csv")
  dirs <- dirs[file.exists(required)]
  if (!length(dirs)) stop("No BookCh11_Hmsc_full_Finland_prediction_* output with newdata predictions was found.")
  dirs[which.max(file.info(dirs)$mtime)]
}

short_species <- function(x) {
  v <- gsub("_", " ", x)
  v <- sub("^([A-Z])[a-z]+ ", "\\1. ", v)
  v
}

read_site_table <- function(path) {
  x <- safe_read_csv(path)
  if (!"site_id" %in% names(x)) names(x)[1] <- "site_id"
  x %>% mutate(site_id = as.character(site_id))
}

source_dir <- if (length(args) >= 2) normalizePath(args[2], winslash = "/", mustWork = TRUE) else latest_hmsc_ch11(project_root)

paths <- list(
  beta = file.path(source_dir, "results", "parameter_estimates_Beta_species_environment_long_main.csv"),
  omega = file.path(source_dir, "results", "parameter_estimates_species_associations_long_main_random_level_1.csv"),
  predictions = file.path(source_dir, "standard", "newdata_predictions_long.csv"),
  newcoords = file.path(source_dir, "data", "newcoords.csv"),
  run_summary = file.path(source_dir, "standard", "run_summary.csv"),
  fit_metrics = file.path(source_dir, "standard", "fit_metrics.csv")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing)) stop("Missing required Hmsc worked-example source files: ", paste(missing, collapse = ", "))

beta_all <- safe_read_csv(paths$beta) %>%
  mutate(
    estimate = suppressWarnings(as.numeric(estimate)),
    support = suppressWarnings(as.numeric(support)),
    predictor_clean = recode(
      predictor,
      `poly(clim, degree = 2, raw = TRUE)1` = "clim",
      `poly(clim, degree = 2, raw = TRUE)2` = "clim^2",
      habCo = "habitat: coniferous",
      habOp = "habitat: open",
      habUrb = "habitat: urban",
      habWe = "habitat: wetland",
      .default = predictor
    ),
    species_label = short_species(species)
  )

species_rank <- beta_all %>%
  group_by(species, species_label) %>%
  summarise(max_abs_beta = max(abs(estimate), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_abs_beta)) %>%
  slice_head(n = 25)

selected_species <- species_rank$species
species_levels <- rev(species_rank$species_label)

beta_panel <- beta_all %>%
  filter(species %in% selected_species) %>%
  mutate(
    species_label = factor(species_label, levels = species_levels),
    predictor_clean = factor(predictor_clean, levels = c("habitat: coniferous", "habitat: open", "habitat: urban", "habitat: wetland", "clim", "clim^2"))
  )

beta_clip <- quantile(abs(beta_panel$estimate), probs = 0.95, na.rm = TRUE)
beta_clip <- max(beta_clip, 1)
beta_panel <- beta_panel %>%
  mutate(
    estimate_display = pmax(pmin(estimate, beta_clip), -beta_clip),
    support_class = case_when(
      is.na(support) ~ "not estimated",
      support >= 0.95 ~ ">=0.95",
      support >= 0.80 ~ "0.80-0.95",
      TRUE ~ "<0.80"
    )
  )

omega_all <- safe_read_csv(paths$omega) %>%
  mutate(
    estimate = suppressWarnings(as.numeric(estimate)),
    support = suppressWarnings(as.numeric(support)),
    species_i_label = short_species(species_i),
    species_j_label = short_species(species_j)
  )

omega_panel <- omega_all %>%
  filter(species_i %in% selected_species, species_j %in% selected_species) %>%
  transmute(
    species_i,
    species_j,
    species_i_label,
    species_j_label,
    estimate,
    support,
    association_type,
    scale,
    comparable,
    note
  )

diag_rows <- data.frame(
  species_i = selected_species,
  species_j = selected_species,
  species_i_label = short_species(selected_species),
  species_j_label = short_species(selected_species),
  estimate = 1,
  support = NA_real_,
  association_type = "diagonal",
  scale = "identity",
  comparable = "within_engine_only",
  note = "Diagonal added for heatmap display only.",
  stringsAsFactors = FALSE
)

omega_square <- bind_rows(
  omega_panel,
  omega_panel %>%
    transmute(
      species_i = species_j,
      species_j = species_i,
      species_i_label = species_j_label,
      species_j_label = species_i_label,
      estimate,
      support,
      association_type,
      scale,
      comparable,
      note
    ),
  diag_rows
) %>%
  mutate(
    species_i_label = factor(species_i_label, levels = rev(species_levels)),
    species_j_label = factor(species_j_label, levels = species_levels)
  )

pred_all <- safe_read_csv(paths$predictions) %>%
  mutate(
    site_id = as.character(site_id),
    species = as.character(response_id),
    predicted_mean = suppressWarnings(as.numeric(predicted_mean)),
    predicted_lower = suppressWarnings(as.numeric(predicted_lower)),
    predicted_upper = suppressWarnings(as.numeric(predicted_upper))
  )

pred_species <- pred_all %>%
  group_by(species) %>%
  summarise(
    n_points = n(),
    mean_prediction = mean(predicted_mean, na.rm = TRUE),
    sd_prediction = sd(predicted_mean, na.rm = TRUE),
    min_prediction = min(predicted_mean, na.rm = TRUE),
    max_prediction = max(predicted_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sd_prediction)) %>%
  slice_head(n = 3)

coords <- read_site_table(paths$newcoords) %>%
  mutate(
    x = suppressWarnings(as.numeric(x)),
    y = suppressWarnings(as.numeric(y))
  )

prediction_panel <- pred_all %>%
  filter(species %in% pred_species$species) %>%
  left_join(coords, by = "site_id") %>%
  filter(is.finite(x), is.finite(y), is.finite(predicted_mean)) %>%
  mutate(species_label = factor(short_species(species), levels = short_species(pred_species$species)))

run_summary <- safe_read_csv(paths$run_summary)
fit_metrics <- safe_read_csv(paths$fit_metrics)

readr::write_csv(beta_all, file.path(output_root, "figure6_panelA_beta_all_species_predictors.csv"))
readr::write_csv(beta_panel, file.path(output_root, "figure6_panelA_beta_display_subset.csv"))
readr::write_csv(omega_all, file.path(output_root, "figure6_panelB_omega_all_species_pairs.csv"))
readr::write_csv(omega_square, file.path(output_root, "figure6_panelB_omega_display_matrix.csv"))
readr::write_csv(pred_species, file.path(output_root, "figure6_panelC_selected_prediction_species.csv"))
readr::write_csv(prediction_panel, file.path(output_root, "figure6_panelC_spatial_prediction_points.csv"))
readr::write_csv(run_summary, file.path(output_root, "figure6_hmsc_run_summary.csv"))
readr::write_csv(fit_metrics, file.path(output_root, "figure6_hmsc_fit_metrics.csv"))

theme_paper <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(face = "bold", size = base_size + 1.5),
      plot.subtitle = element_text(size = base_size, color = "grey25"),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey94", color = "grey80", linewidth = 0.3),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 7, 4, 4)
    )
}

pA <- ggplot(beta_panel, aes(predictor_clean, species_label, fill = estimate_display)) +
  geom_tile(color = "white", linewidth = 0.18) +
  scale_fill_gradient2(
    low = "#2B6C9E",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-beta_clip, beta_clip),
    oob = squish,
    name = "Beta"
  ) +
  labs(
    title = "A. Environmental coefficients",
    subtitle = "Hmsc posterior mean Beta; strongest 25 species shown",
    x = NULL,
    y = NULL
  ) +
  theme_paper(base_size = 6.4) +
  theme(
    axis.text.y = element_text(size = 5.1),
    axis.text.x = element_text(size = 5.8),
    legend.position = "right"
  )

pB <- ggplot(omega_square, aes(species_j_label, species_i_label, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.12) +
  scale_fill_gradient2(
    low = "#2B6C9E",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    oob = squish,
    name = "Omega"
  ) +
  coord_equal() +
  labs(
    title = "B. Species association",
    subtitle = "Hmsc Omega; same 25 species as panel A",
    x = NULL,
    y = NULL
  ) +
  theme_paper(base_size = 6.1) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

pC <- ggplot(prediction_panel, aes(x, y, color = predicted_mean)) +
  geom_point(size = 0.22, alpha = 0.92) +
  facet_wrap(~species_label, nrow = 1) +
  coord_equal(expand = FALSE) +
  scale_color_gradientn(
    colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
    limits = c(0, 1),
    oob = squish,
    name = "Predicted\noccurrence"
  ) +
  labs(
    title = "C. Spatial prediction",
    subtitle = "Out-of-sample predictions over 9,454 Finland grid cells",
    x = "Grid x",
    y = "Grid y"
  ) +
  theme_paper(base_size = 6.5) +
  theme(
    axis.text = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 5.8),
    legend.position = "right"
  )

figure6 <- (pA | pB) / pC +
  plot_layout(heights = c(1.05, 0.95), widths = c(1.05, 1)) +
  plot_annotation(
    title = "Hmsc Worked Example: Ecological Model Outputs",
    theme = theme(plot.title = element_text(face = "bold", size = 12, family = "Arial"))
  )

width_in <- 7.4
height_in <- 6.2
ggsave(file.path(output_root, "Figure6_Hmsc_Worked_Example_Ecological_Results.pdf"), figure6, width = width_in, height = height_in, units = "in", device = cairo_pdf)
ggsave(file.path(output_root, "Figure6_Hmsc_Worked_Example_Ecological_Results.png"), figure6, width = width_in, height = height_in, units = "in", dpi = 450, bg = "white")
ggsave(file.path(output_root, "Figure6_Hmsc_Worked_Example_Ecological_Results.tiff"), figure6, width = width_in, height = height_in, units = "in", dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(output_root, "Figure6_Hmsc_Worked_Example_Ecological_Results.svg"), figure6, width = width_in, height = height_in, units = "in", bg = "white")

note <- c(
  "# Figure 6 source note",
  "",
  paste0("Project root: ", project_root),
  paste0("Hmsc worked-example source: ", source_dir),
  "",
  "Core conclusion: JSDM Studio produces ecological model outputs, not only software logs.",
  "Panel A uses Hmsc Beta posterior means from results/parameter_estimates_Beta_species_environment_long_main.csv.",
  "Panel B uses Hmsc Omega residual associations from results/parameter_estimates_species_associations_long_main_random_level_1.csv and should be interpreted within Hmsc only.",
  "Panel C uses standard/newdata_predictions_long.csv joined to data/newcoords.csv. The three mapped species were selected by the largest spatial standard deviation in predicted occurrence.",
  "The effect heatmap displays the 25 species with the largest absolute Beta values for readability; complete source data are exported alongside the figure.",
  "This figure is a worked-example ecological result figure, not a cross-engine ranking figure."
)
writeLines(note, file.path(output_root, "Figure6_Hmsc_Worked_Example_Ecological_Results_source_note.md"))

message("Figure 6 written to: ", output_root)
message("Hmsc worked-example source: ", source_dir)
