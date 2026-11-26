#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

message("Building Figure 6: worked-example ecological outputs")

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
  project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

output_root <- file.path(project_root, "output")
candidate_dirs <- list.dirs(output_root, recursive = FALSE, full.names = TRUE)
candidate_dirs <- candidate_dirs[grepl("BookCh11_Hmsc_full_Finland_prediction_", basename(candidate_dirs))]

has_required <- function(path) {
  all(file.exists(
    file.path(path, "standard", "effects_species_environment.csv"),
    file.path(path, "standard", "associations_species_species.csv"),
    file.path(path, "predictions", "newdata_predicted_mean_main.csv"),
    file.path(path, "data", "newcoords.csv")
  ))
}

candidate_dirs <- candidate_dirs[vapply(candidate_dirs, has_required, logical(1))]
if (length(candidate_dirs) == 0) {
  stop("No complete BookCh11 Hmsc output directory was found.")
}

info <- file.info(candidate_dirs)
run_dir <- candidate_dirs[order(info$mtime, decreasing = TRUE)][[1]]
message("Using source output: ", basename(run_dir))

figure_dir <- file.path(project_root, "manuscript_figures", "Figure6_worked_example")
source_dir <- file.path(figure_dir, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

abbr_species <- function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) >= 2) {
      paste0(substr(p[[1]], 1, 1), ". ", paste(p[-1], collapse = " "))
    } else {
      p[[1]]
    }
  }, character(1))
}

predictor_label <- function(x) {
  out <- x
  out[out == "habCo"] <- "Coniferous\nhabitat"
  out[out == "habOp"] <- "Open\nhabitat"
  out[out == "habUrb"] <- "Urban\nhabitat"
  out[out == "habWe"] <- "Wetland\nhabitat"
  out[out == "poly(clim, degree = 2, raw = TRUE)1"] <- "Climate\nlinear"
  out[out == "poly(clim, degree = 2, raw = TRUE)2"] <- "Climate\nquadratic"
  out
}

nice_engine_title <- "Hmsc book Chapter 11 bird example"

effects <- read_csv(file.path(run_dir, "standard", "effects_species_environment.csv"))
effects$estimate <- suppressWarnings(as.numeric(effects$estimate))
effects <- effects[
  effects$effect_type == "environment" &
    !effects$predictor %in% c("(Intercept)") &
    !grepl("^Gamma_", effects$predictor) &
    is.finite(effects$estimate),
]

predictor_order <- c(
  "habCo",
  "habOp",
  "habUrb",
  "habWe",
  "poly(clim, degree = 2, raw = TRUE)1",
  "poly(clim, degree = 2, raw = TRUE)2"
)
effects <- effects[effects$predictor %in% predictor_order, ]

species_score <- aggregate(abs(effects$estimate), list(species = effects$species), sum, na.rm = TRUE)
names(species_score)[2] <- "abs_beta_sum"
species_score <- species_score[order(species_score$abs_beta_sum, decreasing = TRUE), ]
selected_species <- head(species_score$species, 12)

effect_matrix <- reshape(
  effects[effects$species %in% selected_species, c("species", "predictor", "estimate")],
  idvar = "species",
  timevar = "predictor",
  direction = "wide"
)
row.names(effect_matrix) <- effect_matrix$species
effect_values <- as.matrix(effect_matrix[, setdiff(names(effect_matrix), "species"), drop = FALSE])
effect_values[is.na(effect_values)] <- 0
if (nrow(effect_values) > 1) {
  species_order <- row.names(effect_values)[hclust(dist(effect_values))$order]
} else {
  species_order <- selected_species
}

effects_sel <- effects[effects$species %in% selected_species, ]
effects_sel$species_label <- abbr_species(effects_sel$species)
species_label_map <- setNames(abbr_species(species_order), species_order)
effects_sel$species_label <- factor(
  effects_sel$species_label,
  levels = rev(unname(species_label_map))
)
effects_sel$predictor_label <- factor(
  predictor_label(effects_sel$predictor),
  levels = predictor_label(predictor_order)
)
effect_limit <- as.numeric(stats::quantile(abs(effects_sel$estimate), probs = 0.95, na.rm = TRUE))
if (!is.finite(effect_limit) || effect_limit <= 0) effect_limit <- max(abs(effects_sel$estimate), na.rm = TRUE)
effects_sel$plot_estimate <- pmax(pmin(effects_sel$estimate, effect_limit), -effect_limit)

write.csv(effects_sel, file.path(source_dir, "figure6_panel_a_environmental_effects_selected.csv"), row.names = FALSE)

assoc <- read_csv(file.path(run_dir, "standard", "associations_species_species.csv"))
assoc$estimate <- suppressWarnings(as.numeric(assoc$estimate))
assoc <- assoc[
  assoc$species_i %in% selected_species &
    assoc$species_j %in% selected_species &
    is.finite(assoc$estimate),
]
assoc_rev <- assoc
names(assoc_rev)[names(assoc_rev) == "species_i"] <- "species_j_tmp"
names(assoc_rev)[names(assoc_rev) == "species_j"] <- "species_i"
names(assoc_rev)[names(assoc_rev) == "species_j_tmp"] <- "species_j"
diag_df <- data.frame(
  engine = "Hmsc",
  species_i = selected_species,
  species_j = selected_species,
  estimate = NA_real_,
  lower = NA_real_,
  upper = NA_real_,
  association_type = "Omega",
  scale = "engine_specific",
  comparable = FALSE,
  note = "Diagonal omitted"
)
assoc_plot <- rbind(assoc, assoc_rev[, names(assoc)], diag_df[, names(assoc)])
assoc_plot$species_i_label <- factor(abbr_species(assoc_plot$species_i), levels = unname(species_label_map))
assoc_plot$species_j_label <- factor(abbr_species(assoc_plot$species_j), levels = rev(unname(species_label_map)))
write.csv(assoc_plot, file.path(source_dir, "figure6_panel_b_hmsc_omega_selected.csv"), row.names = FALSE)

coords <- read_csv(file.path(run_dir, "data", "newcoords.csv"))
names(coords)[1] <- "site_id"
coords$x <- suppressWarnings(as.numeric(coords$x))
coords$y <- suppressWarnings(as.numeric(coords$y))

pred <- read_csv(file.path(run_dir, "predictions", "newdata_predicted_mean_main.csv"))
names(pred)[1] <- "site_id"
species_cols <- setdiff(names(pred), "site_id")
pred[species_cols] <- lapply(pred[species_cols], function(v) suppressWarnings(as.numeric(v)))
pred_sds <- vapply(pred[species_cols], stats::sd, numeric(1), na.rm = TRUE)
pred_species <- names(sort(pred_sds, decreasing = TRUE))[seq_len(min(3, length(pred_sds)))]

pred_long <- do.call(rbind, lapply(pred_species, function(sp) {
  data.frame(
    site_id = pred$site_id,
    species = sp,
    probability = pred[[sp]],
    stringsAsFactors = FALSE
  )
}))
pred_long <- merge(pred_long, coords[, c("site_id", "x", "y")], by = "site_id", all.x = TRUE)
pred_long$species_label <- factor(abbr_species(pred_long$species), levels = abbr_species(pred_species))
write.csv(pred_long, file.path(source_dir, "figure6_panel_c_spatial_predictions_selected.csv"), row.names = FALSE)

writeLines(
  c(
    paste0("source_output=", normalizePath(run_dir, winslash = "/", mustWork = TRUE)),
    paste0("selected_effect_species=", paste(selected_species, collapse = ",")),
    paste0("selected_prediction_species=", paste(pred_species, collapse = ",")),
    paste0("effect_color_limit_95pct_abs=", signif(effect_limit, 6)),
    "association_panel=Hmsc Omega only; engine-specific, not cross-engine ranked",
    "prediction_panel=posterior mean occurrence probability at new spatial prediction sites"
  ),
  con = file.path(source_dir, "figure6_source_notes.txt")
)

base_theme <- theme_classic(base_size = 10.5, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text = element_text(colour = "black", size = 8.8),
    axis.title = element_text(colour = "black", size = 10.5),
    legend.title = element_text(size = 9.5, face = "bold"),
    legend.text = element_text(size = 8.8),
    strip.text = element_text(size = 10.5, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9.5, colour = "grey25", hjust = 0),
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0, 1),
    panel.grid = element_blank()
  )

p_a <- ggplot(effects_sel, aes(x = predictor_label, y = species_label, fill = plot_estimate)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2b6aa8",
    mid = "white",
    high = "#b83a36",
    midpoint = 0,
    limits = c(-effect_limit, effect_limit),
    oob = scales::squish,
    name = "Beta\nmean"
  ) +
  labs(
    title = "Environmental coefficients",
    subtitle = "Species-specific Hmsc Beta estimates",
    x = NULL,
    y = NULL
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "right"
  )

p_b <- ggplot(assoc_plot, aes(x = species_i_label, y = species_j_label, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#2b6aa8",
    mid = "white",
    high = "#b83a36",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey92",
    oob = scales::squish,
    name = "Omega"
  ) +
  coord_fixed() +
  labs(
    title = "Residual association",
    subtitle = "Hmsc Omega; engine-specific estimand",
    x = NULL,
    y = NULL
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.5),
    axis.text.y = element_text(size = 7.5),
    legend.position = "right"
  )

p_c <- ggplot(pred_long, aes(x = x, y = y, colour = probability)) +
  geom_point(size = 0.45, alpha = 0.95) +
  facet_wrap(~ species_label, nrow = 1) +
  coord_equal() +
  scale_colour_gradientn(
    colours = c("#edf6f9", "#a7d8c9", "#41a6a6", "#1d5fa7", "#17356f"),
    limits = c(0, 1),
    oob = scales::squish,
    name = "Predicted\nprobability"
  ) +
  labs(
    title = "Spatial prediction",
    subtitle = "Posterior mean occurrence probability at new Finland grid points",
    x = NULL,
    y = NULL
  ) +
  base_theme +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.35),
    strip.background = element_rect(fill = "grey95", colour = "grey70", linewidth = 0.3)
  )

fig <- (p_a | p_b) / p_c +
  plot_layout(heights = c(1.15, 1.0), widths = c(1.2, 1.0), guides = "collect") +
  plot_annotation(
    title = "Worked example ecological outputs",
    subtitle = "JSDM Studio generates interpretable model results, not only run logs",
    caption = paste0(
      nice_engine_title,
      ". Association panel shows Hmsc Omega only; it is not a cross-engine ranking metric."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 17, face = "bold", family = "Arial", hjust = 0),
      plot.subtitle = element_text(size = 11, colour = "grey25", family = "Arial", hjust = 0),
      plot.caption = element_text(size = 8.8, colour = "grey30", family = "Arial", hjust = 0),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

out_base <- file.path(figure_dir, "Figure6_worked_example_ecological_outputs")
width_in <- 7.4
height_in <- 6.65

ggsave(paste0(out_base, ".pdf"), fig, width = width_in, height = height_in, device = cairo_pdf)
svglite::svglite(paste0(out_base, ".svg"), width = width_in, height = height_in)
print(fig)
dev.off()
ragg::agg_png(paste0(out_base, ".png"), width = width_in, height = height_in, units = "in", res = 450)
print(fig)
dev.off()
ragg::agg_tiff(paste0(out_base, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw")
print(fig)
dev.off()

manifest <- data.frame(
  file = basename(paste0(out_base, c(".pdf", ".svg", ".png", ".tiff"))),
  path = normalizePath(paste0(out_base, c(".pdf", ".svg", ".png", ".tiff")), winslash = "/", mustWork = TRUE),
  bytes = file.info(paste0(out_base, c(".pdf", ".svg", ".png", ".tiff")))$size,
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(figure_dir, "Figure6_manifest.csv"), row.names = FALSE)

message("Figure 6 files written to: ", normalizePath(figure_dir, winslash = "/", mustWork = TRUE))
