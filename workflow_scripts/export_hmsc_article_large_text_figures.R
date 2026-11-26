#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else file.path("output", "BookCh11_Hmsc_full_Finland_prediction_20260827_155248")
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
result_dir <- file.path(base_dir, "results")
pred_dir <- file.path(base_dir, "predictions")
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(result_dir, "article_figures_large_text")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
unlink(Sys.glob(file.path(out_dir, "*large_text*")), force = TRUE)
unlink(Sys.glob(file.path(out_dir, "Fig_*")), force = TRUE)

theme_set(
  theme_classic(base_size = 18, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.7, colour = "black"),
      axis.ticks = element_line(linewidth = 0.7, colour = "black"),
      axis.text = element_text(colour = "black", size = 17),
      axis.title = element_text(colour = "black", size = 20),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold"),
      plot.title = element_text(size = 24, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 17, colour = "#4d5b6a"),
      panel.grid = element_blank()
    )
)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

read_csv_safe <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

first_col_to_id <- function(d, id_name = "site_id") {
  if (ncol(d) > 0 && (names(d)[1] == "" || names(d)[1] == "X")) names(d)[1] <- id_name
  d
}

pretty_term <- function(x, max_chars = 28) {
  out <- as.character(x)
  out <- gsub("`", "", out, fixed = TRUE)
  out <- gsub("^\\(Intercept\\)$", "Intercept", out)
  out <- gsub("^poly\\(([^,]+),[^)]*\\)([0-9]+)$", "\\1 poly\\2", out)
  out <- gsub("^poly\\(([^,]+),.*\\)$", "\\1 poly", out)
  out <- gsub("^I\\((.*)\\)$", "\\1", out)
  out <- gsub("^Random: spatial[_ ]unit$", "Spatial RE", out)
  out <- gsub("^Random: spatial unit$", "Spatial RE", out)
  out <- gsub("^Random: sample$", "Sample RE", out)
  out <- gsub("^hab$", "Habitat", out)
  out <- gsub("^clim$", "Climate", out)
  out <- gsub("_", " ", out, fixed = TRUE)
  out <- gsub("\\s+", " ", out)
  too_long <- nchar(out) > max_chars
  if (any(too_long)) out[too_long] <- abbreviate(out[too_long], minlength = max_chars, strict = TRUE)
  out
}

long_to_matrix <- function(d, row_col, col_col, value_col) {
  rows <- unique(as.character(d[[row_col]]))
  cols <- unique(as.character(d[[col_col]]))
  m <- matrix(NA_real_, nrow = length(rows), ncol = length(cols), dimnames = list(rows, cols))
  idx <- cbind(match(as.character(d[[row_col]]), rows), match(as.character(d[[col_col]]), cols))
  m[idx] <- as.numeric(d[[value_col]])
  m
}

save_base_plot <- function(draw_fun, stem, width = 10, height = 7, pointsize = 20, png_scale = 450) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  tif_path <- file.path(out_dir, paste0(stem, ".tiff"))
  grDevices::pdf(pdf_path, width = width, height = height, pointsize = pointsize, family = "sans")
  draw_fun()
  grDevices::dev.off()
  grDevices::png(png_path, width = width, height = height, units = "in", res = png_scale,
                 pointsize = pointsize, bg = "white", type = if (capabilities("cairo")) "cairo" else "windows")
  draw_fun()
  grDevices::dev.off()
  grDevices::tiff(tif_path, width = width, height = height, units = "in", res = png_scale,
                  pointsize = pointsize, bg = "white", compression = "lzw",
                  type = if (capabilities("cairo")) "cairo" else "windows")
  draw_fun()
  grDevices::dev.off()
  invisible(c(pdf = pdf_path, png = png_path, tiff = tif_path))
}

save_gg <- function(p, stem, width = 8, height = 6, dpi = 450) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  tif_path <- file.path(out_dir, paste0(stem, ".tiff"))
  svg_path <- file.path(out_dir, paste0(stem, ".svg"))
  ggsave(pdf_path, p, width = width, height = height, units = "in", device = cairo_pdf, family = "sans")
  if (requireNamespace("svglite", quietly = TRUE)) ggsave(svg_path, p, width = width, height = height, units = "in", device = svglite::svglite)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(png_path, p, width = width, height = height, units = "in", dpi = dpi, device = ragg::agg_png, bg = "white")
    ggsave(tif_path, p, width = width, height = height, units = "in", dpi = dpi, device = ragg::agg_tiff, compression = "lzw", bg = "white")
  } else {
    ggsave(png_path, p, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  }
  invisible(c(pdf = pdf_path, png = png_path, tiff = tif_path, svg = svg_path))
}

clip_cap <- function(x) {
  vals <- x[is.finite(x)]
  cap <- suppressWarnings(as.numeric(quantile(abs(vals), 0.98, na.rm = TRUE)))
  if (!is.finite(cap) || cap <= 0) cap <- max(abs(vals), na.rm = TRUE)
  if (!is.finite(cap) || cap <= 0) cap <- 1
  cap
}

beta_long <- read_csv_safe(file.path(result_dir, "parameter_estimates_Beta_species_environment_long_main.csv"))
beta <- long_to_matrix(beta_long, "predictor", "species", "estimate")
beta_support <- if ("support" %in% names(beta_long)) long_to_matrix(beta_long, "predictor", "species", "support") else NULL
beta_keep <- !tolower(rownames(beta)) %in% c("(intercept)", "intercept")
beta <- beta[beta_keep, , drop = FALSE]
if (!is.null(beta_support)) beta_support <- beta_support[beta_keep, , drop = FALSE]
beta_cap <- clip_cap(beta)

tree_path <- file.path(data_dir, "phylogeny.nwk")
tree <- NULL
if (file.exists(tree_path) && requireNamespace("ape", quietly = TRUE)) {
  tree <- tryCatch(ape::read.tree(tree_path), error = function(e) NULL)
}

draw_beta <- function() {
  beta_plot <- matrix(pmax(-beta_cap, pmin(beta_cap, beta)),
                      nrow = nrow(beta), ncol = ncol(beta), dimnames = dimnames(beta))
  species_order <- colnames(beta_plot)
  tree_note <- "phylogeny order"
  tree_available <- !is.null(tree)
  if (isTRUE(tree_available)) {
    keep <- intersect(tree$tip.label, colnames(beta_plot))
    if (length(keep) >= 2) {
      t2 <- ape::keep.tip(tree, keep)
      species_order <- t2$tip.label
    } else {
      tree_available <- FALSE
      tree_note <- "input order"
    }
  }
  if (!isTRUE(tree_available) && is.matrix(beta_plot) && ncol(beta_plot) > 2) {
    hc <- tryCatch(hclust(dist(t(beta_plot))), error = function(e) NULL)
    if (!is.null(hc)) {
      species_order <- hc$labels[hc$order]
      tree_note <- "Beta order"
    }
  }
  beta_plot <- beta_plot[, species_order, drop = FALSE]
  sup_plot <- if (!is.null(beta_support)) beta_support[, species_order, drop = FALSE] else NULL
  pal <- grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(121)
  graphics::layout(matrix(c(1, 2, 3), nrow = 1), widths = c(2.1, 4.8, 0.55))
  op <- par(oma = c(0.3, 0.3, 2.4, 0.8), xpd = NA)
  on.exit({ par(op); layout(1) }, add = TRUE)
  if (isTRUE(tree_available)) {
    par(mar = c(4.8, 0.1, 1.0, 0.1))
    t2 <- ape::keep.tip(tree, intersect(tree$tip.label, species_order))
    ape::plot.phylo(t2, direction = "rightwards", show.tip.label = FALSE,
                    edge.color = "#7f7f7f", edge.width = 1.6, no.margin = FALSE)
    mtext("Phylogeny", side = 3, line = 0.1, cex = 1.05, font = 2)
  } else {
    par(mar = c(4.8, 0.1, 1.0, 0.1))
    plot.new()
    text(0.5, 0.5, tree_note, cex = 1.1, col = "#4d5b6a")
  }
  par(mar = c(4.8, 2.0, 1.0, 0.4))
  image(seq_len(nrow(beta_plot)), seq_len(ncol(beta_plot)), beta_plot,
        axes = FALSE, xlab = "", ylab = "", zlim = c(-beta_cap, beta_cap), col = pal, useRaster = TRUE)
  axis(1, at = seq_len(nrow(beta_plot)), labels = pretty_term(rownames(beta_plot), 20),
       las = 2, cex.axis = 1.05, lwd = 1.2)
  axis(2, at = c(1, ncol(beta_plot)), labels = c("Species", ""), tick = FALSE, cex.axis = 0.95, las = 2)
  box(lwd = 1.2, col = "#404040")
  abline(v = seq(0.5, nrow(beta_plot) + 0.5, 1), col = adjustcolor("white", 0.55), lwd = 0.8)
  if (!is.null(sup_plot)) {
    supported <- which(sup_plot >= 0.95, arr.ind = TRUE)
    if (nrow(supported) > 0) points(supported[, 1], supported[, 2], pch = 16, cex = 0.55, col = adjustcolor("black", 0.70))
  }
  mtext("Predictor", side = 1, line = 4.1, cex = 1.15)
  mtext("Beta effects", side = 3, line = 0.1, cex = 1.2, font = 2)
  par(mar = c(4.8, 0.1, 1.0, 2.8))
  legend_y <- seq(-beta_cap, beta_cap, length.out = length(pal))
  image(1, legend_y, matrix(legend_y, nrow = 1), axes = FALSE, xlab = "", ylab = "",
        zlim = c(-beta_cap, beta_cap), col = pal, useRaster = TRUE)
  axis(4, at = pretty(c(-beta_cap, beta_cap), n = 5), labels = round(pretty(c(-beta_cap, beta_cap), n = 5), 1),
       cex.axis = 0.95, lwd = 1.2)
  mtext("Beta", side = 3, line = 0.1, cex = 1.0, font = 2)
  mtext("Species-environment response", outer = TRUE, side = 3, line = 0.9, adj = 0, cex = 1.15, font = 2)
}

gamma_long <- read_csv_safe(file.path(result_dir, "parameter_estimates_trait_environment_gamma_long_main.csv"))
gamma_long <- gamma_long[!tolower(gamma_long$predictor) %in% c("(intercept)", "intercept") &
                           !tolower(gamma_long$trait) %in% c("(intercept)", "intercept"), , drop = FALSE]
gamma <- long_to_matrix(gamma_long, "predictor", "trait", "estimate")
gamma_support <- if ("support" %in% names(gamma_long)) long_to_matrix(gamma_long, "predictor", "trait", "support") else NULL
gamma_cap <- clip_cap(gamma)
gamma_plot_df <- gamma_long
gamma_plot_df$predictor_label <- pretty_term(gamma_plot_df$predictor, 20)
gamma_plot_df$trait_label <- pretty_term(gamma_plot_df$trait, 18)
gamma_plot_df$predictor_label <- factor(gamma_plot_df$predictor_label, levels = rev(unique(gamma_plot_df$predictor_label)))
gamma_plot_df$trait_label <- factor(gamma_plot_df$trait_label, levels = unique(gamma_plot_df$trait_label))
gamma_plot_df$supported <- if ("support" %in% names(gamma_plot_df)) gamma_plot_df$support >= 0.95 else FALSE
p_gamma <- ggplot(gamma_plot_df, aes(x = trait_label, y = predictor_label, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_point(data = gamma_plot_df[gamma_plot_df$supported %in% TRUE, , drop = FALSE],
             shape = 16, size = 4.0, colour = "black") +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
                       limits = c(-gamma_cap, gamma_cap), oob = squish, name = "Gamma") +
  labs(title = "Trait-environment response", x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 19, colour = "black"),
    axis.text.y = element_text(size = 19, colour = "black"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "#404040", fill = NA, linewidth = 0.8),
    plot.margin = margin(18, 24, 18, 24)
  )

draw_gamma <- function() {
  g <- matrix(pmax(-gamma_cap, pmin(gamma_cap, gamma)),
              nrow = nrow(gamma), ncol = ncol(gamma), dimnames = dimnames(gamma))
  pal <- grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(121)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(5.2, 0.55))
  op <- par(oma = c(0.3, 1.2, 2.4, 1.0), xpd = NA)
  on.exit({ par(op); layout(1) }, add = TRUE)
  par(mar = c(4.1, 6.3, 0.8, 0.3))
  g2 <- g[nrow(g):1, , drop = FALSE]
  image(seq_len(ncol(g2)), seq_len(nrow(g2)), t(g2), axes = FALSE, xlab = "", ylab = "",
        zlim = c(-gamma_cap, gamma_cap), col = pal, useRaster = TRUE)
  axis(1, at = seq_len(ncol(g2)), labels = pretty_term(colnames(g2), 18), las = 2, cex.axis = 0.92, lwd = 1.2)
  axis(2, at = seq_len(nrow(g2)), labels = pretty_term(rev(rownames(g)), 20), las = 2, cex.axis = 0.92, lwd = 1.2)
  box(lwd = 1.2, col = "#404040")
  abline(v = seq(0.5, ncol(g2) + 0.5, 1), col = adjustcolor("white", 0.55), lwd = 0.8)
  abline(h = seq(0.5, nrow(g2) + 0.5, 1), col = adjustcolor("white", 0.55), lwd = 0.8)
  if (!is.null(gamma_support)) {
    s <- gamma_support[nrow(gamma_support):1, , drop = FALSE]
    supported <- which(s >= 0.95, arr.ind = TRUE)
    if (nrow(supported) > 0) points(supported[, 2], supported[, 1], pch = 16, cex = 1.0, col = adjustcolor("black", 0.70))
  }
  mtext("Trait-environment response", outer = TRUE, side = 3, line = 0.9, adj = 0, cex = 1.15, font = 2)
  par(mar = c(4.1, 0.1, 0.8, 3.3))
  legend_y <- seq(-gamma_cap, gamma_cap, length.out = length(pal))
  image(1, legend_y, matrix(legend_y, nrow = 1), axes = FALSE, xlab = "", ylab = "",
        zlim = c(-gamma_cap, gamma_cap), col = pal, useRaster = TRUE)
  axis(4, at = pretty(c(-gamma_cap, gamma_cap), n = 5), labels = round(pretty(c(-gamma_cap, gamma_cap), n = 5), 1),
       cex.axis = 0.78, lwd = 1.2)
  mtext("Gamma", side = 3, line = 0.1, cex = 1.0, font = 2)
}

rho_path <- file.path(result_dir, "parameter_estimates_phylogenetic_signal_rho_draws_main.csv")
rho <- if (file.exists(rho_path)) read_csv_safe(rho_path) else data.frame(rho = numeric())
p_rho <- ggplot(rho, aes(x = rho)) +
  geom_histogram(aes(y = after_stat(density)), bins = 12, fill = "#8fbdd3", colour = "white", linewidth = 0.5) +
  geom_vline(xintercept = median(rho$rho, na.rm = TRUE), linewidth = 1.4, colour = "#b2182b") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Phylogenetic signal", x = "rho", y = "Density") +
  theme(plot.margin = margin(12, 18, 12, 12))

assoc_long <- read_csv_safe(file.path(result_dir, "parameter_estimates_species_associations_long_main_random_level_1.csv"))
assoc_species <- sort(unique(c(assoc_long$species_i, assoc_long$species_j)))
assoc_mat <- matrix(NA_real_, length(assoc_species), length(assoc_species), dimnames = list(assoc_species, assoc_species))
assoc_mat[cbind(match(assoc_long$species_i, assoc_species), match(assoc_long$species_j, assoc_species))] <- assoc_long$estimate
assoc_mat[cbind(match(assoc_long$species_j, assoc_species), match(assoc_long$species_i, assoc_species))] <- assoc_long$estimate
assoc_df <- expand.grid(species_i = assoc_species, species_j = assoc_species, stringsAsFactors = FALSE)
assoc_df$estimate <- as.numeric(assoc_mat[cbind(match(assoc_df$species_i, assoc_species), match(assoc_df$species_j, assoc_species))])
p_assoc <- ggplot(assoc_df, aes(species_i, species_j, fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
                       limits = c(-1, 1), oob = squish, name = "Omega") +
  coord_fixed(expand = FALSE) +
  labs(title = "Residual associations", x = "Species", y = "Species") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "#404040", fill = NA, linewidth = 0.8),
    plot.margin = margin(12, 18, 12, 12)
  )

vp_long <- read_csv_safe(file.path(result_dir, "parameter_estimates_variance_partitioning_long_main.csv"))
vp_long$predictor_label <- pretty_term(vp_long$predictor, 22)
vp_mean <- aggregate(variance_proportion ~ predictor_label, vp_long, mean)
vp_mean <- vp_mean[order(vp_mean$variance_proportion), , drop = FALSE]
p_vp_mean <- ggplot(vp_mean, aes(x = variance_proportion, y = reorder(predictor_label, variance_proportion), fill = variance_proportion)) +
  geom_col(width = 0.72) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
  scale_fill_gradient(low = "#8ccfbf", high = "#2b6f9e", guide = "none") +
  labs(title = "Mean variance partitioning", x = "Mean proportion", y = NULL)
top_species <- names(sort(tapply(vp_long$variance_proportion, vp_long$species, max, na.rm = TRUE), decreasing = TRUE))[seq_len(min(18, length(unique(vp_long$species))))]
vp_top <- vp_long[vp_long$species %in% top_species, , drop = FALSE]
vp_top$species_label <- factor(pretty_term(vp_top$species, 18), levels = pretty_term(top_species, 18))
p_vp_species <- ggplot(vp_top, aes(x = species_label, y = variance_proportion, fill = predictor_label)) +
  geom_col(width = 0.82) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_brewer(palette = "Set2", name = "Predictor") +
  labs(title = "Species-level partitioning", x = "Species", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 14),
        legend.position = "bottom")
p_vp <- p_vp_species / p_vp_mean + plot_layout(heights = c(1.1, 0.9))

newcoords <- first_col_to_id(read_csv_safe(file.path(data_dir, "newcoords.csv")))
rich <- read_csv_safe(file.path(pred_dir, "newdata_predicted_species_richness_main.csv"))
focal <- read_csv_safe(file.path(pred_dir, "newdata_predicted_focal_species_main.csv"))
map_df <- merge(newcoords, rich, by = "site_id", all.x = TRUE)
map_df <- merge(map_df, focal, by = "site_id", all.x = TRUE)
p_map_rich <- ggplot(map_df, aes(x = x, y = y, colour = predicted_species_richness)) +
  geom_point(size = 0.75, alpha = 0.95) +
  coord_equal() +
  scale_colour_viridis_c(option = "C", name = "Richness") +
  labs(title = "Predicted richness", x = "x", y = "y") +
  theme(axis.text = element_text(size = 14))
p_map_focal <- ggplot(map_df, aes(x = x, y = y, colour = predicted_probability)) +
  geom_point(size = 0.75, alpha = 0.95) +
  coord_equal() +
  scale_colour_viridis_c(option = "C", name = "Probability") +
  labs(title = gsub("_", " ", unique(na.omit(map_df$species))[1] %||% "Focal species"), x = "x", y = "y") +
  theme(axis.text = element_text(size = 14))
p_map <- p_map_rich + p_map_focal + plot_annotation(title = "Full-grid prediction")

grad <- read_csv_safe(file.path(pred_dir, "gradient_response_curves_main_clim.csv"))
grad$panel <- ifelse(grad$curve_type == "summed_response", "Richness",
                     ifelse(grad$curve_type == "focal_species", gsub("_", " ", grad$response), paste("Community", grad$trait)))
p_grad <- ggplot(grad, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#9e9ac8", alpha = 0.28, na.rm = TRUE) +
  geom_line(linewidth = 1.3, colour = "#3f3b78") +
  facet_wrap(~panel, scales = "free_y", ncol = 3) +
  labs(title = "Climate response", x = "clim", y = "Prediction") +
  theme(plot.margin = margin(12, 18, 12, 12))

save_base_plot(draw_beta, "F01_Beta_phylo_large_text", width = 13.0, height = 7.6, pointsize = 20)
save_gg(p_gamma, "F02_Gamma_trait_env_large_text", width = 8.0, height = 6.2)
save_gg(p_rho, "F03_rho_large_text", width = 7.2, height = 5.2)
save_gg(p_assoc, "F04_Omega_assoc_large_text", width = 7.2, height = 6.4)
save_gg(p_vp, "F05_variance_partitioning_large_text", width = 10.5, height = 8.6)
save_gg(p_map, "F06_Finland_prediction_large_text", width = 11.5, height = 6.4)
save_gg(p_grad, "F07_climate_response_large_text", width = 11.5, height = 5.4)

param_pdf <- file.path(out_dir, "parameter_estimates_article_large_text.pdf")
grDevices::pdf(param_pdf, width = 11.5, height = 7.0, pointsize = 22, family = "sans", onefile = TRUE)
draw_beta()
print(p_gamma)
print(p_rho)
print(p_assoc)
print(p_vp)
grDevices::dev.off()

pred_pdf <- file.path(out_dir, "predictions_article_large_text.pdf")
grDevices::pdf(pred_pdf, width = 11.5, height = 6.4, pointsize = 22, family = "sans", onefile = TRUE)
print(p_map)
print(p_grad)
grDevices::dev.off()

manifest <- data.frame(
  figure = c(
    "F01_Beta_phylo_large_text",
    "F02_Gamma_trait_env_large_text",
    "F03_rho_large_text",
    "F04_Omega_assoc_large_text",
    "F05_variance_partitioning_large_text",
    "F06_Finland_prediction_large_text",
    "F07_climate_response_large_text",
    "parameter_estimates_article_large_text",
    "predictions_article_large_text"
  ),
  scientific_role = c(
    "Species by environment Beta effects ordered by phylogeny",
    "Trait by environment Gamma effects",
    "Phylogenetic signal rho posterior summary",
    "Species residual association matrix",
    "Variance partitioning by predictor and species",
    "Full-grid Finland prediction maps",
    "Climate response curves for richness, focal species and trait-weighted prediction",
    "Multi-page parameter-estimate figure set",
    "Multi-page prediction figure set"
  ),
  text_policy = "Large text, short labels, article-shrink ready",
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(out_dir, "ARTICLE_FIGURE_MANIFEST.csv"), row.names = FALSE)
writeLines(c(
  "# HMSC article figures with large text",
  "",
  "These files were regenerated from the fitted HMSC Chapter 11 output, not enlarged from the original PDFs.",
  "Text and axis labels are intentionally large and concise so the figures remain readable after shrinking in a manuscript layout.",
  "Matrix-style species figures keep species names in the source CSV tables rather than printing all names on the panel.",
  "",
  "Recommended use:",
  "- Use PDF or SVG for editable vector text.",
  "- Use TIFF for journal upload when raster files are requested.",
  "- Use PNG for quick manuscript drafting or slides.",
  "",
  "Black dots in Beta/Gamma heatmaps indicate posterior support >= 0.95 when available."
), file.path(out_dir, "README_article_figures_large_text.md"))

cat("ARTICLE_FIGURES_LARGE_TEXT_OK\n")
cat(out_dir, "\n")
