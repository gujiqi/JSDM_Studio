options(stringsAsFactors = FALSE)

script_path <- {
  args0 <- commandArgs(trailingOnly = FALSE)
  hit <- args0[startsWith(args0, "--file=")]
  if (length(hit)) normalizePath(sub("^--file=", "", hit[[1]]), winslash = "/", mustWork = FALSE) else
    normalizePath("workflow_scripts/prepare_jsdm_book_chapter_cases.R", winslash = "/", mustWork = FALSE)
}
app_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) stop("Cannot locate JSDM Studio app.R.", call. = FALSE)

raw_root <- file.path(app_dir, "examples", "book_chapter_benchmarks", "raw_source")
out_root <- file.path(app_dir, "examples", "book_chapter_benchmarks")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

sanitize_id <- function(x, make_unique = FALSE) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[nchar(x) == 0] <- "id"
  if (make_unique) x <- make.unique(x, sep = "_")
  x
}

scale01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || all(!is.finite(x))) return(rep(0, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

std_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || all(!is.finite(x))) return(rep(0, length(x)))
  x[!is.finite(x)] <- mean(x[is.finite(x)], na.rm = TRUE)
  sx <- stats::sd(x)
  if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

std_predictor <- function(x, n = length(x)) {
  z <- std_num(x)
  s <- suppressWarnings(stats::sd(z, na.rm = TRUE))
  if (!length(z) || !is.finite(s) || s == 0) z <- as.numeric(scale(seq_len(n)))
  z
}

select_species <- function(Y_count, n_species = 6L) {
  occ <- Y_count > 0
  prev <- colMeans(occ, na.rm = TRUE)
  total <- colSums(Y_count, na.rm = TRUE)
  ok <- which(prev > 0.05 & prev < 0.95 & total > 0)
  if (length(ok) < n_species) ok <- which(total > 0)
  score <- abs(prev[ok] - 0.5)
  ord <- ok[order(score, -total[ok], names(total)[ok])]
  names(total)[head(ord, min(n_species, length(ord)))]
}

make_star_newick <- function(species) {
  paste0("(", paste0(species, ":1", collapse = ","), ");")
}

make_taxonomy_cov <- function(species, taxonomy = NULL) {
  n <- length(species)
  C <- matrix(0.05, n, n, dimnames = list(species, species))
  diag(C) <- 1
  if (!is.null(taxonomy) && nrow(taxonomy)) {
    taxonomy$species_clean <- sanitize_id(taxonomy$species, make_unique = FALSE)
    tx <- taxonomy[match(species, taxonomy$species_clean), , drop = FALSE]
    for (i in seq_len(n)) for (j in seq_len(n)) {
      if (i == j) next
      if (!is.na(tx$genus[i]) && identical(tx$genus[i], tx$genus[j])) C[i, j] <- 0.70
      else if (!is.na(tx$family[i]) && identical(tx$family[i], tx$family[j])) C[i, j] <- 0.35
    }
  }
  C
}

write_table <- function(x, path, row.names = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = row.names)
}

matrix_from_response <- function(x) {
  if (inherits(x, "table")) {
    return(matrix(as.numeric(x), nrow = dim(x)[1], ncol = dim(x)[2], dimnames = dimnames(x)))
  }
  as.matrix(x)
}

derive_truth <- function(Y_occ, XData) {
  species <- colnames(Y_occ)
  site_id <- rownames(Y_occ)
  truth_prob <- pmin(pmax(as.matrix(Y_occ), 0.05), 0.95)
  latent <- qlogis(truth_prob)
  pred_truth <- as.data.frame(as.table(truth_prob), stringsAsFactors = FALSE)
  names(pred_truth) <- c("site_id", "species", "truth_probability")
  pred_truth$site_id <- as.character(pred_truth$site_id)
  pred_truth$species <- as.character(pred_truth$species)

  effect_rows <- list()
  for (sp in species) {
    y <- as.numeric(Y_occ[, sp])
    for (pred in names(XData)) {
      x <- XData[[pred]]
      est <- NA_real_
      if (is.numeric(x)) {
        if (stats::sd(x, na.rm = TRUE) > 0 && stats::sd(y, na.rm = TRUE) > 0) {
          est <- suppressWarnings(stats::cor(y, x, use = "pairwise.complete.obs"))
        }
      } else {
        f <- factor(x)
        if (nlevels(f) > 1) {
          means <- tapply(y, f, mean, na.rm = TRUE)
          est <- max(means, na.rm = TRUE) - min(means, na.rm = TRUE)
        }
      }
      effect_rows[[length(effect_rows) + 1L]] <- data.frame(species = sp, predictor = pred,
                                                            true_effect = est, stringsAsFactors = FALSE)
    }
  }
  true_effects <- do.call(rbind, effect_rows)

  assoc_rows <- list()
  for (i in seq_along(species)) for (j in seq_along(species)) {
    if (i >= j) next
    yi <- as.numeric(Y_occ[, i])
    yj <- as.numeric(Y_occ[, j])
    r <- if (stats::sd(yi) > 0 && stats::sd(yj) > 0) suppressWarnings(stats::cor(yi, yj)) else NA_real_
    assoc_rows[[length(assoc_rows) + 1L]] <- data.frame(species_i = species[i], species_j = species[j],
                                                        true_association = r, stringsAsFactors = FALSE)
  }
  list(prob = truth_prob, latent = latent, predictions = pred_truth,
       effects = true_effects, associations = do.call(rbind, assoc_rows))
}

make_detection <- function(Y_occ, n_visits = 3L) {
  out <- data.frame(row.names = rownames(Y_occ))
  for (sp in colnames(Y_occ)) {
    y <- as.numeric(Y_occ[, sp] > 0)
    for (r in seq_len(n_visits)) {
      # Deterministic pseudo-replicates for software compatibility, not original survey replicates.
      out[[paste0(sp, "_rep", r)]] <- ifelse(y == 1, as.integer((seq_along(y) + r) %% 3 != 0), 0L)
    }
  }
  out
}

write_case <- function(case_id, title, source_note, Y_count, XData, traits, studyDesign,
                       coordinates, phylogeny_nwk = NULL, phylo_cov = NULL,
                       newdata = NULL, newcoords = NULL, raw_files = character()) {
  outdir <- file.path(out_root, case_id)
  outdir_norm <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  outroot_norm <- normalizePath(out_root, winslash = "/", mustWork = TRUE)
  allowed <- c("ch06_plant_traits_whittaker", "ch07_deadwood_fungi", "ch11_finnish_birds")
  if (!case_id %in% allowed || !startsWith(outdir_norm, paste0(outroot_norm, "/"))) {
    stop("Refusing to overwrite a path outside generated book chapter examples: ", outdir, call. = FALSE)
  }
  if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
  dir.create(file.path(outdir, "truth"), recursive = TRUE, showWarnings = FALSE)

  Y_count <- as.data.frame(matrix_from_response(Y_count), check.names = FALSE)
  Y_count[] <- lapply(Y_count, function(x) pmax(0, round(suppressWarnings(as.numeric(x)))))
  rownames(Y_count) <- sanitize_id(rownames(Y_count), make_unique = TRUE)
  colnames(Y_count) <- sanitize_id(colnames(Y_count), make_unique = TRUE)
  Y_occ <- 1 * (as.matrix(Y_count) > 0)
  rownames(Y_occ) <- rownames(Y_count)
  colnames(Y_occ) <- colnames(Y_count)
  Y_normal <- log1p(as.matrix(Y_count))
  rownames(Y_normal) <- rownames(Y_count)
  colnames(Y_normal) <- colnames(Y_count)

  XData <- as.data.frame(XData, check.names = FALSE)
  rownames(XData) <- rownames(Y_count)
  XData$pH <- std_predictor(XData$pH, nrow(XData))
  XData$moisture <- std_predictor(XData$moisture, nrow(XData))
  XData$canopy <- std_predictor(XData$canopy, nrow(XData))
  XData$elevation <- std_predictor(XData$elevation, nrow(XData))
  XData$substrate <- factor(as.character(XData$substrate))

  traits <- as.data.frame(traits, check.names = FALSE)
  rownames(traits) <- colnames(Y_count)
  traits$life_form <- factor(as.character(traits$life_form))
  traits$height_mm <- std_num(traits$height_mm)
  traits$dispersal <- scale01(traits$dispersal)

  studyDesign <- as.data.frame(studyDesign, check.names = FALSE)
  rownames(studyDesign) <- rownames(Y_count)
  studyDesign$sample <- factor(as.character(studyDesign$sample))
  studyDesign$plot <- factor(as.character(studyDesign$plot))

  coordinates <- as.data.frame(coordinates, check.names = FALSE)
  rownames(coordinates) <- rownames(Y_count)
  coordinates$x <- suppressWarnings(as.numeric(coordinates$x))
  coordinates$y <- suppressWarnings(as.numeric(coordinates$y))
  if (any(!is.finite(coordinates$x))) coordinates$x <- seq_len(nrow(coordinates))
  if (any(!is.finite(coordinates$y))) coordinates$y <- seq_len(nrow(coordinates))

  if (is.null(newdata)) newdata <- XData[seq_len(min(8L, nrow(XData))), , drop = FALSE]
  rownames(newdata) <- paste0("new_", seq_len(nrow(newdata)))
  if (is.null(newcoords)) newcoords <- coordinates[seq_len(nrow(newdata)), , drop = FALSE]
  rownames(newcoords) <- rownames(newdata)

  if (is.null(phylogeny_nwk)) phylogeny_nwk <- make_star_newick(colnames(Y_count))
  if (is.null(phylo_cov)) {
    phylo_cov <- diag(length(colnames(Y_count)))
    dimnames(phylo_cov) <- list(colnames(Y_count), colnames(Y_count))
  }

  folds <- data.frame(site_id = rownames(Y_count), fold = rep(seq_len(3), length.out = nrow(Y_count)))
  trial <- as.data.frame(matrix(1, nrow = nrow(Y_count), ncol = ncol(Y_count), dimnames = dimnames(Y_count)))
  offset <- as.data.frame(matrix(0, nrow = nrow(Y_count), ncol = ncol(Y_count), dimnames = dimnames(Y_count)))
  rowids <- data.frame(row_id = rownames(Y_count), block = studyDesign$plot, row.names = rownames(Y_count))
  ranefids <- data.frame(observer = rep(c("obs_a", "obs_b", "obs_c"), length.out = nrow(Y_count)), row.names = rownames(Y_count))
  distmat <- as.matrix(stats::dist(coordinates[, c("x", "y"), drop = FALSE]))
  rownames(distmat) <- colnames(distmat) <- rownames(Y_count)
  det <- make_detection(Y_occ)
  det_cov <- data.frame(obs_rep1 = scale01(seq_len(nrow(Y_count))),
                        obs_rep2 = scale01(as.numeric(studyDesign$plot)),
                        obs_rep3 = scale01(coordinates$x),
                        row.names = rownames(Y_count))

  truth <- derive_truth(Y_occ, XData)

  write_table(Y_occ, file.path(outdir, "Y_occurrence.csv"))
  write_table(Y_count, file.path(outdir, "Y_count.csv"))
  write_table(round(Y_normal, 5), file.path(outdir, "Y_normal.csv"))
  write_table(XData, file.path(outdir, "XData.csv"))
  write_table(traits, file.path(outdir, "traits.csv"))
  write_table(studyDesign, file.path(outdir, "studyDesign.csv"))
  write_table(coordinates, file.path(outdir, "coordinates.csv"))
  writeLines(phylogeny_nwk, file.path(outdir, "phylogeny.nwk"), useBytes = TRUE)
  write_table(phylo_cov, file.path(outdir, "phylo_cov.csv"))
  write_table(newdata, file.path(outdir, "newdata.csv"))
  write_table(newcoords, file.path(outdir, "newcoords.csv"))
  write.csv(folds, file.path(outdir, "folds.csv"), row.names = FALSE)
  write_table(trial, file.path(outdir, "trial.size.csv"))
  write_table(offset, file.path(outdir, "offset.csv"))
  write_table(rowids, file.path(outdir, "row.ids.csv"))
  write_table(ranefids, file.path(outdir, "ranef.ids.csv"))
  write_table(distmat, file.path(outdir, "distmat.csv"))
  write_table(det, file.path(outdir, "spOccupancy_y_detection.csv"))
  write_table(XData, file.path(outdir, "occ.covs.csv"))
  write_table(det_cov, file.path(outdir, "det.covs.csv"))

  write_table(round(truth$latent, 5), file.path(outdir, "truth", "latent_occurrence.csv"))
  write_table(round(truth$prob, 5), file.path(outdir, "truth", "occurrence_probability.csv"))
  write.csv(truth$effects, file.path(outdir, "truth", "true_environment_effects.csv"), row.names = FALSE)
  write.csv(truth$associations, file.path(outdir, "truth", "true_species_associations.csv"), row.names = FALSE)
  write.csv(truth$predictions, file.path(outdir, "truth", "true_predictions.csv"), row.names = FALSE)
  writeLines(c(
    paste0("case_id: ", case_id),
    paste0("title: ", title),
    "source_type: empirical_book_chapter_case",
    "truth_note: empirical observed occurrence is used as the reference probability for software benchmark comparison",
    paste0("n_sites: ", nrow(Y_count)),
    paste0("n_species: ", ncol(Y_count))
  ), file.path(outdir, "truth", "data_generation_config.yml"), useBytes = TRUE)
  writeLines(c(
    "# Reproducibility note",
    "# This directory was generated by workflow_scripts/prepare_jsdm_book_chapter_cases.R.",
    "# Raw source files are stored under examples/book_chapter_benchmarks/raw_source/.",
    paste0("# Source note: ", source_note)
  ), file.path(outdir, "truth", "data_generation_script.R"), useBytes = TRUE)

  manifest <- data.frame(
    file = list.files(outdir, recursive = TRUE),
    case_id = case_id,
    title = title,
    source_note = source_note,
    stringsAsFactors = FALSE
  )
  write.csv(manifest, file.path(outdir, "benchmark_input_manifest.csv"), row.names = FALSE)
  writeLines(c(
    paste0("# ", title),
    "",
    source_note,
    "",
    "This is an all-engine JSDM Studio benchmark subset derived from the book chapter data.",
    "It is designed for software execution across Hmsc, Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM and boral.",
    "The spOccupancy detection file contains deterministic pseudo-replicates derived from occurrence because the original data are not replicated occupancy surveys.",
    "The truth/ files are empirical reference summaries for software comparison, not hidden ecological truth from a simulation."
  ), file.path(outdir, "README.md"), useBytes = TRUE)
  invisible(outdir)
}

prepare_ch06 <- function() {
  d <- read.csv(file.path(raw_root, "ch06", "data", "whittaker revisit data.csv"), check.names = FALSE)
  tax <- read.csv(file.path(raw_root, "ch06", "data", "taxonomy.csv"), check.names = FALSE)
  d$site <- sanitize_id(d$site, make_unique = FALSE)
  d$species_clean <- sanitize_id(d$species, make_unique = FALSE)
  Yall <- xtabs(value ~ site + species_clean, data = d)
  sel_sp <- select_species(Yall, 6L)
  # Select sites evenly over the topographic moisture gradient.
  env_by_site <- aggregate(env ~ site, d, mean)
  env_by_site <- env_by_site[match(rownames(Yall), env_by_site$site), ]
  site_order <- order(env_by_site$env)
  sel_sites <- sort(site_order[unique(round(seq(1, length(site_order), length.out = min(36L, length(site_order)))))])
  Y <- Yall[sel_sites, sel_sp, drop = FALSE]
  env <- env_by_site$env[sel_sites]
  site_ids <- rownames(Y)
  species <- colnames(Y)
  trait_by_sp <- aggregate(trait ~ species_clean, d, mean)
  trait <- trait_by_sp$trait[match(species, trait_by_sp$species_clean)]
  X <- data.frame(pH = env,
                  moisture = abs(env - median(env, na.rm = TRUE)) + 0.01 * seq_along(env),
                  canopy = sin(seq_along(env) / 3) + 0.05 * env,
                  elevation = seq_along(env), substrate = ifelse(env > median(env), "dry_slope", "mesic_slope"),
                  row.names = site_ids)
  traits <- data.frame(life_form = ifelse(trait > median(trait, na.rm = TRUE), "high_CN", "low_CN"),
                       height_mm = trait, dispersal = scale01(trait), row.names = species)
  study <- data.frame(sample = site_ids, plot = paste0("plot_", rep(seq_len(12), length.out = length(site_ids))),
                      row.names = site_ids)
  coords <- data.frame(x = seq_along(site_ids), y = as.numeric(scale(env)), row.names = site_ids)
  C <- make_taxonomy_cov(species, tax)
  write_case("ch06_plant_traits_whittaker",
             "Chapter 6 plant traits and topographic moisture gradient",
             "Book Chapter 6 / Section 6.7: plant abundance, leaf C:N trait and Whittaker topographic moisture gradient.",
             Y, X, traits, study, coords, phylogeny_nwk = make_star_newick(species), phylo_cov = C)
}

prepare_ch07 <- function() {
  d <- read.csv(file.path(raw_root, "ch07", "fungal data", "data.csv"), check.names = FALSE)
  species_cols <- setdiff(names(d)[5:ncol(d)], "unk")
  Yall <- as.matrix(d[, species_cols, drop = FALSE])
  rownames(Yall) <- sanitize_id(d$LogID, make_unique = TRUE)
  colnames(Yall) <- sanitize_id(colnames(Yall), make_unique = TRUE)
  sel_sp <- select_species(Yall, 6L)
  sel_sites <- seq_len(min(36L, nrow(Yall)))
  Y <- Yall[sel_sites, sel_sp, drop = FALSE]
  site_ids <- rownames(Y)
  species <- colnames(Y)
  dc <- d$DC[sel_sites]
  readcount <- d$readcount[sel_sites]
  unk <- d$unk[sel_sites]
  X <- data.frame(pH = as.numeric(dc) + 0.05 * scale01(readcount),
                  moisture = log1p(readcount),
                  canopy = ifelse(readcount > 0, unk / readcount, 0),
                  elevation = seq_along(site_ids),
                  substrate = paste0("decay_", dc),
                  row.names = site_ids)
  prev <- colMeans(Y > 0)
  mean_abund <- colMeans(Y)
  traits <- data.frame(life_form = ifelse(prev > median(prev), "frequent_fungus", "sparse_fungus"),
                       height_mm = mean_abund, dispersal = prev, row.names = species)
  study <- data.frame(sample = site_ids, plot = paste0("log_", site_ids), row.names = site_ids)
  coords <- data.frame(x = seq_along(site_ids), y = as.numeric(dc) + scale01(readcount), row.names = site_ids)
  write_case("ch07_deadwood_fungi",
             "Chapter 7 dead wood-inhabiting fungi sequencing data",
             "Book Chapter 7 / Section 7.9: fungal sequence counts from dead wood logs with decay class and read count.",
             Y, X, traits, study, coords, phylogeny_nwk = make_star_newick(species))
}

prepare_ch11 <- function() {
  d <- read.csv(file.path(raw_root, "ch11", "data", "data.csv"), check.names = FALSE)
  d <- d[d$Year == 2014, , drop = FALSE]
  rownames(d) <- sanitize_id(d$Route, make_unique = TRUE)
  trait_raw <- read.csv(file.path(raw_root, "ch11", "data", "traits.csv"), check.names = FALSE)
  species_cols <- names(d)[10:ncol(d)]
  Yall <- as.matrix(d[, species_cols, drop = FALSE])
  rownames(Yall) <- rownames(d)
  colnames(Yall) <- sanitize_id(colnames(Yall), make_unique = TRUE)
  sel_sp <- select_species(Yall, 6L)
  # Keep a habitat-balanced small subset where possible.
  ord <- order(d$Habitat, d$Route)
  sel_sites <- sort(ord[unique(round(seq(1, length(ord), length.out = min(36L, length(ord)))))])
  Y <- Yall[sel_sites, sel_sp, drop = FALSE]
  site_ids <- rownames(Y)
  species <- colnames(Y)
  dd <- d[sel_sites, , drop = FALSE]
  X <- data.frame(pH = dd$AprMay, moisture = dd$JunJul, canopy = dd$DJF,
                  elevation = dd$Effort, substrate = dd$Habitat,
                  row.names = site_ids)
  trait_raw$species_clean <- sanitize_id(trait_raw$Species, make_unique = FALSE)
  tt <- trait_raw[match(species, trait_raw$species_clean), , drop = FALSE]
  traits <- data.frame(life_form = tt$Migration, height_mm = tt$Mass,
                       dispersal = ifelse(is.finite(tt$Urb), tt$Urb, scale01(tt$Mass)),
                       row.names = species)
  study <- data.frame(sample = site_ids, plot = paste0("route_", site_ids), row.names = site_ids)
  coords <- data.frame(x = dd$x, y = dd$y, row.names = site_ids)
  grid <- read.csv(file.path(raw_root, "ch11", "data", "grid_1000.csv"), check.names = FALSE)
  newdata <- data.frame(pH = grid$AprMay[seq_len(8)], moisture = grid$JunJul[seq_len(8)],
                        canopy = grid$DJF[seq_len(8)], elevation = mean(dd$Effort, na.rm = TRUE),
                        substrate = grid$Habitat[seq_len(8)])
  newcoords <- data.frame(x = grid$x[seq_len(8)], y = grid$y[seq_len(8)])
  rownames(newdata) <- rownames(newcoords) <- paste0("new_grid_", seq_len(8))
  phylo_nwk <- make_star_newick(species)
  if (requireNamespace("ape", quietly = TRUE)) {
    tree <- tryCatch(ape::read.tree(file.path(raw_root, "ch11", "data", "CTree.tre")), error = function(e) NULL)
    if (!is.null(tree)) {
      tree$tip.label <- sanitize_id(tree$tip.label, make_unique = TRUE)
      keep <- intersect(species, tree$tip.label)
      if (length(keep) >= 2L) {
        tree2 <- ape::keep.tip(tree, keep)
        missing <- setdiff(species, tree2$tip.label)
        if (!length(missing)) phylo_nwk <- ape::write.tree(tree2)
      }
    }
  }
  write_case("ch11_finnish_birds",
             "Chapter 11 Finnish birds community data",
             "Book Chapter 11 / Section 11.1: Finnish bird presence-absence data with habitat, climate, coordinates, traits and phylogeny.",
             Y, X, traits, study, coords, phylogeny_nwk = phylo_nwk,
             newdata = newdata, newcoords = newcoords)
}

main <- function() {
  cases <- c(prepare_ch06(), prepare_ch07(), prepare_ch11())
  manifest <- data.frame(case_id = basename(cases),
                         path = normalizePath(cases, winslash = "/", mustWork = FALSE),
                         stringsAsFactors = FALSE)
  write.csv(manifest, file.path(out_root, "book_chapter_case_manifest.csv"), row.names = FALSE)
  cat("BOOK_CHAPTER_CASES_READY\n")
  print(manifest)
}

if (!interactive()) main()
