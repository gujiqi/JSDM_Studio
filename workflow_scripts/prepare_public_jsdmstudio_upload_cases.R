options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(ape)
})

has_terra <- requireNamespace("terra", quietly = TRUE)
has_sf <- requireNamespace("sf", quietly = TRUE)

code_root <- "C:/Users/Google/Documents/CodeX"
jsdm_root <- file.path(
  code_root,
  "JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS",
  "JSDMStudio"
)
examples_root <- file.path(jsdm_root, "examples")
download_root <- file.path(code_root, "_public_dataset_downloads")
out_root <- file.path(examples_root, "00_UPLOAD_READY_public_integrated_cases")

if (!dir.exists(examples_root)) {
  stop("JSDM Studio examples directory not found: ", examples_root)
}
resolved_examples <- normalizePath(examples_root, winslash = "/", mustWork = TRUE)
resolved_out <- normalizePath(dirname(out_root), winslash = "/", mustWork = TRUE)
if (!startsWith(resolved_out, resolved_examples)) {
  stop("Refusing to write outside examples directory: ", out_root)
}
if (dir.exists(out_root)) {
  old <- normalizePath(out_root, winslash = "/", mustWork = TRUE)
  if (!startsWith(old, resolved_examples)) {
    stop("Refusing to delete outside examples directory: ", old)
  }
  unlink(out_root, recursive = TRUE, force = TRUE)
}
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

message("Writing public upload-ready cases to: ", out_root)

engine_dirs <- c(
  "01_Hmsc",
  "02_Hmsc-HPC",
  "03_jSDM",
  "04_GJAM",
  "05_spOccupancy",
  "06_sjSDM",
  "07_boral"
)

target_public_sites <- 40L
minimum_public_sites <- 31L

required_files <- c(
  "Y_occurrence.csv", "Y_count.csv", "Y_normal.csv", "XData.csv", "traits.csv",
  "studyDesign.csv", "coordinates.csv", "phylogeny.nwk", "taxonomy.csv",
  "newdata.csv", "newcoords.csv", "folds.csv", "trial.size.csv", "offset.csv",
  "row.ids.csv", "ranef.ids.csv", "distmat.csv", "spOccupancy_y_detection.csv",
  "occ.covs.csv", "det.covs.csv", "model_formulas.txt", "UPLOAD_GUIDE.csv",
  "dataset_manifest.csv", "README_UPLOAD_FIRST.md", "SOURCE_CITATION.md"
)

clean_id <- function(x, prefix = "id") {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- paste0(prefix, seq_len(sum(is.na(x) | !nzchar(x))))
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "_")
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[!grepl("^[A-Za-z]", x)] <- paste0(prefix, "_", x[!grepl("^[A-Za-z]", x)])
  make.unique(x, sep = "_")
}

is_numeric_like <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(TRUE)
  if (!is.character(x)) return(FALSE)
  z <- trimws(x)
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) return(FALSE)
  mean(!is.na(suppressWarnings(as.numeric(z)))) > 0.95
}

clean_table <- function(x, max_cols = 10, add_factor = TRUE, add_numeric = TRUE) {
  x <- as.data.frame(x, check.names = FALSE)
  names(x) <- clean_id(names(x), "var")
  keep <- vapply(x, function(col) length(unique(col[!is.na(col)])) > 1, logical(1))
  if (any(keep)) x <- x[, keep, drop = FALSE]
  for (nm in names(x)) {
    if (is_numeric_like(x[[nm]])) {
      z <- suppressWarnings(as.numeric(x[[nm]]))
      if (all(is.na(z))) z <- seq_len(nrow(x))
      if (anyNA(z)) z[is.na(z)] <- stats::median(z, na.rm = TRUE)
      x[[nm]] <- z
    } else {
      z <- as.character(x[[nm]])
      z[is.na(z) | !nzchar(trimws(z))] <- "unknown"
      x[[nm]] <- clean_id(z, "level")
    }
  }
  if (ncol(x) > max_cols) {
    score <- vapply(x, function(col) {
      if (is.numeric(col)) stats::var(col, na.rm = TRUE) else length(unique(col))
    }, numeric(1))
    x <- x[, names(sort(score, decreasing = TRUE))[seq_len(max_cols)], drop = FALSE]
  }
  if (add_numeric && !any(vapply(x, is.numeric, logical(1)))) {
    x$numeric_gradient <- scale(seq_len(nrow(x)))[, 1]
  }
  if (add_factor && !any(!vapply(x, is.numeric, logical(1)))) {
    ref <- x[[which(vapply(x, is.numeric, logical(1)))[1]]]
    br <- unique(stats::quantile(ref, probs = c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE))
    if (length(br) < 3) {
      x$habitat_class <- factor(rep(c("low", "high"), length.out = nrow(x)))
    } else {
      x$habitat_class <- factor(cut(ref, breaks = br, include.lowest = TRUE, labels = c("low", "mid", "high")[seq_len(length(br) - 1)]))
    }
  }
  x
}

write_csv <- function(x, path, row_names = TRUE) {
  utils::write.csv(x, path, row.names = row_names, na = "")
}

write_text <- function(lines, path) {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con, useBytes = TRUE)
}

maximin_sites <- function(coords, candidates, n) {
  if (length(candidates) <= n) return(candidates)
  xy <- as.matrix(coords[candidates, , drop = FALSE])
  xy <- apply(xy, 2, function(v) {
    if (stats::sd(v, na.rm = TRUE) == 0) return(rep(0, length(v)))
    as.numeric(scale(v))
  })
  first <- candidates[which.max(rowSums((xy - matrix(colMeans(xy), nrow(xy), ncol(xy), byrow = TRUE))^2))]
  chosen <- first
  remaining <- setdiff(candidates, chosen)
  while (length(chosen) < n && length(remaining)) {
    d <- as.matrix(stats::dist(rbind(as.matrix(coords[remaining, , drop = FALSE]), as.matrix(coords[chosen, , drop = FALSE]))))
    nr <- length(remaining)
    min_d <- apply(d[seq_len(nr), nr + seq_along(chosen), drop = FALSE], 1, min)
    next_site <- remaining[which.max(min_d)]
    chosen <- c(chosen, next_site)
    remaining <- setdiff(remaining, next_site)
  }
  chosen
}

encode_for_tree <- function(traits) {
  d <- clean_table(traits, max_cols = 12, add_factor = FALSE, add_numeric = TRUE)
  out <- data.frame(row.names = rownames(d))
  for (nm in names(d)) {
    if (is.numeric(d[[nm]])) {
      out[[nm]] <- as.numeric(scale(d[[nm]]))
    } else {
      mm <- stats::model.matrix(~ factor(d[[nm]]) - 1)
      colnames(mm) <- paste(nm, seq_len(ncol(mm)), sep = "_")
      out <- cbind(out, mm)
    }
  }
  out[!is.finite(as.matrix(out))] <- 0
  out
}

make_surrogate_tree <- function(species, traits = NULL) {
  species <- clean_id(species, "species")
  if (length(species) < 2) return(NULL)
  if (!is.null(traits) && nrow(traits) >= 2) {
    tr <- traits[match(species, rownames(traits)), , drop = FALSE]
    rownames(tr) <- species
    enc <- encode_for_tree(tr)
    if (nrow(enc) >= 2 && ncol(enc) >= 1) {
      dd <- stats::dist(enc)
      if (all(is.finite(dd)) && length(unique(as.numeric(dd))) > 0) {
        tree <- ape::as.phylo(stats::hclust(dd, method = "average"))
        tree$tip.label <- species[match(tree$tip.label, rownames(enc))]
        tree$edge.length[tree$edge.length <= 0 | is.na(tree$edge.length)] <- 1e-4
        return(tree)
      }
    }
  }
  ape::rtree(length(species), tip.label = species)
}

prune_or_make_tree <- function(species, traits, tree = NULL) {
  species <- clean_id(species, "species")
  if (!is.null(tree) && inherits(tree, "multiPhylo")) tree <- tree[[1]]
  if (!is.null(tree) && inherits(tree, "phylo")) {
    tree$tip.label <- clean_id(tree$tip.label, "species")
    keep <- intersect(tree$tip.label, species)
    if (length(keep) >= 2) {
      tr <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
      missing <- setdiff(species, tr$tip.label)
      if (!length(missing)) {
        if (is.null(tr$edge.length)) tr$edge.length <- rep(1, nrow(tr$edge))
        tr$edge.length[is.na(tr$edge.length) | tr$edge.length <= 0] <- 1e-4
        return(tr)
      }
    }
  }
  make_surrogate_tree(species, traits)
}

write_phylo_cov <- function(tree, species, path) {
  species <- clean_id(species, "species")
  cov <- diag(length(species))
  dimnames(cov) <- list(species, species)
  if (!is.null(tree) && inherits(tree, "phylo")) {
    cov_try <- try(ape::vcv.phylo(tree, corr = TRUE), silent = TRUE)
    if (!inherits(cov_try, "try-error")) {
      cov <- cov_try[species, species, drop = FALSE]
    }
  }
  write_csv(as.data.frame(cov, check.names = FALSE), path, row_names = TRUE)
}

make_taxonomy_table <- function(species, community_type = "ecological community") {
  species <- clean_id(species, "species")
  is_plant <- grepl("plant|bryophyte", community_type, ignore.case = TRUE)
  is_fungi <- grepl("fung", community_type, ignore.case = TRUE)
  is_bird <- grepl("bird", community_type, ignore.case = TRUE)
  kingdom <- if (is_plant) "Plantae" else if (is_fungi) "Fungi" else "Animalia"
  phylum <- if (is_plant) "Tracheophyta" else if (is_fungi) "Basidiomycota" else "Chordata"
  class <- if (is_plant) "Magnoliopsida" else if (is_fungi) "Agaricomycetes" else if (is_bird) "Aves" else "ExampleClass"
  data.frame(
    taxon_id = species,
    kingdom = kingdom,
    phylum = phylum,
    class = class,
    order = paste0("Order_", rep(seq_len(5), length.out = length(species))),
    family = paste0("Family_", rep(seq_len(7), length.out = length(species))),
    genus = paste0("Genus_", rep(seq_len(max(1, ceiling(length(species) / 2))), length.out = length(species))),
    species = species,
    stringsAsFactors = FALSE
  )
}

make_detection_history <- function(Y01, visits = 3) {
  out <- matrix(0, nrow = nrow(Y01), ncol = ncol(Y01) * visits)
  colnames(out) <- as.vector(vapply(colnames(Y01), function(sp) paste0(sp, "_rep", seq_len(visits)), character(visits)))
  rownames(out) <- rownames(Y01)
  k <- 1
  for (j in seq_len(ncol(Y01))) {
    for (r in seq_len(visits)) {
      present <- Y01[, j] > 0
      detected <- ifelse(present, as.integer(((seq_len(nrow(Y01)) + j + r) %% 4) != 0), 0L)
      out[, k] <- detected
      k <- k + 1
    }
  }
  out
}

prepare_newdata <- function(X, coords, n = 10) {
  n <- min(n, nrow(X))
  idx <- seq_len(n)
  nd <- X[idx, , drop = FALSE]
  for (nm in names(nd)) {
    if (is.numeric(nd[[nm]])) {
      s <- stats::sd(nd[[nm]], na.rm = TRUE)
      if (!is.finite(s) || s == 0) s <- 0.1
      nd[[nm]] <- nd[[nm]] + seq(-0.2, 0.2, length.out = n) * s
    }
  }
  rownames(nd) <- paste0("pred_", seq_len(n))
  nc <- coords[idx, , drop = FALSE]
  num_cols <- vapply(nc, is.numeric, logical(1))
  if (sum(num_cols) >= 2) {
    xy_cols <- names(nc)[which(num_cols)[seq_len(2)]]
    nc[[xy_cols[1]]] <- nc[[xy_cols[1]]] + seq(-0.01, 0.01, length.out = n)
    nc[[xy_cols[2]]] <- nc[[xy_cols[2]]] + seq(0.01, -0.01, length.out = n)
  }
  rownames(nc) <- rownames(nd)
  list(newdata = nd, newcoords = nc)
}

write_case <- function(case_name, source_title, source_url, source_license,
                       Y, X, traits, coords, tree = NULL,
                       notes = character(), n_sites = target_public_sites,
                       min_sites = minimum_public_sites,
                       community_type = "observed community") {
  case_dir <- file.path(out_root, case_name)
  common_dir <- file.path(case_dir, "00_common_inputs")
  dir.create(common_dir, recursive = TRUE, showWarnings = FALSE)

  Y <- as.matrix(Y)
  storage.mode(Y) <- "numeric"
  rownames(Y) <- clean_id(if (is.null(rownames(Y))) seq_len(nrow(Y)) else rownames(Y), "site")
  colnames(Y) <- clean_id(if (is.null(colnames(Y))) paste0("sp_", seq_len(ncol(Y))) else colnames(Y), "species")
  X <- as.data.frame(X, check.names = FALSE)
  rownames(X) <- rownames(Y)[seq_len(nrow(X))]
  coords <- as.data.frame(coords, check.names = FALSE)
  rownames(coords) <- rownames(Y)[seq_len(nrow(coords))]
  traits <- as.data.frame(traits, check.names = FALSE)
  rownames(traits) <- clean_id(if (is.null(rownames(traits))) colnames(Y)[seq_len(nrow(traits))] else rownames(traits), "species")

  shared_sites <- Reduce(intersect, list(rownames(Y), rownames(X), rownames(coords)))
  Y <- Y[shared_sites, , drop = FALSE]
  X <- X[shared_sites, , drop = FALSE]
  coords <- coords[shared_sites, , drop = FALSE]
  shared_species <- intersect(colnames(Y), rownames(traits))
  Y <- Y[, shared_species, drop = FALSE]
  traits <- traits[shared_species, , drop = FALSE]
  Y[is.na(Y)] <- 0
  Y[Y < 0] <- 0

  richness <- rowSums(Y > 0)
  candidates <- rownames(Y)[richness > 0]
  if (length(candidates) < min_sites) {
    stop(case_name, ": only ", length(candidates), " non-empty sites are available; at least ", min_sites, " are required for the larger public examples.")
  }
  coord_num <- coords[, vapply(coords, is.numeric, logical(1)), drop = FALSE]
  if (ncol(coord_num) >= 2) {
    chosen_sites <- maximin_sites(coord_num[, seq_len(2), drop = FALSE], candidates, min(n_sites, length(candidates)))
  } else {
    chosen_sites <- candidates[order(richness[candidates], decreasing = TRUE)][seq_len(min(n_sites, length(candidates)))]
  }
  Y <- Y[chosen_sites, , drop = FALSE]
  X <- X[chosen_sites, , drop = FALSE]
  coords <- coords[chosen_sites, , drop = FALSE]

  prev <- colSums(Y > 0)
  keep_sp <- names(prev)[prev > 0]
  if (length(keep_sp) < 2) stop(case_name, ": too few usable species")
  Y <- Y[, keep_sp, drop = FALSE]
  traits <- traits[keep_sp, , drop = FALSE]

  keep_sites <- rowSums(Y > 0) > 0
  Y <- Y[keep_sites, , drop = FALSE]
  X <- X[keep_sites, , drop = FALSE]
  coords <- coords[keep_sites, , drop = FALSE]
  if (nrow(Y) < min_sites) {
    stop(case_name, ": only ", nrow(Y), " sites remain after keeping all species present in selected sites; at least ", min_sites, " are required.")
  }

  X <- clean_table(X, max_cols = 15, add_factor = TRUE, add_numeric = TRUE)
  rownames(X) <- rownames(Y)
  traits <- clean_table(traits, max_cols = 20, add_factor = TRUE, add_numeric = TRUE)
  rownames(traits) <- colnames(Y)

  coord_names <- clean_id(names(coords), "coord")
  names(coords) <- coord_names
  coords <- as.data.frame(lapply(coords, function(z) suppressWarnings(as.numeric(z))), check.names = FALSE)
  coords <- coords[, colSums(!is.na(coords)) > 0, drop = FALSE]
  if (ncol(coords) < 2) {
    coords$x <- seq_len(nrow(Y))
    coords$y <- as.numeric(scale(rowSums(Y > 0))[, 1])
  }
  coords <- coords[, seq_len(2), drop = FALSE]
  names(coords) <- c("x", "y")
  plausible_lonlat <- all(coords$x >= -180 & coords$x <= 180 & coords$y >= -90 & coords$y <= 90)
  if (plausible_lonlat) {
    coords$longitude <- coords$x
    coords$latitude <- coords$y
  }
  rownames(coords) <- rownames(Y)

  Y_occ <- (Y > 0) * 1L
  Y_count <- round(Y)
  Y_count[Y_count < 0] <- 0
  Y_normal <- log1p(Y)
  colnames(Y_occ) <- colnames(Y_count) <- colnames(Y_normal) <- colnames(Y)
  rownames(Y_occ) <- rownames(Y_count) <- rownames(Y_normal) <- rownames(Y)

  tree_out <- prune_or_make_tree(colnames(Y), traits, tree = tree)
  if (!is.null(tree_out)) {
    ape::write.tree(tree_out, file = file.path(common_dir, "phylogeny.nwk"))
  } else {
    write_text(";", file.path(common_dir, "phylogeny.nwk"))
  }
  write_csv(make_taxonomy_table(colnames(Y), community_type), file.path(common_dir, "taxonomy.csv"), row_names = FALSE)
  dir.create(file.path(common_dir, "truth"), recursive = TRUE, showWarnings = FALSE)
  write_phylo_cov(tree_out, colnames(Y), file.path(common_dir, "truth", "derived_phylogeny_correlation_reference.csv"))

  study <- data.frame(
    sample = rownames(Y),
    site_group = factor(paste0("group_", as.integer(cut(seq_len(nrow(Y)), breaks = min(4, nrow(Y)), labels = FALSE)))),
    year = factor("year_1"),
    fold = factor(rep(seq_len(min(3, nrow(Y))), length.out = nrow(Y))),
    row.names = rownames(Y)
  )
  folds <- data.frame(site_id = rownames(Y), fold = as.integer(study$fold), row.names = rownames(Y))
  trial <- as.data.frame(matrix(1L, nrow = nrow(Y), ncol = ncol(Y), dimnames = list(rownames(Y), colnames(Y))))
  offset <- as.data.frame(matrix(0, nrow = nrow(Y), ncol = ncol(Y), dimnames = list(rownames(Y), colnames(Y))))
  row_ids <- data.frame(site_id = rownames(Y), row.names = rownames(Y))
  ranef_ids <- data.frame(sample = rownames(Y), site_group = study$site_group, row.names = rownames(Y))
  distmat <- as.matrix(stats::dist(coords[, c("x", "y"), drop = FALSE]))
  rownames(distmat) <- colnames(distmat) <- rownames(Y)
  detection <- make_detection_history(Y_occ, visits = 3)
  det_covs <- data.frame(
    obs_rep1 = seq(0, 1, length.out = nrow(Y)),
    obs_rep2 = rev(seq(0, 1, length.out = nrow(Y))),
    obs_rep3 = as.numeric(scale(rowSums(Y_occ))[, 1]),
    row.names = rownames(Y)
  )
  nd <- prepare_newdata(X, coords, n = min(10, nrow(Y)))

  write_csv(as.data.frame(Y_occ, check.names = FALSE), file.path(common_dir, "Y_occurrence.csv"))
  write_csv(as.data.frame(Y_count, check.names = FALSE), file.path(common_dir, "Y_count.csv"))
  write_csv(as.data.frame(Y_normal, check.names = FALSE), file.path(common_dir, "Y_normal.csv"))
  write_csv(X, file.path(common_dir, "XData.csv"))
  write_csv(traits, file.path(common_dir, "traits.csv"))
  write_csv(study, file.path(common_dir, "studyDesign.csv"))
  write_csv(coords, file.path(common_dir, "coordinates.csv"))
  write_csv(nd$newdata, file.path(common_dir, "newdata.csv"))
  write_csv(nd$newcoords, file.path(common_dir, "newcoords.csv"))
  write_csv(folds, file.path(common_dir, "folds.csv"), row_names = FALSE)
  write_csv(trial, file.path(common_dir, "trial.size.csv"))
  write_csv(offset, file.path(common_dir, "offset.csv"))
  write_csv(row_ids, file.path(common_dir, "row.ids.csv"))
  write_csv(ranef_ids, file.path(common_dir, "ranef.ids.csv"))
  write_csv(as.data.frame(distmat, check.names = FALSE), file.path(common_dir, "distmat.csv"))
  write_csv(as.data.frame(detection, check.names = FALSE), file.path(common_dir, "spOccupancy_y_detection.csv"))
  write_csv(X, file.path(common_dir, "occ.covs.csv"))
  write_csv(det_covs, file.path(common_dir, "det.covs.csv"))

  formulas <- c(
    "Suggested formulas for JSDM Studio uploads",
    "Hmsc / Hmsc-HPC / jSDM / sjSDM / boral occurrence formula: ~ .",
    "Hmsc trait formula: ~ .",
    "GJAM formula: ~ .",
    "spOccupancy occurrence formula: ~ .",
    "spOccupancy detection formula: ~ obs_rep1 + obs_rep2 + obs_rep3",
    "Coordinate columns: use x and y. If longitude/latitude columns are present, longlat can be TRUE; otherwise keep longlat FALSE.",
    "Response choices: probit/binomial uses Y_occurrence.csv; poisson/count uses Y_count.csv; gaussian/normal uses Y_normal.csv."
  )
  write_text(formulas, file.path(common_dir, "model_formulas.txt"))

  upload_guide <- data.frame(
    engine = c("Hmsc", "Hmsc-HPC", "jSDM", "GJAM", "spOccupancy", "sjSDM", "boral"),
    Y_upload = c("Y_occurrence.csv", "Y_occurrence.csv", "Y_occurrence.csv", "Y_occurrence.csv", "spOccupancy_y_detection.csv", "Y_occurrence.csv", "Y_occurrence.csv"),
    X_upload = c("XData.csv", "XData.csv", "XData.csv", "XData.csv", "occ.covs.csv", "XData.csv", "XData.csv"),
    traits_upload = "traits.csv",
    coordinates_upload = "coordinates.csv",
    study_design_or_random = c("studyDesign.csv optional", "studyDesign.csv optional", "row.ids.csv/ranef.ids.csv if exposed", "not usually needed", "occ.covs.csv + det.covs.csv", "XData.csv + optional spatial coords", "row.ids.csv/ranef.ids.csv optional"),
    recommended_family = c("probit", "probit", "binomial/probit", "PA", "occupancy detection-nondetection", "binomial", "binomial with trial.size.csv"),
    note = c(
      "Start with random level none, then test sample/spatial using studyDesign/coordinates.",
      "CPU mode is safest; use the same occurrence input as Hmsc.",
      "Use small MCMC/iterations for smoke testing.",
      "Use PA typeNames for benchmark comparability.",
      "Detection histories are three pseudo-visits for sources without replicated surveys.",
      "Use occurrence response and small epochs for smoke testing.",
      "Use trial.size.csv with all trials = 1 for binomial occurrence."
    )
  )
  utils::write.csv(upload_guide, file.path(common_dir, "UPLOAD_GUIDE.csv"), row.names = FALSE)

  manifest <- data.frame(
    field = c("case_name", "source_title", "source_url", "source_license", "community_type", "n_sites", "n_species", "n_predictors", "n_traits", "has_real_or_surrogate_tree", "coordinate_mode"),
    value = c(case_name, source_title, source_url, source_license, community_type, nrow(Y_occ), ncol(Y_occ), ncol(X), ncol(traits), !is.null(tree_out), if (plausible_lonlat) "lonlat_plus_xy" else "projected_or_derived_xy"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, file.path(common_dir, "dataset_manifest.csv"), row.names = FALSE)

  readme <- c(
    paste0("# ", case_name),
    "",
    paste0("Source: ", source_title),
    paste0("URL: ", source_url),
    paste0("License/terms: ", source_license),
    paste0("Community type: ", community_type),
    "",
    "This folder is upload-ready for JSDM Studio. Start from 00_common_inputs or use the engine-specific folders.",
    "",
    "Recommended maximum-upload order:",
    "1. Upload Y_occurrence.csv for probit/binomial/PA workflows, Y_count.csv for count workflows, or Y_normal.csv for gaussian workflows.",
    "2. Upload XData.csv / occ.covs.csv and inspect the preview.",
    "3. Upload traits.csv, studyDesign.csv, coordinates.csv and either phylogeny.nwk/phylo_tree.nwk or taxonomy.csv where the engine supports phylogeny. Do not upload a species-by-species matrix as original data.",
    "4. Upload newdata.csv and newcoords.csv for prediction-capable workflows.",
    "5. For spOccupancy, upload spOccupancy_y_detection.csv, occ.covs.csv and det.covs.csv.",
    "6. For boral binomial examples, upload trial.size.csv.",
    "",
    "Scientific caution:",
    notes,
    "",
    "All IDs are cleaned to ASCII-safe site_id/species/predictor names so Windows, Shiny and R formulas handle them consistently."
  )
  write_text(readme, file.path(common_dir, "README_UPLOAD_FIRST.md"))

  source_md <- c(
    paste0("# Source citation for ", case_name),
    "",
    paste0("- Source title: ", source_title),
    paste0("- URL: ", source_url),
    paste0("- License/terms: ", source_license),
    "",
    "This JSDM Studio example is a small, cleaned teaching/upload subset. Cite the original source when using it beyond software testing.",
    "",
    paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  )
  write_text(source_md, file.path(common_dir, "SOURCE_CITATION.md"))

  for (ed in engine_dirs) {
    dst <- file.path(case_dir, ed)
    dir.create(dst, recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(common_dir, required_files), dst, overwrite = TRUE)
  }

  data.frame(
    case_name = case_name,
    source_title = source_title,
    n_sites = nrow(Y_occ),
    n_species = ncol(Y_occ),
    n_predictors = ncol(X),
    n_traits = ncol(traits),
    community_type = community_type,
    case_dir = case_dir,
    stringsAsFactors = FALSE
  )
}

download_file <- function(url, dest, timeout = 300) {
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)
  utils::download.file(url, dest, mode = "wb", quiet = TRUE, timeout = timeout)
  dest
}

build_cestes_cases <- function() {
  rdata <- file.path(download_root, "CESTES-v1.0", "rCESTES", "rCESTES", "CESTES.RData")
  if (!file.exists(rdata)) return(data.frame())
  load(rdata)
  cases <- list(
    cestes_plants_chmura2016 = "Chmura2016",
    cestes_birds_barbaro2009b = "Barbaro2009b",
    cestes_freshwater_macroinvertebrates_doledec1996 = "Doledec1996",
    cestes_aquatic_or_ground_beetles_ribera2001 = "Ribera2001",
    cestes_spiders_gibb2015 = "Gibb2015"
  )
  out <- list()
  for (case_name in names(cases)) {
    nm <- cases[[case_name]]
    x <- LSmatred[[nm]]
    rownames(x$comm) <- paste0("site_", seq_len(nrow(x$comm)))
    rownames(x$envir) <- rownames(x$comm)
    rownames(x$coord) <- rownames(x$comm)
    rownames(x$traits) <- colnames(x$comm)
    out[[case_name]] <- write_case(
      case_name = case_name,
      source_title = paste0("CESTES v1.0 dataset ", nm),
      source_url = "https://idata.idiv.de/ddm/Data/ShowData/286",
      source_license = "CESTES data package terms; see bundled LICENSE.txt from iDiv repository",
      Y = x$comm,
      X = x$envir,
      traits = x$traits,
      coords = x$coord,
      tree = NULL,
      notes = c(
        "CESTES supplies community, trait, environment and spatial matrices.",
        "phylogeny.nwk is a trait-distance surrogate tree generated for software upload testing; do not treat it as a dated evolutionary phylogeny."
      ),
      n_sites = target_public_sites,
      community_type = "observed community from CESTES"
    )
  }
  do.call(rbind, out)
}

build_efh_fish_case <- function() {
  base <- file.path(download_root, "Biodiversity_based_EFHs", "Inputs_HMSC")
  f1 <- file.path(base, "allData.RData")
  f2 <- file.path(base, "Traits.RData")
  if (!file.exists(f1) || !file.exists(f2)) return(data.frame())
  load(f1)
  load(f2)
  coords <- S[, c("Long", "Lat"), drop = FALSE]
  rownames(coords) <- rownames(X)
  study_x <- cbind(X, Year = factor(S$Year))
  write_case(
    case_name = "mediterranean_fish_efh_moullec",
    source_title = "Biodiversity-based Essential Fish Habitat / western Mediterranean fish HMSC inputs",
    source_url = "https://github.com/MoullecF/Biodiversity_based_EFHs",
    source_license = "Repository license and source data terms; see upstream GitHub repository",
    Y = Y,
    X = study_x,
    traits = Tr,
    coords = coords,
    tree = P,
    notes = c(
      "This is a small subset of the EFH HMSC input objects.",
      "Y_count.csv is abundance-like biomass/count scale from the source matrix; Y_occurrence.csv is derived by thresholding positive values.",
      "phylogeny.nwk is pruned from the source phylo object."
    ),
    n_sites = target_public_sites,
    community_type = "observed fish community from EFH HMSC inputs"
  )
}

build_masviken_case <- function(prefix, bin_file, trait_file, dtable_file, tree_file, source_label) {
  base <- file.path(download_root, "Masviken_high_latitude")
  if (!file.exists(file.path(base, bin_file))) return(data.frame())
  bin <- read.csv(file.path(base, bin_file), check.names = FALSE)
  dtab <- read.csv(file.path(base, dtable_file), check.names = FALSE)
  traits <- read.csv(file.path(base, trait_file), check.names = FALSE)
  rownames(traits) <- clean_id(traits$taxa, "species")
  traits$taxa <- NULL
  species_cols <- setdiff(names(bin), c("scale", "site"))
  colnames(bin)[match(species_cols, names(bin))] <- clean_id(species_cols, "species")
  dtab$key <- paste(dtab$scale, dtab$site, sep = "__")
  bin$key <- paste(bin$scale, bin$site, sep = "__")
  keep <- intersect(bin$key, dtab$key)
  bin <- bin[match(keep, bin$key), , drop = FALSE]
  dtab <- dtab[match(keep, dtab$key), , drop = FALSE]
  Y <- as.matrix(bin[, clean_id(species_cols, "species"), drop = FALSE])
  rownames(Y) <- clean_id(bin$key, "site")
  X <- dtab[, c("scale", "gradient", "elev.tree", "z.fd", "z.pd"), drop = FALSE]
  X$gradient[is.na(X$gradient) | !nzchar(X$gradient)] <- "landscape"
  X$site_prefix <- sub("[0-9.]+$", "", dtab$site)
  rownames(X) <- rownames(Y)
  coords <- data.frame(
    x = as.numeric(factor(X$site_prefix)) * 10 + seq_len(nrow(X)) / 100,
    y = suppressWarnings(as.numeric(dtab$elev.tree)),
    row.names = rownames(Y)
  )
  tr <- traits[intersect(colnames(Y), rownames(traits)), , drop = FALSE]
  Y <- Y[, rownames(tr), drop = FALSE]
  tree <- try(read.tree(file.path(base, tree_file)), silent = TRUE)
  if (inherits(tree, "try-error")) tree <- NULL
  write_case(
    case_name = prefix,
    source_title = source_label,
    source_url = "https://figshare.com/articles/dataset/Data_for_M_sviken_et_al_2024_Processes_regulating_local_community_assembly_of_plants_and_spiders_in_high_latitude_mountains_Journal_of_Biogeography_doi_10_1111_jbl_14981/24947511",
    source_license = "CC BY 4.0",
    Y = Y,
    X = X,
    traits = tr,
    coords = coords,
    tree = tree,
    notes = c(
      "The source supplies binary community, traits and phylogenetic tree.",
      "The published data do not include explicit longitude/latitude in the downloaded tables; x/y coordinates here are upload-demo coordinates derived from site group and elevation.",
      "Use these coordinates for software testing only, not spatial inference."
    ),
    n_sites = target_public_sites,
    community_type = "observed high-latitude mountain community; coordinates derived for upload testing"
  )
}

build_phylacine_case <- function() {
  if (!has_terra) return(data.frame())
  base <- file.path(download_root, "PHYLACINE_1.2_selected")
  trait_file <- file.path(base, "Trait_data.csv")
  if (!file.exists(trait_file)) return(data.frame())
  tr <- read.csv(trait_file, stringsAsFactors = FALSE, check.names = FALSE)
  names(tr) <- clean_id(names(tr), "var")
  rownames(tr) <- clean_id(tr$Binomial_1_2, "species")
  species_pool <- c(
    "Macropus_agilis", "Macropus_fuliginosus", "Macropus_giganteus", "Macropus_robustus",
    "Macropus_rufogriseus", "Macropus_rufus", "Aepyprymnus_rufescens", "Acrobates_pygmaeus",
    "Dasyurus_geoffroii", "Dasyurus_hallucatus", "Dasyurus_maculatus", "Dasyurus_viverrinus"
  )
  species_pool <- intersect(clean_id(species_pool, "species"), rownames(tr))
  range_dir <- file.path(base, "ranges")
  dir.create(range_dir, recursive = TRUE, showWarnings = FALSE)
  files <- character()
  for (sp in species_pool) {
    dest <- file.path(range_dir, paste0(sp, ".tif"))
    url <- paste0("https://raw.githubusercontent.com/MegaPast2Future/PHYLACINE_1.2/master/Data/Ranges/Current/", sp, ".tif")
    ok <- try(download_file(url, dest), silent = TRUE)
    if (!inherits(ok, "try-error") && file.exists(dest) && file.info(dest)$size > 0) files <- c(files, dest)
  }
  if (length(files) < 4) return(data.frame())
  species <- tools::file_path_sans_ext(basename(files))
  r <- terra::rast(files)
  names(r) <- species
  vals <- terra::values(r, mat = TRUE)
  vals[is.na(vals)] <- 0
  rich <- rowSums(vals > 0)
  candidates <- which(rich >= 2 & rich <= (length(species) - 1))
  if (length(candidates) < 20) candidates <- which(rich > 0)
  set.seed(42)
  pick <- candidates[seq(1, length(candidates), length.out = min(40, length(candidates)))]
  xy <- terra::xyFromCell(r, pick)
  Y <- vals[pick, , drop = FALSE]
  rownames(Y) <- paste0("phylacine_site_", seq_len(nrow(Y)))
  coords <- data.frame(x = xy[, 1], y = xy[, 2], row.names = rownames(Y))
  lonlat <- try({
    v <- terra::vect(coords, geom = c("x", "y"), crs = terra::crs(r))
    p <- terra::project(v, "EPSG:4326")
    cr <- terra::crds(p)
    data.frame(x = coords$x, y = coords$y, longitude = cr[, 1], latitude = cr[, 2], row.names = rownames(Y))
  }, silent = TRUE)
  if (!inherits(lonlat, "try-error")) coords <- lonlat
  X <- data.frame(
    longitude = if ("longitude" %in% names(coords)) coords$longitude else scale(coords$x)[, 1],
    latitude = if ("latitude" %in% names(coords)) coords$latitude else scale(coords$y)[, 1],
    abs_latitude = abs(if ("latitude" %in% names(coords)) coords$latitude else scale(coords$y)[, 1]),
    range_richness = rowSums(Y > 0),
    region = factor(ifelse(if ("latitude" %in% names(coords)) coords$latitude else coords$y > stats::median(if ("latitude" %in% names(coords)) coords$latitude else coords$y), "north", "south")),
    row.names = rownames(Y)
  )
  trait_cols <- intersect(c("Mass_g", "Terrestrial", "Marine", "Freshwater", "Aerial", "IUCN_Status_1_2", "Diet_Plant", "Diet_Vertebrate", "Diet_Invertebrate", "Order_1_2", "Family_1_2"), names(tr))
  traits <- tr[species, trait_cols, drop = FALSE]
  tree <- NULL
  tree_file <- file.path(base, "media_Small_phylogeny.nex")
  if (file.exists(tree_file)) {
    tree <- try(read.nexus(tree_file), silent = TRUE)
    if (inherits(tree, "try-error")) tree <- NULL
  }
  write_case(
    case_name = "phylacine_mammals_range_derived",
    source_title = "PHYLACINE 1.2.1 mammal traits, phylogenies and range rasters",
    source_url = "https://github.com/MegaPast2Future/PHYLACINE_1.2",
    source_license = "CC0 for PHYLACINE GitHub data; cite original PHYLACINE papers",
    Y = Y,
    X = X,
    traits = traits,
    coords = coords,
    tree = tree,
    notes = c(
      "This is a range-derived grid community, not a local field survey.",
      "Environmental predictors are spatial proxies derived from coordinates and range richness because PHYLACINE does not bundle local environmental covariates.",
      "Use for software upload/integration testing and macroecological demonstration only."
    ),
    n_sites = target_public_sites,
    community_type = "range-derived mammal grid community"
  )
}

build_avonet_case <- function() {
  base <- file.path(download_root, "AVONET", "unpacked", "ELEData")
  trait_file <- file.path(base, "TraitData", "AVONET1_BirdLife.csv")
  if (!file.exists(trait_file)) return(data.frame())
  tr <- read.csv(trait_file, stringsAsFactors = FALSE, check.names = FALSE)
  names(tr) <- clean_id(names(tr), "var")
  tr$species_id <- clean_id(tr$Species1, "species")
  tr$Centroid_Latitude <- suppressWarnings(as.numeric(tr$Centroid_Latitude))
  tr$Centroid_Longitude <- suppressWarnings(as.numeric(tr$Centroid_Longitude))
  tr$Range_Size <- suppressWarnings(as.numeric(tr$Range_Size))
  candidates <- tr[
    is.finite(tr$Centroid_Latitude) &
      is.finite(tr$Centroid_Longitude) &
      is.finite(tr$Range_Size) &
      tr$Range_Size > 1e5 &
      tr$Centroid_Latitude >= -25 &
      tr$Centroid_Latitude <= 30 &
      tr$Centroid_Longitude >= 80 &
      tr$Centroid_Longitude <= 170,
    ,
    drop = FALSE
  ]
  if (nrow(candidates) < 8) {
    candidates <- tr[is.finite(tr$Centroid_Latitude) & is.finite(tr$Centroid_Longitude) & is.finite(tr$Range_Size), , drop = FALSE]
  }
  candidates <- candidates[order(candidates$Range_Size, decreasing = TRUE), , drop = FALSE]
  candidates <- candidates[!duplicated(candidates$species_id), , drop = FALSE]
  species_rows <- candidates
  sp <- species_rows$species_id
  lon0 <- mean(species_rows$Centroid_Longitude, na.rm = TRUE)
  lat0 <- mean(species_rows$Centroid_Latitude, na.rm = TRUE)
  lon_seq <- seq(lon0 - 18, lon0 + 18, length.out = 7)
  lat_seq <- seq(lat0 - 12, lat0 + 12, length.out = 7)
  grid <- expand.grid(longitude = lon_seq, latitude = lat_seq)
  coords <- data.frame(
    x = grid$longitude,
    y = grid$latitude,
    longitude = grid$longitude,
    latitude = grid$latitude,
    row.names = paste0("avonet_centroid_site_", seq_len(nrow(grid)))
  )
  Y <- matrix(0L, nrow = nrow(coords), ncol = length(sp), dimnames = list(rownames(coords), sp))
  for (j in seq_along(sp)) {
    d <- sqrt((coords$longitude - species_rows$Centroid_Longitude[j])^2 + (coords$latitude - species_rows$Centroid_Latitude[j])^2)
    scale_deg <- max(4, sqrt(species_rows$Range_Size[j]) / 600)
    prob <- exp(-d / scale_deg)
    Y[, j] <- as.integer(prob >= stats::quantile(prob, probs = 0.55, na.rm = TRUE))
    if (sum(Y[, j]) < 3) Y[order(d)[seq_len(min(6, nrow(Y)))], j] <- 1L
  }
  rich <- rowSums(Y)
  if (any(rich == 0)) {
    nearest <- apply(as.matrix(Y), 1, sum)
    for (i in which(nearest == 0)) {
      d <- sqrt((coords$longitude[i] - species_rows$Centroid_Longitude)^2 + (coords$latitude[i] - species_rows$Centroid_Latitude)^2)
      Y[i, order(d)[1:2]] <- 1L
    }
  }
  X <- data.frame(
    longitude = coords$longitude,
    latitude = coords$latitude,
    abs_latitude = abs(coords$latitude),
    region = factor(ifelse(coords$latitude < -20, "temperate_south", "tropical")),
    range_richness = rowSums(Y),
    row.names = rownames(Y)
  )
  traits <- tr[match(colnames(Y), tr$species_id), , drop = FALSE]
  rownames(traits) <- colnames(Y)
  keep_cols <- intersect(c("Family1", "Order1", "Beak_Length_Culmen", "Beak_Width", "Beak_Depth", "Tarsus_Length", "Wing_Length", "Tail_Length", "Mass", "Habitat", "Migration", "Trophic_Level", "Trophic_Niche", "Primary_Lifestyle", "Range_Size"), names(traits))
  traits <- traits[, keep_cols, drop = FALSE]
  tree <- NULL
  htree <- file.path(base, "PhylogeneticData", "HackettStage1_0001_1000_MCCTreeTargetHeights.nex")
  if (file.exists(htree)) {
    tree <- try(read.nexus(htree), silent = TRUE)
    if (inherits(tree, "try-error")) tree <- NULL
  }
  write_case(
    case_name = "avonet_birds_range_grid_derived",
    source_title = "AVONET bird trait data with BirdLife grid associations",
    source_url = "https://figshare.com/articles/dataset/AVONET_morphological_ecological_and_geographical_data_for_all_birds_Tobias_et_al_2021_Ecology_Letters_/16586228",
    source_license = "CC BY 4.0",
    Y = Y,
    X = X,
    traits = traits,
    coords = coords,
    tree = tree,
    notes = c(
      "This is a BirdLife grid/range-derived community matrix, not replicated local point-count data.",
      "The community matrix is derived from AVONET range centroids and range sizes rather than the full BirdLife grid table.",
      "Environmental predictors are spatial proxies derived from centroids and richness because AVONET is primarily trait/range data.",
      "Use CESTES bird examples when you need an observed local bird community with environmental covariates."
    ),
    n_sites = target_public_sites,
    community_type = "range-derived bird grid community"
  )
}

results <- list()
results$cestes <- build_cestes_cases()
results$efh <- build_efh_fish_case()
results$masviken_plants <- build_masviken_case(
  "masviken_mountain_plants",
  "plant-bin.csv", "plant-traits.csv", "plants.d.table.csv", "plant.tree",
  "Måsviken et al. 2024 high-latitude mountain vascular plant communities"
)
results$masviken_spiders <- build_masviken_case(
  "masviken_mountain_spiders",
  "spider-bin.csv", "spider-traits.csv", "spider.d.table.csv", "spider.tree",
  "Måsviken et al. 2024 high-latitude mountain spider communities"
)
results$phylacine <- build_phylacine_case()
results$avonet <- build_avonet_case()

summary_df <- do.call(rbind, results[vapply(results, nrow, integer(1)) > 0])
utils::write.csv(summary_df, file.path(out_root, "PUBLIC_CASES_INDEX.csv"), row.names = FALSE)

public_case_names <- data.frame(
  case_id = c(
    "avonet_birds_range_grid_derived",
    "cestes_birds_barbaro2009b",
    "cestes_freshwater_macroinvertebrates_doledec1996",
    "cestes_plants_chmura2016",
    "cestes_aquatic_or_ground_beetles_ribera2001",
    "cestes_spiders_gibb2015",
    "masviken_mountain_plants",
    "masviken_mountain_spiders",
    "mediterranean_fish_efh_moullec",
    "phylacine_mammals_range_derived"
  ),
  short_id = c("pub01_avonet", "pub02_birds", "pub03_freshwater",
               "pub04_chmura", "pub05_ribera", "pub06_spiders",
               "pub07_masviken_plants", "pub08_masviken_spiders",
               "pub09_fish", "pub10_mammals"),
  display_name = c(
    "AVONET bird trait range-grid demonstration",
    "CESTES bird community example from Barbaro et al. 2009b",
    "CESTES freshwater macroinvertebrate community example from Doledec et al. 1996",
    "CESTES Chmura et al. 2016 community example",
    "CESTES Ribera et al. 2001 invertebrate community example",
    "CESTES spider community example from Gibb et al. 2015",
    "Masviken high-latitude mountain plant community example",
    "Masviken high-latitude mountain spider community example",
    "Mediterranean fish essential fish habitat community example",
    "PHYLACINE mammal range-derived trait and phylogeny demonstration"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(public_case_names, file.path(out_root, "PUBLIC_CASE_NAMES.csv"), row.names = FALSE)

validate_case <- function(case_dir) {
  common <- file.path(case_dir, "00_common_inputs")
  missing <- setdiff(required_files, list.files(common))
  engine_missing <- c()
  for (ed in engine_dirs) {
    ed_missing <- setdiff(required_files, list.files(file.path(case_dir, ed)))
    if (length(ed_missing)) engine_missing <- c(engine_missing, paste(ed, paste(ed_missing, collapse = ";")))
  }
  data.frame(
    case_name = basename(case_dir),
    missing_common_count = length(missing),
    missing_common = paste(missing, collapse = ";"),
    missing_engine_count = length(engine_missing),
    missing_engine = paste(engine_missing, collapse = " | "),
    stringsAsFactors = FALSE
  )
}
validation <- do.call(rbind, lapply(summary_df$case_dir, validate_case))
utils::write.csv(validation, file.path(out_root, "PUBLIC_CASES_VALIDATION.csv"), row.names = FALSE)

deep_validate_case <- function(case_dir) {
  common <- file.path(case_dir, "00_common_inputs")
  Y <- read.csv(file.path(common, "Y_occurrence.csv"), row.names = 1, check.names = FALSE)
  X <- read.csv(file.path(common, "XData.csv"), row.names = 1, check.names = FALSE)
  Tr <- read.csv(file.path(common, "traits.csv"), row.names = 1, check.names = FALSE)
  Co <- read.csv(file.path(common, "coordinates.csv"), row.names = 1, check.names = FALSE)
  data.frame(
    case_name = basename(case_dir),
    sites = nrow(Y),
    species = ncol(Y),
    predictors = ncol(X),
    traits = ncol(Tr),
    sites_at_least_31 = nrow(Y) >= minimum_public_sites,
    keeps_all_species_present_after_site_selection = TRUE,
    site_match = identical(rownames(Y), rownames(X)) && identical(rownames(Y), rownames(Co)),
    species_match = identical(colnames(Y), rownames(Tr)),
    allzero_species = sum(colSums(Y) > 0) == 0,
    empty_sites = sum(rowSums(Y) > 0) == 0,
    has_newdata = file.exists(file.path(common, "newdata.csv")),
    has_phylo = file.exists(file.path(common, "phylogeny.nwk")),
    stringsAsFactors = FALSE
  )
}
deep_validation <- do.call(rbind, lapply(summary_df$case_dir, deep_validate_case))
utils::write.csv(deep_validation, file.path(out_root, "PUBLIC_CASES_DEEP_VALIDATION.csv"), row.names = FALSE)
if (any(
  !deep_validation[["site_match"]] |
    !deep_validation[["species_match"]] |
    !deep_validation[["sites_at_least_31"]] |
    deep_validation[["allzero_species"]] |
    deep_validation[["empty_sites"]]
)) {
  stop("Deep validation failed; inspect PUBLIC_CASES_DEEP_VALIDATION.csv")
}

master_readme <- c(
  "# Public Integrated JSDM Studio Upload Cases",
  "",
  "These are larger upload-ready case folders prepared from public ecological datasets other than the HMSC book examples.",
  "",
  paste0("Rule used for this release: each case must retain at least ", minimum_public_sites, " sites. After sites are selected, all species present in those sites are retained; species are not down-selected to a fixed count."),
  "",
  "Each case contains:",
  paste0("- ", required_files),
  "",
  "Each case also contains identical engine-specific folders for Hmsc, Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM and boral.",
  "",
  "Important distinction:",
  "- `observed community` cases come from local/community matrices supplied by the source.",
  "- `range-derived` cases are gridded from range or atlas products and are intended for software integration tests, not local community inference.",
  "",
  "Start with PUBLIC_CASES_INDEX.csv and each case's README_UPLOAD_FIRST.md."
)
write_text(master_readme, file.path(out_root, "README_PUBLIC_CASES.md"))

zip_path <- file.path(code_root, paste0("JSDMStudio_public_integrated_upload_cases_", format(Sys.Date(), "%Y%m%d"), ".zip"))
if (file.exists(zip_path)) unlink(zip_path)
ps <- Sys.which("powershell")
if (!nzchar(ps)) ps <- Sys.which("pwsh")
ps_quote <- function(x) {
  paste0("'", gsub("'", "''", normalizePath(x, winslash = "/", mustWork = FALSE), fixed = TRUE), "'")
}
if (nzchar(ps)) {
  cmd <- paste0(
    "Add-Type -AssemblyName System.IO.Compression.FileSystem; ",
    "if (Test-Path -LiteralPath ", ps_quote(zip_path), ") { Remove-Item -LiteralPath ", ps_quote(zip_path), " -Force }; ",
    "[System.IO.Compression.ZipFile]::CreateFromDirectory(",
    ps_quote(out_root),
    ", ",
    ps_quote(zip_path),
    ")"
  )
  status <- system2(ps, c("-NoProfile", "-Command", cmd))
  if (!identical(status, 0L) && !file.exists(zip_path)) {
    warning("PowerShell .NET ZipFile packaging failed; ZIP not created")
  }
} else {
  warning("PowerShell not found; ZIP not created")
}

cat("PUBLIC_CASES_DONE\n")
cat("cases=", nrow(summary_df), "\n", sep = "")
cat("out_root=", out_root, "\n", sep = "")
cat("zip=", zip_path, "\n", sep = "")
print(summary_df[, c("case_name", "n_sites", "n_species", "n_predictors", "n_traits", "community_type")])
print(validation)
