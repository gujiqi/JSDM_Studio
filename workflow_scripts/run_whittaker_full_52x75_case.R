# Rebuild and fit the complete Whittaker vascular-plant case used in Chapter 6.
#
# This script deliberately keeps the original 52 sites and 75 species. It creates
# auditable trait and phylogeny inputs, then runs occurrence and abundance Hmsc
# workflows through the same JSDM Studio adapter used by the GUI.

options(stringsAsFactors = FALSE)

if (.Platform$OS.type == "windows") {
  Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
}

args0 <- commandArgs(trailingOnly = FALSE)
script_arg <- args0[startsWith(args0, "--file=")]
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
} else {
  normalizePath("workflow_scripts/run_whittaker_full_52x75_case.R", winslash = "/", mustWork = FALSE)
}
app_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) stop("Cannot locate JSDM Studio app.R.", call. = FALSE)
setwd(app_dir)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.atomic(x) && is.na(x))) y else x
}

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(
    mode = "full",
    occurrence_samples = "1000",
    occurrence_transient = "2500",
    occurrence_thin = "5",
    occurrence_chains = "4",
    abundance_samples = "1000",
    abundance_transient = "2500",
    abundance_thin = "5",
    abundance_chains = "4",
    abundance_distribution = "lognormal poisson",
    skip_occurrence = "false",
    skip_abundance = "false",
    cv = "true",
    output_root = ""
  )
  for (a in args) {
    if (!startsWith(a, "--")) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    out[[kv[[1]]]] <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "true"
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x)) return(default)
  tolower(trimws(as.character(x)[1])) %in% c("true", "t", "1", "yes", "y")
}

as_int <- function(x, default, minimum = 0L) {
  z <- suppressWarnings(as.integer(x))
  if (!length(z) || !is.finite(z)) z <- as.integer(default)
  max(as.integer(minimum), z)
}

write_csv <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = row.names, na = "", fileEncoding = "UTF-8")
}

copy_tree <- function(src, dst) {
  files <- list.files(src, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  for (f in files) {
    rel <- substring(normalizePath(f, winslash = "/", mustWork = FALSE),
                     nchar(normalizePath(src, winslash = "/", mustWork = FALSE)) + 2L)
    target <- file.path(dst, rel)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(f, target, overwrite = TRUE)) stop("Could not copy source file: ", f, call. = FALSE)
  }
  invisible(files)
}

mode_value <- function(x) {
  x <- x[!is.na(x) & nzchar(trimws(as.character(x)))]
  if (!length(x)) return(NA_character_)
  tab <- sort(table(as.character(x)), decreasing = TRUE)
  sort(names(tab)[tab == max(tab)])[1]
}

canonical_binomial <- function(x) {
  vapply(strsplit(trimws(as.character(x)), "[[:space:]]+"), function(z) {
    if (length(z) < 2) stop("Species name does not contain a binomial: ", paste(z, collapse = " "), call. = FALSE)
    paste(z[1:2], collapse = " ")
  }, character(1))
}

species_id <- function(x) gsub(" ", "_", canonical_binomial(x), fixed = TRUE)

make_matrix <- function(dat, site_order, species_order, value_col) {
  out <- matrix(0, nrow = length(site_order), ncol = length(species_order),
                dimnames = list(site_order, species_order))
  ri <- match(as.character(dat[["site"]]), site_order)
  species_values <- as.character(dat[["species"]])
  ci <- match(species_values, species_order)
  if (anyNA(ci)) ci <- match(species_id(species_values), species_order)
  if (anyNA(ri) || anyNA(ci)) stop("Could not align raw site/species identifiers.", call. = FALSE)
  key <- paste(ri, ci, sep = ":")
  if (anyDuplicated(key)) stop("Raw Whittaker data contain duplicate site-species rows.", call. = FALSE)
  out[cbind(ri, ci)] <- as.numeric(dat[[value_col]])
  out
}

constant_by_group <- function(value, group, expected, label) {
  split_values <- split(value, group)
  bad <- names(split_values)[vapply(split_values, function(z) length(unique(z[is.finite(z)])) != 1L, logical(1))]
  if (length(bad)) stop(label, " is not constant within: ", paste(bad, collapse = ", "), call. = FALSE)
  out <- vapply(expected, function(g) unique(split_values[[g]][is.finite(split_values[[g]])])[1], numeric(1))
  if (any(!is.finite(out))) stop(label, " contains missing or non-finite values.", call. = FALSE)
  out
}

impute_categorical_trait <- function(values, taxonomy, tree, trait_name) {
  species <- taxonomy[["species_id"]]
  names(values) <- species
  original <- values
  method <- ifelse(is.na(values) | !nzchar(values), NA_character_, "global_database_species_match")
  donor <- ifelse(is.na(method), NA_character_, species)
  donor_distance <- ifelse(is.na(method), NA_real_, 0)

  for (rank in c("genus", "family")) {
    missing <- which(is.na(values) | !nzchar(values))
    if (!length(missing)) break
    for (i in missing) {
      peers <- which(taxonomy[[rank]] == taxonomy[[rank]][i] & !(is.na(values) | !nzchar(values)))
      if (!length(peers)) next
      values[i] <- mode_value(values[peers])
      chosen <- peers[values[peers] == values[i]][1]
      method[i] <- paste0(rank, "_mode_imputation")
      donor[i] <- species[chosen]
    }
  }

  missing <- which(is.na(values) | !nzchar(values))
  if (length(missing)) {
    distances <- ape::cophenetic.phylo(tree)
    observed <- which(!(is.na(values) | !nzchar(values)))
    if (!length(observed)) stop("No observed values are available for trait ", trait_name, ".", call. = FALSE)
    for (i in missing) {
      d <- distances[species[i], species[observed]]
      min_d <- min(d, na.rm = TRUE)
      nearest <- observed[which(d == min_d)]
      values[i] <- mode_value(values[nearest])
      chosen <- nearest[values[nearest] == values[i]][1]
      method[i] <- "phylogenetic_nearest_neighbour_imputation"
      donor[i] <- species[chosen]
      donor_distance[i] <- min_d
    }
  }

  list(
    values = values,
    audit = data.frame(
      species = species,
      trait = trait_name,
      original_value = original,
      final_value = values,
      method = method,
      donor_species = donor,
      patristic_distance = donor_distance,
      stringsAsFactors = FALSE
    )
  )
}

prepare_case <- function() {
  if (!requireNamespace("ape", quietly = TRUE)) stop("Package ape is required.", call. = FALSE)
  if (!requireNamespace("V.PhyloMaker2", quietly = TRUE)) stop("Package V.PhyloMaker2 is required.", call. = FALSE)

  raw_root <- file.path(app_dir, "examples", "book_chapter_benchmarks", "raw_source", "ch06")
  raw_data_file <- file.path(raw_root, "data", "whittaker revisit data.csv")
  raw_taxonomy_file <- file.path(raw_root, "data", "taxonomy.csv")
  if (!file.exists(raw_data_file) || !file.exists(raw_taxonomy_file)) {
    stop("Whittaker Chapter 6 source files were not found.", call. = FALSE)
  }

  case_dir <- file.path(app_dir, "examples", "00_ACTUAL_UPLOAD_DATA",
                        "01_vascular_plants__ch06_whittaker_full_52x75")
  original_dir <- file.path(case_dir, "00_original_source")
  common_dir <- file.path(case_dir, "01_common_inputs")
  occurrence_dir <- file.path(case_dir, "02_Hmsc_occurrence")
  abundance_dir <- file.path(case_dir, "03_Hmsc_abundance")
  log1p_dir <- file.path(case_dir, "04_Hmsc_log1p_normal_sensitivity")
  for (d in c(original_dir, common_dir, occurrence_dir, abundance_dir, log1p_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  copy_tree(raw_root, original_dir)

  raw <- read.csv(raw_data_file, check.names = FALSE, stringsAsFactors = FALSE)
  taxonomy_raw <- read.csv(raw_taxonomy_file, check.names = FALSE, stringsAsFactors = FALSE)
  required_raw <- c("site", "species", "trait", "env", "value")
  if (!all(required_raw %in% names(raw))) stop("Raw Whittaker table has an unexpected schema.", call. = FALSE)
  if (!all(c("species", "family", "genus") %in% names(taxonomy_raw))) stop("Taxonomy table has an unexpected schema.", call. = FALSE)

  site_raw <- unique(as.character(raw[["site"]]))
  site_num <- suppressWarnings(as.numeric(site_raw))
  site_raw <- if (all(is.finite(site_num))) site_raw[order(site_num)] else sort(site_raw)
  sites <- sprintf("site_%03d", seq_along(site_raw))
  site_map <- data.frame(original_site = site_raw, site_id = sites, stringsAsFactors = FALSE)
  raw[["site"]] <- site_map[["site_id"]][match(as.character(raw[["site"]]), site_map[["original_site"]])]

  taxonomy <- taxonomy_raw
  taxonomy[["original_species"]] <- taxonomy[["species"]]
  taxonomy[["canonical_binomial"]] <- canonical_binomial(taxonomy[["species"]])
  taxonomy[["species_id"]] <- gsub(" ", "_", taxonomy[["canonical_binomial"]], fixed = TRUE)
  if (anyDuplicated(taxonomy[["species_id"]])) stop("Canonical binomial names are not unique.", call. = FALSE)
  taxonomy <- taxonomy[, c("species_id", "original_species", "canonical_binomial", "family", "genus", "epithet")]
  species <- taxonomy[["species_id"]]

  raw_map <- setNames(taxonomy[["species_id"]], taxonomy[["original_species"]])
  if (length(setdiff(unique(raw[["species"]]), names(raw_map)))) stop("Raw response species are missing from taxonomy.", call. = FALSE)
  raw[["species"]] <- unname(raw_map[raw[["species"]]])
  if (nrow(raw) != length(sites) * length(species)) stop("The full site-by-species grid is incomplete.", call. = FALSE)

  Y_count <- make_matrix(raw, sites, species, "value")
  if (any(Y_count < 0) || any(abs(Y_count - round(Y_count)) > 1e-8)) {
    stop("Whittaker abundance values must be non-negative integers.", call. = FALSE)
  }
  Y_occurrence <- 1L * (Y_count > 0)
  Y_log1p <- log1p(Y_count)
  env <- constant_by_group(as.numeric(raw[["env"]]), raw[["site"]], sites, "TMG")
  leaf_cn <- constant_by_group(as.numeric(raw[["trait"]]), raw[["species"]], species, "leaf C:N")

  env_mean <- mean(env)
  env_sd <- stats::sd(env)
  tmg_z <- as.numeric((env - env_mean) / env_sd)
  XData <- data.frame(TMG = env, TMG_z = tmg_z, TMG_z2 = tmg_z^2,
                      row.names = sites, check.names = FALSE)
  study <- data.frame(sample = sites, row.names = sites, stringsAsFactors = FALSE)

  phylo_input <- data.frame(species = taxonomy[["species_id"]], genus = taxonomy[["genus"]],
                            family = taxonomy[["family"]], stringsAsFactors = FALSE)
  phylo_result <- V.PhyloMaker2::phylo.maker(
    sp.list = phylo_input,
    tree = V.PhyloMaker2::GBOTB.extended.TPL,
    nodes = V.PhyloMaker2::nodes.info.1.TPL,
    output.sp.list = TRUE,
    scenarios = "S3"
  )
  tree <- phylo_result[["scenario.3"]]
  phylo_status <- phylo_result[["species.list"]]
  if (length(tree[["tip.label"]]) != length(species) || !setequal(tree[["tip.label"]], species)) {
    stop("V.PhyloMaker2 output tips do not match the response species.", call. = FALSE)
  }
  tree <- ape::keep.tip(tree, species)
  phylo_cov <- ape::vcv.phylo(tree, corr = TRUE)
  phylo_cov <- phylo_cov[species, species, drop = FALSE]

  growth_url <- "https://raw.githubusercontent.com/traitecoevo/growthform/master/growth_form.csv"
  growth_snapshot <- file.path(original_dir, "external_traits", "growth_form_v0.2.3_snapshot.csv")
  dir.create(dirname(growth_snapshot), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(growth_snapshot)) {
    utils::download.file(growth_url, growth_snapshot, mode = "wb", quiet = FALSE)
  }
  growth <- read.csv(growth_snapshot, check.names = FALSE, stringsAsFactors = FALSE)
  growth_match <- growth[match(taxonomy[["canonical_binomial"]], growth[["sp"]]), , drop = FALSE]
  growth_match[["requested_species"]] <- taxonomy[["canonical_binomial"]]
  growth_match[["species_id"]] <- species
  growth_match <- growth_match[, c("species_id", "requested_species", "sp", "support", "woodiness", "nfixing")]

  support_imp <- impute_categorical_trait(as.character(growth_match[["support"]]), taxonomy, tree, "growth_support")
  wood_imp <- impute_categorical_trait(as.character(growth_match[["woodiness"]]), taxonomy, tree, "woodiness")
  support_labels <- c(F = "free_standing", C = "climber", E = "epiphyte", P = "parasite",
                      M = "mycoheterotroph", A = "aquatic", H = "hemiepiphyte")
  wood_labels <- c(W = "woody", H = "herbaceous")
  support_final <- unname(support_labels[support_imp[["values"]]])
  wood_final <- unname(wood_labels[wood_imp[["values"]]])
  if (anyNA(support_final) || anyNA(wood_final)) stop("Trait code translation failed.", call. = FALSE)

  TrData <- data.frame(
    leaf_CN = as.numeric(leaf_cn),
    woodiness = factor(wood_final, levels = c("herbaceous", "woody")),
    growth_support = factor(support_final,
                            levels = unique(c("free_standing", "climber", "epiphyte", "parasite",
                                              "mycoheterotroph", "aquatic", "hemiepiphyte")[
                                                c("free_standing", "climber", "epiphyte", "parasite",
                                                  "mycoheterotroph", "aquatic", "hemiepiphyte") %in% support_final
                                              ])),
    row.names = species,
    check.names = FALSE
  )

  leaf_audit <- data.frame(
    species = species,
    trait = "leaf_CN",
    original_value = as.character(leaf_cn),
    final_value = as.character(leaf_cn),
    method = "observed_in_Whittaker_source_data",
    donor_species = species,
    patristic_distance = 0,
    stringsAsFactors = FALSE
  )
  support_audit <- support_imp[["audit"]]
  support_audit[["original_value"]] <- unname(support_labels[support_audit[["original_value"]]])
  support_audit[["final_value"]] <- unname(support_labels[support_audit[["final_value"]]])
  wood_audit <- wood_imp[["audit"]]
  wood_audit[["original_value"]] <- unname(wood_labels[wood_audit[["original_value"]]])
  wood_audit[["final_value"]] <- unname(wood_labels[wood_audit[["final_value"]]])
  trait_audit <- rbind(leaf_audit, support_audit, wood_audit)

  pred_tmg <- seq(min(env), max(env), length.out = 101L)
  pred_z <- (pred_tmg - env_mean) / env_sd
  pred_ids <- sprintf("gradient_%03d", seq_along(pred_tmg))
  newdata <- data.frame(TMG = pred_tmg, TMG_z = pred_z, TMG_z2 = pred_z^2,
                        row.names = pred_ids, check.names = FALSE)
  newstudy <- data.frame(sample = pred_ids, row.names = pred_ids, stringsAsFactors = FALSE)

  name_map <- taxonomy[, c("species_id", "original_species", "canonical_binomial", "family", "genus")]
  phylo_audit <- merge(name_map, phylo_status, by.x = "species_id", by.y = "species", all.x = TRUE, sort = FALSE)
  phylo_audit <- phylo_audit[match(species, phylo_audit[["species_id"]]), , drop = FALSE]
  phylo_audit[["response_column_match"]] <- phylo_audit[["species_id"]] == species
  phylo_audit[["trait_row_match"]] <- phylo_audit[["species_id"]] == rownames(TrData)
  phylo_audit[["tree_tip_match"]] <- phylo_audit[["species_id"]] %in% tree[["tip.label"]]

  validation <- data.frame(
    check = c("sites", "species", "raw_long_rows", "Y_count_dimensions", "Y_occurrence_dimensions",
              "all_zero_species_count", "ubiquitous_species_count", "finite_environment", "three_traits",
              "no_missing_traits_after_imputation", "tree_tip_count", "tree_Y_alignment",
              "trait_Y_alignment", "prediction_rows"),
    observed = c(length(sites), length(species), nrow(raw), paste(dim(Y_count), collapse = "x"),
                 paste(dim(Y_occurrence), collapse = "x"), sum(colSums(Y_occurrence) == 0),
                 sum(colSums(Y_occurrence) == nrow(Y_occurrence)), all(is.finite(as.matrix(XData))),
                 ncol(TrData), sum(is.na(TrData)), length(tree[["tip.label"]]),
                 identical(colnames(Y_count), species) && setequal(tree[["tip.label"]], species),
                 identical(rownames(TrData), species), nrow(newdata)),
    expected = c(52, 75, 3900, "52x75", "52x75", 0, 1, TRUE, 3, 0, 75, TRUE, TRUE, 101),
    stringsAsFactors = FALSE
  )
  validation[["passed"]] <- as.character(validation[["observed"]]) == as.character(validation[["expected"]])
  if (!all(validation[["passed"]])) {
    stop("Prepared Whittaker data failed validation: ",
         paste(validation[["check"]][!validation[["passed"]]], collapse = ", "), call. = FALSE)
  }

  write_csv(raw, file.path(common_dir, "whittaker_long_full.csv"))
  write_csv(site_map, file.path(common_dir, "site_id_map.csv"))
  write_csv(name_map, file.path(common_dir, "species_name_map.csv"))
  write_csv(Y_count, file.path(common_dir, "Y_count.csv"), row.names = TRUE)
  write_csv(Y_occurrence, file.path(common_dir, "Y_occurrence.csv"), row.names = TRUE)
  write_csv(Y_log1p, file.path(common_dir, "Y_log1p_abundance.csv"), row.names = TRUE)
  write_csv(XData, file.path(common_dir, "XData.csv"), row.names = TRUE)
  write_csv(TrData, file.path(common_dir, "TrData.csv"), row.names = TRUE)
  write_csv(study, file.path(common_dir, "studyDesign.csv"), row.names = TRUE)
  write_csv(taxonomy, file.path(common_dir, "taxonomy.csv"))
  write_csv(growth_match, file.path(common_dir, "growthform_species_match_raw.csv"))
  write_csv(trait_audit, file.path(common_dir, "trait_provenance_and_imputation.csv"))
  write_csv(phylo_audit, file.path(common_dir, "phylogeny_match_and_alignment.csv"))
  write_csv(phylo_cov, file.path(common_dir, "phylogeny_correlation_derived.csv"), row.names = TRUE)
  write_csv(newdata, file.path(common_dir, "newdata_TMG_gradient.csv"), row.names = TRUE)
  write_csv(newstudy, file.path(common_dir, "prediction_studyDesign.csv"), row.names = TRUE)
  write_csv(validation, file.path(common_dir, "data_validation_report.csv"))
  ape::write.tree(tree, file = file.path(common_dir, "phylogeny_VPhyloMaker2_S3.nwk"))

  data_dictionary <- data.frame(
    file = c("Y.csv", "XData.csv", "TrData.csv", "studyDesign.csv", "phylogeny.nwk", "newdata.csv",
             "prediction_studyDesign.csv", "taxonomy.csv"),
    rows = c(52, 52, 75, 52, NA, 101, 101, 75),
    columns = c(75, 3, 3, 1, NA, 3, 1, 6),
    role = c("site-by-species response", "TMG predictor plus standardized linear and derived quadratic columns; the primary Chapter 6 model uses the linear TMG_z term",
             "leaf C:N plus woodiness and growth-support traits", "site random-effect identifiers",
             "dated vascular-plant phylogeny aligned to Y columns", "TMG prediction gradient",
             "new random-level identifiers for prediction", "original and canonical taxonomy"),
    stringsAsFactors = FALSE
  )
  write_csv(data_dictionary, file.path(common_dir, "DATA_DICTIONARY.csv"))

  writeLines(c(
    "# Complete Whittaker vascular-plant case (52 sites x 75 species)",
    "",
    "This folder preserves the complete Chapter 6 Whittaker dataset. No sites or species were selected out.",
    "",
    "## Responses",
    "",
    "- Y_occurrence.csv: 0/1 transformation of the original abundance values.",
    "- Y_count.csv: original non-negative integer plant abundance values.",
    "- Y_log1p_abundance.csv: deterministic log1p transformation used only for the all-species Gaussian sensitivity model.",
    "",
    "## Environment",
    "",
    "TMG is the observed topographic moisture gradient; larger values are topographically drier sites.",
    "TMG_z and TMG_z2 are deterministic standardized linear and derived quadratic columns. The primary model uses TMG_z only, matching the Chapter 6 linear structure; TMG_z2 is retained for transparent sensitivity analyses.",
    "",
    "## Traits",
    "",
    "leaf_CN is measured in the original Whittaker source data. woodiness and growth_support are matched",
    "to the versioned global growth-form database. Missing categorical values are filled first by genus mode,",
    "then family mode, then nearest patristic neighbour on the V.PhyloMaker2 tree. Every value and method is",
    "listed in trait_provenance_and_imputation.csv; imputed values must not be described as direct measurements.",
    "",
    "## Phylogeny",
    "",
    "phylogeny_VPhyloMaker2_S3.nwk was generated with V.PhyloMaker2 scenario S3 against GBOTB.extended.TPL.",
    "Infraspecific names were mapped to their canonical binomials for backbone matching; the complete map and",
    "the prune/bind status are in phylogeny_match_and_alignment.csv. All 75 tree tips match Y and TrData.",
    "",
    "## Important interpretation",
    "",
    "The original Whittaker data contain one measured trait (leaf C:N) and no geographic coordinates.",
    "The added categorical traits and generated phylogeny are transparent analytical enrichments, not original",
    "field measurements. The prediction grid is an environmental-gradient prediction, not a geographic map."
  ), file.path(case_dir, "README_FULL_WHITTAKER_CASE.md"), useBytes = TRUE)

  writeLines(c(
    "# Data sources and citations",
    "",
    "- Miller JED, Damschen EI, Ives AR (2018). Functional traits and community composition:",
    "  a comparison among community-weighted means, weighted correlations, and multilevel models.",
    "  Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210x.13119",
    "- Dryad source data: https://doi.org/10.5061/dryad.7gj0s3b",
    "- Taseski et al. (2019), global growth form database. https://doi.org/10.1002/ecy.2614",
    "- Zanne et al. (2014), woodiness source. https://doi.org/10.1038/nature12872",
    "- Growth-form snapshot URL: https://raw.githubusercontent.com/traitecoevo/growthform/master/growth_form.csv",
    "- Jin & Qian (2022), V.PhyloMaker2. https://doi.org/10.1016/j.pld.2022.05.005",
    "- Smith & Brown (2018), seed-plant megaphylogeny. https://doi.org/10.1002/ajb2.1019"
  ), file.path(case_dir, "DATA_SOURCES_AND_CITATIONS.md"), useBytes = TRUE)

  write_engine_folder <- function(path, Y) {
    write_csv(Y, file.path(path, "Y.csv"), row.names = TRUE)
    write_csv(XData, file.path(path, "XData.csv"), row.names = TRUE)
    write_csv(TrData, file.path(path, "TrData.csv"), row.names = TRUE)
    write_csv(study, file.path(path, "studyDesign.csv"), row.names = TRUE)
    write_csv(newdata, file.path(path, "newdata.csv"), row.names = TRUE)
    write_csv(newstudy, file.path(path, "prediction_studyDesign.csv"), row.names = TRUE)
    write_csv(taxonomy, file.path(path, "taxonomy.csv"))
    write_csv(trait_audit, file.path(path, "trait_provenance_and_imputation.csv"))
    write_csv(phylo_audit, file.path(path, "phylogeny_match_and_alignment.csv"))
    write_csv(validation, file.path(path, "data_validation_report.csv"))
    ape::write.tree(tree, file = file.path(path, "phylogeny.nwk"))
  }
  write_engine_folder(occurrence_dir, Y_occurrence)
  write_engine_folder(abundance_dir, Y_count)
  write_engine_folder(log1p_dir, Y_log1p)

  list(
    case_dir = case_dir,
    original_dir = original_dir,
    common_dir = common_dir,
    occurrence_dir = occurrence_dir,
    abundance_dir = abundance_dir,
    log1p_dir = log1p_dir,
    Y_count = as.data.frame(Y_count, check.names = FALSE),
    Y_occurrence = as.data.frame(Y_occurrence, check.names = FALSE),
    Y_log1p = as.data.frame(Y_log1p, check.names = FALSE),
    X = XData,
    traits = TrData,
    study = study,
    taxonomy = taxonomy,
    newdata = newdata,
    newstudy = newstudy,
    newick = paste(readLines(file.path(common_dir, "phylogeny_VPhyloMaker2_S3.nwk"), warn = FALSE), collapse = "\n"),
    validation = validation
  )
}

load_runner <- function() {
  path <- file.path(app_dir, "workflow_scripts", "final_case_study_runner.R")
  txt <- readLines(path, warn = FALSE)
  txt <- txt[!grepl("^\\s*if \\(!interactive\\(\\) && identical\\(sys.nframe\\(\\), 0L\\)\\) run_final_cases\\(\\)\\s*$", txt)]
  e <- new.env(parent = globalenv())
  eval(parse(text = paste(txt, collapse = "\n")), envir = e)
  e
}

make_cfg <- function(runner, dat, distribution, samples, transient, thin, chains, do_cv = TRUE) {
  cfg <- runner$configure_hmsc(dat, list(
    hmsc_samples = as.character(samples),
    hmsc_transient = as.character(transient),
    hmsc_thin = as.character(thin),
    hmsc_chains = as.character(chains),
    hmsc_large_samples = as.character(samples),
    hmsc_large_transient = as.character(transient)
  ), large = FALSE, spatial = FALSE, full_bird = FALSE)
  cfg$project_name <- paste0("Whittaker_full_52x75_", gsub("[^A-Za-z0-9]+", "_", distribution))
  cfg$model$distr <- distribution
  # Match the published Chapter 6 model: one linear topographic-moisture
  # gradient. The three requested traits remain in the hierarchical model.
  cfg$model$XFormula <- "~ TMG_z"
  cfg$model$TrFormula <- "~ leaf_CN + woodiness + growth_support"
  cfg$model$use_traits <- TRUE
  cfg$model$use_phylogeny <- TRUE
  cfg$data$phylogeny <- "phylogeny.nwk"
  # Follow the published Chapter 6 model structure. With one unique level per
  # site, a site random level competes strongly with Beta and produced poor
  # species-level PSRF in the full 52 x 75 case.
  cfg$model$random_mode <- "none"
  cfg$model$random_effect_column <- "sample"
  cfg$model$nfMin <- 1L
  cfg$model$nfMax <- 4L
  cfg$model$seed <- 20260918L
  cfg$mcmc$samples <- as.integer(samples)
  cfg$mcmc$transient <- as.integer(transient)
  cfg$mcmc$thin <- as.integer(thin)
  cfg$mcmc$nChains <- as.integer(chains)
  cfg$mcmc$nParallel <- 1L
  cfg$mcmc$verbose <- 0L
  cfg$mcmc$preset <- "Whittaker Chapter 6 structure, extended MCMC"
  # The book code uses sampleMcmc() defaults. Fixed-effect GLM initialisation
  # can be unstable for sparse lognormal-Poisson species, so abundance chains
  # intentionally start from independent random/default values.
  if (identical(distribution, "lognormal poisson")) cfg$mcmc$initPar <- "random"
  cfg$outputs$cv <- isTRUE(do_cv)
  cfg$outputs$diagnostics <- TRUE
  cfg$outputs$parameters <- TRUE
  cfg$outputs$variance_partitioning <- TRUE
  cfg$outputs$omega <- TRUE
  cfg$outputs$gradients <- TRUE
  cfg$outputs$plotting$plotTree <- TRUE
  cfg$outputs$plotting$show.sp.names.beta <- FALSE
  cfg$outputs$plotting$show.sp.names.omega <- FALSE
  cfg$outputs$predictions$nfolds <- 5L
  cfg$outputs$predictions$partition_column <- "sample"
  cfg$outputs$predictions$species.list <- paste(colnames(dat$Y)[seq_len(min(6L, ncol(dat$Y)))], collapse = ",")
  cfg$outputs$predictions$trait.list <- "leaf_CN,woodiness,growth_support"
  cfg$outputs$predictions$env.list <- "TMG_z"
  cfg$outputs$convergence$showOmega <- FALSE
  cfg$outputs$convergence$maxOmega <- 5L
  cfg$outputs$convergence$showAlpha <- FALSE
  cfg
}

run_model <- function(runner, root, prep, response, distribution, case_id,
                      samples, transient, thin, chains, upload_dir, do_cv = TRUE) {
  dat <- list(
    Y = response,
    X = prep[["X"]],
    traits = prep[["traits"]],
    study = prep[["study"]],
    coords = NULL,
    taxonomy = prep[["taxonomy"]],
    newdata = prep[["newdata"]],
    newstudy = prep[["newstudy"]],
    newcoords = NULL,
    newick = prep[["newick"]],
    upload_engine_folder = upload_dir,
    case_input_root = prep[["case_dir"]]
  )
  dat <- runner$make_reference_truth(dat, "whittaker_observed_reference_not_simulated_truth")
  cfg <- make_cfg(runner, dat, distribution, samples, transient, thin, chains, do_cv = do_cv)
  run_warnings <- character()
  result <- withCallingHandlers(
    runner$run_hmsc_with_prediction(
      root, case_id, dat, cfg,
      paste0("Complete Whittaker case: 52 sites, 75 species, three traits, V.PhyloMaker2 S3 tree; ",
             distribution, " response; extended MCMC with ", chains, " chains, ", samples,
             " retained samples, transient=", transient, ", thin=", thin, ".")
    ),
    warning = function(w) {
      run_warnings <<- c(run_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(run_warnings) && dir.exists(result$output_folder[[1]])) {
    warning_file <- file.path(result$output_folder[[1]], "diagnostics", "R_runtime_warnings.txt")
    dir.create(dirname(warning_file), recursive = TRUE, showWarnings = FALSE)
    writeLines(unique(run_warnings), warning_file, useBytes = TRUE)
    result$warnings <- paste(unique(run_warnings), collapse = " | ")
  }
  model_path <- file.path(result$output_folder[[1]], "models", "hmsc_model_main.rds")
  if (identical(result$status[[1]], "fitted") && file.exists(model_path)) {
    audit_env <- new.env(parent = globalenv())
    sys.source(file.path(app_dir, "workflow_scripts", "audit_hmsc_convergence.R"), envir = audit_env)
    tryCatch(
      audit_env$audit_hmsc_convergence(model_path, file.path(result$output_folder[[1]], "diagnostics")),
      error = function(e) writeLines(conditionMessage(e), file.path(result$output_folder[[1]], "diagnostics", "MCMC_parameter_audit_error.txt"), useBytes = TRUE)
    )
  }
  runner$refresh_result(
    result, dat,
    "CASE_STUDY_01_VascularPlants_WhittakerFull",
    "Do species show trait- and phylogeny-structured responses along the topographic moisture gradient?",
    paste0("Primary Hmsc ", distribution, " model")
  )
}

main <- function() {
  cli <- parse_cli()
  prep <- prepare_case()
  cat("PREPARED_CASE_DIR=", normalizePath(prep[["case_dir"]], winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("DATA_DIMENSIONS=52x75\n")
  cat("TRAITS=leaf_CN,woodiness,growth_support\n")
  cat("PHYLOGENY_TIPS=75\n")
  if (tolower(cli[["mode"]]) %in% c("prepare", "prepare-only", "inputs")) return(invisible(prep))

  runner <- load_runner()
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  root <- trimws(cli[["output_root"]] %||% "")
  if (!nzchar(root)) root <- file.path(app_dir, "output", paste0("CASE_STUDY_01_VascularPlants_WhittakerFull_", stamp))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)

  results <- list()
  do_cv <- as_bool(cli[["cv"]], TRUE)
  if (!as_bool(cli[["skip_occurrence"]], FALSE)) {
    results[[length(results) + 1L]] <- run_model(
      runner, root, prep, prep[["Y_occurrence"]], "probit",
      "Whittaker_full_52x75_occurrence_probit",
      as_int(cli[["occurrence_samples"]], 1000L, 1L),
      as_int(cli[["occurrence_transient"]], 2500L, 0L),
      as_int(cli[["occurrence_thin"]], 5L, 1L),
      as_int(cli[["occurrence_chains"]], 4L, 1L),
      prep[["occurrence_dir"]], do_cv = do_cv
    )
  }
  if (!as_bool(cli[["skip_abundance"]], FALSE)) {
    abundance_distribution <- tolower(trimws(as.character(cli[["abundance_distribution"]] %||% "lognormal poisson")))
    use_log1p <- abundance_distribution %in% c("normal", "gaussian", "log1p normal", "log1p_normal")
    abundance_response <- if (use_log1p) prep[["Y_log1p"]] else prep[["Y_count"]]
    abundance_distr <- if (use_log1p) "normal" else "lognormal poisson"
    abundance_case_id <- if (use_log1p) {
      "Whittaker_full_52x75_log1p_abundance_normal_sensitivity"
    } else {
      "Whittaker_full_52x75_abundance_lognormal_poisson"
    }
    abundance_upload_dir <- if (use_log1p) prep[["log1p_dir"]] else prep[["abundance_dir"]]
    results[[length(results) + 1L]] <- run_model(
      runner, root, prep, abundance_response, abundance_distr,
      abundance_case_id,
      as_int(cli[["abundance_samples"]], 1000L, 1L),
      as_int(cli[["abundance_transient"]], 2500L, 0L),
      as_int(cli[["abundance_thin"]], 5L, 1L),
      as_int(cli[["abundance_chains"]], 4L, 1L),
      abundance_upload_dir, do_cv = do_cv
    )
  }
  if (!length(results)) stop("Both occurrence and abundance runs were skipped.", call. = FALSE)
  summary <- do.call(rbind, results)
  rownames(summary) <- NULL
  write_csv(summary, file.path(root, "WHITTAKER_FULL_RUN_SUMMARY.csv"))
  write_csv(prep[["validation"]], file.path(root, "WHITTAKER_FULL_INPUT_VALIDATION.csv"))

  report <- c(
    "# Complete Whittaker vascular-plant analysis",
    "",
    paste0("Generated: ", Sys.time()),
    "",
    "- Input: all 52 sites and all 75 species.",
    "- Measured trait: leaf C:N from the source dataset.",
    "- Added public traits: woodiness and growth support, with auditable hierarchical/phylogenetic imputation.",
    "- Phylogeny: V.PhyloMaker2 scenario S3; all tree tips, response columns and trait rows aligned.",
    paste0("- Models executed: ", paste(summary$case_id, collapse = ", "), "."),
    "- Environmental formula: a linear standardized TMG effect, matching the Chapter 6 model structure.",
    "",
    "## Run summary",
    "",
    paste(capture.output(print(summary[, c("case_id", "status", "file_count", "missing_count", "runtime_seconds", "output_folder", "zip")], row.names = FALSE)), collapse = "\n")
  )
  writeLines(report, file.path(root, "WHITTAKER_FULL_ANALYSIS_REPORT.md"), useBytes = TRUE)

  master_zip_path <- file.path(dirname(root), paste0(basename(root), "_ALL_OUTPUTS.zip"))
  master_zip <- runner$runner_env$make_zip_file(root, master_zip_path)
  sha256 <- if (requireNamespace("openssl", quietly = TRUE)) {
    paste(as.character(openssl::sha256(file(master_zip))), collapse = "")
  } else NA_character_
  writeLines(c(paste0("file=", normalizePath(master_zip, winslash = "/", mustWork = FALSE)),
               paste0("size_bytes=", file.info(master_zip)[["size"]]),
               paste0("sha256=", sha256)),
             file.path(root, "WHITTAKER_FULL_MASTER_ZIP_CHECKSUM.txt"), useBytes = TRUE)
  cat("OUTPUT_ROOT=", normalizePath(root, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("MASTER_ZIP=", normalizePath(master_zip, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("MASTER_ZIP_SHA256=", sha256, "\n", sep = "")
  invisible(summary)
}

if (!interactive() && identical(sys.nframe(), 0L)) main()
