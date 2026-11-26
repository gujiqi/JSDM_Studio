options(stringsAsFactors = FALSE)

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_arg <- args0[startsWith(args0, file_arg)]
script_path <- if (length(script_arg)) sub(file_arg, "", script_arg[[1]]) else "workflow_scripts/summarize_public_case_data_types.R"
app_dir <- normalizePath(file.path(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) {
  app_dir <- "C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio"
}
root <- file.path(app_dir, "examples", "00_UPLOAD_READY_public_integrated_cases")
names_path <- file.path(root, "PUBLIC_CASE_NAMES.csv")
if (!file.exists(names_path)) stop("Missing PUBLIC_CASE_NAMES.csv", call. = FALSE)

read_table <- function(path, row.names = 1) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(row.names) && ncol(x) >= row.names) {
    rn <- x[[row.names]]
    if (!anyDuplicated(rn) && all(nzchar(as.character(rn)))) {
      rownames(x) <- as.character(rn)
      x <- x[-row.names]
    }
  }
  x
}

is_numeric_like_col <- function(z) {
  if (is.numeric(z) || is.integer(z)) return(TRUE)
  zz <- as.character(z)
  zz <- zz[!is.na(zz) & nzchar(trimws(zz))]
  if (!length(zz)) return(FALSE)
  mean(!is.na(suppressWarnings(as.numeric(zz)))) > 0.95
}

has_numeric_like <- function(x) any(vapply(x, is_numeric_like_col, logical(1)))
has_categorical_like <- function(x) any(!vapply(x, is_numeric_like_col, logical(1)))

core_files <- c(
  "Y_occurrence.csv", "Y_count.csv", "Y_normal.csv", "XData.csv", "traits.csv",
  "studyDesign.csv", "coordinates.csv", "phylogeny.nwk", "taxonomy.csv",
  "newdata.csv", "newcoords.csv", "folds.csv", "trial.size.csv", "offset.csv",
  "row.ids.csv", "ranef.ids.csv", "distmat.csv", "spOccupancy_y_detection.csv",
  "occ.covs.csv", "det.covs.csv", "model_formulas.txt", "UPLOAD_GUIDE.csv",
  "dataset_manifest.csv", "README_UPLOAD_FIRST.md", "SOURCE_CITATION.md"
)

nm <- read.csv(names_path, check.names = FALSE, stringsAsFactors = FALSE)
rows <- lapply(seq_len(nrow(nm)), function(i) {
  id <- nm$case_id[[i]]
  d <- file.path(root, id, "00_common_inputs")
  y <- read_table(file.path(d, "Y_occurrence.csv"))
  x <- read_table(file.path(d, "XData.csv"))
  tr <- read_table(file.path(d, "traits.csv"))
  co <- read_table(file.path(d, "coordinates.csv"))
  det <- read_table(file.path(d, "spOccupancy_y_detection.csv"))
  nd <- read_table(file.path(d, "newdata.csv"))
  nc <- read_table(file.path(d, "newcoords.csv"))
  tax_path <- file.path(d, "taxonomy.csv")
  tax <- if (file.exists(tax_path)) read_table(tax_path, row.names = NULL) else data.frame()
  study <- read_table(file.path(d, "studyDesign.csv"))
  distmat <- read_table(file.path(d, "distmat.csv"))
  data.frame(
    case_id = id,
    short_id = nm$short_id[[i]],
    display_name = nm$display_name[[i]],
    sites = nrow(y),
    species = ncol(y),
    has_occurrence_Y = file.exists(file.path(d, "Y_occurrence.csv")),
    has_count_Y = file.exists(file.path(d, "Y_count.csv")),
    has_normal_Y = file.exists(file.path(d, "Y_normal.csv")),
    env_cols = ncol(x),
    env_has_numeric = has_numeric_like(x),
    env_has_categorical = has_categorical_like(x),
    trait_cols = ncol(tr),
    traits_has_numeric = has_numeric_like(tr),
    traits_has_categorical = has_categorical_like(tr),
    study_cols = ncol(study),
    study_has_sample = "sample" %in% names(study),
    coord_cols = ncol(co),
    coords_have_x_y = all(c("x", "y") %in% names(co)),
    coords_have_lon_lat = all(c("longitude", "latitude") %in% names(co)),
    has_phylogeny_nwk = file.exists(file.path(d, "phylogeny.nwk")) && file.info(file.path(d, "phylogeny.nwk"))$size > 0,
    has_taxonomy_csv = file.exists(file.path(d, "taxonomy.csv")),
    taxonomy_rank_cols = sum(tolower(names(tax)) %in% c("kingdom", "phylum", "division", "class", "order", "family", "genus", "species")),
    taxonomy_matches_species = any(tolower(names(tax)) %in% c("species", "species_name", "species_id", "taxon", "taxon_id", "taxon_name")) &&
      nrow(tax) == ncol(y),
    has_derived_phylogeny_reference = file.exists(file.path(d, "truth", "derived_phylogeny_correlation_reference.csv")),
    has_newdata = nrow(nd) > 0,
    has_newcoords = nrow(nc) > 0,
    detection_columns = ncol(det),
    inferred_detection_visits = if (ncol(y) > 0) ncol(det) / ncol(y) else NA_real_,
    has_trial_size = file.exists(file.path(d, "trial.size.csv")),
    has_offset = file.exists(file.path(d, "offset.csv")),
    has_row_ids = file.exists(file.path(d, "row.ids.csv")),
    has_ranef_ids = file.exists(file.path(d, "ranef.ids.csv")),
    distmat_square = nrow(distmat) == ncol(distmat) && nrow(distmat) == nrow(y),
    has_folds = file.exists(file.path(d, "folds.csv")),
    all_core_files = all(file.exists(file.path(d, core_files))),
    stringsAsFactors = FALSE
  )
})

out <- do.call(rbind, rows)
write.csv(out, file.path(root, "PUBLIC_CASE_DATA_TYPE_COVERAGE.csv"), row.names = FALSE)
print(out[, c("short_id", "sites", "species", "env_cols", "env_has_numeric",
              "env_has_categorical", "trait_cols", "traits_has_numeric",
              "traits_has_categorical", "coords_have_x_y", "coords_have_lon_lat",
              "has_phylogeny_nwk", "phylo_cov_matches_species",
              "inferred_detection_visits", "has_newdata", "has_trial_size",
              "distmat_square", "all_core_files")], row.names = FALSE)
