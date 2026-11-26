args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_root <- file.path(project_root, "examples", "00_ACTUAL_UPLOAD_DATA")
if (!dir.exists(data_root)) {
  stop("Cannot find examples/00_ACTUAL_UPLOAD_DATA under project root: ", project_root)
}

read_table_safely <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
}

row_id_cols <- c("", "...1", "id", "ID", "row", "row_id", "unit_id", "site_id", "sample_id", "plot_id")
response_id_cols <- c(row_id_cols, "site", "sample", "plot", "unit")
predictor_id_cols <- c(row_id_cols, "site", "sample", "plot", "unit")
trait_id_cols <- c(row_id_cols, "species", "species_id", "taxon", "taxon_id")
design_id_cols <- c("", "...1", "id", "ID", "row", "row_id")
coord_id_cols <- c(row_id_cols, "site", "sample", "plot", "unit")
general_id_cols <- unique(c(response_id_cols, trait_id_cols))

id_cols_for_role <- function(role) {
  switch(
    role,
    species = response_id_cols,
    response = response_id_cols,
    predictors = predictor_id_cols,
    traits = trait_id_cols,
    design = design_id_cols,
    coords = coord_id_cols,
    matrix = trait_id_cols,
    columns = predictor_id_cols,
    general_id_cols
  )
}

count_data_columns <- function(dat, id_cols = general_id_cols) {
  if (is.null(dat)) return(NA_integer_)
  cols <- names(dat)
  sum(!tolower(cols) %in% tolower(id_cols))
}

dim_text <- function(path, role = "columns") {
  dat <- read_table_safely(path)
  if (is.null(dat)) return("-")
  id_cols <- id_cols_for_role(role)
  paste0(nrow(dat), " x ", count_data_columns(dat, id_cols), " ", role)
}

has_file <- function(...) file.exists(file.path(...))

data_columns <- function(path, role = "columns") {
  dat <- read_table_safely(path)
  if (is.null(dat)) return(character(0))
  cols <- names(dat)
  cols[!tolower(cols) %in% tolower(id_cols_for_role(role))]
}

collapse_or_dash <- function(x) {
  x <- x[nzchar(x)]
  if (length(x) == 0) "-" else paste(x, collapse = ", ")
}

extract_case_name <- function(folder_name) {
  sub("^[0-9]+_[^_]+(?:_[^_]+)*__", "", folder_name)
}

case_dirs <- list.dirs(data_root, full.names = TRUE, recursive = FALSE)
case_dirs <- case_dirs[grepl("^[0-9]+_", basename(case_dirs))]
case_dirs <- case_dirs[order(basename(case_dirs))]

taxon_map <- c(
  "01_vascular_plants__ch06_plant_traits_whittaker" = "Vascular plants",
  "02_deadwood_fungi__ch07_deadwood_fungi" = "Deadwood fungi",
  "03_birds__ch11_finnish_birds" = "Birds",
  "04_bryophytes__tibet_bryophyte_30site_spatial_example" = "Bryophytes",
  "05_freshwater_macroinvertebrates__cestes_freshwater_macroinvertebrates_doledec1996" = "Freshwater macroinvertebrates",
  "06_aquatic_ground_beetles__cestes_aquatic_or_ground_beetles_ribera2001" = "Aquatic or ground beetles",
  "07_spiders__cestes_spiders_gibb2015" = "Spiders",
  "08_marine_fish__mediterranean_fish_efh_moullec" = "Marine fish",
  "09_mammals__phylacine_mammals_range_derived" = "Mammals",
  "10_synthetic_benchmark__universal_benchmark" = "Synthetic benchmark"
)

publication_case_map <- c(
  "01_vascular_plants__ch06_plant_traits_whittaker" = "Book Ch. 6 plant trait-environment case",
  "02_deadwood_fungi__ch07_deadwood_fungi" = "Book Ch. 7 deadwood fungal community case",
  "03_birds__ch11_finnish_birds" = "Book Ch. 11 Finnish bird prediction case",
  "04_bryophytes__tibet_bryophyte_30site_spatial_example" = "Tibet bryophyte spatial community case",
  "05_freshwater_macroinvertebrates__cestes_freshwater_macroinvertebrates_doledec1996" = "Freshwater macroinvertebrate community case",
  "06_aquatic_ground_beetles__cestes_aquatic_or_ground_beetles_ribera2001" = "Aquatic/ground beetle trait case",
  "07_spiders__cestes_spiders_gibb2015" = "Spider community trait case",
  "08_marine_fish__mediterranean_fish_efh_moullec" = "Mediterranean fish EFH community case",
  "09_mammals__phylacine_mammals_range_derived" = "PHYLACINE mammal range-derived case",
  "10_synthetic_benchmark__universal_benchmark" = "Universal synthetic benchmark"
)

rows <- lapply(case_dirs, function(case_dir) {
  folder <- basename(case_dir)
  y_path <- file.path(case_dir, "Hmsc", "Y.csv")
  y <- read_table_safely(y_path)
  sites <- if (is.null(y)) NA_integer_ else nrow(y)
  species <- if (is.null(y)) NA_integer_ else count_data_columns(y, response_id_cols)
  sp_y <- read_table_safely(file.path(case_dir, "spOccupancy", "y.csv"))
  det_cols <- if (is.null(sp_y)) NA_integer_ else count_data_columns(sp_y, response_id_cols)
  visits <- if (is.na(det_cols) || is.na(species) || species == 0) NA_integer_ else as.integer(det_cols / species)
  has_normal <- has_file(case_dir, "GJAM", "Y.csv") && !grepl("bryophyte", folder)
  response_text <- paste0(
    "community Y: ", sites, " x ", species,
    "; count Y: ", if (has_file(case_dir, "boral", "Y.csv")) paste0(sites, " x ", species) else "-",
    "; normal Y: ", if (has_normal) paste0(sites, " x ", species) else "-"
  )
  data.frame(
    taxonomic_group = unname(taxon_map[folder]),
    retained_case = unname(publication_case_map[folder]),
    case_folder = folder,
    sites = sites,
    species = species,
    visits = ifelse(is.na(visits), "-", visits),
    core_response_matrices = response_text,
    environmental_covariates = dim_text(file.path(case_dir, "Hmsc", "XData.csv"), "predictors"),
    traits = dim_text(file.path(case_dir, "Hmsc", "TrData.csv"), "traits"),
    spatial_and_design = paste(
      paste0("studyDesign ", dim_text(file.path(case_dir, "Hmsc", "studyDesign.csv"), "design")),
      paste0("coordinates ", dim_text(file.path(case_dir, "Hmsc", "coordinates.csv"), "coords")),
      sep = "; "
    ),
    phylogeny = if (has_file(case_dir, "Hmsc", "phylogeny.nwk")) paste0(species, " tips; phylo_cov ", species, " x ", species) else "-",
    prediction_uploads = paste(
      paste0("newdata ", dim_text(file.path(case_dir, "Hmsc", "newdata.csv"), "predictors")),
      paste0("newcoords ", dim_text(file.path(case_dir, "Hmsc", "newcoords.csv"), "coords")),
      sep = "; "
    ),
    spoccupancy_detection = if (is.null(sp_y)) "-" else paste0("y.csv ", sites, " x ", det_cols, " detection columns = ", species, " species x ", visits, " visits; det.covs ", dim_text(file.path(case_dir, "spOccupancy", "det.covs.csv"), "columns")),
    advanced_branch_files = "typeNames/specByTrait for GJAM; trials/long format for jSDM; folds for CV; trial.size/offset/row.ids/ranef.ids/distmat for boral",
    recommended_use = if (grepl("ch11", folder)) {
      "Best for prediction-map demonstrations; upload Hmsc/Hmsc-HPC prediction files and spatial files."
    } else if (grepl("bryophyte", folder)) {
      "Best real spatial high-diversity example; use Hmsc/Hmsc-HPC/spOccupancy/sjSDM spatial branches carefully."
    } else if (grepl("universal", folder)) {
      "Best software benchmark; use for all-model smoke tests and unified comparison."
    } else {
      "Best as a compact cross-engine upload example with traits, spatial files and prediction files."
    },
    stringsAsFactors = FALSE
  )
})

summary_table <- do.call(rbind, rows)

case_short_map <- c(
  "01_vascular_plants__ch06_plant_traits_whittaker" = "Plant trait-environment",
  "02_deadwood_fungi__ch07_deadwood_fungi" = "Deadwood fungi",
  "03_birds__ch11_finnish_birds" = "Finnish birds",
  "04_bryophytes__tibet_bryophyte_30site_spatial_example" = "Tibetan bryophytes",
  "05_freshwater_macroinvertebrates__cestes_freshwater_macroinvertebrates_doledec1996" = "Doledec et al. case",
  "06_aquatic_ground_beetles__cestes_aquatic_or_ground_beetles_ribera2001" = "Ribera et al. case",
  "07_spiders__cestes_spiders_gibb2015" = "Gibb et al. case",
  "08_marine_fish__mediterranean_fish_efh_moullec" = "Mediterranean fish",
  "09_mammals__phylacine_mammals_range_derived" = "PHYLACINE-derived case",
  "10_synthetic_benchmark__universal_benchmark" = "Universal benchmark"
)

primary_use_map <- c(
  "01_vascular_plants__ch06_plant_traits_whittaker" = "Cross-engine",
  "02_deadwood_fungi__ch07_deadwood_fungi" = "Cross-engine",
  "03_birds__ch11_finnish_birds" = "Prediction",
  "04_bryophytes__tibet_bryophyte_30site_spatial_example" = "Spatial / high diversity",
  "05_freshwater_macroinvertebrates__cestes_freshwater_macroinvertebrates_doledec1996" = "Cross-engine",
  "06_aquatic_ground_beetles__cestes_aquatic_or_ground_beetles_ribera2001" = "Trait-rich",
  "07_spiders__cestes_spiders_gibb2015" = "Trait-rich",
  "08_marine_fish__mediterranean_fish_efh_moullec" = "High diversity",
  "09_mammals__phylacine_mammals_range_derived" = "Cross-engine",
  "10_synthetic_benchmark__universal_benchmark" = "Software validation"
)

concise_rows <- lapply(case_dirs, function(case_dir) {
  folder <- basename(case_dir)
  y <- read_table_safely(file.path(case_dir, "Hmsc", "Y.csv"))
  sites <- if (is.null(y)) NA_integer_ else nrow(y)
  species <- if (is.null(y)) NA_integer_ else count_data_columns(y, response_id_cols)
  x <- read_table_safely(file.path(case_dir, "Hmsc", "XData.csv"))
  traits <- read_table_safely(file.path(case_dir, "Hmsc", "TrData.csv"))
  pred <- read_table_safely(file.path(case_dir, "Hmsc", "newdata.csv"))
  has_normal <- has_file(case_dir, "GJAM", "Y.csv") && !grepl("bryophyte", folder)
  env_n <- if (is.null(x)) NA_integer_ else count_data_columns(x, predictor_id_cols)
  trait_n <- if (is.null(traits)) NA_integer_ else count_data_columns(traits, trait_id_cols)
  pred_n <- if (is.null(pred)) NA_integer_ else nrow(pred)
  design_cols <- data_columns(file.path(case_dir, "Hmsc", "studyDesign.csv"), "design")
  coord_cols <- data_columns(file.path(case_dir, "Hmsc", "coordinates.csv"), "coords")
  coord_type <- if (any(grepl("lon|lat|longitude|latitude", coord_cols, ignore.case = TRUE))) {
    "geographic or projected site coordinates"
  } else {
    "projected/local 2-D site coordinates"
  }
  data.frame(
    Taxonomic_group = unname(taxon_map[folder]),
    Case = unname(case_short_map[folder]),
    Sites = sites,
    Species = species,
    Response_uploads = if (has_normal) {
      "B/C/G via Y.csv; detection via spOccupancy/y.csv"
    } else {
      "B/C via Y.csv; detection via spOccupancy/y.csv"
    },
    Env_uploads = paste0(env_n, " vars: XData.csv, env.csv, occ.covs.csv"),
    Trait_uploads = if (is.null(traits)) "-" else paste0(trait_n, " traits: TrData.csv, traits.csv, trait_data.csv, specByTrait.csv"),
    Spatial_design_type_and_uploads = paste0(
      "Grouping/random-effect factor table: studyDesign.csv (",
      collapse_or_dash(design_cols),
      "); ",
      coord_type,
      ": coordinates.csv, coords.csv, spatial.csv (",
      collapse_or_dash(coord_cols),
      "); boral row.ids.csv/ranef.ids.csv when row or random effects are enabled"
    ),
    Phylogeny_type_and_uploads = if (has_file(case_dir, "Hmsc", "phylogeny.nwk")) {
      paste0(
        "Newick tree with branch lengths: phylogeny.nwk/phylo_tree.nwk (",
        species,
        " tips); taxonomy rank table: taxonomy.csv (kingdom/phylum/class/order/family/genus/species). Species-by-species correlation is derived internally and saved as output."
      )
    } else {
      "-"
    },
    Prediction_uploads = if (is.null(pred)) {
      "-"
    } else {
      paste0(pred_n, " sites: newdata.csv, newcoords.csv, prediction_studyDesign.csv, new_spatial.csv, prediction_ids.csv")
    },
    Other_original_inputs = "typeNames.csv, traitTypes.csv, trial.size.csv, offset.csv, distmat.csv, folds.csv, det.covs.csv",
    Primary_use = unname(primary_use_map[folder]),
    stringsAsFactors = FALSE
  )
})

concise_case_table <- do.call(rbind, concise_rows)

workflow_guide <- data.frame(
  workflow = c("Hmsc", "Hmsc-HPC", "jSDM", "GJAM", "spOccupancy", "sjSDM", "boral"),
  minimum_first_run_uploads = c(
    "Y.csv + XData.csv",
    "Y.csv + XData.csv",
    "Y.csv + XData.csv",
    "Y.csv + XData.csv + typeNames.csv",
    "y.csv + occ.covs.csv + det.covs.csv",
    "Y.csv + env.csv",
    "Y.csv + XData.csv"
  ),
  maximum_branch_uploads = c(
    "Y.csv, XData.csv, TrData.csv, studyDesign.csv, coordinates.csv, phylogeny.nwk, newdata.csv, newcoords.csv, prediction_studyDesign.csv",
    "Y.csv, XData.csv, traits.csv, studyDesign.csv, coordinates.csv, taxonomy.csv or phylo_tree.nwk, newdata.csv",
    "Y.csv, XData.csv, trait_data.csv, trials.csv, long_format.csv, newdata.csv, prediction_ids.csv",
    "Y.csv, XData.csv, typeNames.csv, specByTrait.csv, traitTypes.csv, newdata.csv",
    "y.csv, occ.covs.csv, det.covs.csv, coords.csv, species.csv, integrated_sources.csv, newdata.csv, newcoords.csv, folds.csv",
    "Y.csv, env.csv, spatial.csv, traits.csv, newdata.csv, new_spatial.csv, species_groups.csv, folds.csv",
    "Y.csv, XData.csv, traits.csv, trial.size.csv, offset.csv, row.ids.csv, ranef.ids.csv, distmat.csv, newdata.csv"
  ),
  when_to_use_optional_files = c(
    "Traits, phylogeny, sample random effects, spatial random effects, or prediction are enabled.",
    "CPU pyhmsc/HDF5 branch needs traits, Newick or taxonomy-derived phylogeny, random or spatial effects, or prediction.",
    "Binomial trials, trait summaries, long-format checks, or new-site prediction are enabled.",
    "Mixed response scales, trait response summaries, or prediction are enabled. typeNames is central.",
    "Detection modelling, spatial occupancy, integrated sources, folds, or new-site prediction are enabled.",
    "Spatial neural terms, traits, folds, groups, or new-site prediction are enabled.",
    "Binomial trial sizes, offsets, row/random effects, distance matrices, traits, or prediction are enabled."
  ),
  note = c(
    "Do not upload all optional files for the minimal run; choose the matching branch first.",
    "Not a full replacement for Hmsc-R; CPU workflow is guarded and should report model_defined/check_failed honestly when native support is absent.",
    "Use occurrence/binomial-compatible Y for the main cross-engine case.",
    "Do not compare GJAM observation-scale parameters directly with Hmsc/jSDM link-scale parameters.",
    "Detection effects must be interpreted separately from occurrence effects.",
    "PyTorch/reticulate environment is the main dependency risk.",
    "JAGS/rjags/R2jags is the main dependency risk."
  ),
  stringsAsFactors = FALSE
)

safe_write_csv <- function(dat, path) {
  tryCatch(
    {
      utils::write.csv(dat, path, row.names = FALSE, fileEncoding = "UTF-8")
      message("WROTE: ", path)
      path
    },
    error = function(e) {
      fallback <- sub("\\.csv$", paste0("_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"), path)
      utils::write.csv(dat, fallback, row.names = FALSE, fileEncoding = "UTF-8")
      message("LOCKED, WROTE FALLBACK: ", fallback)
      fallback
    }
  )
}

safe_write_csv(summary_table, file.path(data_root, "ACTUAL_UPLOAD_SIMPLIFIED_PUBLICATION_TABLE.csv"))
safe_write_csv(workflow_guide, file.path(data_root, "ACTUAL_UPLOAD_WORKFLOW_MIN_MAX_GUIDE.csv"))
safe_write_csv(concise_case_table, file.path(data_root, "ACTUAL_UPLOAD_CONCISE_SOFTWARE_INPUT_TABLE.csv"))

escape_md <- function(x) {
  x <- gsub("\\|", "\\\\|", x)
  x <- gsub("\n", "<br>", x, fixed = TRUE)
  x
}

write_markdown_table <- function(dat, path, title, intro) {
  lines <- c(
    paste0("# ", title),
    "",
    intro,
    "",
    paste0("| ", paste(names(dat), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |")
  )
  for (i in seq_len(nrow(dat))) {
    lines <- c(lines, paste0("| ", paste(vapply(dat[i, ], function(x) escape_md(as.character(x)), character(1)), collapse = " | "), " |"))
  }
  writeLines(lines, path, useBytes = TRUE)
}

write_markdown_table(
  summary_table,
  file.path(data_root, "ACTUAL_UPLOAD_SIMPLIFIED_PUBLICATION_TABLE.md"),
  "Simplified Actual Upload Data Table",
  "This table is designed for papers, manuals and reviewers. It lists the real JSDM Studio upload datasets by retained taxonomic case. Optional branch files are summarized instead of being treated as required inputs."
)

write_markdown_table(
  workflow_guide,
  file.path(data_root, "ACTUAL_UPLOAD_WORKFLOW_MIN_MAX_GUIDE.md"),
  "Workflow Minimum and Maximum Upload Guide",
  "Use the minimum files for a first smoke test. Use maximum branch uploads only after enabling the corresponding model branch in the GUI."
)

write_markdown_table(
  concise_case_table,
  file.path(data_root, "ACTUAL_UPLOAD_CONCISE_SOFTWARE_INPUT_TABLE.md"),
  "Concise Software Input Table",
  "B = binary occurrence, C = count, G = Gaussian/continuous. The upload columns list the actual JSDM Studio file names; optional files should be used only when the matching workflow branch is enabled."
)

readme_path <- file.path(data_root, "README_ACTUAL_UPLOAD_DATA.md")
readme_lines <- c(
  "# Actual Upload Data for JSDM Studio",
  "",
  "This folder contains cleaned, deduplicated example data that can be uploaded directly into JSDM Studio workflow panels.",
  "",
  "Most files are not required for every run. Use `ACTUAL_UPLOAD_WORKFLOW_MIN_MAX_GUIDE.md` to decide which files belong to a minimal first run and which files belong only to optional branches.",
  "",
  "Main reviewer-facing tables:",
  "",
  "- `ACTUAL_UPLOAD_SIMPLIFIED_PUBLICATION_TABLE.md`: compact case-by-case table for manuscripts and manuals.",
  "- `ACTUAL_UPLOAD_SIMPLIFIED_PUBLICATION_TABLE.csv`: CSV version of the same table.",
  "- `ACTUAL_UPLOAD_CONCISE_SOFTWARE_INPUT_TABLE.md`: concise paper-style table with actual GUI upload file names.",
  "- `ACTUAL_UPLOAD_CONCISE_SOFTWARE_INPUT_TABLE.csv`: CSV version of the concise paper-style table.",
  "- `ACTUAL_UPLOAD_WORKFLOW_MIN_MAX_GUIDE.md`: workflow-specific minimum and maximum upload guide.",
  "- `ACTUAL_UPLOAD_WORKFLOW_MIN_MAX_GUIDE.csv`: CSV version of the workflow guide.",
  "- `ACTUAL_UPLOAD_PUBLICATION_TABLE.csv`: full technical table with all per-workflow upload details.",
  "- `ACTUAL_UPLOAD_DATA_INDEX.csv`: long technical inventory.",
  "",
  "Folder convention:",
  "",
  "`<case>/<workflow>/UPLOAD_MAP.csv` tells which file in that workflow folder should be selected in each JSDM Studio upload box.",
  "",
  "For example, the Hmsc workflow uses `Y.csv`, `XData.csv`, optional `TrData.csv`, optional `studyDesign.csv`, optional `coordinates.csv`, optional `phylogeny.nwk`, and optional prediction files. The spOccupancy workflow uses detection-nondetection files such as `y.csv`, `occ.covs.csv` and `det.covs.csv`; it should not be treated as the same upload contract as Hmsc."
)
writeLines(readme_lines, readme_path, useBytes = TRUE)

message("WROTE:")
message(file.path(data_root, "ACTUAL_UPLOAD_SIMPLIFIED_PUBLICATION_TABLE.md"))
message(file.path(data_root, "ACTUAL_UPLOAD_WORKFLOW_MIN_MAX_GUIDE.md"))
message(file.path(data_root, "ACTUAL_UPLOAD_CONCISE_SOFTWARE_INPUT_TABLE.md"))
