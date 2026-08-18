options(stringsAsFactors = FALSE)

script_path <- {
  args0 <- commandArgs(trailingOnly = FALSE)
  hit <- args0[startsWith(args0, "--file=")]
  if (length(hit)) normalizePath(sub("^--file=", "", hit[[1]]), winslash = "/", mustWork = FALSE) else
    normalizePath("workflow_scripts/make_book_chapter_upload_ready_examples.R", winslash = "/", mustWork = FALSE)
}
app_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) stop("Cannot locate JSDM Studio app.R.", call. = FALSE)

case_root <- file.path(app_dir, "examples", "book_chapter_benchmarks")
out_root <- file.path(app_dir, "examples", "00_UPLOAD_READY_book_chapter_cases")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

cases <- c(
  ch06_plant_traits_whittaker = "Chapter 6 plant traits and Whittaker topographic moisture gradient",
  ch07_deadwood_fungi = "Chapter 7 dead wood-inhabiting fungi sequencing data",
  ch11_finnish_birds = "Chapter 11 Finnish birds community data"
)

read_case_csv <- function(case_id, name, row.names = 1) {
  read.csv(file.path(case_root, case_id, name), row.names = row.names,
           check.names = FALSE, stringsAsFactors = FALSE)
}

write_csv <- function(x, path, row.names = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = row.names)
}

copy_file <- function(src, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  file.copy(src, dest, overwrite = TRUE)
}

make_long_format <- function(Y, X) {
  grid <- as.data.frame(as.table(as.matrix(Y)), stringsAsFactors = FALSE)
  names(grid) <- c("site", "species", "presence")
  grid$site <- as.character(grid$site)
  grid$species <- as.character(grid$species)
  X2 <- X
  X2$site <- rownames(X2)
  merge(grid, X2, by = "site", all.x = TRUE, sort = FALSE)
}

write_engine_guide <- function(case_dir, case_id, title) {
  guide <- data.frame(
    engine = c(
      rep("Hmsc", 6), rep("Hmsc-HPC", 7), rep("jSDM", 7), rep("GJAM", 9),
      rep("spOccupancy", 9), rep("sjSDM", 5), rep("boral", 9)
    ),
    ui_upload_slot = c(
      "Y.csv - response matrix", "XData.csv - environmental predictors", "TrData.csv - traits / response attributes",
      "studyDesign.csv - grouping / random-effect design", "coordinates.csv - spatial coordinates", "phylogeny / taxonomy file",
      "Y.csv", "XData.csv", "traits.csv", "studyDesign.csv", "coordinates.csv", "phylo_cov.csv", "newdata.csv",
      "Y.csv", "XData.csv", "trait_data.csv", "long_format.csv", "trials.csv", "newdata.csv", "prediction_ids.csv",
      "Y.csv", "XData.csv", "typeNames.csv", "specByTrait.csv", "traitTypes.csv", "effort.csv", "newdata.csv", "holdoutIndex.csv", "censor.csv",
      "y.csv", "occ.covs.csv", "det.covs.csv", "coords.csv", "species.csv", "integrated_sources.csv", "newdata.csv", "newcoords.csv", "folds.csv",
      "Y.csv", "env.csv", "spatial.csv", "traits.csv", "newdata.csv",
      "Y.csv", "XData.csv", "traits.csv", "row.ids.csv", "ranef.ids.csv", "distmat.csv", "offset.csv", "newdata.csv", "trial.size.csv"
    ),
    file_to_upload = c(
      "01_Hmsc/Y.csv", "01_Hmsc/XData.csv", "01_Hmsc/TrData.csv", "01_Hmsc/studyDesign.csv", "01_Hmsc/coordinates.csv", "01_Hmsc/phylogeny.nwk",
      "02_Hmsc-HPC/Y.csv", "02_Hmsc-HPC/XData.csv", "02_Hmsc-HPC/traits.csv", "02_Hmsc-HPC/studyDesign.csv", "02_Hmsc-HPC/coordinates.csv", "02_Hmsc-HPC/phylo_cov.csv", "02_Hmsc-HPC/newdata.csv",
      "03_jSDM/Y.csv", "03_jSDM/XData.csv", "03_jSDM/trait_data.csv", "03_jSDM/long_format.csv", "03_jSDM/trials.csv", "03_jSDM/newdata.csv", "03_jSDM/prediction_ids.csv",
      "04_GJAM/Y.csv", "04_GJAM/XData.csv", "04_GJAM/typeNames.csv", "04_GJAM/specByTrait.csv", "04_GJAM/traitTypes.csv", "04_GJAM/effort.csv", "04_GJAM/newdata.csv", "04_GJAM/holdoutIndex.csv", "04_GJAM/censor.csv",
      "05_spOccupancy/y.csv", "05_spOccupancy/occ.covs.csv", "05_spOccupancy/det.covs.csv", "05_spOccupancy/coords.csv", "05_spOccupancy/species.csv", "05_spOccupancy/integrated_sources.csv", "05_spOccupancy/newdata.csv", "05_spOccupancy/newcoords.csv", "05_spOccupancy/folds.csv",
      "06_sjSDM/Y.csv", "06_sjSDM/env.csv", "06_sjSDM/spatial.csv", "06_sjSDM/traits.csv", "06_sjSDM/newdata.csv",
      "07_boral/Y.csv", "07_boral/XData.csv", "07_boral/traits.csv", "07_boral/row.ids.csv", "07_boral/ranef.ids.csv", "07_boral/distmat.csv", "07_boral/offset.csv", "07_boral/newdata.csv", "07_boral/trial.size.csv"
    ),
    required_for_default = c(
      TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
      TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
      TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
      TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      TRUE, TRUE, FALSE, FALSE, FALSE,
      TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
    ),
    note = c(
      "Use Y occurrence for probit.", "Predictors are standardized and include substrate factor.", "Enable Use TrData / traits before using this file.",
      "Use random effect mode = sample and grouping column = plot.", "Use only for spatial random effect modes.", "Use only if phylogeny is enabled.",
      "Use occurrence/probit-compatible Y.", "Predictors for Hmsc-HPC.", "Optional traits.", "Optional random-effect design.", "Optional spatial input.", "Optional phylogenetic covariance.", "Optional prediction data.",
      "Default binomial_probit uses Y/X; traits optional.", "Site predictors.", "Enable traits if used.", "Only for long-format branch.", "Only for binomial_logit branch.", "Optional prediction data.", "Optional prediction ID filter.",
      "Use PA typeNames.", "Predictors.", "Upload for response types.", "Optional trait conversion.", "Optional trait types.", "Optional effort/offset information.", "Optional prediction data.", "Optional holdout rows.", "Leave unused unless censoring is enabled.",
      "Replicated detection/non-detection matrix.", "Occurrence covariates.", "Detection covariates.", "Coordinates.", "Species metadata.", "Optional integrated model metadata.", "Optional prediction covariates.", "Optional prediction coordinates.", "Optional folds.",
      "Occurrence response.", "Environmental module.", "Spatial module uses coordinates as simple spatial predictors.", "Species traits metadata.", "Optional prediction data.",
      "Occurrence/binomial response.", "Predictors.", "Optional traits.", "Optional row IDs.", "Optional random-effect IDs.", "Optional distance matrix.", "Optional offset.", "Optional prediction data.", "Binomial trial size matrix."
    ),
    stringsAsFactors = FALSE
  )
  write.csv(guide, file.path(case_dir, "UPLOAD_GUIDE.csv"), row.names = FALSE)
  writeLines(c(
    paste0("# ", title),
    "",
    "Open the engine folder that matches the JSDM Studio workflow panel.",
    "Upload the files using the exact GUI slot names listed in UPLOAD_GUIDE.csv.",
    "",
    "Fastest default tests:",
    "- Hmsc: upload only Y.csv and XData.csv, keep probit, random level none.",
    "- Hmsc-HPC: upload Y.csv and XData.csv, CPU quick settings.",
    "- jSDM: upload Y.csv and XData.csv; trait_data.csv is optional.",
    "- GJAM: upload Y.csv, XData.csv and typeNames.csv; keep PA response type.",
    "- spOccupancy: upload y.csv, occ.covs.csv and det.covs.csv; coords optional for spatial models.",
    "- sjSDM: upload Y.csv and env.csv; spatial.csv optional.",
    "- boral: upload Y.csv and XData.csv; trial.size.csv optional for binomial.",
    "",
    "For full benchmark coverage, upload every file in that engine folder to the matching optional slot."
  ), file.path(case_dir, "README_UPLOAD_FIRST.md"), useBytes = TRUE)
}

for (case_id in names(cases)) {
  src <- file.path(case_root, case_id)
  if (!dir.exists(src)) stop("Missing generated book case: ", src, call. = FALSE)
  case_dir <- file.path(out_root, case_id)
  dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)

  Y <- read_case_csv(case_id, "Y_occurrence.csv")
  X <- read_case_csv(case_id, "XData.csv")
  Tr <- read_case_csv(case_id, "traits.csv")
  study <- read_case_csv(case_id, "studyDesign.csv")
  coords <- read_case_csv(case_id, "coordinates.csv")
  newdata <- read_case_csv(case_id, "newdata.csv")
  newcoords <- read_case_csv(case_id, "newcoords.csv")
  trial <- read_case_csv(case_id, "trial.size.csv")
  species <- colnames(Y)
  sites <- rownames(Y)

  common <- file.path(case_dir, "00_common_inputs")
  dir.create(common, recursive = TRUE, showWarnings = FALSE)
  for (f in c("Y_occurrence.csv", "Y_count.csv", "Y_normal.csv", "XData.csv", "traits.csv",
              "studyDesign.csv", "coordinates.csv", "phylogeny.nwk", "phylo_cov.csv",
              "newdata.csv", "newcoords.csv", "folds.csv", "trial.size.csv", "offset.csv",
              "row.ids.csv", "ranef.ids.csv", "distmat.csv", "spOccupancy_y_detection.csv",
              "occ.covs.csv", "det.covs.csv")) {
    copy_file(file.path(src, f), file.path(common, f))
  }

  hmsc <- file.path(case_dir, "01_Hmsc")
  write_csv(Y, file.path(hmsc, "Y.csv"))
  write_csv(X, file.path(hmsc, "XData.csv"))
  write_csv(Tr, file.path(hmsc, "TrData.csv"))
  write_csv(study, file.path(hmsc, "studyDesign.csv"))
  write_csv(coords, file.path(hmsc, "coordinates.csv"))
  copy_file(file.path(src, "phylogeny.nwk"), file.path(hmsc, "phylogeny.nwk"))
  copy_file(file.path(src, "phylo_cov.csv"), file.path(hmsc, "phylo_cov.csv"))

  hpc <- file.path(case_dir, "02_Hmsc-HPC")
  for (nm in c("Y.csv", "XData.csv", "traits.csv", "studyDesign.csv", "coordinates.csv", "phylo_cov.csv", "newdata.csv")) dir.create(hpc, recursive = TRUE, showWarnings = FALSE)
  write_csv(Y, file.path(hpc, "Y.csv"))
  write_csv(X, file.path(hpc, "XData.csv"))
  write_csv(Tr, file.path(hpc, "traits.csv"))
  write_csv(study, file.path(hpc, "studyDesign.csv"))
  write_csv(coords, file.path(hpc, "coordinates.csv"))
  copy_file(file.path(src, "phylo_cov.csv"), file.path(hpc, "phylo_cov.csv"))
  write_csv(newdata, file.path(hpc, "newdata.csv"))

  jsdm <- file.path(case_dir, "03_jSDM")
  write_csv(Y, file.path(jsdm, "Y.csv"))
  write_csv(X, file.path(jsdm, "XData.csv"))
  write_csv(Tr, file.path(jsdm, "trait_data.csv"))
  write_csv(make_long_format(Y, X), file.path(jsdm, "long_format.csv"), row.names = FALSE)
  write_csv(trial, file.path(jsdm, "trials.csv"))
  write_csv(newdata, file.path(jsdm, "newdata.csv"))
  write_csv(data.frame(Id_sites = "all", Id_species = "all"), file.path(jsdm, "prediction_ids.csv"), row.names = FALSE)

  gjam <- file.path(case_dir, "04_GJAM")
  write_csv(Y, file.path(gjam, "Y.csv"))
  write_csv(X, file.path(gjam, "XData.csv"))
  write_csv(data.frame(response = species, typeName = "PA"), file.path(gjam, "typeNames.csv"), row.names = FALSE)
  write_csv(Tr, file.path(gjam, "specByTrait.csv"))
  write_csv(data.frame(trait = names(Tr), typeName = c("CAT", rep("CON", ncol(Tr) - 1))), file.path(gjam, "traitTypes.csv"), row.names = FALSE)
  write_csv(data.frame(effort = rep(1, length(sites)), row.names = sites), file.path(gjam, "effort.csv"))
  write_csv(newdata, file.path(gjam, "newdata.csv"))
  write_csv(data.frame(row = integer(), note = character()), file.path(gjam, "holdoutIndex.csv"), row.names = FALSE)
  write_csv(data.frame(response = character(), value = numeric(), lower = numeric(), upper = numeric()), file.path(gjam, "censor.csv"), row.names = FALSE)

  spocc <- file.path(case_dir, "05_spOccupancy")
  copy_file(file.path(src, "spOccupancy_y_detection.csv"), file.path(spocc, "y.csv"))
  copy_file(file.path(src, "occ.covs.csv"), file.path(spocc, "occ.covs.csv"))
  copy_file(file.path(src, "det.covs.csv"), file.path(spocc, "det.covs.csv"))
  write_csv(coords, file.path(spocc, "coords.csv"))
  write_csv(data.frame(species = species), file.path(spocc, "species.csv"), row.names = FALSE)
  write_csv(data.frame(source = "book_chapter_upload_ready", rows = length(sites)), file.path(spocc, "integrated_sources.csv"), row.names = FALSE)
  write_csv(newdata, file.path(spocc, "newdata.csv"))
  write_csv(newcoords, file.path(spocc, "newcoords.csv"))
  copy_file(file.path(src, "folds.csv"), file.path(spocc, "folds.csv"))

  sjsdm <- file.path(case_dir, "06_sjSDM")
  write_csv(Y, file.path(sjsdm, "Y.csv"))
  write_csv(X, file.path(sjsdm, "env.csv"))
  write_csv(coords, file.path(sjsdm, "spatial.csv"))
  write_csv(Tr, file.path(sjsdm, "traits.csv"))
  write_csv(newdata, file.path(sjsdm, "newdata.csv"))

  boral <- file.path(case_dir, "07_boral")
  write_csv(Y, file.path(boral, "Y.csv"))
  write_csv(X, file.path(boral, "XData.csv"))
  write_csv(Tr, file.path(boral, "traits.csv"))
  for (f in c("row.ids.csv", "ranef.ids.csv", "distmat.csv", "offset.csv", "trial.size.csv")) copy_file(file.path(src, f), file.path(boral, f))
  write_csv(newdata, file.path(boral, "newdata.csv"))

  write_engine_guide(case_dir, case_id, unname(cases[[case_id]]))
}

summary <- data.frame(
  case_id = names(cases),
  title = unname(cases),
  upload_ready_folder = normalizePath(file.path(out_root, names(cases)), winslash = "/", mustWork = FALSE),
  stringsAsFactors = FALSE
)
write.csv(summary, file.path(out_root, "UPLOAD_READY_CASE_INDEX.csv"), row.names = FALSE)

writeLines(c(
  "# JSDM Studio Upload-Ready Book Chapter Cases",
  "",
  "This folder is the clean upload entry point. Use this instead of browsing the larger examples directory.",
  "",
  "Choose one case folder:",
  "- ch06_plant_traits_whittaker",
  "- ch07_deadwood_fungi",
  "- ch11_finnish_birds",
  "",
  "Inside each case folder, choose the R package / workflow folder:",
  "- 01_Hmsc",
  "- 02_Hmsc-HPC",
  "- 03_jSDM",
  "- 04_GJAM",
  "- 05_spOccupancy",
  "- 06_sjSDM",
  "- 07_boral",
  "",
  "Then upload the files with the same names to the matching JSDM Studio workflow panel.",
  "Open each case's README_UPLOAD_FIRST.md and UPLOAD_GUIDE.csv for exact slot-by-slot mapping.",
  "",
  "The benchmark runs already completed successfully for all three cases and all seven engines. See output/book_chapter_benchmark_run_audit_20260818.csv."
), file.path(out_root, "README_START_HERE.md"), useBytes = TRUE)

zipfile <- file.path(app_dir, "output", "JSDMStudio_book_chapter_examples_UPLOAD_READY_20260818.zip")
dir.create(dirname(zipfile), recursive = TRUE, showWarnings = FALSE)
if (file.exists(zipfile)) unlink(zipfile)
rel_files <- list.files(out_root, recursive = TRUE, all.files = FALSE, no.. = TRUE)
if (requireNamespace("zip", quietly = TRUE)) {
  zip::zipr(zipfile = zipfile, files = rel_files, root = out_root, recurse = FALSE, mode = "mirror")
} else {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(out_root)
  utils::zip(zipfile, rel_files, flags = "-r9Xq")
}
cat("UPLOAD_READY_FOLDER=", normalizePath(out_root, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("UPLOAD_READY_ZIP=", normalizePath(zipfile, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("UPLOAD_READY_FILES=", length(list.files(out_root, recursive = TRUE)), "\n", sep = "")
