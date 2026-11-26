options(stringsAsFactors = FALSE)

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_arg <- args0[startsWith(args0, file_arg)]
script_path <- if (length(script_arg)) sub(file_arg, "", script_arg[[1]]) else "workflow_scripts/summarize_public_case_dimensions.R"
app_dir <- normalizePath(file.path(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) {
  app_dir <- "C:/Users/Google/Documents/CodeX/JSDMStudio_FINAL_COMPLETE_MAX_20260811_HMSCSTUDIO_INSPIRED_EXPLANATIONS/JSDMStudio"
}
root <- file.path(app_dir, "examples", "00_UPLOAD_READY_public_integrated_cases")
names_path <- file.path(root, "PUBLIC_CASE_NAMES.csv")
if (!file.exists(names_path)) stop("Missing PUBLIC_CASE_NAMES.csv", call. = FALSE)

read_dim <- function(path, row.names = 1) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(row.names) && ncol(x) >= row.names) {
    rn <- x[[row.names]]
    if (!anyDuplicated(rn) && all(nzchar(as.character(rn)))) x <- x[-row.names]
  }
  dim(x)
}

nm <- read.csv(names_path, check.names = FALSE, stringsAsFactors = FALSE)
rows <- lapply(seq_len(nrow(nm)), function(i) {
  id <- nm$case_id[[i]]
  d <- file.path(root, id, "00_common_inputs")
  y <- read_dim(file.path(d, "Y_occurrence.csv"))
  x <- read_dim(file.path(d, "XData.csv"))
  tr <- read_dim(file.path(d, "traits.csv"))
  sd <- read_dim(file.path(d, "studyDesign.csv"))
  co <- read_dim(file.path(d, "coordinates.csv"))
  nd <- read_dim(file.path(d, "newdata.csv"))
  det <- read_dim(file.path(d, "spOccupancy_y_detection.csv"))
  files <- list.files(d, full.names = TRUE)
  size_kb <- round(sum(file.info(files)$size, na.rm = TRUE) / 1024, 1)
  visits <- if (y[[2]] > 0) det[[2]] / y[[2]] else NA_real_
  data.frame(
    case_id = id,
    short_id = nm$short_id[[i]],
    display_name = nm$display_name[[i]],
    sites = y[[1]],
    species = y[[2]],
    X_predictors = x[[2]],
    trait_cols = tr[[2]],
    study_cols = sd[[2]],
    coord_cols = co[[2]],
    detection_columns = det[[2]],
    inferred_visits = visits,
    newdata_rows = nd[[1]],
    common_files = length(files),
    common_size_kb = size_kb,
    sites_gt_30 = y[[1]] > 30,
    stringsAsFactors = FALSE
  )
})

out <- do.call(rbind, rows)
write.csv(out, file.path(root, "PUBLIC_CASE_DIMENSIONS.csv"), row.names = FALSE)
print(out[, c("short_id", "sites", "species", "X_predictors", "trait_cols",
              "detection_columns", "inferred_visits", "newdata_rows",
              "common_size_kb", "sites_gt_30")], row.names = FALSE)
if (any(!out$sites_gt_30)) stop("At least one case has <=30 sites.", call. = FALSE)
