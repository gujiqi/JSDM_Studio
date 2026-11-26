
# ============================================================
# JSDM Studio: engine workflows and local workbench
# ============================================================

options(shiny.maxRequestSize = 1024^3)
studio_version <- "0.4.8"

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(yaml)
  library(htmltools)
})

fileInput <- function(...) {
  # Shiny's distant negative offsets create horizontal scrolling in RTL layouts.
  htmltools::tagQuery(shiny::fileInput(...))$find(".shiny-input-file")$
    removeAttrs("style")$addClass("studio-file-native")$allTags()
}

safe_require <- function(pkg) {
  suppressWarnings(suppressMessages(requireNamespace(pkg, quietly = TRUE)))
}
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (length(x) == 1 && is.atomic(x) && is.na(x)) return(y)
  x
}

neutralize_windows_r_locale_env <- function() {
  if (.Platform$OS.type == "windows") {
    Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
    # Reticulate reads Python's redirected discovery output in UTF-8 on Windows R.
    Sys.setenv(PYTHONIOENCODING = "utf-8")
  }
  invisible(TRUE)
}
neutralize_windows_r_locale_env()

app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
dir.create(file.path(app_dir, "output"), showWarnings = FALSE, recursive = TRUE)
source(file.path(app_dir, "R", "parameter_reference.R"), local = TRUE)
source(file.path(app_dir, "R", "project_recommendations.R"), local = TRUE)
source(file.path(app_dir, "R", "studio_preferences.R"), local = TRUE, encoding = "UTF-8")

# ------------------------------------------------------------
# Safe utility functions
# ------------------------------------------------------------

timestamp_id <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

read_studio_csv <- function(path, id_columns = c("", "...1", "row.names", "rowname", "row_id", "id",
                                                "site", "site_id", "sample", "sample_id", "taxon", "taxon_id",
                                                "species", "species_id"), text_columns = character()) {
  dat <- read.csv(path, check.names = FALSE, colClasses = "character",
                  row.names = NULL, fileEncoding = "UTF-8-BOM")
  ids <- NULL
  id_column <- if (ncol(dat) > 1L && names(dat)[1] %in% id_columns) names(dat)[1] else NULL
  if (!is.null(id_column)) {
    ids <- dat[[1L]]
    if (anyNA(ids) || any(!nzchar(ids)) || any(ids != trimws(ids)) || anyDuplicated(ids)) {
      stop("The first ID column must contain unique, nonempty IDs without surrounding spaces.", call. = FALSE)
    }
    dat <- dat[-1L]
  }
  if (!ncol(dat) || any(!nzchar(names(dat))) || anyDuplicated(names(dat))) {
    stop("Data columns must have nonempty, unique names.", call. = FALSE)
  }
  values <- setdiff(names(dat), text_columns)
  dat[values] <- lapply(dat[values], utils::type.convert, as.is = TRUE)
  if (!is.null(ids)) {
    rownames(dat) <- ids
    attr(dat, "studio_explicit_ids") <- TRUE
  }
  dat
}

read_csv_safe <- function(path, ...) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  tryCatch(read_studio_csv(path, ...), error = function(e) NULL)
}

studio_labelled_rows <- function(dat) {
  if (is.null(dat)) return(FALSE)
  ids <- rownames(dat)
  isTRUE(attr(dat, "studio_explicit_ids")) ||
    (!is.null(ids) && !identical(ids, as.character(seq_len(nrow(dat)))))
}

studio_row_order_issue <- function(reference, dat, label) {
  if (!studio_labelled_rows(reference) || !studio_labelled_rows(dat) ||
      nrow(reference) != nrow(dat) || identical(rownames(reference), rownames(dat))) return(character())
  if (setequal(rownames(reference), rownames(dat))) {
    paste(label, "has the right IDs in a different order. Reorder it to match the response rows.")
  } else {
    paste(label, "row IDs do not match the response row IDs.")
  }
}

studio_trait_order_issue <- function(Y, traits) {
  if (is.null(Y) || !studio_labelled_rows(traits) || nrow(traits) != ncol(Y) ||
      identical(rownames(traits), colnames(Y))) return(character())
  if (setequal(rownames(traits), colnames(Y))) {
    "Reorder the trait rows to match the species columns of Y."
  } else {
    "Trait row IDs do not match the species columns of Y."
  }
}

write_lines <- function(x, path) {
  writeLines(as.character(x), con = path, useBytes = TRUE)
}

engine_output_dirs <- function(engine = NULL) {
  core <- c(
    "", "inputs", "data", "models", "results", "tables", "plots",
    "diagnostics", "report", "predictions", "reproducible_script",
    "standard", "workflow_scripts"
  )
  key <- tolower(gsub("[^a-z0-9]+", "", engine %||% ""))
  extra <- character()
  if (grepl("hmschpc", key)) {
    extra <- c("samples", "compile", "hdf5", "python_logs")
  } else if (grepl("jsdm", key) && !grepl("sjsdm", key)) {
    extra <- c("mcmc")
  } else if (grepl("gjam", key)) {
    extra <- c("chains", "sensitivity", "ordination", "missing_data")
  } else if (grepl("spoccupancy|spocc", key)) {
    extra <- c("samples", "spatial", "model_assessment")
  } else if (grepl("sjsdm", key)) {
    extra <- c("anova", "internal_structure", "importance", "weights", "residuals")
  } else if (grepl("boral", key)) {
    extra <- c("jags", "mcmc", "ordination", "residuals", "random_effects", "variable_selection")
  } else if (grepl("hmsc", key)) {
    extra <- c("samples")
  }
  unique(c(core, extra))
}

make_engine_run_dir <- function(engine, project_name = "JSDMWorkbench") {
  clean_name <- gsub("[^A-Za-z0-9_-]+", "_", project_name %||% "JSDMWorkbench")
  engine_clean <- gsub("[^A-Za-z0-9_-]+", "_", engine)
  outdir <- tempfile(paste0(clean_name, "_", engine_clean, "_", timestamp_id(), "_"),
                     tmpdir = file.path(app_dir, "output"))
  if (!dir.create(outdir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create run output folder: ", outdir, call. = FALSE)
  }
  for (d in engine_output_dirs(engine)) {
    dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
  }
  outdir
}

copy_upload <- function(fileinfo, outdir, target_name) {
  if (is.null(fileinfo) || is.null(fileinfo$datapath) || !file.exists(fileinfo$datapath)) return(NULL)
  dest_inputs <- file.path(outdir, "inputs", target_name)
  dest_data <- file.path(outdir, "data", target_name)
  dir.create(dirname(dest_inputs), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(dest_data), recursive = TRUE, showWarnings = FALSE)
  file.copy(fileinfo$datapath, dest_inputs, overwrite = TRUE)
  file.copy(fileinfo$datapath, dest_data, overwrite = TRUE)
  dest_inputs
}

make_zip <- function(outdir) {
  if (is.null(outdir) || !dir.exists(outdir)) stop("Cannot create ZIP because output folder does not exist.", call. = FALSE)
  if (exists("ensure_output_contract", mode = "function")) {
    ensure_output_contract(outdir)
  }
  zipfile <- paste0(outdir, ".zip")
  if (file.exists(zipfile)) unlink(zipfile)
  manifest <- file.path(outdir, "standard", "output_manifest.csv")
  dir.create(dirname(manifest), recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(file = list.files(outdir, recursive = TRUE), stringsAsFactors = FALSE),
            manifest, row.names = FALSE)
  rel_files <- list.files(outdir, recursive = TRUE, all.files = FALSE, no.. = TRUE)
  if (!length(rel_files)) stop("Cannot create ZIP because output folder has no files.", call. = FALSE)
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(
      zipfile = zipfile,
      files = rel_files,
      recurse = FALSE,
      root = outdir,
      mode = "mirror"
    )
  } else {
    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(outdir)
    utils::zip(zipfile, files = rel_files, flags = "-r9Xq")
  }
  if (!file.exists(zipfile) || file.info(zipfile)$size <= 0) {
    stop("ZIP creation failed or produced an empty file.", call. = FALSE)
  }
  zipfile
}

write_session_info <- function(outdir) {
  dir.create(file.path(outdir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    writeLines(capture.output(utils::sessionInfo()), file.path(outdir, "diagnostics", "session_info.txt"))
  }, error = function(e) {
    writeLines(paste("sessionInfo() failed:", conditionMessage(e)),
               file.path(outdir, "diagnostics", "session_info.txt"))
  })
}

write_data_check_messages <- function(outdir, messages = character()) {
  dir.create(file.path(outdir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
  if (is.null(messages) || length(messages) == 0) messages <- "No data-check messages were recorded."
  write.csv(data.frame(message = as.character(messages), stringsAsFactors = FALSE),
            file.path(outdir, "diagnostics", "data_check_messages.csv"), row.names = FALSE)
}

contract_status_values <- function() c("fitted", "model_defined", "check_failed", "fit_failed")

normalize_engine_status <- function(status) {
  if (is.null(status)) status <- list(status = "fit_failed")
  is_list <- is.list(status)
  raw <- if (is_list) status$status else status
  raw <- as.character(raw %||% "fit_failed")[1]
  if (!nzchar(raw) || is.na(raw)) raw <- "fit_failed"
  key <- tolower(trimws(raw))
  msg <- if (is_list) paste(c(status$warnings %||% character(), status$errors %||% character()), collapse = "; ") else ""
  dependency_problem <- grepl("not available|not installed|missing|required package|dependency|could not be loaded", msg, ignore.case = TRUE)
  mapped <- if (key %in% contract_status_values()) {
    key
  } else if (key %in% c("not_run", "waiting", "not_started", "missing", "missing_status", "not_ready_or_unknown", "not_ready", "stopped")) {
    "check_failed"
  } else if (key %in% c("scaffold_only", "ready")) {
    if (dependency_problem) "check_failed" else "model_defined"
  } else if (key %in% c("ready_for_real_fit", "model_ready", "compiled", "compile_only")) {
    "model_defined"
  } else if (grepl("scaffold|defined|compile|ready", key)) {
    if (dependency_problem) "check_failed" else "model_defined"
  } else if (grepl("fail|error|unknown", key)) {
    "fit_failed"
  } else if (grepl("complete|success|done", key)) {
    "model_defined"
  } else {
    "fit_failed"
  }
  if (!is_list) {
    status <- list(status = raw, warnings = character(), errors = character())
  }
  if (!identical(mapped, raw) && is.null(status$raw_status)) status$raw_status <- raw
  status$status <- mapped
  if (is.null(status$warnings) || length(status$warnings) == 0) status$warnings <- character()
  if (is.null(status$errors) || length(status$errors) == 0) status$errors <- character()
  if (identical(mapped, "model_defined") && length(status$warnings) == 0) {
    status$warnings <- "Model boundary, scaffold or executable script was defined; no fitted posterior or performance evidence is claimed."
  }
  if (identical(mapped, "check_failed") && length(status$errors) == 0) {
    status$errors <- "Workflow did not reach fitting. Inspect diagnostics/data_check_messages.csv and used_config.yml."
  }
  if (identical(mapped, "fit_failed") && length(status$errors) == 0) {
    status$errors <- "Fitting or post-processing failed, or the workflow ended without a recognized success status."
  }
  status
}

write_engine_status <- function(outdir, status) {
  status <- normalize_engine_status(status)
  dir.create(file.path(outdir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  if (is.null(status$warnings) || length(status$warnings) == 0) status$warnings <- character()
  if (is.null(status$errors) || length(status$errors) == 0) status$errors <- character()
  default_caveats <- c(
      "A fitted software-test run means the workflow and output contract executed; publication inference still requires adequate sampler settings and convergence diagnostics.",
      "Association parameters are engine-specific and should be compared only where standard tables explicitly mark them as comparable."
  )
  status$scientific_caveats <- unique(c(
    status$scientific_caveats %||% character(),
    status$risk_notes %||% character(),
    default_caveats
  ))
  status$risk_notes <- NULL
  operational_risks <- unique(as.character(status$operational_risks %||% character()))
  status$operational_risks <- operational_risks
  if (length(operational_risks) == 0) {
    status$operational_risks <- character()
  }
  writeLines(jsonlite::toJSON(status, pretty = TRUE, auto_unbox = TRUE),
             file.path(outdir, "diagnostics", "engine_status.json"))
  write.csv(data.frame(
    engine = status$engine %||% "unknown",
    status = status$status %||% "unknown",
    message = paste(c(status$warnings, status$errors), collapse = "; "),
    stringsAsFactors = FALSE
  ), file.path(outdir, "tables", "engine_status.csv"), row.names = FALSE)
  if (length(status$errors) > 0) {
    err_name <- paste0(gsub("[^A-Za-z0-9_-]+", "_", status$engine %||% "engine"), "_error.txt")
    writeLines(as.character(status$errors), file.path(outdir, "diagnostics", err_name))
  }
  write.csv(data.frame(
    engine = status$engine %||% "unknown",
    caveat = c("software_test_vs_publication_inference", "engine_specific_association_scales",
               rep("engine_or_case_specific_caveat", max(0L, length(status$scientific_caveats) - 2L))),
    control = c("Increase sampler settings and inspect diagnostics before ecological inference.",
                "Use comparable flags and notes in standard association tables.",
                rep("Read engine-specific standard tables, reports and diagnostics before interpretation.",
                    max(0L, length(status$scientific_caveats) - 2L))),
    note = as.character(status$scientific_caveats),
    stringsAsFactors = FALSE
  ), file.path(outdir, "diagnostics", "scientific_caveats.csv"), row.names = FALSE)
  if (length(operational_risks) > 0) {
    write.csv(data.frame(
      engine = status$engine %||% "unknown",
      risk = seq_along(operational_risks),
      control = "Treat this as an operational risk and inspect diagnostics before trusting the output.",
      note = operational_risks,
      stringsAsFactors = FALSE
    ), file.path(outdir, "diagnostics", "risk_flags.csv"), row.names = FALSE)
  }
  write_session_info(outdir)
  invisible(status)
}

write_standard_outputs <- function(outdir, engine, status, Y = NULL, X = NULL) {
  status <- normalize_engine_status(status)
  status_value <- status$status %||% "fit_failed"
  dir.create(file.path(outdir, "standard"), recursive = TRUE, showWarnings = FALSE)
  species_names <- if (!is.null(Y) && !is.null(colnames(Y))) colnames(Y) else NA_character_
  site_names <- if (!is.null(Y) && !is.null(rownames(Y))) rownames(Y) else NA_character_
  predictor_names <- if (!is.null(X) && !is.null(colnames(X))) colnames(X) else NA_character_
  write.csv(data.frame(
    run_id = basename(outdir),
    engine = engine,
    status = status_value,
    n_sites = if (is.null(Y)) NA_integer_ else nrow(Y),
    n_responses = if (is.null(Y)) NA_integer_ else ncol(Y),
    n_predictors = if (is.null(X)) NA_integer_ else ncol(X),
    stringsAsFactors = FALSE
  ), file.path(outdir, "standard", "run_summary.csv"), row.names = FALSE)
  standard_tables <- list(
    effects_long = data.frame(engine=engine, response_id=species_names[1], predictor=predictor_names[1], direction=NA_character_, estimate=NA_real_, lower=NA_real_, upper=NA_real_, notes=status_value, stringsAsFactors=FALSE),
    predictions_long = data.frame(engine=engine, site_id=NA_character_, response_id=NA_character_, observed=NA_real_, predicted_mean=NA_real_, predicted_lower=NA_real_, predicted_upper=NA_real_, stringsAsFactors=FALSE),
    associations_long = data.frame(engine=engine, response_1=NA_character_, response_2=NA_character_, association_type=NA_character_, estimate=NA_real_, comparable_level=status_value, stringsAsFactors=FALSE),
    fit_metrics = data.frame(engine=engine, metric=NA_character_, response_id=NA_character_, value=NA_real_, notes=status_value, stringsAsFactors=FALSE),
    diagnostics_long = data.frame(engine=engine, diagnostic="engine_status", status=status_value, value=NA_real_, note=paste(c(status$warnings, status$errors), collapse="; "), stringsAsFactors=FALSE),
    effects_species_environment = data.frame(engine=engine, species=species_names[1], predictor=predictor_names[1], estimate=NA_real_, lower=NA_real_, upper=NA_real_, statistic=NA_real_, p_or_support=NA_real_, effect_type="environment", scale="engine_specific", comparable=FALSE, note=status_value, stringsAsFactors=FALSE),
    predictions_site_species = data.frame(engine=engine, site_id=site_names[1], species=species_names[1], observed=NA_real_, predicted=NA_real_, truth_probability=NA_real_, residual=NA_real_, prediction_scale="engine_specific", stringsAsFactors=FALSE),
    associations_species_species = data.frame(engine=engine, species_i=NA_character_, species_j=NA_character_, estimate=NA_real_, lower=NA_real_, upper=NA_real_, association_type="engine_specific", scale="engine_specific", comparable=FALSE, note=status_value, stringsAsFactors=FALSE)
  )
  for (nm in names(standard_tables)) {
    write.csv(standard_tables[[nm]], file.path(outdir, "standard", paste0(nm, ".csv")), row.names = FALSE)
  }
  write.csv(data.frame(file = list.files(outdir, recursive = TRUE), stringsAsFactors = FALSE),
    file.path(outdir, "standard", "output_manifest.csv"), row.names = FALSE)
}

read_engine_status_safe <- function(outdir, engine = NULL) {
  status_file <- file.path(outdir, "diagnostics", "engine_status.json")
  if (file.exists(status_file) && requireNamespace("jsonlite", quietly = TRUE)) {
    st <- tryCatch(jsonlite::fromJSON(status_file, simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(st)) return(normalize_engine_status(st))
  }
  run_file <- file.path(outdir, "standard", "run_summary.csv")
  if (file.exists(run_file)) {
    run <- tryCatch(read.csv(run_file, check.names = FALSE), error = function(e) NULL)
    if (!is.null(run) && nrow(run) > 0) {
      return(normalize_engine_status(list(
        engine = engine %||% run$engine[1] %||% "unknown",
        status = run$status[1] %||% "fit_failed",
        warnings = character(),
        errors = character()
      )))
    }
  }
  normalize_engine_status(list(engine = engine %||% basename(outdir), status = "check_failed",
                               errors = "No engine_status.json or standard/run_summary.csv was available."))
}

standard_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
}

ensure_output_contract <- function(outdir, engine = NULL, status = NULL, Y = NULL, X = NULL) {
  if (is.null(outdir) || !dir.exists(outdir)) return(invisible(FALSE))
  status <- normalize_engine_status(status %||% read_engine_status_safe(outdir, engine))
  engine <- engine %||% status$engine %||% basename(outdir)
  status$engine <- engine
  for (d in engine_output_dirs(engine)) dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(file.path(outdir, "used_config.yml"))) {
    if (requireNamespace("yaml", quietly = TRUE)) {
      try(yaml::write_yaml(list(engine = engine, status = status$status, note = "Minimal config written by output-contract repair."), file.path(outdir, "used_config.yml")), silent = TRUE)
    } else {
      writeLines(c(paste0("engine: ", engine), paste0("status: ", status$status)), file.path(outdir, "used_config.yml"))
    }
  }
  if (!file.exists(file.path(outdir, "diagnostics", "data_check_messages.csv"))) {
    write_data_check_messages(outdir, "No data-check messages were recorded before output-contract repair.")
  }
  status <- write_engine_status(outdir, status)
  std <- file.path(outdir, "standard")
  dir.create(std, recursive = TRUE, showWarnings = FALSE)
  placeholder <- function(nm, dat) {
    p <- file.path(std, paste0(nm, ".csv"))
    if (!file.exists(p)) write.csv(dat, p, row.names = FALSE)
  }
  placeholder("run_summary", data.frame(run_id = basename(outdir), engine = engine, status = status$status,
                                        n_sites = if (is.null(Y)) NA_integer_ else nrow(Y),
                                        n_responses = if (is.null(Y)) NA_integer_ else ncol(Y),
                                        n_predictors = if (is.null(X)) NA_integer_ else ncol(X),
                                        stringsAsFactors = FALSE))
  run_summary <- standard_read_csv(file.path(std, "run_summary.csv"))
  if (!is.null(run_summary) && nrow(run_summary) > 0) {
    if (!"engine" %in% names(run_summary)) run_summary$engine <- engine
    if (!"status" %in% names(run_summary)) run_summary$status <- status$status
    run_summary$engine[1] <- engine
    run_summary$status[1] <- status$status
    write.csv(run_summary, file.path(std, "run_summary.csv"), row.names = FALSE)
  }
  placeholder("effects_long", data.frame(engine=engine, response_id=NA_character_, predictor=NA_character_, direction=NA_character_, estimate=NA_real_, lower=NA_real_, upper=NA_real_, notes=status$status, stringsAsFactors=FALSE))
  placeholder("predictions_long", data.frame(engine=engine, site_id=NA_character_, response_id=NA_character_, observed=NA_real_, predicted_mean=NA_real_, predicted_lower=NA_real_, predicted_upper=NA_real_, stringsAsFactors=FALSE))
  placeholder("associations_long", data.frame(engine=engine, response_1=NA_character_, response_2=NA_character_, association_type=NA_character_, estimate=NA_real_, comparable_level=status$status, stringsAsFactors=FALSE))
  placeholder("fit_metrics", data.frame(engine=engine, metric=NA_character_, response_id=NA_character_, value=NA_real_, notes=status$status, stringsAsFactors=FALSE))
  diag_note <- paste(c(status$warnings %||% character(), status$errors %||% character()), collapse = "; ")
  write.csv(data.frame(engine=engine, diagnostic="engine_status", status=status$status,
                       value=NA_real_, note=diag_note, stringsAsFactors=FALSE),
            file.path(std, "diagnostics_long.csv"), row.names = FALSE)
  eff <- standard_read_csv(file.path(std, "effects_long.csv"))
  if (!is.null(eff) && nrow(eff) > 0) {
    species <- eff$response_id %||% eff$species %||% NA_character_
    pred <- eff$predictor %||% NA_character_
    write.csv(data.frame(engine=engine, species=species, predictor=pred,
                         estimate=suppressWarnings(as.numeric(eff$estimate %||% NA_real_)),
                         lower=suppressWarnings(as.numeric(eff$lower %||% NA_real_)),
                         upper=suppressWarnings(as.numeric(eff$upper %||% NA_real_)),
                         statistic=NA_real_, p_or_support=NA_real_,
                         effect_type=eff$effect_type %||% if (identical(engine, "spOccupancy")) "occurrence_or_detection_check_effect_type" else "environment",
                         scale=eff$scale %||% "engine_specific",
                         comparable=FALSE,
                         note=eff$notes %||% status$status,
                         stringsAsFactors=FALSE),
              file.path(std, "effects_species_environment.csv"), row.names = FALSE)
  }
  pred <- standard_read_csv(file.path(std, "predictions_long.csv"))
  if (!is.null(pred) && nrow(pred) > 0) {
    predicted <- pred$predicted_mean %||% pred$predicted %||% NA_real_
    observed <- pred$observed %||% NA_real_
    write.csv(data.frame(engine=engine, site_id=pred$site_id %||% NA_character_,
                         species=pred$response_id %||% pred$species %||% NA_character_,
                         observed=suppressWarnings(as.numeric(observed)),
                         predicted=suppressWarnings(as.numeric(predicted)),
                         truth_probability=NA_real_,
                         residual=if (identical(engine, "spOccupancy")) NA_real_ else suppressWarnings(as.numeric(observed)) - suppressWarnings(as.numeric(predicted)),
                         prediction_scale=pred$prediction_scale %||% if (identical(engine, "spOccupancy")) "occupancy_probability" else "engine_specific",
                         stringsAsFactors=FALSE),
              file.path(std, "predictions_site_species.csv"), row.names = FALSE)
  }
  assoc <- standard_read_csv(file.path(std, "associations_long.csv"))
  if (!is.null(assoc) && nrow(assoc) > 0) {
    write.csv(data.frame(engine=engine,
                         species_i=assoc$response_1 %||% assoc$species_i %||% NA_character_,
                         species_j=assoc$response_2 %||% assoc$species_j %||% NA_character_,
                         estimate=suppressWarnings(as.numeric(assoc$estimate %||% NA_real_)),
                         lower=NA_real_, upper=NA_real_,
                         association_type=assoc$association_type %||% "engine_specific",
                         scale="engine_specific",
                         comparable=FALSE,
                         note=paste("Associations are engine-specific; raw numeric values are not interchangeable across engines. Status:", status$status),
                         stringsAsFactors=FALSE),
              file.path(std, "associations_species_species.csv"), row.names = FALSE)
  }
  write.csv(data.frame(file = list.files(outdir, recursive = TRUE), stringsAsFactors = FALSE),
            file.path(std, "output_manifest.csv"), row.names = FALSE)
  write_folder_readmes(outdir, engine)
  invisible(TRUE)
}

write_reproducible_stub <- function(outdir, engine) {
  dir.create(file.path(outdir, "reproducible_script"), recursive = TRUE, showWarnings = FALSE)
  engine_id <- gsub("[^A-Za-z0-9_]+", "_", engine)
  writeLines(c(
    paste0("# Reproducible ", engine, " analysis/scaffold script generated by JSDMWorkbench"),
    "# This script is executable. In model_defined or check_failed states it rebuilds exported manifest tables,",
    "# diagnostics and reproducibility index from used_config.yml and data/.",
    "args <- commandArgs(trailingOnly = FALSE)",
    "file_arg <- '--file='",
    "script_arg <- args[startsWith(args, file_arg)]",
    "if (length(script_arg) > 0) {",
    "  script_path <- sub(file_arg, '', script_arg[[1]])",
    "  if (!file.exists(script_path)) script_path <- gsub('~+~', ' ', script_path, fixed = TRUE)",
    "  setwd(dirname(dirname(normalizePath(script_path, winslash = '/', mustWork = TRUE))))",
    "}",
    "dir.create('tables', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('diagnostics', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('standard', showWarnings = FALSE, recursive = TRUE)",
    "`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x",
    "cfg <- if (requireNamespace('yaml', quietly = TRUE) && file.exists('used_config.yml')) yaml::read_yaml('used_config.yml') else list()",
    "status_file <- file.path('diagnostics', 'engine_status.json')",
    "input_files <- if (dir.exists('data')) list.files('data', full.names = FALSE) else character()",
    "manifest <- data.frame(file = list.files('.', recursive = TRUE), stringsAsFactors = FALSE)",
    "write.csv(manifest, file.path('standard', 'reproduced_output_manifest.csv'), row.names = FALSE)",
    paste0("write.csv(data.frame(engine = cfg$engine %||% '", engine, "', input_file = input_files, stringsAsFactors = FALSE), file.path('tables', 'reproduced_input_manifest.csv'), row.names = FALSE)"),
    "if (file.exists(status_file)) cat(readLines(status_file), sep = '\\n')",
    "message('Reproducibility scaffold refreshed. Production fitting requires the corresponding engine package and adapter.')"
  ), file.path(outdir, "reproducible_script", paste0("run_this_", engine_id, "_analysis.R")))
  writeLines(c(
    paste0("# Workflow wrapper for ", engine),
    paste0("source(file.path('reproducible_script', 'run_this_", engine_id, "_analysis.R'))")
  ), file.path(outdir, "workflow_scripts", paste0("run_", engine_id, "_workflow.R")))
}

write_sjsdm_reproducible_script <- function(outdir) {
  dir.create(file.path(outdir, "reproducible_script"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "workflow_scripts"), recursive = TRUE, showWarnings = FALSE)
  script <- c(
    "# Reproducible sjSDM 1.0.7 analysis script generated by JSDM Studio",
    "# Run from the output folder with: Rscript reproducible_script/run_this_sjSDM_analysis.R",
    "args <- commandArgs(trailingOnly = FALSE)",
    "file_arg <- '--file='",
    "script_arg <- args[startsWith(args, file_arg)]",
    "if (length(script_arg) > 0) {",
    "  script_path <- sub(file_arg, '', script_arg[[1]])",
    "  if (!file.exists(script_path)) script_path <- gsub('~+~', ' ', script_path, fixed = TRUE)",
    "  setwd(dirname(dirname(normalizePath(script_path, winslash = '/', mustWork = TRUE))))",
    "}",
    "dir.create('models', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('tables', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('predictions', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('anova', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('internal_structure', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('importance', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('weights', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('plots', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('diagnostics', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('standard', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('results', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('report', showWarnings = FALSE, recursive = TRUE)",
    "dir.create('residuals', showWarnings = FALSE, recursive = TRUE)",
    "`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (is.atomic(x) && length(x) == 1L && is.na(x))) y else x",
    "configure_sjsdm_python <- function() {",
    "  configured <- Sys.getenv('RETICULATE_PYTHON', unset = '')",
    "  if (nzchar(configured) && file.exists(configured)) return(normalizePath(configured, winslash = '/', mustWork = FALSE))",
    "  local_appdata <- Sys.getenv('LOCALAPPDATA', unset = '')",
    "  shared_python <- file.path(local_appdata, 'JSDMStudio', 'hmschpc-python310', 'Scripts', 'python.exe')",
    "  if (.Platform$OS.type != 'windows' || !nzchar(local_appdata) || !file.exists(shared_python)) {",
    "    stop('sjSDM needs a Python environment with PyTorch, torch_optimizer, pyro and madgrad. Use the supplied platform launcher or set RETICULATE_PYTHON to that environment.', call. = FALSE)",
    "  }",
    "  Sys.setenv(RETICULATE_PYTHON = normalizePath(shared_python, winslash = '/', mustWork = FALSE))",
    "  Sys.getenv('RETICULATE_PYTHON')",
    "}",
    "as_bool <- function(x, default = FALSE) { if (is.null(x) || length(x) == 0 || is.na(x)) return(default); if (is.logical(x)) return(isTRUE(x)); tolower(as.character(x)) %in% c('true','t','1','yes','y') }",
    "as_num <- function(x, default) { z <- suppressWarnings(as.numeric(x)); if (!length(z) || !is.finite(z[1])) default else z[1] }",
    "as_int <- function(x, default) as.integer(round(as_num(x, default)))",
    "read_csv_safe <- function(path) {",
    "  if (!file.exists(path)) return(NULL)",
    "  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)",
    "  if (ncol(dat) > 1) {",
    "    first <- dat[[1]]",
    "    first_name <- names(dat)[1] %||% ''",
    "    first_chr <- as.character(first)",
    "    first_num <- suppressWarnings(as.numeric(first_chr))",
    "    sequence_index <- all(!is.na(first_num)) && identical(as.integer(first_num), seq_len(length(first_num)))",
    "    row_id_name <- first_name %in% c('', 'X', 'X.1', '...1', 'row.names', 'rowname', 'row_id', 'id', 'site', 'site_id', 'sample', 'sample_id', 'taxon', 'taxon_id', 'species', 'species_id')",
    "    row_id_text <- !all(!is.na(first_num))",
    "    if (!anyDuplicated(first_chr) && (row_id_name || sequence_index || row_id_text)) { dat <- dat[-1]; rownames(dat) <- make.unique(first_chr) }",
    "  }",
    "  dat",
    "}",
    "clean_df <- function(df) {",
    "  if (is.null(df)) return(NULL)",
    "  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)",
    "  for (nm in names(df)) {",
    "    if (is.character(df[[nm]])) {",
    "      x <- trimws(df[[nm]]); x[x == ''] <- NA",
    "      nx <- suppressWarnings(as.numeric(x))",
    "      df[[nm]] <- if (all(is.na(x) | !is.na(nx))) nx else factor(x)",
    "    }",
    "  }",
    "  df",
    "}",
    "clean_new_df <- function(df, template = NULL) {",
    "  if (is.null(df)) return(NULL)",
    "  out <- clean_df(df)",
    "  adj <- list()",
    "  add_adj <- function(column, row_id, original, replacement, reason) adj[[length(adj)+1]] <<- data.frame(column=column, row_id=row_id, original=as.character(original), replacement=as.character(replacement), reason=reason, stringsAsFactors=FALSE)",
    "  if (!is.null(template)) {",
    "    for (nm in intersect(names(out), names(template))) {",
    "      if (is.factor(template[[nm]])) {",
    "        raw <- trimws(as.character(out[[nm]])); raw[raw == ''] <- NA",
    "        lv <- levels(template[[nm]]); repl <- lv[1]",
    "        bad <- is.na(raw) | !(raw %in% lv)",
    "        if (any(bad)) { for (rid in rownames(out)[bad]) add_adj(nm, rid, out[[nm]][rownames(out)==rid][1], repl, 'unseen_or_missing_factor_level_replaced_by_training_reference'); raw[bad] <- repl }",
    "        out[[nm]] <- factor(raw, levels = lv)",
    "      } else if (is.numeric(template[[nm]])) {",
    "        vals <- suppressWarnings(as.numeric(out[[nm]])); bad <- !is.finite(vals)",
    "        if (any(bad)) { repl <- suppressWarnings(stats::median(template[[nm]], na.rm=TRUE)); if (!is.finite(repl)) repl <- 0; for (rid in rownames(out)[bad]) add_adj(nm, rid, out[[nm]][rownames(out)==rid][1], repl, 'non_numeric_or_missing_value_imputed_from_training_median'); vals[bad] <- repl }",
    "        out[[nm]] <- vals",
    "      }",
    "    }",
    "  }",
    "  dir.create('predictions', showWarnings=FALSE, recursive=TRUE)",
    "  if (length(adj)) write.csv(do.call(rbind, adj), file.path('predictions','sjSDM_newdata_adjustments.csv'), row.names=FALSE) else write.csv(data.frame(column=NA_character_, row_id=NA_character_, original=NA_character_, replacement=NA_character_, reason='no_newdata_adjustments_needed'), file.path('predictions','sjSDM_newdata_adjustments.csv'), row.names=FALSE)",
    "  out",
    "}",
    "num_matrix <- function(x) {",
    "  x <- as.matrix(x)",
    "  storage.mode(x) <- 'numeric'",
    "  x",
    "}",
    "family_from_cfg <- function(x) {",
    "  x <- tolower(as.character(x %||% 'binomial_probit'))",
    "  if (x == 'binomial_probit') return(stats::binomial('probit'))",
    "  if (x == 'binomial_logit') return(stats::binomial('logit'))",
    "  if (x == 'poisson_log') return(stats::poisson('log'))",
    "  if (x %in% c('nbinom','negative_binomial_log','negative.binomial')) return('nbinom')",
    "  if (x == 'gaussian_identity') return(stats::gaussian('identity'))",
    "  stop('Unsupported sjSDM family in used_config.yml: ', x, call. = FALSE)",
    "}",
    "hidden_vec <- function(x) {",
    "  z <- suppressWarnings(as.integer(trimws(unlist(strsplit(as.character(x %||% '10,10,10'), ',')))))",
    "  if (!length(z) || anyNA(z) || any(z <= 0)) stop('DNN hidden layers must be comma-separated positive integers.', call. = FALSE)",
    "  if (length(z) == 1) z <- rep(z, 2)",
    "  z",
    "}",
    "optimizer_from_cfg <- function(name, weight_decay) {",
    "  name <- as.character(name %||% 'Adamax')",
    "  wd <- as_num(weight_decay, 0.002)",
    "  switch(name,",
    "    Adamax = sjSDM::Adamax(weight_decay = wd),",
    "    RMSprop = sjSDM::RMSprop(weight_decay = wd),",
    "    SGD = sjSDM::SGD(weight_decay = wd),",
    "    AccSGD = sjSDM::AccSGD(weight_decay = wd),",
    "    AdaBound = sjSDM::AdaBound(weight_decay = wd),",
    "    madgrad = sjSDM::madgrad(weight_decay = wd),",
    "    stop('Unsupported optimizer for sjSDM 1.0.7: ', name, call. = FALSE)",
    "  )",
    "}",
    "write_status <- function(status, errors = character(), warnings = character()) {",
    "  payload <- list(engine='sjSDM', status=status, errors=as.character(errors), warnings=as.character(warnings), time=as.character(Sys.time()))",
    "  if (requireNamespace('jsonlite', quietly = TRUE)) writeLines(jsonlite::toJSON(payload, pretty=TRUE, auto_unbox=TRUE), file.path('diagnostics','engine_status.json'))",
    "  write.csv(data.frame(engine='sjSDM', status=status, message=paste(c(warnings, errors), collapse='; ')), file.path('tables','engine_status.csv'), row.names=FALSE)",
    "}",
    "step_warnings <- character()",
    "safe_step <- function(name, expr) {",
    "  tryCatch(withCallingHandlers(force(expr), warning = function(w) { step_warnings <<- c(step_warnings, paste(name, conditionMessage(w), sep=': ')); invokeRestart('muffleWarning') }), error = function(e) { msg <- conditionMessage(e); step_warnings <<- c(step_warnings, paste(name, msg, sep=': ')); writeLines(msg, file.path('diagnostics', paste0('sjSDM_', gsub('[^A-Za-z0-9_]+', '_', name), '_error.txt'))); NULL })",
    "}",
    "coef_block <- function(cf, model, idx = 1, col_names = NULL) {",
    "  z <- if (is.list(cf)) { if (!is.null(names(cf)) && 'env' %in% names(cf) && idx == 1) cf$env else if (!is.null(names(cf)) && 'spatial' %in% names(cf) && idx == 2) cf$spatial else cf[[idx]] } else cf",
    "  if (is.list(z) && length(z) == 1) z <- z[[1]]",
    "  z <- as.matrix(z)",
    "  if (!is.null(model$species) && nrow(z) == length(model$species)) rownames(z) <- model$species",
    "  if (!is.null(col_names) && ncol(z) == length(col_names)) colnames(z) <- col_names else if (!is.null(model$names) && ncol(z) == length(model$names)) colnames(z) <- model$names",
    "  z",
    "}",
    "name_species_matrix <- function(mat, model) {",
    "  mat <- num_matrix(mat)",
    "  if (!is.null(model$species) && nrow(mat) == length(model$species)) rownames(mat) <- model$species",
    "  if (!is.null(model$species) && ncol(mat) == length(model$species)) colnames(mat) <- model$species",
    "  mat",
    "}",
    "write_effects_long <- function(mat) {",
    "  tab <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)",
    "  names(tab) <- c('response_id','predictor','estimate')",
    "  tab$engine <- 'sjSDM'; tab$direction <- ifelse(tab$estimate > 0, 'positive', ifelse(tab$estimate < 0, 'negative', 'zero')); tab$lower <- NA_real_; tab$upper <- NA_real_; tab$notes <- 'coef.sjSDM environmental/spatial coefficient'",
    "  write.csv(tab[, c('engine','response_id','predictor','direction','estimate','lower','upper','notes')], file.path('standard','effects_long.csv'), row.names = FALSE)",
    "}",
    "write_dnn_coef_artifacts <- function(obj, component) {",
    "  dir.create('weights', showWarnings = FALSE, recursive = TRUE)",
    "  pieces <- if (is.list(obj)) obj else list(weights = obj)",
    "  piece_names <- names(pieces); if (is.null(piece_names) || length(piece_names) != length(pieces) || any(!nzchar(piece_names))) piece_names <- paste0('part_', seq_along(pieces))",
    "  count_values <- function(z) { vals <- tryCatch(unlist(z, recursive = TRUE, use.names = FALSE), error = function(e) NULL); if (is.null(vals)) NA_integer_ else length(vals) }",
    "  saveRDS(obj, file.path('weights', paste0(component, '_dnn_coef_weights.rds')))",
    "  manifest <- data.frame(component = component, part = piece_names, class = vapply(pieces, function(z) paste(class(z), collapse='|'), character(1)), length = vapply(pieces, count_values, integer(1)), stringsAsFactors = FALSE)",
    "  write.csv(manifest, file.path('weights', paste0(component, '_dnn_coef_weights_manifest.csv')), row.names = FALSE)",
    "  write.csv(data.frame(component = component, coefficient_type = 'DNN_weight_object', status = 'saved_as_rds', note = 'DNN weights are not a species-by-predictor coefficient matrix; use weights/*_dnn_coef_weights.rds for reproducibility.', stringsAsFactors = FALSE), file.path('tables', paste0('coef_', component, '.csv')), row.names = FALSE)",
    "  invisible(NULL)",
    "}",
    "write_predictions_long <- function(pred, Yobs = NULL, row_ids = NULL, out_file = file.path('standard','predictions_long.csv')) {",
    "  if (!is.null(row_ids) && length(row_ids) == nrow(pred)) rownames(pred) <- row_ids",
    "  else if (is.null(rownames(pred)) || all(grepl('^[0-9]+$', rownames(pred)))) rownames(pred) <- if (!is.null(Yobs) && nrow(Yobs) == nrow(pred)) rownames(Yobs) %||% paste0('site_', seq_len(nrow(pred))) else paste0('pred_site_', seq_len(nrow(pred)))",
    "  if (is.null(colnames(pred)) || all(grepl('^V[0-9]+$', colnames(pred)))) colnames(pred) <- model$species %||% colnames(Y) %||% paste0('sp_', seq_len(ncol(pred)))",
    "  tab <- as.data.frame(as.table(as.matrix(pred)), stringsAsFactors = FALSE)",
    "  names(tab) <- c('site_id','response_id','predicted_mean')",
    "  obs <- if (!is.null(Yobs) && all(dim(Yobs) == dim(pred))) as.vector(as.matrix(Yobs)) else NA_real_",
    "  tab <- data.frame(engine='sjSDM', site_id=tab$site_id, response_id=tab$response_id, observed=obs, predicted_mean=tab$predicted_mean, predicted_lower=NA_real_, predicted_upper=NA_real_, stringsAsFactors = FALSE)",
    "  write.csv(tab, out_file, row.names = FALSE)",
    "}",
    "write_associations_long <- function(mat, association_type) {",
    "  mat <- as.matrix(mat); ids <- which(upper.tri(mat, diag = FALSE), arr.ind = TRUE)",
    "  if (is.null(rownames(mat))) rownames(mat) <- model$species %||% paste0('sp_', seq_len(nrow(mat)))",
    "  if (is.null(colnames(mat))) colnames(mat) <- model$species %||% paste0('sp_', seq_len(ncol(mat)))",
    "  tab <- data.frame(engine='sjSDM', response_1=rownames(mat)[ids[,1]], response_2=colnames(mat)[ids[,2]], association_type=association_type, estimate=mat[ids], comparable_level='fitted', stringsAsFactors = FALSE)",
    "  write.csv(tab, file.path('standard','associations_long.csv'), row.names = FALSE)",
    "}",
    "nested_to_long <- function(x) {",
    "  rows <- list()",
    "  add <- function(component, item, value) rows[[length(rows) + 1L]] <<- data.frame(component = component, item = item, value = as.character(value), stringsAsFactors = FALSE)",
    "  walk <- function(obj, path = 'root') {",
    "    if (is.null(obj)) return(NULL)",
    "    if (is.function(obj) || is.environment(obj) || typeof(obj) %in% c('externalptr','builtin','special')) return(NULL)",
    "    if (is.data.frame(obj) || is.matrix(obj)) {",
    "      m <- as.data.frame(obj, check.names = FALSE)",
    "      rn <- rownames(m)",
    "      for (i in seq_len(nrow(m))) for (j in names(m)) add(path, paste0(if (!is.null(rn) && nzchar(rn[i])) rn[i] else i, ':', j), m[[j]][i])",
    "      return(NULL)",
    "    }",
    "    if (is.list(obj)) {",
    "      nms <- names(obj); if (is.null(nms)) nms <- as.character(seq_along(obj))",
    "      for (i in seq_along(obj)) walk(obj[[i]], paste(path, nms[[i]], sep = '.'))",
    "      return(NULL)",
    "    }",
    "    if (!is.atomic(obj)) { add(path, 'value', paste(deparse(obj), collapse = ' ')); return(NULL) }",
    "    vals <- as.vector(obj); nms <- names(vals); if (is.null(nms)) nms <- as.character(seq_along(vals))",
    "    for (i in seq_along(vals)) add(path, nms[[i]], vals[[i]])",
    "  }",
    "  walk(x)",
    "  if (length(rows) == 0) data.frame(component = character(), item = character(), value = character(), stringsAsFactors = FALSE) else do.call(rbind, rows)",
    "}",
    "write_runtime_diagnostics <- function() {",
    "  writeLines(capture.output(sessionInfo()), file.path('diagnostics','session_info.txt'))",
    "  backend_lines <- c(paste('sjSDM version:', as.character(utils::packageVersion('sjSDM'))), paste('RETICULATE_PYTHON:', Sys.getenv('RETICULATE_PYTHON', unset = '')))",
    "  backend_lines <- c(backend_lines, paste('reticulate namespace available:', requireNamespace('reticulate', quietly = TRUE)))",
    "  if (requireNamespace('reticulate', quietly = TRUE)) {",
    "    backend_lines <- c(backend_lines, tryCatch({",
    "      cfg <- reticulate::py_config(); reticulate::py_run_string('import torch, torch_optimizer, pyro, madgrad')",
    "      c(paste('Python:', cfg$python), paste('Python version:', cfg$version), paste('PyTorch version:', reticulate::py_eval('torch.__version__')))",
    "    }, error = function(e) paste('PyTorch backend check failed:', conditionMessage(e))))",
    "  }",
    "  writeLines(backend_lines, file.path('diagnostics','torch_diagnostic.txt'))",
    "  if (!file.exists(file.path('diagnostics','data_check_messages.csv'))) write.csv(data.frame(Message = 'Reproducible script executed directly; GUI data-check messages were not supplied.', stringsAsFactors = FALSE), file.path('diagnostics','data_check_messages.csv'), row.names = FALSE)",
    "}",
    "write_input_manifest <- function(Y, env_dat = NULL, spatial_dat = NULL, traits_dat = NULL, groups_dat = NULL, folds_dat = NULL) {",
    "  row <- function(file, obj, role, used) data.frame(file = file, rows = if (is.null(obj)) NA_integer_ else nrow(obj), columns = if (is.null(obj)) NA_integer_ else ncol(obj), role = role, used_in_fit = used, status = if (is.null(obj)) 'not_uploaded' else 'readable', stringsAsFactors = FALSE)",
    "  manifest <- rbind(row('Y.csv', as.data.frame(Y), 'response matrix', TRUE), row('env.csv', env_dat, 'environmental module predictors', TRUE), row('spatial.csv', spatial_dat, 'spatial module predictors/eigenvectors', !is.null(spatial_dat)), row('traits.csv', traits_dat, 'species traits for post-hoc interpretation/assembly plots', FALSE), row('species_groups.csv', groups_dat, 'species grouping metadata for reporting/plots', FALSE), row('folds.csv', folds_dat, 'uploaded CV fold metadata; sjSDM_cv uses integer CV folds', FALSE))",
    "  write.csv(manifest, file.path('tables','sjSDM_input_manifest.csv'), row.names = FALSE)",
    "}",
    "write_run_artifacts <- function(status, warnings = character(), errors = character()) {",
    "  run_tab <- data.frame(engine = 'sjSDM', status = status, n_sites = if (exists('Y')) nrow(Y) else NA_integer_, n_responses = if (exists('Y')) ncol(Y) else NA_integer_, n_env_predictors = if (exists('env_dat')) ncol(env_dat) else NA_integer_, n_spatial_predictors = if (exists('spatial_dat') && !is.null(spatial_dat)) ncol(spatial_dat) else 0L, stringsAsFactors = FALSE)",
    "  write.csv(run_tab, file.path('results','sjSDM_run_summary.csv'), row.names = FALSE)",
    "  writeLines(c('sjSDM results', '=============', paste('Status:', status), '', 'Core outputs:', '- models/sjSDM_model.rds', '- tables/coef_environment.csv and optional coef_spatial.csv', '- tables/covariance_matrix.csv and correlation_matrix.csv', '- predictions/predictions.csv', '- standard/run_summary.csv, effects_long.csv, predictions_long.csv and associations_long.csv', '', 'Diagnostics:', '- diagnostics/engine_status.json', '- diagnostics/session_info.txt', '- diagnostics/torch_diagnostic.txt', '', 'Warnings/errors:', paste(c(warnings, errors), collapse = '; ')), file.path('results','README_sjSDM_results.txt'))",
    "}",
    "write_simple_report <- function(status, warnings = character(), errors = character()) {",
    "  esc <- function(x) gsub('&', '&amp;', gsub('<', '&lt;', gsub('>', '&gt;', as.character(x))))",
    "  files <- list.files('.', recursive = TRUE)",
    "  cfg_text <- if (exists('cfg', inherits = TRUE)) paste(capture.output(str(cfg, max.level = 3)), collapse = '\\n') else 'Configuration was not loaded before the script failed.'",
    "  html <- c('<!doctype html><html><head><meta charset=utf-8><title>sjSDM report</title></head><body>', '<h1>sjSDM analysis report</h1>', paste0('<p><b>Status:</b> ', esc(status), '</p>'), paste0('<p><b>Generated:</b> ', esc(Sys.time()), '</p>'), '<h2>Configuration</h2><pre>', esc(cfg_text), '</pre>', '<h2>Warnings and errors</h2><pre>', esc(paste(c(warnings, errors), collapse = '\\n')), '</pre>', '<h2>Output files</h2><pre>', esc(paste(files, collapse = '\\n')), '</pre>', '</body></html>')",
    "  writeLines(html, file.path('report','sjSDM_report.html'))",
    "}",
    "fill_empty_output_dirs <- function() {",
    "  dirs <- c('inputs','data','models','results','tables','plots','diagnostics','report','predictions','reproducible_script','standard','workflow_scripts','anova','internal_structure','importance','weights','residuals')",
    "  for (d in dirs) if (dir.exists(d) && length(list.files(d, all.files = FALSE, no.. = TRUE)) == 0) writeLines(c('sjSDM output folder', 'No file was produced in this folder for the selected settings. See diagnostics/engine_status.json and results/README_sjSDM_results.txt.'), file.path(d, 'README.txt'))",
    "}",
    "find_assembly_pred <- function(name, env_dat = NULL, spatial_dat = NULL, traits_dat = NULL) {",
    "  name <- trimws(as.character(name %||% ''))",
    "  if (!nzchar(name)) return(list(pred = NULL, response = 'sites', source = 'default'))",
    "  if (!is.null(env_dat) && name %in% names(env_dat)) return(list(pred = env_dat[[name]], response = 'sites', source = paste0('env:', name)))",
    "  if (!is.null(spatial_dat) && name %in% names(spatial_dat)) return(list(pred = spatial_dat[[name]], response = 'sites', source = paste0('spatial:', name)))",
    "  if (!is.null(traits_dat) && name %in% names(traits_dat)) return(list(pred = traits_dat[[name]], response = 'species', source = paste0('traits:', name)))",
    "  stop('assembly predictor was requested but not found in env, spatial or traits data: ', name, call. = FALSE)",
    "}",
    "tryCatch({",
    "  if (!requireNamespace('yaml', quietly = TRUE)) stop('Package yaml is required for this reproducible script.', call. = FALSE)",
    "  if (!requireNamespace('sjSDM', quietly = TRUE)) stop('Package sjSDM is required for real fitting.', call. = FALSE)",
    "  if (!requireNamespace('reticulate', quietly = TRUE)) stop('Package reticulate is required for the sjSDM Python/PyTorch backend.', call. = FALSE)",
    "  sjsdm_python <- configure_sjsdm_python()",
    "  suppressWarnings(suppressPackageStartupMessages(library(sjSDM)))",
    "  if (!isTRUE(sjSDM::is_torch_available())) stop('sjSDM cannot import PyTorch from RETICULATE_PYTHON: ', sjsdm_python, call. = FALSE)",
    "  reticulate::py_run_string('import torch, torch_optimizer, pyro, madgrad')",
    "  cfg <- yaml::read_yaml('used_config.yml')",
    "  write_runtime_diagnostics()",
    "  Y_raw <- read_csv_safe(file.path('data','Y.csv'))",
    "  if (is.null(Y_raw)) stop('data/Y.csv is required.', call. = FALSE)",
    "  Y <- as.matrix(data.frame(lapply(Y_raw, function(z) suppressWarnings(as.numeric(z))), check.names = FALSE))",
    "  rownames(Y) <- rownames(Y_raw); colnames(Y) <- colnames(Y_raw)",
    "  if (anyNA(Y)) stop('sjSDM training Y must be numeric and complete. Missing values are only appropriate for conditional prediction, not fitting.', call. = FALSE)",
    "  env_raw <- read_csv_safe(file.path('data','env.csv'))",
    "  env_kind <- tolower(as.character(cfg$model$env_model %||% 'linear'))",
    "  if (env_kind %in% c('none','intercept-only','intercept_only')) { env_dat <- data.frame(intercept = rep(1, nrow(Y))); env_formula <- ~ 1 } else { env_dat <- clean_df(env_raw); if (is.null(env_dat)) stop('env.csv is required unless env module is intercept-only.', call. = FALSE); env_formula <- as.formula(cfg$model$env_formula %||% '~ .') }",
    "  if (env_kind == 'dnn') {",
    "    env_obj <- sjSDM::DNN(data = env_dat, formula = env_formula, hidden = hidden_vec(cfg$dnn_optimizer$hidden), activation = cfg$dnn_optimizer$activation %||% 'selu', bias = as_bool(cfg$dnn_optimizer$bias, TRUE), lambda = as_num(cfg$regularization_biotic$env_lambda, 0), alpha = as_num(cfg$regularization_biotic$env_alpha, 0.5), dropout = as_num(cfg$dnn_optimizer$dropout, 0))",
    "  } else {",
    "    env_obj <- sjSDM::linear(data = env_dat, formula = env_formula, lambda = as_num(cfg$regularization_biotic$env_lambda, 0), alpha = as_num(cfg$regularization_biotic$env_alpha, 0.5))",
    "  }",
    "  spatial_raw <- read_csv_safe(file.path('data','spatial.csv'))",
    "  spatial_kind <- tolower(as.character(cfg$model$spatial_model %||% 'none'))",
    "  if (spatial_kind %in% c('spatial eigenvectors','spatial_eigenvectors')) spatial_kind <- 'eigenvectors'",
    "  spatial_obj <- NULL",
    "  spatial_dat <- NULL",
    "  if (!(spatial_kind %in% c('none',''))) {",
    "    if (is.null(spatial_raw)) stop('spatial.csv is required for a sjSDM spatial module.', call. = FALSE)",
    "    spatial_dat <- clean_df(spatial_raw)",
    "    if (as_bool(cfg$spatial_anova_metacommunity$generateSpatialEV, FALSE)) {",
    "      coords <- as.matrix(data.frame(lapply(spatial_dat, function(z) suppressWarnings(as.numeric(z))), check.names = FALSE))",
    "      if (ncol(coords) < 2 || any(!is.finite(coords[, 1:2, drop=FALSE]))) stop('generateSpatialEV needs at least two numeric coordinate columns.', call. = FALSE)",
    "      spatial_dat <- as.data.frame(sjSDM::generateSpatialEV(coords[, 1:2, drop=FALSE], threshold = as_num(cfg$spatial_anova_metacommunity$spatial_ev_threshold, 0)))",
    "      k <- min(ncol(spatial_dat), as_int(cfg$spatial_anova_metacommunity$spatial_ev_k, ncol(spatial_dat)))",
    "      spatial_dat <- spatial_dat[, seq_len(k), drop = FALSE]",
    "    }",
    "    spatial_formula <- as.formula(cfg$model$spatial_formula %||% '~ 0 + .')",
    "    if (spatial_kind == 'dnn') spatial_obj <- sjSDM::DNN(data = spatial_dat, formula = spatial_formula, hidden = hidden_vec(cfg$dnn_optimizer$hidden), activation = cfg$dnn_optimizer$activation %||% 'selu', bias = as_bool(cfg$dnn_optimizer$bias, TRUE), lambda = as_num(cfg$regularization_biotic$spatial_lambda, 0), alpha = as_num(cfg$regularization_biotic$spatial_alpha, 0.5), dropout = as_num(cfg$dnn_optimizer$dropout, 0))",
    "    else spatial_obj <- sjSDM::linear(data = spatial_dat, formula = spatial_formula, lambda = as_num(cfg$regularization_biotic$spatial_lambda, 0), alpha = as_num(cfg$regularization_biotic$spatial_alpha, 0.5))",
    "  }",
    "  traits_dat <- clean_df(read_csv_safe(file.path('data','traits.csv')))",
    "  groups_dat <- clean_df(read_csv_safe(file.path('data','species_groups.csv')))",
    "  folds_dat <- read_csv_safe(file.path('data','folds.csv'))",
    "  if (!is.null(traits_dat)) write.csv(traits_dat, file.path('tables','traits_metadata.csv'))",
    "  if (!is.null(groups_dat)) write.csv(groups_dat, file.path('tables','species_groups.csv'), row.names = FALSE)",
    "  if (!is.null(folds_dat)) write.csv(folds_dat, file.path('tables','folds_uploaded.csv'), row.names = FALSE)",
    "  write_input_manifest(Y, env_dat, spatial_dat, traits_dat, groups_dat, folds_dat)",
    "  if (!as_bool(cfg$spatial_anova_metacommunity$include_space_in_anova, TRUE) && !is.null(spatial_obj)) {",
    "    msg <- 'include_space_in_anova=FALSE noted: sjSDM::anova has no API flag to remove a fitted spatial module, so the fitted spatial fraction is still reported. To omit space from ANOVA, fit spatial module = none.'",
    "    step_warnings <- c(step_warnings, msg)",
    "    writeLines(msg, file.path('anova','space_fraction_note.txt'))",
    "  }",
    "  control <- sjSDM::sjSDMControl(optimizer = optimizer_from_cfg(cfg$dnn_optimizer$optimizer, cfg$dnn_optimizer$weight_decay), scheduler = as_int(cfg$control$scheduler, 0), lr_reduce_factor = as_num(cfg$control$lr_reduce_factor, 0.99), early_stopping_training = as_int(cfg$control$early_stopping_training, 0), mixed = as_bool(cfg$control$mixed, FALSE))",
    "  biotic_df <- suppressWarnings(as.numeric(cfg$regularization_biotic$biotic_df)); if (!length(biotic_df) || !is.finite(biotic_df)) biotic_df <- NULL else biotic_df <- as.integer(biotic_df)",
    "  biotic_obj <- sjSDM::bioticStruct(df = biotic_df, lambda = as_num(cfg$regularization_biotic$biotic_lambda, 0), alpha = as_num(cfg$regularization_biotic$biotic_alpha, 0.5), on_diag = as_bool(cfg$regularization_biotic$on_diag, FALSE), reg_on_Cov = as_bool(cfg$regularization_biotic$reg_on_Cov, TRUE), inverse = as_bool(cfg$regularization_biotic$inverse, FALSE))",
    "  model <- sjSDM::sjSDM(Y = Y, env = env_obj, biotic = biotic_obj, spatial = spatial_obj, family = family_from_cfg(cfg$model$family), iter = as_int(cfg$model$iter, 100), step_size = min(as_int(cfg$model$step_size, 50), nrow(Y)), learning_rate = as_num(cfg$dnn_optimizer$learning_rate, 0.003), se = as_bool(cfg$model$se, FALSE), sampling = as_int(cfg$model$sampling, 5000), parallel = as_int(cfg$model$parallel, 0), control = control, device = as.character(cfg$dnn_optimizer$device %||% 'cpu'), dtype = as.character(cfg$model$dtype %||% 'float32'), seed = as_int(cfg$model$seed, 1234), verbose = as_bool(cfg$model$verbose, TRUE))",
    "  saveRDS(model, file.path('models','sjSDM_model.rds'))",
    "  if (file.exists(file.path('data','pretrained_weights.rds'))) safe_step('setWeights', { pretrained <- readRDS(file.path('data','pretrained_weights.rds')); model <<- sjSDM::setWeights(model, pretrained); saveRDS(model, file.path('models','sjSDM_model_after_setWeights.rds')); write.csv(data.frame(status='applied', source='data/pretrained_weights.rds'), file.path('weights','pretrained_weights_status.csv'), row.names = FALSE) })",
    "  if (as_bool(cfg$regularization_biotic$tune_regularization, FALSE) && as_int(cfg$regularization_biotic$cv_k, 0) >= 2 && as_int(cfg$regularization_biotic$tune_steps, 0) > 0) safe_step('sjSDM_cv', { cv_res <- tryCatch(sjSDM::sjSDM_cv(Y = Y, env = env_obj, biotic = biotic_obj, spatial = spatial_obj, tune = 'random', CV = as_int(cfg$regularization_biotic$cv_k, 2), tune_steps = as_int(cfg$regularization_biotic$tune_steps, 1), device = as.character(cfg$dnn_optimizer$device %||% 'cpu'), sampling = max(100L, min(as_int(cfg$model$sampling, 5000), 1000L)), family = family_from_cfg(cfg$model$family), iter = as_int(cfg$model$iter, 100), step_size = min(as_int(cfg$model$step_size, 50), nrow(Y)), learning_rate = as_num(cfg$dnn_optimizer$learning_rate, 0.003), parallel = as_int(cfg$model$parallel, 0), control = control, dtype = as.character(cfg$model$dtype %||% 'float32'), seed = as_int(cfg$model$seed, 1234)), error = function(e) e); if (inherits(cv_res, 'error')) { msg <- conditionMessage(cv_res); step_warnings <<- c(step_warnings, paste('sjSDM_cv optional tuning failed; main fit still completed:', msg)); writeLines(msg, file.path('diagnostics','sjSDM_cv_error.txt')); write.csv(data.frame(status='sjSDM_cv_failed_main_fit_completed', error=msg, CV=as_int(cfg$regularization_biotic$cv_k, 2), tune_steps=as_int(cfg$regularization_biotic$tune_steps, 1), stringsAsFactors = FALSE), file.path('tables','sjSDM_cv_result_long.csv'), row.names = FALSE) } else { saveRDS(cv_res, file.path('tables','sjSDM_cv_result.rds')); write.csv(nested_to_long(cv_res), file.path('tables','sjSDM_cv_result_long.csv'), row.names = FALSE) } })",
    "  cf <- safe_step('coef', coef(model))",
    "  env_coef <- NULL",
    "  if (!is.null(cf) && identical(env_kind, 'dnn') && as_bool(cfg$outputs$coef, TRUE)) {",
    "    dnn_obj <- if (is.list(cf) && !is.null(names(cf)) && 'env' %in% names(cf)) cf$env else if (is.list(cf)) cf[[1]] else cf",
    "    safe_step('coef_environment_dnn_weights', write_dnn_coef_artifacts(dnn_obj, 'environment'))",
    "    write.csv(data.frame(engine='sjSDM', response_id=NA_character_, predictor=NA_character_, direction=NA_character_, estimate=NA_real_, lower=NA_real_, upper=NA_real_, notes='environment DNN weights saved in weights/environment_dnn_coef_weights.rds; no direct species-by-predictor coefficient matrix is available', stringsAsFactors = FALSE), file.path('standard','effects_long.csv'), row.names = FALSE)",
    "  } else if (!is.null(cf)) {",
    "    env_coef <- safe_step('coef_environment', coef_block(cf, model, 1))",
    "  }",
    "  if (!is.null(env_coef) && as_bool(cfg$outputs$coef, TRUE)) { write.csv(env_coef, file.path('tables','coef_environment.csv')); write_effects_long(env_coef) }",
    "  if (!is.null(spatial_obj) && !identical(spatial_kind, 'dnn') && !is.null(cf) && is.list(cf) && length(cf) > 1 && as_bool(cfg$outputs$coef, TRUE)) safe_step('coef_spatial', write.csv(coef_block(cf, model, 2, col_names = names(spatial_dat)), file.path('tables','coef_spatial.csv')))",
    "  if (!is.null(spatial_obj) && identical(spatial_kind, 'dnn') && !is.null(cf) && is.list(cf)) { spatial_dnn_obj <- if (!is.null(names(cf)) && 'spatial' %in% names(cf)) cf$spatial else if (length(cf) > 1) cf[[2]] else NULL; if (!is.null(spatial_dnn_obj)) safe_step('coef_spatial_dnn_weights', write_dnn_coef_artifacts(spatial_dnn_obj, 'spatial')) }",
    "  cov_mat <- safe_step('covariance', name_species_matrix(sjSDM::getCov(model), model)); if (!is.null(cov_mat) && as_bool(cfg$outputs$covariance_correlation, TRUE)) { write.csv(cov_mat, file.path('tables','covariance_matrix.csv')); write_associations_long(cov_mat, 'covariance') }",
    "  cor_mat <- safe_step('correlation', name_species_matrix(sjSDM::getCor(model), model)); if (!is.null(cor_mat) && as_bool(cfg$outputs$covariance_correlation, TRUE)) write.csv(cor_mat, file.path('tables','correlation_matrix.csv'))",
    "  if (as_bool(cfg$outputs$weights, TRUE)) safe_step('weights', saveRDS(sjSDM::getWeights(model), file.path('weights','model_weights.rds')))",
    "  if (as_bool(cfg$model$se, FALSE)) safe_step('standard_errors', { se_obj <- sjSDM::getSe(model, step_size = min(as_int(cfg$model$step_size, 50), nrow(Y)), parallel = as_int(cfg$model$parallel, 0)); saveRDS(se_obj, file.path('tables','standard_errors.rds')); se_long <- nested_to_long(se_obj); write.csv(se_long, file.path('tables','standard_errors.csv'), row.names = FALSE); write.csv(data.frame(component='not_available', item='p_value', value=NA_real_, note='sjSDM 1.0.7 getSe returns standard-error objects; p-values are not a stable exported API output.', stringsAsFactors = FALSE), file.path('tables','p_values.csv'), row.names = FALSE) })",
    "  if (as_bool(cfg$outputs$predict, TRUE)) { pred_note <- data.frame(status='not_attempted', reason='', fallback='', stringsAsFactors=FALSE); tryCatch({ pred_train <- predict(model, type = 'link'); write.csv(pred_train, file.path('predictions','predictions.csv')); write_predictions_long(pred_train, Y); nd <- read_csv_safe(file.path('data','newdata.csv')); spn <- read_csv_safe(file.path('data','new_spatial.csv')); nd_clean <- clean_new_df(nd, env_dat); spn_clean <- clean_new_df(spn, spatial_dat); if (!is.null(nd_clean) || !is.null(spn_clean)) { if (!is.null(spatial_obj) && !is.null(nd_clean) && is.null(spn_clean)) stop('new_spatial.csv is required when predicting new rows from a spatial sjSDM model.', call. = FALSE); if (!is.null(spatial_obj) && !is.null(spn_clean) && isTRUE(as_bool(cfg$spatial_anova_metacommunity$generateSpatialEV, FALSE)) && !identical(names(spn_clean), names(spatial_dat))) stop('Prediction with generateSpatialEV needs precomputed new_spatial eigenvectors with the same columns as the training spatial design; raw new coordinates cannot be projected safely by this smoke-test adapter.', call. = FALSE); pred_new <- predict(model, newdata = nd_clean, SP = spn_clean, type = 'link'); write.csv(pred_new, file.path('predictions','newdata_predictions.csv')); write_predictions_long(pred_new, NULL, row_ids = rownames(nd_clean), out_file = file.path('predictions','newdata_predictions_long.csv')) }; pred_note <- data.frame(status='saved_training_probability_and_optional_newdata_probability', reason='', fallback='', stringsAsFactors=FALSE) }, error=function(e) { msg <- conditionMessage(e); writeLines(msg, file.path('diagnostics','sjSDM_predict_note.txt')); pred <- predict(model, type='link'); write.csv(pred, file.path('predictions','predictions.csv')); write_predictions_long(pred, Y); pred_note <<- data.frame(status='saved_training_probability_newdata_skipped', reason=msg, fallback='predictions/predictions.csv', stringsAsFactors=FALSE) }); write.csv(pred_note, file.path('predictions','sjSDM_predict_manifest.csv'), row.names=FALSE) }",
    "  if (as_bool(cfg$outputs$Rsquared, TRUE)) safe_step('Rsquared', { r2 <- sjSDM::Rsquared(model, verbose = as_bool(cfg$model$verbose, TRUE)); r2_num <- as.numeric(r2); r2_names <- names(r2_num); if (is.null(r2_names) || length(r2_names) != length(r2_num) || any(!nzchar(r2_names))) r2_names <- paste0('Rsquared_', seq_along(r2_num)); r2_tab <- data.frame(metric=r2_names, value=r2_num, stringsAsFactors = FALSE); write.csv(r2_tab, file.path('tables','Rsquared.csv'), row.names = FALSE); write.csv(r2_tab, file.path('tables','Rsquared_total.csv'), row.names = FALSE); write.csv(data.frame(response_id = model$species %||% colnames(Y), metric = 'Rsquared_species_not_exported_by_sjSDM_1.0.7', value = NA_real_, stringsAsFactors = FALSE), file.path('tables','Rsquared_species.csv'), row.names = FALSE); write.csv(data.frame(site_id = rownames(Y) %||% paste0('site_', seq_len(nrow(Y))), metric = 'Rsquared_sites_not_exported_by_sjSDM_1.0.7', value = NA_real_, stringsAsFactors = FALSE), file.path('tables','Rsquared_sites.csv'), row.names = FALSE); write.csv(data.frame(engine='sjSDM', metric=r2_names, response_id=NA_character_, value=r2_num, notes='sjSDM::Rsquared'), file.path('standard','fit_metrics.csv'), row.names = FALSE) })",
    "  if (as_bool(cfg$outputs$residuals, TRUE)) safe_step('residuals', { res <- residuals(model); write.csv(res, file.path('tables','residuals.csv')); write.csv(res, file.path('residuals','residuals.csv')) })",
    "  if (as_bool(cfg$spatial_anova_metacommunity$do_anova, TRUE)) safe_step('anova', {",
    "    an <- anova(model, samples = as_int(cfg$spatial_anova_metacommunity$anova_samples, 5000), verbose = as_bool(cfg$model$verbose, TRUE))",
    "    saveRDS(an, file.path('anova','sjSDM_anova.rds'))",
    "    an_sum <- summary(an, fractions = cfg$spatial_anova_metacommunity$internal_fractions %||% 'proportional')",
    "    write.csv(an_sum, file.path('anova','sjSDM_anova_summary.csv'))",
    "    if (!is.null(an$results)) write.csv(an$results, file.path('anova','sjSDM_anova_results.csv'), row.names = FALSE)",
    "    if (!is.null(an$species)) write.csv(nested_to_long(an$species), file.path('anova','anova_species.csv'), row.names = FALSE)",
    "    if (!is.null(an$sites)) write.csv(nested_to_long(an$sites), file.path('anova','anova_sites.csv'), row.names = FALSE)",
    "    if (as_bool(cfg$outputs$plots, TRUE)) { grDevices::pdf(file.path('anova','anova_plot.pdf'), width = 8, height = 6); vals <- if ('R2 McFadden' %in% names(an_sum)) an_sum[['R2 McFadden']] else an_sum[[1]]; graphics::barplot(as.numeric(vals), names.arg = rownames(an_sum), las = 2, col = '#6FA8DC', main = 'sjSDM ANOVA fractions'); grDevices::dev.off() }",
    "    if (as_bool(cfg$spatial_anova_metacommunity$do_internal, TRUE) && !is.null(spatial_obj)) {",
    "      internal <- sjSDM::internalStructure(an, fractions = cfg$spatial_anova_metacommunity$internal_fractions %||% 'proportional')",
    "      saveRDS(internal, file.path('internal_structure','internal_structure.rds'))",
    "      if (!is.null(internal$internals$Species)) write.csv(internal$internals$Species, file.path('internal_structure','internal_structure_species.csv'))",
    "      if (!is.null(internal$internals$Sites)) write.csv(internal$internals$Sites, file.path('internal_structure','internal_structure_sites.csv'))",
    "      write.csv(nested_to_long(internal), file.path('internal_structure','internal_structure_long.csv'), row.names = FALSE)",
    "      if (as_bool(cfg$outputs$plots, TRUE)) { grDevices::pdf(file.path('internal_structure','internal_structure_plot.pdf'), width = 8, height = 6); tabp <- if (!is.null(internal$internals$Species)) internal$internals$Species else internal$internals$Sites; num_cols <- names(tabp)[vapply(tabp, is.numeric, logical(1))]; if (length(num_cols) > 0) graphics::barplot(as.numeric(tabp[[num_cols[[1]]]]), names.arg = rownames(tabp), las = 2, col = '#8DD3C7', main = 'sjSDM internal structure') else plot.new(); grDevices::dev.off() }",
    "      if (as_bool(cfg$spatial_anova_metacommunity$do_assembly, FALSE)) tryCatch({ ap <- find_assembly_pred(cfg$spatial_anova_metacommunity$assembly_predictor, env_dat, spatial_dat, traits_dat); grDevices::pdf(file.path('internal_structure','assembly_effects.pdf'), width = 8, height = 6); try(sjSDM::plotAssemblyEffects(internal, response = ap$response, pred = ap$pred), silent = TRUE); grDevices::dev.off(); write.csv(data.frame(status='plotted', response=ap$response, predictor_source=ap$source, stringsAsFactors = FALSE), file.path('internal_structure','assembly_effects.csv'), row.names = FALSE) }, error = function(e) { writeLines(conditionMessage(e), file.path('diagnostics','sjSDM_assembly_error.txt')); write.csv(data.frame(status='skipped', reason=conditionMessage(e), stringsAsFactors = FALSE), file.path('internal_structure','assembly_effects.csv'), row.names = FALSE) })",
    "    } else if (as_bool(cfg$spatial_anova_metacommunity$do_internal, TRUE) && is.null(spatial_obj)) {",
    "      msg <- 'internalStructure skipped: sjSDM internalStructure is only supported for spatial models.'",
    "      write.csv(data.frame(status='skipped', reason=msg, stringsAsFactors = FALSE), file.path('internal_structure','internal_structure_skipped.csv'), row.names = FALSE)",
    "    } else {",
    "      write.csv(data.frame(status='skipped', reason='do_internal is FALSE', stringsAsFactors = FALSE), file.path('internal_structure','internal_structure_skipped.csv'), row.names = FALSE)",
    "    }",
    "  }) else write.csv(data.frame(status='skipped', reason='do_anova is FALSE', stringsAsFactors = FALSE), file.path('anova','sjSDM_anova_skipped.csv'), row.names = FALSE)",
    "  if (as_bool(cfg$outputs$importance, TRUE) && !is.null(env_coef) && !is.null(cov_mat)) { if (identical(spatial_kind, 'dnn')) { msg <- 'importance skipped: sjSDM getImportance expects linear spatial coefficients; spatial DNN weights were saved in weights/spatial_dnn_coef_weights.rds instead.'; writeLines(msg, file.path('diagnostics','sjSDM_importance_skipped.txt')); write.csv(data.frame(status='skipped', reason=msg), file.path('importance','importance_skipped.csv'), row.names = FALSE); write.csv(data.frame(component='spatial_DNN', item='importance', value=NA, note=msg, stringsAsFactors = FALSE), file.path('importance','importance_summary.csv'), row.names = FALSE) } else safe_step('importance', { imp_fun <- getFromNamespace('getImportance', 'sjSDM'); sp_coef <- if (!is.null(spatial_obj) && !is.null(cf) && is.list(cf) && length(cf) > 1) coef_block(cf, model, 2, col_names = names(spatial_dat)) else NULL; x_mat <- num_matrix(model$data$X); sp_mat <- if (!is.null(sp_coef) && !is.null(spatial_dat)) num_matrix(spatial_dat) else NULL; imp <- imp_fun(beta = t(num_matrix(env_coef)), sp = if (!is.null(sp_coef)) t(num_matrix(sp_coef)) else NULL, association = num_matrix(cov_mat), covX = stats::cov(x_mat), covSP = if (!is.null(sp_mat)) stats::cov(sp_mat) else NULL); saveRDS(imp, file.path('importance','importance.rds')); if (is.list(imp) && !is.null(imp$env)) write.csv(imp$env, file.path('importance','env_importance.csv')); if (is.list(imp) && !is.null(imp$spatial)) write.csv(imp$spatial, file.path('importance','spatial_importance.csv')); if (is.list(imp) && !is.null(imp$biotic)) write.csv(data.frame(response_id = names(imp$biotic) %||% model$species %||% seq_along(imp$biotic), importance = as.numeric(imp$biotic)), file.path('importance','biotic_importance.csv'), row.names = FALSE); write.csv(nested_to_long(imp), file.path('importance','importance_summary.csv'), row.names = FALSE); if (as_bool(cfg$outputs$plots, TRUE)) { grDevices::pdf(file.path('importance','importance_plot.pdf'), width = 8, height = 6); try(plot(imp), silent = TRUE); grDevices::dev.off() } }) } else if (as_bool(cfg$outputs$importance, TRUE)) write.csv(data.frame(status='skipped', reason='importance requires coefficients and covariance outputs', stringsAsFactors = FALSE), file.path('importance','importance_skipped.csv'), row.names = FALSE)",
    "  if (as_bool(cfg$outputs$plots, TRUE)) safe_step('plot', { pdf(file.path('plots','sjSDM_plot.pdf'), width = 8, height = 6); on.exit(dev.off(), add = TRUE); if (!is.null(env_coef)) { m_env <- num_matrix(env_coef); op <- par(mar = c(7, 7, 3, 2)); image(num_matrix(t(m_env[nrow(m_env):1, , drop = FALSE])), axes = FALSE, col = hcl.colors(64, 'Viridis')); axis(1, at = seq(0, 1, length.out = ncol(m_env)), labels = colnames(m_env), las = 2, cex.axis = 0.7); axis(2, at = seq(0, 1, length.out = nrow(m_env)), labels = rev(rownames(m_env)), las = 2, cex.axis = 0.8); title('sjSDM environmental coefficients'); par(op) }; if (!is.null(cor_mat)) { m_cor <- num_matrix(cor_mat); op <- par(mar = c(6, 6, 3, 2)); image(num_matrix(t(m_cor[nrow(m_cor):1, , drop = FALSE])), axes = FALSE, col = hcl.colors(64, 'Blue-Red 3')); axis(1, at = seq(0, 1, length.out = ncol(m_cor)), labels = colnames(m_cor), las = 2); axis(2, at = seq(0, 1, length.out = nrow(m_cor)), labels = rev(rownames(m_cor)), las = 2); title('sjSDM species correlation'); par(op) } })",
    "  write.csv(data.frame(run_id=basename(getwd()), engine='sjSDM', status='fitted', n_sites=nrow(Y), n_responses=ncol(Y), n_predictors=ncol(env_dat), stringsAsFactors=FALSE), file.path('standard','run_summary.csv'), row.names = FALSE)",
    "  write_run_artifacts('fitted', warnings = step_warnings)",
    "  write_status('fitted', warnings = step_warnings)",
    "  write_simple_report('fitted', warnings = step_warnings)",
    "  fill_empty_output_dirs()",
    "  write.csv(data.frame(file = list.files('.', recursive = TRUE), stringsAsFactors = FALSE), file.path('standard','output_manifest.csv'), row.names = FALSE)",
    "}, error = function(e) {",
    "  msg <- conditionMessage(e)",
    "  writeLines(msg, file.path('diagnostics','sjSDM_reproducible_error.txt'))",
    "  write_run_artifacts('fit_failed', errors = msg)",
    "  write_status('fit_failed', errors = msg)",
    "  write_simple_report('fit_failed', errors = msg)",
    "  fill_empty_output_dirs()",
    "  stop(e)",
    "})"
  )
  writeLines(script, file.path(outdir, "reproducible_script", "run_this_sjSDM_analysis.R"))
  writeLines(c(
    "# Workflow wrapper for sjSDM",
    "source(file.path('reproducible_script', 'run_this_sjSDM_analysis.R'))"
  ), file.path(outdir, "workflow_scripts", "run_sjSDM_workflow.R"))
}

write_folder_readmes <- function(outdir, engine) {
  dirs <- list.dirs(outdir, recursive = TRUE, full.names = TRUE)
  for (d in dirs) {
    if (identical(normalizePath(d, winslash = "/", mustWork = FALSE), normalizePath(outdir, winslash = "/", mustWork = FALSE))) next
    if (length(list.files(d, all.files = FALSE, no.. = TRUE)) == 0) {
      writeLines(c(
        paste0("JSDMWorkbench ", engine, " output folder"),
        "This folder is intentionally present for the selected engine workflow.",
        "No engine-specific file was produced here for the current settings/status."
      ), file.path(d, "README.txt"))
    }
  }
}

write_engine_scaffold_outputs <- function(outdir, engine, cfg, status, Y = NULL, X = NULL) {
  status <- normalize_engine_status(status)
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "models"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "results"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "predictions"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
  engine_clean <- gsub("[^A-Za-z0-9_]+", "_", engine)
  status_value <- status$status %||% "unknown"
  n_sites <- if (is.null(Y)) NA_integer_ else nrow(Y)
  n_responses <- if (is.null(Y)) NA_integer_ else ncol(Y)
  n_predictors <- if (is.null(X)) NA_integer_ else ncol(X)
  write.csv(data.frame(
    engine = engine, status = status_value, n_sites = n_sites,
    n_responses = n_responses, n_predictors = n_predictors,
    production_fit_connected = identical(status_value, "fitted"),
    stringsAsFactors = FALSE
  ), file.path(outdir, "tables", "model_run_specification.csv"), row.names = FALSE)
  write.csv(data.frame(
    object = paste0(engine, " fitted model"),
    expected_path = file.path("models", paste0(engine_clean, "_model.rds")),
    status = if (identical(status_value, "fitted")) "available_when_written_by_adapter" else "not_fitted_model_defined_or_failed",
    note = "The GUI exported settings, diagnostics and executable scripts. Statistical model objects require the production engine adapter.",
    stringsAsFactors = FALSE
  ), file.path(outdir, "models", "model_object_manifest.csv"), row.names = FALSE)
  writeLines(c(
    paste0(engine, " model object"),
    "====================",
    paste0("Status: ", status_value),
    "No fitted model object is saved unless the production fitting adapter completes successfully."
  ), file.path(outdir, "models", paste0(engine_clean, "_MODEL_STATUS.txt")))
  write.csv(data.frame(
    output_type = c("prediction_table", "uncertainty_table"),
    expected_file = c("predictions.csv", "prediction_uncertainty.csv"),
    status = status_value,
    note = "Scaffold row. Real values are written only after production fitting/prediction.",
    stringsAsFactors = FALSE
  ), file.path(outdir, "predictions", "prediction_manifest.csv"), row.names = FALSE)
  write.csv(data.frame(
    plot = c("diagnostic_plot", "parameter_plot", "prediction_plot"),
    expected_file = c("diagnostic.pdf", "parameters.pdf", "predictions.pdf"),
    status = status_value,
    stringsAsFactors = FALSE
  ), file.path(outdir, "plots", "plot_manifest.csv"), row.names = FALSE)
  write.csv(data.frame(
    result = c("configuration", "data_check", "engine_status", "standard_outputs", "reproducible_script"),
    file = c("used_config.yml", "diagnostics/data_check_messages.csv", "diagnostics/engine_status.json",
             "standard/run_summary.csv", paste0("reproducible_script/run_this_", engine_clean, "_analysis.R")),
    required_for_review = TRUE,
    stringsAsFactors = FALSE
  ), file.path(outdir, "results", paste0(engine_clean, "_result_index.csv")), row.names = FALSE)

  if (identical(engine, "jSDM")) {
    dir.create(file.path(outdir, "mcmc"), recursive = TRUE, showWarnings = FALSE)
    write.csv(data.frame(parameter_block = c("sp", "gamma", "latent", "alpha", "V_alpha", "V", "Deviance"),
                         expected_object = paste0("mcmc_", c("sp", "gamma", "latent", "alpha", "V_alpha", "V", "Deviance"), ".rds"),
                         status = status_value), file.path(outdir, "mcmc", "mcmc_manifest.csv"), row.names = FALSE)
    write.csv(data.frame(species = colnames(Y %||% data.frame(sp_1 = numeric())), predictor = "intercept", estimate = NA_real_, lower = NA_real_, upper = NA_real_, status = status_value),
              file.path(outdir, "tables", "beta_summary.csv"), row.names = FALSE)
    write.csv(data.frame(response_1 = "sp_1", response_2 = "sp_2", correlation_type = "residual", estimate = NA_real_, status = status_value),
              file.path(outdir, "tables", "residual_cor_scaffold.csv"), row.names = FALSE)
  } else if (identical(engine, "GJAM")) {
    for (d in c("chains", "sensitivity", "ordination", "missing_data")) dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
    write.csv(data.frame(chain = c("bgibbs", "fgibbs", "sgibbs", "ygibbs"), expected_file = paste0(c("bgibbs", "fgibbs", "sgibbs", "ygibbs"), ".rds"), status = status_value),
              file.path(outdir, "chains", "chain_manifest.csv"), row.names = FALSE)
    write.csv(data.frame(parameter = c("betaMu", "betaSe", "corMu", "sigMu", "fmatrix"), expected_table = paste0(c("betaMu", "betaSe", "corMu", "sigMu", "fmatrix"), ".csv"), status = status_value),
              file.path(outdir, "tables", "gjam_parameter_manifest.csv"), row.names = FALSE)
    write.csv(data.frame(component = c("gjamSensitivity", "IIE", "inverse_prediction"), status = status_value),
              file.path(outdir, "sensitivity", "sensitivity_manifest.csv"), row.names = FALSE)
    write.csv(data.frame(axis = NA_integer_, score = NA_real_, status = status_value),
              file.path(outdir, "ordination", "ordination_scores.csv"), row.names = FALSE)
  } else if (identical(engine, "spOccupancy")) {
    for (d in c("samples", "spatial", "model_assessment")) dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
    write.csv(data.frame(block = c("occurrence_beta", "detection_alpha", "latent_z", "spatial_w", "factor_loadings"),
                         expected_file = c("beta_samples.rds", "alpha_samples.rds", "z_samples.rds", "w_samples.rds", "lambda_samples.rds"),
                         status = status_value), file.path(outdir, "samples", "posterior_sample_manifest.csv"), row.names = FALSE)
    write.csv(data.frame(parameter = c("occurrence", "detection"), mean = NA_real_, lower = NA_real_, upper = NA_real_, Rhat = NA_real_, ESS = NA_real_, status = status_value),
              file.path(outdir, "tables", "summary_beta_alpha_scaffold.csv"), row.names = FALSE)
    write.csv(data.frame(metric = c("PPC", "WAIC", "kfold"), value = NA_real_, status = status_value),
              file.path(outdir, "model_assessment", "assessment_manifest.csv"), row.names = FALSE)
    write.csv(data.frame(spatial_output = c("spatial_random_effects", "spatial_parameters", "SVC_surfaces"), status = status_value),
              file.path(outdir, "spatial", "spatial_manifest.csv"), row.names = FALSE)
  } else if (identical(engine, "sjSDM")) {
    for (d in c("anova", "internal_structure", "importance", "weights", "residuals")) dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
    write.csv(data.frame(component = c("environment", "spatial", "biotic"), coefficient = NA_real_, status = status_value),
              file.path(outdir, "tables", "coef_scaffold.csv"), row.names = FALSE)
    write.csv(data.frame(response_1 = "sp_1", response_2 = "sp_2", covariance = NA_real_, correlation = NA_real_, status = status_value),
              file.path(outdir, "tables", "covariance_correlation_scaffold.csv"), row.names = FALSE)
    write.csv(data.frame(fraction = c("environment", "space", "association", "residual"), value = NA_real_, status = status_value),
              file.path(outdir, "anova", "sjSDM_anova_results.csv"), row.names = FALSE)
    write.csv(data.frame(predictor = colnames(X %||% data.frame(env_1 = numeric())), importance = NA_real_, status = status_value),
              file.path(outdir, "importance", "importance_summary.csv"), row.names = FALSE)
    write.csv(data.frame(weight_object = c("env_weights", "spatial_weights", "model_weights"), status = status_value),
              file.path(outdir, "weights", "weights_manifest.csv"), row.names = FALSE)
    write.csv(data.frame(site_id = NA_character_, response_id = NA_character_, residual = NA_real_, status = status_value),
              file.path(outdir, "residuals", "residuals.csv"), row.names = FALSE)
  } else if (identical(engine, "boral")) {
    for (d in c("jags", "mcmc", "ordination", "residuals", "random_effects", "variable_selection")) dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
    writeLines(c("model {", "  # Placeholder JAGS model file for boral production adapter.", "  # Real JAGS code is generated by boral during fitting.", "}"),
               file.path(outdir, "jags", "jagsboralmodel.txt"))
    write.csv(data.frame(parameter = c("beta", "theta", "lv", "sigma"), median = NA_real_, lower = NA_real_, upper = NA_real_, status = status_value),
              file.path(outdir, "tables", "summary_boral.csv"), row.names = FALSE)
    write.csv(data.frame(axis = NA_integer_, site_score = NA_real_, status = status_value),
              file.path(outdir, "ordination", "lv_scores.csv"), row.names = FALSE)
    write.csv(data.frame(response_1 = "sp_1", response_2 = "sp_2", residual_correlation = NA_real_, status = status_value),
              file.path(outdir, "tables", "residual_cor.csv"), row.names = FALSE)
    write.csv(data.frame(effect = "row_effect", estimate = NA_real_, status = status_value),
              file.path(outdir, "random_effects", "ranef_predictions.csv"), row.names = FALSE)
    write.csv(data.frame(predictor = colnames(X %||% data.frame(env_1 = numeric())), posterior_inclusion_probability = NA_real_, status = status_value),
              file.path(outdir, "variable_selection", "ssvs_posterior_probabilities.csv"), row.names = FALSE)
    write.csv(data.frame(metric = c("DIC", "logLik_conditional", "logLik_marginal"), value = NA_real_, status = status_value),
              file.path(outdir, "tables", "boral_fit_metrics.csv"), row.names = FALSE)
  }
  write_folder_readmes(outdir, engine)
}

create_not_run_zip <- function(engine, project_name = "JSDMWorkbench") {
  outdir <- make_engine_run_dir(paste0(engine, "_not_run"), project_name)
  status <- list(
    engine = engine,
    status = "check_failed",
    runtime_seconds = 0,
    warnings = character(),
    errors = paste0(engine, " has not been run in this session. This diagnostic ZIP was created instead of returning a fake ZIP.")
  )
  write_engine_status(outdir, status)
  write_data_check_messages(outdir, "No run was started before download.")
  write_standard_outputs(outdir, engine, status)
  write_reproducible_stub(outdir, engine)
  writeLines(c(
    paste(engine, "was not run"),
    "====================",
    "This is a real diagnostic ZIP created by the download handler.",
    "Run the engine workflow first to generate model outputs.",
    "Start with diagnostics/engine_status.json and diagnostics/data_check_messages.csv."
  ), file.path(outdir, "results", paste0("README_", gsub("[^A-Za-z0-9]+", "_", engine), "_NOT_RUN.txt")))
  make_zip(outdir)
}

create_check_failed_zip <- function(engine, project_name, cfg = list(), check = NULL,
                                    files = list(), Y = NULL, X = NULL, extra = NULL) {
  outdir <- make_engine_run_dir(paste0(engine, "_check_failed"), project_name)
  tryCatch(yaml::write_yaml(cfg, file.path(outdir, "used_config.yml")), error = function(e) NULL)
  if (length(files) > 0) {
    for (nm in names(files)) copy_upload(files[[nm]], outdir, nm)
  }
  tryCatch(write.csv(data_summary(Y, X, extra = extra), file.path(outdir, "tables", "data_summary.csv"), row.names = FALSE),
           error = function(e) writeLines(conditionMessage(e), file.path(outdir, "diagnostics", paste0(gsub("[^A-Za-z0-9]+", "_", engine), "_data_summary_error.txt"))))
  msgs <- if (!is.null(check) && !is.null(check$messages)) check$messages else paste0(engine, " check failed before messages were created.")
  write_data_check_messages(outdir, msgs)
  status <- list(
    engine = engine,
    status = "check_failed",
    runtime_seconds = 0,
    warnings = msgs,
    errors = paste0(engine, " workflow stopped because input data or settings did not pass validation.")
  )
  write_engine_status(outdir, status)
  write_standard_outputs(outdir, engine, status, Y, X)
  write_reproducible_stub(outdir, engine)
  writeLines(c(
    paste(engine, "check failed"),
    "====================",
    "The application created this diagnostic ZIP instead of returning an empty or fake ZIP.",
    "Open diagnostics/data_check_messages.csv first.",
    "Open diagnostics/engine_status.json for machine-readable status."
  ), file.path(outdir, "results", paste0("README_", gsub("[^A-Za-z0-9]+", "_", engine), "_CHECK_FAILED.txt")))
  list(outdir = outdir, zip = make_zip(outdir), status = status)
}

copy_zip_to_download <- function(zip_path, outdir, engine, file, project_name = "JSDMWorkbench") {
  if ((is.null(zip_path) || !file.exists(zip_path)) && !is.null(outdir) && dir.exists(outdir)) {
    zip_path <- tryCatch(make_zip(outdir), error = function(e) {
      writeLines(conditionMessage(e), file.path(outdir, "diagnostics", paste0(gsub("[^A-Za-z0-9]+", "_", engine), "_zip_error.txt")))
      NULL
    })
  }
  if (is.null(zip_path) || !file.exists(zip_path)) {
    zip_path <- create_not_run_zip(engine, project_name)
  }
  ok <- file.copy(zip_path, file, overwrite = TRUE)
  if (!isTRUE(ok) || !file.exists(file) || file.info(file)$size <= 0) {
    stop(paste0("Download ZIP copy failed for ", engine, "."), call. = FALSE)
  }
  invisible(file)
}

data_summary <- function(Y, X, Tr = NULL, study = NULL, coord = NULL, extra = NULL) {
  add <- function(file, obj, role, required = "No") {
    data.frame(
      File = file,
      Required = required,
      Role = role,
      Rows = if (is.null(obj)) NA_integer_ else nrow(obj),
      Columns = if (is.null(obj)) NA_integer_ else ncol(obj),
      Status = if (is.null(obj)) "Not uploaded or unreadable" else "Readable",
      stringsAsFactors = FALSE
    )
  }
  out <- list(
    add("Y.csv", Y, "Response matrix", "Yes"),
    add("XData.csv", X, "Predictor / covariate table", "Yes"),
    add("TrData.csv", Tr, "Traits or response attributes", "No"),
    add("studyDesign.csv", study, "Grouping variables / random effects", "No"),
    add("coordinates.csv", coord, "Spatial coordinates", "No")
  )
  if (!is.null(extra)) out <- c(out, extra)
  do.call(rbind, out)
}









hmsc_factor_preserve_order <- function(x) {
  x <- as.character(x)
  x[trimws(x) == ""] <- NA
  factor(x, levels = unique(x[!is.na(x)]))
}

hmsc_clean_dataframe_types <- function(df, role = "XData") {
  if (is.null(df)) return(NULL)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  role_key <- tolower(gsub("[^a-z]", "", role %||% ""))
  study_role <- role_key %in% c("studydesign", "study")
  for (nm in names(df)) {
    if (study_role) {
      df[[nm]] <- hmsc_factor_preserve_order(df[[nm]])
      next
    }
    if (is.character(df[[nm]])) {
      x <- trimws(df[[nm]])
      x[x == ""] <- NA
      numeric_x <- suppressWarnings(as.numeric(x))
      if (all(is.na(x) | !is.na(numeric_x))) {
        df[[nm]] <- numeric_x
      } else {
        df[[nm]] <- as.factor(x)
      }
    }
    if (is.logical(df[[nm]])) {
      df[[nm]] <- as.factor(df[[nm]])
    }
  }
  df
}

hmsc_normalize_spatial_method <- function(spatial_method = "NNGP") {
  method <- trimws(as.character(spatial_method %||% "NNGP"))
  method_key <- toupper(method)
  if (method_key %in% c("FULL", "GP", "GAUSSIAN PROCESS", "GAUSSIAN_PROCESS")) return("Full")
  if (method_key %in% c("NNGP", "NEAREST NEIGHBOUR GP", "NEAREST_NEIGHBOR_GP")) return("NNGP")
  if (method_key %in% c("GPP", "GAUSSIAN PREDICTIVE PROCESS", "GAUSSIAN_PREDICTIVE_PROCESS")) return("GPP")
  "NNGP"
}

hmsc_normalize_random_mode <- function(random_mode = "none", spatial_method = "NNGP") {
  mode <- tolower(trimws(as.character(random_mode %||% "none")))
  mode_key <- gsub("[^a-z0-9]+", "_", mode)
  mode_key <- gsub("^_|_$", "", mode_key)
  method <- hmsc_normalize_spatial_method(spatial_method)
  if (mode_key %in% c("", "none", "no_random_level", "no_random_effect")) {
    return(list(type = "none", mode = "none", spatial_method = method))
  }
  if (mode_key %in% c("sample", "group", "grouping", "sample_random_effect", "group_random_effect")) {
    return(list(type = "sample", mode = "sample", spatial_method = method))
  }
  if (mode_key %in% c("spatial_full", "full_spatial", "spatial_gp", "full")) {
    return(list(type = "spatial", mode = "spatial_full", spatial_method = "Full"))
  }
  if (mode_key %in% c("spatial_nngp", "nngp_spatial", "nngp")) {
    return(list(type = "spatial", mode = "spatial_nngp", spatial_method = "NNGP"))
  }
  if (mode_key %in% c("spatial_gpp", "gpp_spatial", "gpp")) {
    return(list(type = "spatial", mode = "spatial_gpp", spatial_method = "GPP"))
  }
  if (mode_key %in% c("spatial", "spatial_sdata")) {
    return(list(type = "spatial", mode = paste0("spatial_", tolower(method)), spatial_method = method))
  }
  list(type = "none", mode = "none", spatial_method = method)
}

hmsc_unit_ids <- function(n, preferred = NULL, prefix = "unit") {
  preferred <- as.character(preferred %||% character())
  if (length(preferred) == n && !anyDuplicated(preferred) && all(nzchar(preferred))) {
    return(preferred)
  }
  sprintf("%s_%03d", prefix, seq_len(n))
}

hmsc_prepare_study_design <- function(studyDesign = NULL, n, units = NULL, ensure_col = "sample") {
  units <- hmsc_unit_ids(n, units, prefix = "sample")
  if (is.null(studyDesign)) {
    studyDesign <- data.frame(sample = units, stringsAsFactors = FALSE)
    rownames(studyDesign) <- units
  } else {
    studyDesign <- as.data.frame(studyDesign, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(studyDesign) != n) {
      stop("studyDesign.csv rows must match Y rows.", call. = FALSE)
    }
    if (is.null(rownames(studyDesign)) || any(!nzchar(rownames(studyDesign)))) {
      rownames(studyDesign) <- units
    }
  }
  if (nzchar(ensure_col %||% "") && !(ensure_col %in% names(studyDesign))) {
    studyDesign[[ensure_col]] <- units
  }
  hmsc_clean_dataframe_types(studyDesign, role = "studyDesign")
}

hmsc_select_coordinate_columns <- function(coord, lon_col = "longitude", lat_col = "latitude") {
  if (is.null(coord)) return(list(ok = FALSE, lon_col = NA_character_, lat_col = NA_character_, message = "coordinates.csv is missing."))
  nm <- names(coord)
  find_col <- function(x) {
    x <- trimws(as.character(x %||% ""))
    if (!nzchar(x)) return(NA_character_)
    hit <- nm[tolower(nm) == tolower(x)]
    if (length(hit) > 0) hit[1] else NA_character_
  }
  requested <- c(find_col(lon_col), find_col(lat_col))
  if (all(!is.na(requested))) {
    return(list(ok = TRUE, lon_col = requested[1], lat_col = requested[2], message = NULL))
  }
  candidates <- list(c("longitude", "latitude"), c("lon", "lat"), c("x", "y"), c("easting", "northing"))
  for (pair in candidates) {
    hit <- c(find_col(pair[1]), find_col(pair[2]))
    if (all(!is.na(hit))) {
      return(list(ok = TRUE, lon_col = hit[1], lat_col = hit[2],
                  message = paste0("Coordinate columns auto-detected as ", hit[1], ", ", hit[2], ".")))
    }
  }
  numeric_cols <- nm[vapply(coord, function(z) {
    zz <- suppressWarnings(as.numeric(z))
    all(is.na(z) | !is.na(zz))
  }, logical(1))]
  if (length(numeric_cols) >= 2) {
    return(list(ok = TRUE, lon_col = numeric_cols[1], lat_col = numeric_cols[2],
                message = paste0("Coordinate columns auto-detected as first numeric columns: ", numeric_cols[1], ", ", numeric_cols[2], ".")))
  }
  list(ok = FALSE, lon_col = NA_character_, lat_col = NA_character_,
       message = paste0("Coordinate columns were not found. Requested ", lon_col, ", ", lat_col,
                        "; available columns: ", paste(nm, collapse = ", "), "."))
}

hmsc_prepare_spatial_matrix <- function(coord, lon_col, lat_col, units) {
  coord <- as.data.frame(coord, stringsAsFactors = FALSE, check.names = FALSE)
  selected <- hmsc_select_coordinate_columns(coord, lon_col, lat_col)
  if (!isTRUE(selected$ok)) stop(selected$message, call. = FALSE)
  xy <- as.matrix(coord[, c(selected$lon_col, selected$lat_col), drop = FALSE])
  suppressWarnings(storage.mode(xy) <- "numeric")
  if (any(!is.finite(xy))) {
    stop("Spatial coordinate columns must be numeric and finite.", call. = FALSE)
  }
  rownames(xy) <- units
  attr(xy, "coordinate_message") <- selected$message
  xy
}

hmsc_make_gpp_knots <- function(xy, sKnot_file = "") {
  sKnot_file <- trimws(as.character(sKnot_file %||% ""))
  candidate_paths <- character()
  if (nzchar(sKnot_file)) {
    candidate_paths <- c(sKnot_file, file.path("data", sKnot_file), file.path("inputs", sKnot_file))
  }
  for (path in candidate_paths) {
    if (file.exists(path)) {
      knots <- read.csv(path, row.names = 1, check.names = FALSE)
      knots <- as.data.frame(knots[, seq_len(min(2, ncol(knots))), drop = FALSE])
      knots[] <- lapply(knots, function(z) suppressWarnings(as.numeric(z)))
      if (nrow(knots) > 0 && ncol(knots) == 2 && all(is.finite(as.matrix(knots)))) return(knots)
    }
  }
  n_knots <- max(2L, min(10L, floor(sqrt(nrow(xy)))))
  knots <- tryCatch(constructKnots(as.data.frame(xy), nKnots = n_knots), error = function(e) NULL)
  if (is.null(knots) || nrow(knots) == 0) {
    stop("GPP spatial model needs knot locations; automatic constructKnots() failed. Choose Spatial Full/NNGP or provide a valid sKnot file.", call. = FALSE)
  }
  as.data.frame(knots)
}

hmsc_optional_number <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || !is.finite(x[1])) return(NULL)
  x[1]
}

hmsc_parse_name_list <- function(x) {
  vals <- trimws(unlist(strsplit(as.character(x %||% ""), "[,;]+")))
  vals[nzchar(vals)]
}

hmsc_flatten_config <- function(x, prefix = "") {
  rows <- list()
  for (nm in names(x)) {
    key <- if (nzchar(prefix)) paste(prefix, nm, sep = ".") else nm
    val <- x[[nm]]
    if (is.list(val) && !inherits(val, "data.frame")) {
      rows <- c(rows, hmsc_flatten_config(val, key))
    } else {
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = key,
        value = paste(as.character(val), collapse = ","),
        stringsAsFactors = FALSE
      )
    }
  }
  rows
}

hmsc_parse_alphapw <- function(x) {
  x <- trimws(as.character(x %||% ""))
  if (!nzchar(x)) return(NULL)
  vals <- suppressWarnings(as.numeric(unlist(strsplit(x, "[,;[:space:]]+"))))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(NULL)
  if (length(vals) %% 2 != 0) {
    stop("alphapw must be supplied as numeric distance,probability pairs.", call. = FALSE)
  }
  matrix(vals, ncol = 2, byrow = TRUE)
}

hmsc_apply_random_level_priors <- function(rL, cfg) {
  pri <- cfg$model$priors %||% list()
  args <- list(rL = rL, setDefault = isTRUE(pri$setDefault))
  for (nm in c("a1", "b1", "a2", "b2")) {
    val <- hmsc_optional_number(pri[[nm]])
    if (!is.null(val)) args[[nm]] <- val
  }
  alphapw <- hmsc_parse_alphapw(pri$alphapw %||% "")
  if (!is.null(alphapw)) args$alphapw <- alphapw
  nfMax <- hmsc_optional_number(cfg$model$nfMax)
  nfMin <- hmsc_optional_number(cfg$model$nfMin)
  if (!is.null(nfMax)) args$nfMax <- as.integer(nfMax)
  if (!is.null(nfMin)) {
    if (!is.null(rL$xDim) && is.finite(rL$xDim) && rL$xDim > 0 && nfMin < 2) nfMin <- 2
    args$nfMin <- as.integer(nfMin)
  }
  tryCatch(do.call(setPriors, args), error = function(e) {
    stop(paste("Invalid HMSC random-level prior settings:", conditionMessage(e)), call. = FALSE)
  })
}

hmsc_normalize_initPar <- function(x) {
  x <- trimws(tolower(as.character(x %||% "")))
  if (identical(x, "fixed effects")) return("fixed effects")
  NULL
}

hmsc_normalize_omega_order <- function(x) {
  key <- trimws(as.character(x %||% "original"))
  key_lower <- tolower(key)
  if (key_lower %in% c("alphabetical", "alpha", "alphabet")) return("alphabet")
  if (key_lower == "aoe") return("AOE")
  if (key_lower == "fpc") return("FPC")
  if (key_lower == "hclust") return("hclust")
  "original"
}

hmsc_quiet_optional <- function(expr, fallback = NULL) {
  msg_file <- tempfile()
  con <- file(msg_file, open = "wt")
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    close(con)
    unlink(msg_file)
  }, add = TRUE)
  tryCatch(suppressWarnings(suppressMessages(force(expr))), error = function(e) fallback)
}

hmsc_build_mcmc_updater <- function(cfg) {
  up <- cfg$mcmc$updater %||% list()
  mode_info <- hmsc_normalize_random_mode(cfg$model$random_mode %||% "none", cfg$model$spatial_method %||% "NNGP")
  advanced_type <- gsub("[^a-z0-9]+", "_", tolower(trimws(as.character(cfg$model$random_level_type %||% "none"))))
  advanced_type <- gsub("^_|_$", "", advanced_type)
  gamma_eta <- isTRUE(up$GammaEta)
  if (advanced_type == "covariate_dependent_xdata" || mode_info$spatial_method %in% c("GPP", "NNGP")) {
    gamma_eta <- FALSE
  }
  list(
    GammaEta = gamma_eta,
    Gamma2 = FALSE
  )
}

hmsc_find_workflow_file <- function(file_name) {
  file_name <- trimws(as.character(file_name %||% ""))
  if (!nzchar(file_name)) return(NULL)
  candidates <- c(file_name, file.path("data", file_name), file.path("inputs", file_name))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) return(NULL)
  candidates[1]
}

hmsc_read_optional_matrix <- function(file_name, expected_rows = NULL, expected_cols = NULL, label = "matrix") {
  path <- hmsc_find_workflow_file(file_name)
  if (is.null(path)) return(NULL)
  mat <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (any(!is.finite(mat))) stop(paste0(label, " must be a numeric matrix: ", path), call. = FALSE)
  if (!is.null(expected_rows) && nrow(mat) != expected_rows) {
    stop(paste0(label, " row count must be ", expected_rows, "."), call. = FALSE)
  }
  if (!is.null(expected_cols) && ncol(mat) != expected_cols) {
    stop(paste0(label, " column count must be ", expected_cols, "."), call. = FALSE)
  }
  mat
}

hmsc_read_optional_dataframe <- function(file_name, expected_rows = NULL, label = "data frame") {
  path <- hmsc_find_workflow_file(file_name)
  if (is.null(path)) return(NULL)
  dat <- read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(expected_rows) && nrow(dat) != expected_rows) {
    stop(paste0(label, " row count must be ", expected_rows, "."), call. = FALSE)
  }
  hmsc_clean_dataframe_types(dat, role = label)
}

hmsc_read_square_numeric_matrix <- function(file_name, label = "matrix") {
  path <- hmsc_find_workflow_file(file_name)
  if (is.null(path)) return(NULL)
  mat <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (any(!is.finite(mat))) stop(paste0(label, " must be numeric and finite: ", path), call. = FALSE)
  if (nrow(mat) != ncol(mat)) stop(paste0(label, " must be a square matrix."), call. = FALSE)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop(paste0(label, " must have row names and column names."), call. = FALSE)
  }
  mat
}

hmsc_align_random_level_matrix <- function(mat, levels, label = "random-level matrix") {
  levels <- as.character(levels)
  if (all(levels %in% rownames(mat)) && all(levels %in% colnames(mat))) {
    mat <- mat[levels, levels, drop = FALSE]
  } else if (nrow(mat) == length(levels)) {
    rownames(mat) <- levels
    colnames(mat) <- levels
  } else {
    stop(paste0(label, " dimensions must match the random-level units."), call. = FALSE)
  }
  mat
}

hmsc_align_random_level_xdata <- function(xData, levels, label = "random-level xData") {
  levels <- as.character(levels)
  if (all(levels %in% rownames(xData))) {
    xData <- xData[levels, , drop = FALSE]
  } else if (nrow(xData) == length(levels)) {
    rownames(xData) <- levels
  } else {
    stop(paste0(label, " rows must match the random-level units."), call. = FALSE)
  }
  xData
}

hmsc_make_positive_definite_correlation <- function(C, eps = 1e-6) {
  C <- as.matrix(C)
  suppressWarnings(storage.mode(C) <- "numeric")
  C[!is.finite(C)] <- 0
  C <- (C + t(C)) / 2
  diag(C) <- 1
  ev <- eigen(C, symmetric = TRUE)
  if (min(ev$values, na.rm = TRUE) <= eps) {
    C <- ev$vectors %*% diag(pmax(ev$values, eps), nrow = length(ev$values)) %*% t(ev$vectors)
    C <- (C + t(C)) / 2
    d <- sqrt(pmax(diag(C), eps))
    C <- C / outer(d, d)
    diag(C) <- 1
  }
  C
}

hmsc_taxonomy_to_C <- function(taxonomy, species_names) {
  taxonomy <- as.data.frame(taxonomy, stringsAsFactors = FALSE, check.names = FALSE)
  nm_low <- tolower(names(taxonomy))
  species_col <- names(taxonomy)[nm_low %in% c("species", "species_name", "species_id", "taxon", "taxon_id", "taxon_name", "response", "response_id")]
  rownames_match <- !is.null(rownames(taxonomy)) && length(setdiff(species_names, rownames(taxonomy))) == 0
  if (!rownames_match && length(species_col) > 0) {
    rownames(taxonomy) <- as.character(taxonomy[[species_col[1]]])
  }
  if (is.null(rownames(taxonomy)) || any(!nzchar(rownames(taxonomy)))) {
    stop("Taxonomy file must have species names either as row names or in a species/taxon column.", call. = FALSE)
  }
  missing_species <- setdiff(species_names, rownames(taxonomy))
  if (length(missing_species) > 0) {
    stop(paste0("Taxonomy file is missing species: ", paste(missing_species, collapse = ", ")), call. = FALSE)
  }
  taxonomy <- taxonomy[species_names, , drop = FALSE]
  rank_order <- c("kingdom", "phylum", "division", "class", "order", "family", "genus", "species")
  rank_cols <- names(taxonomy)[tolower(names(taxonomy)) %in% rank_order]
  rank_cols <- rank_cols[order(match(tolower(rank_cols), rank_order))]
  rank_cols <- setdiff(rank_cols, species_col[1] %||% character())
  if (length(rank_cols) == 0) {
    stop("Taxonomy file must contain at least one taxonomic-rank column such as kingdom, phylum, class, order, family, genus or species.", call. = FALSE)
  }
  ranks <- taxonomy[, rank_cols, drop = FALSE]
  ranks[] <- lapply(ranks, function(z) {
    z <- trimws(as.character(z))
    z[z == ""] <- NA_character_
    z
  })
  n <- length(species_names)
  C <- matrix(0, n, n, dimnames = list(species_names, species_names))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      comparable <- !is.na(unlist(ranks[i, , drop = TRUE])) & !is.na(unlist(ranks[j, , drop = TRUE]))
      if (!any(comparable)) {
        C[i, j] <- if (i == j) 1 else 0
      } else {
        C[i, j] <- mean(unlist(ranks[i, comparable, drop = TRUE]) == unlist(ranks[j, comparable, drop = TRUE]))
      }
    }
  }
  hmsc_make_positive_definite_correlation(C)
}

hmsc_read_phylogeny_or_taxonomy <- function(file_name, species_names) {
  path <- hmsc_find_workflow_file(file_name)
  if (is.null(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv", "tsv")) {
    taxonomy <- if (ext == "tsv") {
      read.delim(path, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
    } else {
      read.csv(path, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
    }
    return(list(C = hmsc_taxonomy_to_C(taxonomy, species_names), source = "taxonomy_table", path = path))
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("The ape package is required to read Newick phylogeny files.", call. = FALSE)
  }
  tree <- tryCatch(ape::read.tree(path), error = function(e) NULL)
  if (is.null(tree) || is.null(tree$tip.label)) {
    stop("Phylogeny file is neither a readable Newick tree nor a supported taxonomy CSV/TSV table.", call. = FALSE)
  }
  missing_species <- setdiff(species_names, tree$tip.label)
  if (length(missing_species) > 0) {
    stop(paste0("Phylogeny tree is missing species: ", paste(missing_species, collapse = ", ")), call. = FALSE)
  }
  if (length(setdiff(tree$tip.label, species_names)) > 0) {
    tree <- ape::keep.tip(tree, species_names)
  }
  list(phyloTree = tree, source = "newick_tree", path = path)
}

hmsc_write_derived_phylogeny_outputs <- function(outdir, phy = NULL, species_names = NULL, Cmat = NULL,
                                                source = NA_character_, input_file = NA_character_) {
  if (is.null(species_names) || !length(species_names)) return(invisible(FALSE))
  derived <- NULL
  mode <- "none"
  note <- "No phylogeny/taxonomy input was used."
  if (!is.null(Cmat)) {
    derived <- hmsc_make_positive_definite_correlation(Cmat)
    mode <- "derived_correlation_matrix"
    note <- "Species-by-species correlation matrix was derived internally from the accepted phylogeny/taxonomy input."
  } else if (!is.null(phy) && !is.null(phy$C)) {
    derived <- hmsc_make_positive_definite_correlation(phy$C)
    mode <- "derived_from_taxonomy_table"
    note <- "Taxonomic ranks were converted internally to a positive-definite species-by-species correlation matrix."
  } else if (!is.null(phy) && !is.null(phy$phyloTree)) {
    derived <- tryCatch({
      tree <- phy$phyloTree
      if (length(setdiff(species_names, tree$tip.label)) == 0) {
        if (length(setdiff(tree$tip.label, species_names)) > 0) tree <- ape::keep.tip(tree, species_names)
        C <- ape::vcv.phylo(tree, corr = TRUE)
        C <- C[species_names, species_names, drop = FALSE]
        hmsc_make_positive_definite_correlation(C)
      } else {
        NULL
      }
    }, error = function(e) NULL)
    mode <- "derived_from_newick_tree"
    note <- "Newick tree was passed to Hmsc and an internal reviewer-facing species correlation matrix was also exported."
  }
  if (is.null(derived)) return(invisible(FALSE))
  rownames(derived) <- species_names
  colnames(derived) <- species_names
  for (subdir in c("data", "tables", "standard")) {
    dir.create(file.path(outdir, subdir), recursive = TRUE, showWarnings = FALSE)
    write.csv(derived, file.path(outdir, subdir, "derived_phylogeny_correlation.csv"))
  }
  write.csv(
    data.frame(
      source = as.character(source %||% NA_character_),
      input_file = as.character(input_file %||% NA_character_),
      derived_file = "derived_phylogeny_correlation.csv",
      mode = mode,
      user_upload_contract = "Upload a Newick tree or taxonomy rank table; do not upload a species-by-species covariance/correlation matrix.",
      note = note,
      stringsAsFactors = FALSE
    ),
    file.path(outdir, "tables", "phylogeny_input_contract.csv"),
    row.names = FALSE
  )
  invisible(TRUE)
}

hmsc_read_design_csv <- function(path, optional_id = c("", "row_id", "site_id")) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  if (ncol(dat) && names(dat)[1] %in% optional_id) {
    ids <- dat[[1]]
    if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
      stop("The first column of the study design must contain unique, nonempty row IDs.", call. = FALSE)
    }
    dat <- dat[-1]
    rownames(dat) <- ids
  }
  if (any(!nzchar(names(dat))) || anyDuplicated(names(dat))) {
    stop("Study design columns must have unique, nonempty names.", call. = FALSE)
  }
  dat
}

hmsc_align_rows <- function(dat, ids, label) {
  if (is.null(dat)) return(NULL)
  if (nrow(dat) != length(ids)) stop(label, " must have one row per response or prediction row.", call. = FALSE)
  current <- rownames(dat)
  implicit <- as.character(seq_along(ids))
  if (!is.null(current) && !identical(current, implicit) && !identical(ids, implicit)) {
    if (!setequal(current, ids)) stop(label, " IDs do not match the response or prediction IDs.", call. = FALSE)
    dat <- dat[ids, , drop = FALSE]
  }
  dat
}

hmsc_group_factor <- function(study, col, label = "studyDesign.csv", training = TRUE) {
  if (is.null(study) || !col %in% names(study)) {
    stop(label, " is missing grouping column '", col, "'.", call. = FALSE)
  }
  groups <- hmsc_factor_preserve_order(trimws(as.character(study[[col]])))
  if (anyNA(groups)) stop(label, " has missing group IDs in '", col, "'.", call. = FALSE)
  if (training && nlevels(groups) < 2L) stop("Grouping column '", col, "' needs at least two distinct units.", call. = FALSE)
  groups
}

hmsc_build_random_level <- function(Ymat, studyDesign = NULL, coord = NULL, cfg) {
  if (!is.null(rownames(Ymat))) {
    studyDesign <- hmsc_align_rows(studyDesign, rownames(Ymat), "studyDesign.csv")
    coord <- hmsc_align_rows(coord, rownames(Ymat), "coordinates.csv")
  }
  primary <- hmsc_build_primary_random_level(Ymat, studyDesign, coord, cfg)
  columns <- as.character(cfg$model$additional_random_levels %||% character())
  if (anyNA(columns) || any(!nzchar(columns)) || anyDuplicated(c(names(primary$ranLevels), columns))) {
    stop("Select one distinct studyDesign.csv column for each grouping level.", call. = FALSE)
  }
  if (length(columns)) {
    if (is.null(primary$studyDesign)) primary$studyDesign <- studyDesign
    primary$summary <- primary$summary[!is.na(primary$summary$random_level), , drop = FALSE]
    if (identical(primary$mode_info$type, "none")) {
      primary$mode_info <- list(type = "sample", mode = "grouping", spatial_method = NA_character_)
    }
    for (col in columns) {
      groups <- hmsc_group_factor(studyDesign, col)
      primary$studyDesign[[col]] <- groups
      rL <- Hmsc::HmscRandomLevel(units = levels(groups))
      primary$ranLevels[[col]] <- hmsc_apply_random_level_priors(rL, cfg)
      primary$summary <- rbind(primary$summary,
        data.frame(random_mode = "grouping", random_level = col, spatial_method = NA_character_, units = nlevels(groups)))
    }
  }
  chosen <- trimws(strsplit(as.character(cfg$model$ranLevelsUsed %||% ""), ",", fixed = TRUE)[[1]])
  chosen <- chosen[nzchar(chosen)]
  if (length(chosen)) {
    missing <- setdiff(chosen, names(primary$ranLevels))
    if (length(missing)) stop("ranLevelsUsed contains undefined levels: ", paste(missing, collapse = ", "), call. = FALSE)
    primary$ranLevels <- primary$ranLevels[unique(chosen)]
    primary$summary <- primary$summary[match(names(primary$ranLevels), primary$summary$random_level), , drop = FALSE]
  }
  if (length(primary$ranLevels)) {
    primary$studyDesign <- hmsc_clean_dataframe_types(primary$studyDesign, role = "studyDesign")
  }
  primary
}

hmsc_build_primary_random_level <- function(Ymat, studyDesign = NULL, coord = NULL, cfg) {
  mode_info <- hmsc_normalize_random_mode(cfg$model$random_mode %||% "none", cfg$model$spatial_method %||% "NNGP")
  n <- nrow(Ymat)
  y_units <- hmsc_unit_ids(n, rownames(Ymat), prefix = "sample")
  warnings <- character()
  advanced_type <- gsub("[^a-z0-9]+", "_", tolower(trimws(as.character(cfg$model$random_level_type %||% "none"))))
  advanced_type <- gsub("^_|_$", "", advanced_type)
  advanced_enabled <- !(advanced_type %in% c("", "none"))
  if (advanced_enabled && mode_info$type != "none") {
    warnings <- c(warnings, "Advanced HmscRandomLevel type was ignored because a primary random-level design was selected.")
  }
  if (advanced_enabled && mode_info$type == "none") {
    unit_col <- trimws(as.character(cfg$model$units_column %||% "sample"))
    if (!nzchar(unit_col)) unit_col <- "sample"
    if (advanced_type == "n_only") {
      N <- as.integer(cfg$model$random_N %||% 1L)
      if (!is.finite(N) || N < 1) stop("Advanced N-only random level requires N >= 1.", call. = FALSE)
      level_ids <- as.character(seq_len(N))
      studyDesign <- hmsc_prepare_study_design(studyDesign, n, units = y_units, ensure_col = unit_col)
      studyDesign[[unit_col]] <- factor(level_ids[((seq_len(n) - 1L) %% N) + 1L], levels = level_ids)
      rL <- HmscRandomLevel(N = N)
      rL <- hmsc_apply_random_level_priors(rL, cfg)
      return(list(studyDesign = studyDesign, ranLevels = setNames(list(rL), unit_col),
                  mode_info = list(type = "advanced", mode = "advanced_N_only", spatial_method = NA_character_),
                  warnings = warnings,
                  summary = data.frame(random_mode = "advanced_N_only", random_level = unit_col, spatial_method = NA_character_, units = N)))
    }
    if (advanced_type == "spatial_sdata") {
      if (is.null(coord)) stop("Advanced spatial sData random level requires coordinates.csv.", call. = FALSE)
      if (nrow(coord) != n) stop("coordinates.csv rows must match Y rows.", call. = FALSE)
      spatial_units <- hmsc_unit_ids(n, rownames(coord) %||% y_units, prefix = "spatial")
      studyDesign <- hmsc_prepare_study_design(studyDesign, n, units = y_units, ensure_col = unit_col)
      studyDesign[[unit_col]] <- factor(spatial_units, levels = spatial_units)
      xy <- hmsc_prepare_spatial_matrix(coord, cfg$model$lon_col %||% "longitude", cfg$model$lat_col %||% "latitude", spatial_units)
      coord_message <- attr(xy, "coordinate_message")
      if (!is.null(coord_message)) warnings <- c(warnings, coord_message)
      method <- hmsc_normalize_spatial_method(cfg$model$sMethod %||% "NNGP")
      longlat <- isTRUE(cfg$model$longlat)
      if (method == "NNGP") {
        nn <- as.integer(cfg$model$nNeighbours %||% 10)
        if (nn < 1) stop("NNGP requires nNeighbours >= 1.", call. = FALSE)
        if (nn >= nrow(xy)) {
          warnings <- c(warnings, paste0("NNGP nNeighbours was reduced from ", nn, " to ", max(1L, nrow(xy) - 1L), " because it must be smaller than the number of spatial units."))
          nn <- max(1L, nrow(xy) - 1L)
        }
        rL <- HmscRandomLevel(sData = xy, sMethod = "NNGP", nNeighbours = nn, longlat = longlat)
      } else if (method == "GPP") {
        knots <- hmsc_make_gpp_knots(xy, cfg$model$sKnot_file %||% "")
        rL <- HmscRandomLevel(sData = xy, sMethod = "GPP", sKnot = knots, longlat = longlat)
      } else {
        rL <- HmscRandomLevel(sData = xy, sMethod = "Full", longlat = longlat)
      }
      rL <- hmsc_apply_random_level_priors(rL, cfg)
      return(list(studyDesign = studyDesign, ranLevels = setNames(list(rL), unit_col),
                  mode_info = list(type = "advanced", mode = "advanced_spatial_sData", spatial_method = method),
                  warnings = warnings,
                  summary = data.frame(random_mode = "advanced_spatial_sData", random_level = unit_col, spatial_method = method, units = nrow(xy))))
    }
    studyDesign <- hmsc_prepare_study_design(studyDesign, n, units = y_units, ensure_col = unit_col)
    studyDesign[[unit_col]] <- hmsc_factor_preserve_order(studyDesign[[unit_col]])
    level_ids <- levels(studyDesign[[unit_col]])
    if (advanced_type == "unstructured_units") {
      rL <- HmscRandomLevel(units = level_ids)
      advanced_mode <- "advanced_units"
    } else if (advanced_type == "distance_matrix_distmat") {
      distMat <- hmsc_read_square_numeric_matrix(cfg$model$distMat_file %||% "", label = "distMat")
      if (is.null(distMat)) stop("Advanced distance-matrix random level requires distMat_file.", call. = FALSE)
      distMat <- hmsc_align_random_level_matrix(distMat, level_ids, label = "distMat")
      rL <- HmscRandomLevel(distMat = distMat)
      advanced_mode <- "advanced_distMat"
    } else if (advanced_type == "covariate_dependent_xdata") {
      xData <- hmsc_read_optional_dataframe(cfg$model$xData_file %||% "", label = "random-level xData")
      if (is.null(xData)) stop("Advanced covariate-dependent random level requires xData_file.", call. = FALSE)
      xData <- hmsc_align_random_level_xdata(xData, level_ids, label = "random-level xData")
      rL <- HmscRandomLevel(xData = xData)
      advanced_mode <- "advanced_xData"
    } else {
      stop(paste0("Unsupported advanced HmscRandomLevel type: ", cfg$model$random_level_type), call. = FALSE)
    }
    rL <- hmsc_apply_random_level_priors(rL, cfg)
    return(list(studyDesign = studyDesign, ranLevels = setNames(list(rL), unit_col),
                mode_info = list(type = "advanced", mode = advanced_mode, spatial_method = NA_character_),
                warnings = warnings,
                summary = data.frame(random_mode = advanced_mode, random_level = unit_col, spatial_method = NA_character_, units = length(level_ids))))
  }
  if (mode_info$type == "none") {
    return(list(studyDesign = NULL, ranLevels = NULL, mode_info = mode_info, warnings = warnings,
                summary = data.frame(random_mode = "none", random_level = NA_character_, spatial_method = NA_character_, units = 0)))
  }
  if (mode_info$type == "sample") {
    col <- trimws(as.character(cfg$model$random_effect_column %||% "sample"))
    if (!nzchar(col)) col <- "sample"
    groups <- hmsc_group_factor(studyDesign, col)
    studyDesign <- hmsc_prepare_study_design(studyDesign, n, units = y_units, ensure_col = col)
    studyDesign[[col]] <- groups
    rL <- HmscRandomLevel(units = levels(studyDesign[[col]]))
    rL <- hmsc_apply_random_level_priors(rL, cfg)
    return(list(studyDesign = studyDesign, ranLevels = setNames(list(rL), col), mode_info = mode_info, warnings = warnings,
                summary = data.frame(random_mode = "sample", random_level = col, spatial_method = NA_character_, units = nlevels(studyDesign[[col]]))))
  }
  if (mode_info$type == "spatial") {
    if (is.null(coord)) stop("Spatial random effect requires coordinates.csv.", call. = FALSE)
    if (nrow(coord) != n) stop("coordinates.csv rows must match Y rows.", call. = FALSE)
    if (!is.null(studyDesign) && nrow(studyDesign) != n) stop("studyDesign.csv rows must match Y rows.", call. = FALSE)
    spatial_col <- "spatial_unit"
    spatial_units <- hmsc_unit_ids(n, rownames(coord) %||% y_units, prefix = "spatial")
    studyDesign <- hmsc_prepare_study_design(studyDesign, n, units = y_units, ensure_col = spatial_col)
    studyDesign[[spatial_col]] <- factor(spatial_units, levels = spatial_units)
    xy <- hmsc_prepare_spatial_matrix(coord, cfg$model$lon_col %||% "longitude", cfg$model$lat_col %||% "latitude", spatial_units)
    coord_message <- attr(xy, "coordinate_message")
    if (!is.null(coord_message)) warnings <- c(warnings, coord_message)
    method <- mode_info$spatial_method
    longlat <- isTRUE(cfg$model$longlat)
    if (method == "NNGP") {
      nn <- as.integer(cfg$model$nNeighbours %||% 10)
      if (nn < 1) stop("NNGP requires nNeighbours >= 1.", call. = FALSE)
      if (nn >= nrow(xy)) {
        warnings <- c(warnings, paste0("NNGP nNeighbours was reduced from ", nn, " to ", max(1L, nrow(xy) - 1L), " because it must be smaller than the number of spatial units."))
        nn <- max(1L, nrow(xy) - 1L)
      }
      rL <- HmscRandomLevel(sData = xy, sMethod = "NNGP", nNeighbours = nn, longlat = longlat)
    } else if (method == "GPP") {
      knots <- hmsc_make_gpp_knots(xy, cfg$model$sKnot_file %||% "")
      rL <- HmscRandomLevel(sData = xy, sMethod = "GPP", sKnot = knots, longlat = longlat)
      if (dir.exists("data")) {
        suppressWarnings(try(write.csv(knots, file.path("data", "generated_sKnot_GPP.csv"), row.names = TRUE), silent = TRUE))
      }
    } else {
      rL <- HmscRandomLevel(sData = xy, sMethod = "Full", longlat = longlat)
    }
    rL <- hmsc_apply_random_level_priors(rL, cfg)
    return(list(studyDesign = studyDesign, ranLevels = setNames(list(rL), spatial_col), mode_info = mode_info, warnings = warnings,
                summary = data.frame(random_mode = mode_info$mode, random_level = spatial_col, spatial_method = method, units = nrow(xy))))
  }
  stop("Unsupported HMSC random effect mode.", call. = FALSE)
}

clean_predictor_types <- function(df, role = "predictors") {
  if (is.null(df)) return(NULL)
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  for (nm in names(df)) {
    if (is.character(df[[nm]])) {
      x <- trimws(df[[nm]])
      x[x == ""] <- NA
      numeric_x <- suppressWarnings(as.numeric(x))
      if (all(is.na(x) | !is.na(numeric_x))) {
        df[[nm]] <- numeric_x
      } else {
        df[[nm]] <- as.factor(x)
      }
    } else if (is.logical(df[[nm]])) {
      df[[nm]] <- as.factor(df[[nm]])
    }
  }
  df
}

hmsc_prediction_site_ids <- function(df, prefix = "new_site") {
  ids <- rownames(df)
  if (is.null(ids) || length(ids) != nrow(df) || any(!nzchar(ids))) {
    ids <- paste0(prefix, "_", seq_len(nrow(df)))
  }
  make.unique(as.character(ids))
}

hmsc_prepare_prediction_xdata <- function(newdata, training_xdata, XFormula = "~ .") {
  if (is.null(newdata)) {
    return(list(XData = NULL, site_ids = character(), messages = character()))
  }
  if (is.null(training_xdata)) stop("HMSC newdata prediction requires fitted training XData.", call. = FALSE)
  nd <- as.data.frame(newdata, stringsAsFactors = FALSE, check.names = FALSE)
  tx <- hmsc_clean_dataframe_types(training_xdata, role = "XData")
  site_ids <- hmsc_prediction_site_ids(nd, "new_site")
  id_cols <- intersect(names(nd), c("site_id", "sample_id", "row_id", "id"))
  if (length(id_cols) > 0) {
    candidate_ids <- as.character(nd[[id_cols[1]]])
    if (length(candidate_ids) == nrow(nd) && all(nzchar(candidate_ids))) {
      site_ids <- make.unique(candidate_ids)
    }
    if (!(id_cols[1] %in% names(tx))) nd[[id_cols[1]]] <- NULL
  }
  rownames(nd) <- site_ids
  f <- tryCatch(as.formula(XFormula %||% "~ ."), error = function(e) NULL)
  if (is.null(f)) stop("HMSC XFormula is invalid, so newdata cannot be checked for prediction.", call. = FALSE)
  formula_txt <- paste(deparse(f), collapse = "")
  required <- if (identical(trimws(formula_txt), "~." ) || identical(trimws(formula_txt), "~ .")) names(tx) else all.vars(f)
  missing <- setdiff(required, names(nd))
  if (length(missing) > 0) {
    stop(paste0("HMSC newdata.csv is missing predictor column(s) required by XFormula: ", paste(missing, collapse = ", ")), call. = FALSE)
  }
  keep <- intersect(names(tx), names(nd))
  nd <- nd[, keep, drop = FALSE]
  for (nm in names(nd)) {
    if (is.factor(tx[[nm]])) {
      vals <- as.character(nd[[nm]])
      bad <- setdiff(unique(vals[!is.na(vals)]), levels(tx[[nm]]))
      if (length(bad) > 0) {
        stop(paste0("HMSC newdata.csv column '", nm, "' has unseen categorical level(s): ", paste(bad, collapse = ", "),
                    ". Use the same factor levels as XData.csv."), call. = FALSE)
      }
      nd[[nm]] <- factor(vals, levels = levels(tx[[nm]]))
    } else if (is.numeric(tx[[nm]]) || is.integer(tx[[nm]])) {
      nd[[nm]] <- suppressWarnings(as.numeric(nd[[nm]]))
    } else {
      nd[[nm]] <- tx[[nm]][0]
    }
  }
  if (any(is.na(nd))) stop("HMSC newdata.csv contains NA after type conversion; Hmsc::predict() does not allow NA in XData.", call. = FALSE)
  list(
    XData = nd,
    site_ids = rownames(nd),
    messages = c(
      paste0("HMSC newdata.csv will be used for out-of-sample prediction: ", nrow(nd), " prediction rows."),
      paste0("HMSC newdata predictors matched to training XData columns: ", paste(names(nd), collapse = ", "))
    )
  )
}

hmsc_prepare_prediction_study_design <- function(new_study, model, new_site_ids) {
  if (is.null(model$rLNames) || length(model$rLNames) == 0) return(NULL)
  sd <- hmsc_align_rows(new_study, new_site_ids, "prediction_studyDesign.csv")
  if (is.null(sd)) sd <- data.frame(row.names = new_site_ids)
  rownames(sd) <- new_site_ids
  for (nm in model$rLNames) {
    if (!nm %in% names(sd) && (model$ranLevels[[nm]]$sDim %||% 0) > 0) sd[[nm]] <- new_site_ids
    sd[[nm]] <- hmsc_group_factor(sd, nm, "prediction_studyDesign.csv", training = FALSE)
  }
  sd[, model$rLNames, drop = FALSE]
}

hmsc_prepare_prediction_ranlevels <- function(model, ranLevels, newcoord, cfg, new_site_ids, studyDesign) {
  if (is.null(model$rLNames) || length(model$rLNames) == 0 || is.null(ranLevels)) return(ranLevels)
  has_spatial <- vapply(ranLevels[model$rLNames], function(rL) isTRUE((rL$sDim %||% 0) > 0), logical(1))
  if (!any(has_spatial)) return(ranLevels)
  newcoord <- hmsc_align_rows(newcoord, new_site_ids, "newcoords.csv")
  out <- ranLevels
  for (nm in names(out)) {
    rL <- out[[nm]]
    if (!isTRUE((rL$sDim %||% 0) > 0)) next
    units <- as.character(studyDesign[[nm]])
    is_old <- units %in% levels(model$dfPi[, nm])
    if (is.null(newcoord) && all(is_old)) next
    if (is.null(newcoord)) stop("newcoords.csv is required for new spatial units.", call. = FALSE)
    if (is.null(rL$s)) stop("This fitted spatial random level does not store coordinates for new units.", call. = FALSE)
    is_longlat <- inherits(rL$s, "SpatialPoints")
    s_old <- if (is_longlat) sp::coordinates(rL$s) else as.matrix(rL$s)
    xy <- hmsc_prepare_spatial_matrix(newcoord, cfg$model$lon_col %||% "longitude",
                                     cfg$model$lat_col %||% "latitude", units)
    if (ncol(s_old) != ncol(xy)) stop("Training and prediction coordinate dimensions differ.", call. = FALSE)
    colnames(xy) <- colnames(s_old)
    if (is_longlat && any(abs(xy) > rep(c(180, 90), each = nrow(xy)))) {
      stop("Prediction longitude and latitude must be valid decimal degrees.", call. = FALSE)
    }
    first <- match(units, units)
    if (!isTRUE(all.equal(unname(xy), unname(xy[first, , drop = FALSE]), check.attributes = FALSE))) {
      stop("Each spatial unit must have consistent coordinates in newcoords.csv.", call. = FALSE)
    }
    if (any(is_old) && !isTRUE(all.equal(unname(xy[is_old, , drop = FALSE]),
                                       unname(s_old[units[is_old], , drop = FALSE]), check.attributes = FALSE))) {
      stop("An existing spatial unit was assigned different coordinates.", call. = FALSE)
    }
    added <- xy[!is_old & !duplicated(units), , drop = FALSE]
    s_all <- rbind(s_old, added)
    rL$s <- if (is_longlat) sp::SpatialPoints(s_all, proj4string = rL$s@proj4string) else s_all
    out[[nm]] <- rL
  }
  out
}

hmsc_predict_new_sites <- function(model, XData, new_study, newcoord, cfg, ids, nParallel = 1) {
  spatial <- any(vapply(model$ranLevels[model$rLNames], function(rl) (rl$sDim %||% 0) > 0, logical(1)))
  sd <- hmsc_prepare_prediction_study_design(new_study, model, ids)
  rl <- hmsc_prepare_prediction_ranlevels(model, model$ranLevels, newcoord, cfg, ids, sd)
  seed <- as.integer(cfg$model$seed %||% 123L)
  withr::with_seed(seed, stats::predict(model, XData = XData, studyDesign = sd, ranLevels = rl,
    expected = TRUE, predictEtaMean = !spatial, predictEtaMeanField = spatial, nParallel = nParallel))
}

hmsc_summarise_prediction_object <- function(pred, site_ids = NULL, species_ids = NULL) {
  if (is.null(pred)) return(NULL)
  if (is.list(pred) && length(pred) > 0 && all(vapply(pred, function(x) is.matrix(x) || is.data.frame(x), logical(1)))) {
    mats <- lapply(pred, function(x) as.matrix(x))
    nr <- nrow(mats[[1]]); nc <- ncol(mats[[1]])
    arr <- array(unlist(lapply(mats, as.numeric)), dim = c(nr, nc, length(mats)))
    dimnames(arr) <- list(rownames(mats[[1]]), colnames(mats[[1]]), names(pred) %||% paste0("sample_", seq_along(mats)))
    pred <- arr
  }
  if (length(dim(pred)) == 3) {
    d <- dim(pred)
    site_n <- length(site_ids %||% character())
    sp_n <- length(species_ids %||% character())
    margin <- c(1, 2)
    if (site_n > 0 && sp_n > 0) {
      if (d[1] == site_n && d[2] == sp_n) margin <- c(1, 2)
      else if (d[2] == site_n && d[3] == sp_n) margin <- c(2, 3)
      else if (d[1] == sp_n && d[2] == site_n) margin <- c(2, 1)
    }
    pred_mean <- apply(pred, margin, mean, na.rm = TRUE)
    pred_lower <- apply(pred, margin, stats::quantile, probs = 0.025, na.rm = TRUE)
    pred_upper <- apply(pred, margin, stats::quantile, probs = 0.975, na.rm = TRUE)
  } else {
    if (is.null(dim(pred))) {
      if (is.null(site_ids) || is.null(species_ids) || length(pred) != length(site_ids) * length(species_ids)) {
        stop("HMSC prediction vector cannot be reshaped because site/species dimensions are unknown.", call. = FALSE)
      }
      pred_mean <- matrix(as.numeric(pred), nrow = length(site_ids), ncol = length(species_ids),
                          dimnames = list(site_ids, species_ids))
    } else {
      pred_mean <- as.matrix(pred)
    }
    pred_lower <- pred_mean
    pred_upper <- pred_mean
  }
  if (!is.null(site_ids) && nrow(pred_mean) == length(site_ids)) {
    rownames(pred_mean) <- rownames(pred_lower) <- rownames(pred_upper) <- site_ids
  }
  if (!is.null(species_ids) && ncol(pred_mean) == length(species_ids)) {
    colnames(pred_mean) <- colnames(pred_lower) <- colnames(pred_upper) <- species_ids
  }
  list(mean = pred_mean, lower = pred_lower, upper = pred_upper)
}

hmsc_plot_prediction_heatmap <- function(pred_mean, main = "HMSC predicted values",
                                         observed = NULL, value_label = "predicted mean") {
  mat <- suppressWarnings(apply(as.matrix(pred_mean), c(1, 2), as.numeric))
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = nrow(as.matrix(pred_mean)))
  rownames(mat) <- rownames(pred_mean) %||% paste0("site_", seq_len(nrow(mat)))
  colnames(mat) <- colnames(pred_mean) %||% paste0("sp_", seq_len(ncol(mat)))
  finite_vals <- mat[is.finite(mat)]
  if (length(finite_vals) == 0) {
    plot.new()
    title(main = main)
    text(0.5, 0.5, paste("No finite", value_label, "values were available."))
    return(invisible(FALSE))
  }
  pick_ticks <- function(n, max_ticks = 14L) {
    if (n <= max_ticks) return(seq_len(n))
    unique(pmax(1L, pmin(n, round(seq(1, n, length.out = max_ticks)))))
  }
  zlim <- range(finite_vals, na.rm = TRUE)
  if (isTRUE(all(zlim >= 0 & zlim <= 1))) zlim <- c(0, 1)
  pal <- grDevices::hcl.colors(101, "YlGnBu", rev = FALSE)
  op <- par(mar = c(7.5, 8.5, 4, 5))
  on.exit(par(op), add = TRUE)
  image(
    x = seq_len(ncol(mat)), y = seq_len(nrow(mat)),
    z = t(mat[nrow(mat):1, , drop = FALSE]),
    col = pal, zlim = zlim, axes = FALSE,
    xlab = "Species", ylab = "Site / prediction point", main = main
  )
  x_at <- pick_ticks(ncol(mat), max_ticks = 18L)
  axis(1, at = x_at, labels = abbreviate(colnames(mat)[x_at], minlength = 8),
       las = 2, cex.axis = 0.7)
  y_labels <- rev(rownames(mat))
  y_at <- seq_along(y_labels)
  if (length(y_at) > 30) y_at <- pick_ticks(length(y_labels), max_ticks = 20L)
  y_labels <- y_labels[y_at]
  axis(2, at = y_at, labels = y_labels, las = 2, cex.axis = 0.65)
  if (nrow(mat) > length(y_at)) {
    mtext(paste0(nrow(mat), " prediction rows; axis labels sampled"), side = 2, line = 7.2, cex = 0.65)
  }
  box()
  mtext(value_label, side = 4, line = 2.5)
  if (!is.null(observed)) {
    obs <- suppressWarnings(apply(as.matrix(observed), c(1, 2), as.numeric))
    if (all(dim(obs) == dim(mat))) {
      residual <- obs - mat
      res_vals <- residual[is.finite(residual)]
      if (length(res_vals) > 0) {
        plot(jitter(as.numeric(mat)), jitter(as.numeric(obs)), pch = 16,
             col = grDevices::adjustcolor("#1f6f78", alpha.f = 0.45),
             xlab = value_label, ylab = "observed", main = paste(main, "- observed vs predicted"))
        abline(0, 1, col = "#b35b1f", lwd = 2)
        legend("topleft", bty = "n",
               legend = paste0("mean residual = ", round(mean(res_vals, na.rm = TRUE), 3)))
      }
    }
  }
  invisible(TRUE)
}

hmsc_plot_newdata_spatial_prediction <- function(pred_mean, newcoord, cfg, site_ids,
                                                 main = "HMSC new-site spatial prediction") {
  if (is.null(newcoord) || is.null(pred_mean)) return(invisible(FALSE))
  xy <- tryCatch(
    hmsc_prepare_spatial_matrix(newcoord, cfg$model$lon_col %||% "longitude",
                                cfg$model$lat_col %||% "latitude", site_ids),
    error = function(e) NULL
  )
  if (is.null(xy) || nrow(xy) == 0) return(invisible(FALSE))
  mat <- suppressWarnings(apply(as.matrix(pred_mean), c(1, 2), as.numeric))
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = nrow(xy))
  if (nrow(mat) != nrow(xy)) return(invisible(FALSE))
  score <- rowMeans(mat, na.rm = TRUE)
  score_label <- "mean across species"
  finite_vals <- score[is.finite(score)]
  if (length(finite_vals) == 0) return(invisible(FALSE))
  if (diff(range(finite_vals, na.rm = TRUE)) == 0 && ncol(mat) > 1) {
    species_sd <- apply(mat, 2, stats::sd, na.rm = TRUE)
    if (any(is.finite(species_sd)) && max(species_sd, na.rm = TRUE) > 0) {
      best_species <- which.max(species_sd)
      score <- mat[, best_species]
      score_label <- paste0("species: ", colnames(mat)[best_species] %||% paste0("sp_", best_species))
      finite_vals <- score[is.finite(score)]
    }
  }
  pal <- grDevices::hcl.colors(101, "YlGnBu", rev = FALSE)
  vr <- range(finite_vals, na.rm = TRUE)
  if (diff(vr) == 0) {
    idx <- rep(51L, length(score))
  } else {
    idx <- pmax(1L, pmin(101L, round(1 + 100 * (score - vr[1]) / diff(vr))))
  }
  point_cex <- if (nrow(xy) > 5000) 0.45 else if (nrow(xy) > 1000) 0.65 else if (nrow(xy) > 250) 1.0 else 2.4
  point_pch <- if (nrow(xy) > 1000) 16 else 21
  if (point_pch == 16) {
    plot(xy[, 1], xy[, 2], pch = point_pch, col = pal[idx], cex = point_cex, asp = 1,
         xlab = colnames(xy)[1] %||% "x", ylab = colnames(xy)[2] %||% "y",
         main = paste0(main, "\n", score_label))
  } else {
    plot(xy[, 1], xy[, 2], pch = point_pch, bg = pal[idx], col = "#1b4d59", cex = point_cex, asp = 1,
         xlab = colnames(xy)[1] %||% "x", ylab = colnames(xy)[2] %||% "y",
         main = paste0(main, "\n", score_label))
  }
  if (nrow(xy) <= 250) text(xy[, 1], xy[, 2], labels = seq_len(nrow(xy)), cex = 0.65, pos = 3)
  legend("topright", bty = "n",
         legend = c(score_label,
                    paste0("points = ", nrow(xy)),
                    paste0("mean = ", round(mean(finite_vals), 3)),
                    paste0("range = ", paste(round(range(finite_vals), 3), collapse = "-"))))
  invisible(TRUE)
}

hmsc_plot_newdata_spatial_score <- function(score, newcoord, cfg, site_ids,
                                            main = "HMSC new-site spatial score",
                                            score_label = "score") {
  if (is.null(newcoord) || is.null(score)) return(invisible(FALSE))
  xy <- tryCatch(
    hmsc_prepare_spatial_matrix(newcoord, cfg$model$lon_col %||% "longitude",
                                cfg$model$lat_col %||% "latitude", site_ids),
    error = function(e) NULL
  )
  if (is.null(xy) || nrow(xy) == 0 || length(score) != nrow(xy)) return(invisible(FALSE))
  score <- suppressWarnings(as.numeric(score))
  finite_vals <- score[is.finite(score)]
  if (length(finite_vals) == 0) return(invisible(FALSE))
  pal <- grDevices::hcl.colors(101, "YlGnBu", rev = FALSE)
  vr <- range(finite_vals, na.rm = TRUE)
  idx <- if (diff(vr) == 0) rep(51L, length(score)) else pmax(1L, pmin(101L, round(1 + 100 * (score - vr[1]) / diff(vr))))
  point_cex <- if (nrow(xy) > 5000) 0.45 else if (nrow(xy) > 1000) 0.65 else if (nrow(xy) > 250) 1.0 else 2.4
  point_pch <- if (nrow(xy) > 1000) 16 else 21
  if (point_pch == 16) {
    plot(xy[, 1], xy[, 2], pch = point_pch, col = pal[idx], cex = point_cex, asp = 1,
         xlab = colnames(xy)[1] %||% "x", ylab = colnames(xy)[2] %||% "y",
         main = paste0(main, "\n", score_label))
  } else {
    plot(xy[, 1], xy[, 2], pch = point_pch, bg = pal[idx], col = "#1b4d59", cex = point_cex, asp = 1,
         xlab = colnames(xy)[1] %||% "x", ylab = colnames(xy)[2] %||% "y",
         main = paste0(main, "\n", score_label))
  }
  if (nrow(xy) <= 250) text(xy[, 1], xy[, 2], labels = seq_len(nrow(xy)), cex = 0.65, pos = 3)
  legend("topright", bty = "n",
         legend = c(score_label,
                    paste0("points = ", nrow(xy)),
                    paste0("mean = ", round(mean(finite_vals), 3)),
                    paste0("range = ", paste(round(range(finite_vals), 3), collapse = "-"))))
  invisible(TRUE)
}

hmsc_prepare_beta_heatmap_data <- function(m, postBeta, include_intercept = FALSE) {
  if (is.null(postBeta) || is.null(postBeta$mean)) return(NULL)
  beta <- as.matrix(postBeta$mean)
  cov_names <- rownames(beta)
  if (is.null(cov_names) || length(cov_names) != nrow(beta) || any(!nzchar(cov_names))) {
    cov_names <- m$covNames %||% colnames(m$X) %||% paste0("covariate_", seq_len(nrow(beta)))
  }
  sp_names <- colnames(beta)
  if (is.null(sp_names) || length(sp_names) != ncol(beta) || any(!nzchar(sp_names))) {
    sp_names <- m$spNames %||% paste0("species_", seq_len(ncol(beta)))
  }
  rownames(beta) <- make.unique(as.character(cov_names))
  colnames(beta) <- make.unique(as.character(sp_names))
  support <- NULL
  if (!is.null(postBeta$support)) {
    support <- as.matrix(postBeta$support)
    if (identical(dim(support), dim(beta))) {
      rownames(support) <- rownames(beta)
      colnames(support) <- colnames(beta)
    } else {
      support <- NULL
    }
  }
  if (!isTRUE(include_intercept) && nrow(beta) > 1) {
    keep <- !tolower(rownames(beta)) %in% c("(intercept)", "intercept")
    if (any(keep)) {
      beta <- beta[keep, , drop = FALSE]
      if (!is.null(support)) support <- support[keep, , drop = FALSE]
    }
  }
  tree <- m$phyloTree %||% NULL
  list(beta = beta, support = support, tree = tree)
}

hmsc_write_beta_heatmap_long <- function(beta, support, path, model_name = "main") {
  if (is.null(beta) || nrow(beta) == 0 || ncol(beta) == 0) return(invisible(FALSE))
  grid <- expand.grid(predictor = rownames(beta), species = colnames(beta), stringsAsFactors = FALSE)
  grid$model <- model_name
  grid$estimate <- as.numeric(beta[cbind(match(grid$predictor, rownames(beta)), match(grid$species, colnames(beta)))])
  grid$support <- if (!is.null(support)) {
    as.numeric(support[cbind(match(grid$predictor, rownames(support)), match(grid$species, colnames(support)))])
  } else {
    NA_real_
  }
  grid <- grid[, c("model", "species", "predictor", "estimate", "support")]
  write.csv(grid, path, row.names = FALSE)
  invisible(TRUE)
}

hmsc_draw_phylo_beta_heatmap <- function(m, postBeta, model_name = "main",
                                         include_intercept = FALSE,
                                         support_level = 0.95,
                                         show_species_names = TRUE) {
  dat <- hmsc_prepare_beta_heatmap_data(m, postBeta, include_intercept = include_intercept)
  if (is.null(dat) || is.null(dat$beta) || nrow(dat$beta) == 0 || ncol(dat$beta) == 0) {
    plot.new()
    title(main = paste0(model_name, ": phylogeny + Beta heatmap"))
    text(0.5, 0.5, "Beta estimates were not available.")
    return(invisible(FALSE))
  }
  beta <- dat$beta
  support <- dat$support
  tree <- dat$tree
  tree_available <- !is.null(tree) && requireNamespace("ape", quietly = TRUE)
  tip_order <- colnames(beta)
  tree_note <- "species order: input order"
  if (tree_available) {
    keep_tips <- intersect(tree$tip.label, colnames(beta))
    if (length(keep_tips) >= 2) {
      tree <- ape::keep.tip(tree, keep_tips)
      beta <- beta[, tree$tip.label, drop = FALSE]
      if (!is.null(support)) support <- support[, tree$tip.label, drop = FALSE]
      tree_note <- "species order: phylogeny"
    } else {
      tree_available <- FALSE
      tree_note <- "phylogeny skipped: fewer than two matching tips"
    }
  }
  if (!tree_available && ncol(beta) > 2) {
    hc <- tryCatch(stats::hclust(stats::dist(t(beta))), error = function(e) NULL)
    if (!is.null(hc)) {
      beta <- beta[, hc$labels[hc$order], drop = FALSE]
      if (!is.null(support)) support <- support[, colnames(beta), drop = FALSE]
      tree_note <- "species order: Beta similarity"
    }
  }
  finite_beta <- beta[is.finite(beta)]
  if (length(finite_beta) == 0) {
    plot.new()
    title(main = paste0(model_name, ": phylogeny + Beta heatmap"))
    text(0.5, 0.5, "No finite Beta estimates were available.")
    return(invisible(FALSE))
  }
  cap <- suppressWarnings(as.numeric(stats::quantile(abs(finite_beta), probs = 0.98, na.rm = TRUE)))
  if (!is.finite(cap) || cap <= 0) cap <- max(abs(finite_beta), na.rm = TRUE)
  if (!is.finite(cap) || cap <= 0) cap <- 1
  beta_clip <- matrix(pmax(-cap, pmin(cap, beta)), nrow = nrow(beta), ncol = ncol(beta),
                      dimnames = dimnames(beta))
  pal <- grDevices::colorRampPalette(c("#2b6f9e", "#f7f7f5", "#b84a3a"))(101)
  graphics::layout(matrix(c(1, 2, 3), nrow = 1), widths = c(if (tree_available) 2.25 else 0.25, 3.4, 0.45))
  on.exit(graphics::layout(1), add = TRUE)
  op <- par(oma = c(0, 0, 2.8, 0), xpd = NA)
  on.exit(par(op), add = TRUE)
  n_sp <- ncol(beta_clip)
  n_cov <- nrow(beta_clip)
  if (tree_available) {
    par(mar = c(7, 0.2, 2.5, 0.1))
    ape::plot.phylo(tree, direction = "rightwards", show.tip.label = FALSE,
                    no.margin = FALSE, edge.color = "#8a8078", edge.width = 0.8)
    last_plot <- get("last_plot.phylo", envir = get(".PlotPhyloEnv", envir = asNamespace("ape")))
    tip_order <- tree$tip.label[order(last_plot$yy[seq_len(ape::Ntip(tree))])]
    beta_clip <- beta_clip[, tip_order, drop = FALSE]
    beta <- beta[, tip_order, drop = FALSE]
    if (!is.null(support)) support <- support[, tip_order, drop = FALSE]
    title("phylogeny", line = 0.4, cex.main = 0.9)
  } else {
    par(mar = c(7, 0.2, 2.5, 0.1))
    plot.new()
    text(0.5, 0.5, "No\nphylogeny", cex = 0.75, col = "#6f6f6f")
  }
  par(mar = c(7, 0.4, 2.5, if (isTRUE(show_species_names) && n_sp <= 80) 5.0 else 0.6))
  image(x = seq_len(n_cov), y = seq_len(n_sp), z = beta_clip,
        axes = FALSE, xlab = "", ylab = "", col = pal,
        zlim = c(-cap, cap), useRaster = TRUE)
  axis(1, at = seq_len(n_cov), labels = rownames(beta_clip), las = 2, cex.axis = 0.62)
  if (isTRUE(show_species_names) && n_sp <= 80) {
    axis(4, at = seq_len(n_sp), labels = colnames(beta_clip), las = 2, cex.axis = 0.42, tick = FALSE, line = -0.2)
  }
  box(col = "#9a9a9a")
  abline(v = seq(0.5, n_cov + 0.5, by = 1), col = grDevices::adjustcolor("white", alpha.f = 0.55), lwd = 0.35)
  if (!is.null(support) && identical(dim(support), dim(beta_clip))) {
    supported <- which(support >= support_level, arr.ind = TRUE)
    if (nrow(supported) > 0) {
      points(supported[, 1], supported[, 2], pch = 16, cex = 0.22, col = grDevices::adjustcolor("black", alpha.f = 0.65))
    }
  }
  title("environmental effects (Beta)", line = 0.4, cex.main = 0.9)
  par(mar = c(7, 0.1, 2.5, 2.2))
  legend_y <- seq(-cap, cap, length.out = length(pal))
  image(x = 1, y = legend_y, z = matrix(legend_y, nrow = 1), axes = FALSE,
        xlab = "", ylab = "", col = pal, zlim = c(-cap, cap), useRaster = TRUE)
  axis(4, at = pretty(c(-cap, cap), n = 5), labels = round(pretty(c(-cap, cap), n = 5), 2), cex.axis = 0.65)
  title("Beta", line = 0.4, cex.main = 0.8)
  mtext(paste0(model_name, ": phylogeny-ordered species-environment Beta heatmap"),
        outer = TRUE, side = 3, line = 1.1, cex = 1.0, font = 2)
  mtext(paste0("Posterior mean Beta; colors clipped at |", round(cap, 2),
               "| for display. Black dots: support >= ", support_level,
               ". ", tree_note, "."),
        outer = TRUE, side = 3, line = -0.3, cex = 0.68)
  invisible(TRUE)
}

hmsc_save_phylo_beta_heatmap <- function(m, postBeta, model_name = "main", result_dir = "results",
                                         include_intercept = FALSE,
                                         support_level = 0.95,
                                         show_species_names = TRUE) {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(result_dir, paste0("parameter_estimates_phylo_beta_heatmap_", model_name))
  ok <- TRUE
  tryCatch({
    grDevices::pdf(paste0(base, ".pdf"), width = 11, height = 7.5)
    hmsc_draw_phylo_beta_heatmap(m, postBeta, model_name = model_name,
                                 include_intercept = include_intercept,
                                 support_level = support_level,
                                 show_species_names = show_species_names)
    grDevices::dev.off()
  }, error = function(e) {
    ok <<- FALSE
    try(grDevices::dev.off(), silent = TRUE)
  })
  tryCatch({
    grDevices::png(paste0(base, ".png"), width = 2600, height = 1700, res = 220)
    hmsc_draw_phylo_beta_heatmap(m, postBeta, model_name = model_name,
                                 include_intercept = include_intercept,
                                 support_level = support_level,
                                 show_species_names = show_species_names)
    grDevices::dev.off()
  }, error = function(e) {
    ok <<- FALSE
    try(grDevices::dev.off(), silent = TRUE)
  })
  invisible(ok)
}

hmsc_prepare_gamma_heatmap_data <- function(m, postGamma, include_intercepts = TRUE) {
  if (is.null(postGamma) || is.null(postGamma$mean)) return(NULL)
  gamma <- as.matrix(postGamma$mean)
  if (nrow(gamma) == 0 || ncol(gamma) == 0) return(NULL)
  predictor_names <- rownames(gamma)
  if (is.null(predictor_names) || length(predictor_names) != nrow(gamma) || any(!nzchar(predictor_names))) {
    predictor_names <- m$covNames %||% NULL
  }
  if (is.null(predictor_names) || length(predictor_names) != nrow(gamma) || any(!nzchar(predictor_names))) {
    predictor_names <- tryCatch(colnames(stats::model.matrix(m$XFormula, data = m$XData)), error = function(e) NULL)
  }
  if (is.null(predictor_names) || length(predictor_names) != nrow(gamma) || any(!nzchar(predictor_names))) {
    predictor_names <- paste0("predictor_term_", seq_len(nrow(gamma)))
  }
  trait_names <- colnames(gamma)
  if (is.null(trait_names) || length(trait_names) != ncol(gamma) || any(!nzchar(trait_names))) {
    trait_names <- m$trNames %||% NULL
  }
  if (is.null(trait_names) || length(trait_names) != ncol(gamma) || any(!nzchar(trait_names))) {
    trait_names <- tryCatch(colnames(stats::model.matrix(m$TrFormula, data = m$TrData)), error = function(e) NULL)
  }
  if (is.null(trait_names) || length(trait_names) != ncol(gamma) || any(!nzchar(trait_names))) {
    trait_names <- paste0("trait_term_", seq_len(ncol(gamma)))
  }
  rownames(gamma) <- make.unique(as.character(predictor_names))
  colnames(gamma) <- make.unique(as.character(trait_names))
  support <- NULL
  if (!is.null(postGamma$support)) {
    support <- as.matrix(postGamma$support)
    if (identical(dim(support), dim(gamma))) {
      rownames(support) <- rownames(gamma)
      colnames(support) <- colnames(gamma)
    } else {
      support <- NULL
    }
  }
  gamma_all <- gamma
  support_all <- support
  if (!isTRUE(include_intercepts)) {
    keep_rows <- !tolower(rownames(gamma)) %in% c("(intercept)", "intercept")
    keep_cols <- !tolower(colnames(gamma)) %in% c("(intercept)", "intercept")
    if (any(keep_rows) && any(keep_cols)) {
      gamma <- gamma[keep_rows, keep_cols, drop = FALSE]
      if (!is.null(support)) support <- support[keep_rows, keep_cols, drop = FALSE]
    }
  }
  list(gamma = gamma, support = support, gamma_all = gamma_all, support_all = support_all)
}

hmsc_write_gamma_heatmap_long <- function(gamma, support, path, model_name = "main") {
  if (is.null(gamma) || nrow(gamma) == 0 || ncol(gamma) == 0) return(invisible(FALSE))
  grid <- expand.grid(predictor = rownames(gamma), trait = colnames(gamma), stringsAsFactors = FALSE)
  grid$model <- model_name
  grid$estimate <- as.numeric(gamma[cbind(match(grid$predictor, rownames(gamma)), match(grid$trait, colnames(gamma)))])
  grid$support <- if (!is.null(support)) {
    as.numeric(support[cbind(match(grid$predictor, rownames(support)), match(grid$trait, colnames(support)))])
  } else {
    NA_real_
  }
  grid$effect_type <- "trait_environment_gamma"
  grid$scale <- "HMSC Gamma posterior mean"
  grid$note <- "Gamma links environmental predictor terms to trait terms; it is a trait-mediated environmental response, not a species-specific Beta coefficient."
  grid <- grid[, c("model", "predictor", "trait", "estimate", "support", "effect_type", "scale", "note")]
  write.csv(grid, path, row.names = FALSE)
  invisible(TRUE)
}

hmsc_pretty_hmsc_terms <- function(x, max_chars = 28) {
  out <- as.character(x)
  out <- gsub("`", "", out, fixed = TRUE)
  out <- gsub("^\\(Intercept\\)$", "Intercept", out)
  out <- gsub("^poly\\(([^,]+),[^)]*\\)([0-9]+)$", "\\1 poly\\2", out)
  out <- gsub("^I\\((.*)\\)$", "\\1", out)
  out <- gsub("\\s+", " ", out)
  too_long <- nchar(out) > max_chars
  if (any(too_long)) out[too_long] <- abbreviate(out[too_long], minlength = max_chars, strict = TRUE)
  out
}

hmsc_draw_trait_environment_gamma_heatmap <- function(m, postGamma, model_name = "main",
                                                     include_intercepts = FALSE,
                                                     support_level = 0.95) {
  dat <- hmsc_prepare_gamma_heatmap_data(m, postGamma, include_intercepts = include_intercepts)
  if (is.null(dat) || is.null(dat$gamma) || nrow(dat$gamma) == 0 || ncol(dat$gamma) == 0) {
    plot.new()
    title(main = paste0(model_name, ": trait-environment response"))
    text(0.5, 0.5, "Gamma estimates were not available.")
    return(invisible(FALSE))
  }
  gamma <- dat$gamma
  support <- dat$support
  finite_gamma <- gamma[is.finite(gamma)]
  if (length(finite_gamma) == 0) {
    plot.new()
    title(main = paste0(model_name, ": trait-environment response"))
    text(0.5, 0.5, "No finite Gamma estimates were available.")
    return(invisible(FALSE))
  }
  cap <- suppressWarnings(as.numeric(stats::quantile(abs(finite_gamma), probs = 0.98, na.rm = TRUE)))
  if (!is.finite(cap) || cap <= 0) cap <- max(abs(finite_gamma), na.rm = TRUE)
  if (!is.finite(cap) || cap <= 0) cap <- 1
  gamma_clip <- matrix(pmax(-cap, pmin(cap, gamma)), nrow = nrow(gamma), ncol = ncol(gamma),
                       dimnames = dimnames(gamma))
  pal <- grDevices::colorRampPalette(c("#2b6f9e", "#f8f8f5", "#b84a3a"))(101)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  layout(matrix(c(1, 2), nrow = 1), widths = c(5.5, 0.55))
  on.exit(layout(1), add = TRUE)
  par(oma = c(0, 0, 2.6, 0), mar = c(8.5, 7.0, 2.4, 0.4), xpd = NA)
  n_pred <- nrow(gamma_clip)
  n_trait <- ncol(gamma_clip)
  gamma_plot <- gamma_clip[n_pred:1, , drop = FALSE]
  image(x = seq_len(n_trait), y = seq_len(n_pred), z = t(gamma_plot),
        axes = FALSE, xlab = "", ylab = "", col = pal,
        zlim = c(-cap, cap), useRaster = TRUE)
  axis(1, at = seq_len(n_trait), labels = hmsc_pretty_hmsc_terms(colnames(gamma_clip), 24), las = 2, cex.axis = 0.72)
  axis(2, at = seq_len(n_pred), labels = hmsc_pretty_hmsc_terms(rev(rownames(gamma_clip)), 28), las = 2, cex.axis = 0.72)
  box(col = "#9a9a9a")
  abline(v = seq(0.5, n_trait + 0.5, by = 1), col = grDevices::adjustcolor("white", alpha.f = 0.55), lwd = 0.35)
  abline(h = seq(0.5, n_pred + 0.5, by = 1), col = grDevices::adjustcolor("white", alpha.f = 0.55), lwd = 0.35)
  if (!is.null(support) && identical(dim(support), dim(gamma_clip))) {
    support_plot <- support[n_pred:1, , drop = FALSE]
    supported <- which(support_plot >= support_level, arr.ind = TRUE)
    if (nrow(supported) > 0) {
      points(supported[, 2], supported[, 1], pch = 16, cex = 0.72,
             col = grDevices::adjustcolor("black", alpha.f = 0.68))
    }
  }
  title("trait-mediated environmental response (Gamma)", line = 0.3, cex.main = 0.95)
  par(mar = c(8.5, 0.2, 2.4, 2.2))
  legend_y <- seq(-cap, cap, length.out = length(pal))
  image(x = 1, y = legend_y, z = matrix(legend_y, nrow = 1), axes = FALSE,
        xlab = "", ylab = "", col = pal, zlim = c(-cap, cap), useRaster = TRUE)
  axis(4, at = pretty(c(-cap, cap), n = 5), labels = round(pretty(c(-cap, cap), n = 5), 2), cex.axis = 0.68)
  title("Gamma", line = 0.3, cex.main = 0.82)
  mtext(paste0(model_name, ": trait-environment response heatmap"),
        outer = TRUE, side = 3, line = 1.0, cex = 1.0, font = 2)
  mtext(paste0("Posterior mean Gamma; colors clipped at |", round(cap, 2),
               "| for display. Black dots: support >= ", support_level,
               ". Intercepts ", if (isTRUE(include_intercepts)) "shown" else "hidden", "."),
        outer = TRUE, side = 3, line = -0.35, cex = 0.68)
  invisible(TRUE)
}

hmsc_save_trait_environment_gamma_heatmap <- function(m, postGamma, model_name = "main", result_dir = "results",
                                                      include_intercepts = FALSE,
                                                      support_level = 0.95) {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(result_dir, paste0("parameter_estimates_trait_environment_gamma_heatmap_", model_name))
  ok <- TRUE
  tryCatch({
    grDevices::pdf(paste0(base, ".pdf"), width = 8.8, height = 7.2)
    hmsc_draw_trait_environment_gamma_heatmap(m, postGamma, model_name = model_name,
                                             include_intercepts = include_intercepts,
                                             support_level = support_level)
    grDevices::dev.off()
  }, error = function(e) {
    ok <<- FALSE
    try(grDevices::dev.off(), silent = TRUE)
  })
  tryCatch({
    grDevices::png(paste0(base, ".png"), width = 2200, height = 1800, res = 220)
    hmsc_draw_trait_environment_gamma_heatmap(m, postGamma, model_name = model_name,
                                             include_intercepts = include_intercepts,
                                             support_level = support_level)
    grDevices::dev.off()
  }, error = function(e) {
    ok <<- FALSE
    try(grDevices::dev.off(), silent = TRUE)
  })
  invisible(ok)
}

hmsc_plot_species_order <- function(m, species_names, value_matrix = NULL) {
  species_names <- make.unique(as.character(species_names))
  tree <- m$phyloTree %||% NULL
  if (!is.null(tree) && requireNamespace("ape", quietly = TRUE)) {
    keep <- intersect(tree$tip.label, species_names)
    if (length(keep) >= 2) {
      tree <- ape::keep.tip(tree, keep)
      return(list(order = tree$tip.label, tree = tree, note = "phylogeny"))
    }
  }
  if (!is.null(value_matrix) && length(species_names) > 2) {
    mat <- as.matrix(value_matrix)
    if (ncol(mat) == length(species_names)) {
      colnames(mat) <- species_names
      hc <- tryCatch(stats::hclust(stats::dist(t(mat))), error = function(e) NULL)
      if (!is.null(hc)) return(list(order = hc$labels[hc$order], tree = NULL, note = "estimate similarity"))
    } else if (nrow(mat) == length(species_names)) {
      rownames(mat) <- species_names
      hc <- tryCatch(stats::hclust(stats::dist(mat)), error = function(e) NULL)
      if (!is.null(hc)) return(list(order = hc$labels[hc$order], tree = NULL, note = "estimate similarity"))
    }
  }
  list(order = species_names, tree = NULL, note = "input order")
}

hmsc_extract_mcmc_numeric <- function(x) {
  values <- numeric()
  walk <- function(z) {
    if (inherits(z, "mcmc.list") || inherits(z, "mcmc")) {
      values <<- c(values, suppressWarnings(as.numeric(as.matrix(z))))
    } else if (is.list(z)) {
      for (ii in seq_along(z)) walk(z[[ii]])
    } else if (is.numeric(z)) {
      values <<- c(values, suppressWarnings(as.numeric(z)))
    }
  }
  walk(x)
  values[is.finite(values)]
}

hmsc_prepare_rho_signal_data <- function(m, model_name = "main") {
  draws <- numeric()
  if (requireNamespace("coda", quietly = TRUE)) {
    mpost <- tryCatch(convertToCodaObject(m, spNamesNumbers = c(TRUE, FALSE),
                                          covNamesNumbers = c(TRUE, FALSE)),
                      error = function(e) NULL)
    if (!is.null(mpost) && !is.null(mpost$Rho)) draws <- hmsc_extract_mcmc_numeric(mpost$Rho)
  }
  post <- tryCatch(getPostEstimate(m, parName = "rho"), error = function(e) NULL)
  mean_est <- if (length(draws) > 0) mean(draws, na.rm = TRUE) else suppressWarnings(as.numeric(post$mean %||% NA_real_))[1]
  med_est <- if (length(draws) > 0) stats::median(draws, na.rm = TRUE) else mean_est
  lower <- if (length(draws) > 1) stats::quantile(draws, 0.025, na.rm = TRUE, names = FALSE) else NA_real_
  upper <- if (length(draws) > 1) stats::quantile(draws, 0.975, na.rm = TRUE, names = FALSE) else NA_real_
  support <- suppressWarnings(as.numeric(post$support %||% NA_real_))[1]
  support_neg <- suppressWarnings(as.numeric(post$supportNeg %||% NA_real_))[1]
  list(
    draws = data.frame(model = model_name, draw_id = seq_along(draws), rho = draws, stringsAsFactors = FALSE),
    summary = data.frame(model = model_name, parameter = "rho", mean = mean_est, median = med_est,
                         lower = lower, upper = upper, support_positive = support,
                         support_negative = support_neg, n_draws = length(draws),
                         note = "Hmsc phylogenetic signal rho; near zero means weak phylogenetic signal in residual species responses.",
                         stringsAsFactors = FALSE)
  )
}

hmsc_draw_rho_signal_curve <- function(m, model_name = "main") {
  dat <- hmsc_prepare_rho_signal_data(m, model_name)
  s <- dat$summary
  draws <- dat$draws$rho
  op <- par(mar = c(4.2, 4.6, 3.6, 1.2))
  on.exit(par(op), add = TRUE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  axis(1); axis(2, las = 1); box(col = "#8b8b8b")
  title(main = paste0(model_name, ": phylogenetic signal (rho)"), xlab = "rho", ylab = "relative posterior density")
  if (length(draws) >= 3 && length(unique(round(draws, 8))) >= 3) {
    den <- stats::density(pmin(1, pmax(0, draws)), from = 0, to = 1, na.rm = TRUE)
    y <- den$y / max(den$y, na.rm = TRUE)
    polygon(c(den$x, rev(den$x)), c(rep(0, length(y)), rev(y)),
            col = grDevices::adjustcolor("#7aa6d6", alpha.f = 0.36), border = NA)
    lines(den$x, y, col = "#325f91", lwd = 2)
  } else if (length(draws) > 0) {
    h <- hist(pmin(1, pmax(0, draws)), breaks = seq(0, 1, length.out = 21), plot = FALSE)
    y <- if (max(h$counts) > 0) h$counts / max(h$counts) else h$counts
    rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], y,
         col = grDevices::adjustcolor("#7aa6d6", alpha.f = 0.45), border = "white")
  } else {
    text(0.5, 0.5, "rho posterior draws were not available;\nshowing posterior summary only.", cex = 0.85, col = "#5a5a5a")
  }
  if (is.finite(s$median)) abline(v = s$median, col = "#1b4d59", lwd = 2)
  if (is.finite(s$lower) && is.finite(s$upper)) segments(s$lower, 0.08, s$upper, 0.08, lwd = 4, col = grDevices::adjustcolor("#1b4d59", alpha.f = 0.55))
  legend("topright", bty = "n", cex = 0.82,
         legend = c(paste0("median = ", round(s$median, 3)),
                    paste0("mean = ", round(s$mean, 3)),
                    paste0("95% interval = ", ifelse(is.finite(s$lower), round(s$lower, 3), "NA"),
                           "-", ifelse(is.finite(s$upper), round(s$upper, 3), "NA")),
                    paste0("support positive = ", ifelse(is.finite(s$support_positive), round(s$support_positive, 3), "NA"))))
  invisible(TRUE)
}

hmsc_save_rho_signal_curve <- function(m, model_name = "main", result_dir = "results") {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  dat <- hmsc_prepare_rho_signal_data(m, model_name)
  write.csv(dat$summary, file.path(result_dir, paste0("parameter_estimates_phylogenetic_signal_rho_", model_name, ".csv")), row.names = FALSE)
  if (nrow(dat$draws) > 0) write.csv(dat$draws, file.path(result_dir, paste0("parameter_estimates_phylogenetic_signal_rho_draws_", model_name, ".csv")), row.names = FALSE)
  base <- file.path(result_dir, paste0("parameter_estimates_phylogenetic_signal_rho_", model_name))
  ok <- TRUE
  tryCatch({ grDevices::pdf(paste0(base, ".pdf"), width = 6.2, height = 4.4); hmsc_draw_rho_signal_curve(m, model_name); grDevices::dev.off() },
           error = function(e) { ok <<- FALSE; try(grDevices::dev.off(), silent = TRUE) })
  tryCatch({ grDevices::png(paste0(base, ".png"), width = 1600, height = 1100, res = 240); hmsc_draw_rho_signal_curve(m, model_name); grDevices::dev.off() },
           error = function(e) { ok <<- FALSE; try(grDevices::dev.off(), silent = TRUE) })
  invisible(ok)
}

hmsc_prepare_vp_book_tables <- function(VP, model_name = "main") {
  if (is.null(VP) || is.null(VP$vals)) return(NULL)
  vals <- as.matrix(VP$vals)
  if (nrow(vals) == 0 || ncol(vals) == 0) return(NULL)
  rownames(vals) <- rownames(vals) %||% paste0("predictor_", seq_len(nrow(vals)))
  colnames(vals) <- colnames(vals) %||% paste0("species_", seq_len(ncol(vals)))
  grid <- expand.grid(predictor = rownames(vals), species = colnames(vals), stringsAsFactors = FALSE)
  long <- data.frame(model = model_name, species = grid$species, predictor = grid$predictor,
                     variance_proportion = as.numeric(vals[cbind(match(grid$predictor, rownames(vals)), match(grid$species, colnames(vals)))]),
                     stringsAsFactors = FALSE)
  means <- data.frame(model = model_name, predictor = rownames(vals),
                      mean_variance_proportion = rowMeans(vals, na.rm = TRUE),
                      stringsAsFactors = FALSE)
  means <- means[order(means$mean_variance_proportion, decreasing = TRUE), , drop = FALSE]
  list(vals = vals, long = long, means = means)
}

hmsc_draw_vp_bookstyle <- function(VP, model_name = "main", max_species_labels = 24) {
  dat <- hmsc_prepare_vp_book_tables(VP, model_name)
  if (is.null(dat)) {
    plot.new(); title(main = paste0(model_name, ": variance partitioning")); text(0.5, 0.5, "Variance partitioning was not available.")
    return(invisible(FALSE))
  }
  vals <- dat$vals
  vals[!is.finite(vals)] <- 0
  if (nrow(vals) > 1) vals <- vals[order(rowMeans(vals, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]
  species_order <- order(colSums(vals * seq_len(nrow(vals)), na.rm = TRUE), colnames(vals))
  vals <- vals[, species_order, drop = FALSE]
  pal <- grDevices::hcl.colors(max(3, nrow(vals)), "Spectral", rev = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(2.6, 1.25))
  on.exit(graphics::layout(1), add = TRUE)
  op <- par(oma = c(0, 0, 2.5, 0))
  on.exit(par(op), add = TRUE)
  par(mar = c(6.5, 4.5, 2.2, 0.6))
  bp <- barplot(vals, col = pal, border = NA, space = 0.08, ylim = c(0, 1),
                ylab = "Variance proportion", main = "Variance partitioning", axes = FALSE)
  axis(2, las = 1)
  box(col = "#8b8b8b")
  nsp <- ncol(vals)
  lab_count <- if (nsp <= max_species_labels) nsp else min(14L, max_species_labels)
  lab_idx <- unique(round(seq(1, nsp, length.out = lab_count)))
  lab <- colnames(vals)[lab_idx]
  if (nsp > max_species_labels) lab <- abbreviate(lab, minlength = 8)
  axis(1, at = bp[lab_idx], labels = lab, las = 2, cex.axis = if (nsp > 40) 0.36 else 0.55, tick = FALSE)
  if (nsp > max_species_labels) mtext(paste0(nsp, " species; x-axis labels sampled"), side = 1, line = 5.4, cex = 0.58)
  legend("topright", bty = "n", fill = pal, border = NA, cex = if (nrow(vals) > 7) 0.46 else 0.56,
         legend = paste0(rownames(vals), " (mean = ", round(rowMeans(vals, na.rm = TRUE), 3), ")"))
  par(mar = c(4.6, 9.4, 2.2, 1.0))
  means <- dat$means[order(dat$means$mean_variance_proportion), , drop = FALSE]
  bar_cols <- grDevices::hcl.colors(max(3, nrow(means)), "viridis")
  barplot(means$mean_variance_proportion, horiz = TRUE, names.arg = means$predictor,
          las = 1, col = bar_cols, border = NA, xlab = "Mean proportion", main = "Mean VP per predictor", cex.names = 0.66)
  box(col = "#8b8b8b")
  mtext(paste0(model_name, ": Hmsc variance partitioning by species and predictor"),
        outer = TRUE, side = 3, line = 0.8, cex = 1.0, font = 2)
  invisible(TRUE)
}

hmsc_save_vp_bookstyle <- function(VP, model_name = "main", result_dir = "results") {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  dat <- hmsc_prepare_vp_book_tables(VP, model_name)
  if (is.null(dat)) return(invisible(FALSE))
  write.csv(dat$long, file.path(result_dir, paste0("parameter_estimates_variance_partitioning_long_", model_name, ".csv")), row.names = FALSE)
  write.csv(dat$means, file.path(result_dir, paste0("parameter_estimates_variance_partitioning_mean_by_predictor_", model_name, ".csv")), row.names = FALSE)
  base <- file.path(result_dir, paste0("parameter_estimates_variance_partitioning_bookstyle_", model_name))
  ok <- TRUE
  tryCatch({ grDevices::pdf(paste0(base, ".pdf"), width = 11, height = 7.2); hmsc_draw_vp_bookstyle(VP, model_name); grDevices::dev.off() },
           error = function(e) { ok <<- FALSE; try(grDevices::dev.off(), silent = TRUE) })
  tryCatch({ grDevices::png(paste0(base, ".png"), width = 2600, height = 1700, res = 220); hmsc_draw_vp_bookstyle(VP, model_name); grDevices::dev.off() },
           error = function(e) { ok <<- FALSE; try(grDevices::dev.off(), silent = TRUE) })
  invisible(ok)
}

hmsc_write_association_long <- function(mat, support = NULL, path, model_name = "main", level_name = "random_level") {
  if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 2) return(invisible(FALSE))
  mat <- as.matrix(mat)
  species <- rownames(mat) %||% colnames(mat) %||% paste0("species_", seq_len(nrow(mat)))
  rownames(mat) <- colnames(mat) <- make.unique(as.character(species))
  if (!is.null(support) && identical(dim(support), dim(mat))) {
    support <- as.matrix(support)
    rownames(support) <- colnames(support) <- rownames(mat)
  }
  idx <- which(upper.tri(mat), arr.ind = TRUE)
  out <- data.frame(model = model_name, random_level = level_name,
                    species_i = rownames(mat)[idx[, 1]], species_j = colnames(mat)[idx[, 2]],
                    estimate = as.numeric(mat[idx]),
                    support = if (!is.null(support)) as.numeric(support[idx]) else NA_real_,
                    association_type = "Hmsc residual association (Omega correlation)",
                    scale = "latent residual correlation",
                    comparable = "within_engine_only",
                    note = "Do not directly equate this Omega correlation with associations from other JSDM engines.",
                    stringsAsFactors = FALSE)
  write.csv(out, path, row.names = FALSE)
  invisible(TRUE)
}

hmsc_draw_association_heatmap <- function(m, mat, support = NULL, model_name = "main", level_name = "random_level",
                                          support_level = 0.95, max_species = 80,
                                          show_species_names = TRUE) {
  if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 2) {
    plot.new(); title(main = paste0(model_name, ": species associations")); text(0.5, 0.5, "Association matrix was not available.")
    return(invisible(FALSE))
  }
  mat <- as.matrix(mat)
  species <- rownames(mat) %||% colnames(mat) %||% paste0("species_", seq_len(nrow(mat)))
  rownames(mat) <- colnames(mat) <- make.unique(as.character(species))
  if (!is.null(support) && identical(dim(support), dim(mat))) {
    support <- as.matrix(support)
    rownames(support) <- colnames(support) <- rownames(mat)
  } else {
    support <- NULL
  }
  ord <- hmsc_plot_species_order(m, rownames(mat), mat)$order
  ord <- intersect(ord, rownames(mat))
  if (length(ord) > max_species) ord <- ord[seq_len(max_species)]
  mat <- mat[ord, ord, drop = FALSE]
  if (!is.null(support)) support <- support[ord, ord, drop = FALSE]
  diag(mat) <- 1
  mat[!is.finite(mat)] <- NA_real_
  n <- nrow(mat)
  pal <- grDevices::colorRampPalette(c("#4b306a", "#2b7b7c", "#f0d35b"))(101)
  op <- par(mar = c(if (show_species_names && n <= 80) 8.5 else 4.5,
                    if (show_species_names && n <= 80) 8.5 else 4.5, 4.2, 5.0))
  on.exit(par(op), add = TRUE)
  image(x = seq_len(n), y = seq_len(n), z = t(mat[n:1, , drop = FALSE]),
        axes = FALSE, xlab = "", ylab = "", col = pal, zlim = c(-1, 1), useRaster = TRUE,
        main = paste0(model_name, ": species association matrix"))
  if (show_species_names && n <= 80) {
    axis(1, at = seq_len(n), labels = colnames(mat), las = 2, cex.axis = if (n > 50) 0.38 else 0.55, tick = FALSE)
    axis(2, at = seq_len(n), labels = rev(rownames(mat)), las = 2, cex.axis = if (n > 50) 0.38 else 0.55, tick = FALSE)
  }
  box(col = "#8b8b8b")
  if (!is.null(support)) {
    idx <- which(support >= support_level & row(support) != col(support), arr.ind = TRUE)
    if (nrow(idx) > 0) points(idx[, 2], n - idx[, 1] + 1, pch = 16, cex = 0.18, col = grDevices::adjustcolor("black", alpha.f = 0.65))
  }
  pretty_vals <- pretty(c(-1, 1), n = 5)
  usr <- par("usr")
  x0 <- usr[2] + 0.7
  yseq <- seq(usr[3], usr[4], length.out = length(pal))
  rect(x0, yseq[-length(yseq)], x0 + 0.35, yseq[-1], col = pal[-length(pal)], border = NA, xpd = NA)
  axis(4, at = approx(x = seq(-1, 1, length.out = length(yseq)), y = yseq, xout = pretty_vals)$y,
       labels = pretty_vals, las = 1, cex.axis = 0.72, line = 2.2)
  mtext("Omega correlation", side = 4, line = 4.0, cex = 0.78)
  mtext(paste0("Random level: ", level_name, ". Dots: support >= ", support_level,
               ". Association is Hmsc-specific latent residual correlation."),
        side = 3, line = 0.35, cex = 0.7)
  invisible(TRUE)
}

hmsc_save_association_heatmap <- function(m, mat, support = NULL, model_name = "main", level_name = "random_level", result_dir = "results",
                                          support_level = 0.95, show_species_names = TRUE) {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  safe_level <- gsub("[^A-Za-z0-9_-]+", "_", level_name)
  hmsc_write_association_long(mat, support, file.path(result_dir, paste0("parameter_estimates_species_associations_long_", model_name, "_", safe_level, ".csv")),
                              model_name = model_name, level_name = level_name)
  base <- file.path(result_dir, paste0("parameter_estimates_species_association_matrix_", model_name, "_", safe_level))
  ok <- TRUE
  tryCatch({ grDevices::pdf(paste0(base, ".pdf"), width = 9.2, height = 8.2); hmsc_draw_association_heatmap(m, mat, support, model_name, level_name, support_level, show_species_names = show_species_names); grDevices::dev.off() },
           error = function(e) { ok <<- FALSE; try(grDevices::dev.off(), silent = TRUE) })
  tryCatch({ grDevices::png(paste0(base, ".png"), width = 2100, height = 1900, res = 220); hmsc_draw_association_heatmap(m, mat, support, model_name, level_name, support_level, show_species_names = show_species_names); grDevices::dev.off() },
           error = function(e) { ok <<- FALSE; try(grDevices::dev.off(), silent = TRUE) })
  invisible(ok)
}

hmsc_gradient_axis <- function(Gradient, covariate) {
  xdat <- Gradient$XDataNew %||% NULL
  if (is.null(xdat) || !(covariate %in% colnames(xdat))) {
    return(list(x = seq_len(nrow(xdat %||% data.frame(row = seq_len(1)))), labels = NULL, is_numeric = TRUE))
  }
  raw <- xdat[[covariate]]
  x <- suppressWarnings(as.numeric(as.character(raw)))
  is_num <- all(is.finite(x))
  if (!is_num) {
    labs <- as.character(raw)
    x <- seq_along(labs)
  } else {
    labs <- NULL
  }
  list(x = x, labels = labs, is_numeric = is_num)
}

hmsc_gradient_trait_matrix <- function(TrData, species, trait_list = character()) {
  if (is.null(TrData) || nrow(TrData) == 0 || ncol(TrData) == 0) return(NULL)
  trait_data <- TrData
  if (all(species %in% rownames(trait_data))) {
    trait_data <- trait_data[species, , drop = FALSE]
  } else if (nrow(trait_data) == length(species)) {
    rownames(trait_data) <- species
  } else {
    return(NULL)
  }
  numeric_traits <- vapply(trait_data, is.numeric, logical(1))
  if (length(trait_list) > 0) numeric_traits <- numeric_traits & colnames(trait_data) %in% trait_list
  if (!any(numeric_traits)) return(NULL)
  mat <- as.matrix(trait_data[, numeric_traits, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat
}

hmsc_gradient_response_table <- function(m, Gradient, predY, covariate, model_name = "main",
                                         TrData = NULL, species_list = character(), trait_list = character()) {
  species <- m$spNames %||% colnames(m$Y) %||% paste0("species_", seq_len(ncol(m$Y)))
  axis <- hmsc_gradient_axis(Gradient, covariate)
  ps <- hmsc_summarise_prediction_object(predY, site_ids = paste0("gradient_", seq_along(axis$x)), species_ids = species)
  pm <- ps$mean; pl <- ps$lower; pu <- ps$upper
  if (nrow(pm) != length(axis$x)) return(NULL)
  species_filter <- intersect(species_list, colnames(pm))
  focal <- if (length(species_filter) > 0) species_filter[1] else {
    ssd <- suppressWarnings(apply(pm, 2, stats::sd, na.rm = TRUE))
    ssd[!is.finite(ssd)] <- -Inf
    colnames(pm)[which.max(ssd)]
  }
  rows <- list()
  rows[[length(rows) + 1]] <- data.frame(model = model_name, covariate = covariate, curve_type = "summed_response",
                                         response = "all_species", trait = NA_character_, x = axis$x,
                                         x_label = if (is.null(axis$labels)) as.character(round(axis$x, 6)) else axis$labels,
                                         mean = rowSums(pm, na.rm = TRUE),
                                         lower = rowSums(pl, na.rm = TRUE),
                                         upper = rowSums(pu, na.rm = TRUE),
                                         stringsAsFactors = FALSE)
  if (!is.null(focal) && length(focal) == 1 && nzchar(focal)) {
    rows[[length(rows) + 1]] <- data.frame(model = model_name, covariate = covariate, curve_type = "focal_species",
                                           response = focal, trait = NA_character_, x = axis$x,
                                           x_label = if (is.null(axis$labels)) as.character(round(axis$x, 6)) else axis$labels,
                                           mean = as.numeric(pm[, focal]), lower = as.numeric(pl[, focal]), upper = as.numeric(pu[, focal]),
                                           stringsAsFactors = FALSE)
  }
  trait_mat <- hmsc_gradient_trait_matrix(TrData, colnames(pm), trait_list)
  if (!is.null(trait_mat) && ncol(trait_mat) > 0) {
    denom <- pmax(rowSums(pm, na.rm = TRUE), .Machine$double.eps)
    cwm <- sweep(as.matrix(pm) %*% trait_mat, 1, denom, "/")
    trait_nm <- colnames(cwm)[1]
    rows[[length(rows) + 1]] <- data.frame(model = model_name, covariate = covariate, curve_type = "community_weighted_trait",
                                           response = "predicted_community", trait = trait_nm, x = axis$x,
                                           x_label = if (is.null(axis$labels)) as.character(round(axis$x, 6)) else axis$labels,
                                           mean = as.numeric(cwm[, trait_nm]), lower = NA_real_, upper = NA_real_,
                                           stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

hmsc_panel_curve <- function(tab, title, ylab, ribbon = TRUE) {
  if (is.null(tab) || nrow(tab) == 0) {
    plot.new(); title(main = title); text(0.5, 0.5, "No curve data available.")
    return(invisible(FALSE))
  }
  x <- suppressWarnings(as.numeric(tab$x))
  y <- suppressWarnings(as.numeric(tab$mean))
  lo <- suppressWarnings(as.numeric(tab$lower))
  hi <- suppressWarnings(as.numeric(tab$upper))
  finite <- is.finite(x) & is.finite(y)
  if (!any(finite)) {
    plot.new(); title(main = title); text(0.5, 0.5, "No finite curve values.")
    return(invisible(FALSE))
  }
  ord <- order(x)
  x <- x[ord]; y <- y[ord]; lo <- lo[ord]; hi <- hi[ord]
  ylim <- range(c(y, lo, hi), na.rm = TRUE)
  plot(x, y, type = "n", xlab = unique(tab$covariate)[1], ylab = ylab, main = title, ylim = ylim)
  if (isTRUE(ribbon) && all(is.finite(lo)) && all(is.finite(hi))) {
    polygon(c(x, rev(x)), c(lo, rev(hi)), col = grDevices::adjustcolor("#6d6aa8", alpha.f = 0.24), border = NA)
  }
  lines(x, y, col = "#25245f", lwd = 2.2)
  points(x, y, pch = 16, cex = 0.65, col = grDevices::adjustcolor("#25245f", alpha.f = 0.7))
  invisible(TRUE)
}

hmsc_draw_gradient_response_panel <- function(m, Gradient, predY, covariate, model_name = "main",
                                              TrData = NULL, species_list = character(), trait_list = character()) {
  tab <- hmsc_gradient_response_table(m, Gradient, predY, covariate, model_name, TrData, species_list, trait_list)
  if (is.null(tab) || nrow(tab) == 0) {
    plot.new(); title(main = paste0(model_name, ": gradient response")); text(0.5, 0.5, "Gradient predictions were not available.")
    return(invisible(FALSE))
  }
  op <- par(mfrow = c(2, 2), oma = c(0, 0, 2.8, 0), mar = c(4.2, 4.5, 2.6, 1.0))
  on.exit(par(op), add = TRUE)
  hmsc_panel_curve(tab[tab$curve_type == "summed_response", , drop = FALSE],
                   "Summed community response", "summed predicted response")
  focal_tab <- tab[tab$curve_type == "focal_species", , drop = FALSE]
  hmsc_panel_curve(focal_tab, if (nrow(focal_tab) > 0) unique(focal_tab$response)[1] else "Focal species response",
                   "predicted response")
  trait_tab <- tab[tab$curve_type == "community_weighted_trait", , drop = FALSE]
  hmsc_panel_curve(trait_tab, if (nrow(trait_tab) > 0) paste0("CWM trait: ", unique(trait_tab$trait)[1]) else "Community-weighted trait",
                   "community-weighted mean", ribbon = FALSE)
  plot.new()
  text(0, 0.85, "Interpretation", adj = 0, font = 2)
  text(0, 0.68, paste0("Covariate: ", covariate), adj = 0)
  text(0, 0.52, "Lines show posterior mean predictions.", adj = 0)
  text(0, 0.36, "Bands show approximate 95% pointwise intervals when available.", adj = 0)
  text(0, 0.20, "These are quick-test diagnostic curves unless MCMC settings are publication scale.", adj = 0)
  mtext(paste0(model_name, ": Hmsc gradient response curves"), outer = TRUE, side = 3, line = 1.0, cex = 1.0, font = 2)
  invisible(TRUE)
}

hmsc_save_gradient_response_panel <- function(m, Gradient, predY, covariate, model_name = "main", result_dir = "predictions",
                                              TrData = NULL, species_list = character(), trait_list = character()) {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  safe_cov <- gsub("[^A-Za-z0-9_-]+", "_", covariate)
  tab <- hmsc_gradient_response_table(m, Gradient, predY, covariate, model_name, TrData, species_list, trait_list)
  if (is.null(tab) || nrow(tab) == 0) return(invisible(FALSE))
  write.csv(tab, file.path(result_dir, paste0("gradient_response_curves_", model_name, "_", safe_cov, ".csv")), row.names = FALSE)
  base <- file.path(result_dir, paste0("gradient_response_curves_", model_name, "_", safe_cov))
  ok <- TRUE
  tryCatch({ grDevices::pdf(paste0(base, ".pdf"), width = 9.6, height = 7.2); hmsc_draw_gradient_response_panel(m, Gradient, predY, covariate, model_name, TrData, species_list, trait_list); grDevices::dev.off() },
           error = function(e) { ok <<- FALSE; try(grDevices::dev.off(), silent = TRUE) })
  tryCatch({ grDevices::png(paste0(base, ".png"), width = 2200, height = 1650, res = 220); hmsc_draw_gradient_response_panel(m, Gradient, predY, covariate, model_name, TrData, species_list, trait_list); grDevices::dev.off() },
           error = function(e) { ok <<- FALSE; try(grDevices::dev.off(), silent = TRUE) })
  invisible(ok)
}

numeric_matrix_check <- function(dat, label = "Y") {
  if (is.null(dat)) return(list(matrix = NULL, ok = FALSE, messages = paste0(label, " is missing.")))
  raw <- as.matrix(dat)
  numeric_values <- suppressWarnings(as.numeric(raw))
  mat <- matrix(numeric_values, nrow = nrow(raw), ncol = ncol(raw), dimnames = dimnames(raw))
  raw_chr <- trimws(as.character(raw))
  bad <- !is.na(raw_chr) & nzchar(raw_chr) & is.na(numeric_values)
  messages <- character()
  if (any(bad)) {
    messages <- c(messages, paste0(label, " contains non-numeric response values that cannot be converted."))
  }
  list(matrix = mat, ok = !any(bad), messages = messages)
}

response_family_messages <- function(Ymat, family, engine = "engine") {
  ok <- TRUE
  msg <- character()
  if (is.null(Ymat)) return(list(ok = FALSE, messages = paste0(engine, " response matrix is missing.")))
  fam <- tolower(family %||% "")
  is_binary_family <- fam %in% c("probit", "binomial", "binomial_probit", "binomial_logit", "binomial_probit_sp_constrained", "pa")
  is_count_family <- fam %in% c("poisson", "poisson_log", "nbinom", "negative.binomial", "negative_binomial", "negative_binomial_log", "ztpoisson", "ztnegative.binomial", "da", "cc")
  is_gaussian_family <- fam %in% c("normal", "gaussian", "gaussian_identity", "con")
  is_beta_family <- fam %in% c("beta", "fc")
  is_ordinal_family <- fam %in% c("ordinal", "oc", "cat")
  is_positive_family <- fam %in% c("lnormal", "lognormal", "gamma", "exponential")
  if (any(!is.finite(Ymat), na.rm = TRUE)) msg <- c(msg, paste0(engine, " response matrix contains non-finite values."))
  if (is_binary_family && any(!(Ymat %in% c(0, 1, NA)))) {
    ok <- FALSE
    msg <- c(msg, paste0(engine, " binary/probit/binomial responses must be 0/1/NA."))
  }
  if (is_count_family) {
    if (any(Ymat < 0, na.rm = TRUE)) {
      ok <- FALSE
      msg <- c(msg, paste0(engine, " count responses must be non-negative."))
    }
    if (any(abs(Ymat - round(Ymat)) > 1e-8, na.rm = TRUE)) {
      ok <- FALSE
      msg <- c(msg, paste0(engine, " count responses must be integer-valued."))
    }
    if (fam %in% c("ztpoisson", "ztnegative.binomial") && any(Ymat < 1, na.rm = TRUE)) {
      ok <- FALSE
      msg <- c(msg, paste0(engine, " zero-truncated count responses must be positive integers."))
    }
  }
  if (is_positive_family && any(Ymat <= 0, na.rm = TRUE)) {
    ok <- FALSE
    msg <- c(msg, paste0(engine, " ", fam, " responses must be positive."))
  }
  if (is_beta_family && any(Ymat <= 0 | Ymat >= 1, na.rm = TRUE)) {
    ok <- FALSE
    msg <- c(msg, paste0(engine, " beta/fraction responses must be strictly between 0 and 1."))
  }
  if (is_ordinal_family && any(Ymat < 1 | abs(Ymat - round(Ymat)) > 1e-8, na.rm = TRUE)) {
    ok <- FALSE
    msg <- c(msg, paste0(engine, " ordinal responses must be positive integer levels starting at 1."))
  }
  if (is_gaussian_family && any(is.na(Ymat))) {
    msg <- c(msg, paste0(engine, " gaussian/normal responses contain NA; check whether the engine can handle missing responses."))
  }
  list(ok = ok, messages = msg)
}

validate_one_sided_formula <- function(formula_text, data, label = "formula") {
  ok <- TRUE
  msg <- character()
  formula_text <- paste(as.character(formula_text %||% ""), collapse = " ")
  if (is.null(data) || !nzchar(trimws(formula_text))) return(list(ok = ok, messages = msg))
  f <- tryCatch(as.formula(formula_text), error = function(e) e)
  if (inherits(f, "error")) {
    return(list(ok = FALSE, messages = paste0(label, " is not a valid R formula: ", conditionMessage(f))))
  }
  vars <- setdiff(all.vars(f), ".")
  missing <- setdiff(vars, colnames(data))
  if (length(missing) > 0) {
    ok <- FALSE
    msg <- c(msg, paste0(label, " variables missing from data: ", paste(missing, collapse = ", ")))
  }
  list(ok = ok, messages = msg)
}

run_hmsc_s1s7_pipeline <- function(outdir, cfg, Y, XData, TrData = NULL, studyDesign = NULL, coord = NULL,
                                   newdata = NULL, newStudyDesign = NULL, newcoord = NULL,
                                   log_fun = message) {
  start_time <- Sys.time()
  result <- list(status = "not_started", warnings = character(), errors = character(), runtime_seconds = 0)
  safe_log <- function(...) {
    txt <- paste(...)
    try(log_fun(txt), silent = TRUE)
  }
  fail_now <- function(msg) {
    result$status <<- "fit_failed"
    result$errors <<- c(result$errors, msg)
    safe_log("HMSC S1-S7 error:", msg)
    try(writeLines(msg, file.path(outdir, "diagnostics", "HMSC_S1S7_error.txt")), silent = TRUE)
    invisible(NULL)
  }

  for (d in c("data", "models", "results", "plots", "tables", "predictions", "diagnostics", "workflow_scripts", "reproducible_script", "standard", "report")) {
    dir.create(file.path(outdir, d), showWarnings = FALSE, recursive = TRUE)
  }

  parameter_audit <- tryCatch(do.call(rbind, hmsc_flatten_config(cfg)), error = function(e) NULL)
  if (!is.null(parameter_audit)) {
    parameter_audit$status <- "read_by_server"
    parameter_audit$status[grepl("^mcmc.updater.(Beta|Gamma|Omega)$", parameter_audit$parameter)] <- "documented_Hmsc_sampleMcmc_has_no_direct_argument"
    try(write.csv(parameter_audit, file.path(outdir, "tables", "HMSC_parameter_audit.csv"), row.names = FALSE), silent = TRUE)
  }

  # Write input objects to data/ so the exported S1-S7 scripts are truly reproducible.
  try(write.csv(Y, file.path(outdir, "data", "Y.csv"), row.names = TRUE), silent = TRUE)
  try(write.csv(XData, file.path(outdir, "data", "XData.csv"), row.names = TRUE), silent = TRUE)
  if (!is.null(TrData)) try(write.csv(TrData, file.path(outdir, "data", "TrData.csv"), row.names = TRUE), silent = TRUE)
  if (!is.null(studyDesign)) try(write.csv(studyDesign, file.path(outdir, "data", "studyDesign.csv"), row.names = TRUE), silent = TRUE)
  if (!is.null(coord)) try(write.csv(coord, file.path(outdir, "data", "coordinates.csv"), row.names = TRUE), silent = TRUE)
  if (!is.null(newdata)) try(write.csv(newdata, file.path(outdir, "data", "newdata.csv"), row.names = TRUE), silent = TRUE)
  if (!is.null(newStudyDesign)) try(write.csv(newStudyDesign, file.path(outdir, "data", "prediction_studyDesign.csv"), row.names = TRUE), silent = TRUE)
  if (!is.null(newcoord)) try(write.csv(newcoord, file.path(outdir, "data", "newcoords.csv"), row.names = TRUE), silent = TRUE)

  if (!requireNamespace("Hmsc", quietly = TRUE)) {
    fail_now("Hmsc package is not installed. Install Hmsc before running the real S1-S7 workflow.")
    result$runtime_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    writeLines(jsonlite::toJSON(result, pretty = TRUE, auto_unbox = TRUE), file.path(outdir, "diagnostics", "HMSC_S1S7_pipeline_status.json"))
    writeLines(jsonlite::toJSON(list(engine = "Hmsc", status = result$status, warnings = result$warnings, errors = result$errors, runtime_seconds = result$runtime_seconds), pretty = TRUE, auto_unbox = TRUE),
               file.path(outdir, "diagnostics", "engine_status.json"))
    writeLines(capture.output(sessionInfo()), file.path(outdir, "diagnostics", "session_info.txt"))
    return(result)
  }

  # Optional packages used by plotting/result scripts.
  for (pkg in c("coda", "corrplot", "colorspace", "writexl", "vioplot", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      result$warnings <- c(result$warnings, paste0("Optional package not installed: ", pkg))
    }
  }

  # Generate executable S1-S7 scripts in the output folder. These scripts use relative paths only.
  generate_hmsc_s1s7_scripts <- function(outdir, cfg) {
    random_helpers <- c("%||%", "hmsc_factor_preserve_order", "hmsc_clean_dataframe_types",
      "hmsc_normalize_spatial_method", "hmsc_normalize_random_mode", "hmsc_unit_ids",
      "hmsc_prepare_study_design", "hmsc_select_coordinate_columns", "hmsc_prepare_spatial_matrix",
      "hmsc_make_gpp_knots", "hmsc_optional_number", "hmsc_parse_alphapw",
      "hmsc_apply_random_level_priors", "hmsc_find_workflow_file", "hmsc_read_optional_dataframe",
      "hmsc_read_square_numeric_matrix", "hmsc_align_random_level_matrix", "hmsc_align_random_level_xdata",
      "hmsc_read_design_csv", "hmsc_align_rows", "hmsc_group_factor", "hmsc_build_random_level",
      "hmsc_build_primary_random_level", "hmsc_prepare_prediction_study_design",
      "hmsc_prepare_prediction_ranlevels", "hmsc_predict_new_sites")
    dump(random_helpers, file = file.path(outdir, "workflow_scripts", "HMSC_random_level_helpers.R"),
         envir = environment(hmsc_build_random_level))
    samples <- as.integer(cfg$mcmc$samples %||% 50)
    thin <- as.integer(cfg$mcmc$thin %||% 1)
    transient <- as.integer(cfg$mcmc$transient %||% 25)
    nChains <- as.integer(cfg$mcmc$nChains %||% 2)
    nParallel <- as.integer(cfg$mcmc$nParallel %||% 1)
    nfolds <- as.integer(cfg$outputs$predictions$nfolds %||% cfg$outputs$nfolds %||% 2)
    distr <- cfg$model$distr %||% "probit"
    XFormula <- cfg$model$XFormula %||% "~ ."
    TrFormula <- cfg$model$TrFormula %||% "~ ."
    mode_info <- hmsc_normalize_random_mode(cfg$model$random_mode %||% "none", cfg$model$spatial_method %||% "NNGP")
    lon_col <- cfg$model$lon_col %||% "longitude"
    lat_col <- cfg$model$lat_col %||% "latitude"
    spatial_method <- mode_info$spatial_method
    use_traits <- isTRUE(cfg$model$use_traits)
    use_phylogeny <- isTRUE(cfg$model$use_phylogeny)
    phylogeny_file <- cfg$data$phylogeny %||% ""
    C_file <- cfg$model$C_file %||% ""
    initPar_text <- trimws(as.character(cfg$mcmc$initPar %||% "fixed effects"))
    alignPost <- isTRUE(cfg$mcmc$alignPost %||% TRUE)
    fromPrior <- isTRUE(cfg$mcmc$sample_prior %||% FALSE)
    pool_chains <- isTRUE(cfg$mcmc$pool_chains %||% FALSE)
    updater_GammaEta <- isTRUE(cfg$mcmc$updater$GammaEta %||% TRUE)
    compute_waic <- isTRUE(cfg$outputs$waic %||% FALSE)
    XScale <- isTRUE(cfg$model$XScale %||% TRUE)
    TrScale <- isTRUE(cfg$model$TrScale %||% TRUE)
    YScale <- isTRUE(cfg$model$YScale %||% FALSE)
    truncateNumberOfFactors <- isTRUE(cfg$model$truncateNumberOfFactors %||% TRUE)
    Loff_file <- cfg$model$Loff_file %||% ""
    use_XRRR <- isTRUE(cfg$model$use_XRRR %||% FALSE)
    XRRR_file <- cfg$model$XRRR_file %||% ""
    XRRRFormula <- cfg$model$XRRRFormula %||% "~ ."
    XRRRScale <- isTRUE(cfg$model$XRRRScale %||% TRUE)
    ncRRR <- as.integer(cfg$model$ncRRR %||% 2)
    advanced_random_type <- cfg$model$random_level_type %||% "none"
    advanced_random_key <- gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(trimws(as.character(advanced_random_type)))))
    if (advanced_random_key == "covariate_dependent_xdata") updater_GammaEta <- FALSE
    if (spatial_method %in% c("GPP", "NNGP")) updater_GammaEta <- FALSE
    if (identical(spatial_method, "GPP") && isTRUE(alignPost)) {
      alignPost <- FALSE
      result$scientific_caveats <- c(result$scientific_caveats %||% character(), "Hmsc alignPost was disabled for GPP spatial quick fitting because alignPosterior can return NA with tiny smoke-test posterior samples.")
    }
    partition_column <- cfg$outputs$predictions$partition_column %||% ""
    env_list <- paste(hmsc_parse_name_list(cfg$outputs$predictions$env.list %||% ""), collapse = ",")
    species_list <- paste(hmsc_parse_name_list(cfg$outputs$predictions$species.list %||% ""), collapse = ",")
    trait_list <- paste(hmsc_parse_name_list(cfg$outputs$predictions$trait.list %||% ""), collapse = ",")
    computeSAIR <- isTRUE(cfg$outputs$predictions$computeSAIR %||% FALSE)

    s1 <- c(
      "# S1_define_models.R -- generated by JSDMWorkbench",
      paste0("set.seed(", as.integer(cfg$model$seed %||% 1234L), ")"),
      "library(Hmsc)",
      "localDir <- '.'",
      "dataDir <- file.path(localDir, 'data')",
      "modelDir <- file.path(localDir, 'models')",
      "if (!dir.exists(modelDir)) dir.create(modelDir, recursive = TRUE)",
      "Y <- read.csv(file.path(dataDir, 'Y.csv'), row.names = 1, check.names = FALSE)",
      "XData <- read.csv(file.path(dataDir, 'XData.csv'), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)",
      "factor_preserve_order <- function(x) { x <- as.character(x); x[trimws(x) == ''] <- NA; factor(x, levels = unique(x[!is.na(x)])) }",
      "clean_types <- function(df) { for (nm in names(df)) { if (is.character(df[[nm]])) { x <- trimws(df[[nm]]); x[x == ''] <- NA; nx <- suppressWarnings(as.numeric(x)); df[[nm]] <- if (all(is.na(x) | !is.na(nx))) nx else factor_preserve_order(x) }; if (is.logical(df[[nm]])) df[[nm]] <- factor_preserve_order(df[[nm]]) }; df }",
      "select_coord_cols <- function(coord, lon_col, lat_col) { nm <- names(coord); find_col <- function(x) { hit <- nm[tolower(nm) == tolower(x)]; if (length(hit) > 0) hit[1] else NA_character_ }; req <- c(find_col(lon_col), find_col(lat_col)); if (all(!is.na(req))) return(req); for (pair in list(c('longitude','latitude'), c('lon','lat'), c('x','y'), c('easting','northing'))) { hit <- c(find_col(pair[1]), find_col(pair[2])); if (all(!is.na(hit))) return(hit) }; num <- nm[vapply(coord, function(z) all(is.na(z) | !is.na(suppressWarnings(as.numeric(z)))), logical(1))]; if (length(num) >= 2) return(num[1:2]); stop('Coordinate columns not found in coordinates.csv.') }",
      "find_workflow_file <- function(file_name) { file_name <- trimws(as.character(file_name)); if (!nzchar(file_name)) return(NULL); candidates <- c(file_name, file.path('data', file_name), file.path('inputs', file_name)); candidates <- candidates[file.exists(candidates)]; if (length(candidates) == 0) NULL else candidates[1] }",
      "read_numeric_matrix <- function(file_name, expected_rows=NULL, expected_cols=NULL, label='matrix') { path <- find_workflow_file(file_name); if (is.null(path)) return(NULL); mat <- as.matrix(read.csv(path, row.names=1, check.names=FALSE)); storage.mode(mat) <- 'numeric'; if (any(!is.finite(mat))) stop(paste0(label, ' must be numeric and finite.')); if (!is.null(expected_rows) && nrow(mat) != expected_rows) stop(paste0(label, ' row count mismatch.')); if (!is.null(expected_cols) && ncol(mat) != expected_cols) stop(paste0(label, ' column count mismatch.')); mat }",
      "read_optional_dataframe <- function(file_name, expected_rows=NULL, label='data frame') { path <- find_workflow_file(file_name); if (is.null(path)) return(NULL); dat <- read.csv(path, row.names=1, check.names=FALSE, stringsAsFactors=FALSE); if (!is.null(expected_rows) && nrow(dat) != expected_rows) stop(paste0(label, ' row count mismatch.')); clean_types(dat) }",
      "make_pd_corr <- function(C, eps=1e-6) { C <- as.matrix(C); storage.mode(C) <- 'numeric'; C[!is.finite(C)] <- 0; C <- (C + t(C))/2; diag(C) <- 1; ev <- eigen(C, symmetric=TRUE); if (min(ev$values, na.rm=TRUE) <= eps) { C <- ev$vectors %*% diag(pmax(ev$values, eps), nrow=length(ev$values)) %*% t(ev$vectors); C <- (C + t(C))/2; d <- sqrt(pmax(diag(C), eps)); C <- C / outer(d, d); diag(C) <- 1 }; C }",
      "taxonomy_to_C <- function(taxonomy, species_names) { taxonomy <- as.data.frame(taxonomy, stringsAsFactors=FALSE, check.names=FALSE); nm_low <- tolower(names(taxonomy)); species_col <- names(taxonomy)[nm_low %in% c('species','species_name','species_id','taxon','taxon_id','taxon_name','response','response_id')]; rownames_match <- !is.null(rownames(taxonomy)) && length(setdiff(species_names, rownames(taxonomy))) == 0; if (!rownames_match && length(species_col) > 0) rownames(taxonomy) <- as.character(taxonomy[[species_col[1]]]); missing_species <- setdiff(species_names, rownames(taxonomy)); if (length(missing_species) > 0) stop(paste0('Taxonomy file is missing species: ', paste(missing_species, collapse=', '))); taxonomy <- taxonomy[species_names,,drop=FALSE]; rank_order <- c('kingdom','phylum','division','class','order','family','genus','species'); rank_cols <- names(taxonomy)[tolower(names(taxonomy)) %in% rank_order]; rank_cols <- rank_cols[order(match(tolower(rank_cols), rank_order))]; if (length(species_col) > 0) rank_cols <- setdiff(rank_cols, species_col[1]); if (length(rank_cols) == 0) stop('Taxonomy file must contain taxonomic-rank columns.'); ranks <- taxonomy[,rank_cols,drop=FALSE]; ranks[] <- lapply(ranks, function(z) { z <- trimws(as.character(z)); z[z==''] <- NA_character_; z }); n <- length(species_names); C <- matrix(0, n, n, dimnames=list(species_names, species_names)); for (i in seq_len(n)) for (j in seq_len(n)) { comparable <- !is.na(unlist(ranks[i,,drop=TRUE])) & !is.na(unlist(ranks[j,,drop=TRUE])); C[i,j] <- if (!any(comparable)) as.numeric(i == j) else mean(unlist(ranks[i,comparable,drop=TRUE]) == unlist(ranks[j,comparable,drop=TRUE])) }; make_pd_corr(C) }",
      "read_optional_matrix <- function(file_name, species_names) { mat <- read_numeric_matrix(file_name, label='C matrix'); if (is.null(mat)) return(NULL); if (all(species_names %in% rownames(mat)) && all(species_names %in% colnames(mat))) mat <- mat[species_names, species_names, drop=FALSE] else { rownames(mat) <- species_names; colnames(mat) <- species_names }; make_pd_corr(mat) }",
      "read_phylogeny_or_taxonomy <- function(file_name, species_names) { path <- find_workflow_file(file_name); if (is.null(path)) return(NULL); ext <- tolower(tools::file_ext(path)); if (ext %in% c('csv','tsv')) { taxonomy <- if (ext == 'tsv') read.delim(path, row.names=NULL, check.names=FALSE, stringsAsFactors=FALSE) else read.csv(path, row.names=NULL, check.names=FALSE, stringsAsFactors=FALSE); return(list(C=taxonomy_to_C(taxonomy, species_names), source='taxonomy_table', path=path)) }; if (!requireNamespace('ape', quietly=TRUE)) stop('ape is required to read Newick phylogeny files.'); tree <- ape::read.tree(path); missing_species <- setdiff(species_names, tree$tip.label); if (length(missing_species) > 0) stop(paste0('Phylogeny tree is missing species: ', paste(missing_species, collapse=', '))); extra <- setdiff(tree$tip.label, species_names); if (length(extra) > 0) tree <- ape::keep.tip(tree, species_names); list(phyloTree=tree, source='newick_tree', path=path) }",
      "write_derived_phylogeny <- function(phy=NULL, species_names, Cmat=NULL, source=NA_character_, input_file=NA_character_) { derived <- NULL; if (!is.null(Cmat)) derived <- make_pd_corr(Cmat) else if (!is.null(phy) && !is.null(phy$C)) derived <- make_pd_corr(phy$C) else if (!is.null(phy) && !is.null(phy$phyloTree) && requireNamespace('ape', quietly=TRUE)) derived <- tryCatch({ tree <- phy$phyloTree; if (length(setdiff(tree$tip.label, species_names)) > 0) tree <- ape::keep.tip(tree, species_names); C <- ape::vcv.phylo(tree, corr=TRUE); make_pd_corr(C[species_names, species_names, drop=FALSE]) }, error=function(e) NULL); if (is.null(derived)) return(invisible(FALSE)); rownames(derived) <- species_names; colnames(derived) <- species_names; for (d in c('data','tables','standard')) { dir.create(d, showWarnings=FALSE, recursive=TRUE); write.csv(derived, file.path(d, 'derived_phylogeny_correlation.csv')) }; write.csv(data.frame(source=source, input_file=input_file, derived_file='derived_phylogeny_correlation.csv', user_upload_contract='Upload a Newick tree or taxonomy rank table; do not upload a species-by-species covariance/correlation matrix.', stringsAsFactors=FALSE), file.path('tables','phylogeny_input_contract.csv'), row.names=FALSE); invisible(TRUE) }",
      "XData <- clean_types(XData)",
      "Y <- as.matrix(Y); storage.mode(Y) <- 'numeric'",
      paste0("XFormula <- as.formula('", gsub("'", "\\\\'", XFormula), "')"),
      paste0("distr <- '", distr, "'"),
      paste0("XScale <- ", if (XScale) "TRUE" else "FALSE"),
      paste0("YScale <- ", if (YScale) "TRUE" else "FALSE"),
      paste0("truncateNumberOfFactors <- ", if (truncateNumberOfFactors) "TRUE" else "FALSE"),
      paste0("Loff_file <- '", gsub("'", "\\\\'", Loff_file), "'"),
      "model_args <- list(Y = Y, XData = XData, XFormula = XFormula, distr = distr, XScale = XScale, YScale = YScale, truncateNumberOfFactors = truncateNumberOfFactors)",
      "Loff <- read_numeric_matrix(Loff_file, expected_rows=nrow(Y), expected_cols=ncol(Y), label='Loff offset'); if (!is.null(Loff)) model_args$Loff <- Loff",
      if (use_traits) paste0("TrData <- read.csv(file.path(dataDir, 'TrData.csv'), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)\nTrData <- clean_types(TrData)\nmodel_args$TrData <- TrData\nmodel_args$TrFormula <- as.formula('", gsub("'", "\\\\'", TrFormula), "')\nmodel_args$TrScale <- ", if (TrScale) "TRUE" else "FALSE") else "# traits not used",
      paste0("use_XRRR <- ", if (use_XRRR) "TRUE" else "FALSE"),
      paste0("XRRR_file <- '", gsub("'", "\\\\'", XRRR_file), "'"),
      paste0("XRRRFormula <- as.formula('", gsub("'", "\\\\'", XRRRFormula), "')"),
      paste0("XRRRScale <- ", if (XRRRScale) "TRUE" else "FALSE"),
      paste0("ncRRR <- ", ncRRR),
      "if (use_XRRR) { XRRRData <- read_optional_dataframe(XRRR_file, expected_rows=nrow(Y), label='XRRRData'); if (is.null(XRRRData)) stop('Use XRRR is enabled, but XRRRData file was not found.'); model_args$XRRRData <- XRRRData; model_args$XRRRFormula <- XRRRFormula; model_args$XRRRScale <- XRRRScale; model_args$ncRRR <- ncRRR }",
      paste0("use_phylogeny <- ", if (use_phylogeny) "TRUE" else "FALSE"),
      paste0("phylogeny_file <- '", gsub("'", "\\\\'", phylogeny_file), "'"),
      paste0("C_file <- '", gsub("'", "\\\\'", C_file), "'"),
      "if (use_phylogeny) {",
      "  Cmat <- read_optional_matrix(C_file, colnames(Y))",
      "  if (!is.null(Cmat)) {",
      "    model_args$C <- Cmat",
      "    write.csv(data.frame(enabled=TRUE, source='C_matrix', file=C_file, mode='C'), file.path('tables', 'S1_phylogeny_summary.csv'), row.names=FALSE)",
      "    write_derived_phylogeny(species_names=colnames(Y), Cmat=Cmat, source='deprecated_internal_C_matrix', input_file=C_file)",
      "  } else {",
      "    phy <- read_phylogeny_or_taxonomy(phylogeny_file, colnames(Y))",
      "    if (is.null(phy)) stop('Use phylogeny/taxonomy is enabled, but no Newick tree or taxonomy rank table was found.')",
      "    if (!is.null(phy$C)) { model_args$C <- phy$C; mode <- 'C' } else { model_args$phyloTree <- phy$phyloTree; mode <- 'phyloTree' }",
      "    write.csv(data.frame(enabled=TRUE, source=phy$source, file=phy$path, mode=mode), file.path('tables', 'S1_phylogeny_summary.csv'), row.names=FALSE)",
      "    write_derived_phylogeny(phy=phy, species_names=colnames(Y), source=phy$source, input_file=phy$path)",
      "  }",
      "} else {",
      "  write.csv(data.frame(enabled=FALSE, source=NA_character_, file=NA_character_, mode=NA_character_), file.path('tables', 'S1_phylogeny_summary.csv'), row.names=FALSE)",
      "}",
      "source(file.path('workflow_scripts', 'HMSC_random_level_helpers.R'))",
      paste0("random_cfg <- list(model = ", paste(capture.output(dput(cfg$model)), collapse = "\n"), ")"),
      "study <- hmsc_read_design_csv(file.path(dataDir, 'studyDesign.csv'))",
      "coord <- hmsc_read_design_csv(file.path(dataDir, 'coordinates.csv'))",
      "random_level <- hmsc_build_random_level(Y, study, coord, random_cfg)",
      "if (length(random_level$ranLevels)) { model_args$studyDesign <- random_level$studyDesign; model_args$ranLevels <- random_level$ranLevels }",
      "write.csv(random_level$summary, file.path('tables', 'S1_random_level_summary.csv'), row.names = FALSE)",
      "m <- do.call(Hmsc, model_args)",
      "models <- list(main = m)",
      "save(models, file = file.path(modelDir, 'unfitted_models.RData'))"
    )
    writeLines(s1, file.path(outdir, "workflow_scripts", "S1_define_models.R"))

    s2 <- c(
      "# S2_fit_models.R -- generated by JSDMWorkbench",
      paste0("set.seed(", as.integer(cfg$model$seed %||% 1234L), ")"),
      "library(Hmsc)",
      "modelDir <- file.path('.', 'models')",
      "load(file.path(modelDir, 'unfitted_models.RData'))",
      paste0("samples <- ", samples),
      paste0("thin <- ", thin),
      paste0("transient <- ", transient),
      paste0("nChains <- ", nChains),
      paste0("nParallel <- ", nParallel),
      paste0("verbose <- ", as.integer(cfg$mcmc$verbose %||% 10)),
      paste0("initPar <- ", if (tolower(initPar_text) == "fixed effects") "'fixed effects'" else "NULL"),
      paste0("alignPost <- ", if (alignPost) "TRUE" else "FALSE"),
      paste0("fromPrior <- ", if (fromPrior) "TRUE" else "FALSE"),
      paste0("pool_chains <- ", if (pool_chains) "TRUE" else "FALSE"),
      paste0("updater <- list(GammaEta = ", if (updater_GammaEta) "TRUE" else "FALSE", ", Gamma2 = FALSE)"),
      "for (mi in seq_along(models)) {",
      "  models[[mi]] <- sampleMcmc(models[[mi]], samples = samples, thin = thin, transient = transient, initPar = initPar, verbose = verbose, adaptNf = rep(max(1L, ceiling(0.4 * samples * thin)), models[[mi]]$nr), nChains = nChains, nParallel = nParallel, updater = updater, fromPrior = fromPrior, alignPost = alignPost, engine = 'R')",
      "  if (pool_chains && exists('poolMcmcChains', mode='function')) { pooled <- tryCatch(poolMcmcChains(models[[mi]]$postList), error=function(e) NULL); if (!is.null(pooled)) models[[mi]]$postList <- pooled }",
      "}",
      "save(models, file = file.path(modelDir, paste0('models_thin_', thin, '_samples_', samples, '_chains_', nChains, '.Rdata')))"
    )
    writeLines(s2, file.path(outdir, "workflow_scripts", "S2_fit_models.R"))

    s3 <- c(
      "# S3_evaluate_convergence.R -- generated by JSDMWorkbench",
      "set.seed(1)",
      "library(Hmsc); library(coda)",
      "modelDir <- file.path('.', 'models'); resultDir <- file.path('.', 'results'); tableDir <- file.path('.', 'tables'); predDir <- file.path('.', 'predictions')",
      "if (!dir.exists(resultDir)) dir.create(resultDir, recursive = TRUE)",
      paste0("samples <- ", samples), paste0("thin <- ", thin), paste0("nChains <- ", nChains),
      paste0("do_diagnostics <- ", if (isTRUE(cfg$outputs$diagnostics %||% TRUE)) "TRUE" else "FALSE"),
      paste0("showBeta <- ", if (isTRUE(cfg$outputs$convergence$showBeta %||% TRUE)) "TRUE" else "FALSE"),
      paste0("showGamma <- ", if (isTRUE(cfg$outputs$convergence$showGamma %||% TRUE)) "TRUE" else "FALSE"),
      paste0("showOmega <- ", if (isTRUE(cfg$outputs$convergence$showOmega %||% TRUE)) "TRUE" else "FALSE"),
      paste0("showRho <- ", if (isTRUE(cfg$outputs$convergence$showRho %||% TRUE)) "TRUE" else "FALSE"),
      paste0("showAlpha <- ", if (isTRUE(cfg$outputs$convergence$showAlpha %||% TRUE)) "TRUE" else "FALSE"),
      paste0("doEffective <- ", if (isTRUE(cfg$outputs$convergence$effectiveSize %||% TRUE)) "TRUE" else "FALSE"),
      paste0("doGelman <- ", if (isTRUE(cfg$outputs$convergence$gelmanPSRF %||% TRUE)) "TRUE" else "FALSE"),
      "load(file.path(modelDir, paste0('models_thin_', thin, '_samples_', samples, '_chains_', nChains, '.Rdata')))",
      "txt <- file.path(resultDir, 'MCMC_convergence.txt'); cat('MCMC convergence diagnostics\\n\\n', file = txt)",
      "pdf(file.path(resultDir, 'MCMC_convergence.pdf'))",
      "extract_mcmc_lists <- function(x) { out <- list(); walk <- function(z, nm) { if (inherits(z, 'mcmc.list')) out[[nm]] <<- z else if (is.list(z)) for (ii in seq_along(z)) walk(z[[ii]], paste0(nm, '_', ii)) }; walk(x, 'value'); out }",
      "if (!do_diagnostics) { plot.new(); title(main='MCMC diagnostics disabled'); text(0.5,0.5,'MCMC diagnostics disabled by output settings.'); cat('MCMC diagnostics disabled by output settings.\\n', file=txt, append=TRUE) } else",
      "for (j in seq_along(models)) {",
      "  mpost <- tryCatch(convertToCodaObject(models[[j]], spNamesNumbers = c(TRUE, FALSE), covNamesNumbers = c(TRUE, FALSE)), error = function(e) NULL)",
      "  targets <- c(Beta=showBeta, Gamma=showGamma, Omega=showOmega, Rho=showRho, Alpha=showAlpha)",
      "  cat('\\nModel: ', names(models)[j], '\\n', file = txt, append = TRUE)",
      "  if (is.null(mpost)) { cat('convertToCodaObject unavailable for this model.\\n', file=txt, append=TRUE); next }",
      "  for (parName in names(targets)[targets]) {",
      "    if (is.null(mpost[[parName]])) { cat(parName, 'not available.\\n', file=txt, append=TRUE); next }",
      "    lists <- extract_mcmc_lists(mpost[[parName]])",
      "    if (length(lists) == 0) { cat(parName, 'has no mcmc.list entries.\\n', file=txt, append=TRUE); next }",
      "    for (lnm in names(lists)) {",
      "      obj <- lists[[lnm]]",
      "      if (doGelman && nChains > 1) { psrf <- tryCatch(coda::gelman.diag(obj, multivariate=FALSE)$psrf[,1], error=function(e) NULL); if (!is.null(psrf)) { pn <- names(psrf); if (is.null(pn) || length(pn) != length(psrf)) pn <- paste0('par_', seq_along(psrf)); write.csv(data.frame(parameter=pn, PSRF=as.numeric(psrf)), file.path(tableDir, paste0('S3_', parName, '_PSRF_', lnm, '_', names(models)[j], '.csv')), row.names=FALSE); finite <- psrf[is.finite(psrf)]; if (length(finite)>0) hist(finite, main=paste('PSRF', parName, names(models)[j]), xlab='PSRF') else { plot.new(); title(main=paste('PSRF', parName, names(models)[j])); text(0.5,0.5,'No finite PSRF values available.') } } }",
      "      if (doEffective) { ess <- tryCatch(coda::effectiveSize(obj), error=function(e) NULL); if (!is.null(ess)) { en <- names(ess); if (is.null(en) || length(en) != length(ess)) en <- paste0('par_', seq_along(ess)); write.csv(data.frame(parameter=en, effectiveSize=as.numeric(ess)), file.path(tableDir, paste0('S3_', parName, '_effectiveSize_', lnm, '_', names(models)[j], '.csv')), row.names=FALSE) } }",
      "    }",
      "  }",
      "}",
      "dev.off()"
    )
    writeLines(s3, file.path(outdir, "workflow_scripts", "S3_evaluate_convergence.R"))

    s4 <- c(
      "# S4_compute_model_fit.R -- generated by JSDMWorkbench",
      "set.seed(1)",
      "library(Hmsc)",
      "modelDir <- file.path('.', 'models'); tableDir <- file.path('.', 'tables'); predDir <- file.path('.', 'predictions')",
      paste0("samples <- ", samples), paste0("thin <- ", thin), paste0("nChains <- ", nChains), paste0("nParallel <- ", nParallel), paste0("nfolds <- ", nfolds),
      paste0("compute_waic <- ", if (compute_waic) "TRUE" else "FALSE"),
      paste0("do_predicted <- ", if (isTRUE(cfg$outputs$predicted %||% TRUE)) "TRUE" else "FALSE"),
      paste0("do_fit <- ", if (isTRUE(cfg$outputs$fit %||% TRUE)) "TRUE" else "FALSE"),
      paste0("do_cv <- ", if (isTRUE(cfg$outputs$cv %||% TRUE)) "TRUE" else "FALSE"),
      paste0("partition_column <- '", gsub("'", "\\\\'", partition_column), "'"),
      "load(file.path(modelDir, paste0('models_thin_', thin, '_samples_', samples, '_chains_', nChains, '.Rdata')))",
      "MF <- list(); MFCV <- list(); WAIC <- list()",
      "for (mi in seq_along(models)) {",
      "  MF[mi] <- list(NULL); MFCV[mi] <- list(NULL); WAIC[mi] <- list('disabled')",
      "  preds <- if (do_predicted || do_fit || do_cv) tryCatch(computePredictedValues(models[[mi]]), error = function(e) { message('computePredictedValues failed for ', names(models)[mi], ': ', e$message); NULL }) else NULL",
      "  if (do_predicted && !is.null(preds)) saveRDS(preds, file.path(predDir, paste0('predicted_values_', names(models)[mi], '.rds')))",
      "  MF[mi] <- list(if (do_fit && !is.null(preds)) tryCatch(evaluateModelFit(hM = models[[mi]], predY = preds), error = function(e) list(fit_error = e$message)) else if (do_fit) list(fit_error = 'computePredictedValues failed') else NULL)",
      "  partition <- tryCatch({ args <- list(hM=models[[mi]], nfolds=nfolds); if (nzchar(partition_column) && !is.null(models[[mi]]$studyDesign) && partition_column %in% names(models[[mi]]$studyDesign)) args$column <- partition_column; do.call(createPartition, args) }, error = function(e) NULL)",
      "  if (do_cv && !is.null(partition) && !is.null(preds)) {",
      "    preds_cv <- tryCatch(computePredictedValues(models[[mi]], partition = partition, nParallel = nParallel), error = function(e) NULL)",
      "    MFCV[mi] <- list(if (!is.null(preds_cv)) tryCatch(evaluateModelFit(hM = models[[mi]], predY = preds_cv), error = function(e) list(fit_error = e$message)) else NULL)",
      "  }",
      "  WAIC[mi] <- list(if (compute_waic) tryCatch(computeWAIC(models[[mi]]), error = function(e) e$message) else 'disabled')",
      "}",
      "names(MF) <- names(models); names(MFCV) <- names(models); names(WAIC) <- names(models)",
      "save(MF, MFCV, WAIC, file = file.path(modelDir, paste0('MF_thin_', thin, '_samples_', samples, '_chains_', nChains, '_nfolds_', nfolds, '.Rdata')))"
    )
    writeLines(s4, file.path(outdir, "workflow_scripts", "S4_compute_model_fit.R"))

    s5 <- c(
      "# S5_show_model_fit.R -- generated by JSDMWorkbench",
      "library(Hmsc)",
      "modelDir <- file.path('.', 'models'); resultDir <- file.path('.', 'results')",
      paste0("samples <- ", samples), paste0("thin <- ", thin), paste0("nChains <- ", nChains), paste0("nfolds <- ", nfolds),
      "load(file.path(modelDir, paste0('MF_thin_', thin, '_samples_', samples, '_chains_', nChains, '_nfolds_', nfolds, '.Rdata')))",
      "pdf(file.path(resultDir, paste0('model_fit_nfolds_', nfolds, '.pdf')))",
      "safe_fit_plot <- function(x, y, xlim, ylim, xlab, ylab, main) { ok <- is.finite(x) & is.finite(y); if (any(ok)) { plot(x[ok], y[ok], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main); abline(0,1) } else { plot.new(); title(main=main); text(0.5, 0.5, 'No finite fit metrics for this quick-test run.') } }",
      "for (j in seq_along(MF)) {",
      "  cMF <- MF[[j]]; cMFCV <- MFCV[[j]]",
      "  if (!is.null(cMF$TjurR2) && !is.null(cMFCV$TjurR2)) safe_fit_plot(cMF$TjurR2, cMFCV$TjurR2, c(-1,1), c(-1,1), 'explanatory power', 'predictive power', paste0(names(MF)[j], ': Tjur R2'))",
      "  if (!is.null(cMF$R2) && !is.null(cMFCV$R2)) safe_fit_plot(cMF$R2, cMFCV$R2, c(-1,1), c(-1,1), 'explanatory power', 'predictive power', paste0(names(MF)[j], ': R2'))",
      "  if (!is.null(cMF$AUC) && !is.null(cMFCV$AUC)) safe_fit_plot(cMF$AUC, cMFCV$AUC, c(0,1), c(0,1), 'explanatory power', 'predictive power', paste0(names(MF)[j], ': AUC'))",
      "}",
      "dev.off()"
    )
    writeLines(s5, file.path(outdir, "workflow_scripts", "S5_show_model_fit.R"))

    phylo_beta_helper <- c(
      "# HMSC_phylo_beta_heatmap_helpers.R -- generated by JSDM Studio",
      "prep_beta_heatmap <- function(m, postBeta, include_intercept=FALSE) {",
      "  if (is.null(postBeta) || is.null(postBeta$mean)) return(NULL)",
      "  beta <- as.matrix(postBeta$mean)",
      "  rn <- rownames(beta); if (is.null(rn) || length(rn) != nrow(beta) || any(!nzchar(rn))) rn <- if (!is.null(m$covNames)) m$covNames else if (!is.null(m$X)) colnames(m$X) else paste0('covariate_', seq_len(nrow(beta)))",
      "  cn <- colnames(beta); if (is.null(cn) || length(cn) != ncol(beta) || any(!nzchar(cn))) cn <- if (!is.null(m$spNames)) m$spNames else paste0('species_', seq_len(ncol(beta)))",
      "  rownames(beta) <- make.unique(as.character(rn)); colnames(beta) <- make.unique(as.character(cn))",
      "  support <- NULL; if (!is.null(postBeta$support)) { support <- as.matrix(postBeta$support); if (identical(dim(support), dim(beta))) { rownames(support) <- rownames(beta); colnames(support) <- colnames(beta) } else support <- NULL }",
      "  if (!isTRUE(include_intercept) && nrow(beta) > 1) { keep <- !tolower(rownames(beta)) %in% c('(intercept)', 'intercept'); if (any(keep)) { beta <- beta[keep,,drop=FALSE]; if (!is.null(support)) support <- support[keep,,drop=FALSE] } }",
      "  list(beta=beta, support=support, tree=m$phyloTree)",
      "}",
      "write_beta_heatmap_long <- function(beta, support, path, model_name='main') {",
      "  if (is.null(beta) || nrow(beta)==0 || ncol(beta)==0) return(FALSE)",
      "  grid <- expand.grid(predictor=rownames(beta), species=colnames(beta), stringsAsFactors=FALSE); grid$model <- model_name",
      "  grid$estimate <- as.numeric(beta[cbind(match(grid$predictor, rownames(beta)), match(grid$species, colnames(beta)))])",
      "  grid$support <- if (!is.null(support)) as.numeric(support[cbind(match(grid$predictor, rownames(support)), match(grid$species, colnames(support)))]) else NA_real_",
      "  write.csv(grid[, c('model','species','predictor','estimate','support')], path, row.names=FALSE); TRUE",
      "}",
      "draw_phylo_beta_heatmap <- function(m, postBeta, model_name='main', include_intercept=FALSE, support_level=0.95, show_species_names=TRUE) {",
      "  dat <- prep_beta_heatmap(m, postBeta, include_intercept); if (is.null(dat)) { plot.new(); title(main=paste0(model_name, ': phylogeny + Beta heatmap')); text(0.5,0.5,'Beta estimates unavailable'); return(FALSE) }",
      "  beta <- dat$beta; support <- dat$support; tree <- dat$tree; tree_available <- !is.null(tree) && requireNamespace('ape', quietly=TRUE); tree_note <- 'species order: input order'",
      "  if (tree_available) { keep <- intersect(tree$tip.label, colnames(beta)); if (length(keep) >= 2) { tree <- ape::keep.tip(tree, keep); beta <- beta[, tree$tip.label, drop=FALSE]; if (!is.null(support)) support <- support[, tree$tip.label, drop=FALSE]; tree_note <- 'species order: phylogeny' } else tree_available <- FALSE }",
      "  if (!tree_available && ncol(beta) > 2) { hc <- tryCatch(stats::hclust(stats::dist(t(beta))), error=function(e) NULL); if (!is.null(hc)) { beta <- beta[, hc$labels[hc$order], drop=FALSE]; if (!is.null(support)) support <- support[, colnames(beta), drop=FALSE]; tree_note <- 'species order: Beta similarity' } }",
      "  finite_beta <- beta[is.finite(beta)]; if (length(finite_beta)==0) { plot.new(); title(main=paste0(model_name, ': phylogeny + Beta heatmap')); text(0.5,0.5,'No finite Beta estimates'); return(FALSE) }",
      "  cap <- suppressWarnings(as.numeric(stats::quantile(abs(finite_beta), 0.98, na.rm=TRUE))); if (!is.finite(cap) || cap <= 0) cap <- max(abs(finite_beta), na.rm=TRUE); if (!is.finite(cap) || cap <= 0) cap <- 1",
      "  beta_clip <- matrix(pmax(-cap, pmin(cap, beta)), nrow=nrow(beta), ncol=ncol(beta), dimnames=dimnames(beta)); pal <- grDevices::colorRampPalette(c('#2b6f9e','#f7f7f5','#b84a3a'))(101); graphics::layout(matrix(c(1,2,3), nrow=1), widths=c(if (tree_available) 2.25 else 0.25, 3.4, 0.45)); on.exit(graphics::layout(1), add=TRUE); op <- par(oma=c(0,0,2.8,0), xpd=NA); on.exit(par(op), add=TRUE)",
      "  n_sp <- ncol(beta_clip); n_cov <- nrow(beta_clip)",
      "  if (tree_available) { par(mar=c(7,0.2,2.5,0.1)); ape::plot.phylo(tree, direction='rightwards', show.tip.label=FALSE, no.margin=FALSE, edge.color='#8a8078', edge.width=0.8); lp <- get('last_plot.phylo', envir=get('.PlotPhyloEnv', envir=asNamespace('ape'))); tip_order <- tree$tip.label[order(lp$yy[seq_len(ape::Ntip(tree))])]; beta_clip <- beta_clip[, tip_order, drop=FALSE]; if (!is.null(support)) support <- support[, tip_order, drop=FALSE]; title('phylogeny', line=0.4, cex.main=0.9) } else { par(mar=c(7,0.2,2.5,0.1)); plot.new(); text(0.5,0.5,'No\\nphylogeny', cex=0.75, col='#6f6f6f') }",
      "  par(mar=c(7,0.4,2.5, if (isTRUE(show_species_names) && n_sp <= 80) 5.0 else 0.6)); image(x=seq_len(n_cov), y=seq_len(n_sp), z=beta_clip, axes=FALSE, xlab='', ylab='', col=pal, zlim=c(-cap,cap), useRaster=TRUE); axis(1, at=seq_len(n_cov), labels=rownames(beta_clip), las=2, cex.axis=0.62); if (isTRUE(show_species_names) && n_sp <= 80) axis(4, at=seq_len(n_sp), labels=colnames(beta_clip), las=2, cex.axis=0.42, tick=FALSE, line=-0.2); box(col='#9a9a9a'); abline(v=seq(0.5,n_cov+0.5,by=1), col=grDevices::adjustcolor('white', alpha.f=0.55), lwd=0.35)",
      "  if (!is.null(support) && identical(dim(support), dim(beta_clip))) { supported <- which(support >= support_level, arr.ind=TRUE); if (nrow(supported) > 0) points(supported[,1], supported[,2], pch=16, cex=0.22, col=grDevices::adjustcolor('black', alpha.f=0.65)) }",
      "  title('environmental effects (Beta)', line=0.4, cex.main=0.9); par(mar=c(7,0.1,2.5,2.2)); ly <- seq(-cap,cap,length.out=length(pal)); image(x=1, y=ly, z=matrix(ly,nrow=1), axes=FALSE, xlab='', ylab='', col=pal, zlim=c(-cap,cap), useRaster=TRUE); axis(4, at=pretty(c(-cap,cap), n=5), labels=round(pretty(c(-cap,cap), n=5),2), cex.axis=0.65); title('Beta', line=0.4, cex.main=0.8)",
      "  mtext(paste0(model_name, ': phylogeny-ordered species-environment Beta heatmap'), outer=TRUE, side=3, line=1.1, cex=1.0, font=2); mtext(paste0('Posterior mean Beta; colors clipped at |', round(cap,2), '| for display. Black dots: support >= ', support_level, '. ', tree_note, '.'), outer=TRUE, side=3, line=-0.3, cex=0.68); TRUE",
      "}",
      "save_phylo_beta_heatmap <- function(m, postBeta, model_name='main', result_dir='results') {",
      "  base <- file.path(result_dir, paste0('parameter_estimates_phylo_beta_heatmap_', model_name)); grDevices::pdf(paste0(base,'.pdf'), width=11, height=7.5); draw_phylo_beta_heatmap(m, postBeta, model_name=model_name, show_species_names=TRUE); grDevices::dev.off(); grDevices::png(paste0(base,'.png'), width=2600, height=1700, res=220); draw_phylo_beta_heatmap(m, postBeta, model_name=model_name, show_species_names=TRUE); grDevices::dev.off(); TRUE",
      "}"
    )
    writeLines(phylo_beta_helper, file.path(outdir, "workflow_scripts", "HMSC_phylo_beta_heatmap_helpers.R"))
    bookstyle_helper_functions <- c(
      "hmsc_summarise_prediction_object",
      "hmsc_plot_species_order",
      "hmsc_extract_mcmc_numeric",
      "hmsc_prepare_gamma_heatmap_data",
      "hmsc_write_gamma_heatmap_long",
      "hmsc_pretty_hmsc_terms",
      "hmsc_draw_trait_environment_gamma_heatmap",
      "hmsc_save_trait_environment_gamma_heatmap",
      "hmsc_prepare_rho_signal_data",
      "hmsc_draw_rho_signal_curve",
      "hmsc_save_rho_signal_curve",
      "hmsc_prepare_vp_book_tables",
      "hmsc_draw_vp_bookstyle",
      "hmsc_save_vp_bookstyle",
      "hmsc_write_association_long",
      "hmsc_draw_association_heatmap",
      "hmsc_save_association_heatmap",
      "hmsc_gradient_axis",
      "hmsc_gradient_trait_matrix",
      "hmsc_gradient_response_table",
      "hmsc_panel_curve",
      "hmsc_draw_gradient_response_panel",
      "hmsc_save_gradient_response_panel"
    )
    bookstyle_helper <- c(
      "# HMSC_bookstyle_plot_helpers.R -- generated by JSDM Studio",
      "`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x",
      capture.output(dump(bookstyle_helper_functions, file = ""))
    )
    writeLines(bookstyle_helper, file.path(outdir, "workflow_scripts", "HMSC_bookstyle_plot_helpers.R"))

    s6 <- c(
      "# S6_show_parameter_estimates.R -- generated by JSDMWorkbench",
      "library(Hmsc)",
      "if (requireNamespace('corrplot', quietly = TRUE)) library(corrplot)",
      "modelDir <- file.path('.', 'models'); resultDir <- file.path('.', 'results')",
      paste0("samples <- ", samples), paste0("thin <- ", thin), paste0("nChains <- ", nChains),
      paste0("do_params <- ", if (isTRUE(cfg$outputs$parameters %||% TRUE)) "TRUE" else "FALSE"),
      paste0("do_vp <- ", if (isTRUE(cfg$outputs$variance_partitioning %||% TRUE)) "TRUE" else "FALSE"),
      paste0("do_omega <- ", if (isTRUE(cfg$outputs$omega %||% TRUE)) "TRUE" else "FALSE"),
      paste0("do_rho <- ", if (isTRUE(cfg$outputs$convergence$showRho %||% TRUE)) "TRUE" else "FALSE"),
      paste0("plot_beta <- ", if (isTRUE(cfg$outputs$plotting$plotBeta %||% TRUE)) "TRUE" else "FALSE"),
      paste0("plot_gamma <- ", if (isTRUE(cfg$outputs$plotting$plotGamma %||% TRUE)) "TRUE" else "FALSE"),
      paste0("show_sp_beta <- ", if (isTRUE(cfg$outputs$plotting$show.sp.names.beta %||% FALSE)) "TRUE" else "FALSE"),
      paste0("plot_tree <- ", if (isTRUE(cfg$outputs$plotting$plotTree %||% FALSE)) "TRUE" else "FALSE"),
      paste0("show_sp_omega <- ", if (isTRUE(cfg$outputs$plotting$show.sp.names.omega %||% TRUE)) "TRUE" else "FALSE"),
      paste0("omega_order <- '", gsub("'", "\\\\'", hmsc_normalize_omega_order(cfg$outputs$plotting$omega.order %||% "original")), "'"),
      paste0("maxOmega <- ", as.integer(cfg$outputs$convergence$maxOmega %||% 10)),
      "quiet_optional <- function(expr, fallback=NULL) { tf <- tempfile(); con <- file(tf, open='wt'); sink(con, type='message'); on.exit({ sink(type='message'); close(con); unlink(tf) }, add=TRUE); tryCatch(suppressWarnings(suppressMessages(force(expr))), error=function(e) fallback) }",
      "helper <- file.path('workflow_scripts', 'HMSC_phylo_beta_heatmap_helpers.R'); if (file.exists(helper)) source(helper)",
      "book_helper <- file.path('workflow_scripts', 'HMSC_bookstyle_plot_helpers.R'); if (file.exists(book_helper)) source(book_helper)",
      "load(file.path(modelDir, paste0('models_thin_', thin, '_samples_', samples, '_chains_', nChains, '.Rdata')))",
      "txt <- file.path(resultDir, 'parameter_estimates.txt'); cat('Parameter estimates\\n\\n', file=txt)",
      "pdf(file.path(resultDir, 'parameter_estimates.pdf'))",
      "for (j in seq_along(models)) {",
      "  m <- models[[j]]",
      "  VP <- if (do_params && do_vp) quiet_optional(computeVariancePartitioning(m)) else NULL",
      "  if (!is.null(VP)) { write.csv(VP$vals, file.path(resultDir, paste0('parameter_estimates_VP_', names(models)[j], '.csv'))); if (exists('hmsc_draw_vp_bookstyle')) quiet_optional(hmsc_draw_vp_bookstyle(VP, model_name=names(models)[j])); if (exists('hmsc_save_vp_bookstyle')) quiet_optional(hmsc_save_vp_bookstyle(VP, model_name=names(models)[j], result_dir=resultDir)); quiet_optional(plotVariancePartitioning(hM=m, VP=VP, main=paste0('Variance partitioning: ', names(models)[j]))) }",
      "  if (do_params && do_rho && exists('hmsc_draw_rho_signal_curve')) { quiet_optional(hmsc_draw_rho_signal_curve(m, model_name=names(models)[j])); if (exists('hmsc_save_rho_signal_curve')) quiet_optional(hmsc_save_rho_signal_curve(m, model_name=names(models)[j], result_dir=resultDir)) }",
      "  postBeta <- if (do_params) quiet_optional(getPostEstimate(m, parName='Beta')) else NULL",
      "  if (!is.null(postBeta)) {",
      "    write.csv(as.data.frame(t(postBeta$mean)), file.path(resultDir, paste0('parameter_estimates_Beta_mean_', names(models)[j], '.csv')))",
      "    if (exists('prep_beta_heatmap')) { bh <- prep_beta_heatmap(m, postBeta, include_intercept=FALSE); if (!is.null(bh)) write_beta_heatmap_long(bh$beta, bh$support, file.path(resultDir, paste0('parameter_estimates_Beta_species_environment_long_', names(models)[j], '.csv')), model_name=names(models)[j]); if (plot_beta) { quiet_optional(draw_phylo_beta_heatmap(m, postBeta, model_name=names(models)[j], include_intercept=FALSE, support_level=0.95, show_species_names=TRUE)); quiet_optional(save_phylo_beta_heatmap(m, postBeta, model_name=names(models)[j], result_dir=resultDir)) } }",
      "    if (plot_beta) quiet_optional(plotBeta(m, post=postBeta, supportLevel=0.95, param='Sign', plotTree=plot_tree, spNamesNumbers=c(TRUE, show_sp_beta)))",
      "  }",
      "  postGamma <- if (do_params) quiet_optional(getPostEstimate(m, parName='Gamma')) else NULL",
      "  if (!is.null(postGamma)) {",
      "    write.csv(as.data.frame(t(postGamma$mean)), file.path(resultDir, paste0('parameter_estimates_Gamma_mean_', names(models)[j], '.csv')))",
      "    if (!is.null(postGamma$support)) write.csv(as.data.frame(t(postGamma$support)), file.path(resultDir, paste0('parameter_estimates_Gamma_support_', names(models)[j], '.csv')))",
      "    if (exists('hmsc_prepare_gamma_heatmap_data')) { gh_all <- hmsc_prepare_gamma_heatmap_data(m, postGamma, include_intercepts=TRUE); if (!is.null(gh_all)) quiet_optional(hmsc_write_gamma_heatmap_long(gh_all$gamma_all, gh_all$support_all, file.path(resultDir, paste0('parameter_estimates_trait_environment_gamma_long_', names(models)[j], '.csv')), model_name=names(models)[j])) }",
      "    if (exists('hmsc_save_trait_environment_gamma_heatmap')) quiet_optional(hmsc_save_trait_environment_gamma_heatmap(m, postGamma, model_name=names(models)[j], result_dir=resultDir, include_intercepts=FALSE, support_level=0.95))",
      "    if (plot_gamma && exists('hmsc_draw_trait_environment_gamma_heatmap')) quiet_optional(hmsc_draw_trait_environment_gamma_heatmap(m, postGamma, model_name=names(models)[j], include_intercepts=FALSE, support_level=0.95))",
      "    if (plot_gamma) quiet_optional(plotGamma(m, post=postGamma, supportLevel=0.95, param='Sign'))",
      "  }",
      "  OmegaCor <- if (do_params && do_omega) quiet_optional(computeAssociations(m)) else NULL",
      "  if (!is.null(OmegaCor)) for (r in seq_along(OmegaCor)) { level_name <- names(OmegaCor)[r]; if (is.null(level_name) || length(level_name)==0 || is.na(level_name) || !nzchar(level_name)) level_name <- paste0('random_level_', r); write.csv(OmegaCor[[r]]$mean, file.path(resultDir, paste0('parameter_estimates_Omega_mean_', names(models)[j], '_', level_name, '.csv'))); if (exists('hmsc_draw_association_heatmap')) quiet_optional(hmsc_draw_association_heatmap(m, OmegaCor[[r]]$mean, OmegaCor[[r]]$support, model_name=names(models)[j], level_name=level_name, show_species_names=show_sp_omega)); if (exists('hmsc_save_association_heatmap')) quiet_optional(hmsc_save_association_heatmap(m, OmegaCor[[r]]$mean, OmegaCor[[r]]$support, model_name=names(models)[j], level_name=level_name, result_dir=resultDir, show_species_names=show_sp_omega)); if (requireNamespace('corrplot', quietly = TRUE)) { mat <- OmegaCor[[r]]$mean; if (is.finite(maxOmega) && nrow(mat) > maxOmega && maxOmega > 1) mat <- mat[seq_len(maxOmega), seq_len(maxOmega), drop=FALSE]; tl_pos <- if (show_sp_omega) 'lt' else 'n'; quiet_optional(corrplot::corrplot(mat, method='color', order=omega_order, tl.pos=tl_pos)) } }",
      "}",
      "dev.off()"
    )
    writeLines(s6, file.path(outdir, "workflow_scripts", "S6_show_parameter_estimates.R"))

    s7 <- c(
      "# S7_make_predictions.R -- generated by JSDMWorkbench",
      "library(Hmsc)",
      "modelDir <- file.path('.', 'models'); resultDir <- file.path('.', 'results'); predDir <- file.path('.', 'predictions')",
      paste0("samples <- ", samples), paste0("thin <- ", thin), paste0("nChains <- ", nChains),
      paste0("env_list <- '", gsub("'", "\\\\'", env_list), "'"),
      paste0("species_list <- '", gsub("'", "\\\\'", species_list), "'"),
      paste0("trait_list <- '", gsub("'", "\\\\'", trait_list), "'"),
      paste0("computeSAIR <- ", if (computeSAIR) "TRUE" else "FALSE"),
      paste0("do_gradients <- ", if (isTRUE(cfg$outputs$gradients %||% FALSE)) "TRUE" else "FALSE"),
      paste0("nParallel <- ", nParallel),
      paste0("lon_col <- '", gsub("'", "\\\\'", lon_col), "'"),
      paste0("lat_col <- '", gsub("'", "\\\\'", lat_col), "'"),
      "quiet_optional <- function(expr, fallback=NULL) { tf <- tempfile(); con <- file(tf, open='wt'); sink(con, type='message'); on.exit({ sink(type='message'); close(con); unlink(tf) }, add=TRUE); tryCatch(suppressWarnings(suppressMessages(force(expr))), error=function(e) fallback) }",
      "book_helper <- file.path('workflow_scripts', 'HMSC_bookstyle_plot_helpers.R'); if (file.exists(book_helper)) source(book_helper)",
      "summarise_prediction <- function(pred, site_ids=NULL, species_ids=NULL) {",
      "  if (is.list(pred) && length(pred)>0 && all(vapply(pred, function(x) is.matrix(x) || is.data.frame(x), logical(1)))) {",
      "    mats <- lapply(pred, as.matrix); arr <- array(unlist(lapply(mats, as.numeric)), dim=c(nrow(mats[[1]]), ncol(mats[[1]]), length(mats)))",
      "    dimnames(arr) <- list(rownames(mats[[1]]), colnames(mats[[1]]), names(pred)); out <- apply(arr, c(1,2), mean, na.rm=TRUE)",
      "  } else if (length(dim(pred)) == 3) out <- apply(pred, c(1,2), mean, na.rm=TRUE) else out <- as.matrix(pred)",
      "  if (!is.null(site_ids) && nrow(out) == length(site_ids)) rownames(out) <- site_ids",
      "  if (!is.null(species_ids) && ncol(out) == length(species_ids)) colnames(out) <- species_ids",
      "  out",
      "}",
      "plot_prediction_heatmap <- function(pm, main, observed=NULL, value_label='predicted mean') {",
      "  mat <- suppressWarnings(apply(as.matrix(pm), c(1,2), as.numeric)); if (is.null(dim(mat))) mat <- matrix(mat, nrow=nrow(pm))",
      "  if (!is.null(rownames(pm))) rownames(mat) <- rownames(pm); if (!is.null(colnames(pm))) colnames(mat) <- colnames(pm)",
      "  pick_ticks <- function(n, max_ticks=14L) { if (n <= max_ticks) seq_len(n) else unique(pmax(1L, pmin(n, round(seq(1, n, length.out=max_ticks))))) }",
      "  pal <- grDevices::hcl.colors(101, 'YlGnBu', rev=TRUE); vr <- range(mat, finite=TRUE, na.rm=TRUE)",
      "  op <- par(no.readonly=TRUE); on.exit(par(op), add=TRUE); par(mar=c(8,7,4,7))",
      "  image(t(mat[nrow(mat):1,,drop=FALSE]), axes=FALSE, col=pal, zlim=if (all(is.finite(vr)) && diff(vr)>0) vr else NULL, main=main, xlab='Species', ylab='Site / prediction point')",
      "  x_idx <- pick_ticks(ncol(mat), 18L); y_idx <- pick_ticks(nrow(mat), 20L)",
      "  axis(1, at=if (ncol(mat)==1) 0.5 else (x_idx-1)/(ncol(mat)-1), labels=abbreviate(colnames(mat)[x_idx], minlength=10), las=2, cex.axis=0.75)",
      "  axis(2, at=if (nrow(mat)==1) 0.5 else (y_idx-1)/(nrow(mat)-1), labels=rev(rownames(mat))[y_idx], las=2, cex.axis=0.75)",
      "  if (nrow(mat) > length(y_idx)) mtext(paste0(nrow(mat), ' prediction rows; axis labels sampled'), side=2, line=5.2, cex=0.65)",
      "  mtext(value_label, side=4, line=4)",
      "  if (!is.null(observed) && all(dim(observed) == dim(mat))) { plot(as.numeric(mat), as.numeric(observed), pch=16, col=rgb(0.1,0.45,0.55,0.45), xlab=value_label, ylab='observed', main=paste0(main, ' - observed vs predicted')); abline(0,1,col='gray50',lty=2); abline(h=0.5,col='gray85',lty=3); abline(v=0.5,col='gray85',lty=3); legend('topleft', bty='n', legend=paste0('mean residual = ', round(mean(as.numeric(observed)-as.numeric(mat), na.rm=TRUE), 3))) }",
      "}",
      "prepare_spatial_xy <- function(dat, site_ids, lon_col, lat_col) {",
      "  if (ncol(dat) > 1 && names(dat)[1] %in% c('', 'X', 'row.names', 'rowname', 'row_id', 'id', 'site_id', 'sample_id')) { rownames(dat) <- make.unique(as.character(dat[[1]])); dat <- dat[-1] }",
      "  if (!all(c(lon_col, lat_col) %in% colnames(dat)) && all(c('x','y') %in% colnames(dat))) { lon_col <- 'x'; lat_col <- 'y' }",
      "  if (!all(c(lon_col, lat_col) %in% colnames(dat))) stop('newcoords.csv must include spatial coordinate columns.', call.=FALSE)",
      "  xy <- as.matrix(dat[, c(lon_col, lat_col), drop=FALSE]); mode(xy) <- 'numeric'",
      "  if (!is.null(site_ids) && length(site_ids) == nrow(xy)) rownames(xy) <- site_ids",
      "  xy",
      "}",
      "plot_newdata_spatial <- function(pm, coord, site_ids, lon_col, lat_col, main) {",
      "  xy <- tryCatch(prepare_spatial_xy(coord, site_ids, lon_col, lat_col), error=function(e) NULL); if (is.null(xy) || nrow(xy) != nrow(pm)) return(FALSE)",
      "  mat <- suppressWarnings(apply(as.matrix(pm), c(1,2), as.numeric)); score <- rowMeans(mat, na.rm=TRUE); score_label <- 'mean across species'; finite_vals <- score[is.finite(score)]",
      "  if (length(finite_vals)==0) return(FALSE)",
      "  if (diff(range(finite_vals, na.rm=TRUE)) == 0 && ncol(mat) > 1) { ssd <- apply(mat, 2, stats::sd, na.rm=TRUE); if (any(is.finite(ssd)) && max(ssd, na.rm=TRUE) > 0) { best <- which.max(ssd); score <- mat[, best]; score_label <- paste0('species: ', colnames(mat)[best]); finite_vals <- score[is.finite(score)] } }",
      "  pal <- grDevices::hcl.colors(101, 'YlGnBu', rev=FALSE); vr <- range(finite_vals, na.rm=TRUE); idx <- if (diff(vr)==0) rep(51L, length(score)) else pmax(1L, pmin(101L, round(1 + 100 * (score - vr[1]) / diff(vr))))",
      "  point_cex <- if (nrow(xy) > 5000) 0.45 else if (nrow(xy) > 1000) 0.65 else if (nrow(xy) > 250) 1.0 else 2.4; point_pch <- if (nrow(xy) > 1000) 16 else 21",
      "  if (point_pch == 16) plot(xy[,1], xy[,2], pch=point_pch, col=pal[idx], cex=point_cex, asp=1, xlab=colnames(xy)[1], ylab=colnames(xy)[2], main=paste0(main, '\\n', score_label)) else plot(xy[,1], xy[,2], pch=point_pch, bg=pal[idx], col='#1b4d59', cex=point_cex, asp=1, xlab=colnames(xy)[1], ylab=colnames(xy)[2], main=paste0(main, '\\n', score_label))",
      "  if (nrow(xy) <= 250) text(xy[,1], xy[,2], labels=seq_len(nrow(xy)), cex=0.65, pos=3); legend('topright', bty='n', legend=c(score_label, paste0('points = ', nrow(xy)), paste0('mean = ', round(mean(finite_vals), 3)), paste0('range = ', paste(round(range(finite_vals), 3), collapse='-')))); TRUE",
      "}",
      "plot_newdata_score <- function(score, coord, site_ids, lon_col, lat_col, main, score_label='score') {",
      "  xy <- tryCatch(prepare_spatial_xy(coord, site_ids, lon_col, lat_col), error=function(e) NULL); score <- suppressWarnings(as.numeric(score)); if (is.null(xy) || nrow(xy) != length(score)) return(FALSE)",
      "  finite_vals <- score[is.finite(score)]; if (length(finite_vals)==0) return(FALSE)",
      "  pal <- grDevices::hcl.colors(101, 'YlGnBu', rev=FALSE); vr <- range(finite_vals, na.rm=TRUE); idx <- if (diff(vr)==0) rep(51L, length(score)) else pmax(1L, pmin(101L, round(1 + 100 * (score - vr[1]) / diff(vr))))",
      "  point_cex <- if (nrow(xy) > 5000) 0.45 else if (nrow(xy) > 1000) 0.65 else if (nrow(xy) > 250) 1.0 else 2.4; point_pch <- if (nrow(xy) > 1000) 16 else 21",
      "  if (point_pch == 16) plot(xy[,1], xy[,2], pch=point_pch, col=pal[idx], cex=point_cex, asp=1, xlab=colnames(xy)[1], ylab=colnames(xy)[2], main=paste0(main, '\\n', score_label)) else plot(xy[,1], xy[,2], pch=point_pch, bg=pal[idx], col='#1b4d59', cex=point_cex, asp=1, xlab=colnames(xy)[1], ylab=colnames(xy)[2], main=paste0(main, '\\n', score_label))",
      "  if (nrow(xy) <= 250) text(xy[,1], xy[,2], labels=seq_len(nrow(xy)), cex=0.65, pos=3); legend('topright', bty='n', legend=c(score_label, paste0('points = ', nrow(xy)), paste0('mean = ', round(mean(finite_vals), 3)), paste0('range = ', paste(round(range(finite_vals), 3), collapse='-')))); TRUE",
      "}",
      "load(file.path(modelDir, paste0('models_thin_', thin, '_samples_', samples, '_chains_', nChains, '.Rdata')))",
      "prediction_plot_count <- 0L",
      "pdf(file.path(resultDir, 'predictions.pdf'))",
      "for (j in seq_along(models)) {",
      "  m <- models[[j]]",
      "  covariates <- intersect(all.vars(m$XFormula), colnames(m$XData))",
      "  if (length(covariates) == 0) covariates <- colnames(m$XData)",
      "  env_filter <- trimws(unlist(strsplit(env_list, '[,;]+'))); env_filter <- env_filter[nzchar(env_filter)]; if (length(env_filter) > 0) { keep <- intersect(covariates, env_filter); if (length(keep) > 0) covariates <- keep }",
      "  write.csv(data.frame(model=names(models)[j], species_list=species_list, trait_list=trait_list, env_list=paste(covariates, collapse=','), computeSAIR=computeSAIR), file.path('tables', paste0('S7_prediction_settings_', names(models)[j], '.csv')), row.names=FALSE)",
      "  if (!do_gradients) {",
      "    pred_file <- file.path(resultDir, paste0('predicted_values_', names(models)[j], '.rds'))",
      "    pred_obj <- if (file.exists(pred_file)) readRDS(pred_file) else quiet_optional(computePredictedValues(m, expected=TRUE), NULL)",
      "    if (!is.null(pred_obj)) { pm_fit <- summarise_prediction(pred_obj, rownames(m$Y), m$spNames); write.csv(pm_fit, file.path(predDir, paste0('fitted_predicted_mean_', names(models)[j], '.csv'))); plot_prediction_heatmap(pm_fit, main=paste0(names(models)[j], ': fitted prediction mean'), observed=m$Y, value_label='fitted predicted mean'); prediction_plot_count <- prediction_plot_count + 1L } else { plot.new(); title(main=paste0(names(models)[j], ': no fitted predictions available')) }",
      "    next",
      "  }",
      "  for (covariate in covariates) {",
      "    Gradient <- quiet_optional(constructGradient(m, focalVariable = covariate))",
      "    if (!is.null(Gradient)) { predY <- quiet_optional(predict(m, Gradient=Gradient, expected=TRUE)); if (!is.null(predY)) { saveRDS(predY, file.path(predDir, paste0('gradient_prediction_', names(models)[j], '_', covariate, '.rds'))); sp_vec <- trimws(unlist(strsplit(species_list, '[,;]+'))); sp_vec <- sp_vec[nzchar(sp_vec)]; tr_vec <- trimws(unlist(strsplit(trait_list, '[,;]+'))); tr_vec <- tr_vec[nzchar(tr_vec)]; trait_data <- if (file.exists(file.path('data', 'TrData.csv'))) read.csv(file.path('data', 'TrData.csv'), row.names=1, check.names=FALSE, stringsAsFactors=FALSE) else NULL; if (exists('hmsc_draw_gradient_response_panel') && isTRUE(quiet_optional(hmsc_draw_gradient_response_panel(m, Gradient, predY, covariate, model_name=names(models)[j], TrData=trait_data, species_list=sp_vec, trait_list=tr_vec), FALSE))) prediction_plot_count <- prediction_plot_count + 1L; if (exists('hmsc_save_gradient_response_panel')) quiet_optional(hmsc_save_gradient_response_panel(m, Gradient, predY, covariate, model_name=names(models)[j], result_dir=predDir, TrData=trait_data, species_list=sp_vec, trait_list=tr_vec)); quiet_optional(plotGradient(m, Gradient, pred=predY, yshow=0, measure='S', showData=TRUE, main=paste0(names(models)[j], ': ', covariate))) } }",
      "  }",
      "}",
      "newdata_file <- file.path('data', 'newdata.csv')",
      "if (file.exists(newdata_file)) {",
      "  nd <- read.csv(newdata_file, check.names=FALSE, stringsAsFactors=FALSE)",
      "  if (ncol(nd) > 1 && names(nd)[1] %in% c('', 'X', 'row.names', 'rowname', 'row_id', 'id', 'site_id', 'sample_id')) { rownames(nd) <- make.unique(as.character(nd[[1]])); nd <- nd[-1] }",
      "  new_rows <- list()",
      "  for (j in seq_along(models)) {",
      "    m <- models[[j]]",
      "    ndj <- nd",
      "    keep <- intersect(colnames(m$XData), colnames(ndj)); ndj <- ndj[, keep, drop=FALSE]",
      "    miss <- setdiff(all.vars(m$XFormula), colnames(ndj)); if (length(miss) > 0 && paste(deparse(m$XFormula), collapse='') != '~.') { writeLines(paste('newdata missing columns:', paste(miss, collapse=',')), file.path('diagnostics', paste0('S7_newdata_prediction_error_', names(models)[j], '.txt'))); next }",
      "    for (nm in colnames(ndj)) if (is.factor(m$XData[[nm]])) ndj[[nm]] <- factor(as.character(ndj[[nm]]), levels=levels(m$XData[[nm]])) else if (is.numeric(m$XData[[nm]])) ndj[[nm]] <- suppressWarnings(as.numeric(ndj[[nm]]))",
      "    if (any(is.na(ndj))) { writeLines('newdata contains NA or unseen factor levels after matching training XData.', file.path('diagnostics', paste0('S7_newdata_prediction_error_', names(models)[j], '.txt'))); next }",
      "    source(file.path('workflow_scripts', 'HMSC_random_level_helpers.R'))",
      paste0("    random_cfg <- list(model = ", paste(capture.output(dput(cfg$model)), collapse = "\n"), ")"),
      "    pred_sd <- hmsc_read_design_csv(file.path('data', 'prediction_studyDesign.csv'))",
      "    pred_coord <- hmsc_read_design_csv(file.path('data', 'newcoords.csv'))",
      "    pred_new <- tryCatch(hmsc_predict_new_sites(m, ndj, pred_sd, pred_coord, random_cfg, rownames(ndj), nParallel), error=function(e) e)",
      "    if (inherits(pred_new, 'error')) { writeLines(conditionMessage(pred_new), file.path('diagnostics', paste0('S7_newdata_prediction_error_', names(models)[j], '.txt'))); next }",
      "    pm <- if (is.list(pred_new) && length(pred_new)>0 && all(vapply(pred_new, function(x) is.matrix(x) || is.data.frame(x), logical(1)))) { mats <- lapply(pred_new, as.matrix); arr <- array(unlist(lapply(mats, as.numeric)), dim=c(nrow(mats[[1]]), ncol(mats[[1]]), length(mats))); dimnames(arr) <- list(rownames(mats[[1]]), colnames(mats[[1]]), names(pred_new)); apply(arr, c(1,2), mean, na.rm=TRUE) } else if (length(dim(pred_new)) == 3) apply(pred_new, c(1,2), mean, na.rm=TRUE) else as.matrix(pred_new)",
      "    rownames(pm) <- rownames(ndj); if (is.null(colnames(pm))) colnames(pm) <- m$spNames",
      "    write.csv(pm, file.path(predDir, paste0('newdata_predicted_mean_', names(models)[j], '.csv')))",
      "    plot_prediction_heatmap(pm, main=paste0(names(models)[j], ': newdata predicted mean'), observed=NULL, value_label='newdata predicted mean'); prediction_plot_count <- prediction_plot_count + 1L",
      "    if (file.exists(file.path('data', 'newcoords.csv'))) { coord_new <- read.csv(file.path('data', 'newcoords.csv'), check.names=FALSE, stringsAsFactors=FALSE); if (isTRUE(plot_newdata_spatial(pm, coord_new, rownames(pm), lon_col, lat_col, paste0(names(models)[j], ': new-site mean prediction')))) prediction_plot_count <- prediction_plot_count + 1L }",
      "    site_ids <- rownames(pm); if (is.null(site_ids) || length(site_ids) != nrow(pm)) site_ids <- paste0('newsite_', seq_len(nrow(pm)))",
      "    richness <- rowSums(pm, na.rm=TRUE); write.csv(data.frame(site_id=site_ids, predicted_species_richness=as.numeric(richness), stringsAsFactors=FALSE), file.path(predDir, paste0('newdata_predicted_species_richness_', names(models)[j], '.csv')), row.names=FALSE)",
      "    focal_values <- NULL; focal_species <- 'Corvus_monedula'; if (!focal_species %in% colnames(pm)) { ssd <- suppressWarnings(apply(pm, 2, stats::sd, na.rm=TRUE)); ssd[!is.finite(ssd)] <- -Inf; focal_species <- colnames(pm)[which.max(ssd)] }",
      "    if (!is.null(focal_species) && length(focal_species)==1L && !is.na(focal_species) && nzchar(focal_species)) { focal_values <- pm[, focal_species]; write.csv(data.frame(site_id=site_ids, species=focal_species, predicted_probability=as.numeric(focal_values), stringsAsFactors=FALSE), file.path(predDir, paste0('newdata_predicted_focal_species_', names(models)[j], '.csv')), row.names=FALSE) }",
      "    if (file.exists(file.path('data', 'newcoords.csv'))) { coord_new <- read.csv(file.path('data', 'newcoords.csv'), check.names=FALSE, stringsAsFactors=FALSE); if (isTRUE(plot_newdata_score(richness, coord_new, site_ids, lon_col, lat_col, paste0(names(models)[j], ': new-site predicted species richness'), 'predicted species richness'))) prediction_plot_count <- prediction_plot_count + 1L; if (exists('focal_values') && isTRUE(plot_newdata_score(focal_values, coord_new, site_ids, lon_col, lat_col, paste0(names(models)[j], ': new-site focal species prediction'), paste0('predicted probability: ', focal_species)))) prediction_plot_count <- prediction_plot_count + 1L }",
      "    trait_file <- file.path('data', 'TrData.csv')",
      "    if (file.exists(trait_file)) { tdat <- read.csv(trait_file, row.names=1, check.names=FALSE, stringsAsFactors=FALSE); if (all(colnames(pm) %in% rownames(tdat))) tdat <- tdat[colnames(pm),,drop=FALSE] else if (nrow(tdat) == ncol(pm)) rownames(tdat) <- colnames(pm) else tdat <- NULL; if (!is.null(tdat)) { numeric_traits <- vapply(tdat, is.numeric, logical(1)); if (any(numeric_traits)) { tm <- as.matrix(tdat[, numeric_traits, drop=FALSE]); storage.mode(tm) <- 'numeric'; cwm <- sweep(as.matrix(pm) %*% tm, 1, pmax(rowSums(pm, na.rm=TRUE), .Machine$double.eps), '/'); rownames(cwm) <- site_ids; write.csv(cwm, file.path(predDir, paste0('newdata_community_weighted_mean_traits_', names(models)[j], '.csv'))); if (file.exists(file.path('data', 'newcoords.csv'))) { coord_new <- read.csv(file.path('data', 'newcoords.csv'), check.names=FALSE, stringsAsFactors=FALSE); trait_nm <- colnames(cwm)[1]; if (isTRUE(plot_newdata_score(cwm[, trait_nm], coord_new, site_ids, lon_col, lat_col, paste0(names(models)[j], ': new-site CWM trait prediction'), paste0('community weighted mean: ', trait_nm)))) prediction_plot_count <- prediction_plot_count + 1L } } } }",
      "    grid <- expand.grid(site_id=site_ids, response_id=colnames(pm), stringsAsFactors=FALSE)",
      "    new_rows[[length(new_rows)+1]] <- data.frame(engine='Hmsc', site_id=grid$site_id, response_id=grid$response_id, observed=NA_real_, predicted_mean=as.numeric(as.vector(pm)), predicted_lower=NA_real_, predicted_upper=NA_real_, prediction_context='newdata', stringsAsFactors=FALSE)",
      "  }",
      "  write.csv(if (length(new_rows)>0) do.call(rbind, new_rows) else data.frame(status='no_newdata_predictions_saved'), file.path('standard', 'newdata_predictions_long.csv'), row.names=FALSE)",
      "}",
      "if (length(list.files(predDir, recursive=FALSE)) == 0) write.csv(data.frame(status='no_prediction_files', reason=ifelse(do_gradients, 'no gradient predictions produced in this quick-test run', 'prediction and gradient outputs disabled by config')), file.path(predDir, 'prediction_manifest.csv'), row.names=FALSE)",
      "if (prediction_plot_count == 0L) { plot.new(); title(main='No fitted, gradient or newdata prediction plots were produced') }",
      "dev.off()"
    )
    writeLines(s7, file.path(outdir, "workflow_scripts", "S7_make_predictions.R"))

    s8 <- c(
      "# S8_standardize_outputs.R -- generated by JSDMWorkbench",
      "library(Hmsc)",
      "modelDir <- file.path('.', 'models'); resultDir <- file.path('.', 'results'); predDir <- file.path('.', 'predictions')",
      "standardDir <- file.path('.', 'standard'); plotDir <- file.path('.', 'plots'); reportDir <- file.path('.', 'report')",
      "dir.create(standardDir, recursive = TRUE, showWarnings = FALSE); dir.create(plotDir, recursive = TRUE, showWarnings = FALSE); dir.create(reportDir, recursive = TRUE, showWarnings = FALSE)",
      paste0("samples <- ", samples), paste0("thin <- ", thin), paste0("nChains <- ", nChains), paste0("nfolds <- ", nfolds),
      paste0("do_params <- ", if (isTRUE(cfg$outputs$parameters %||% TRUE)) "TRUE" else "FALSE"),
      "load(file.path(modelDir, paste0('models_thin_', thin, '_samples_', samples, '_chains_', nChains, '.Rdata')))",
      "mf_file <- file.path(modelDir, paste0('MF_thin_', thin, '_samples_', samples, '_chains_', nChains, '_nfolds_', nfolds, '.Rdata')); if (file.exists(mf_file)) load(mf_file)",
      "Y <- as.matrix(read.csv(file.path('data', 'Y.csv'), row.names = 1, check.names = FALSE)); storage.mode(Y) <- 'numeric'",
      "XData <- read.csv(file.path('data', 'XData.csv'), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)",
      "fit_rows <- list(); if (exists('MF')) for (i in seq_along(MF)) { x <- MF[[i]]; fit_rows[[length(fit_rows)+1]] <- data.frame(engine='Hmsc', metric=c('mean_TjurR2','mean_R2','mean_AUC','mean_SR2'), response_id=NA_character_, value=c(if (!is.null(x$TjurR2)) mean(x$TjurR2, na.rm=TRUE) else NA, if (!is.null(x$R2)) mean(x$R2, na.rm=TRUE) else NA, if (!is.null(x$AUC)) mean(x$AUC, na.rm=TRUE) else NA, if (!is.null(x$SR2)) mean(x$SR2, na.rm=TRUE) else NA), notes=names(MF)[i], stringsAsFactors=FALSE) }",
      "write.csv(if (length(fit_rows)>0) do.call(rbind, fit_rows) else data.frame(engine='Hmsc', metric=NA_character_, response_id=NA_character_, value=NA_real_, notes='fit not computed'), file.path(standardDir, 'fit_metrics.csv'), row.names=FALSE)",
      "effect_rows <- list(); if (do_params) for (j in seq_along(models)) { m <- models[[j]]; postBeta <- tryCatch(getPostEstimate(m, parName='Beta'), error=function(e) NULL); if (!is.null(postBeta)) { mat <- postBeta$mean; predictors <- if (!is.null(m$covNames)) m$covNames else paste0('Beta_', seq_len(nrow(mat))); species <- if (!is.null(colnames(mat))) colnames(mat) else m$spNames; grid <- expand.grid(predictor=predictors, response_id=species, stringsAsFactors=FALSE); effect_rows[[length(effect_rows)+1]] <- data.frame(engine='Hmsc', response_id=grid$response_id, predictor=grid$predictor, direction=ifelse(as.vector(mat)>=0,'positive','negative'), estimate=as.numeric(as.vector(mat)), lower=NA_real_, upper=NA_real_, notes=paste0('model=', names(models)[j], '; parameter=Beta'), stringsAsFactors=FALSE) } }",
      "write.csv(if (length(effect_rows)>0) do.call(rbind, effect_rows) else data.frame(engine='Hmsc', response_id=NA_character_, predictor=NA_character_, direction=NA_character_, estimate=NA_real_, lower=NA_real_, upper=NA_real_, notes='no effects available'), file.path(standardDir, 'effects_long.csv'), row.names=FALSE)",
      "pred_rows <- list(); for (j in seq_along(models)) { pf <- file.path(predDir, paste0('predicted_values_', names(models)[j], '.rds')); if (!file.exists(pf)) next; arr <- readRDS(pf); if (length(dim(arr)) == 3) { pm <- apply(arr, c(1,2), mean, na.rm=TRUE); pl <- apply(arr, c(1,2), quantile, probs=0.025, na.rm=TRUE); pu <- apply(arr, c(1,2), quantile, probs=0.975, na.rm=TRUE) } else { pm <- as.matrix(arr); pl <- pm; pu <- pm }; sites <- if (!is.null(rownames(pm))) rownames(pm) else rownames(Y); species <- if (!is.null(colnames(pm))) colnames(pm) else colnames(Y); grid <- expand.grid(site_id=sites, response_id=species, stringsAsFactors=FALSE); obs <- Y[cbind(match(grid$site_id, rownames(Y)), match(grid$response_id, colnames(Y)))]; pred_rows[[length(pred_rows)+1]] <- data.frame(engine='Hmsc', site_id=grid$site_id, response_id=grid$response_id, observed=as.numeric(obs), predicted_mean=as.numeric(as.vector(pm)), predicted_lower=as.numeric(as.vector(pl)), predicted_upper=as.numeric(as.vector(pu)), prediction_context='training', stringsAsFactors=FALSE) }",
      "npf <- file.path(standardDir, 'newdata_predictions_long.csv'); if (file.exists(npf)) { np <- tryCatch(read.csv(npf, check.names=FALSE, stringsAsFactors=FALSE), error=function(e) NULL); if (!is.null(np) && all(c('site_id','response_id','predicted_mean') %in% names(np))) pred_rows[[length(pred_rows)+1]] <- np }",
      "write.csv(if (length(pred_rows)>0) do.call(rbind, pred_rows) else data.frame(engine='Hmsc', site_id=NA_character_, response_id=NA_character_, observed=NA_real_, predicted_mean=NA_real_, predicted_lower=NA_real_, predicted_upper=NA_real_, prediction_context=NA_character_), file.path(standardDir, 'predictions_long.csv'), row.names=FALSE)",
      "assoc_rows <- list(); for (j in seq_along(models)) { omega <- tryCatch(computeAssociations(models[[j]]), error=function(e) NULL); if (is.null(omega)) next; for (r in seq_along(omega)) { mat <- omega[[r]]$mean; if (is.null(mat) || nrow(mat) < 2) next; idx <- which(upper.tri(mat), arr.ind=TRUE); assoc_rows[[length(assoc_rows)+1]] <- data.frame(engine='Hmsc', response_1=rownames(mat)[idx[,1]], response_2=colnames(mat)[idx[,2]], association_type='Omega', estimate=as.numeric(mat[idx]), comparable_level=paste0(names(models)[j], ':', names(omega)[r]), stringsAsFactors=FALSE) } }",
      "write.csv(if (length(assoc_rows)>0) do.call(rbind, assoc_rows) else data.frame(engine='Hmsc', response_1=NA_character_, response_2=NA_character_, association_type='Omega', estimate=NA_real_, comparable_level='no_random_level_or_association_disabled'), file.path(standardDir, 'associations_long.csv'), row.names=FALSE)",
      "plot_sources <- c(MCMC_convergence=file.path(resultDir,'MCMC_convergence.pdf'), model_fit=file.path(resultDir, paste0('model_fit_nfolds_', nfolds, '.pdf')), parameter_estimates=file.path(resultDir,'parameter_estimates.pdf'), predictions=file.path(resultDir,'predictions.pdf'))",
      "plot_manifest <- data.frame(plot=names(plot_sources), source=unname(plot_sources), copied_to=NA_character_, stringsAsFactors=FALSE); for (i in seq_along(plot_sources)) if (file.exists(plot_sources[i])) { dest <- file.path(plotDir, paste0(names(plot_sources)[i], '.pdf')); file.copy(plot_sources[i], dest, overwrite=TRUE); plot_manifest$copied_to[i] <- dest }; write.csv(plot_manifest, file.path(plotDir, 'plot_manifest.csv'), row.names=FALSE)",
      "write.csv(data.frame(run_id=basename(getwd()), engine='Hmsc', status='fitted', n_sites=nrow(Y), n_responses=ncol(Y), n_predictors=ncol(XData), samples=samples, thin=thin, nChains=nChains, stringsAsFactors=FALSE), file.path(standardDir, 'run_summary.csv'), row.names=FALSE)",
      "writeLines(c('<!DOCTYPE html><html><head><meta charset=\"UTF-8\"><title>Hmsc Report</title></head><body><h1>Hmsc workflow report</h1><p>Status: fitted</p><p>Generated by the executable S1-S7 scripts.</p></body></html>'), file.path(reportDir, 'Hmsc_report.html'))",
      "write.csv(data.frame(file=list.files('.', recursive=TRUE), stringsAsFactors=FALSE), file.path(standardDir, 'output_manifest.csv'), row.names=FALSE)"
    )
    writeLines(s8, file.path(outdir, "workflow_scripts", "S8_standardize_outputs.R"))

    writeLines(c(
      "# run_all_HMSC_S1_to_S7.R -- generated by JSDMWorkbench",
      "source('workflow_scripts/S1_define_models.R')",
      "source('workflow_scripts/S2_fit_models.R')",
      "source('workflow_scripts/S3_evaluate_convergence.R')",
      "source('workflow_scripts/S4_compute_model_fit.R')",
      "source('workflow_scripts/S5_show_model_fit.R')",
      "source('workflow_scripts/S6_show_parameter_estimates.R')",
      "source('workflow_scripts/S7_make_predictions.R')",
      "source('workflow_scripts/S8_standardize_outputs.R')"
    ), file.path(outdir, "reproducible_script", "run_all_HMSC_S1_to_S7.R"))
    writeLines(c(
      "# run_this_HMSC_analysis.R -- generated by JSDMWorkbench",
      "args <- commandArgs(trailingOnly = FALSE)",
      "file_arg <- '--file='",
      "script_arg <- args[startsWith(args, file_arg)]",
      "if (length(script_arg) > 0) {",
      "  script_path <- sub(file_arg, '', script_arg[[1]])",
      "  if (!file.exists(script_path)) script_path <- gsub('~+~', ' ', script_path, fixed = TRUE)",
      "  setwd(dirname(dirname(normalizePath(script_path, winslash = '/', mustWork = TRUE))))",
      "}",
      "source(file.path('reproducible_script', 'run_all_HMSC_S1_to_S7.R'))"
    ), file.path(outdir, "reproducible_script", "run_this_HMSC_analysis.R"))
  }

  generate_hmsc_s1s7_scripts(outdir, cfg)

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(outdir)

  tryCatch({
    library(Hmsc)
    set.seed(as.integer(cfg$model$seed %||% 1234))

    # S1 define models --------------------------------------------------------
    result$status <- "running_S1"
    safe_log("HMSC S1: defining model")
    Ymat <- as.matrix(Y)
    suppressWarnings(storage.mode(Ymat) <- "numeric")
    XData <- hmsc_clean_dataframe_types(XData, role = "XData")
    XFormula <- as.formula(cfg$model$XFormula %||% "~ .")
    distr <- cfg$model$distr %||% "probit"
    model_args <- list(Y = Ymat, XData = XData, XFormula = XFormula, distr = distr)
    model_args$XScale <- isTRUE(cfg$model$XScale)
    model_args$YScale <- isTRUE(cfg$model$YScale)
    model_args$truncateNumberOfFactors <- isTRUE(cfg$model$truncateNumberOfFactors)
    phylo_summary <- data.frame(enabled = isTRUE(cfg$model$use_phylogeny), source = NA_character_, file = NA_character_, mode = NA_character_, stringsAsFactors = FALSE)
    Loff <- hmsc_read_optional_matrix(cfg$model$Loff_file %||% "", expected_rows = nrow(Ymat), expected_cols = ncol(Ymat), label = "Loff offset")
    if (!is.null(Loff)) model_args$Loff <- Loff
    if (isTRUE(cfg$model$use_XRRR)) {
      XRRRData <- hmsc_read_optional_dataframe(cfg$model$XRRR_file %||% "", expected_rows = nrow(Ymat), label = "XRRRData")
      if (is.null(XRRRData)) stop("Use XRRR is enabled, but XRRRData file was not found in data/ or inputs/.", call. = FALSE)
      model_args$XRRRData <- XRRRData
      model_args$XRRRFormula <- as.formula(cfg$model$XRRRFormula %||% "~ .")
      model_args$XRRRScale <- isTRUE(cfg$model$XRRRScale)
      model_args$ncRRR <- as.integer(cfg$model$ncRRR %||% 2)
    }
    if (isTRUE(cfg$model$use_phylogeny)) {
      Cmat <- hmsc_read_optional_matrix(cfg$model$C_file %||% "", expected_rows = ncol(Ymat), expected_cols = ncol(Ymat), label = "deprecated internal C matrix")
      if (!is.null(Cmat)) {
        if (!is.null(rownames(Cmat)) && !is.null(colnames(Cmat)) && all(colnames(Ymat) %in% rownames(Cmat)) && all(colnames(Ymat) %in% colnames(Cmat))) {
          Cmat <- Cmat[colnames(Ymat), colnames(Ymat), drop = FALSE]
        } else {
          rownames(Cmat) <- colnames(Ymat)
          colnames(Cmat) <- colnames(Ymat)
        }
        Cmat <- hmsc_make_positive_definite_correlation(Cmat)
        model_args$C <- Cmat
        phylo_summary <- data.frame(enabled = TRUE, source = "deprecated_internal_C_matrix", file = cfg$model$C_file %||% "", mode = "C", stringsAsFactors = FALSE)
        hmsc_write_derived_phylogeny_outputs(outdir, species_names = colnames(Ymat), Cmat = Cmat,
                                             source = "deprecated_internal_C_matrix", input_file = cfg$model$C_file %||% "")
      } else {
        phy <- hmsc_read_phylogeny_or_taxonomy(cfg$data$phylogeny %||% "", colnames(Ymat))
        if (is.null(phy)) {
          stop("Use phylogeny/taxonomy is enabled, but no Newick tree or taxonomy rank table was found.", call. = FALSE)
        }
        if (!is.null(phy$C)) {
          model_args$C <- phy$C
          phylo_summary <- data.frame(enabled = TRUE, source = phy$source, file = phy$path, mode = "C", stringsAsFactors = FALSE)
        } else {
          model_args$phyloTree <- phy$phyloTree
          phylo_summary <- data.frame(enabled = TRUE, source = phy$source, file = phy$path, mode = "phyloTree", stringsAsFactors = FALSE)
        }
        hmsc_write_derived_phylogeny_outputs(outdir, phy = phy, species_names = colnames(Ymat),
                                             source = phy$source, input_file = phy$path)
      }
    }

    if (isTRUE(cfg$model$use_traits) && !is.null(TrData)) {
      TrData <- hmsc_clean_dataframe_types(TrData, role = "TrData")
      model_args$TrData <- TrData
      model_args$TrFormula <- as.formula(cfg$model$TrFormula %||% "~ .")
      model_args$TrScale <- isTRUE(cfg$model$TrScale)
    }

    random_level <- hmsc_build_random_level(Ymat, studyDesign = studyDesign, coord = coord, cfg = cfg)
    result$warnings <- c(result$warnings, random_level$warnings)
    random_mode <- random_level$mode_info$mode
    if (!is.null(random_level$studyDesign) && !is.null(random_level$ranLevels)) {
      model_args$studyDesign <- random_level$studyDesign
      model_args$ranLevels <- random_level$ranLevels
    }

    m <- do.call(Hmsc, model_args)
    models <- list(main = m)
    save(models, file = file.path("models", "unfitted_models.RData"))
    write.csv(data.frame(model = names(models), distr = distr, XFormula = deparse(XFormula), random_mode = random_mode, n_sites = nrow(Ymat), n_species = ncol(Ymat)),
              file.path("tables", "S1_defined_models.csv"), row.names = FALSE)
    write.csv(random_level$summary, file.path("tables", "S1_random_level_summary.csv"), row.names = FALSE)
    write.csv(phylo_summary, file.path("tables", "S1_phylogeny_summary.csv"), row.names = FALSE)
    parameter_audit <- do.call(rbind, hmsc_flatten_config(cfg))
    parameter_audit$status <- "read_by_server"
    parameter_audit$status[grepl("^mcmc.updater.(Beta|Gamma|Omega)$", parameter_audit$parameter)] <- "documented_Hmsc_sampleMcmc_has_no_direct_argument"
    write.csv(parameter_audit, file.path("tables", "HMSC_parameter_audit.csv"), row.names = FALSE)

    # S2 fit models -----------------------------------------------------------
    result$status <- "running_S2"
    safe_log("HMSC S2: fitting model")
    samples <- as.integer(cfg$mcmc$samples %||% 50)
    thin <- as.integer(cfg$mcmc$thin %||% 1)
    nChains <- as.integer(cfg$mcmc$nChains %||% 2)
    nParallel <- as.integer(cfg$mcmc$nParallel %||% 1)
    transient <- as.integer(cfg$mcmc$transient %||% max(1, ceiling(0.5 * samples * thin)))
    verbose <- as.integer(cfg$mcmc$verbose %||% 10)
    initPar <- hmsc_normalize_initPar(cfg$mcmc$initPar %||% "fixed effects")
    alignPost <- isTRUE(cfg$mcmc$alignPost %||% TRUE)
    if (identical(random_mode, "spatial_gpp") && isTRUE(alignPost)) {
      alignPost <- FALSE
      result$scientific_caveats <- c(result$scientific_caveats %||% character(), "Hmsc alignPost was disabled for GPP spatial fitting because alignPosterior can fail with NA latent-factor summaries in quick/small posterior runs.")
    }
    fromPrior <- isTRUE(cfg$mcmc$sample_prior %||% FALSE)
    updater <- hmsc_build_mcmc_updater(cfg)
    set.seed(as.integer(cfg$model$seed %||% 1234L))
    fitted <- list()
    for (mi in seq_along(models)) {
      safe_log("HMSC S2: sampleMcmc for model", names(models)[mi])
      fitted[[mi]] <- sampleMcmc(models[[mi]], samples = samples, thin = thin, transient = transient,
                                 initPar = initPar, verbose = verbose,
                                 adaptNf = rep(max(1L, ceiling(0.4 * samples * thin)), models[[mi]]$nr),
                                 nChains = nChains, nParallel = nParallel, updater = updater,
                                 fromPrior = fromPrior, alignPost = alignPost, engine = "R")
      if (isTRUE(cfg$mcmc$pool_chains %||% FALSE) && exists("poolMcmcChains", mode = "function")) {
        pooled <- tryCatch(poolMcmcChains(fitted[[mi]]$postList), error = function(e) {
          result$warnings <<- c(result$warnings, paste0("pool_chains failed for ", names(models)[mi], ": ", conditionMessage(e)))
          NULL
        })
        if (!is.null(pooled)) fitted[[mi]]$postList <- pooled
      }
    }
    names(fitted) <- names(models)
    models <- fitted
    model_file <- file.path("models", paste0("models_thin_", thin, "_samples_", samples, "_chains_", nChains, ".Rdata"))
    save(models, file = model_file)
    if (isTRUE(cfg$outputs$save_model %||% TRUE)) {
      saveRDS(models[[1]], file.path("models", "hmsc_model_main.rds"))
    }
    write.csv(data.frame(model_file = model_file, samples = samples, thin = thin, transient = transient, nChains = nChains, nParallel = nParallel,
                         initPar = cfg$mcmc$initPar %||% NA_character_, alignPost = alignPost, fromPrior = fromPrior,
                         pool_chains = isTRUE(cfg$mcmc$pool_chains %||% FALSE),
                         updater_GammaEta = updater$GammaEta, updater_Gamma2 = updater$Gamma2),
              file.path("tables", "S2_fit_models.csv"), row.names = FALSE)

    # S3 convergence ----------------------------------------------------------
    result$status <- "running_S3"
    safe_log("HMSC S3: convergence diagnostics")
    conv_txt <- file.path("results", "MCMC_convergence.txt")
    cat("MCMC convergence diagnostics\n\n", file = conv_txt)
    pdf(file.path("results", "MCMC_convergence.pdf"), width = 8, height = 6)
    do_diag <- isTRUE(cfg$outputs$diagnostics %||% TRUE)
    extract_mcmc_lists <- function(x) {
      out <- list()
      walk <- function(z, nm) {
        if (inherits(z, "mcmc.list")) {
          out[[nm]] <<- z
        } else if (is.list(z)) {
          for (ii in seq_along(z)) walk(z[[ii]], paste0(nm, "_", ii))
        }
      }
      walk(x, "value")
      out
    }
    if (!do_diag) {
      plot.new(); title(main = "MCMC diagnostics disabled"); text(0.5, 0.5, "MCMC diagnostics disabled by output settings.")
      cat("MCMC diagnostics disabled by output settings.\n", file = conv_txt, append = TRUE)
    } else {
      conv <- cfg$outputs$convergence %||% list()
      targets <- c(Beta = isTRUE(conv$showBeta %||% TRUE),
                   Gamma = isTRUE(conv$showGamma %||% TRUE),
                   Omega = isTRUE(conv$showOmega %||% TRUE),
                   Rho = isTRUE(conv$showRho %||% TRUE),
                   Alpha = isTRUE(conv$showAlpha %||% TRUE))
      for (j in seq_along(models)) {
        cat("\nModel: ", names(models)[j], "\n", file = conv_txt, append = TRUE)
        mpost <- tryCatch(convertToCodaObject(models[[j]], spNamesNumbers = c(TRUE, FALSE), covNamesNumbers = c(TRUE, FALSE)), error = function(e) {
          result$warnings <<- c(result$warnings, paste0("convertToCodaObject failed for ", names(models)[j], ": ", conditionMessage(e)))
          NULL
        })
        if (is.null(mpost) || !requireNamespace("coda", quietly = TRUE)) {
          cat("Coda diagnostics unavailable for this model.\n", file = conv_txt, append = TRUE)
          next
        }
        for (par_name in names(targets)[targets]) {
          if (is.null(mpost[[par_name]])) {
            cat(par_name, "not available.\n", file = conv_txt, append = TRUE)
            next
          }
          lists <- extract_mcmc_lists(mpost[[par_name]])
          if (length(lists) == 0) {
            cat(par_name, "has no mcmc.list entries.\n", file = conv_txt, append = TRUE)
            next
          }
          for (lnm in names(lists)) {
            obj <- lists[[lnm]]
            if (isTRUE(conv$gelmanPSRF %||% TRUE) && nChains > 1) {
              psrf <- tryCatch(coda::gelman.diag(obj, multivariate = FALSE)$psrf[,1], error = function(e) NULL)
              if (!is.null(psrf)) {
                cat("\n", par_name, lnm, "PSRF summary:\n", file = conv_txt, append = TRUE)
                cat(capture.output(summary(psrf)), sep = "\n", file = conv_txt, append = TRUE)
                psrf_finite <- psrf[is.finite(psrf)]
                if (length(psrf_finite) > 0) {
                  hist(psrf_finite, main = paste("PSRF", par_name, names(models)[j]), xlab = "PSRF")
                } else {
                  plot.new(); title(main = paste("PSRF", par_name, names(models)[j])); text(0.5, 0.5, "No finite PSRF values available.")
                }
                psrf_names <- names(psrf)
                if (is.null(psrf_names) || length(psrf_names) != length(psrf)) psrf_names <- paste0("par_", seq_along(psrf))
                write.csv(data.frame(parameter = psrf_names, PSRF = as.numeric(psrf)),
                          file.path("tables", paste0("S3_", par_name, "_PSRF_", lnm, "_", names(models)[j], ".csv")), row.names = FALSE)
              }
            }
            if (isTRUE(conv$effectiveSize %||% TRUE)) {
              ess <- tryCatch(coda::effectiveSize(obj), error = function(e) NULL)
              if (!is.null(ess)) {
                ess_names <- names(ess)
                if (is.null(ess_names) || length(ess_names) != length(ess)) ess_names <- paste0("par_", seq_along(ess))
                write.csv(data.frame(parameter = ess_names, effectiveSize = as.numeric(ess)),
                          file.path("tables", paste0("S3_", par_name, "_effectiveSize_", lnm, "_", names(models)[j], ".csv")), row.names = FALSE)
              }
            }
          }
        }
      }
    }
    dev.off()

    # S4 compute model fit ----------------------------------------------------
    result$status <- "running_S4"
    safe_log("HMSC S4: computing model fit")
    MF <- list(); MFCV <- list(); WAIC <- list(); PRED <- list()
    nfolds <- as.integer(cfg$outputs$predictions$nfolds %||% cfg$outputs$nfolds %||% 2)
    do_predicted <- isTRUE(cfg$outputs$predicted %||% TRUE)
    do_fit <- isTRUE(cfg$outputs$fit %||% TRUE)
    do_cv <- isTRUE(cfg$outputs$cv %||% TRUE)
    do_waic <- isTRUE(cfg$outputs$waic %||% FALSE)
    for (mi in seq_along(models)) {
      m <- models[[mi]]
      PRED[mi] <- list(NULL)
      MF[mi] <- list(NULL)
      MFCV[mi] <- list(NULL)
      WAIC[mi] <- list("disabled")
      preds <- NULL
      if (do_predicted || do_fit || do_cv) {
        preds <- tryCatch(computePredictedValues(m, nParallel = nParallel), error = function(e) {
          result$warnings <<- c(result$warnings, paste0("computePredictedValues failed for ", names(models)[mi], ": ", conditionMessage(e)))
          NULL
        })
        if (!is.null(preds)) {
          PRED[mi] <- list(preds)
          if (do_predicted) saveRDS(preds, file.path("predictions", paste0("predicted_values_", names(models)[mi], ".rds")))
        } else {
          PRED[mi] <- list(NULL)
        }
      }
      if (do_fit && !is.null(preds)) {
        MF[mi] <- list(tryCatch(evaluateModelFit(hM = m, predY = preds), error = function(e) {
          msg <- conditionMessage(e)
          write.csv(data.frame(model = names(models)[mi], status = "skipped", reason = msg, stringsAsFactors = FALSE),
                    file.path("diagnostics", paste0("evaluateModelFit_skipped_", names(models)[mi], ".csv")), row.names = FALSE)
          if (grepl("No control observation", msg, ignore.case = TRUE)) {
            result$scientific_caveats <<- c(result$scientific_caveats %||% character(), paste0("evaluateModelFit skipped for ", names(models)[mi], " because the quick/small test data have no control observation; fitted predictions and standard prediction tables were still exported."))
          } else {
            result$warnings <<- c(result$warnings, paste0("evaluateModelFit failed for ", names(models)[mi], ": ", msg))
          }
          list(fit_error = msg)
        }))
      } else {
        MF[mi] <- list(if (do_fit) list(fit_error = "computePredictedValues failed or predictions disabled") else NULL)
      }
      WAIC[mi] <- list(if (do_waic) tryCatch(computeWAIC(m), error = function(e) e$message) else "disabled")
      if (do_cv && !is.null(preds)) {
        partition_col <- trimws(as.character(cfg$outputs$predictions$partition_column %||% ""))
        partition <- tryCatch({
          args <- list(hM = m, nfolds = nfolds)
          if (nzchar(partition_col) && !is.null(m$studyDesign) && partition_col %in% names(m$studyDesign)) args$column <- partition_col
          do.call(createPartition, args)
        }, error = function(e) {
          result$warnings <<- c(result$warnings, paste0("createPartition failed for ", names(models)[mi], ": ", conditionMessage(e)))
          NULL
        })
        if (!is.null(partition)) {
          preds_cv <- tryCatch(computePredictedValues(m, partition = partition, nParallel = nParallel), error = function(e) {
            result$warnings <<- c(result$warnings, paste0("cross-validation computePredictedValues failed for ", names(models)[mi], ": ", conditionMessage(e)))
            NULL
          })
          MFCV[mi] <- list(if (!is.null(preds_cv)) {
            tryCatch(evaluateModelFit(hM = m, predY = preds_cv), error = function(e) {
              msg <- conditionMessage(e)
              write.csv(data.frame(model = names(models)[mi], status = "skipped", reason = msg, stringsAsFactors = FALSE),
                        file.path("diagnostics", paste0("cross_validation_evaluateModelFit_skipped_", names(models)[mi], ".csv")), row.names = FALSE)
              if (grepl("No control observation", msg, ignore.case = TRUE)) {
                result$scientific_caveats <<- c(result$scientific_caveats %||% character(), paste0("cross-validation evaluateModelFit skipped for ", names(models)[mi], " because the quick/small partition has no control observation."))
              } else {
                result$warnings <<- c(result$warnings, paste0("cross-validation evaluateModelFit failed for ", names(models)[mi], ": ", msg))
              }
              list(fit_error = msg)
            })
          } else NULL)
        } else {
          MFCV[mi] <- list(NULL)
        }
      } else {
        MFCV[mi] <- list(NULL)
      }
    }
    names(MF) <- names(models); names(MFCV) <- names(models); names(WAIC) <- names(models); names(PRED) <- names(models)
    mf_file <- file.path("models", paste0("MF_thin_", thin, "_samples_", samples, "_chains_", nChains, "_nfolds_", nfolds, ".Rdata"))
    save(MF, MFCV, WAIC, file = mf_file)
    mf_rows <- do.call(rbind, lapply(seq_along(MF), function(i) {
      x <- MF[[i]]
      data.frame(model = names(MF)[i],
                 mean_TjurR2 = if (!is.null(x$TjurR2)) mean(x$TjurR2, na.rm = TRUE) else NA,
                 mean_R2 = if (!is.null(x$R2)) mean(x$R2, na.rm = TRUE) else NA,
                 mean_AUC = if (!is.null(x$AUC)) mean(x$AUC, na.rm = TRUE) else NA,
                 mean_SR2 = if (!is.null(x$SR2)) mean(x$SR2, na.rm = TRUE) else NA)
    }))
    write.csv(mf_rows, file.path("tables", "S4_model_fit_summary.csv"), row.names = FALSE)

    # S5 show model fit -------------------------------------------------------
    result$status <- "running_S5"
    safe_log("HMSC S5: plotting model fit")
    pdf(file.path("results", paste0("model_fit_nfolds_", nfolds, ".pdf")), width = 8, height = 6)
    for (j in seq_along(MF)) {
      cMF <- MF[[j]]; cMFCV <- MFCV[[j]]
      if (!is.null(cMF$TjurR2) && !is.null(cMFCV$TjurR2)) { plot(cMF$TjurR2, cMFCV$TjurR2, xlim=c(-1,1), ylim=c(-1,1), xlab="explanatory power", ylab="predictive power", main=paste0(names(MF)[j], ": Tjur R2")); abline(0,1) }
      if (!is.null(cMF$R2) && !is.null(cMFCV$R2)) { plot(cMF$R2, cMFCV$R2, xlim=c(-1,1), ylim=c(-1,1), xlab="explanatory power", ylab="predictive power", main=paste0(names(MF)[j], ": R2")); abline(0,1) }
      if (!is.null(cMF$AUC) && !is.null(cMFCV$AUC)) { plot(cMF$AUC, cMFCV$AUC, xlim=c(0,1), ylim=c(0,1), xlab="explanatory power", ylab="predictive power", main=paste0(names(MF)[j], ": AUC")); abline(0,1); abline(v=.5); abline(h=.5) }
      if (!is.null(cMF$SR2) && !is.null(cMFCV$SR2)) { plot(cMF$SR2, cMFCV$SR2, xlim=c(-1,1), ylim=c(-1,1), xlab="explanatory power", ylab="predictive power", main=paste0(names(MF)[j], ": SR2")); abline(0,1) }
    }
    dev.off()

    # S6 parameter estimates --------------------------------------------------
    result$status <- "running_S6"
    safe_log("HMSC S6: parameter estimates")
    text.file <- file.path("results", "parameter_estimates.txt")
    cat("Parameter estimates\n\n", file=text.file)
    pdf(file.path("results", "parameter_estimates.pdf"), width=9, height=7)
    do_params <- isTRUE(cfg$outputs$parameters %||% TRUE)
    do_vp <- isTRUE(cfg$outputs$variance_partitioning %||% TRUE)
    do_omega <- isTRUE(cfg$outputs$omega %||% TRUE)
    for (j in seq_along(models)) {
      m <- models[[j]]
      cat("\n", names(models)[j], "\n", file=text.file, append=TRUE)
      VP <- if (do_params && do_vp) hmsc_quiet_optional(computeVariancePartitioning(m)) else NULL
      if (!is.null(VP)) {
        write.csv(VP$vals, file.path("results", paste0("parameter_estimates_VP_", names(models)[j], ".csv")))
        hmsc_quiet_optional(hmsc_draw_vp_bookstyle(VP, model_name = names(models)[j]))
        hmsc_quiet_optional(hmsc_save_vp_bookstyle(VP, model_name = names(models)[j], result_dir = "results"))
        VP_plot <- VP
        if (isTRUE(cfg$outputs$plotting$var.part.order.explained %||% TRUE) && !is.null(VP_plot$vals) && is.matrix(VP_plot$vals)) {
          ord <- order(colMeans(VP_plot$vals, na.rm = TRUE), decreasing = TRUE)
          VP_plot$vals <- VP_plot$vals[, ord, drop = FALSE]
        }
        hmsc_quiet_optional(plotVariancePartitioning(hM=m, VP=VP_plot, main=paste0("Proportion of explained variance, ", names(models)[j])))
      }
      if (do_params && isTRUE(cfg$outputs$convergence$showRho %||% TRUE)) {
        hmsc_quiet_optional(hmsc_draw_rho_signal_curve(m, model_name = names(models)[j]))
        hmsc_quiet_optional(hmsc_save_rho_signal_curve(m, model_name = names(models)[j], result_dir = "results"))
      }
      postBeta <- if (do_params) hmsc_quiet_optional(getPostEstimate(m, parName="Beta")) else NULL
      if (!is.null(postBeta)) {
        write.csv(as.data.frame(t(postBeta$mean)), file.path("results", paste0("parameter_estimates_Beta_mean_", names(models)[j], ".csv")))
        if (isTRUE(cfg$outputs$beta_support %||% TRUE)) write.csv(as.data.frame(t(postBeta$support)), file.path("results", paste0("parameter_estimates_Beta_support_", names(models)[j], ".csv")))
        beta_heatmap_data <- hmsc_prepare_beta_heatmap_data(m, postBeta, include_intercept = FALSE)
        if (!is.null(beta_heatmap_data)) {
          hmsc_write_beta_heatmap_long(
            beta_heatmap_data$beta, beta_heatmap_data$support,
            file.path("results", paste0("parameter_estimates_Beta_species_environment_long_", names(models)[j], ".csv")),
            model_name = names(models)[j]
          )
        }
        if (isTRUE(cfg$outputs$plotting$plotBeta %||% TRUE)) {
          show_sp_beta <- isTRUE(cfg$outputs$plotting$show.sp.names.beta %||% FALSE)
          hmsc_quiet_optional(hmsc_draw_phylo_beta_heatmap(
            m, postBeta, model_name = names(models)[j],
            include_intercept = FALSE,
            support_level = 0.95,
            show_species_names = TRUE
          ))
          hmsc_quiet_optional(hmsc_save_phylo_beta_heatmap(
            m, postBeta, model_name = names(models)[j], result_dir = "results",
            include_intercept = FALSE,
            support_level = 0.95,
            show_species_names = TRUE
          ))
          hmsc_quiet_optional(plotBeta(m, post=postBeta, supportLevel=0.95, param="Sign",
                                       plotTree=isTRUE(cfg$outputs$plotting$plotTree %||% FALSE),
                                       covNamesNumbers=c(TRUE,FALSE), spNamesNumbers=c(TRUE, show_sp_beta)))
        }
      }
      postGamma <- if (do_params) hmsc_quiet_optional(getPostEstimate(m, parName="Gamma")) else NULL
      if (!is.null(postGamma)) {
        write.csv(as.data.frame(t(postGamma$mean)), file.path("results", paste0("parameter_estimates_Gamma_mean_", names(models)[j], ".csv")))
        if (isTRUE(cfg$outputs$gamma_support %||% TRUE)) write.csv(as.data.frame(t(postGamma$support)), file.path("results", paste0("parameter_estimates_Gamma_support_", names(models)[j], ".csv")))
        gamma_heatmap_data <- hmsc_prepare_gamma_heatmap_data(m, postGamma, include_intercepts = TRUE)
        if (!is.null(gamma_heatmap_data)) {
          hmsc_write_gamma_heatmap_long(
            gamma_heatmap_data$gamma_all, gamma_heatmap_data$support_all,
            file.path("results", paste0("parameter_estimates_trait_environment_gamma_long_", names(models)[j], ".csv")),
            model_name = names(models)[j]
          )
        }
        hmsc_quiet_optional(hmsc_save_trait_environment_gamma_heatmap(
          m, postGamma, model_name = names(models)[j], result_dir = "results",
          include_intercepts = FALSE,
          support_level = 0.95
        ))
        if (isTRUE(cfg$outputs$plotting$plotGamma %||% TRUE)) {
          hmsc_quiet_optional(hmsc_draw_trait_environment_gamma_heatmap(
            m, postGamma, model_name = names(models)[j],
            include_intercepts = FALSE,
            support_level = 0.95
          ))
          hmsc_quiet_optional(plotGamma(m, post=postGamma, supportLevel=0.95, param="Sign", covNamesNumbers=c(TRUE,FALSE), trNamesNumbers=c(TRUE,FALSE)))
        }
      }
      OmegaCor <- if (do_params && do_omega) hmsc_quiet_optional(computeAssociations(m)) else NULL
      if (!is.null(OmegaCor)) {
        for (r in seq_along(OmegaCor)) {
          omega_level_name <- names(OmegaCor)[r]
          if (is.null(omega_level_name) || length(omega_level_name) == 0 || is.na(omega_level_name) || !nzchar(omega_level_name)) {
            omega_level_name <- paste0("random_level_", r)
          }
          write.csv(OmegaCor[[r]]$mean, file.path("results", paste0("parameter_estimates_Omega_mean_", names(models)[j], "_", omega_level_name, ".csv")))
          if (isTRUE(cfg$outputs$omega_support %||% TRUE)) write.csv(OmegaCor[[r]]$support, file.path("results", paste0("parameter_estimates_Omega_support_", names(models)[j], "_", omega_level_name, ".csv")))
          hmsc_quiet_optional(hmsc_draw_association_heatmap(
            m, OmegaCor[[r]]$mean, OmegaCor[[r]]$support,
            model_name = names(models)[j], level_name = omega_level_name,
            support_level = 0.95,
            show_species_names = isTRUE(cfg$outputs$plotting$show.sp.names.omega %||% TRUE)
          ))
          hmsc_quiet_optional(hmsc_save_association_heatmap(
            m, OmegaCor[[r]]$mean, OmegaCor[[r]]$support,
            model_name = names(models)[j], level_name = omega_level_name,
            result_dir = "results",
            support_level = 0.95,
            show_species_names = isTRUE(cfg$outputs$plotting$show.sp.names.omega %||% TRUE)
          ))
          if (requireNamespace("corrplot", quietly = TRUE)) {
            omega_mat <- OmegaCor[[r]]$mean
            max_omega <- as.integer(cfg$outputs$convergence$maxOmega %||% nrow(omega_mat))
            if (is.finite(max_omega) && max_omega > 1 && nrow(omega_mat) > max_omega) omega_mat <- omega_mat[seq_len(max_omega), seq_len(max_omega), drop = FALSE]
            omega_order <- hmsc_normalize_omega_order(cfg$outputs$plotting$omega.order %||% "original")
            tl_pos <- if (isTRUE(cfg$outputs$plotting$show.sp.names.omega %||% TRUE)) "lt" else "n"
            hmsc_quiet_optional(corrplot::corrplot(omega_mat, method="color", order=omega_order, tl.pos=tl_pos, main=paste0("Associations: ", names(models)[j])))
          }
        }
      }
    }
    dev.off()

    # S7 predictions ----------------------------------------------------------
    result$status <- "running_S7"
    safe_log("HMSC S7: predictions")
    pdf(file.path("results", "predictions.pdf"), width=9, height=7)
    prediction_plot_count <- 0L
    new_prediction_rows <- list()
    new_prediction_manifest <- data.frame(
      status = if (is.null(newdata)) "not_requested" else if (isTRUE(cfg$outputs$predicted %||% TRUE)) "pending" else "skipped",
      reason = if (is.null(newdata)) "newdata.csv was not uploaded" else if (isTRUE(cfg$outputs$predicted %||% TRUE)) "" else "computePredictedValues / prediction output is disabled",
      rows = if (is.null(newdata)) 0L else nrow(newdata),
      stringsAsFactors = FALSE
    )
    for (j in seq_along(models)) {
      m <- models[[j]]
      covariates <- tryCatch(intersect(all.vars(m$XFormula), colnames(m$XData)), error=function(e) character())
      if (length(covariates) == 0) covariates <- tryCatch(colnames(m$XData), error=function(e) character())
      env_filter <- hmsc_parse_name_list(cfg$outputs$predictions$env.list %||% "")
      if (length(env_filter) > 0) {
        keep <- intersect(covariates, env_filter)
        if (length(keep) > 0) {
          covariates <- keep
        } else {
          result$warnings <- c(result$warnings, paste0("env.list did not match any environmental covariates for ", names(models)[j], "."))
        }
      }
      write.csv(data.frame(
        model = names(models)[j],
        species_list = paste(hmsc_parse_name_list(cfg$outputs$predictions$species.list %||% ""), collapse = ","),
        trait_list = paste(hmsc_parse_name_list(cfg$outputs$predictions$trait.list %||% ""), collapse = ","),
        env_list = paste(covariates, collapse = ","),
        computeSAIR = isTRUE(cfg$outputs$predictions$computeSAIR %||% FALSE),
        stringsAsFactors = FALSE
      ), file.path("tables", paste0("S7_prediction_settings_", names(models)[j], ".csv")), row.names = FALSE)
      if (isTRUE(cfg$outputs$gradients %||% FALSE) && length(covariates) > 0) {
        for (covariate in covariates) {
          Gradient <- hmsc_quiet_optional(constructGradient(m, focalVariable=covariate))
          if (!is.null(Gradient)) {
            predY <- hmsc_quiet_optional(predict(m, Gradient=Gradient, expected=TRUE))
            if (!is.null(predY)) {
              saveRDS(predY, file.path("predictions", paste0("gradient_prediction_", names(models)[j], "_", covariate, ".rds")))
              gradient_species <- hmsc_parse_name_list(cfg$outputs$predictions$species.list %||% "")
              gradient_traits <- hmsc_parse_name_list(cfg$outputs$predictions$trait.list %||% "")
              if (isTRUE(hmsc_quiet_optional(hmsc_draw_gradient_response_panel(
                m, Gradient, predY, covariate,
                model_name = names(models)[j],
                TrData = TrData,
                species_list = gradient_species,
                trait_list = gradient_traits
              ), FALSE))) {
                prediction_plot_count <- prediction_plot_count + 1L
              }
              hmsc_quiet_optional(hmsc_save_gradient_response_panel(
                m, Gradient, predY, covariate,
                model_name = names(models)[j],
                result_dir = "predictions",
                TrData = TrData,
                species_list = gradient_species,
                trait_list = gradient_traits
              ))
              hmsc_quiet_optional(plotGradient(m, Gradient, pred=predY, yshow=0, measure="S", showData=TRUE, main=paste0(names(models)[j], ": ", covariate)))
            }
          }
        }
      } else {
        plotted_training <- FALSE
        if (!is.null(PRED[[j]])) {
          train_ps <- tryCatch(
            hmsc_summarise_prediction_object(
              PRED[[j]],
              site_ids = rownames(Ymat) %||% paste0("site_", seq_len(nrow(Ymat))),
              species_ids = colnames(Ymat) %||% paste0("sp_", seq_len(ncol(Ymat)))
            ),
            error = function(e) NULL
          )
          if (!is.null(train_ps) && !is.null(train_ps$mean)) {
            hmsc_plot_prediction_heatmap(
              train_ps$mean,
              main = paste0(names(models)[j], ": fitted prediction mean"),
              observed = Ymat,
              value_label = "fitted predicted mean"
            )
            prediction_plot_count <- prediction_plot_count + 1L
            plotted_training <- TRUE
          }
        }
        if (!plotted_training) {
          plot.new()
          title(main = paste0(names(models)[j], ": no fitted prediction plot"))
          text(0.5, 0.5, "No fitted predicted values were available. Check predictions/ and diagnostics/.")
        }
      }
    }
    if (!is.null(newdata) && isTRUE(cfg$outputs$predicted %||% TRUE)) {
      pred_prep <- tryCatch(hmsc_prepare_prediction_xdata(newdata, XData, cfg$model$XFormula %||% "~ ."), error = function(e) e)
      if (inherits(pred_prep, "error")) {
        msg <- conditionMessage(pred_prep)
        result$warnings <- c(result$warnings, paste0("HMSC newdata prediction skipped: ", msg))
        writeLines(msg, file.path("diagnostics", "HMSC_newdata_prediction_error.txt"))
        new_prediction_manifest <- data.frame(status = "skipped", reason = msg, rows = nrow(newdata), stringsAsFactors = FALSE)
      } else {
        for (j in seq_along(models)) {
          m <- models[[j]]
          pred_new <- tryCatch(hmsc_predict_new_sites(m, pred_prep$XData, newStudyDesign,
            newcoord, cfg, pred_prep$site_ids, nParallel), error = function(e) e)
          if (inherits(pred_new, "error")) {
            msg <- conditionMessage(pred_new)
            result$warnings <- c(result$warnings, paste0("HMSC newdata prediction failed for ", names(models)[j], ": ", msg))
            writeLines(msg, file.path("diagnostics", paste0("HMSC_newdata_prediction_error_", names(models)[j], ".txt")))
            next
          }
          species <- m$spNames %||% colnames(Ymat) %||% paste0("sp_", seq_len(ncol(Ymat)))
          ps <- hmsc_summarise_prediction_object(pred_new, site_ids = pred_prep$site_ids, species_ids = species)
          pred_mean <- ps$mean
          pred_lower <- ps$lower
          pred_upper <- ps$upper
          species <- colnames(pred_mean) %||% species
          write.csv(pred_mean, file.path("predictions", paste0("newdata_predicted_mean_", names(models)[j], ".csv")))
          write.csv(pred_lower, file.path("predictions", paste0("newdata_predicted_lower_", names(models)[j], ".csv")))
          write.csv(pred_upper, file.path("predictions", paste0("newdata_predicted_upper_", names(models)[j], ".csv")))
          saveRDS(pred_new, file.path("predictions", paste0("newdata_predicted_values_", names(models)[j], ".rds")))
          hmsc_plot_prediction_heatmap(
            pred_mean,
            main = paste0(names(models)[j], ": newdata predicted mean"),
            observed = NULL,
            value_label = "newdata predicted mean"
          )
          prediction_plot_count <- prediction_plot_count + 1L
          if (!is.null(newcoord)) {
            if (isTRUE(hmsc_plot_newdata_spatial_prediction(
              pred_mean, newcoord, cfg, pred_prep$site_ids,
              main = paste0(names(models)[j], ": new-site mean prediction")
            ))) {
              prediction_plot_count <- prediction_plot_count + 1L
            }
          }
          pred_site_ids <- rownames(pred_mean)
          if (is.null(pred_site_ids) || length(pred_site_ids) != nrow(pred_mean)) {
            pred_site_ids <- pred_prep$site_ids
          }
          if (is.null(pred_site_ids) || length(pred_site_ids) != nrow(pred_mean)) {
            pred_site_ids <- paste0("newsite_", seq_len(nrow(pred_mean)))
          }
          richness <- rowSums(pred_mean, na.rm = TRUE)
          write.csv(
            data.frame(site_id = pred_site_ids,
                       predicted_species_richness = as.numeric(richness),
                       stringsAsFactors = FALSE),
            file.path("predictions", paste0("newdata_predicted_species_richness_", names(models)[j], ".csv")),
            row.names = FALSE
          )
          if (!is.null(newcoord) && isTRUE(hmsc_plot_newdata_spatial_score(
            richness, newcoord, cfg, pred_site_ids,
            main = paste0(names(models)[j], ": new-site predicted species richness"),
            score_label = "predicted species richness"
          ))) {
            prediction_plot_count <- prediction_plot_count + 1L
          }
          focal_species <- "Corvus_monedula"
          if (!focal_species %in% colnames(pred_mean)) {
            species_sd <- suppressWarnings(apply(pred_mean, 2, stats::sd, na.rm = TRUE))
            species_sd[!is.finite(species_sd)] <- -Inf
            focal_species <- colnames(pred_mean)[which.max(species_sd)]
          }
          if (!is.null(focal_species) && length(focal_species) == 1L && !is.na(focal_species) && nzchar(focal_species)) {
            focal_values <- pred_mean[, focal_species]
            write.csv(
              data.frame(site_id = pred_site_ids,
                         species = focal_species,
                         predicted_probability = as.numeric(focal_values),
                         stringsAsFactors = FALSE),
              file.path("predictions", paste0("newdata_predicted_focal_species_", names(models)[j], ".csv")),
              row.names = FALSE
            )
            if (!is.null(newcoord) && isTRUE(hmsc_plot_newdata_spatial_score(
              focal_values, newcoord, cfg, pred_site_ids,
              main = paste0(names(models)[j], ": new-site focal species prediction"),
              score_label = paste0("predicted probability: ", focal_species)
            ))) {
              prediction_plot_count <- prediction_plot_count + 1L
            }
          }
          if (!is.null(TrData) && nrow(TrData) > 0 && ncol(TrData) > 0) {
            trait_data <- TrData
            if (all(colnames(pred_mean) %in% rownames(trait_data))) {
              trait_data <- trait_data[colnames(pred_mean), , drop = FALSE]
            } else if (nrow(trait_data) == ncol(pred_mean)) {
              rownames(trait_data) <- colnames(pred_mean)
            } else {
              trait_data <- NULL
            }
            if (!is.null(trait_data)) {
              numeric_traits <- vapply(trait_data, is.numeric, logical(1))
              if (any(numeric_traits)) {
                trait_matrix <- as.matrix(trait_data[, numeric_traits, drop = FALSE])
                storage.mode(trait_matrix) <- "numeric"
                denom <- pmax(rowSums(pred_mean, na.rm = TRUE), .Machine$double.eps)
                cwm <- sweep(as.matrix(pred_mean) %*% trait_matrix, 1, denom, "/")
                rownames(cwm) <- pred_site_ids
                write.csv(cwm, file.path("predictions", paste0("newdata_community_weighted_mean_traits_", names(models)[j], ".csv")))
                cwm_trait <- colnames(cwm)[1]
                if (!is.null(newcoord) && isTRUE(hmsc_plot_newdata_spatial_score(
                  cwm[, cwm_trait], newcoord, cfg, pred_site_ids,
                  main = paste0(names(models)[j], ": new-site CWM trait prediction"),
                  score_label = paste0("community weighted mean: ", cwm_trait)
                ))) {
                  prediction_plot_count <- prediction_plot_count + 1L
                }
              }
            }
          }
          grid <- expand.grid(site_id = pred_prep$site_ids, response_id = species, stringsAsFactors = FALSE)
          new_prediction_rows[[length(new_prediction_rows) + 1]] <- data.frame(
            engine = "Hmsc", site_id = grid$site_id, response_id = grid$response_id,
            observed = NA_real_, predicted_mean = as.numeric(as.vector(pred_mean)),
            predicted_lower = as.numeric(as.vector(pred_lower)), predicted_upper = as.numeric(as.vector(pred_upper)),
            prediction_context = "newdata",
            stringsAsFactors = FALSE)
          new_prediction_manifest <- data.frame(status = "saved", reason = paste(pred_prep$messages, collapse = "; "), rows = nrow(pred_prep$XData), stringsAsFactors = FALSE)
        }
      }
    }
    write.csv(new_prediction_manifest, file.path("predictions", "newdata_prediction_manifest.csv"), row.names = FALSE)
    if (prediction_plot_count == 0L) {
      plot.new()
      title(main = "HMSC predictions")
      text(0.5, 0.5, "No fitted or newdata prediction plots were generated. Check diagnostics/ for details.")
    }
    dev.off()
    if (length(list.files("predictions", recursive = FALSE)) == 0) {
      write.csv(
        data.frame(
          status = "no_prediction_files",
          reason = if (isTRUE(cfg$outputs$gradients %||% FALSE)) "no gradient predictions produced in this quick-test run" else "prediction and gradient outputs disabled by config",
          stringsAsFactors = FALSE
        ),
        file.path("predictions", "prediction_manifest.csv"),
        row.names = FALSE
      )
    }

    # Standard comparison outputs --------------------------------------------
    standard_fit <- do.call(rbind, lapply(seq_len(nrow(mf_rows)), function(i) {
      data.frame(engine = "Hmsc",
                 metric = c("mean_TjurR2", "mean_R2", "mean_AUC", "mean_SR2"),
                 response_id = NA_character_,
                 value = as.numeric(mf_rows[i, c("mean_TjurR2", "mean_R2", "mean_AUC", "mean_SR2")]),
                 notes = mf_rows$model[i],
                 stringsAsFactors = FALSE)
    }))
    write.csv(standard_fit, file.path("standard", "fit_metrics.csv"), row.names = FALSE)

    effect_rows <- list()
    if (isTRUE(cfg$outputs$parameters %||% TRUE)) {
      for (j in seq_along(models)) {
        m <- models[[j]]
        postBeta <- tryCatch(getPostEstimate(m, parName = "Beta"), error = function(e) NULL)
        if (!is.null(postBeta) && !is.null(postBeta$mean)) {
          mat <- postBeta$mean
          predictors <- m$covNames %||% rownames(mat) %||% paste0("Beta_", seq_len(nrow(mat)))
          species <- colnames(mat) %||% m$spNames %||% paste0("sp_", seq_len(ncol(mat)))
          grid <- expand.grid(predictor = predictors, response_id = species, stringsAsFactors = FALSE)
          effect_rows[[length(effect_rows) + 1]] <- data.frame(
            engine = "Hmsc", response_id = grid$response_id, predictor = grid$predictor,
            direction = ifelse(as.vector(mat) >= 0, "positive", "negative"),
            estimate = as.numeric(as.vector(mat)), lower = NA_real_, upper = NA_real_,
            notes = paste0("model=", names(models)[j], "; parameter=Beta; posterior_support=", round(as.numeric(as.vector(postBeta$support %||% matrix(NA_real_, nrow(mat), ncol(mat)))), 4)),
            stringsAsFactors = FALSE)
        }
        postGamma <- tryCatch(getPostEstimate(m, parName = "Gamma"), error = function(e) NULL)
        if (!is.null(postGamma) && !is.null(postGamma$mean)) {
          mat <- postGamma$mean
          predictors <- paste0("Gamma_", seq_len(nrow(mat)))
          traits <- colnames(mat) %||% paste0("trait_", seq_len(ncol(mat)))
          grid <- expand.grid(predictor = predictors, response_id = paste0("trait:", traits), stringsAsFactors = FALSE)
          effect_rows[[length(effect_rows) + 1]] <- data.frame(
            engine = "Hmsc", response_id = grid$response_id, predictor = grid$predictor,
            direction = ifelse(as.vector(mat) >= 0, "positive", "negative"),
            estimate = as.numeric(as.vector(mat)), lower = NA_real_, upper = NA_real_,
            notes = paste0("model=", names(models)[j], "; parameter=Gamma"),
            stringsAsFactors = FALSE)
        }
      }
    }
    effects_long <- if (length(effect_rows) > 0) do.call(rbind, effect_rows) else data.frame(engine="Hmsc", response_id=NA_character_, predictor=NA_character_, direction=NA_character_, estimate=NA_real_, lower=NA_real_, upper=NA_real_, notes="No Beta/Gamma estimates available.", stringsAsFactors=FALSE)
    write.csv(effects_long, file.path("standard", "effects_long.csv"), row.names = FALSE)

    prediction_rows <- list()
    for (j in seq_along(PRED)) {
      arr <- PRED[[j]]
      if (is.null(arr)) next
      sites <- rownames(Ymat) %||% paste0("site_", seq_len(nrow(Ymat)))
      species <- colnames(Ymat) %||% paste0("sp_", seq_len(ncol(Ymat)))
      ps <- hmsc_summarise_prediction_object(arr, site_ids = sites, species_ids = species)
      pred_mean <- ps$mean
      pred_lower <- ps$lower
      pred_upper <- ps$upper
      sites <- rownames(pred_mean) %||% sites
      species <- colnames(pred_mean) %||% species
      grid <- expand.grid(site_id = sites, response_id = species, stringsAsFactors = FALSE)
      obs_row <- match(grid$site_id, rownames(Ymat) %||% sites)
      obs_col <- match(grid$response_id, colnames(Ymat) %||% species)
      obs <- Ymat[cbind(obs_row, obs_col)]
      prediction_rows[[length(prediction_rows) + 1]] <- data.frame(
        engine = "Hmsc", site_id = grid$site_id, response_id = grid$response_id,
        observed = as.numeric(obs), predicted_mean = as.numeric(as.vector(pred_mean)),
        predicted_lower = as.numeric(as.vector(pred_lower)), predicted_upper = as.numeric(as.vector(pred_upper)),
        prediction_context = "training",
        stringsAsFactors = FALSE)
    }
    if (length(new_prediction_rows) > 0) {
      write.csv(do.call(rbind, new_prediction_rows), file.path("standard", "newdata_predictions_long.csv"), row.names = FALSE)
      prediction_rows <- c(prediction_rows, new_prediction_rows)
    }
    predictions_long <- if (length(prediction_rows) > 0) do.call(rbind, prediction_rows) else data.frame(engine="Hmsc", site_id=NA_character_, response_id=NA_character_, observed=NA_real_, predicted_mean=NA_real_, predicted_lower=NA_real_, predicted_upper=NA_real_, prediction_context=NA_character_, stringsAsFactors=FALSE)
    write.csv(predictions_long, file.path("standard", "predictions_long.csv"), row.names = FALSE)

    association_rows <- list()
    if (isTRUE(cfg$outputs$omega %||% TRUE)) {
      for (j in seq_along(models)) {
        omega <- tryCatch(computeAssociations(models[[j]]), error = function(e) NULL)
        if (!is.null(omega)) {
          for (r in seq_along(omega)) {
            mat <- omega[[r]]$mean
            if (is.null(mat) || nrow(mat) < 2) next
            idx <- which(upper.tri(mat), arr.ind = TRUE)
            association_rows[[length(association_rows) + 1]] <- data.frame(
              engine = "Hmsc",
              response_1 = rownames(mat)[idx[, 1]] %||% paste0("sp_", idx[, 1]),
              response_2 = colnames(mat)[idx[, 2]] %||% paste0("sp_", idx[, 2]),
              association_type = "Omega",
              estimate = as.numeric(mat[idx]),
              comparable_level = paste0(names(models)[j], ":", names(omega)[r]),
              stringsAsFactors = FALSE)
          }
        }
      }
    }
    associations_long <- if (length(association_rows) > 0) do.call(rbind, association_rows) else data.frame(engine="Hmsc", response_1=NA_character_, response_2=NA_character_, association_type="Omega", estimate=NA_real_, comparable_level="no_random_level_or_association_disabled", stringsAsFactors=FALSE)
    write.csv(associations_long, file.path("standard", "associations_long.csv"), row.names = FALSE)

    plot_sources <- c(
      MCMC_convergence = file.path("results", "MCMC_convergence.pdf"),
      model_fit = file.path("results", paste0("model_fit_nfolds_", nfolds, ".pdf")),
      parameter_estimates = file.path("results", "parameter_estimates.pdf"),
      predictions = file.path("results", "predictions.pdf")
    )
    plot_manifest <- data.frame(plot = names(plot_sources), source = unname(plot_sources), copied_to = NA_character_, stringsAsFactors = FALSE)
    for (i in seq_along(plot_sources)) {
      if (file.exists(plot_sources[i])) {
        dest <- file.path("plots", paste0(names(plot_sources)[i], ".pdf"))
        file.copy(plot_sources[i], dest, overwrite = TRUE)
        plot_manifest$copied_to[i] <- dest
      }
    }
    write.csv(plot_manifest, file.path("plots", "plot_manifest.csv"), row.names = FALSE)

    write.csv(data.frame(run_id=basename(outdir), engine="Hmsc", status="fitted", n_sites=nrow(Ymat), n_responses=ncol(Ymat), n_predictors=ncol(XData), samples=samples, thin=thin, transient=transient, nChains=nChains, nParallel=nParallel, random_mode=random_mode, stringsAsFactors=FALSE),
              file.path("standard", "run_summary.csv"), row.names = FALSE)
    make_html_report(outdir, "Hmsc", cfg, "fitted")
    write.csv(data.frame(file=list.files(".", recursive=TRUE), stringsAsFactors=FALSE), file.path("standard", "output_manifest.csv"), row.names=FALSE)
    result$status <- "fitted"
  }, error = function(e) {
    result$status <<- "fit_failed"
    result$errors <<- c(result$errors, conditionMessage(e))
    safe_log("HMSC S1-S7 failed:", conditionMessage(e))
    try(writeLines(conditionMessage(e), file.path("diagnostics", "HMSC_S1S7_error.txt")), silent = TRUE)
  })

  result$runtime_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  diag_dir <- if (dir.exists("diagnostics")) "diagnostics" else file.path(outdir, "diagnostics")
  writeLines(jsonlite::toJSON(result, pretty = TRUE, auto_unbox = TRUE), file.path(diag_dir, "HMSC_S1S7_pipeline_status.json"))
  writeLines(jsonlite::toJSON(list(engine = "Hmsc", status = result$status, warnings = result$warnings, errors = result$errors, runtime_seconds = result$runtime_seconds),
                              pretty = TRUE, auto_unbox = TRUE),
             file.path(diag_dir, "engine_status.json"))
  writeLines(capture.output(sessionInfo()), file.path(diag_dir, "session_info.txt"))
  if (!file.exists(file.path(outdir, "report", "Hmsc_report.html"))) {
    try(make_html_report(outdir, "Hmsc", cfg, paste(result$status, paste(c(result$warnings, result$errors), collapse = "; "))), silent = TRUE)
  }
  result
}

validate_hmsc <- function(Y = NULL, X = NULL, Tr = NULL, study = NULL, coord = NULL,
                          newdata = NULL, newstudy = NULL, newcoord = NULL,
                          distr = "probit", XFormula = "~ .", TrFormula = "~ .",
                          use_traits = FALSE, use_phylogeny = FALSE,
                          random_mode = "none", random_effect_column = "sample",
                          spatial_method = "NNGP", nNeighbours = 10,
                          lon_col = "longitude", lat_col = "latitude",
                          samples = 50, transient = 25, thin = 1, nChains = 2, nParallel = 1,
                          nfMin = 1, nfMax = 10, nfolds = 2, partition_column = "sample",
                          additional_random_levels = character()) {
  msg <- c(studio_row_order_issue(Y, X, "XData.csv"),
           if (use_traits) studio_trait_order_issue(Y, Tr))
  ok <- !length(msg)
  mode_info <- hmsc_normalize_random_mode(random_mode, spatial_method)
  random_type <- mode_info$type
  spatial_method <- mode_info$spatial_method
  random_effect_column <- trimws(random_effect_column %||% "")
  if (!nzchar(random_effect_column)) random_effect_column <- "sample"
  grouping <- c(if (random_type == "sample") random_effect_column, as.character(additional_random_levels))
  grouping_check <- tryCatch({
    if (anyNA(grouping) || any(!nzchar(grouping)) || anyDuplicated(grouping) ||
        (random_type == "spatial" && "spatial_unit" %in% grouping)) {
      stop("Select one distinct studyDesign.csv column for each grouping level.", call. = FALSE)
    }
    if (!is.null(Y) && !is.null(rownames(Y))) {
      hmsc_align_rows(study, rownames(Y), "studyDesign.csv")
      hmsc_align_rows(coord, rownames(Y), "coordinates.csv")
    }
    for (col in grouping) {
      hmsc_group_factor(study, col)
      if (!is.null(newdata)) {
        nd_ids <- hmsc_prediction_site_ids(newdata)
        hmsc_group_factor(hmsc_align_rows(newstudy, nd_ids, "prediction_studyDesign.csv"),
                          col, "prediction_studyDesign.csv", training = FALSE)
      }
    }
    NULL
  }, error = conditionMessage)
  if (!is.null(grouping_check)) {
    ok <- FALSE
    msg <- c(msg, grouping_check)
  } else if (length(grouping)) {
    msg <- c(msg, paste("Grouping levels:", paste(grouping, collapse = ", ")))
  }
  if (is.null(Y)) {
    ok <- FALSE
    msg <- c(msg, "Y.csv is required. Rows must be sampling units/sites and columns must be species or response variables.")
  } else {
    ycheck <- numeric_matrix_check(Y, "HMSC Y")
    if (!isTRUE(ycheck$ok)) ok <- FALSE
    msg <- c(msg, ycheck$messages)
    fam_check <- response_family_messages(ycheck$matrix, distr, "HMSC")
    if (!isTRUE(fam_check$ok)) ok <- FALSE
    msg <- c(msg, fam_check$messages)
  }
  if (is.null(X)) {
    ok <- FALSE
    msg <- c(msg, "XData.csv is required. Rows must match Y rows and columns are environmental predictors.")
  } else if (!is.null(Y) && nrow(X) != nrow(Y)) {
    ok <- FALSE
    msg <- c(msg, paste0("Y and XData row counts must match. Y rows=", nrow(Y), ", XData rows=", nrow(X), "."))
  }
  # formula variable checks, best effort
  if (!is.null(X)) {
    X_clean <- clean_predictor_types(X, "XData")
    char_cols <- names(X)[vapply(X, is.character, logical(1))]
    factor_cols <- names(X_clean)[vapply(X_clean, is.factor, logical(1))]
    numeric_cols <- names(X_clean)[vapply(X_clean, is.numeric, logical(1))]
    if (length(char_cols) > 0) {
      msg <- c(msg, paste0("Character predictors detected in XData.csv. Numeric-looking columns will be converted to numeric; categorical columns will be converted to factors before HMSC fitting: ", paste(char_cols, collapse=", ")))
    }
    if (length(factor_cols) > 0) msg <- c(msg, paste0("Categorical XData predictors for HMSC: ", paste(factor_cols, collapse=", ")))
    if (length(numeric_cols) > 0) msg <- c(msg, paste0("Numeric XData predictors for HMSC: ", paste(numeric_cols, collapse=", ")))
  }
  if (!is.null(X) && nzchar(trimws(XFormula))) {
    f <- tryCatch(as.formula(XFormula), error = function(e) NULL)
    if (is.null(f)) {
      ok <- FALSE
      msg <- c(msg, paste0("XFormula is not a valid R formula: ", XFormula))
    } else {
      vars <- all.vars(f)
      missing <- setdiff(vars, colnames(X))
      if (length(missing) > 0 && XFormula != "~ .") {
        ok <- FALSE
        msg <- c(msg, paste0("XFormula variables missing from XData.csv: ", paste(missing, collapse=", ")))
      }
    }
  }
  if (!is.null(newdata)) {
    if (nrow(newdata) < 1) {
      ok <- FALSE
      msg <- c(msg, "newdata.csv has no rows. Remove it or provide at least one prediction site.")
    }
    if (is.null(X)) {
      ok <- FALSE
      msg <- c(msg, "newdata.csv was uploaded but XData.csv is missing; HMSC prediction needs the training XData schema.")
    } else {
      prep <- tryCatch(hmsc_prepare_prediction_xdata(newdata, X, XFormula), error = function(e) e)
      if (inherits(prep, "error")) {
        ok <- FALSE
        msg <- c(msg, conditionMessage(prep))
      } else {
        msg <- c(msg, prep$messages)
      }
    }
    if (random_type == "sample") {
      if (is.null(newstudy)) {
        ok <- FALSE
        msg <- c(msg, "HMSC newdata prediction with Sample/grouping random level requires prediction_studyDesign.csv.")
      } else if (nrow(newstudy) != nrow(newdata)) {
        ok <- FALSE
        msg <- c(msg, "prediction_studyDesign.csv rows must match newdata.csv rows.")
      } else if (!(random_effect_column %in% colnames(newstudy))) {
        ok <- FALSE
        msg <- c(msg, paste0("prediction_studyDesign.csv is missing grouping column '", random_effect_column, "'."))
      }
    }
    if (random_type == "spatial") {
      if (is.null(newcoord)) {
        ok <- FALSE
        msg <- c(msg, paste0("HMSC newdata prediction with spatial random level '", mode_info$mode, "' requires newcoords.csv."))
      } else if (nrow(newcoord) != nrow(newdata)) {
        ok <- FALSE
        msg <- c(msg, "newcoords.csv rows must match newdata.csv rows.")
      } else {
        coord_new_sel <- hmsc_select_coordinate_columns(newcoord, lon_col, lat_col)
        if (!isTRUE(coord_new_sel$ok)) {
          ok <- FALSE
          msg <- c(msg, paste0("newcoords.csv: ", coord_new_sel$message))
        } else if (!is.null(coord_new_sel$message)) {
          msg <- c(msg, paste0("newcoords.csv: ", coord_new_sel$message))
        }
      }
    }
  } else {
    if (!is.null(newstudy)) msg <- c(msg, "prediction_studyDesign.csv was uploaded but newdata.csv is missing; it will be ignored.")
    if (!is.null(newcoord)) msg <- c(msg, "newcoords.csv was uploaded but newdata.csv is missing; it will be ignored.")
  }
  if (isTRUE(use_traits)) {
    if (is.null(Tr)) {
      ok <- FALSE
      msg <- c(msg, "Use TrData / traits is enabled, but TrData.csv is missing.")
    } else if (!is.null(Y) && nrow(Tr) != ncol(Y)) {
      ok <- FALSE
      msg <- c(msg, paste0("TrData rows must match Y columns/species. TrData rows=", nrow(Tr), ", Y columns=", ncol(Y), "."))
    } else {
      tr_names <- rownames(Tr)
      y_names <- colnames(Y)
      has_named_traits <- !is.null(tr_names) && length(tr_names) == nrow(Tr) &&
        any(nzchar(tr_names)) && !identical(tr_names, as.character(seq_len(nrow(Tr))))
      if (!is.null(Y) && !is.null(y_names) && isTRUE(has_named_traits)) {
        missing_in_traits <- setdiff(y_names, tr_names)
        extra_in_traits <- setdiff(tr_names, y_names)
        if (length(missing_in_traits) > 0 || length(extra_in_traits) > 0) {
          ok <- FALSE
          msg <- c(msg, paste0("TrData row names must match Y species column names exactly. Missing in TrData: ",
                               paste(missing_in_traits %||% character(), collapse=", "),
                               "; extra in TrData: ", paste(extra_in_traits %||% character(), collapse=", "), "."))
        } else if (!identical(tr_names, y_names)) {
          ok <- FALSE
          msg <- c(msg, "TrData row order must match Y species column order. Reorder TrData rows before fitting so trait-mediated Gamma estimates map to the correct species.")
        } else {
          msg <- c(msg, "TrData row names exactly match Y species column names.")
        }
      }
      Tr_clean <- clean_predictor_types(Tr, "TrData")
      factor_cols <- names(Tr_clean)[vapply(Tr_clean, is.factor, logical(1))]
      if (length(factor_cols) > 0) msg <- c(msg, paste0("Categorical TrData variables for HMSC: ", paste(factor_cols, collapse=", ")))
      ftr <- tryCatch(as.formula(TrFormula), error = function(e) NULL)
      if (is.null(ftr)) {
        ok <- FALSE
        msg <- c(msg, paste0("TrFormula is not a valid R formula: ", TrFormula))
      } else {
        vars <- all.vars(ftr)
        missing <- setdiff(vars, colnames(Tr))
        if (length(missing) > 0 && TrFormula != "~ .") {
          ok <- FALSE
          msg <- c(msg, paste0("TrFormula variables missing from TrData.csv: ", paste(missing, collapse=", ")))
        }
      }
    }
  }
  if (isTRUE(use_phylogeny)) {
    msg <- c(msg, "Phylogeny is enabled. Make sure species names in the tree/matrix match Y column names exactly.")
  }
  if (!is.null(study) && !is.null(Y) && nrow(study) != nrow(Y)) {
    ok <- FALSE
    msg <- c(msg, "studyDesign.csv rows must match Y rows.")
  } else if (!is.null(study)) {
    msg <- c(msg, "All studyDesign.csv columns will be converted to factors before HMSC fitting, as required by Hmsc.")
  }
  if (random_type == "sample") {
    if (!nzchar(trimws(random_effect_column %||% ""))) random_effect_column <- "sample"
    if (is.null(study)) {
      ok <- FALSE
      msg <- c(msg, "Random effect mode 'sample' requires studyDesign.csv.")
    } else if (!(random_effect_column %in% colnames(study))) {
      ok <- FALSE
      msg <- c(msg, paste0("Grouping column '", random_effect_column, "' is not in studyDesign.csv. Available: ", paste(colnames(study), collapse=", ")))
    } else {
      msg <- c(msg, paste0("Sample random level will use grouping column '", random_effect_column, "'."))
    }
  }
  if (random_type == "spatial") {
    if (is.null(coord)) {
      ok <- FALSE
      msg <- c(msg, paste0("Spatial random level '", mode_info$mode, "' requires coordinates.csv."))
    } else {
      coord_sel <- hmsc_select_coordinate_columns(coord, lon_col, lat_col)
      if (!isTRUE(coord_sel$ok)) {
        ok <- FALSE
        msg <- c(msg, coord_sel$message)
      } else {
        if (!is.null(coord_sel$message)) msg <- c(msg, coord_sel$message)
        xy <- as.matrix(coord[, c(coord_sel$lon_col, coord_sel$lat_col), drop = FALSE])
        suppressWarnings(storage.mode(xy) <- "numeric")
        if (any(!is.finite(xy))) {
          ok <- FALSE
          msg <- c(msg, "Spatial coordinate columns must be numeric and finite.")
        } else if (any(duplicated(as.data.frame(xy, check.names = FALSE)))) {
          ok <- FALSE
          msg <- c(msg, "Spatial coordinate rows contain exact duplicates. Remove duplicate coordinates or aggregate duplicated sites before fitting a spatial HMSC random level.")
        }
      }
      if (!is.null(Y) && nrow(coord) != nrow(Y)) {
        ok <- FALSE
        msg <- c(msg, "coordinates.csv rows must match Y rows.")
      }
      if (spatial_method == "NNGP" && nNeighbours < 1) {
        ok <- FALSE
        msg <- c(msg, "NNGP requires nNeighbours >= 1.")
      }
      if (spatial_method == "GPP") {
        msg <- c(msg, "GPP spatial mode will use sKnot_file if present; otherwise it will generate knots with Hmsc::constructKnots().")
      }
    }
  }
  if (samples < 1 || transient < 0 || thin < 1 || nChains < 1 || nParallel < 1) {
    ok <- FALSE
    msg <- c(msg, "MCMC settings invalid: samples >= 1, transient >= 0, thin >= 1, nChains >= 1, nParallel >= 1 are required.")
  }
  if (nParallel > nChains) msg <- c(msg, "nParallel is larger than nChains. It will not speed up beyond the number of chains.")
  if (nfMin < 0 || nfMax < 1 || nfMin > nfMax) {
    ok <- FALSE
    msg <- c(msg, "Random-level latent factor prior must satisfy 0 <= nfMin <= nfMax.")
  }
  if (nfolds < 2) msg <- c(msg, "Cross-validation nfolds < 2. Cross-validation outputs should be skipped or set nfolds >= 2.")
  if (length(msg) == 0) msg <- "Hmsc data and settings check passed."
  list(ok = ok, messages = msg)
}

validate_boral <- function(Y = NULL, X = NULL, traits = NULL, rowids = NULL, ranefids = NULL,
                           distmat = NULL, offset = NULL, family = "negative.binomial",
                            family_text = "", num.lv = 2, lv.type = "independent",
                            row.eff = "none", n.burnin = 10000, n.iteration = 40000,
                            n.thin = 30, trial.size = 1, calc.ics = FALSE,
                            use_traits = FALSE, which.traits = "", use_ssvs = FALSE,
                            ssvs.index = "", save.model = FALSE, do.fit = TRUE,
                            formula.X = "~ .") {
  msg <- c(studio_row_order_issue(Y, X, "XData.csv"),
           if (use_traits) studio_trait_order_issue(Y, traits),
           studio_row_order_issue(Y, offset, "offset.csv"))
  ok <- !length(msg)
  allowed <- c("binomial", "poisson", "negative.binomial", "normal", "tweedie", "exponential",
               "gamma", "lnormal", "beta", "ordinal", "ztpoisson", "ztnegative.binomial")
  if (is.null(Y)) {
    ok <- FALSE
    msg <- c(msg, "Y.csv is required. Rows should be sites/samples and columns should be species/responses.")
  }
  fam_vec <- NULL
  if (nzchar(trimws(family_text))) {
    fam_vec <- tolower(trimws(unlist(strsplit(family_text, ","))))
  } else {
    fam_vec <- tolower(trimws(family))
  }
  fam_vec[fam_vec == "lognormal"] <- "lnormal"
  if (identical(lv.type, "power.exponential")) lv.type <- "powered.exponential"
  if (!is.null(Y) && length(fam_vec) == 1) fam_vec <- rep(fam_vec, ncol(Y))
  if (!is.null(Y) && length(fam_vec) != ncol(Y)) {
    ok <- FALSE
    msg <- c(msg, paste0("Family vector length must be 1 or equal to Y columns. Current length=", length(fam_vec), ", Y columns=", ncol(Y), "."))
  }
  bad_fam <- setdiff(fam_vec, allowed)
  if (length(bad_fam) > 0) {
    ok <- FALSE
    msg <- c(msg, paste0("Invalid boral family value(s): ", paste(unique(bad_fam), collapse=", "), ". Allowed: ", paste(allowed, collapse=", ")))
  }
  if (!is.null(Y)) {
    ycheck <- numeric_matrix_check(Y, "boral Y")
    if (!isTRUE(ycheck$ok)) ok <- FALSE
    msg <- c(msg, ycheck$messages)
    Ym <- ycheck$matrix
    for (fam in unique(fam_vec)) {
      cols <- fam_vec == fam
      if (identical(fam, "binomial")) {
        trial_num <- suppressWarnings(as.numeric(trial.size))
        if (!is.finite(trial_num) || trial_num < 1) {
          ok <- FALSE
          msg <- c(msg, "boral binomial trial.size must be a positive number.")
        } else {
          vals <- Ym[, cols, drop = FALSE]
          if (any(vals < 0 | abs(vals - round(vals)) > 1e-8 | vals > trial_num, na.rm = TRUE)) {
            ok <- FALSE
            msg <- c(msg, paste0("boral binomial responses must be integer successes between 0 and trial.size=", trial_num, "."))
          }
        }
      } else {
        fam_check <- response_family_messages(Ym[, cols, drop = FALSE], fam, paste("boral", fam))
        if (!isTRUE(fam_check$ok)) ok <- FALSE
        msg <- c(msg, fam_check$messages)
      }
    }
  }
  if (!is.null(X)) {
    X_clean <- clean_predictor_types(X, "boral XData")
    factor_cols <- names(X_clean)[vapply(X_clean, is.factor, logical(1))]
    if (length(factor_cols) > 0) msg <- c(msg, paste0("boral XData categorical predictors will be factors: ", paste(factor_cols, collapse=", ")))
    fcheck <- validate_one_sided_formula(formula.X, X_clean, "boral formula.X")
    if (!isTRUE(fcheck$ok)) ok <- FALSE
    msg <- c(msg, fcheck$messages)
  }
  if (!is.null(traits)) {
    traits_clean <- clean_predictor_types(traits, "boral traits")
    factor_cols <- names(traits_clean)[vapply(traits_clean, is.factor, logical(1))]
    if (length(factor_cols) > 0) msg <- c(msg, paste0("boral traits categorical variables will be factors: ", paste(factor_cols, collapse=", ")))
  }
  if (!is.null(rowids)) msg <- c(msg, "boral row.ids will be treated as grouping factors.")
  if (!is.null(ranefids)) msg <- c(msg, "boral ranef.ids will be treated as grouping factors.")
  if (!is.null(X) && !is.null(Y) && nrow(X) != nrow(Y)) {
    ok <- FALSE
    msg <- c(msg, paste0("XData rows must match Y rows. X rows=", nrow(X), ", Y rows=", nrow(Y), "."))
  }
  if (isTRUE(use_traits) && is.null(traits)) {
    ok <- FALSE
    msg <- c(msg, "Use traits is enabled, but traits.csv is not uploaded.")
  }
  if (!is.null(traits) && !is.null(Y) && nrow(traits) != ncol(Y)) {
    ok <- FALSE
    msg <- c(msg, paste0("traits rows should match Y columns/species. traits rows=", nrow(traits), ", Y columns=", ncol(Y), "."))
  }
  if (!is.null(rowids) && !is.null(Y) && nrow(rowids) != nrow(Y)) {
    ok <- FALSE
    msg <- c(msg, paste0("row.ids rows should match Y rows. row.ids rows=", nrow(rowids), ", Y rows=", nrow(Y), "."))
  }
  if (!is.null(ranefids) && !is.null(Y) && nrow(ranefids) != nrow(Y)) {
    ok <- FALSE
    msg <- c(msg, paste0("ranef.ids rows should match Y rows. ranef.ids rows=", nrow(ranefids), ", Y rows=", nrow(Y), "."))
  }
  if (!is.null(offset) && !is.null(Y)) {
    if (!(nrow(offset) == nrow(Y) && ncol(offset) == ncol(Y))) {
      ok <- FALSE
      msg <- c(msg, paste0("offset matrix should have same dimensions as Y. offset=", nrow(offset), "x", ncol(offset), ", Y=", nrow(Y), "x", ncol(Y), "."))
    }
  }
  spatial_lvs <- c("exponential","squared.exponential","powered.exponential","spherical")
  if (lv.type %in% spatial_lvs && is.null(distmat)) {
    ok <- FALSE
    msg <- c(msg, paste0("lv.type='", lv.type, "' requires a distance matrix distmat.csv."))
  }
  if (!is.null(distmat) && !is.null(Y) && (nrow(distmat) != nrow(Y) || ncol(distmat) != nrow(Y))) {
    ok <- FALSE
    msg <- c(msg, paste0("distmat.csv should be a square n_sites x n_sites matrix. distmat=", nrow(distmat), "x", ncol(distmat), ", Y rows=", nrow(Y), "."))
  }
  if (num.lv < 0) {
    ok <- FALSE
    msg <- c(msg, "num.lv must be >= 0.")
  }
  if (num.lv == 0 && lv.type != "independent") {
    msg <- c(msg, "lv.type is ignored when num.lv = 0.")
  }
  if (row.eff != "none" && is.null(rowids)) {
    msg <- c(msg, "row.eff is not 'none' and row.ids.csv is missing. boral can auto-create row-specific IDs, but uploading row.ids.csv is recommended for reproducibility.")
  }
  if (row.eff == "random" && !is.null(rowids)) {
    msg <- c(msg, "row.eff='random' and row.ids are supplied. They will be converted to grouping-factor integer IDs before fitting.")
  }
  if (n.burnin < 0 || n.iteration <= 0 || n.thin <= 0) {
    ok <- FALSE
    msg <- c(msg, "n.burnin must be >=0; n.iteration and n.thin must be >0.")
  }
  if (n.burnin >= n.iteration) {
    msg <- c(msg, "n.burnin is >= n.iteration. In boral this may leave little/no posterior sampling depending on control interpretation. Check settings.")
  }
  if (n.iteration %% n.thin != 0) {
    msg <- c(msg, "n.iteration should usually be divisible by n.thin for clean posterior storage.")
  }
  if (trial.size < 1) {
    ok <- FALSE
    msg <- c(msg, "trial.size must be >= 1.")
  }
  if (isTRUE(use_traits) && !nzchar(trimws(which.traits))) {
    msg <- c(msg, "Traits are enabled but which.traits is empty. Adapter will use a default all-traits-for-all-coefficients rule unless specified.")
  }
  if (isTRUE(use_ssvs) && !nzchar(trimws(ssvs.index))) {
    msg <- c(msg, "SSVS is enabled but ssvs.index is empty. Adapter will need a default or explicit selection scheme.")
  }
  if (isTRUE(save.model)) {
    msg <- c(msg, "save.model=TRUE can be memory-consuming because raw JAGS model/samples can be returned.")
  }
  if (!isTRUE(do.fit)) {
    msg <- c(msg, "do.fit=FALSE selected: boral will build the model object/script but not run MCMC.")
  }
  if (length(msg) == 0) msg <- "boral basic data and setting checks passed."
  list(ok = ok, messages = msg, family_vec = fam_vec)
}

validate_sjsdm <- function(Y=NULL, env=NULL, spatial=NULL, traits=NULL, newdata=NULL,
                           family="binomial_probit", env_model="linear", spatial_model="none",
                           env_formula="~ .", spatial_formula="~ 0 + .",
                           biotic_lambda=0, biotic_alpha=0.5, iter=100, sampling=5000,
                           learning_rate=0.003, step_size=50, parallel=0, cv_k=0,
                           dnn_hidden="10,10,10", dropout=0, device="cpu", dtype="float32",
                           verbose=TRUE, optimizer="Adamax", weight_decay=0.002,
                           scheduler=0, lr_reduce_factor=0.99, early_stopping_training=0,
                           mixed=FALSE, generate_spatial_ev=FALSE, spatial_ev_threshold=0,
                           anova_samples=5000, tune_steps=0) {
  msg <- c(studio_row_order_issue(Y, env, "env.csv"),
           studio_row_order_issue(Y, spatial, "spatial.csv"))
  ok <- !length(msg)
  family_key <- tolower(trimws(as.character(family %||% "binomial_probit")))
  if (family_key %in% c("negative_binomial_log", "negative.binomial", "negative_binomial")) family_key <- "nbinom"
  env_model_key <- tolower(trimws(as.character(env_model %||% "linear")))
  if (env_model_key %in% c("none", "intercept", "intercept only", "intercept_only")) env_model_key <- "intercept_only"
  spatial_model_key <- tolower(trimws(as.character(spatial_model %||% "none")))
  if (spatial_model_key %in% c("spatial eigenvectors", "spatial_eigenvectors", "eigenvectors")) spatial_model_key <- "eigenvectors"
  if (is.null(Y)) { ok <- FALSE; msg <- c(msg, "Y.csv is required. Rows should be sites/samples and columns should be species/OTUs/responses.") }
  if (is.null(env) && env_model_key != "intercept_only") {
    ok <- FALSE
    msg <- c(msg, "env.csv / XData.csv is required unless env module = intercept-only. sjSDM always needs an environmental design matrix, even if it contains only an intercept.")
  }
  if (is.null(env) && env_model_key == "intercept_only") msg <- c(msg, "No env.csv uploaded and env module is intercept-only; the reproducible script will build a one-column intercept design.")
  if (!is.null(Y) && !is.null(env) && nrow(Y) != nrow(env)) { ok <- FALSE; msg <- c(msg, paste0("Y and env rows must match. Y rows=", nrow(Y), ", env rows=", nrow(env), ".")) }
  if (!is.null(Y)) {
    ycheck <- numeric_matrix_check(Y, "sjSDM Y")
    if (!isTRUE(ycheck$ok)) ok <- FALSE
    msg <- c(msg, ycheck$messages)
    fam_check <- response_family_messages(ycheck$matrix, family_key, "sjSDM")
    if (!isTRUE(fam_check$ok)) ok <- FALSE
    msg <- c(msg, fam_check$messages)
    if (anyNA(ycheck$matrix)) {
      ok <- FALSE
      msg <- c(msg, "sjSDM needs a complete training response matrix. Missing responses are not zeros; prepare the analysis data before fitting.")
    }
  }
  if (!is.null(env)) {
    env_clean <- clean_predictor_types(env, "sjSDM env")
    factor_cols <- names(env_clean)[vapply(env_clean, is.factor, logical(1))]
    if (length(factor_cols) > 0) msg <- c(msg, paste0("sjSDM categorical env predictors will be factors/model-matrix encoded by the production adapter: ", paste(factor_cols, collapse=", ")))
    if (env_model_key != "intercept_only") {
      fcheck <- validate_one_sided_formula(env_formula, env_clean, "sjSDM env formula")
      if (!isTRUE(fcheck$ok)) ok <- FALSE
      msg <- c(msg, fcheck$messages)
    }
  }
  if (isTRUE(generate_spatial_ev) && spatial_model_key == "none") { ok <- FALSE; msg <- c(msg, "generateSpatialEV is enabled, but spatial module is none. Choose spatial eigenvectors or linear spatial module so the generated eigenvectors are actually used.") }
  if ((spatial_model_key != "none" || isTRUE(generate_spatial_ev)) && is.null(spatial)) { ok <- FALSE; msg <- c(msg, "Spatial model/eigenvector generation selected but spatial.csv / coordinates / spatial predictors were not uploaded.") }
  if (!is.null(spatial) && !is.null(Y) && nrow(spatial) != nrow(Y)) { ok <- FALSE; msg <- c(msg, paste0("Spatial rows should match Y rows. spatial rows=", nrow(spatial), ", Y rows=", nrow(Y), ".")) }
  if (!is.null(spatial) && (spatial_model_key != "none" || isTRUE(generate_spatial_ev))) {
    spatial_clean <- clean_predictor_types(spatial, "sjSDM spatial")
    if (isTRUE(generate_spatial_ev)) {
      smat <- suppressWarnings(as.matrix(data.frame(lapply(spatial_clean, as.numeric))))
      if (ncol(smat) < 2 || any(!is.finite(smat[, 1:2, drop=FALSE]))) {
        ok <- FALSE
        msg <- c(msg, "generateSpatialEV requires at least two numeric coordinate columns in spatial.csv.")
      }
    } else {
      sfcheck <- validate_one_sided_formula(spatial_formula, spatial_clean, "sjSDM spatial formula")
      if (!isTRUE(sfcheck$ok)) ok <- FALSE
      msg <- c(msg, sfcheck$messages)
    }
  }
  if (!is.null(newdata) && !is.null(env) && env_model_key != "intercept_only") {
    miss_new <- setdiff(all.vars(tryCatch(as.formula(env_formula), error=function(e) ~ .)), colnames(newdata))
    if (length(miss_new) > 0 && trimws(env_formula) != "~ .") {
      ok <- FALSE
      msg <- c(msg, paste0("newdata.csv is missing env formula variable(s): ", paste(miss_new, collapse=", ")))
    }
  }
  if (!is.null(traits) && !is.null(Y) && nrow(traits) != ncol(Y)) msg <- c(msg, paste0("Trait rows should usually match Y columns/species. trait rows=", nrow(traits), ", Y columns=", ncol(Y), "."))
  if (biotic_lambda < 0 || biotic_alpha < 0 || biotic_alpha > 1) { ok <- FALSE; msg <- c(msg, "bioticStruct requires lambda >= 0 and alpha between 0 and 1.") }
  if (iter <= 0 || sampling <= 0 || learning_rate <= 0 || step_size <= 0 || parallel < 0) { ok <- FALSE; msg <- c(msg, "iter, sampling, learning_rate and step_size must be positive; parallel must be >= 0.") }
  if (dropout < 0 || dropout >= 1) { ok <- FALSE; msg <- c(msg, "DNN dropout must be in [0, 1).") }
  if (env_model_key == "dnn" || spatial_model_key == "dnn") {
    h <- suppressWarnings(as.integer(trimws(unlist(strsplit(dnn_hidden, ",")))))
    if (length(h)==0 || any(is.na(h)) || any(h <= 0)) { ok <- FALSE; msg <- c(msg, "DNN hidden layers must be comma-separated positive integers, e.g. 10,10,10.") }
  }
  if (!(family_key %in% c("binomial_probit", "binomial_logit", "poisson_log", "nbinom", "gaussian_identity"))) { ok <- FALSE; msg <- c(msg, "Unsupported sjSDM family. Use binomial_probit, binomial_logit, poisson_log, nbinom, or gaussian_identity.") }
  if (!(device %in% c("cpu", "gpu"))) { ok <- FALSE; msg <- c(msg, "sjSDM device must be 'cpu' or 'gpu' for version 1.0.7. Use gpu only when torch CUDA is available.") }
  if (device == "gpu") msg <- c(msg, "GPU selected. The reproducible script will use device='gpu'; sjSDM/PyTorch CUDA must be installed and checked before fitting.")
  if (!(dtype %in% c("float32", "float64"))) { ok <- FALSE; msg <- c(msg, "dtype must be float32 or float64.") }
  if (!(optimizer %in% c("Adamax", "RMSprop", "SGD", "AccSGD", "AdaBound", "madgrad"))) { ok <- FALSE; msg <- c(msg, "Unsupported sjSDM optimizer. DiffGrad is present in source but is not exported by sjSDM 1.0.7, so it is not exposed here.") }
  if (weight_decay < 0) { ok <- FALSE; msg <- c(msg, "optimizer weight_decay must be >= 0.") }
  if (scheduler < 0 || early_stopping_training < 0) { ok <- FALSE; msg <- c(msg, "sjSDMControl scheduler and early_stopping_training must be >= 0. Use 0 to disable.") }
  if (lr_reduce_factor <= 0 || lr_reduce_factor >= 1) { ok <- FALSE; msg <- c(msg, "sjSDMControl lr_reduce_factor must be between 0 and 1.") }
  if (isTRUE(mixed) && device != "gpu") msg <- c(msg, "mixed half-precision training is intended for newer GPUs; keep it off on CPU.")
  if (anova_samples < 1000) msg <- c(msg, "ANOVA / variation-partitioning samples are low; larger samples are recommended for stable results.")
  if (tune_steps > 0 && cv_k < 2) msg <- c(msg, "Regularization tuning requested but CV folds < 2. Cross-validation tuning needs at least 2 folds.")
  if (length(msg)==0) msg <- "sjSDM basic data and setting checks passed."
  list(ok=ok, messages=msg)
}

validate_spoccupancy <- function(y = NULL, occ = NULL, det = NULL, coords = NULL, species = NULL,
                                 integrated = NULL, model_type = "PGOcc",
                                 n.batch = 1000, batch.length = 25, n.burn = 0, n.thin = 1,
                                 n.chains = 1, n.factors = 3, NNGP = FALSE, n.neighbors = 15,
                                 k.fold = 0, svc.cols = "", occ.formula = "~ 1", det.formula = "~ 1") {
  ok <- TRUE
  msg <- c()
  spatial_models <- c("spPGOcc", "spMsPGOcc", "sfJSDM", "sfMsPGOcc", "spIntPGOcc", "stPGOcc", "stMsPGOcc",
                      "svcPGBinom", "svcPGOcc", "svcTPGBinom", "svcTPGOcc", "svcMsPGOcc", "svcTMsPGOcc")
  ms_models <- c("msPGOcc", "spMsPGOcc", "lfJSDM", "sfJSDM", "lfMsPGOcc", "sfMsPGOcc",
                 "tMsPGOcc", "stMsPGOcc", "svcMsPGOcc", "svcTMsPGOcc")
  int_models <- c("intPGOcc", "spIntPGOcc", "stIntPGOcc", "intMsPGOcc")
  lf_models <- c("lfJSDM", "sfJSDM", "lfMsPGOcc", "sfMsPGOcc")
  svc_models <- c("svcPGBinom", "svcPGOcc", "svcTPGBinom", "svcTPGOcc", "svcMsPGOcc", "svcTMsPGOcc")
  if (is.null(y)) {
    ok <- FALSE
    msg <- c(msg, "Detection/nondetection y data are required. For occupancy models this is usually sites x replicates, or species x sites x replicates for multi-species models.")
  } else {
    ycheck <- numeric_matrix_check(y, "spOccupancy y")
    if (!isTRUE(ycheck$ok)) ok <- FALSE
    msg <- c(msg, ycheck$messages)
    fam_check <- response_family_messages(ycheck$matrix, "binomial", "spOccupancy")
    if (!isTRUE(fam_check$ok)) ok <- FALSE
    msg <- c(msg, fam_check$messages)
  }
  if (is.null(occ)) {
    msg <- c(msg, "occ.covs / occurrence covariates were not uploaded. The model can only use intercept-only occurrence if the adapter allows it.")
  } else if (!is.null(y)) {
    if (nrow(occ) != nrow(y) && !(model_type %in% ms_models)) {
      msg <- c(msg, paste0("For single-species models, occ.covs rows should match y rows. occ rows=", nrow(occ), ", y rows=", nrow(y), "."))
    }
    occ_clean <- clean_predictor_types(occ, "spOccupancy occ.covs")
    factor_cols <- names(occ_clean)[vapply(occ_clean, is.factor, logical(1))]
    if (length(factor_cols) > 0) msg <- c(msg, paste0("spOccupancy occurrence categorical covariates will be factors/model-matrix encoded by the production adapter: ", paste(factor_cols, collapse=", ")))
    fcheck <- validate_one_sided_formula(occ.formula, occ_clean, "spOccupancy occ.formula")
    if (!isTRUE(fcheck$ok)) ok <- FALSE
    msg <- c(msg, fcheck$messages)
  }
  if (is.null(det)) {
    msg <- c(msg, "det.covs / detection covariates were not uploaded. The model can only use intercept-only detection if the adapter allows it.")
  } else {
    det_clean <- clean_predictor_types(det, "spOccupancy det.covs")
    factor_cols <- names(det_clean)[vapply(det_clean, is.factor, logical(1))]
    if (length(factor_cols) > 0) msg <- c(msg, paste0("spOccupancy detection categorical covariates will be factors/model-matrix encoded by the production adapter: ", paste(factor_cols, collapse=", ")))
    fcheck <- validate_one_sided_formula(det.formula, det_clean, "spOccupancy det.formula")
    if (!isTRUE(fcheck$ok)) ok <- FALSE
    msg <- c(msg, fcheck$messages)
  }
  if (model_type %in% spatial_models && is.null(coords)) {
    ok <- FALSE
    msg <- c(msg, paste0(model_type, " is spatial and requires coords.csv."))
  }
  if (model_type %in% ms_models && is.null(species)) {
    msg <- c(msg, "Multi-species models benefit from species names/traits metadata. If omitted, species names will be inferred from y dimnames or column names when possible.")
  }
  if (model_type %in% int_models && is.null(integrated)) {
    ok <- FALSE
    msg <- c(msg, paste0(model_type, " is an integrated model and requires integrated data-source information."))
  }
  if (n.batch <= 0 || batch.length <= 0) {
    ok <- FALSE
    msg <- c(msg, "n.batch and batch.length must both be positive.")
  }
  total <- n.batch * batch.length
  if (n.burn >= total) {
    ok <- FALSE
    msg <- c(msg, paste0("n.burn must be smaller than total iterations n.batch * batch.length = ", total, "."))
  }
  if (n.thin <= 0 || ((total - n.burn) %% n.thin) != 0) {
    msg <- c(msg, "For clean posterior storage, n.thin should divide n.batch * batch.length - n.burn.")
  }
  if (n.chains < 1) {
    ok <- FALSE
    msg <- c(msg, "n.chains must be at least 1.")
  }
  if (model_type %in% lf_models && n.factors < 1) {
    ok <- FALSE
    msg <- c(msg, "Latent-factor models require n.factors >= 1.")
  }
  if (isTRUE(NNGP) && n.neighbors < 1) {
    ok <- FALSE
    msg <- c(msg, "NNGP requires n.neighbors >= 1.")
  }
  if (k.fold < 0) {
    ok <- FALSE
    msg <- c(msg, "k.fold must be zero or positive.")
  }
  if (model_type %in% svc_models && !nzchar(svc.cols)) {
    msg <- c(msg, "SVC model selected but svc.cols is empty. Add covariate indices/names for spatially varying coefficients.")
  }
  if (!grepl("^~", trimws(occ.formula))) msg <- c(msg, "occ.formula should be one-sided and start with ~, e.g. ~ elev + forest.")
  if (!grepl("^~", trimws(det.formula))) msg <- c(msg, "det.formula should be one-sided and start with ~, e.g. ~ day + observer.")
  if (length(msg) == 0) msg <- "spOccupancy basic data and setting checks passed."
  list(ok = ok, messages = msg)
}

validate_gjam <- function(Y, X, type_table = NULL, type_text = "", single_type = "DA",
                          fcgroups = "", ccgroups = "", ng = 2000, burnin = 500,
                          holdoutN = 0, random = "", notStandard = "", formula_text = "~ .") {
  ok <- TRUE
  msg <- c()
  allowed <- c("PA","CON","CA","DA","FC","CC","OC","CAT")
  if (is.null(Y)) { ok <- FALSE; msg <- c(msg, "Y.csv is required and could not be read.") }
  if (is.null(X)) { ok <- FALSE; msg <- c(msg, "XData.csv is required and could not be read.") }
  if (!is.null(Y) && !is.null(X) && nrow(Y) != nrow(X)) {
    ok <- FALSE
    msg <- c(msg, paste0("Y and XData row numbers differ: Y=", nrow(Y), ", XData=", nrow(X), "."))
  }
  if (!is.null(X) && nzchar(random) && !(random %in% colnames(X))) {
    ok <- FALSE
    msg <- c(msg, paste0("random column '", random, "' is not present in XData."))
  }
  if (!is.null(X)) {
    X_clean <- clean_predictor_types(X, "GJAM XData")
    factor_cols <- names(X_clean)[vapply(X_clean, is.factor, logical(1))]
    if (length(factor_cols) > 0) msg <- c(msg, paste0("GJAM categorical predictors detected and will be converted to factors/model matrix columns by the production adapter: ", paste(factor_cols, collapse=", ")))
    fcheck <- validate_one_sided_formula(formula_text, X_clean, "GJAM formula")
    if (!isTRUE(fcheck$ok)) ok <- FALSE
    msg <- c(msg, fcheck$messages)
  }
  if (!is.null(X) && nzchar(notStandard)) {
    ns <- trimws(unlist(strsplit(notStandard, ",")))
    ns <- ns[nzchar(ns)]
    bad <- setdiff(ns, colnames(X))
    if (length(bad) > 0) msg <- c(msg, paste0("notStandard contains columns not found in XData: ", paste(bad, collapse=", ")))
  }
  type_vec <- NULL
  if (nzchar(type_text)) {
    type_vec <- toupper(trimws(unlist(strsplit(type_text, ","))))
  } else if (!is.null(type_table) && ncol(type_table) >= 1) {
    # Look for a column named typeName/type/typeNames; otherwise use last column.
    cn <- tolower(colnames(type_table))
    idx <- which(cn %in% c("typename","typenames","type","types"))[1]
    if (is.na(idx)) idx <- ncol(type_table)
    type_vec <- toupper(trimws(as.character(type_table[[idx]])))
  } else {
    type_vec <- toupper(single_type)
  }
  if (length(type_vec) == 1 && !is.null(Y)) type_vec <- rep(type_vec, ncol(Y))
  if (!is.null(Y) && length(type_vec) != ncol(Y)) {
    ok <- FALSE
    msg <- c(msg, paste0("typeNames length must be 1 or equal to number of Y columns. Current length=", length(type_vec), ", Y columns=", ncol(Y), "."))
  }
  bad_types <- setdiff(type_vec, allowed)
  if (length(bad_types) > 0) {
    ok <- FALSE
    msg <- c(msg, paste0("Invalid GJAM typeNames: ", paste(unique(bad_types), collapse=", "), ". Allowed: ", paste(allowed, collapse=", ")))
  }
  if (!is.null(Y)) {
    ycheck <- numeric_matrix_check(Y, "GJAM Y")
    if (!isTRUE(ycheck$ok)) ok <- FALSE
    msg <- c(msg, ycheck$messages)
    Ym <- ycheck$matrix
    for (typ in unique(type_vec)) {
      cols <- type_vec == typ
      fam_check <- response_family_messages(Ym[, cols, drop = FALSE], typ, paste("GJAM", typ))
      if (!isTRUE(fam_check$ok)) ok <- FALSE
      msg <- c(msg, fam_check$messages)
    }
    if ("FC" %in% type_vec && !nzchar(fcgroups)) {
      msg <- c(msg, "FC fractional composition columns detected; FCgroups should be supplied for composition groups.")
    }
    if ("CC" %in% type_vec && !nzchar(ccgroups)) {
      msg <- c(msg, "CC count composition columns detected; CCgroups should be supplied for composition groups.")
    }
  }
  if (burnin >= ng) {
    ok <- FALSE
    msg <- c(msg, "burnin must be less than ng.")
  }
  if (!is.null(Y) && holdoutN >= nrow(Y)) {
    ok <- FALSE
    msg <- c(msg, "holdoutN must be smaller than the number of observations.")
  }
  if (length(msg) == 0) msg <- "GJAM basic data and setting checks passed."
  list(ok = ok, messages = msg, type_vec = type_vec)
}

validate_basic <- function(Y, X, Tr = NULL, study = NULL, coord = NULL, engine = "Hmsc") {
  ok <- TRUE
  msg <- c()
  if (is.null(Y)) { ok <- FALSE; msg <- c(msg, "Y.csv is required and could not be read.") }
  if (is.null(X)) { ok <- FALSE; msg <- c(msg, "XData.csv is required and could not be read.") }
  if (!is.null(Y) && !is.null(X) && nrow(Y) != nrow(X)) {
    ok <- FALSE
    msg <- c(msg, paste0("Y and XData row numbers differ: Y=", nrow(Y), ", XData=", nrow(X), "."))
  }
  if (!is.null(Y) && anyNA(Y)) msg <- c(msg, "Y contains missing values. Check distribution and engine support.")
  if (!is.null(X) && anyNA(X)) msg <- c(msg, "XData contains missing values. Most model engines fail when predictors contain NA.")
  if (!is.null(Tr) && !is.null(Y) && nrow(Tr) != ncol(Y)) {
    msg <- c(msg, paste0("TrData rows should usually match Y columns. TrData rows=", nrow(Tr), ", Y columns=", ncol(Y), "."))
  }
  if (!is.null(study) && !is.null(Y) && nrow(study) != nrow(Y)) {
    msg <- c(msg, paste0("studyDesign rows should match Y rows. studyDesign rows=", nrow(study), ", Y rows=", nrow(Y), "."))
  }
  if (!is.null(coord) && !is.null(Y) && nrow(coord) != nrow(Y)) {
    msg <- c(msg, paste0("coordinates rows should match Y rows. coordinates rows=", nrow(coord), ", Y rows=", nrow(Y), "."))
  }
  if (engine == "jSDM" && !is.null(Y)) {
    # jSDM has discrete choices, warn about response values.
    if (all(Y %in% c(0, 1), na.rm = TRUE)) {
      msg <- c(msg, "Y appears binary; jSDM binomial probit/logit is suitable.")
    } else if (all(Y >= 0 & abs(Y - round(Y)) < 1e-8, na.rm = TRUE)) {
      msg <- c(msg, "Y appears non-negative integer count data; jSDM poisson_log may be suitable.")
    } else {
      msg <- c(msg, "Y appears continuous or mixed; jSDM_gaussian may be suitable for continuous responses, but mixed response types should not be forced into jSDM.")
    }
  }
  if (length(msg) == 0) msg <- paste0(engine, " basic data checks passed.")
  list(ok = ok, messages = msg)
}

augment_jsdm_check <- function(check, Y, X, model_type = "binomial_probit", site_formula = "~ .") {
  if (is.null(check)) check <- list(ok = TRUE, messages = character())
  if (!is.null(Y)) {
    ycheck <- numeric_matrix_check(Y, "jSDM Y")
    fam_check <- response_family_messages(ycheck$matrix, model_type, "jSDM")
    check$ok <- isTRUE(check$ok) && isTRUE(ycheck$ok) && isTRUE(fam_check$ok)
    check$messages <- c(check$messages, ycheck$messages, fam_check$messages)
  }
  if (!is.null(X)) {
    x_clean <- clean_predictor_types(X, "jSDM XData")
    factor_cols <- names(x_clean)[vapply(x_clean, is.factor, logical(1))]
    if (length(factor_cols) > 0) {
      check$messages <- c(check$messages, paste0("jSDM categorical predictors will be factors/model-matrix encoded by the production adapter: ", paste(factor_cols, collapse=", ")))
    }
    fcheck <- validate_one_sided_formula(site_formula, x_clean, "jSDM site_formula")
    check$ok <- isTRUE(check$ok) && isTRUE(fcheck$ok)
    check$messages <- c(check$messages, fcheck$messages)
  }
  if (length(check$messages) == 0) check$messages <- "jSDM data and settings check passed."
  check
}

augment_formula_check <- function(check, data, formula_text, label) {
  if (is.null(check)) check <- list(ok = TRUE, messages = character())
  if (!is.null(data)) {
    cleaned <- clean_predictor_types(data, label)
    fcheck <- validate_one_sided_formula(formula_text, cleaned, label)
    check$ok <- isTRUE(check$ok) && isTRUE(fcheck$ok)
    check$messages <- c(check$messages, fcheck$messages)
  }
  check
}

make_html_report <- function(outdir, engine, cfg, status_text) {
  report <- file.path(outdir, "report", paste0(engine, "_report.html"))
  html <- paste0(
    "<!DOCTYPE html><html><head><meta charset='UTF-8'><title>", engine, " Report</title>",
    "<style>body{font-family:Segoe UI,Arial,sans-serif;margin:40px;line-height:1.6;color:#14213d}",
    "h1{font-size:34px}h2{color:#1d4ed8}code{background:#eef2ff;padding:2px 6px;border-radius:4px}",
    ".box{border:1px solid #dbe5f0;border-radius:18px;padding:18px;margin:16px 0;background:#fff}",
    "li{margin:6px 0}</style></head><body>",
    "<h1>", engine, " workflow report</h1>",
    "<div class='box'><b>Status:</b> ", htmlEscape(status_text), "</div>",
    "<h2>Project</h2><p>", htmlEscape(cfg$project_name %||% "JSDMWorkbench"), "</p>",
    "<h2>Engine-specific workflow</h2>",
    "<p>This report comes from the ", htmlEscape(engine), " workflow. Different engines have different data requirements, assumptions, parameters and outputs.</p>",
    "<h2>Output structure</h2><ul>",
     "<li><code>used_config.yml</code>: saved configuration</li>",
     "<li><code>inputs/</code>: copied input files</li>",
     "<li><code>data/</code>: reproducible copies used by exported scripts</li>",
     "<li><code>models/</code>: model objects when fitting is enabled and package is available</li>",
    "<li><code>tables/</code>: tabular summaries</li>",
    "<li><code>plots/</code>: figures</li>",
    "<li><code>diagnostics/</code>: checks and logs</li>",
    "</ul>",
    "<h2>Interpretation reminder</h2>",
    "<p>Compare ecological conclusions and prediction behaviour across engines only after checking whether the models were run on compatible data and settings.</p>",
    "</body></html>"
  )
  write_lines(html, report)
  report
}

jsdm_adapter_file <- file.path(app_dir, "R", "jsdm_adapter.R")
if (file.exists(jsdm_adapter_file)) {
  source(jsdm_adapter_file, local = TRUE)
}
gjam_adapter_file <- file.path(app_dir, "R", "gjam_adapter.R")
if (file.exists(gjam_adapter_file)) {
  source(gjam_adapter_file, local = TRUE)
}
spoccupancy_adapter_file <- file.path(app_dir, "R", "spoccupancy_adapter.R")
if (file.exists(spoccupancy_adapter_file)) {
  source(spoccupancy_adapter_file, local = TRUE)
}
boral_adapter_file <- file.path(app_dir, "R", "boral_adapter.R")
if (file.exists(boral_adapter_file)) {
  source(boral_adapter_file, local = TRUE)
}
hmschpc_adapter_file <- file.path(app_dir, "R", "hmschpc_adapter.R")
if (file.exists(hmschpc_adapter_file)) {
  source(hmschpc_adapter_file, local = TRUE)
}
source(file.path(app_dir, "R", "empirical_benchmark.R"), local = TRUE)

if (!exists("hmschpc_blank", mode = "function")) {
  hmschpc_blank <- function(x) {
    is.null(x) || length(x) == 0 || all(is.na(x)) || !nzchar(trimws(as.character(x)[1]))
  }
}

if (!exists("hmschpc_python_source_dir", mode = "function")) {
  hmschpc_python_source_dir <- function() {
    normalizePath(file.path(app_dir, "external_packages", "hmsc-hpc-main"),
                  winslash = "/", mustWork = FALSE)
  }
}

if (!exists("hmschpc_python_bin", mode = "function")) {
  hmschpc_python_bin <- function(cfg = list()) {
    runtime <- cfg$runtime %||% list()
    py <- runtime$python %||% ""
    if (!hmschpc_blank(py)) return(unname(as.character(py)[1]))

    local_appdata <- Sys.getenv("LOCALAPPDATA", unset = "")
    candidates <- c(
      Sys.getenv("RETICULATE_PYTHON", unset = ""),
      file.path(local_appdata, "JSDMStudio", "hmschpc-python310", "Scripts", "python.exe"),
      file.path(local_appdata, "R-MINI~1", "envs", "r-sjsdm", "python.exe"),
      file.path(local_appdata, "r-miniconda", "envs", "r-sjsdm", "python.exe"),
      file.path(local_appdata, "r-miniconda", "python.exe"),
      unname(Sys.which("python")),
      unname(Sys.which("py"))
    )
    candidates <- candidates[nzchar(candidates)]
    candidates <- candidates[file.exists(candidates) | basename(candidates) %in% c("python", "py")]
    if (length(candidates) > 0) return(normalizePath(candidates[1], winslash = "/", mustWork = FALSE))
    "python"
  }
}

# ------------------------------------------------------------
# Theme and UI helpers
# ------------------------------------------------------------

theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#1d4ed8",
  secondary = "#475569",
  success = "#15803d",
  danger = "#b91c1c",
  warning = "#b45309",
  info = "#0e7490",
  base_font = font_collection("system-ui", "Segoe UI", "sans-serif"),
  code_font = font_collection("ui-monospace", "Consolas", "monospace")
)

css <- HTML("
:root{--ink:#0f172a;--muted:#64748b;--line:#dbe5f0;--blue:#1d4ed8;--cyan:#0e7490;--green:#15803d;--amber:#b45309;--red:#b91c1c;}
body{background:linear-gradient(180deg,#f8fbff 0%,#f3f7fb 45%,#eef3f9 100%);color:var(--ink);}
.section-title{font-size:25px;font-weight:900;letter-spacing:-.025em;margin:22px 0 5px;}
.section-subtitle{color:var(--muted);margin-bottom:16px;}
.cardx{
  background:rgba(255,255,255,.94);border:1px solid var(--line);border-radius:24px;padding:20px;margin-bottom:17px;
  box-shadow:0 14px 34px rgba(15,23,42,.075);
}
.cardx h3{font-size:19px;font-weight:900;letter-spacing:-.015em;margin-top:0;}
.cardx h4{font-size:16px;font-weight:850;margin-top:8px;}
.metric-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:14px;margin:18px 0;}
.metric{background:#fff;border:1px solid var(--line);border-radius:22px;padding:18px;box-shadow:0 10px 28px rgba(15,23,42,.06);}
.metric .icon{font-size:26px;margin-bottom:8px}.metric .label{text-transform:uppercase;font-size:11px;letter-spacing:.1em;color:var(--muted);font-weight:850}.metric .value{font-size:22px;font-weight:900;margin-top:3px}
.note{border-left:5px solid var(--blue);background:#eff6ff;padding:14px 16px;border-radius:15px;margin:12px 0;}
.warn{border-left:5px solid var(--amber);background:#fff7ed;padding:14px 16px;border-radius:15px;margin:12px 0;}
.safe{border-left:5px solid var(--green);background:#f0fdf4;padding:14px 16px;border-radius:15px;margin:12px 0;}
.danger{border-left:5px solid var(--red);background:#fef2f2;padding:14px 16px;border-radius:15px;margin:12px 0;}
.engine-card{border:1px solid var(--line);border-radius:24px;padding:20px;background:#fff;min-height:245px;box-shadow:0 14px 32px rgba(15,23,42,.07);}
.engine-card h3{display:flex;align-items:center;gap:8px;}
.engine-tag{font-size:11px;padding:5px 9px;border-radius:999px;background:#eef2ff;color:#3730a3;font-weight:850;}
.param-help{font-size:13px;color:#475569;background:#f8fafc;border:1px solid #e2e8f0;border-radius:13px;padding:10px;margin:7px 0 12px 0;}
.accordion-button{font-weight:850;}
.form-label{font-weight:750;color:#334155;}
.btn{border-radius:14px;font-weight:800;}
pre.logbox{background:#0f172a;color:#d1fae5;padding:16px;border-radius:18px;min-height:240px;max-height:420px;overflow:auto;font-size:12px;border:1px solid #1e293b;}
.small-muted{font-size:13px;color:var(--muted);}
.pill{display:inline-block;padding:5px 9px;border-radius:999px;background:#eef2ff;color:#3730a3;font-size:12px;font-weight:850;margin:2px;}
.nature-ribbon{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px;margin:16px 0 20px 0;}
.visual-card{border:1px solid #cbd5e1;background:linear-gradient(180deg,#fff,#f8fafc);border-radius:18px;padding:14px;display:flex;gap:12px;align-items:flex-start;min-height:92px;}
.visual-card .mark{width:38px;height:38px;border-radius:12px;background:#0f766e;color:white;display:inline-flex;align-items:center;justify-content:center;font-weight:950;box-shadow:0 8px 20px rgba(15,118,110,.22);}
.visual-card b{display:block;margin-bottom:3px;}
.sync-pair{background:#ffffff;border:1px solid #e2e8f0;border-radius:16px;padding:11px 12px;margin-bottom:12px;box-shadow:0 8px 20px rgba(15,23,42,.045);}
.sync-pair .form-group,.sync-pair .mb-3{margin-bottom:6px!important;}
.sync-control-row{display:grid;grid-template-columns:minmax(0,1fr) 118px;gap:12px;align-items:end;}
.sync-slider-stack .form-label{display:block;margin-bottom:8px;}
.sync-slider-row{display:flex;justify-content:space-between;gap:10px;margin-top:4px;}
.sync-range{width:100%;accent-color:#1d4ed8;}
.sync-scale{font-size:11px;color:#64748b;font-weight:750;white-space:nowrap;}
@media (max-width: 700px){.sync-control-row{grid-template-columns:1fr;}.sync-scale{white-space:normal;}}
.file-preview-shell{background:#fff;border:1px solid var(--line);border-radius:20px;padding:16px;box-shadow:0 10px 28px rgba(15,23,42,.06);}
.file-preview-meta{background:#f8fafc;border:1px solid #e2e8f0;border-radius:14px;padding:12px;margin:8px 0 14px 0;}
.mode-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin:10px 0 16px 0;}
.mode-card{border:1px solid #dbe5f0;border-radius:16px;background:#fff;padding:12px;min-height:102px;}
.mode-card b{display:block;color:#0f172a;margin-bottom:5px;}
.mode-card .step-dot{display:inline-flex;width:24px;height:24px;border-radius:50%;align-items:center;justify-content:center;background:#1d4ed8;color:#fff;font-weight:900;margin-right:6px;}
.conditional-box{border:1px solid #bfdbfe;background:#f8fbff;border-radius:18px;padding:14px;margin:12px 0;}
")

sync_script <- HTML("
(function() {
  function asNumber(value) {
    var n = parseFloat(value);
    return isNaN(n) ? null : n;
  }
  function syncRangeToNumber(range) {
    var pair = $(range).closest('.sync-pair');
    var id = pair.data('syncId');
    var target = $('#' + id);
    if (!target.length) return;
    target.val($(range).val());
    target.trigger('change');
  }
  function syncNumberToRange(numberInput) {
    var pair = $(numberInput).closest('.sync-pair');
    var range = pair.find('.sync-range');
    if (!range.length) return;
    var value = asNumber($(numberInput).val());
    if (value === null) return;
    var minValue = asNumber(range.attr('min'));
    var maxValue = asNumber(range.attr('max'));
    if (minValue !== null) value = Math.max(minValue, value);
    if (maxValue !== null) value = Math.min(maxValue, value);
    range.val(value);
  }
  $(document).on('input change', '.sync-range', function() {
    syncRangeToNumber(this);
  });
  $(document).on('input change', '.sync-number input', function() {
    syncNumberToRange(this);
  });
  $(document).on('shiny:bound', function() {
    $('.sync-pair').each(function() {
      var input = $(this).find('.sync-number input');
      if (input.length) syncNumberToRange(input[0]);
    });
  });
})();
")

studio_symbol <- function(kind, class = NULL, mode = "art") {
  span(class = paste(c("studio-icon", class), collapse = " "),
       `data-studio-icon` = kind, `data-icon-mode` = mode, `aria-hidden` = "true")
}

studio_asset <- function(path) paste0(path, "?v=", studio_version)

studio_app_icon <- paste0("studio-app.svg?v=", unname(tools::md5sum(file.path("www", "studio-app.svg"))))

home_platform <- function(name, symbol, mode, style) {
  tags$li(
    class = paste("studio-platform", style),
    div(class = "studio-platform-device", `aria-hidden` = "true",
        div(class = "studio-platform-screen", icon(symbol))),
    tags$strong(name),
    span(class = "studio-platform-mode", mode)
  )
}
home_workflow_step <- function(name, symbol, detail) {
  tags$li(
    class = "studio-flow-step",
    studio_symbol(symbol, "studio-flow-icon"),
    div(tags$h3(name), p(detail))
  )
}
section_header <- function(title, subtitle = NULL, symbol = NULL) {
  tagList(div(class = paste0("section-title", if (!is.null(symbol)) " studio-title"),
              if (!is.null(symbol)) studio_symbol(symbol, "studio-title-icon"), title),
          if (!is.null(subtitle)) div(class = "section-subtitle", subtitle))
}
help_text <- function(x) div(class="param-help", x)
workflow_overview <- function(symbol, title, description, requirement, notes) {
  stages <- list(
    inputs = c("Inputs", "upload"), validation = c("Check data", "check"),
    settings = c("Model", "tuning"), results = c("Run & results", "report")
  )
  tags$header(class = "studio-workflow-overview",
    div(class = "studio-workflow-heading",
      studio_symbol(symbol, "studio-engine-icon"),
      div(h2(title), p(description))),
    div(class = "studio-workflow-guidance",
      p(class = "studio-workflow-requirement", tags$strong("Before fitting"), requirement),
      tags$details(class = "studio-workflow-notes",
        tags$summary(studio_symbol("guide", mode = "glyph"), "Interpretation notes"),
        tags$ul(lapply(notes, tags$li)))),
    tags$nav(class = "studio-workflow-stages", `aria-label` = "Workflow sections",
      lapply(names(stages), function(stage) {
        tags$button(class = "studio-workflow-stage", type = "button", `data-workflow-stage` = stage,
          studio_symbol(stages[[stage]][2L], mode = "glyph"), tags$span(stages[[stage]][1L]))
      }))
  )
}
file_preview_block <- function(prefix, title = "Uploaded data preview") {
  div(class="file-preview-shell",
      h3(title),
      fluidRow(
        column(4,
          uiOutput(paste0(prefix, "_file_preview_selector")),
          div(class="file-preview-meta", verbatimTextOutput(paste0(prefix, "_file_preview_meta")))
        ),
        column(8, DTOutput(paste0(prefix, "_file_preview_table")))
      ),
      div(class="small-muted", "CSV/TSV previews show up to 100 rows; tree/text files show the first lines; binary files are listed but not parsed.")
  )
}
reviewer_lens_card <- function(mark, title, body, checks = NULL, caution = NULL) {
  div(class="cardx",
      h3(title),
      p(class="small-muted", body),
      if (!is.null(checks)) tagList(tags$b("Reviewer checks"), tags$ul(lapply(checks, tags$li))),
      if (!is.null(caution)) div(class="warn", caution)
  )
}
parameter_story_card <- function(title, why, check, pitfall = NULL) {
  div(class="cardx",
      h3(title),
      p(class="small-muted", why),
      tags$b("Before fitting"),
      p(check),
      if (!is.null(pitfall)) div(class="warn", pitfall)
  )
}
software_review_notice <- function() {
  div(class="note",
      tags$b("Scientific review mode: "),
      "The interface now explains each major modelling choice as a data contract, an ecological claim and a reproducibility requirement. Use the text as a checklist before treating any fitted object as an ecological result."
  )
}
parameter_dictionary_data <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- studio_parameter_reference(ui)
    cache
  }
})
plain_num <- function(x) {
  format(as.numeric(x), scientific = FALSE, trim = TRUE, big.mark = ",")
}
synced_numeric_slider <- function(id, label, value, min = 0, max = NULL, step = 1, width = "100%", help = NULL) {
  min_num <- suppressWarnings(as.numeric(min))
  step_num <- suppressWarnings(as.numeric(step))
  value_num <- suppressWarnings(as.numeric(value))
  if (!is.finite(min_num)) min_num <- 0
  if (!is.finite(step_num) || step_num <= 0) step_num <- 1
  if (!is.finite(value_num)) value_num <- min_num
  range_max <- suppressWarnings(as.numeric(max))
  if (is.null(max) || length(max) == 0 || !is.finite(range_max)) {
    range_max <- base::max(value_num * 5, min_num + step_num * 100, min_num + 1, na.rm = TRUE)
  }
  if (!is.finite(range_max) || range_max <= min_num) range_max <- min_num + step_num * 100
  div(
    class = "sync-pair",
    `data-sync-id` = id,
    div(
      class = "sync-control-row",
      div(
        class = "sync-slider-stack",
        tags$label(class = "form-label", `for` = paste0(id, "_range"), label),
        tags$input(id = paste0(id, "_range"), type = "range", class = "sync-range", min = min_num, max = range_max, step = step_num, value = value_num),
        div(class = "sync-slider-row", span(class = "sync-scale", paste0("Min ", plain_num(min_num))), span(class = "sync-scale", paste0("Max ", plain_num(range_max))))
      ),
      div(class = "sync-number", numericInput(id, "Value", value = value, min = min, max = max, step = step, width = width))
    ),
    if (!is.null(help)) help_text(help)
  )
}

# ------------------------------------------------------------
# UI
# ------------------------------------------------------------


# ---- Compatibility helpers for imported FULL_ENGLISH_FIXED Hmsc panel ----
# These helpers are required by the restored Hmsc UI.
if (!exists("guide_block")) {
  guide_block <- function(title, body) {
    tags$details(
      open = NA,
      tags$summary(title),
      div(class = "param-help", HTML(body))
    )
  }
}
if (!exists("eta_card")) {
  eta_card <- function(title, value, detail = "") {
    div(class = "eta-card",
        div(class = "eta-title", title),
        div(class = "eta-value", value),
        if (!is.null(detail) && nzchar(as.character(detail))) div(class = "eta-detail", detail) else NULL)
  }
}
if (!exists("param_box")) {
  param_box <- function(title, ui, help = NULL) {
    div(class = "param-box",
        div(class = "param-title", title),
        ui,
        if (!is.null(help)) div(class = "param-help", help) else NULL)
  }
}
if (!exists("num_value")) {
  num_value <- function(input, id, default = NULL) {
    val <- input[[id]]
    if (is.null(val) || length(val) == 0 || is.na(val)) return(default)
    as.numeric(val)
  }
}
if (!exists("%||%")) {
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
}
if (!exists("safe_text")) {
  safe_text <- function(x, default = "") {
    if (is.null(x) || length(x) == 0 || is.na(x)) default else as.character(x)
  }
}

studio_engine_catalog <- list(
  Hmsc = list("engine-hmsc", "Traits, phylogeny and multilevel community models.",
              "Community responses and site predictors; optional traits, tree/taxonomy and study design.",
              "Species effects, trait effects, phylogenetic signal, random-level associations and predictions."),
  "Hmsc-HPC" = list("engine-hmschpc", "CPU pyhmsc workflow for the supported HMSC subset.",
                    "Community responses and predictors; optional traits, phylogeny and supported iid or Full spatial levels. TensorFlow required.",
                    "Native posterior draws, parameter summaries and supported predictions. NNGP and GPP are not available in this sampler."),
  jSDM = list("engine-jsdm", "Bayesian regression with latent factors.",
              "Site-by-species data, matching site predictors and optional traits for supported functions.",
              "Species coefficients and latent-factor residual correlations; prediction depends on the chosen function."),
  GJAM = list("engine-gjam", "Joint regression for mixed observation scales.",
              "Responses, predictors and a typeName for each response column. Supply effort or censoring only where applicable.",
              "Response summaries, corMu, predictions and supported inverse predictions."),
  spOccupancy = list("engine-spoccupancy", "Occupancy with an explicit detection process.",
                     "Repeated detection/non-detection surveys and occurrence/detection covariates; coordinates for spatial variants.",
                     "Occurrence and detection coefficients, occupancy probabilities and model-specific diagnostic summaries."),
  sjSDM = list("engine-sjsdm", "Regularised community regression and species covariance.",
               "Responses, environmental predictors and optional spatial predictors. Python/PyTorch required.",
               "Environmental effects, fitted covariance and predictions on the selected response scale."),
  boral = list("engine-boral", "Bayesian ordination and correlated-response regression.",
               "Community responses, predictors and optional latent, grouping or trait structure. JAGS required.",
               "Latent scores, species effects, residual correlations and selected fit summaries.")
)

studio_open_page <- function(page, label = "Open workflow") {
  tags$a(href = "#main_nav", class = "btn btn-outline-secondary studio-link", `data-studio-nav` = page,
         label, studio_symbol("arrow", mode = "glyph"))
}

studio_external_link <- function(path, label) {
  tags$a(href = paste0("https://github.com/gujiqi/JSDM_Studio/", path), target = "_blank",
         rel = "noopener noreferrer", class = "btn btn-outline-secondary studio-link",
         label, studio_symbol("external", mode = "glyph"))
}

studio_guide_row <- function(title, symbol, ..., action = NULL) {
  tags$article(class = "studio-guide-row",
    studio_symbol(symbol, "studio-guide-icon"),
    div(h3(title), lapply(list(...), p)),
    action
  )
}

studio_engine_row <- function(engine) {
  entry <- studio_engine_catalog[[engine]]
  page <- if (engine == "sjSDM") "s-jSDM Workflow" else paste(engine, "Workflow")
  studio_guide_row(engine, entry[[1]], entry[[2]], paste("Inputs:", entry[[3]]),
                   paste("Outputs:", entry[[4]]), action = studio_open_page(page))
}

studio_home_engine <- function(engine) {
  page <- if (engine == "sjSDM") "s-jSDM Workflow" else paste(engine, "Workflow")
  tags$li(tags$a(class = "studio-home-engine", href = "#main_nav", `data-studio-nav` = page,
    studio_symbol(studio_engine_catalog[[engine]][[1]], "studio-engine-chip"), span(engine)))
}

studio_project_match_row <- function(row, leading = FALSE) {
  engine <- row$Engine
  page <- if (engine == "sjSDM") "s-jSDM Workflow" else paste(engine, "Workflow")
  label <- if (leading) "Leading match" else switch(row$Status,
    direct = "Match", alternative = "Alternative", provisional = "Provisional",
    needs_input = "Needs input", not_applicable = "Not applicable", unavailable = "Not connected")
  points <- if (is.na(row$Score)) "N/A" else paste(row$Score, "pts")
  tags$article(class = "studio-guide-row studio-match-row", `data-engine` = engine,
    `data-match-state` = row$Status, `data-match-score` = if (!is.na(row$Score)) row$Score,
    studio_symbol(studio_engine_catalog[[engine]][[1]], "studio-guide-icon"),
    div(class = "studio-match-copy",
      h4(engine), span(class = "studio-match-status", label),
      p(row$Why),
      if (nzchar(row$Route)) p(class = "studio-match-route", row$Route),
      if (nzchar(row$Points)) tags$details(class = "studio-match-breakdown",
        tags$summary("Score details"), p(row$Points), p(row$Caution))),
    div(class = "studio-match-action",
      span(class = "studio-match-score", points),
      htmltools::tagAppendAttributes(studio_open_page(page, "Open"),
        `aria-label` = paste("Open", engine, "workflow"), title = paste("Open", engine, "workflow")))
  )
}

ui <- navbarPage(
  title = span(class = "studio-brand",
    tags$img(src = studio_app_icon, width = 31, height = 31, alt = ""), "JSDM Studio"),
  theme = theme,
  header = tagList(tags$head(
    tags$style(css), tags$link(rel = "stylesheet", href = studio_asset("studio-home.css")),
    tags$link(rel = "stylesheet", href = studio_asset("studio-workbench.css")),
    tags$link(rel = "stylesheet", href = studio_asset("studio-icons.css")),
    tags$link(rel = "stylesheet", href = studio_asset("studio-localization.css")),
    tags$link(rel = "icon", href = studio_app_icon, type = "image/svg+xml"),
    tags$script(sync_script), tags$script(src = studio_asset("studio-icons.js"), defer = NA),
    tags$script(src = studio_asset("studio-workbench.js"), defer = NA),
    tags$script(src = studio_asset("studio-localization.js"), defer = NA)
  ), studio_preferences_ui()),
  id = "main_nav",
  selected = "Home",
  collapsible = TRUE,

  tabPanel("Home",
    tags$main(class = "studio-home",
      tags$header(class = "studio-intro",
        div(class = "studio-home-inner studio-intro-grid",
          div(
            tags$h1("JSDM Studio"),
            p(class = "studio-home-lead", "A workbench for joint species distribution modelling."),
            div(class = "studio-home-actions",
              actionButton("home_open_project", "Open Project", class = "studio-home-action studio-home-primary"),
              actionButton("home_open_engine_guide", "Engine Guide", class = "studio-home-action studio-home-secondary")
            )
          ),
          tags$section(class = "studio-platforms", `aria-labelledby` = "studio-platform-title",
            tags$h2(id = "studio-platform-title", "Run locally"),
            tags$ul(class = "studio-platform-list",
              home_platform("Windows", "windows", "Desktop or browser", "studio-platform-win"),
              home_platform("macOS", "apple", "Local browser", "studio-platform-mac"),
              home_platform("Linux", "linux", "Local browser", "studio-platform-linux")
            )
          )
        )
      ),
      tags$section(class = "studio-flow", `aria-labelledby` = "studio-flow-title",
        div(class = "studio-home-inner",
          tags$h2(id = "studio-flow-title", "Analysis workflow"),
          tags$ol(class = "studio-flow-list", role = "list",
            home_workflow_step("Upload", "upload", "Supply the data your engine needs."),
            home_workflow_step("Check", "check", "Review inputs and model settings."),
            home_workflow_step("Fit", "run", "Run the selected model."),
            home_workflow_step("Diagnose", "diagnostic", "Inspect model checks and warnings."),
            home_workflow_step("Export", "export", "Save results, settings and scripts.")
          )
        )
      ),
      tags$section(class = "studio-engine-strip studio-home-inner", `aria-labelledby` = "studio-engine-title",
        tags$h2(id = "studio-engine-title", "Model workflows"),
        tags$ul(class = "studio-engine-list", role = "list",
          lapply(names(studio_engine_catalog), studio_home_engine)
        ),
        p("Inputs, diagnostics and exports depend on the selected engine and settings.")
      )
    )
  ),

  tabPanel("Project",
    fluidPage(
      tags$main(class = "studio-project-layout",
        section_header("Project", "Study design and workflow selection"),
        div(class = "studio-project-grid",
          tags$section(class = "studio-project-section studio-project-summary",
            h3("Workflow matches"),
            p(class = "studio-match-note", "Rule-based match points, not model accuracy. Based on the declared profile, not checked data."),
            uiOutput("project_recommendation_summary")
          ),
          tags$section(class = "studio-project-section studio-project-profile",
            h3("Research question and data"),
            selectInput("question_template", "Question template", choices = names(studio_project_questions)),
            help_text("The quantity or process to estimate. This does not set a model formula."),
            selectInput("project_response_structure", "Response structure",
              choices = c("One response type" = "single response type",
                          "Repeated detection/non-detection surveys" = "detection-nondetection occupancy surveys",
                          "Mixed response types" = "mixed response types", "Composition" = "composition data",
                          "Not specified" = "unknown")),
            help_text("Repeated detection histories contain separate visits or replicates, not one 0/1 value per site."),
            conditionalPanel("input.project_response_structure === 'single response type'",
              selectInput("project_response_family", "Response family", choices = studio_project_families),
              help_text("Choose the form of the response as analysed. Integer counts and genuinely bounded proportions need different models.")),
            numericInput("project_n_responses", "Species / response columns", 20, min = 1, step = 1),
            help_text("Species or response columns, excluding IDs; count a repeatedly surveyed species only once."),
            conditionalPanel("input.question_template === 'Integrated occupancy with multiple data sources'",
              checkboxInput("project_has_integrated_sources", "Multiple aligned survey sources available", FALSE),
              help_text("The connected single-species integration workflow also needs source/prefix metadata.")),
            checkboxInput("project_has_traits", "Species traits available", FALSE),
            help_text("One row per species, matched to the response columns."),
            checkboxInput("project_has_phylogeny", "Tree or taxonomy available", FALSE),
            help_text("A species-aligned Newick tree or taxonomy table for a clearly labelled surrogate."),
            checkboxInput("project_has_spatial", "Spatial coordinates available", FALSE),
            help_text("Training site coordinates with a documented coordinate reference system."),
            checkboxInput("project_need_prediction", "New-site prediction needed", FALSE),
            help_text("Predictions for new sites with the required matching predictors and spatial inputs."),
            tags$details(class = "studio-disclosure", tags$summary("Project metadata"),
              textInput("project_name", "Project name", value = "JSDM_Studio_project"),
              help_text("A name saved with the analysis. Example: Finnish_birds_2014."),
              textAreaInput("project_question", "Scientific question", rows = 2),
              help_text("Records the hypothesis without changing the fitted formula."),
              numericInput("project_n_sites", "Observations / sites", 100, min = 1, step = 1),
              help_text("Declared number of sampling units; not a check of adequate sample size."),
              numericInput("project_n_predictors", "Candidate predictors", 5, min = 0, step = 1),
              help_text("Number of measured explanatory variables before factor and interaction expansion."),
              checkboxInput("project_many_zeros", "Many zero observations", FALSE),
              help_text("Many zeros are not, by themselves, evidence of zero inflation."),
              DTOutput("project_metadata_table"))
          ),
          tags$section(id = "project_matches", class = "studio-project-section studio-project-results",
            uiOutput("project_recommendation_cards"),
            tags$details(class = "studio-disclosure", tags$summary("Scoring rules"),
              p(class = "studio-match-note",
                "+2 compatible response; +6 direct objective or +2 alternative objective; +1 per relevant additional supported feature. Incompatible workflows and missing required inputs are not scored."))
          )
        )
      )
    )
  ),

  tabPanel("Engine Guide",
    tags$main(class = "studio-guide",
      section_header("Engine Guide", "Choose a model based on the observation process and the quantity you need to estimate."),
      div(class = "studio-guide-grid", lapply(names(studio_engine_catalog), studio_engine_row)),
      tags$details(class = "studio-disclosure", tags$summary("What can be compared across engines?"),
        DTOutput("comparison_meaning_table"))
    )
  ),

  tabPanel("Support",
    tags$main(class = "studio-guide",
      section_header("Support", paste("JSDM Studio", studio_version)),
      studio_guide_row("Guides and parameter reference", "guide",
        "Input formats, model interpretation and replay.",
        action = studio_open_page("Guides", "Open guides")),
      studio_guide_row("Updates and release notes", "update",
        "Check published versions and choose an installer for your operating system.",
        action = studio_external_link("releases", "View releases")),
      studio_guide_row("Feedback and bug reports", "feedback",
        "Include the application version, engine, error message and steps to reproduce. Do not post private data or logs.",
        action = studio_external_link("issues/new", "Open feedback form"))
    )
  ),

  # ----------------------------------------------------------
  # Hmsc workflow tab
  # ----------------------------------------------------------
  


  # ----------------------------------------------------------
  # jSDM workflow tab
  # ----------------------------------------------------------

  navbarMenu("Workflows",
  tabPanel("Hmsc", value = "Hmsc Workflow",
    fluidPage(
      workflow_overview("engine-hmsc", "Hmsc workflow",
        "Hierarchical community models with species traits, phylogeny and random levels.",
        "Match the response family to Y. Spatial projections need predictor rows and coordinates for the new sites.",
        c("For a minimal binary test, use Y.csv with 0/1 values and XData.csv, select probit, and leave random effects, traits and phylogeny off.",
          "Beta describes environmental effects on the selected model scale. Gamma describes trait-related differences in those effects when a matched trait table is included.",
          "Omega summarizes residual association at a specified random level. It does not identify direct interactions.",
          "Full, NNGP and GPP use different spatial inputs. Supply training coordinates separately from newcoords.csv; newdata.csv holds the matching predictors for prediction.",
          "A quick run checks the workflow. For inference, choose an adequate chain length and check convergence.")),
      accordion(
        class = "studio-workflow-form", open = "inputs",
        accordion_panel("1. HMSC data upload", value = "inputs",
          h4("Training inputs"),
          fluidRow(
            column(6,
              fileInput("hmsc_Y_file", "Y.csv  - response matrix", accept = ".csv"),
              help_text("Required. Rows = sites/samples, columns = species/responses. Example for default probit: 0/1 presence-absence matrix. Example columns: sp_1, sp_2, sp_3."),
              fileInput("hmsc_X_file", "XData.csv  - environmental predictors", accept = ".csv"),
              help_text("Required. Rows must match Y. Columns are environmental variables. Numeric and categorical columns are both supported. Character categorical columns such as substrate will be converted to factor before HMSC fitting. Example: pH, moisture, canopy, substrate, elevation. Default XFormula = ~ . uses all columns."),
              fileInput("hmsc_Tr_file", "TrData.csv  - traits / response attributes", accept = ".csv"),
              help_text("Optional. Rows must match Y columns/species. Example: life_form, height_mm, reproduction. Use only when 'Use TrData / traits' is checked.")
            ),
            column(6,
              fileInput("hmsc_study_file", "studyDesign.csv  - grouping / random-effect design", accept = ".csv"),
              help_text("Required for grouping levels. One row per Y row; add site_id to align by ID, or keep rows in Y order. Example columns: plot, site, year."),
              fileInput("hmsc_coord_file", "coordinates.csv  - spatial coordinates", accept = ".csv"),
              help_text("Optional. Rows must match Y. Example columns: longitude, latitude or x, y. Required only when Random effect mode = spatial."),
              fileInput("hmsc_phylo_file", "phylogeny.nwk or taxonomy.csv", accept = c(".tre", ".tree", ".nwk", ".txt", ".csv", ".tsv")),
              help_text("Optional. Upload a Newick tree with species tips, or a taxonomy table with species plus kingdom/phylum/class/order/family/genus/species columns. Do not upload a species-by-species covariance/correlation matrix; JSDM Studio derives that matrix internally when needed.")
            )
          ),
          hr(),
          h4("Prediction inputs for new sites"),
          div(class="safe",
              tags$b("Spatial prediction rule: "),
              "Use newdata.csv for the environmental predictors at the new points and newcoords.csv for their longitude/latitude or x/y coordinates. ",
              "Coordinates alone are not enough unless the model formula does not use environmental predictors."),
          fluidRow(
            column(4,
              fileInput("hmsc_newdata_file", "newdata.csv  - prediction environmental predictors", accept = ".csv"),
              help_text("Optional. Rows = new sites or grid cells to predict. Columns must contain the XFormula variables from XData.csv. A site_id/sample_id/id column is allowed and will become prediction site IDs.")
            ),
            column(4,
              fileInput("hmsc_newcoords_file", "newcoords.csv  - prediction longitude/latitude or x/y", accept = ".csv"),
              help_text("Conditional. Required for Chapter-11-style spatial prediction from Spatial Full, Spatial NNGP or Spatial GPP HMSC models. Rows must match newdata.csv; use longitude/latitude or x/y columns matching the coordinate controls.")
            ),
            column(4,
              fileInput("hmsc_newstudy_file", "prediction_studyDesign.csv  - prediction random-effect design", accept = ".csv"),
              help_text("Required for prediction with grouping levels. Supply one row per newdata row and the same grouping columns, for example plot and year. Use site_id to align rows. Match known group labels exactly.")
            )
          )
        ),
        accordion_panel("2. HMSC data check", value = "validation",
          fluidRow(
            column(4,
              actionButton("hmsc_check", "Check HMSC data", class = "btn-primary"),
              div(class="note", tags$b("Check before running."), " This verifies Y/X dimensions, response distribution, formula variables, trait dimensions and random-effect requirements.")
            ),
            column(8, div(class="cardx", h3("HMSC check messages"), verbatimTextOutput("hmsc_check_messages")))
          ),
          div(class="cardx", h3("HMSC data dimensions"), DTOutput("hmsc_data_table")),
          file_preview_block("hmsc", "HMSC uploaded data preview")
        ),
        accordion_panel("3. HMSC core model settings", value = "settings",
          fluidRow(
            column(4, selectInput("hmsc_distr", "distr  - response distribution", choices = c("probit", "poisson", "normal"), selected = "probit"),
                   help_text("Default/example: probit. Use probit for 0/1 presence-absence, poisson for non-negative integer counts, normal for continuous responses.")),
            column(4, textInput("hmsc_XFormula", "XFormula  - environmental formula", value = "~ ."),
                   help_text("Default/example: ~ . uses all XData columns. Included example data support: ~ pH + moisture + canopy. Categorical predictors are allowed and will be treated as factors.")),
            column(4, textInput("hmsc_TrFormula", "TrFormula  - trait formula", value = "~ ."),
                   help_text("Default/example: ~ . uses all trait columns. Ignored unless Use TrData / traits is checked. Example: ~ life_form + height_mm + reproduction."))
          ),
          fluidRow(
            column(3, checkboxInput("hmsc_use_traits", "Use TrData / traits", FALSE),
                   help_text("Default: off. Turn on only after uploading TrData.csv with one row per species. Enables Gamma / trait-mediated interpretation.")),
            column(3, checkboxInput("hmsc_use_phylogeny", "Use phylogeny / taxonomy", FALSE),
                   help_text("Default: off. Turn on only if the phylogeny/taxonomy file matches species names. Useful for evolutionary interpretation.")),
            column(3, checkboxInput("hmsc_save_model", "Save model object", TRUE),
                   help_text("Default: on. Saves fitted Hmsc model objects in models/.")),
            column(3, selectInput("hmsc_preset", "MCMC preset", choices = c("Quick test", "Standard", "Publication", "Custom"), selected = "Quick test"),
                   help_text("Default: Quick test. Quick test checks the workflow only. Publication analyses need longer chains and convergence checking."))
          )
        ),
        accordion_panel("4. HMSC random effects and spatial settings",
          p(class = "section-subtitle", "Combine grouping columns such as plot and year. Add a spatial level only when your study has coordinates for each spatial unit."),
          selectInput("hmsc_random_mode", "Primary random level",
            choices = c("No primary level" = "none", "Grouping" = "sample",
              "Spatial Full" = "spatial_full", "Spatial NNGP" = "spatial_nngp", "Spatial GPP" = "spatial_gpp")),
          help_text("Optional. Add grouping levels below with or without a primary spatial level. Full uses the full spatial covariance; NNGP uses neighbours; GPP uses knots."),
          conditionalPanel("input.hmsc_random_mode == 'sample'",
            textInput("hmsc_random_effect_column", "Primary grouping column", value = "sample"),
            help_text("Use the exact column name in studyDesign.csv. Example: plot.")
          ),
          conditionalPanel(
            condition = "input.hmsc_random_mode == 'spatial_full' || input.hmsc_random_mode == 'spatial_nngp' || input.hmsc_random_mode == 'spatial_gpp'",
            div(class="conditional-box",
              h4("Step 2 - Spatial coordinate settings"),
              fluidRow(
                column(4, textInput("hmsc_lon_col", "Longitude / x column", value = "longitude"),
                       help_text("Used for all spatial modes. If not found, the app auto-detects longitude/latitude, lon/lat, x/y or the first two numeric coordinate columns.")),
                column(4, textInput("hmsc_lat_col", "Latitude / y column", value = "latitude"),
                       help_text("Used for all spatial modes. If your file uses x/y, change this to y or let the app auto-detect.")),
                column(4, checkboxInput("hmsc_longlat", "longlat", FALSE),
                       help_text("Default: FALSE for projected/local coordinates. Use TRUE only when coordinates are longitude/latitude degrees."))
              )
            )
          ),
          conditionalPanel(
            condition = "input.hmsc_random_mode == 'spatial_full'",
            div(class="conditional-box",
                tags$b("Step 3 - Spatial Full Gaussian process"),
                div(class="small-muted", "No neighbor or knot parameter is needed. Use this for small spatial datasets where a full Gaussian process is computationally reasonable."))
          ),
          conditionalPanel(
            condition = "input.hmsc_random_mode == 'spatial_nngp'",
            div(class="conditional-box",
              h4("Step 3 - Spatial NNGP setting"),
              synced_numeric_slider("hmsc_nNeighbours", "NNGP nNeighbours", value = 10, min = 1, max = 1000, step = 1,
                help = "Used only for Spatial NNGP. Common values are 10-20; the upper limit is high enough for advanced tests but large values can be slow.")
            )
          ),
          conditionalPanel(
            condition = "input.hmsc_random_mode == 'spatial_gpp'",
            div(class="conditional-box",
              h4("Step 3 - Spatial GPP / knots setting"),
              textInput("hmsc_sKnot_file", "sKnot / knots file name", value = ""),
              help_text("Optional knots file for GPP. Example: knots.csv. Leave empty for first runs if the adapter can create safe knots from coordinates.")
            )
          ),
          h4("Additional grouping levels"),
          div(id = "hmsc_group_levels"),
          div(class = "studio-level-actions",
            actionButton("hmsc_add_group", "Add grouping level")),
          help_text("Each row defines a separate nonspatial random level from studyDesign.csv. Example: add plot and year for crossed effects. Use unique nested IDs, such as site1_plot2, for plots nested in sites."),
          tags$details(class = "studio-disclosure",
            tags$summary("Advanced random-level settings"),
            div(class="small-muted", "Most users should leave this block unchanged. These controls exist for older configs and custom HmscRandomLevel objects."),
            textInput("hmsc_units_column", "Advanced units column", value = "sample"),
            help_text("Defines the unit names for legacy HmscRandomLevel modes when no primary level is selected."),
            fluidRow(
              column(4, selectInput("hmsc_spatial_method", "Legacy spatial method override", choices = c("Full", "NNGP", "GPP"), selected = "NNGP"),
                     help_text("Kept for older saved configs that still say random_mode = spatial. New runs should choose Spatial Full, Spatial NNGP or Spatial GPP directly above.")),
              column(4, selectInput("hmsc_random_level_type", "Advanced HmscRandomLevel type", choices = c("none", "unstructured units", "spatial sData", "distance matrix distMat", "covariate-dependent xData", "N only"), selected = "none"),
                     help_text("Advanced. Default: none. Usually follow Random level design instead of editing this.")),
              column(4, selectInput("hmsc_sMethod", "Advanced sMethod", choices = c("Full", "NNGP", "GPP"), selected = "NNGP"),
                     help_text("Advanced spatial HmscRandomLevel setting. Default/example: NNGP."))
            ),
            fluidRow(
              column(4, synced_numeric_slider("hmsc_random_N", "N for N-only level", value = 1, min = 1, max = 1000, step = 1,
                     help = "Advanced. Default: 1. Usually not needed if using a grouping column.")),
              column(4, textInput("hmsc_distMat_file", "distMat file name", value = ""),
                     help_text("Advanced optional file name. Example: distMat.csv. Leave empty unless using distance matrix random level.")),
              column(4, textInput("hmsc_xData_file", "random-level xData file name", value = ""),
                     help_text("Advanced optional file. Example: random_xData.csv. Leave empty for ordinary models."))
            ),
            fluidRow(
              column(6, synced_numeric_slider("hmsc_nfMin", "nfMin", value = 1, min = 0, max = 1000, step = 1,
                     help = "Minimum latent factors for random-level prior. Default/example: 1. Keep default unless you understand latent-factor priors.")),
              column(6, synced_numeric_slider("hmsc_nfMax", "nfMax", value = 10, min = 1, max = 10000, step = 1,
                     help = "Maximum latent factors for random-level prior. Default/example: 10. For many species, larger nfMax may be useful but slower."))
            )
          )
        ),
        accordion_panel("5. HMSC constructor advanced options",
          fluidRow(
            column(3, checkboxInput("hmsc_XScale", "XScale", TRUE),
                   help_text("Default: TRUE. Scale environmental predictors internally. Recommended for most analyses.")),
            column(3, checkboxInput("hmsc_TrScale", "TrScale", TRUE),
                   help_text("Default: TRUE. Scale trait predictors internally. Only matters when traits are used.")),
            column(3, checkboxInput("hmsc_YScale", "YScale", FALSE),
                   help_text("Default: FALSE. Use only for selected normal-response workflows. Keep FALSE for probit/poisson.")),
            column(3, checkboxInput("hmsc_truncateNumberOfFactors", "truncateNumberOfFactors", TRUE),
                   help_text("Default: TRUE. Keeps latent-factor dimension controlled. Safe for most users."))
          ),
          fluidRow(
            column(4, textInput("hmsc_Loff_file", "Loff offset file name", value = ""),
                   help_text("Optional offset matrix. Example: offset.csv. Leave empty for ordinary analysis.")),
            column(4, textInput("hmsc_ranLevelsUsed", "ranLevelsUsed", value = ""),
                   help_text("Optional comma-separated random levels to use. Example: sample,plot. Leave empty to use defined random levels.")),
            column(4, div(class="note",
                   tags$b("Derived phylogeny matrix"),
                   div("Upload Newick or taxonomy in Step 1. JSDM Studio derives and exports derived_phylogeny_correlation.csv internally; users should not upload a species-by-species correlation matrix.")))
          ),
          checkboxInput("hmsc_use_XRRR", "Reduced-rank regression (XRRR)", FALSE),
          help_text("Optional. Model the effects of a group of predictors through a smaller set of linear combinations."),
          conditionalPanel("input.hmsc_use_XRRR",
            div(class = "studio-subsettings",
              textInput("hmsc_XRRR_file", "Reduced-rank covariate file", value = ""),
              help_text("Path to a prepared CSV, for example data/XRRRData.csv. Use site IDs in the first column and match the rows of Y."),
              textInput("hmsc_XRRRFormula", "Reduced-rank formula (XRRRFormula)", value = "~ ."),
              help_text("Choose columns from XRRRData, for example ~ temperature + moisture - 1."),
              synced_numeric_slider("hmsc_ncRRR", "Reduced-rank dimension (ncRRR)", value = 2, min = 1, max = 50, step = 1,
                help = "Number of linear combinations of the selected predictors. Example: 2."),
              checkboxInput("hmsc_XRRRScale", "Scale reduced-rank predictors (XRRRScale)", TRUE),
              help_text("Scale the covariates used in the reduced-rank component. Default: on."))
          )
        ),
        accordion_panel("6. HMSC MCMC and priors",
          div(class="warn", "Quick-test MCMC settings are for testing the software only. Publication settings must be much larger and convergence must be checked."),
          fluidRow(
            column(3, synced_numeric_slider("hmsc_samples", "samples", value = 50, min = 1, max = 20000, step = 10,
                   help = "Posterior samples retained per chain. Quick example: 50. Publication example: 1000+.")),
            column(3, synced_numeric_slider("hmsc_transient", "transient", value = 25, min = 0, max = 1000000, step = 10,
                   help = "Burn-in/transient iterations. Quick example: 25. Publication example: 5000 or more depending on thin.")),
            column(3, synced_numeric_slider("hmsc_thin", "thin", value = 1, min = 1, max = 100000, step = 1,
                   help = "Thinning interval. Quick example: 1. Publication example: 10 or more.")),
            column(3, synced_numeric_slider("hmsc_nChains", "nChains", value = 2, min = 1, max = 100, step = 1,
                   help = "Number of MCMC chains. Default/example: 2. Publication example: 4 or more."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("hmsc_nParallel", "nParallel", value = 1, min = 1, max = 256, step = 1,
                   help = "Parallel chains. Default: 1 for Windows safety. Increase only if parallel R works.")),
            column(3, synced_numeric_slider("hmsc_verbose", "verbose", value = 10, min = 0, max = 100000, step = 1,
                   help = "Progress reporting interval. Default/example: 10. Use 0 to suppress output.")),
            column(3, synced_numeric_slider("hmsc_seed", "seed", value = 1234, min = 1, max = 999999, step = 1,
                   help = "Random seed. Default/example: 1234 for reproducibility.")),
            column(3, selectInput("hmsc_initPar", "initPar", choices = c("fixed effects", "random", "current", "default"), selected = "fixed effects"),
                   help_text("Initial values strategy. Default/example: fixed effects, stable for most simple models."))
          ),
          fluidRow(
            column(3, checkboxInput("hmsc_alignPost", "alignPost", TRUE),
                   help_text("Default: TRUE. Align posterior latent variables to reduce label-switching problems.")),
            column(3, checkboxInput("hmsc_sample_prior", "sample_prior", FALSE),
                   help_text("Default: FALSE. Use TRUE only for prior predictive checks or teaching.")),
            column(3, checkboxInput("hmsc_pool_chains", "pool_chains", FALSE),
                   help_text("Default: FALSE. Keep chains separate until convergence is inspected.")),
            column(3, checkboxInput("hmsc_use_default_priors", "Use HMSC default priors", TRUE),
                   help_text("Default: TRUE. Recommended unless reproducing a specific advanced analysis."))
          ),
          fluidRow(
            column(3, checkboxInput("hmsc_updater_GammaEta", "updater$GammaEta", TRUE),
                   help_text("Default: TRUE. Random-level latent-factor updater. Ignored when no random level exists.")),
            column(3, checkboxInput("hmsc_updater_Beta", "updater$Beta", TRUE),
                   help_text("Default: TRUE. Updates environmental fixed effects Beta.")),
            column(3, checkboxInput("hmsc_updater_Gamma", "updater$Gamma", TRUE),
                   help_text("Default: TRUE. Updates trait-response Gamma when traits are used.")),
            column(3, checkboxInput("hmsc_updater_Omega", "updater$Omega", TRUE),
                   help_text("Default: TRUE. Updates association/random-level covariance parameters when present."))
          ),
          fluidRow(
            column(2, numericInput("hmsc_a1", "a1", value = NA, min = 0, step = 0.1),
                   help_text("Advanced shrinkage prior. Default/example: NA = HMSC default.")),
            column(2, numericInput("hmsc_b1", "b1", value = NA, min = 0, step = 0.1),
                   help_text("Advanced shrinkage prior. Default/example: NA = HMSC default.")),
            column(2, numericInput("hmsc_a2", "a2", value = NA, min = 0, step = 0.1),
                   help_text("Advanced shrinkage prior. Default/example: NA = HMSC default.")),
            column(2, numericInput("hmsc_b2", "b2", value = NA, min = 0, step = 0.1),
                   help_text("Advanced shrinkage prior. Default/example: NA = HMSC default.")),
            column(4, textInput("hmsc_alphapw", "alphapw", value = ""),
                   help_text("Optional spatial scale prior grid. Example: 0.1,0.5,1,2. Leave empty for default."))
          )
        ),
        accordion_panel("7. HMSC results, diagnostics and plots",
          fluidRow(
            column(3, checkboxInput("hmsc_out_predicted", "computePredictedValues", TRUE),
                   help_text("Default: TRUE. Generates predicted values for model fit and prediction outputs.")),
            column(3, checkboxInput("hmsc_out_fit", "evaluateModelFit", TRUE),
                   help_text("Default: TRUE. Computes RMSE/R2/AUC/TjurR2 depending on distr.")),
            column(3, checkboxInput("hmsc_out_cv", "Cross-validation", TRUE),
                   help_text("Default: TRUE. Uses createPartition and computePredictedValues with partitions when production adapter is connected.")),
            column(3, checkboxInput("hmsc_out_waic", "computeWAIC", FALSE),
                   help_text("Default: FALSE. WAIC can be slower. Enable when model comparison needs WAIC."))
          ),
          fluidRow(
            column(3, checkboxInput("hmsc_out_diag", "MCMC diagnostics", TRUE),
                   help_text("Default: TRUE. Uses convertToCodaObject, effectiveSize and Gelman PSRF where available.")),
            column(3, checkboxInput("hmsc_out_params", "Parameter estimates", TRUE),
                   help_text("Default: TRUE. Saves Beta, Gamma, Omega/rho/alpha summaries. Also exports book-style Beta species-environment heatmaps and Gamma trait-environment response heatmaps.")),
            column(3, checkboxInput("hmsc_out_vp", "Variance partitioning", TRUE),
                   help_text("Default: TRUE. Uses computeVariancePartitioning and exports both the native plot and a book-style species stacked-bar plus mean-predictor VP figure.")),
            column(3, checkboxInput("hmsc_out_omega", "Associations / Omega", TRUE),
                   help_text("Default: TRUE. Uses computeAssociations and exports a book-style species association heatmap. Omega is a latent residual association, not direct evidence of biotic interaction."))
          ),
          fluidRow(
            column(3, checkboxInput("hmsc_out_gradients", "Environmental gradients", FALSE),
                   help_text("Default: FALSE. Advanced: uses constructGradient/plotGradient and exports book-style response panels: summed community response, focal species response and numeric trait-weighted community response.")),
            column(3, checkboxInput("hmsc_beta_support", "Beta support", TRUE),
                   help_text("Default: TRUE. Saves support / credible interval summaries for environmental effects.")),
            column(3, checkboxInput("hmsc_gamma_support", "Gamma support", TRUE),
                   help_text("Default: TRUE. Saves support summaries for trait-mediated environment-by-trait effects when traits are used.")),
            column(3, checkboxInput("hmsc_omega_support", "Omega support", TRUE),
                   help_text("Default: TRUE. Saves support summaries for species associations."))
          ),
          h4("Convergence diagnostics"),
          fluidRow(
            column(2, checkboxInput("hmsc_conv_beta", "Beta", TRUE), help_text("Default: TRUE. Fixed-effect diagnostics.")),
            column(2, checkboxInput("hmsc_conv_gamma", "Gamma", TRUE), help_text("Default: TRUE. Trait-effect diagnostics.")),
            column(2, checkboxInput("hmsc_conv_omega", "Omega", TRUE), help_text("Default: TRUE. Association diagnostics.")),
            column(2, synced_numeric_slider("hmsc_maxOmega", "maxOmega", value = 10, min = 1, max = 100, step = 1, help = "Default/example: 10. Limits plotted Omega pairs.")),
            column(2, checkboxInput("hmsc_conv_rho", "rho", TRUE), help_text("Default: TRUE. Phylogenetic signal diagnostics if used. HMSC also exports a posterior rho curve/summary figure in results/.")),
            column(2, checkboxInput("hmsc_conv_alpha", "alpha", TRUE), help_text("Default: TRUE. Spatial scale diagnostics if spatial random level is used."))
          ),
          fluidRow(
            column(3, checkboxInput("hmsc_effective_size", "effectiveSize", TRUE), help_text("Default: TRUE. Effective sample size diagnostics.")),
            column(3, checkboxInput("hmsc_gelman_psrf", "Gelman PSRF", TRUE), help_text("Default: TRUE. Potential scale reduction factor for multiple chains.")),
            column(3, checkboxInput("hmsc_plot_beta", "plotBeta", TRUE), help_text("Default: TRUE. Produces Beta effect figures.")),
            column(3, checkboxInput("hmsc_plot_gamma", "plotGamma", TRUE), help_text("Default: TRUE. Produces native Gamma figures and a book-style trait-environment response heatmap when traits exist."))
          ),
          h4("Plotting options"),
          fluidRow(
            column(3, checkboxInput("hmsc_varpart_order_explained", "VP order explained", TRUE), help_text("Default: TRUE. Order variance partitioning by explained fraction.")),
            column(3, checkboxInput("hmsc_varpart_order_raw", "VP order raw", FALSE), help_text("Default: FALSE. Use raw order instead.")),
            column(3, checkboxInput("hmsc_show_sp_names_beta", "Show species names in Beta plots", FALSE), help_text("Default: FALSE for crowded plots.")),
            column(3, checkboxInput("hmsc_plotTree", "plotTree", FALSE), help_text("Default: FALSE. Enable when phylogeny is used."))
          ),
          fluidRow(
            column(3, selectInput("hmsc_omega_order", "Omega order", choices = c("original", "AOE", "FPC", "hclust", "alphabetical" = "alphabet"), selected = "original"), help_text("Default: original. Controls association matrix order.")),
            column(3, checkboxInput("hmsc_show_sp_names_omega", "Show species names in Omega", TRUE), help_text("Default: TRUE. Turn off for very many species.")),
            column(3, textInput("hmsc_species_list", "species.list", value = ""), help_text("Optional comma-separated species for gradient plots. Example: sp_1,sp_2. Empty = adapter default.")),
            column(3, textInput("hmsc_trait_list", "trait.list", value = ""), help_text("Optional traits for gradient interpretation. Example: height_mm,life_form."))
          ),
          fluidRow(
            column(3, textInput("hmsc_env_list", "env.list", value = ""), help_text("Optional environmental gradients. Example: pH,moisture. Empty = adapter default.")),
            column(3, synced_numeric_slider("hmsc_nfolds", "nfolds", value = 2, min = 2, max = 20, step = 1, help = "Default/example: 2 for quick CV. Publication example: 5 or 10.")),
            column(3, textInput("hmsc_partition_column", "partition column", value = "sample"), help_text("Column for grouped CV partition. Example/default: sample. Leave if no studyDesign.")),
            column(3, checkboxInput("hmsc_compute_sair", "computeSAIR", FALSE), help_text("Default: FALSE. Advanced species association interpretation routine."))
          )
        ),
        accordion_panel("8. HMSC output workflow and download", value = "results",
          div(class="note", tags$b("Output workflow follows your S1-S7 structure:"), " S1 define models, S2 fit models, S3 evaluate convergence, S4 compute model fit, S5 show model fit, S6 show parameter estimates, S7 make predictions."),
          fluidRow(
            column(3, checkboxInput("hmsc_out_copy_inputs", "Copy input files", TRUE), help_text("Default: TRUE. Copies Y, XData, TrData, studyDesign, coordinates and phylogeny into inputs/.")),
            column(3, checkboxInput("hmsc_export_scripts", "Export executable R scripts", TRUE), help_text("Default: TRUE. Saves reproducible_script/run_this_HMSC_analysis.R and workflow_scripts/S1-S7 scripts.")),
            column(3, checkboxInput("hmsc_create_standard", "Create standard comparison tables", TRUE), help_text("Default: TRUE. Exports standard comparison tables. An unavailable quantity is not a fitted estimate.")),
            column(3, checkboxInput("hmsc_real_fit", "Run real HMSC S1-S7 results workflow", TRUE), help_text("Default: TRUE. Runs the uploaded S1-S7 style workflow: define models, fit models, convergence, model fit, parameter estimates and predictions. Turn off only if you want scaffold files without fitting."))
          ),
          fluidRow(
            column(4,
              actionButton("hmsc_run", "Run HMSC workflow", class = "btn-success"),
              br(), br(),
              downloadButton("hmsc_download", "Download HMSC ZIP", icon = NULL),
              div(class="warn", "With fitting enabled, Studio runs Hmsc after validating the inputs. Turning fitting off exports reproducible files without fitted estimates.")
            ),
            column(8,
              div(class="cardx", h3("HMSC run summary"), verbatimTextOutput("hmsc_run_summary")),
              div(class="cardx", h3("HMSC log"), tags$pre(class="logbox", textOutput("hmsc_log", inline = TRUE)))
            )
          ),
          div(class="cardx", h3("HMSC output files"), DTOutput("hmsc_files"))
        )
      )
    )
  ),

  tabPanel("Hmsc-HPC", value = "Hmsc-HPC Workflow",
    fluidPage(
      workflow_overview("engine-hmschpc", "Hmsc-HPC workflow",
        "CPU TensorFlow sampling for supported Hmsc models through pyhmsc.",
        "Check the Python runtime and choose an independent-group or Full spatial random intercept when a random level is needed.",
        c("This Studio workflow supports traits and phylogeny, but does not expose all R Hmsc features or fit NNGP and GPP random levels.",
          "A minimal test uses Y.csv and XData.csv with a matching response family, no random level, and traits and phylogeny off.",
          "Associations derived from latent loadings must be interpreted at their random level; do not compare raw loadings with Omega.",
          "Compiled inputs alone are model_defined. A fitted run requires posterior sampling. Missing dependencies should be reported with diagnostics.")),
      accordion(
        class = "studio-workflow-form", open = "inputs",
        accordion_panel("1. Hmsc-HPC data upload", value = "inputs",
          fluidRow(
            column(6,
              fileInput("hmschpc_Y_file", "Y.csv - response matrix", accept = ".csv"),
              help_text("Required. Rows = sites/samples, columns = species/responses. Use 0/1 for probit, non-negative integers for poisson, numeric values for normal."),
              fileInput("hmschpc_X_file", "XData.csv - environmental predictors", accept = ".csv"),
              help_text("Required. Rows must match Y. Numeric and categorical predictors are supported through patsy design matrices."),
              fileInput("hmschpc_traits_file", "traits.csv - species traits", accept = ".csv"),
              help_text("Optional. Rows must match Y species columns. Used when Use traits is checked."),
              fileInput("hmschpc_newdata_file", "newdata.csv - prediction covariates", accept = ".csv"),
              help_text("Optional. If omitted, fitted/training predictions are exported from XData.csv.")
            ),
            column(6,
              fileInput("hmschpc_study_file", "studyDesign.csv - random-level design", accept = ".csv"),
              help_text("Optional unless iid/spatial random level is selected. Must contain the grouping column, such as plot or site."),
              fileInput("hmschpc_coord_file", "coordinates.csv - spatial coordinates", accept = ".csv"),
              help_text("Optional. Required for spatial_full unless coordinate columns already exist in studyDesign.csv."),
              fileInput("hmschpc_taxonomy_file", "taxonomy.csv - taxonomic rank table", accept = c(".csv", ".tsv")),
              help_text("Optional. Table with species plus kingdom/phylum/class/order/family/genus/species columns. JSDM Studio derives the species-by-species covariance/correlation matrix internally."),
              fileInput("hmschpc_phylo_tree_file", "phylo_tree.nwk - Newick phylogeny", accept = c(".nwk", ".tree", ".tre", ".txt")),
              help_text("Optional. Newick tree with tips matching Y species. The tree is passed to pyhmsc where supported and also summarized as a derived correlation output.")
            )
          )
        ),
        accordion_panel("2. Hmsc-HPC data check", value = "validation",
          fluidRow(
            column(4,
              actionButton("hmschpc_check", "Check Hmsc-HPC data", class = "btn-primary"),
              div(class="note", tags$b("Check before running."), " Validates dimensions, response distribution, formula variables, traits, phylogeny and supported random-level settings.")
            ),
            column(8, div(class="cardx", h3("Hmsc-HPC check messages"), verbatimTextOutput("hmschpc_check_messages")))
          ),
          div(class="cardx", h3("Hmsc-HPC data dimensions"), DTOutput("hmschpc_data_table")),
          file_preview_block("hmschpc", "Hmsc-HPC uploaded data preview")
        ),
        accordion_panel("3. Hmsc-HPC model settings", value = "settings",
          fluidRow(
            column(3, selectInput("hmschpc_distribution", "distribution", choices = c("poisson", "probit", "normal", "gaussian"), selected = "poisson"),
                   help_text("Native pyhmsc supports poisson, probit/bernoulli and normal/gaussian response matrices. The GUI normalizes gaussian to normal in model.yaml.")),
            column(3, textInput("hmschpc_XFormula", "X formula", value = "~ ."),
                   help_text("Patsy/R-like one-sided formula. Default ~ . expands to all XData columns.")),
            column(3, checkboxInput("hmschpc_use_traits", "Use traits", FALSE),
                   help_text("Enable only when traits.csv has one row per species.")),
            column(3, textInput("hmschpc_trait_formula", "trait formula", value = "~ ."),
                   help_text("Trait design formula. Default ~ . uses all trait columns."))
          ),
          fluidRow(
            column(4, selectInput("hmschpc_phylogeny_mode", "Phylogeny mode", choices = c("none", "newick", "taxonomy"), selected = "none"),
                   help_text("Choose newick for phylo_tree.nwk or taxonomy for taxonomy.csv. Do not upload a species-by-species covariance/correlation matrix; it is derived and saved internally.")),
            column(4, selectInput("hmschpc_random_mode", "Random level design", choices = c("none", "iid", "spatial_full", "random_slope_iid"), selected = "none"),
                   help_text("Native sampler supports none, iid and spatial_full. random_slope_iid is compile-only unless Run sampler is off.")),
            column(4, textInput("hmschpc_random_name", "Random level name", value = "plot"),
                   help_text("Name written to model.yaml random_levels. Example: plot, site, year."))
          ),
          fluidRow(
            column(3, textInput("hmschpc_random_column", "Grouping column", value = "plot"),
                   help_text("Column in studyDesign.csv used as the random intercept level.")),
            column(3, textInput("hmschpc_coord_x", "x coordinate column", value = "xcoord"),
                   help_text("Used for spatial_full. Can be in studyDesign.csv or coordinates.csv.")),
            column(3, textInput("hmschpc_coord_y", "y coordinate column", value = "ycoord"),
                   help_text("Used for spatial_full. Can be in studyDesign.csv or coordinates.csv.")),
            column(3, textInput("hmschpc_random_slope_formula", "random slope x_formula", value = ""),
                   help_text("Advanced compile-only branch. Example: ~ elevation. Native sampling is guarded by Hmsc-HPC."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("hmschpc_nf", "nf latent factors", value = 1, min = 1, max = 1000, step = 1,
                   help = "Latent factors for iid/spatial random levels. Small values are safest for tests.")),
            column(3, synced_numeric_slider("hmschpc_nfMin", "nfMin", value = 1, min = 1, max = 1000, step = 1,
                   help = "Minimum latent factors in the random-level prior.")),
            column(3, synced_numeric_slider("hmschpc_nfMax", "nfMax", value = 4, min = 1, max = 10000, step = 1,
                   help = "Maximum latent factors in the random-level prior.")),
            column(3, numericInput("hmschpc_alpha", "spatial alpha", value = 1, min = 0, step = 0.1),
                   help_text("Spatial range prior scale for spatial_full. Keep simple for CPU testing."))
          )
        ),
        accordion_panel("4. Hmsc-HPC CPU sampler settings",
          div(class="warn", "Quick CPU settings are only for software testing. Publication analysis needs larger samples and convergence checking."),
          fluidRow(
            column(3, checkboxInput("hmschpc_run_sampler", "Run sampler", TRUE),
                   help_text("If off, the workflow compiles and validates model artifacts only, returning model_defined.")),
            column(3, synced_numeric_slider("hmschpc_samples", "samples", value = 10, min = 1, max = 20000, step = 1,
                   help = "Posterior samples retained per chain.")),
            column(3, synced_numeric_slider("hmschpc_transient", "transient", value = 10, min = 0, max = 1000000, step = 1,
                   help = "Burn-in/transient iterations.")),
            column(3, synced_numeric_slider("hmschpc_thin", "thin", value = 1, min = 1, max = 100000, step = 1,
                   help = "Thinning interval."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("hmschpc_chains", "chains", value = 2, min = 1, max = 100, step = 1,
                   help = "Number of compiled chains.")),
            column(3, textInput("hmschpc_chains_to_run", "chains to run", value = ""),
                   help_text("Optional zero-based chain ids, comma-separated. Empty runs all compiled chains.")),
            column(3, synced_numeric_slider("hmschpc_verbose", "verbose", value = 5, min = 1, max = 100000, step = 1,
                   help = "Progress interval. Must be at least 1; Hmsc-HPC verbose=0 currently fails in TensorFlow.")),
            column(3, synced_numeric_slider("hmschpc_seed", "rngseed", value = 1234, min = 0, max = 999999, step = 1,
                   help = "Random seed passed to hmsc.run_gibbs_sampler."))
          ),
          fluidRow(
            column(3, selectInput("hmschpc_precision", "floating precision", choices = c("64", "32"), selected = "64"),
                   help_text("fp64 is default and conservative; fp32 can be faster but less stable.")),
            column(3, selectInput("hmschpc_tnlib", "truncated normal library", choices = c("tf", "tfd", "scipy"), selected = "tf"),
                   help_text("Hmsc-HPC truncated-normal sampler backend.")),
            column(3, synced_numeric_slider("hmschpc_hmcleapfrog", "hmcleapfrog", value = 10, min = 1, max = 1000, step = 1,
                   help = "HMC conditional updater leapfrog steps.")),
            column(3, synced_numeric_slider("hmschpc_hmcthin", "hmcthin", value = 0, min = 0, max = 100000, step = 1,
                   help = "Iterations between HMC conditional updater calls; 0 disables HMC updater."))
          ),
          fluidRow(
            column(3, checkboxInput("hmschpc_update_beta_eta", "updbe", FALSE),
                   help_text("Joint Beta/Eta update flag. Keep off for first CPU tests.")),
            column(3, checkboxInput("hmschpc_save_eta", "save Eta", TRUE),
                   help_text("Requested Eta export. The CPU HDF5 runner may force this on because the current Hmsc-HPC writer expects Eta entries.")),
            column(3, checkboxInput("hmschpc_eager", "TensorFlow eager", FALSE),
                   help_text("Debug mode. Slower; graph mode is normal.")),
            column(3, checkboxInput("hmschpc_profile", "TensorFlow profile", FALSE),
                   help_text("Advanced profiling flag. Usually off."))
          ),
          fluidRow(
            column(6, textInput("hmschpc_python", "Python executable", value = hmschpc_python_bin(list())),
                   help_text("Default uses the detected Python executable. The checked environment here used Python 3.10.")),
            column(6, textInput("hmschpc_python_source", "Hmsc-HPC source directory", value = file.path(app_dir, "external_packages", "hmsc-hpc-main")),
                   help_text("Bundled source directory used through PYTHONPATH for reproducible runs."))
          )
        ),
        accordion_panel("5. Hmsc-HPC outputs and download", value = "results",
          fluidRow(
            column(3, checkboxInput("hmschpc_out_predictions", "Predictions", TRUE),
                   help_text("Exports fitted/training predictions and standard predictions_long.csv.")),
            column(3, checkboxInput("hmschpc_out_diagnostics", "Diagnostics", TRUE),
                   help_text("Exports validate-init log, sampler log, Rhat/ESS when arviz is available.")),
            column(3, checkboxInput("hmschpc_out_plots", "Plots", TRUE),
                   help_text("Exports Beta heatmap and observed-vs-fitted plot when matplotlib is available.")),
            column(3, checkboxInput("hmschpc_out_zip", "Create ZIP", TRUE),
                   help_text("Creates a real ZIP of the Hmsc-HPC output folder."))
          ),
          fluidRow(
            column(4,
              actionButton("hmschpc_run", "Run Hmsc-HPC workflow", class = "btn-success"),
              br(), br(),
              downloadButton("hmschpc_download", "Download Hmsc-HPC ZIP", icon = NULL),
              div(class="warn", "The CPU workflow writes executable R and Python scripts plus pyhmsc JSON/HDF5 artifacts.")
            ),
            column(8,
              div(class="cardx", h3("Hmsc-HPC run summary"), verbatimTextOutput("hmschpc_run_summary")),
              div(class="cardx", h3("Hmsc-HPC log"), tags$pre(class="logbox", textOutput("hmschpc_log", inline = TRUE)))
            )
          ),
          div(class="cardx", h3("Hmsc-HPC output files"), DTOutput("hmschpc_files"))
        )
      )
    )
  ),

  tabPanel("jSDM", value = "jSDM Workflow",
    fluidPage(
      workflow_overview("engine-jsdm", "jSDM workflow",
        "Bayesian community regression with optional latent factors, site effects and trait terms.",
        "Choose a function that matches the binary, count or continuous response in Y.",
        c("Which latent factors and site effects are fitted determines the meaning of residual associations.",
          "jSDM trait coefficients are model-specific. Do not equate them with Hmsc Gamma.",
          "jSDM separates environmental and residual correlations. Neither should be treated as the same estimand as Hmsc Omega or boral residual correlation without checking model definitions.",
          "Match predictions by site and species IDs before comparing them.")),
      accordion(
        class = "studio-workflow-form", open = "inputs",
        accordion_panel("1. jSDM data upload", value = "inputs",
          fluidRow(
            column(6,
              fileInput("jsdm_Y_file", "Y.csv  - jSDM response matrix", accept = ".csv"),
              help_text("Rows = sites/sampling units. Columns = species/responses. Binary data suit binomial models; non-negative integer data suit poisson; continuous data suit gaussian."),
              fileInput("jsdm_X_file", "XData.csv  - jSDM site_data", accept = ".csv"),
              help_text("Rows must match Y. Columns are predictors used by site_formula."),
              fileInput("jsdm_trait_file", "trait_data.csv  - optional traits", accept = ".csv"),
              help_text("Optional. Rows should match Y columns when trait_formula is used.")
            ),
            column(6,
              fileInput("jsdm_long_file", "long_format.csv  - optional long-format data", accept = ".csv"),
              help_text("Optional. Used only for jSDM_binomial_probit_long_format workflows."),
              fileInput("jsdm_trials_file", "trials.csv  - optional binomial trials", accept = ".csv"),
              help_text("Optional. For binomial logit models with repeated trials. It can be a vector/table with one value per site; each trial count should be >= observed successes. Use 1 for simple presence/absence."),
              fileInput("jsdm_newdata_file", "newdata.csv  - optional prediction covariates", accept = ".csv"),
              help_text("Optional. New site covariates for predict.jSDM. Columns must match the variables used in site_formula/site_data."),
              fileInput("jsdm_prediction_ids_file", "prediction_ids.csv  - optional site/species IDs", accept = ".csv"),
              help_text("Optional. A helper table for selecting Id_sites and Id_species for prediction. If omitted, the adapter can use all species and selected/default sites."),
              div(class="note", tags$b("jSDM data are independent from Hmsc data."), " Uploading Hmsc files does not automatically run jSDM.")
            )
          )
        ),
        accordion_panel("2. jSDM data check", value = "validation",
          actionButton("jsdm_check", "Check jSDM data", class = "btn-primary"),
          br(), br(),
          fluidRow(
            column(5, div(class="cardx", h3("jSDM check messages"), verbatimTextOutput("jsdm_check_messages"))),
            column(7, div(class="cardx", h3("jSDM data dimensions"), DTOutput("jsdm_data_table")))
          ),
          file_preview_block("jsdm", "jSDM uploaded data preview")
        ),
        accordion_panel("3. jSDM model settings", value = "settings",
          fluidRow(
            column(4,
              selectInput("jsdm_model_type", "model_type",
                choices = c("binomial_probit", "binomial_logit", "poisson_log", "gaussian", "binomial_probit_long_format", "binomial_probit_sp_constrained"),
                selected = "binomial_probit"),
              help_text("Choose the jSDM function family. This is not the same as Hmsc distr.")
            ),
            column(4,
              textInput("jsdm_site_formula", "site_formula", value = "~ ."),
              help_text("One-sided formula for site-level predictors in jSDM.")
            ),
            column(4,
              textInput("jsdm_trait_formula", "trait_formula", value = "~ ."),
              help_text("Trait formula if trait_data is supplied.")
            )
          ),
          fluidRow(
            column(4, synced_numeric_slider("jsdm_n_latent", "n_latent", value = 2, min = 0, max = 20, step = 1,
                   help = "Number of latent variables. Residual correlations require latent variables.")),
            column(4, selectInput("jsdm_site_effect", "site_effect", choices = c("none", "fixed", "random"), selected = "random"),
                   help_text("jSDM site effect. This is not identical to Hmsc studyDesign random levels.")),
            column(4, synced_numeric_slider("jsdm_trials", "trials", value = 1, min = 1, max = 1000, step = 1,
                   help = "Number of trials for binomial logit. Use 1 for presence/absence."))
          ),
          hr(),
          h4("jSDM MCMC controls"),
          fluidRow(
            column(3, synced_numeric_slider("jsdm_burnin", "burnin", value = 100, min = 0, max = 50000, step = 100,
                   help = "Burn-in iterations for jSDM.")),
            column(3, synced_numeric_slider("jsdm_mcmc", "mcmc", value = 100, min = 100, max = 100000, step = 100,
                   help = "Gibbs iterations for jSDM. Increase for real analysis.")),
            column(3, synced_numeric_slider("jsdm_thin", "thin", value = 1, min = 1, max = 500, step = 1,
                   help = "Thinning interval.")),
            column(3, synced_numeric_slider("jsdm_seed", "seed", value = 1234, min = 1, max = 999999, step = 1,
              help = "Random seed supplied to the sampler. Example: 1234. Keep the seed with the package version and settings when reproducing a run."))
          ),
          fluidRow(
            column(4, synced_numeric_slider("jsdm_verbose", "verbose", value = 1, min = 0, max = 1, step = 1,
              help = "Controls progress in the sampler log, without changing the fitted model. 1 prints progress; 0 is quiet.")),
            column(4, selectInput("jsdm_preset", "MCMC preset", choices = c("Quick test", "Standard", "Publication", "Custom"), selected = "Quick test"),
              help_text("Sets burn-in, sampling iterations and thinning. Quick test checks execution; Publication is a longer starting point and requires convergence checks. Custom leaves your values unchanged.")),
            column(4, synced_numeric_slider("jsdm_ropt", "ropt", value = 0.44, min = 0.01, max = 0.99, step = 0.01,
                   help = "Target acceptance rate for adaptive Metropolis where relevant."))
          ),
          hr(),
          h4("jSDM starting values"),
          fluidRow(
            column(3, numericInput("jsdm_beta_start", "beta_start", value = 0),
              help_text("Initial species-by-predictor regression coefficients. A scalar is repeated for every coefficient. Example: 0 starts every environmental slope at zero.")),
            column(3, numericInput("jsdm_gamma_start", "gamma_start", value = 0),
              help_text("Initial coefficients linking traits to species-level regression coefficients. Used only when the selected function fits the supplied traits. Example: 0 starts with no trait modulation.")),
            column(3, numericInput("jsdm_lambda_start", "lambda_start", value = 0),
              help_text("Initial species loadings on latent factors; ignored with zero latent variables. The package imposes triangular identification constraints. Example: leave the scalar at 0 for its default initialization.")),
            column(3, numericInput("jsdm_W_start", "W_start", value = 0),
              help_text("Initial site scores on each latent factor; ignored with zero latent variables. Example: 0 starts scores at zero."))
          ),
          fluidRow(
            column(3, numericInput("jsdm_alpha_start", "alpha_start", value = 0),
              help_text("Initial site intercept effects; ignored when site effects are disabled. Example: 0.")),
            column(3, numericInput("jsdm_V_alpha", "V_alpha", value = 1, min = 0.0001),
              help_text("Positive initial variance for random site effects, or fixed Gaussian prior variance for fixed site effects. Ignored with no site effect. Example: 1.")),
            column(3, numericInput("jsdm_V_start", "V_start", value = 1, min = 0.0001),
              help_text("Positive starting value for the residual variance. Gaussian workflow only. Example: 1.")),
            column(3, div(class="safe", "Scalar starts are safest in GUI mode."))
          ),
          hr(),
          h4("jSDM priors"),
          fluidRow(
            column(3, numericInput("jsdm_shape_Valpha", "shape_Valpha", value = 0.5, min = 0.0001),
              help_text("Shape of the inverse-Gamma prior on random site-effect variance. Used only with random site effects. Example: 0.5.")),
            column(3, numericInput("jsdm_rate_Valpha", "rate_Valpha", value = 0.0005, min = 0.000001),
              help_text("Rate argument of the inverse-Gamma prior on random site-effect variance. Used only with random site effects. Example: 0.0005.")),
            column(3, numericInput("jsdm_shape_V", "shape_V", value = 0.5, min = 0.0001),
              help_text("Shape of the inverse-Gamma prior on Gaussian residual variance. Gaussian workflow only. Example: 0.5.")),
            column(3, numericInput("jsdm_rate_V", "rate_V", value = 0.0005, min = 0.000001),
              help_text("Rate argument of the inverse-Gamma prior on Gaussian residual variance. Gaussian workflow only. Example: 0.0005."))
          ),
          fluidRow(
            column(3, numericInput("jsdm_mu_beta", "mu_beta", value = 0),
              help_text("Mean of the Normal prior for species regression coefficients. The scalar applies to all supported coefficients. Example: 0.")),
            column(3, numericInput("jsdm_V_beta", "V_beta", value = 10, min = 0.0001),
              help_text("Variance of the Normal prior for regression coefficients. A positive scalar supplies a common diagonal variance, not a standard deviation. Example: 10.")),
            column(3, numericInput("jsdm_mu_gamma", "mu_gamma", value = 0),
              help_text("Mean of the Normal prior for the coefficients relating traits to species responses; applies to supported trait models. Example: 0.")),
            column(3, numericInput("jsdm_V_gamma", "V_gamma", value = 10, min = 0.0001),
              help_text("Positive variance of the Normal prior for trait-mediated coefficients; applies to supported trait models. Example: 10."))
          ),
          fluidRow(
            column(3, numericInput("jsdm_mu_lambda", "mu_lambda", value = 0),
              help_text("Mean of the Normal prior on free latent-factor loadings. Used with latent factors. Example: 0.")),
            column(3, numericInput("jsdm_V_lambda", "V_lambda", value = 10, min = 0.0001),
              help_text("Variance of the Normal prior on free latent-factor loadings. A positive scalar specifies a common diagonal variance. Example: 10.")),
            column(6, help_text("These priors are jSDM-specific and should not be interpreted as Hmsc priors."))
          )
        ),

        accordion_panel("4. jSDM function-specific settings",
          div(class="note", tags$b("Different jSDM functions require different response objects."), " The GUI keeps one jSDM workflow but records which function-specific options are active."),
          fluidRow(
            column(4, selectInput("jsdm_response_argument", "response argument mapping",
              choices = c("auto from model_type", "presence_data", "count_data", "response_data", "long_format_data"),
              selected = "auto from model_type"),
              help_text("Which jSDM function argument should receive Y. binomial models use presence_data; poisson uses count_data; gaussian uses response_data; long-format uses long-format data.")),
            column(4, checkboxInput("jsdm_constrained_latent", "Use constrained probit latent model", FALSE),
              help_text("For jSDM_binomial_probit_sp_constrained. Requires n_latent > 0 and presence-absence data. It constrains selected lambda loadings to improve convergence of latent axes.")),
            column(4, synced_numeric_slider("jsdm_constrained_nchains", "constrained preliminary chains", value = 2, min = 1, max = 8, step = 1,
              help = "Number of chains/runs to use in the preliminary step for identifying constrained species/loadings, when the constrained probit model is implemented."))
          ),
          fluidRow(
            column(4, textInput("jsdm_long_site_col", "long format site column", value = "site"),
              help_text("For long-format probit data: column identifying the site/sampling unit.")),
            column(4, textInput("jsdm_long_species_col", "long format species column", value = "species"),
              help_text("For long-format probit data: column identifying the species/response.")),
            column(4, textInput("jsdm_long_response_col", "long format response column", value = "presence"),
              help_text("For long-format probit data: binary response column, usually 0/1."))
          ),
          fluidRow(
            column(4, checkboxInput("jsdm_scale_site_data", "Scale continuous site_data", TRUE),
              help_text("Recommended for many jSDM examples: continuous covariates are often centered/scaled before modelling. Categorical/binary variables should not be blindly scaled.")),
            column(4, checkboxInput("jsdm_include_intercept", "Include intercept in site_formula", TRUE),
              help_text("A one-sided formula normally includes an intercept unless removed with ~ 0 + x. Keep TRUE for ordinary models.")),
            column(4, checkboxInput("jsdm_allow_traits", "Allow trait_data / trait_formula", TRUE),
              help_text("If trait_data is supplied, jSDM uses gamma parameters to model how traits explain species-specific beta responses."))
          )
        ),
        accordion_panel("5. jSDM prediction settings",
          div(class="note", tags$b("predict.jSDM outputs can be large."), " Mean predictions are safest. Quantile and posterior predictions can become memory-heavy."),
          fluidRow(
            column(3, checkboxInput("jsdm_do_predict", "Run predict.jSDM", TRUE),
              help_text("Generate predictions from the fitted jSDM object when a production adapter is connected.")),
            column(3, selectInput("jsdm_predict_type", "prediction type", choices = c("mean", "quantile", "posterior"), selected = "mean"),
              help_text("mean = posterior mean prediction; quantile = mean plus quantiles; posterior = full predictive posterior and may be very large.")),
            column(3, textInput("jsdm_predict_probs", "probs for quantile", value = "0.025,0.5,0.975"),
              help_text("Comma-separated probabilities used only when prediction type is quantile.")),
            column(3, synced_numeric_slider("jsdm_predict_max_sites", "max prediction sites", value = 500, min = 1, max = 50000, step = 50,
              help = "Safety limit for GUI prediction to prevent huge prediction objects."))
          ),
          fluidRow(
            column(4, textInput("jsdm_Id_sites", "Id_sites", value = "all"),
              help_text("Sites for prediction. Use 'all', a comma-separated list of row names/indices, or provide prediction_ids.csv.")),
            column(4, textInput("jsdm_Id_species", "Id_species", value = "all"),
              help_text("Species/responses for prediction. Use 'all', a comma-separated list of names/indices, or provide prediction_ids.csv.")),
            column(4, checkboxInput("jsdm_prediction_histograms", "Prediction histograms", TRUE),
              help_text("Create histograms of predicted theta/probabilities/count means/continuous predictions, depending on model type."))
          )
        ),
        accordion_panel("6. jSDM correlation and diagnostic settings",
          fluidRow(
            column(3, synced_numeric_slider("jsdm_cor_prob", "correlation HPD probability", value = 0.95, min = 0.5, max = 0.999, step = 0.01,
              help = "Probability coverage for residual/environmental correlation intervals. 0.95 is the default in helper functions.")),
            column(3, selectInput("jsdm_cor_type", "correlation estimate type", choices = c("mean", "median"), selected = "mean"),
              help_text("Use posterior mean or posterior median as the point estimate for correlation calculations.")),
            column(3, checkboxInput("jsdm_plot_residual_cor", "plot_residual_cor", TRUE),
              help_text("Create a corrplot-style residual correlation figure when n_latent > 0.")),
            column(3, checkboxInput("jsdm_plot_associations", "plot_associations", TRUE),
              help_text("Create association/network-style plots when supported."))
          ),
          fluidRow(
            column(3, checkboxInput("jsdm_diag_beta", "Trace/density beta", TRUE),
              help_text("Trace and density plots for species-specific beta parameters in mcmc.sp.")),
            column(3, checkboxInput("jsdm_diag_lambda", "Trace/density lambda", TRUE),
              help_text("Trace and density plots for latent-variable loadings lambda when n_latent > 0.")),
            column(3, checkboxInput("jsdm_diag_W", "Trace/density W", TRUE),
              help_text("Trace/density plots for latent variables W for selected sites when n_latent > 0.")),
            column(3, checkboxInput("jsdm_diag_alpha", "Trace/density alpha", TRUE),
              help_text("Trace/density plots for site effects alpha when site_effect is fixed or random."))
          ),
          fluidRow(
            column(3, checkboxInput("jsdm_diag_Valpha", "Trace/density V_alpha", TRUE),
              help_text("Trace/density plots for random site-effect variance V_alpha when site_effect='random'.")),
            column(3, checkboxInput("jsdm_diag_V", "Trace/density residual V", TRUE),
              help_text("For gaussian models, trace/density plots for residual variance V.")),
            column(3, checkboxInput("jsdm_diag_deviance", "Trace/density deviance", TRUE),
              help_text("Trace/density plot and summary of mcmc.Deviance.")),
            column(3, checkboxInput("jsdm_coda_summary", "coda summaries", TRUE),
              help_text("Save coda summaries for available mcmc objects: beta/lambda, gamma, alpha, V_alpha, V and Deviance."))
          )
        ),
        accordion_panel("7. jSDM output table settings",
          fluidRow(
            column(4, checkboxInput("jsdm_table_model_spec", "model_spec table", TRUE),
              help_text("Save model_spec components: family/link, number of latent variables, site effect, data dimensions, prior and MCMC settings.")),
            column(4, checkboxInput("jsdm_table_beta", "beta summary table", TRUE),
              help_text("Save posterior summaries of beta parameters from mcmc.sp.")),
            column(4, checkboxInput("jsdm_table_lambda", "lambda summary table", TRUE),
              help_text("Save posterior summaries of lambda factor loadings when n_latent > 0."))
          ),
          fluidRow(
            column(4, checkboxInput("jsdm_table_gamma", "gamma summary table", TRUE),
              help_text("Save posterior summaries of gamma trait-effect parameters when trait_data is used.")),
            column(4, checkboxInput("jsdm_table_alpha", "alpha summary table", TRUE),
              help_text("Save posterior summaries of site effects alpha when site_effect is not none.")),
            column(4, checkboxInput("jsdm_table_predictions", "prediction table", TRUE),
              help_text("Save prediction outputs from predict.jSDM as CSV/RDS depending on size and type."))
          )
        ),
        accordion_panel("8. jSDM output settings",
          div(class="note", tags$b("jSDM output is separated from Hmsc output."), " The downloaded jSDM ZIP has its own folders and report."),
          fluidRow(
            column(4, checkboxInput("jsdm_out_resid", "Residual correlations", TRUE),
                   help_text("Compute get_residual_cor outputs when n_latent > 0: covariance/correlation matrices, lower/upper HPD intervals and significant correlation matrices.")),
            column(4, checkboxInput("jsdm_out_env", "Environmental correlations", TRUE),
                   help_text("Compute get_enviro_cor outputs: correlations due to shared environmental responses to X beta.")),
            column(4, checkboxInput("jsdm_out_trace", "Traceplots and density plots", TRUE),
                   help_text("Generate PDF plots for beta, lambda, W, alpha, V_alpha, V and Deviance when those objects exist."))
          ),
          fluidRow(
            column(4, checkboxInput("jsdm_out_pred", "Predictions", TRUE),
                   help_text("Save predicted theta / response values from predict.jSDM or fitted prediction slots such as theta_latent, logit_theta_latent, probit_theta_latent or Y_pred.")),
            column(4, checkboxInput("jsdm_out_model", "Save jSDM model object", TRUE),
                   help_text("Save the fitted jSDM object as models/jsdm_model.rds when full fitting is enabled.")),
            column(4, checkboxInput("jsdm_out_report", "Generate jSDM HTML report", TRUE),
                   help_text("Create report/jSDM_report.html summarizing data, settings, diagnostics and output files."))
          ),
          fluidRow(
            column(4, checkboxInput("jsdm_out_mcmc_rds", "Save MCMC objects as RDS", TRUE),
                   help_text("Save raw MCMC components such as mcmc.sp, mcmc.gamma, mcmc.latent, mcmc.alpha, mcmc.V_alpha, mcmc.V and mcmc.Deviance when available.")),
            column(4, checkboxInput("jsdm_out_csv_tables", "Save CSV summary tables", TRUE),
                   help_text("Save compact posterior means, medians, HPD intervals and model metadata as CSV tables.")),
            column(4, checkboxInput("jsdm_out_figures", "Save figure PDFs/PNGs", TRUE),
                   help_text("Save traceplots, density plots, correlation plots and prediction histograms in plots/."))
          ),
          hr(),
          checkboxInput("jsdm_real_fit", "Run real jSDM package fit", TRUE),
          help_text("Default: TRUE. Calls the installed jSDM package and exports an executable reproducible script. Turn off only when you intentionally want the output scaffold without fitting."),
          checkboxInput("jsdm_out_copy_inputs", "Copy jSDM input files to output folder", TRUE),
          help_text("Retain original uploads for an export-only run. A fitted analysis always saves its required inputs for replay."),
          checkboxInput("jsdm_out_config", "Save jSDM used_config.yml", TRUE),
          help_text("Save settings with an export-only run. Fitting always records the configuration it used."),
          checkboxInput("jsdm_out_zip", "Create jSDM output ZIP", TRUE),
          help_text("Package the run folder and its diagnostics as a single download. Example: keep enabled for sharing or replay.")
        ),
        accordion_panel("9. jSDM run and results", value = "results",
          fluidRow(
            column(4,
              actionButton("jsdm_run", "Run jSDM workflow", class = "btn-success"),
              br(), br(),
              downloadButton("jsdm_download", "Download jSDM ZIP", icon = NULL),
              div(class="warn", "With real fit enabled, this workflow runs the installed jSDM package, exports an executable R script, and writes diagnostics if fitting fails.")
            ),
            column(8,
              div(class="cardx", h3("jSDM run summary"), verbatimTextOutput("jsdm_run_summary")),
              div(class="cardx", h3("jSDM log"), tags$pre(class="logbox", textOutput("jsdm_log", inline = TRUE)))
            )
          ),
          div(class="cardx", h3("jSDM output files"), DTOutput("jsdm_files"))
        )
      )
    )
  ),


  tabPanel("GJAM", value = "GJAM Workflow",
    fluidPage(
      workflow_overview("engine-gjam", "GJAM workflow",
        "Joint regression for responses that can differ in measurement scale or observation process.",
        "For mixed responses, specify the correct typeName for every response column.",
        c("PA, CON, CA, DA, FC, CC, OC and CAT specify different response types. A successful fit does not validate those choices.",
          "Censoring, effort and compositional groups affect the observation model. Document them where applicable.",
          "Check the output scale before comparing GJAM summaries with link-scale estimates.",
          "The inverse-prediction output reconstructs covariates for the fitted response records.")),
      accordion(
        class = "studio-workflow-form", open = "inputs",
        accordion_panel("1. GJAM data upload", value = "inputs",
          fluidRow(
            column(6,
              fileInput("gjam_Y_file", "Y.csv  - GJAM response matrix / ydata", accept = ".csv"),
              help_text("Required. Rows = observations/sampling units. Columns = response variables such as species abundance, presence/absence, ordinal scores, composition components, categorical attributes, biomass, basal area, microbiome OTUs or ecological condition variables."),
              fileInput("gjam_X_file", "XData.csv  - GJAM predictors / xdata", accept = ".csv"),
              help_text("Required. Rows must match Y. Columns are predictors used in the formula. Factors should be coded as factor/categorical columns in the production adapter."),
              fileInput("gjam_type_file", "typeNames.csv  - response type table", accept = ".csv"),
              help_text("Recommended. One row per Y column with response name and typeName. Allowed types: PA, CON, CA, DA, FC, CC, OC, CAT. If omitted, the single selected typeName is applied to all response columns.")
            ),
            column(6,
              fileInput("gjam_censor_file", "censor.csv  - optional censoring intervals", accept = ".csv"),
              help_text("Optional. Defines censored values/intervals for gjamCensorY. Use for detection limits, interval observations or custom censoring."),
              fileInput("gjam_effort_file", "effort.csv  - optional sampling effort", accept = ".csv"),
              help_text("Optional. Effort can be plot area, search time, sequencing depth or sampling effort, especially for discrete abundance DA."),
              fileInput("gjam_newdata_file", "newdata.csv  - optional prediction predictors", accept = ".csv"),
              help_text("Optional. New predictor table for gjamPredict. Columns should match variables in formula/xdata.")
            )
          ),
          fluidRow(
            column(4,
              fileInput("gjam_trait_spec_file", "specByTrait.csv  - optional species-by-trait table", accept = ".csv"),
              help_text("Optional. Used by traitList / gjamSpec2Trait to translate species responses into plot-by-trait summaries.")
            ),
            column(4,
              fileInput("gjam_trait_types_file", "traitTypes.csv  - optional trait type table", accept = ".csv"),
              help_text("Optional. Describes trait types used in trait analysis.")
            ),
            column(4,
              fileInput("gjam_holdout_file", "holdoutIndex.csv  - optional holdout rows", accept = ".csv"),
              help_text("Optional. Explicit row indices for out-of-sample prediction. If omitted, use holdoutN.")
            )
          )
        ),
        accordion_panel("2. GJAM data check", value = "validation",
          actionButton("gjam_check", "Check GJAM data", class = "btn-primary"),
          br(), br(),
          fluidRow(
            column(5, div(class="cardx", h3("GJAM check messages"), verbatimTextOutput("gjam_check_messages"))),
            column(7, div(class="cardx", h3("GJAM data dimensions"), DTOutput("gjam_data_table")))
          ),
          file_preview_block("gjam", "GJAM uploaded data preview"),
          div(class="note", tags$b("GJAM compatibility check:"), " Y and XData must have the same number of rows. typeNames must be either one value for all responses or one type per Y column. Composition groups must match FC/CC response columns.")
        ),
        accordion_panel("3. GJAM model settings", value = "settings",
          fluidRow(
            column(4, textInput("gjam_formula", "formula", value = "~ ."),
                   help_text("One-sided formula beginning with ~, for example ~ temp + moisture + treatment. The response is not written on the left side because ydata is supplied separately.")),
            column(4, selectInput("gjam_type_single", "single typeName", choices = c("PA", "CON", "CA", "DA", "FC", "CC", "OC", "CAT"), selected = "DA"),
                   help_text("Used when no typeNames.csv is supplied. PA=presence/absence; CON=continuous; CA=continuous abundance with zero censoring; DA=discrete abundance/count; FC=fractional composition; CC=count composition; OC=ordinal counts; CAT=categorical classes.")),
            column(4, textInput("gjam_notStandard", "notStandard", value = ""),
                   help_text("Comma-separated xdata column names that should not be standardized. GJAM standardizes non-factor predictors by default, then transforms output back. Use for variables that must remain on the original computational scale."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("gjam_ng", "ng", value = 2000, min = 100, max = 100000, step = 500,
                   help = "Number of Gibbs steps. Default in the package manual is 2000. Increase for formal analysis.")),
            column(3, synced_numeric_slider("gjam_burnin", "burnin", value = 500, min = 0, max = 50000, step = 100,
                   help = "Initial Gibbs steps discarded. Must be less than ng.")),
            column(3, synced_numeric_slider("gjam_holdoutN", "holdoutN", value = 0, min = 0, max = 10000, step = 1,
                   help = "Number of observations/rows held out for out-of-sample prediction if holdoutIndex is not supplied.")),
            column(3, synced_numeric_slider("gjam_seed", "seed", value = 1234, min = 1, max = 999999, step = 1,
                   help = "Random seed used before fitting for reproducibility."))
          ),
          fluidRow(
            column(4, textInput("gjam_random", "random", value = ""),
                   help_text("Optional xdata column name for random effects. The column should be a factor and each group should have replication.")),
            column(4, checkboxInput("gjam_FULL", "FULL", FALSE),
                   help_text("If TRUE, save full prediction chains in chains$ygibbs. This can be large.")),
            column(4, checkboxInput("gjam_PREDICTX", "PREDICTX", TRUE),
                   help_text("If TRUE, perform inverse prediction of predictors X from Y. Set FALSE to speed up large runs."))
          ),
          fluidRow(
            column(4, checkboxInput("gjam_REDUCT", "REDUCT", FALSE),
                   help_text("If TRUE, force dimension reduction. If FALSE, allow ordinary covariance unless automatic reduction is triggered by many responses.")),
            column(4, synced_numeric_slider("gjam_reduct_N", "reductList$N", value = 20, min = 1, max = 200, step = 1,
                   help = "Dimension-reduction parameter N. Used in reductList = list(N=..., r=...).")),
            column(4, synced_numeric_slider("gjam_reduct_r", "reductList$r", value = 3, min = 1, max = 50, step = 1,
                   help = "Dimension-reduction rank r. Larger values allow more covariance complexity but increase cost."))
          ),
          fluidRow(
            column(4, synced_numeric_slider("gjam_ematAlpha", "ematAlpha", value = 0.5, min = 0.01, max = 0.99, step = 0.01,
                   help = "Probability assigned for conditional and marginal independence in the environmental response matrix ematrix.")),
            column(4, checkboxInput("gjam_USE_CENSOR", "Use censor list", FALSE),
                   help_text("Use censor.csv to create a censor list via gjamCensorY in the production adapter.")),
            column(4, checkboxInput("gjam_USE_EFFORT", "Use effort", FALSE),
                   help_text("Use effort.csv for effort-adjusted discrete abundance models or sampling-effort adjustments."))
          )
        ),
        accordion_panel("4. GJAM response type and composition settings",
          fluidRow(
            column(4, textAreaInput("gjam_typeNames_text", "typeNames manual override", value = "", rows = 3),
                   help_text("Optional comma-separated response types, one per Y column. Leave empty to use typeNames.csv or the single typeName.")),
            column(4, textAreaInput("gjam_FCgroups", "FCgroups", value = "", rows = 3),
                   help_text("For fractional composition responses. Provide comma-separated group IDs of length S, where non-FC columns use 0. Example: 0,0,1,1,1,2,2,2.")),
            column(4, textAreaInput("gjam_CCgroups", "CCgroups", value = "", rows = 3),
                   help_text("For count composition responses. Provide comma-separated group IDs of length S, where non-CC columns use 0."))
          ),
          fluidRow(
            column(4, selectInput("gjam_composition_reference", "composition reference rule", choices = c("last column in each group", "first column in each group", "adapter default"), selected = "adapter default"),
                   help_text("Composition data have a reference component. The production adapter should document which column is treated as reference.")),
            column(4, checkboxInput("gjam_trimY", "Run gjamTrimY before fitting", FALSE),
                   help_text("Trim response matrix and/or aggregate rare response types before fitting. Useful for very sparse community matrices.")),
            column(4, synced_numeric_slider("gjam_trim_minObs", "minimum observations for trimming", value = 5, min = 1, max = 1000, step = 1,
                   help = "Responses observed fewer than this threshold may be dropped or aggregated by gjamTrimY in the adapter."))
          ),
          div(class="warn", "typeNames is the most important GJAM setting. Wrong response types will make interpretation wrong even if the model runs.")
        ),
        accordion_panel("5. GJAM prior and censoring settings",
          fluidRow(
            column(4, checkboxInput("gjam_use_prior_template", "Use gjamPriorTemplate", FALSE),
                   help_text("Create or load a prior template for coefficient constraints. Useful when known signs or bounded effects are required.")),
            column(4, fileInput("gjam_prior_file", "priorTemplate.csv / prior file", accept = ".csv"),
                   help_text("Optional prior template for coefficients. In the production adapter this can be converted to the prior list expected by GJAM.")),
            column(4, selectInput("gjam_prior_mode", "prior mode", choices = c("non-informative/default", "sign-constrained", "custom file"), selected = "non-informative/default"),
                   help_text("Default GJAM examples commonly use weak/non-informative priors. Sign-constrained priors can be useful when effect direction is known."))
          ),
          fluidRow(
            column(4, textInput("gjam_censor_columns", "censor columns", value = ""),
                   help_text("Optional comma-separated columns in Y affected by censoring. Can also be supplied in censor.csv.")),
            column(4, textInput("gjam_censor_values", "censor values", value = ""),
                   help_text("Observed values or labels that should be treated as censored. Can also be supplied in censor.csv.")),
            column(4, textInput("gjam_censor_intervals", "censor intervals", value = ""),
                   help_text("Intervals defining censored latent values. Use censor.csv for complex interval tables."))
          )
        ),
        accordion_panel("6. GJAM prediction, sensitivity and ordination settings",
          fluidRow(
            column(4, checkboxInput("gjam_do_predict", "Run gjamPredict", TRUE),
                   help_text("Generate in-sample and/or out-of-sample predictions, including conditional prediction when newdata or holdouts are provided.")),
            column(4, checkboxInput("gjam_do_sensitivity", "Run gjamSensitivity", TRUE),
                   help_text("Evaluate sensitivity of the joint response to predictors. This is one of GJAM's key ecological outputs.")),
            column(4, checkboxInput("gjam_do_ordination", "Run gjamOrdination", FALSE),
                   help_text("Ordinate the response matrix / fitted response structure for visualization of community gradients."))
          ),
          fluidRow(
            column(4, checkboxInput("gjam_do_conditional", "Run gjamConditionalParameters", FALSE),
                   help_text("Obtain conditional coefficient matrices for conditional relationships among responses.")),
            column(4, checkboxInput("gjam_do_iie", "Run gjamIIE / gjamIIEplot", FALSE),
                   help_text("Evaluate indirect effects and interactions when appropriate.")),
            column(4, checkboxInput("gjam_do_traits", "Run trait analysis / gjamSpec2Trait", FALSE),
                   help_text("Use species-by-trait information to create plot-by-trait matrices and trait-based outputs."))
          ),
          fluidRow(
            column(4, checkboxInput("gjam_missingX", "Allow missing X prediction", TRUE),
                   help_text("GJAM can conditionally predict missing predictor values; this should be explicitly reported.")),
            column(4, checkboxInput("gjam_missingY", "Allow missing Y prediction", TRUE),
                   help_text("GJAM can conditionally predict missing response values; this should be explicitly reported.")),
            column(4, checkboxInput("gjam_inverse_prediction", "Inverse prediction of X", TRUE),
                   help_text("Inverse prediction predicts environmental inputs from the joint response distribution. It is controlled by PREDICTX."))
          )
        ),
        accordion_panel("7. GJAM output settings",
          div(class="note", tags$b("GJAM has its own output archive."), " It is not mixed with Hmsc or jSDM outputs."),
          fluidRow(
            column(4, checkboxInput("gjam_out_model", "Save GJAM model object", TRUE),
                   help_text("Save fitted object as models/gjam_model.rds when production fitting is enabled.")),
            column(4, checkboxInput("gjam_out_chains", "Save MCMC chains", TRUE),
                   help_text("Save chains such as bgibbs, bgibbsUn, fgibbs, fbgibbs, sgibbs and optionally ygibbs.")),
            column(4, checkboxInput("gjam_out_parameters", "Save parameter tables", TRUE),
                   help_text("Save betaMu, betaSe, betaMuUn, betaSeUn, fBetaMu, fBetaSd, corMu, corSe, sigMu, sigSe, fmatrix, fMu, fSe and ematrix when available."))
          ),
          fluidRow(
            column(4, checkboxInput("gjam_out_fit", "Save fit diagnostics", TRUE),
                   help_text("Save fit metrics such as DIC, rmspeAll, rmspeBySpec, xscore and yscore.")),
            column(4, checkboxInput("gjam_out_prediction", "Save predictions", TRUE),
                   help_text("Save prediction outputs, including richness, xpred, ypred and uncertainty summaries when available.")),
            column(4, checkboxInput("gjam_out_missing", "Save missing-data predictions", TRUE),
                   help_text("Save missing predictor/response locations and predicted means/standard errors."))
          ),
          fluidRow(
            column(4, checkboxInput("gjam_out_plots", "Save gjamPlot figures", TRUE),
                   help_text("Save visual summaries from gjamPlot and optional sensitivity/ordination/IIE plots.")),
            column(4, checkboxInput("gjam_out_report", "Generate GJAM HTML report", TRUE),
                   help_text("Create report/GJAM_report.html with inputs, settings, warnings and output explanation.")),
            column(4, checkboxInput("gjam_out_zip", "Create GJAM output ZIP", TRUE),
                   help_text("Package the GJAM output folder as a downloadable ZIP."))
          ),
          hr(),
          checkboxInput("gjam_real_fit", "Run real GJAM package fit", TRUE),
          help_text("Default: TRUE. Calls the installed gjam package, exports an executable reproducible script, and writes real diagnostics if fitting fails."),
          checkboxInput("gjam_out_copy_inputs", "Copy GJAM input files to output folder", TRUE),
          help_text("Retain original uploads for an export-only run. A fitted analysis always saves its required inputs for replay."),
          checkboxInput("gjam_out_config", "Save GJAM used_config.yml", TRUE),
          help_text("Save settings with an export-only run. Fitting always records the configuration it used.")
        ),
        accordion_panel("8. GJAM run and results", value = "results",
          fluidRow(
            column(4,
              actionButton("gjam_run", "Run GJAM workflow", class = "btn-success"),
              br(), br(),
              downloadButton("gjam_download", "Download GJAM ZIP", icon = NULL),
              div(class="warn", "With real fit enabled, this workflow runs the installed gjam package, exports an executable R script, and writes diagnostics if fitting fails.")
            ),
            column(8,
              div(class="cardx", h3("GJAM run summary"), verbatimTextOutput("gjam_run_summary")),
              div(class="cardx", h3("GJAM log"), tags$pre(class="logbox", textOutput("gjam_log", inline = TRUE)))
            )
          ),
          div(class="cardx", h3("GJAM output files"), DTOutput("gjam_files"))
        )
      )
    )
  ),


  tabPanel("spOccupancy", value = "spOccupancy Workflow",
    fluidPage(
      workflow_overview("engine-spoccupancy", "spOccupancy workflow",
        "Single- and multi-species occupancy models that account for imperfect detection.",
        "Use real repeat surveys. Code a surveyed nondetection as 0 and an unobserved visit as NA.",
        c("Occupancy covariates affect presence probability (psi); detection covariates affect detection probability (p) conditional on presence.",
          "Choose a branch that matches the data: spatial, latent-factor, integrated and varying-coefficient models have different inputs.",
          "Do not equate estimated occupancy probability with detection probability or a sampled latent state.",
          "Studio does not expose upstream multi-season functions.")),
      accordion(
        class = "studio-workflow-form", open = "inputs",
        accordion_panel("1. spOccupancy data upload", value = "inputs",
          fluidRow(
            column(6,
              fileInput("spocc_y_file", "y.csv  - detection / nondetection data", accept = ".csv"),
              help_text("Use one row per site and one column per visit for one species. For multiple species use columns such as sp1_rep1, sp1_rep2, sp2_rep1, sp2_rep2. Use 0 for a surveyed nondetection and NA for an unobserved visit."),
              fileInput("spocc_occ_file", "occ.covs.csv  - occurrence covariates", accept = ".csv"),
              help_text("Site-level covariates for the occurrence process psi. Rows should match sites. Example columns: elevation, forest cover, climate, habitat, management."),
              fileInput("spocc_det_file", "det.covs.csv  - detection covariates", accept = ".csv"),
              help_text("Observation-level or site/replicate-level covariates for detection probability p. Examples: observer, date, time, wind, survey method, effort.")
            ),
            column(6,
              fileInput("spocc_coords_file", "coords.csv  - spatial coordinates", accept = ".csv"),
              help_text("Required for spPGOcc, spMsPGOcc, lfMsPGOcc and svcPGOcc. Supply two projected x/y coordinates per site in the same order as y.csv. Do not supply geographic degrees as distances."),
              fileInput("spocc_species_file", "species.csv  - species names / metadata", accept = ".csv"),
              help_text("Optional species list; the adapter uses species names to interpret multi-species visit columns. Extra trait columns are not fitted by this workflow."),
              fileInput("spocc_integrated_file", "integrated_sources.csv  - integrated data-source metadata", accept = ".csv"),
              help_text("Required for integrated occupancy models. Describes multiple detection-nondetection datasets, source IDs, sample units and data-source specific settings.")
            )
          ),
          fluidRow(
            column(4, fileInput("spocc_newdata_file", "newdata.csv  - prediction occurrence covariates", accept = ".csv"),
                   help_text("Optional. Covariates at prediction sites for predict.spOccupancy model methods.")),
            column(4, fileInput("spocc_newcoords_file", "newcoords.csv  - prediction coordinates", accept = ".csv"),
                   help_text("Required with newdata.csv for spatial and latent-factor predictions. Supply the same projected coordinate system as the training sites."))
          )
        ),
        accordion_panel("2. spOccupancy data check", value = "validation",
          actionButton("spocc_check", "Check spOccupancy data", class = "btn-primary"),
          br(), br(),
          fluidRow(
            column(5, div(class="cardx", h3("spOccupancy check messages"), verbatimTextOutput("spocc_check_messages"))),
            column(7, div(class="cardx", h3("spOccupancy data dimensions"), DTOutput("spocc_data_table")))
          ),
          file_preview_block("spocc", "spOccupancy uploaded data preview"),
          div(class="note", tags$b("Occupancy check:"), " spOccupancy needs detection-nondetection data. Spatial models require coords. Integrated models require data-source metadata. Latent-factor models require n.factors >= 1.")
        ),
        accordion_panel("3. spOccupancy model settings", value = "settings",
          fluidRow(
            column(4, selectInput("spocc_model_type", "Occupancy model function",
              choices = c("Single species" = "PGOcc", "Single species, spatial" = "spPGOcc",
                          "Multiple species" = "msPGOcc", "Multiple species, spatial" = "spMsPGOcc",
                          "Multiple species, latent factors" = "lfMsPGOcc",
                          "Single species, integrated surveys" = "intPGOcc",
                          "Single species, spatially varying coefficients" = "svcPGOcc"),
              selected = "PGOcc"),
              help_text("All seven choices model detection separately from occupancy. The displayed function name is saved in used_config.yml. Multi-season and no-detection functions in the R package are not exposed here.")),
            column(4, textInput("spocc_occ_formula", "occ.formula", value = "~ 1"),
              help_text("One-sided occurrence formula for occupancy probability psi. Example: ~ elev + forest + (1 | region). Random intercepts use lme4 syntax in supported functions.")),
            column(4, textInput("spocc_det_formula", "det.formula", value = "~ 1"),
              help_text("One-sided detection formula for detection probability p. Example: ~ day + observer + wind. Imperfect detection is the key distinction of occupancy models."))
          ),
          fluidRow(
            tags$p(class = "small-muted", "The function determines the data layout. Priors and starting values are set by the adapter and recorded in the exported R script.")
          )
        ),
        accordion_panel("4. spOccupancy spatial / latent-factor / SVC settings",
          fluidRow(
            column(3, selectInput("spocc_cov_model", "cov.model", choices = c("exponential", "spherical", "matern", "gaussian"), selected = "exponential"),
              help_text("Spatial covariance function for spatial models. Matern is flexible; exponential is common and stable.")),
            column(3, checkboxInput("spocc_NNGP", "NNGP", TRUE),
              help_text("Use nearest-neighbor Gaussian process for scalable spatial models. Recommended for large spatial datasets.")),
            column(3, synced_numeric_slider("spocc_n_neighbors", "n.neighbors", value = 15, min = 1, max = 50, step = 1,
              help = "Number of neighbors for NNGP. 15 is common; smaller values such as 5 may work in some datasets.")),
            column(3, selectInput("spocc_search_type", "search.type", choices = c("cb", "brute"), selected = "cb"),
              help_text("Neighbor search type for NNGP. 'cb' is the usual fast coordinate-based search."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("spocc_n_factors", "n.factors", value = 3, min = 1, max = 20, step = 1,
              help = "Number of latent factors in lfMsPGOcc. The factors represent residual variation shared by species.")),
            column(3, textInput("spocc_svc_cols", "svc.cols", value = ""),
              help_text("Occurrence design columns whose coefficients vary spatially in svcPGOcc. Example: 1,2 or elev,forest."))
          )
        ),
        accordion_panel("5. spOccupancy MCMC settings",
          fluidRow(
            column(3, synced_numeric_slider("spocc_n_batch", "n.batch", value = 1000, min = 1, max = 50000, step = 100,
              help = "Number of MCMC batches. Total iterations = n.batch * batch.length.")),
            column(3, synced_numeric_slider("spocc_batch_length", "batch.length", value = 25, min = 1, max = 1000, step = 5,
              help = "Number of MCMC iterations per batch. Total iterations = n.batch * batch.length.")),
            column(3, synced_numeric_slider("spocc_n_burn", "n.burn", value = 2500, min = 0, max = 100000, step = 100,
              help = "Burn-in iterations discarded from posterior samples. Must be less than total iterations.")),
            column(3, synced_numeric_slider("spocc_n_thin", "n.thin", value = 1, min = 1, max = 500, step = 1,
              help = "Thinning interval for posterior storage."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("spocc_n_chains", "n.chains", value = 1, min = 1, max = 8, step = 1,
              help = "Number of MCMC chains. More than one chain is recommended for convergence diagnostics.")),
            column(3, synced_numeric_slider("spocc_accept_rate", "accept.rate", value = 0.43, min = 0.01, max = 0.99, step = 0.01,
              help = "Target acceptance rate for adaptive Metropolis updates such as spatial range parameters.")),
            column(3, synced_numeric_slider("spocc_n_report", "n.report", value = 100, min = 1, max = 10000, step = 10,
              help = "Progress reporting interval in batches.")),
            column(3, synced_numeric_slider("spocc_n_omp_threads", "n.omp.threads", value = 1, min = 1, max = 16, step = 1,
              help = "Number of OpenMP threads. Use 1 for safest Windows operation."))
          ),
          fluidRow(
            column(4, checkboxInput("spocc_verbose", "verbose", TRUE),
              help_text("Print progress while running.")),
            column(4, synced_numeric_slider("spocc_seed", "seed", value = 1234, min = 1, max = 999999, step = 1,
              help = "Random seed for the sampler."))
          )
        ),
        accordion_panel("6. spOccupancy assessment, prediction and outputs",
          fluidRow(
            column(4, checkboxInput("spocc_do_ppc", "Posterior predictive check", TRUE),
              help_text("Compare replicated and observed detection histories using ppcOcc for the selected function.")),
            column(4, checkboxInput("spocc_do_waic", "WAIC", TRUE),
              help_text("Compute WAIC where the fitted function supports it. WAIC is not a held-out prediction score.")),
            column(4, checkboxInput("spocc_do_predict", "Predict occupancy at new sites", TRUE),
              help_text("With newdata.csv, export latent occupancy predictions. Training-site detection and occupancy draws are separate outputs. This workflow does not fit held-out folds."))
          ),
          fluidRow(
            column(3, checkboxInput("spocc_save_model", "Save model object", TRUE),
              help_text("Save fitted spOccupancy object as models/spOccupancy_model.rds when production fitting is connected.")),
            column(3, checkboxInput("spocc_save_samples", "Save posterior samples", TRUE),
              help_text("Save coda mcmc samples and parameter arrays to samples/.")),
            column(3, checkboxInput("spocc_out_zip", "Create spOccupancy output ZIP", TRUE),
              help_text("Package the spOccupancy output folder as a downloadable ZIP."))
          ),
          hr(),
          checkboxInput("spocc_real_fit", "Run real spOccupancy package fit", TRUE),
          help_text("Default: TRUE. Runs the installed spOccupancy package and exports a standalone reproducible R script. If disabled, only the scaffold is written."),
          checkboxInput("spocc_out_copy_inputs", "Copy spOccupancy input files to output folder", TRUE),
          help_text("Retain detection histories and covariates for an export-only run. A fitted analysis always saves these inputs for replay."),
          checkboxInput("spocc_out_config", "Save spOccupancy used_config.yml", TRUE),
          help_text("Save settings with an export-only run. Fitting always records the configuration it used."),
          checkboxInput("spocc_out_report", "Generate spOccupancy HTML report", TRUE),
          help_text("Write a browser-readable workflow report with the run status and exported files. Machine-readable diagnostics remain available.")
        ),
        accordion_panel("7. spOccupancy run and results", value = "results",
          fluidRow(
            column(4,
              actionButton("spocc_run", "Run spOccupancy workflow", class = "btn-success"),
              br(), br(),
              downloadButton("spocc_download", "Download spOccupancy ZIP", icon = NULL),
              div(class="warn", "Real fitting uses the installed spOccupancy package. Unsupported advanced functions fail with diagnostics rather than showing a false Completed status.")
            ),
            column(8,
              div(class="cardx", h3("spOccupancy run summary"), verbatimTextOutput("spocc_run_summary")),
              div(class="cardx", h3("spOccupancy log"), tags$pre(class="logbox", textOutput("spocc_log", inline = TRUE)))
            )
          ),
          div(class="cardx", h3("spOccupancy output files"), DTOutput("spocc_files"))
        )
      )
    )
  ),


  tabPanel("sjSDM", value = "s-jSDM Workflow",
    fluidPage(
      workflow_overview("engine-sjsdm", "sjSDM workflow",
        "Regularized environmental, spatial and species association models fitted with PyTorch.",
        "Check the PyTorch environment and record the selected device and regularization settings.",
        c("Penalty and tuning choices affect coefficient magnitudes and which effects shrink toward zero.",
          "Trait inputs support downstream interpretation and are not an Hmsc-style trait hierarchy.",
          "Residual covariance and correlation describe dependence under the fitted model, not direct interactions or the Hmsc Omega estimand.")),
      accordion(
        class = "studio-workflow-form", open = "inputs",
        accordion_panel("1. sjSDM data upload", value = "inputs",
          fluidRow(
            column(6,
              fileInput("sjsdm_Y_file", "Y.csv  - response matrix", accept = ".csv"),
              help_text("Required. Rows = sites/samples. Columns = species, OTUs, ASVs or response variables. Common use: presence/absence eDNA or large community matrices. Binary, count and continuous families are available."),
              fileInput("sjsdm_env_file", "env.csv / XData.csv  - environmental predictors", accept = ".csv"),
              help_text("Recommended. Rows must match Y. Columns are environmental covariates used in the env module. Example: pH, temperature, precipitation, forest cover, treatment."),
              fileInput("sjsdm_spatial_file", "spatial.csv  - coordinates / spatial predictors", accept = ".csv"),
              help_text("Optional. Rows must match Y. Can be raw coordinates, spatial eigenvectors from generateSpatialEV, or other spatial covariates used in the spatial module.")
            ),
            column(6,
              fileInput("sjsdm_traits_file", "traits.csv  - species traits / metadata", accept = ".csv"),
              help_text("Optional. Rows should match Y columns. Used for post-hoc internal structure / assembly-effect interpretation rather than the core sjSDM fitting function."),
              fileInput("sjsdm_newdata_file", "newdata.csv  - prediction environmental covariates", accept = ".csv"),
              help_text("Optional. New environmental predictors for predict.sjSDM. Use with new_spatial.csv when spatial terms are included."),
              fileInput("sjsdm_group_file", "species_groups.csv  - species grouping / colors", accept = ".csv"),
              help_text("Optional. Used for plot.sjSDM grouping. Example columns: species, group.")
            )
          ),
          fluidRow(
            column(4, fileInput("sjsdm_spatial_new_file", "new_spatial.csv  - prediction spatial covariates", accept = ".csv"), help_text("Optional. Spatial predictors/coordinates for prediction sites.")),
            column(4, fileInput("sjsdm_weights_file", "weights.rds / pretrained weights", accept = c(".rds", ".RDS")), help_text("Optional. Pretrained DNN weights for setWeights, or weights exported from a previous model.")),
            column(4, fileInput("sjsdm_cv_file", "folds.csv - validation metadata", accept = ".csv"), help_text("Optional record of your external validation groups. This file is archived; it does not control the folds generated by sjSDM_cv."))
          )
        ),
        accordion_panel("2. sjSDM data check", value = "validation",
          actionButton("sjsdm_check", "Check sjSDM data", class = "btn-primary"), br(), br(),
          fluidRow(
            column(5, div(class="cardx", h3("sjSDM check messages"), verbatimTextOutput("sjsdm_check_messages"))),
            column(7, div(class="cardx", h3("sjSDM data dimensions"), DTOutput("sjsdm_data_table")))
          ),
          file_preview_block("sjsdm", "sjSDM uploaded data preview"),
          div(class="note", tags$b("sjSDM compatibility check:"), " Y and env rows must match. Spatial predictors must match Y rows. Binary/count/continuous families require suitable Y values. GPU/CUDA availability is checked by the production adapter.")
        ),
        accordion_panel("3. sjSDM model settings", value = "settings",
          fluidRow(
            column(4, selectInput("sjsdm_family", "family", choices = c("binomial(probit)" = "binomial_probit", "binomial(logit)" = "binomial_logit", "poisson(log)" = "poisson_log", "nbinom(log, experimental)" = "nbinom", "gaussian(identity)" = "gaussian_identity"), selected = "binomial_probit"), help_text("Response distribution mapped to sjSDM 1.0.7: binomial('probit'/'logit'), poisson('log'), 'nbinom' with log link, or gaussian('identity').")),
            column(4, selectInput("sjsdm_env_model", "env module", choices = c("linear", "DNN", "intercept-only"), selected = "linear"), help_text("Environmental response model. sjSDM always needs an env design. Use intercept-only only when you deliberately want no environmental predictors.")),
            column(4, textInput("sjsdm_env_formula", "env formula", value = "~ ."), help_text("Formula for environmental predictors. Example: ~ pH + moisture + treatment. Use ~ . for all columns; use ~ 0 + . to remove intercept."))
          ),
          fluidRow(
            column(4, selectInput("sjsdm_spatial_model", "spatial module", choices = c("none", "linear", "DNN", "spatial eigenvectors" = "eigenvectors"), selected = "none"), help_text("Spatial module. none = no spatial term; linear = linear spatial predictors; DNN = nonlinear spatial term; spatial eigenvectors = use uploaded or generated eigenvectors.")),
            column(4, textInput("sjsdm_spatial_formula", "spatial formula", value = "~ 0 + ."), help_text("Formula for spatial predictors. Common example: ~ 0 + . for spatial eigenvectors or coordinates without intercept.")),
            column(4, checkboxInput("sjsdm_se", "se = TRUE", FALSE), help_text("Compute post-hoc standard errors / p-values with getSe or built-in se option. Can be slow for large models."))
          ),
          fluidRow(
            column(4, synced_numeric_slider("sjsdm_iter", "iter", value = 100, min = 1, max = 10000, step = 10, help = "Number of training iterations / epochs. Use small values only for testing; increase for real analysis.")),
            column(4, synced_numeric_slider("sjsdm_sampling", "sampling", value = 5000, min = 100, max = 100000, step = 500, help = "Monte Carlo samples for approximating the joint likelihood. Larger values improve stability but increase runtime/memory.")),
            column(4, synced_numeric_slider("sjsdm_seed", "seed", value = 1234, min = 1, max = 999999, step = 1, help = "Random seed for reproducibility."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("sjsdm_step_size", "step_size / batch size", value = 50, min = 1, max = 100000, step = 10, help = "sjSDM step_size. Internally this is passed as the stochastic-gradient batch size. If left too large for small data, reduce it.")),
            column(3, synced_numeric_slider("sjsdm_parallel", "parallel data-loader workers", value = 0, min = 0, max = 64, step = 1, help = "sjSDM parallel argument. 0 is safest on Windows; use more workers only for large datasets.")),
            column(3, selectInput("sjsdm_dtype", "dtype", choices = c("float32", "float64"), selected = "float32"), help_text("Numerical precision. float32 is sjSDM's GPU-friendly default; float64 can be slower and memory-heavy.")),
            column(3, checkboxInput("sjsdm_verbose", "verbose", TRUE), help_text("Print sjSDM fitting progress in the reproducible script."))
          )
        ),
        accordion_panel("4. sjSDM regularization and biotic structure",
          fluidRow(
            column(4, synced_numeric_slider("sjsdm_env_lambda", "env lambda", value = 0, min = 0, max = 1, step = 0.001, help = "Elastic-net regularization strength for environmental coefficients. lambda = 0 means no regularization.")),
            column(4, synced_numeric_slider("sjsdm_env_alpha", "env alpha", value = 0.5, min = 0, max = 1, step = 0.05, help = "Elastic-net mixing for environmental coefficients. alpha=0 pure lasso; alpha=1 pure ridge; 0.5 mixes both.")),
            column(4, synced_numeric_slider("sjsdm_spatial_lambda", "spatial lambda", value = 0, min = 0, max = 1, step = 0.001, help = "Elastic-net regularization strength for spatial module weights/coefficients."))
          ),
          fluidRow(
            column(4, synced_numeric_slider("sjsdm_spatial_alpha", "spatial alpha", value = 0.5, min = 0, max = 1, step = 0.05, help = "Elastic-net lasso/ridge mixing for spatial module.")),
            column(4, synced_numeric_slider("sjsdm_biotic_lambda", "biotic lambda", value = 0.01, min = 0, max = 1, step = 0.001, help = "Regularization strength for species-species covariance/association matrix using bioticStruct. Important for sparse associations.")),
            column(4, synced_numeric_slider("sjsdm_biotic_alpha", "biotic alpha", value = 0.5, min = 0, max = 1, step = 0.05, help = "Elastic-net mixing for biotic associations. alpha=0 lasso; alpha=1 ridge."))
          ),
          fluidRow(
            column(3, numericInput("sjsdm_biotic_df", "biotic df", value = NA, min = 1, step = 1), help_text("Degrees of freedom for covariance parameterization. Leave NA for package default, usually ncol(Y)/2.")),
            column(3, checkboxInput("sjsdm_on_diag", "regularize on_diag", FALSE), help_text("Regularize diagonal covariance entries. Usually FALSE unless you deliberately want diagonal shrinkage.")),
            column(3, checkboxInput("sjsdm_reg_on_Cov", "reg_on_Cov", TRUE), help_text("Regularize the covariance matrix. This is the ordinary setting for bioticStruct.")),
            column(3, checkboxInput("sjsdm_inverse", "regularize inverse covariance", FALSE), help_text("Regularize inverse covariance matrix instead of covariance matrix. Advanced use only."))
          ),
          fluidRow(
            column(4, checkboxInput("sjsdm_tune_regularization", "Tune regularization by CV", FALSE), help_text("Use sjSDM_cv / custom tuning to choose regularization strengths. Recommended when association sparsity matters.")),
            column(4, synced_numeric_slider("sjsdm_tune_steps", "tune random steps", value = 0, min = 0, max = 500, step = 5, help = "Number of random regularization combinations to evaluate. 0 = no tuning.")),
            column(4, synced_numeric_slider("sjsdm_cv_k", "CV folds", value = 0, min = 0, max = 20, step = 1, help = "Number of folds for cross-validation. Use 0 to skip; use >=2 for tuning or prediction validation."))
          )
        ),
        accordion_panel("5. sjSDM DNN and optimizer settings",
          fluidRow(
            column(4, textInput("sjsdm_dnn_hidden", "DNN hidden layers", value = "10,10,10"), help_text("Comma-separated hidden units for DNN modules. Example: 10,10,10 = three hidden layers with 10 neurons each.")),
            column(4, selectInput("sjsdm_activation", "DNN activation", choices = c("selu","relu","leakyrelu","tanh","sigmoid"), selected = "selu"), help_text("Activation function for DNN modules. selu is package default; relu is common; tanh/sigmoid can be useful for bounded nonlinearities.")),
            column(4, synced_numeric_slider("sjsdm_dropout", "dropout", value = 0, min = 0, max = 0.95, step = 0.05, help = "Dropout probability for DNN regularization. 0 = no dropout. Use carefully with ecological sample sizes."))
          ),
          fluidRow(
            column(4, checkboxInput("sjsdm_bias", "DNN bias", TRUE), help_text("Use bias terms in DNN layers.")),
            column(4, selectInput("sjsdm_optimizer", "optimizer", choices = c("Adamax","RMSprop","SGD","AccSGD","AdaBound","madgrad"), selected = "Adamax"), help_text("PyTorch optimizer exposed and exported by sjSDM 1.0.7. Only valid exported optimizers are shown.")),
            column(4, synced_numeric_slider("sjsdm_learning_rate", "learning rate", value = 0.003, min = 1e-6, max = 0.1, step = 0.001, help = "Optimizer learning rate. Too high can diverge; too low can be slow."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("sjsdm_weight_decay", "weight_decay", value = 0.002, min = 0, max = 1, step = 0.001, help = "L2 penalty in optimizer where supported, separate from model-module elastic-net penalties.")),
            column(3, synced_numeric_slider("sjsdm_scheduler", "control$scheduler", value = 0, min = 0, max = 1000, step = 1, help = "sjSDMControl scheduler patience. 0 disables learning-rate scheduling.")),
            column(3, synced_numeric_slider("sjsdm_lr_reduce_factor", "control$lr_reduce_factor", value = 0.99, min = 0.01, max = 0.999, step = 0.01, help = "Learning-rate reduction factor used by sjSDMControl when scheduler is enabled.")),
            column(3, synced_numeric_slider("sjsdm_early_stopping", "control$early_stopping_training", value = 0, min = 0, max = 10000, step = 1, help = "Early stopping patience on training loss. 0 disables early stopping."))
          ),
          fluidRow(
            column(4, selectInput("sjsdm_device", "device", choices = c("cpu","gpu"), selected = "cpu"), help_text("Compute device. sjSDM 1.0.7 accepts 'cpu', 'gpu' or numeric GPU ids. Use gpu only when torch CUDA is available.")),
            column(4, checkboxInput("sjsdm_mixed", "control$mixed half precision", FALSE), help_text("Mixed half-precision training. Only recommended for newer CUDA GPUs; keep off for CPU.")),
            column(4, div(class="note", "Actual sjSDM call: sjSDM(Y, env, biotic, spatial, family, iter, step_size, learning_rate, se, sampling, parallel, control, device, dtype, seed, verbose)."))
          )
        ),
        accordion_panel("6. sjSDM spatial, ANOVA and metacommunity settings",
          fluidRow(
            column(4, checkboxInput("sjsdm_generate_spatial_ev", "generateSpatialEV", FALSE), help_text("Generate spatial eigenvectors from coordinates before fitting. Useful for spatial autocorrelation when raw coordinates are provided.")),
            column(4, synced_numeric_slider("sjsdm_spatial_ev_k", "number of spatial eigenvectors", value = 20, min = 1, max = 500, step = 1, help = "Number of spatial eigenvectors to retain if generated by the adapter.")),
            column(4, synced_numeric_slider("sjsdm_spatial_ev_threshold", "generateSpatialEV threshold", value = 0, min = 0, max = 1, step = 0.01, help = "Threshold argument passed to generateSpatialEV(coords, threshold). 0 keeps the package default behavior."))
          ),
          fluidRow(
            column(4, checkboxInput("sjsdm_include_space_in_anova", "include space fraction in ANOVA", TRUE), help_text("Include spatial module in variation partitioning / ANOVA.")),
            column(8, div(class="note", "For spatial eigenvector filtering, use spatial.csv with at least two numeric coordinate columns, enable generateSpatialEV, then keep spatial formula as ~ 0 + . for the retained eigenvectors."))
          ),
          fluidRow(
            column(4, checkboxInput("sjsdm_do_anova", "Run anova.sjSDM", TRUE), help_text("Compute variation partitioning for environment, space and associations.")),
            column(4, synced_numeric_slider("sjsdm_anova_samples", "ANOVA samples", value = 5000, min = 100, max = 100000, step = 500, help = "Monte Carlo samples for anova.sjSDM. Many samples are recommended for stable results.")),
            column(4, checkboxInput("sjsdm_do_internal", "Run internalStructure", TRUE), help_text("Partition internal metacommunity structure across species and sites after ANOVA."))
          ),
          fluidRow(
            column(4, selectInput("sjsdm_internal_fractions", "internal fractions handling", choices = c("proportional","discard","equal"), selected = "proportional"), help_text("How to handle shared fractions in internalStructure / summary.sjSDManova.")),
            column(4, checkboxInput("sjsdm_do_assembly", "Run plotAssemblyEffects", TRUE), help_text("Plot assembly effects from internalStructure against environmental, spatial, richness or trait predictors.")),
            column(4, textInput("sjsdm_assembly_predictor", "assembly predictor", value = ""), help_text("Optional column name from env/spatial/traits used in plotAssemblyEffects. Empty = use default summaries."))
          )
        ),
        accordion_panel("7. sjSDM prediction, importance and outputs",
          fluidRow(
            column(4, checkboxInput("sjsdm_do_predict", "Run predict.sjSDM", TRUE), help_text("Generate predictions for uploaded newdata or training env data. If spatial module is used, SP/new_spatial should also be supplied.")),
            column(4, checkboxInput("sjsdm_do_rsquared", "Run Rsquared", TRUE), help_text("Compute R-squared / pseudo-R2 metrics for fitted model.")),
            column(4, checkboxInput("sjsdm_do_importance", "Run getImportance / importance", TRUE), help_text("Compute variable importance for environmental/spatial/DNN predictors when supported."))
          ),
          fluidRow(
            column(4, checkboxInput("sjsdm_do_weights", "Save DNN weights", TRUE), help_text("Save DNN weights from getWeights for reproducibility and reuse with setWeights.")),
            column(4, checkboxInput("sjsdm_do_coef", "Save coef", TRUE), help_text("Save environmental and spatial coefficients from coef.sjSDM.")),
            column(4, checkboxInput("sjsdm_do_cov_cor", "Save getCov / getCor", TRUE), help_text("Save species covariance and correlation matrices from getCov and getCor."))
          ),
          fluidRow(
            column(4, checkboxInput("sjsdm_do_residuals", "Save residuals", TRUE), help_text("Save residuals.sjSDM outputs for diagnostic plots.")),
            column(4, checkboxInput("sjsdm_do_plots", "Save plots", TRUE), help_text("Save plot.sjSDM, plot.sjSDManova, plot.sjSDMimportance, plot.sjSDMinternalStructure and assembly-effect plots where requested.")),
            column(4, checkboxInput("sjsdm_out_zip", "Create sjSDM output ZIP", TRUE), help_text("Package the sjSDM output folder as a downloadable ZIP."))
          ),
          hr(),
          checkboxInput("sjsdm_real_fit", "Fit with sjSDM", TRUE),
          help_text("Run the installed sjSDM backend. Clear to export an analysis scaffold without estimates or predictions."),
          checkboxInput("sjsdm_out_copy_inputs", "Copy sjSDM input files to output folder", TRUE),
          help_text("Retain original uploads for an export-only run. A fitted analysis always saves its required inputs for replay."),
          checkboxInput("sjsdm_out_config", "Save sjSDM used_config.yml", TRUE),
          help_text("Save settings with an export-only run. Fitting always records the configuration it used."),
          checkboxInput("sjsdm_out_report", "Generate sjSDM HTML report", TRUE),
          help_text("Write a browser-readable report with the run status and files. Machine-readable diagnostics remain available.")
        ),
        accordion_panel("8. sjSDM run and results", value = "results",
          fluidRow(
            column(4, actionButton("sjsdm_run", "Run sjSDM workflow", class = "btn-success"), br(), br(), downloadButton("sjsdm_download", "Download sjSDM ZIP", icon = NULL), div(class="warn", "If sjSDM is installed and real fit is enabled, the GUI runs the exported reproducible R script and saves the fitted model, coefficients, associations and predictions. If the package is unavailable, a diagnostic scaffold is still exported.")),
            column(8, div(class="cardx", h3("sjSDM run summary"), verbatimTextOutput("sjsdm_run_summary")), div(class="cardx", h3("sjSDM log"), tags$pre(class="logbox", textOutput("sjsdm_log", inline = TRUE))))
          ),
          div(class="cardx", h3("sjSDM output files"), DTOutput("sjsdm_files"))
        )
      )
    )
  ),

  tabPanel("boral", value = "boral Workflow",
    fluidPage(
      workflow_overview("engine-boral", "boral workflow",
        "Bayesian ordination and multivariate regression with latent variables and optional trait effects.",
        "Check that JAGS is available. For binomial responses, supply trial sizes that match the observations.",
        c("The number and type of latent variables determine the meaning of the ordination and residual dependence.",
          "Fourth-corner terms require matching the trait rows to the species and matching XData to Y by site.",
          "Grouping IDs define which rows share a random effect. Report the grouping alongside association results.",
          "A missing JAGS runtime is a dependency failure, not a fitted analysis.")),
      accordion(
        class = "studio-workflow-form", open = "inputs",
        accordion_panel("1. boral data upload", value = "inputs",
          fluidRow(
            column(6,
              fileInput("boral_Y_file", "Y.csv  - response matrix", accept = ".csv"),
              help_text("Required. Rows = sites/samples. Columns = species or response variables. Examples: spider counts, plant abundance, presence/absence, biomass, ordinal cover scores, traits-as-responses."),
              fileInput("boral_X_file", "XData.csv  - covariate matrix X", accept = ".csv"),
              help_text("Optional. Rows must match Y. Environmental predictors, treatments or site covariates. If omitted and num.lv > 0, boral fits pure latent-variable ordination."),
              fileInput("boral_traits_file", "traits.csv  - species traits", accept = ".csv"),
              help_text("Optional. Rows should match Y columns/species. Used for fourth-corner style models where traits explain species-specific environmental responses.")
            ),
            column(6,
              fileInput("boral_rowids_file", "row.ids.csv  - row effect IDs", accept = ".csv"),
              help_text("Optional. IDs for row effects. Used when sites/samples have grouping, repeated measures or structured row effects."),
              fileInput("boral_ranefids_file", "ranef.ids.csv  - response-specific random intercept IDs", accept = ".csv"),
              help_text("Optional. Matrix/data frame with rows matching Y; columns define random-intercept grouping factors for response-specific random effects."),
              fileInput("boral_distmat_file", "distmat.csv  - distance matrix for structured latent variables", accept = ".csv"),
              help_text("Optional. Square n_sites x n_sites distance matrix. Required for lv.type = exponential, squared.exponential, power.exponential or spherical.")
            )
          ),
          fluidRow(
            column(4, fileInput("boral_offset_file", "offset.csv  - offset matrix", accept = ".csv"),
                   help_text("Optional. Offset matrix with the same dimensions as Y. Use for sampling effort, exposure or known offsets.")),
            column(4, fileInput("boral_newdata_file", "newdata.csv  - prediction covariates", accept = ".csv"),
                   help_text("Optional. New X matrix/data frame for predict.boral. Columns should match X/formula.X.")),
            column(4, fileInput("boral_trials_file", "trial.size.csv  - optional binomial trials", accept = ".csv"),
                   help_text("Optional. Binomial trial sizes. Use scalar trial.size in settings for ordinary Bernoulli data."))
          )
        ),
        accordion_panel("2. boral data check", value = "validation",
          actionButton("boral_check", "Check boral data", class = "btn-primary"),
          br(), br(),
          fluidRow(
            column(5, div(class="cardx", h3("boral check messages"), verbatimTextOutput("boral_check_messages"))),
            column(7, div(class="cardx", h3("boral data dimensions"), DTOutput("boral_data_table")))
          ),
          file_preview_block("boral", "boral uploaded data preview"),
          div(class="note", tags$b("boral compatibility check:"), " Y is required. X rows, row.ids, ranef.ids and offset dimensions must match Y. Spatial latent-variable structures require a distance matrix. Traits should have one row per response/species.")
        ),
        accordion_panel("3. boral model settings", value = "settings",
          fluidRow(
            column(4, selectInput("boral_model_mode", "model mode", choices = c("pure latent-variable ordination", "independent response GLMs", "correlated response GLMs", "trait/fourth-corner model", "random-effects model", "SSVS variable-selection model"), selected = "correlated response GLMs"),
                   help_text("High-level workflow label. It helps users understand whether the model uses only latent variables, only covariates, both covariates and latent variables, traits, random effects or SSVS.")),
            column(4, selectInput("boral_family", "family", choices = c("binomial","poisson","negative.binomial","normal","tweedie","exponential","gamma","lnormal","beta","ordinal","ztpoisson","ztnegative.binomial"), selected = "negative.binomial"),
                   help_text("Default response distribution. Can be overridden per column with family vector. negative.binomial is often useful for overdispersed count data.")),
            column(4, textAreaInput("boral_family_text", "family vector override", rows = 2, value = ""),
                   help_text("Optional comma-separated family values, one per Y column. Example: poisson,poisson,binomial,normal. Leave empty to use the selected family for all responses."))
          ),
          fluidRow(
            column(4, synced_numeric_slider("boral_num_lv", "num.lv", value = 2, min = 0, max = 20, step = 1,
                   help = "Number of latent variables. 2 is common for ordination plots; 0 fits independent GLMs when X is present.")),
            column(4, selectInput("boral_lv_type", "lv.control$type", choices = c("independent","exponential","squared.exponential","powered.exponential","spherical"), selected = "independent"),
                   help_text("Latent-variable correlation structure across sites. independent is fastest. Spatial structures require distmat.csv and are slower.")),
            column(4, textInput("boral_model_name", "model.name", value = "jagsboralmodel.txt"),
                   help_text("Name/path for generated JAGS model file. The output ZIP stores this under jags/ when production fitting is connected."))
          ),
          fluidRow(
            column(4, textInput("boral_formula_X", "formula.X", value = "~ ."),
                   help_text("Formula used to construct X from XData when XData is a data frame. Use ~ . for all predictors; leave as ~ . for simple use.")),
            column(4, textInput("boral_X_ind", "X.ind", value = ""),
                   help_text("Optional subset/indices of X columns included in the model. Leave empty to use all columns created by formula.X.")),
            column(4, synced_numeric_slider("boral_trial_size", "trial.size", value = 1, min = 1, max = 1000, step = 1,
                   help = "Number of binomial trials. Use 1 for Bernoulli presence/absence; larger values for binomial counts."))
          ),
          fluidRow(
            column(4, selectInput("boral_row_eff", "row.eff", choices = c("none","fixed","random"), selected = "none"),
                   help_text("Optional row effect. fixed adjusts for site total abundance; random treats row effects as normally distributed. Avoid unnecessary row effects with rich X.")),
            column(4, checkboxInput("boral_use_offset", "Use offset", FALSE),
                   help_text("Use uploaded offset.csv. Must have same dimensions as Y.")),
            column(4, checkboxInput("boral_do_fit", "do.fit", TRUE),
                   help_text("If TRUE, run JAGS MCMC in production. If FALSE, build model/script without fitting."))
          )
        ),
        accordion_panel("4. boral traits, random effects and variable selection",
          fluidRow(
            column(4, checkboxInput("boral_use_traits", "Use traits", FALSE),
                   help_text("Enable trait-mediated responses. Traits explain species-specific intercepts and environmental coefficients, similar to fourth-corner models.")),
            column(4, textInput("boral_which_traits", "which.traits", value = ""),
                   help_text("Advanced. List-like specification of which trait columns explain each coefficient. Empty = adapter default all traits for all coefficients. Example concept: intercept and each X coefficient can use selected traits.")),
            column(4, checkboxInput("boral_traits_intercept_warning", "traits matrix has no intercept", TRUE),
                   help_text("boral adds the trait intercept automatically; uploaded traits.csv should not include an intercept column."))
          ),
          fluidRow(
            column(4, checkboxInput("boral_use_ranef", "Use response-specific random intercepts", FALSE),
                   help_text("Use ranef.ids.csv to include response-specific random intercepts for sampling design or grouping."),
            ),
            column(4, checkboxInput("boral_use_ssvs", "Use SSVS variable selection", FALSE),
                   help_text("Enable stochastic search variable selection on X coefficients and/or trait coefficients.")),
            column(4, textInput("boral_ssvs_index", "prior.control$ssvs.index", value = ""),
                   help_text("SSVS index for X covariates. -1=no SSVS; 0=response-specific selection; positive integers=grouped selection. Example: -1,-1,1,1,1."))
          ),
          fluidRow(
            column(4, textInput("boral_ssvs_traitsindex", "prior.control$ssvs.traitsindex", value = ""),
                   help_text("SSVS index for trait coefficients. Must match which.traits structure in production adapter.")),
            column(4, synced_numeric_slider("boral_ssvs_g", "prior.control$ssvs.g", value = 1e-6, min = 0, max = 0.1, step = 1e-6,
                   help = "Spike-and-slab variance ratio g for SSVS. Very small values create a tight spike near zero.")),
            column(4, checkboxInput("boral_save_model", "save.model", FALSE),
                   help_text("Save raw JAGS model. This can be memory-consuming because it can retain large MCMC objects."))
          )
        ),
        accordion_panel("5. boral MCMC and prior settings",
          fluidRow(
            column(3, synced_numeric_slider("boral_n_burnin", "mcmc.control$n.burnin", value = 10000, min = 0, max = 100000, step = 1000,
                   help = "MCMC burn-in iterations. Package default is large; use smaller values only for testing.")),
            column(3, synced_numeric_slider("boral_n_iteration", "mcmc.control$n.iteration", value = 40000, min = 100, max = 200000, step = 1000,
                   help = "Total MCMC iterations. Increase for publication analyses and check convergence.")),
            column(3, synced_numeric_slider("boral_n_thin", "mcmc.control$n.thin", value = 30, min = 1, max = 1000, step = 1,
                   help = "MCMC thinning interval. Larger thinning reduces storage but does not fix poor mixing.")),
            column(3, synced_numeric_slider("boral_seed", "mcmc.control$seed", value = 1234, min = 1, max = 999999, step = 1,
                   help = "Random seed for reproducible JAGS sampling."))
          ),
          fluidRow(
            column(4, textInput("boral_prior_type", "prior.control$type", value = "normal,normal,normal,uniform"),
                   help_text("Prior family vector for parameters. Default concept: normal priors for coefficients and uniform for non-negative parameters.")),
            column(4, textInput("boral_hypparams", "prior.control$hypparams", value = "10,10,10,30"),
                   help_text("Hyperparameters controlling prior variances/ranges. Default follows boral examples/manual.")),
            column(4, checkboxInput("boral_calc_ics", "calc.ics", FALSE),
                   help_text("Calculate information criteria / log-likelihood outputs where supported. Can increase runtime."))
          ),
          fluidRow(
            column(4, checkboxInput("boral_save_mcmc_samples", "Save MCMC samples", TRUE),
                   help_text("Save get.mcmcsamples outputs to mcmc/ for post-hoc analysis.")),
            column(4, checkboxInput("boral_save_hpd", "Save HPD intervals", TRUE),
                   help_text("Save get.hpdintervals outputs to tables/ for posterior uncertainty.")),
            column(4, checkboxInput("boral_save_dic", "Save DIC / log-likelihood measures", TRUE),
                   help_text("Save get.dic, calc.condlogLik, calc.marglogLik and related model comparison measures where available."))
          )
        ),
        accordion_panel("6. boral diagnostics, ordination and prediction settings",
          fluidRow(
            column(4, checkboxInput("boral_do_summary", "summary.boral", TRUE),
                   help_text("Save model summary including posterior medians, HPD intervals and main parameter estimates.")),
            column(4, checkboxInput("boral_do_residual_plot", "plot.boral residual diagnostics", TRUE),
                   help_text("Create residual diagnostic plots, including Dunn-Smyth residuals when supported.")),
            column(4, checkboxInput("boral_do_lvsplot", "lvsplot ordination", TRUE),
                   help_text("Create latent-variable ordination / biplot. Most useful when num.lv is 1 or 2."))
          ),
          fluidRow(
            column(4, textInput("boral_ind_spp", "lvsplot ind.spp", value = "default"),
                   help_text("Which species to show on ordination biplot. default = package choice; or comma-separated species indices/names.")),
            column(4, checkboxInput("boral_do_ranefsplot", "ranefsplot", FALSE),
                   help_text("Plot response-specific random intercepts if ranef.ids are used.")),
            column(4, checkboxInput("boral_do_coefsplot", "coefsplot", TRUE),
                   help_text("Plot covariate coefficient estimates and HPD intervals."))
          ),
          fluidRow(
            column(4, checkboxInput("boral_do_env_cor", "get.enviro.cor", TRUE),
                   help_text("Compute species correlations induced by shared environmental responses.")),
            column(4, checkboxInput("boral_do_resid_cor", "get.residual.cor", TRUE),
                   help_text("Compute residual correlations accounted for by latent variables.")),
            column(4, checkboxInput("boral_do_varpart", "calc.varpart", FALSE),
                   help_text("Compute variation partitioning when model structure supports it."))
          ),
          fluidRow(
            column(4, checkboxInput("boral_do_predict", "predict.boral", TRUE),
                   help_text("Predict responses at newdata or fitted sites. Save predictions and uncertainty summaries.")),
            column(4, checkboxInput("boral_do_fitted", "fitted.boral", TRUE),
                   help_text("Save fitted values for fitted sites.")),
            column(4, checkboxInput("boral_do_tidy", "tidyboral", TRUE),
                   help_text("Create tidy tables of posterior summaries for downstream use."))
          )
        ),
        accordion_panel("7. boral output settings",
          div(class="note", tags$b("boral output is separate."), " It stores JAGS scripts, MCMC samples, ordination plots, residual diagnostics, correlations and prediction tables in its own ZIP."),
          fluidRow(
            column(4, checkboxInput("boral_out_model", "Save boral model object", TRUE),
                   help_text("Save fitted boral object as models/boral_model.rds when production fitting is enabled.")),
            column(4, checkboxInput("boral_out_jags", "Save JAGS model file", TRUE),
                   help_text("Save generated JAGS model script in jags/ for transparency and reproducibility.")),
            column(4, checkboxInput("boral_out_tables", "Save CSV summary tables", TRUE),
                   help_text("Save posterior summaries, HPD intervals, coefficients, correlations, residuals, fit measures and predictions as CSV."))
          ),
          fluidRow(
            column(4, checkboxInput("boral_out_plots", "Save PDF/PNG plots", TRUE),
                   help_text("Save ordination, residual, coefficient and random-effect plots.")),
            column(4, checkboxInput("boral_out_report", "Generate boral HTML report", TRUE),
                   help_text("Create report/boral_report.html summarizing data, settings and output files.")),
            column(4, checkboxInput("boral_out_zip", "Create boral output ZIP", TRUE),
                   help_text("Package the boral output folder as a downloadable ZIP."))
          ),
          help_text("boral retains the input files and configuration required by its exported R script, including for model-only runs.")
        ),
        accordion_panel("8. boral run and results", value = "results",
          fluidRow(
            column(4,
              actionButton("boral_run", "Run boral workflow", class = "btn-success"),
              br(), br(),
              downloadButton("boral_download", "Download boral ZIP", icon = NULL),
              div(class="warn", "Fitting requires boral and JAGS. With fitting disabled, Studio exports reproducible files without fitted estimates.")
            ),
            column(8,
              div(class="cardx", h3("boral run summary"), verbatimTextOutput("boral_run_summary")),
              div(class="cardx", h3("boral log"), tags$pre(class="logbox", textOutput("boral_log", inline = TRUE)))
            )
          ),
          div(class="cardx", h3("boral output files"), DTOutput("boral_files"))
        )
      )
    )
  )

  ),

  tabPanel("Universal Benchmark",
    fluidPage(
      section_header("Universal Benchmark", "Simulation recovery and matched-data comparison across the seven workflows.", "benchmark"),
      div(class = "studio-compact studio-benchmark-source",
        radioButtons("univ_data_source", "Data source",
          choices = c("Observed data" = "empirical", "Simulate data" = "synthetic",
                      "Demonstration case" = "builtin", "Import benchmark ZIP" = "upload"),
          selected = "empirical", inline = TRUE)),
      accordion(
        open = c("2. Data", "4. Run benchmark"),
        accordion_panel("1. Run settings",
          conditionalPanel("input.univ_data_source != 'empirical'",
          tags$details(class = "studio-disclosure", tags$summary("Advanced benchmark parameters"),
          fluidRow(
            column(3, synced_numeric_slider("univ_n_sites", "observations / sites", value = 120, min = 30, max = 1000, step = 10,
                   help = "Reviewer-scale default for statistical recovery. Use 30-50 for quick software checks, 120+ for effect and prediction recovery, and larger sizes for scaling experiments.")),
            column(3, synced_numeric_slider("univ_n_species", "responses / species", value = 12, min = 5, max = 100, step = 1,
                   help = "Reviewer-scale default community size. Species names are shared across Y, traits, phylogeny and all standard tables.")),
            column(3, synced_numeric_slider("univ_n_visits", "detection visits", value = 3, min = 2, max = 5, step = 1,
                   help = "Replicated visits for the spOccupancy detection-nondetection input.")),
            column(3, synced_numeric_slider("univ_seed", "benchmark seed", value = 20260601, min = 1, max = 99999999, step = 1,
                   help = "Seed for the latent ecological truth and derived model-specific inputs."))
          ),
          checkboxInput("univ_allow_large", "Allow manuscript-scale benchmark sizes", TRUE),
          div(class="note", "When this is off, the command-line runner clamps generated synthetic benchmarks to 30-40 sites and 5-8 species for fast checks. Keep it on for Figure 4 statistical recovery runs."),
          h4("Data-generating structure"),
          selectInput("univ_generating_structure", "Generating structure",
            choices = c("Correlated latent probit (Hmsc / Hmsc-HPC / jSDM / sjSDM link-aligned)" = "probit_latent",
                        "Logistic occupancy with replicated detection (spOccupancy process-aligned)" = "logit_occupancy_detection"),
            selected = "probit_latent"),
          div(class="note", "These are documented Studio simulators, not exact simulations from fitted R or Python package models. GJAM and boral remain cross-structure fits in this benchmark. Run both generating structures and report them separately; neither is an all-engine fairness ranking."),
          h4("Signal settings"),
          fluidRow(
            column(3, selectInput("univ_truth_profile", "truth profile",
              choices = c("recovery", "quick", "association_stress"), selected = "recovery"),
              help_text("recovery is the default manuscript benchmark: low predictor collinearity, stronger non-zero effects and residual association truth that is actually used to generate Y.")),
            column(3, synced_numeric_slider("univ_effect_scale", "effect signal scale", value = 1.60, min = 0.10, max = 5, step = 0.05,
                   help = "Multiplier for true species-by-predictor effects. Increase when validating coefficient recovery; avoid unrealistically large values for ecological examples.")),
            column(3, synced_numeric_slider("univ_association_strength", "association signal scale", value = 1.65, min = 0, max = 5, step = 0.05,
                   help = "Strength of the low-rank residual species association used to generate the response matrix. This is separate from phylogeny/taxonomy.")),
            column(3, synced_numeric_slider("univ_residual_sd", "residual latent SD", value = 1.05, min = 0.01, max = 5, step = 0.05,
                   help = "Residual latent variation added after environmental effects. Larger values make association structure more visible but can reduce pure environmental recovery."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("univ_observation_noise_sd", "normal observation noise SD", value = 0.15, min = 0.001, max = 5, step = 0.01,
                   help = "Observation noise for the continuous Y_normal benchmark response.")),
            column(9, div(class="warn", tags$b("Statistical recovery note: "),
                          "A and C panels are meaningful only when the truth is identifiable. The benchmark now writes predictor_correlation.csv and residual_association_covariance.csv under truth/ so weak recovery can be diagnosed rather than hidden."))
          ),
          checkboxGroupInput("univ_engines", "Engines to run sequentially",
            choices = c("Hmsc" = "Hmsc", "Hmsc-HPC" = "Hmsc-HPC", "jSDM" = "jSDM", "GJAM" = "GJAM", "spOccupancy" = "spOccupancy", "sjSDM" = "sjSDM", "boral" = "boral"),
            selected = c("Hmsc", "Hmsc-HPC", "jSDM", "GJAM", "spOccupancy", "sjSDM", "boral")),
          h4("Hmsc / Hmsc-HPC sampler length for benchmark"),
          div(class="note", "Default values are reviewer-scale software validation settings, not final publication convergence settings. Increase or replicate them before making strong statistical claims."),
          fluidRow(
            column(3, synced_numeric_slider("univ_hmsc_samples", "Hmsc samples", value = 500, min = 1, max = 20000, step = 10,
                   help = "Posterior samples per Hmsc chain for the Universal Benchmark.")),
            column(3, synced_numeric_slider("univ_hmsc_transient", "Hmsc transient", value = 500, min = 0, max = 50000, step = 10,
                   help = "Hmsc burn-in/transient iterations for the Universal Benchmark.")),
            column(3, synced_numeric_slider("univ_hmsc_thin", "Hmsc thin", value = 1, min = 1, max = 1000, step = 1,
                   help = "Hmsc thinning interval for the Universal Benchmark.")),
            column(3, synced_numeric_slider("univ_hmsc_nChains", "Hmsc chains", value = 2, min = 1, max = 8, step = 1,
                   help = "Hmsc chains for the Universal Benchmark. Keep nParallel = 1 for Windows stability."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("univ_hmschpc_samples", "Hmsc-HPC samples", value = 500, min = 1, max = 20000, step = 10,
                   help = "CPU pyhmsc/Hmsc-HPC retained samples for the Universal Benchmark.")),
            column(3, synced_numeric_slider("univ_hmschpc_transient", "Hmsc-HPC transient", value = 500, min = 0, max = 50000, step = 10,
                   help = "CPU pyhmsc/Hmsc-HPC burn-in/transient iterations for the Universal Benchmark.")),
            column(3, synced_numeric_slider("univ_hmschpc_thin", "Hmsc-HPC thin", value = 1, min = 1, max = 1000, step = 1,
                   help = "CPU pyhmsc/Hmsc-HPC thinning interval for the Universal Benchmark.")),
            column(3, synced_numeric_slider("univ_hmschpc_chains", "Hmsc-HPC chains", value = 1, min = 1, max = 8, step = 1,
                   help = "Compiled Hmsc-HPC chain count. One chain is fastest for CPU benchmark checks; increase only when diagnostics are needed."))
          ),
          h4("Other engine sampler length for benchmark"),
          fluidRow(
            column(3, synced_numeric_slider("univ_jsdm_burnin", "jSDM burnin", value = 200, min = 0, max = 50000, step = 50,
                   help = "Burn-in iterations passed to the jSDM benchmark adapter.")),
            column(3, synced_numeric_slider("univ_jsdm_mcmc", "jSDM mcmc", value = 300, min = 1, max = 100000, step = 50,
                   help = "Saved MCMC iterations passed to the jSDM benchmark adapter.")),
            column(3, synced_numeric_slider("univ_jsdm_thin", "jSDM thin", value = 5, min = 1, max = 1000, step = 1,
                   help = "Thinning interval for jSDM benchmark runs.")),
            column(3, synced_numeric_slider("univ_gjam_ng", "GJAM ng", value = 700, min = 100, max = 100000, step = 100,
                   help = "GJAM MCMC generation count for the benchmark. Increase for manuscript diagnostics."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("univ_gjam_burnin", "GJAM burnin", value = 150, min = 0, max = 50000, step = 50,
                   help = "GJAM burn-in for the benchmark.")),
            column(3, synced_numeric_slider("univ_spocc_n_batch", "spOccupancy n.batch", value = 30, min = 1, max = 10000, step = 5,
                   help = "spOccupancy batch count for the benchmark.")),
            column(3, synced_numeric_slider("univ_spocc_batch_length", "spOccupancy batch.length", value = 10, min = 1, max = 10000, step = 5,
                   help = "spOccupancy samples per batch.")),
            column(3, synced_numeric_slider("univ_spocc_n_burn", "spOccupancy n.burn", value = 100, min = 0, max = 100000, step = 50,
                   help = "spOccupancy burn-in samples."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("univ_spocc_n_thin", "spOccupancy n.thin", value = 5, min = 1, max = 1000, step = 1,
                   help = "spOccupancy thinning interval.")),
            column(3, synced_numeric_slider("univ_spocc_n_chains", "spOccupancy n.chains", value = 1, min = 1, max = 8, step = 1,
                   help = "spOccupancy chain count for benchmark validation.")),
            column(3, synced_numeric_slider("univ_sjsdm_iter", "sjSDM iter", value = 100, min = 1, max = 10000, step = 5,
                   help = "sjSDM training epochs/iterations for the benchmark. Universal Benchmark uses a longer default than a smoke test because coefficient recovery is sensitive to optimization length.")),
            column(3, synced_numeric_slider("univ_sjsdm_sampling", "sjSDM sampling", value = 1000, min = 1, max = 100000, step = 50,
                   help = "sjSDM sampling count for uncertainty/biotic structure summaries."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("univ_sjsdm_step_size", "sjSDM step size", value = 32, min = 1, max = 10000, step = 1,
                   help = "sjSDM optimizer step size / batch setting used by the adapter.")),
            column(3, synced_numeric_slider("univ_sjsdm_learning_rate", "sjSDM learning rate", value = 0.010, min = 1e-6, max = 0.1, step = 0.001,
                   help = "sjSDM optimizer learning rate used only by the Universal Benchmark adapter. The default is tuned for recovery tests on standardized predictors.")),
            column(3, synced_numeric_slider("univ_sjsdm_weight_decay", "sjSDM weight_decay", value = 0, min = 0, max = 1, step = 0.001,
                   help = "Optimizer L2 penalty for sjSDM in the Universal Benchmark. 0 avoids extra shrinkage when the benchmark is judging effect recovery.")),
            column(3, synced_numeric_slider("univ_sjsdm_biotic_lambda", "sjSDM biotic lambda", value = 0.001, min = 0, max = 1, step = 0.001,
                   help = "Regularization strength for the sjSDM species covariance structure. A small value stabilizes association estimates without erasing the simulated signal."))
          ),
          fluidRow(
            column(3, synced_numeric_slider("univ_boral_n_burnin", "boral n.burnin", value = 250, min = 0, max = 100000, step = 50,
                   help = "boral/JAGS burn-in for benchmark validation.")),
            column(3, synced_numeric_slider("univ_boral_n_iteration", "boral n.iteration", value = 750, min = 1, max = 200000, step = 50,
                   help = "boral/JAGS total MCMC iterations for benchmark validation.")),
            column(3, synced_numeric_slider("univ_boral_n_thin", "boral n.thin", value = 10, min = 1, max = 1000, step = 1,
                   help = "boral/JAGS thinning interval."))
          ),
          div(class="note", "Simulation settings are saved with the run. Short chains or optimisation runs are software checks and do not establish inferential adequacy.")
          )),
          conditionalPanel("input.univ_data_source == 'empirical'",
            selectInput("univ_empirical_profile", "Run preset",
              choices = c("Standard" = "standard", "Quick workflow check" = "quick")),
            help_text("All seven engines are attempted in sequence. Quick uses short fits for a smoke test; Standard uses longer fixed per-engine settings. Inspect convergence and uncertainty before inference."))
        ),
        accordion_panel("2. Data",
          conditionalPanel("input.univ_data_source == 'empirical'",
            p("Observed 0/1 incidence, shared environmental predictors and a nonspatial model baseline. Use the individual workflows for traits, phylogeny, spatial structure, counts or continuous responses."),
            fileInput("univ_empirical_y", "Y.csv / Observed species records", accept = ".csv"),
            help_text("Required. First column site_id, followed by species columns with 0/1 values. Example: site_01, 1, 0. Every species needs both 0 and 1 among training sites."),
            fileInput("univ_empirical_x", "XData.csv / Site predictors", accept = ".csv"),
            help_text("Required. First column site_id, then predictors such as TMG and habitat. IDs must match Y.csv. Numeric predictors are scaled from training data; text predictors are treated as factors."),
            tags$details(class = "studio-disclosure", tags$summary("Repeated surveys and new-site prediction"),
              fileInput("univ_empirical_visits", "visits.csv / Repeated detections", accept = ".csv"),
              help_text("Optional; required for spOccupancy. Columns: site_id, visit_id, and the same species as Y.csv. One row per actual visit, with 0/1 survey outcomes; at least two visits per site. Absence of a survey is not a zero. A 1 in Y means detected on at least one visit."),
              fileInput("univ_empirical_new_x", "newdata.csv / New-site predictors", accept = ".csv"),
              help_text("Optional. New, disjoint site_id values and the same predictors as XData.csv. Supply only levels seen in training."),
              fileInput("univ_empirical_test_y", "test_y.csv / Independent test observations", accept = ".csv"),
              help_text("Optional. Independent 0/1 records for the newdata sites and training species. Without this file, new-site predictions are projections and receive no held-out scores.")
            ),
            actionButton("univ_empirical_check", "Check observed data")),
          conditionalPanel("input.univ_data_source == 'synthetic'",
            actionButton("univ_generate", "Generate preview data", class = "btn-primary"),
            help_text("Preview the simulated data from the settings in section 1. A run saves its own inputs.")),
          conditionalPanel("input.univ_data_source == 'builtin'",
            uiOutput("univ_builtin_case_ui"),
            actionButton("univ_use_builtin", "Load demonstration case", class = "btn-primary"),
            help_text("Demonstration cases exercise the software contract. They are not independent empirical analyses.")),
          conditionalPanel("input.univ_data_source == 'upload'",
            fileInput("univ_data_zip", "Complete simulation benchmark ZIP", accept = ".zip"),
            help_text("For a previously generated Studio simulation contract. Use Observed data for real CSV files."),
            actionButton("univ_import_zip", "Import benchmark ZIP", class = "btn-primary"),
            tags$details(class = "studio-disclosure", tags$summary("Required ZIP contents"), DTOutput("univ_required_files"))),
          verbatimTextOutput("univ_data_summary"),
          tags$details(class = "studio-disclosure", tags$summary("Input file inventory"), DTOutput("univ_input_files"))
        ),
        accordion_panel("3. Engine mapping and dependency preflight",
          fluidRow(
            column(4,
              actionButton("univ_preflight", "Run dependency preflight", class = "btn-primary"),
              div(class="note", "Checks R packages, Python path and JAGS-related packages before running engines.")
            ),
            column(8,
              div(class="cardx", h3("Dependency preflight"), DTOutput("univ_preflight_table"))
            )
          )
        ),
        accordion_panel("4. Run benchmark",
          fluidRow(
            column(4,
              actionButton("univ_run_all", "Run benchmark", class = "btn-success"),
              conditionalPanel("input.univ_data_source == 'synthetic'",
                actionButton("univ_run_matrix", "Run both simulation structures", class = "btn-outline-secondary")),
              conditionalPanel("input.univ_data_source != 'empirical'",
                actionButton("univ_compare", "Rebuild simulation comparisons", class = "btn-outline-secondary"))
            ),
            column(8,
              div(class="cardx", h3("Universal benchmark run summary"), verbatimTextOutput("univ_run_summary")),
              conditionalPanel("input.univ_data_source == 'synthetic'",
                tags$details(class = "studio-disclosure", tags$summary("Simulation matrix"), verbatimTextOutput("univ_matrix_summary"))),
              tags$details(class = "studio-disclosure", tags$summary("Run log"),
                tags$pre(class="logbox", textOutput("univ_log", inline = TRUE)))
            )
          ),
          div(class="cardx", h3("Engine status matrix"), DTOutput("univ_status_table"))
        ),
        accordion_panel("5. Unified standard results",
          div(class="cardx", h3("Predictor -> species effects"), DTOutput("univ_effects_table")),
          div(class="cardx", h3("Site -> species predictions"), DTOutput("univ_predictions_table")),
          div(class="cardx", h3("Species association summaries"), DTOutput("univ_associations_table")),
          conditionalPanel("input.univ_data_source == 'empirical'",
            div(class="cardx", h3("Detection-process effects"), DTOutput("univ_detection_table")))
        ),
        accordion_panel("6. Benchmark comparison",
          div(class="cardx", h3("Comparison summary"), DTOutput("univ_compare_table")),
          div(class="cardx", h3("Comparable / non-comparable rules"), DTOutput("univ_comparable_table")),
          div(class="note", "Compare predictions only when the target and validation observations match. Association parameters remain model-specific.")
        ),
        accordion_panel("7. Download",
          fluidRow(
            column(4,
              downloadButton("univ_download", "Download all-model benchmark ZIP", icon = NULL),
              div(class="note", "Includes original inputs, per-engine outputs, model diagnostics, configs and reproducible scripts. Truth tables are included only for simulations.")
            ),
            column(8,
              div(class="cardx", h3("Master output files"), DTOutput("univ_files"))
            )
          )
        )
      )
    )
  ),

  tabPanel("Compare Models",
    tags$main(class = "studio-guide",
      section_header("Compare Models", "Read-only audit of selected run folders"),
      accordion(
        open = c("Run folders", "Run evidence"),
        accordion_panel("Run folders",
          p(class = "small-muted", "Leave a path blank to use the latest run for that engine in this session. Original files are not changed."),
          textInput("compare_hmsc_dir", "Hmsc run folder"),
          help_text("Folder containing used_config.yml and the Hmsc diagnostics and standard summaries."),
          textInput("compare_hmschpc_dir", "Hmsc-HPC run folder"),
          help_text("Select a Hmsc-HPC run, not its parent directory or the ZIP file."),
          textInput("compare_jsdm_dir", "jSDM run folder"),
          help_text("Select an extracted or local jSDM run that contains its status files."),
          textInput("compare_gjam_dir", "GJAM run folder"),
          help_text("Select the GJAM run whose response types and fit you want to inspect."),
          textInput("compare_spocc_dir", "spOccupancy run folder"),
          help_text("Select a spOccupancy run. Occupancy probabilities differ from observed incidence."),
          textInput("compare_sjsdm_dir", "sjSDM run folder"),
          help_text("Select a sjSDM run. Its covariance estimates are specific to the fitted model."),
          textInput("compare_boral_dir", "boral run folder"),
          help_text("Select a boral run. Missing values in summary tables are not results."),
          checkboxGroupInput("compare_engines", "Engines to include",
            choices = names(studio_engine_catalog), selected = names(studio_engine_catalog)),
          help_text("Include only engines whose run evidence you need to inspect."),
          actionButton("compare_run", "Inspect runs", class = "btn-primary"),
          downloadButton("compare_download", "Export full audit CSV", icon = NULL)
        ),
        accordion_panel("Run evidence", DTOutput("compare_table"), verbatimTextOutput("compare_notes")),
        accordion_panel("Workflow capabilities", DTOutput("compare_outputs_matrix"))
      )
    )
  ),

  tabPanel("Guides",
    tags$main(class = "studio-guide",
      section_header("Guides", "Data contracts, interpretation and replay"),
      accordion(
        open = "Input files and identifiers",
        accordion_panel("Input files and identifiers",
          tags$ul(
            tags$li(tags$b("Responses and predictors. "), "Use one row per observation. With IDs in both tables, use the same ID order. Without IDs, Studio pairs rows by position. Add site_id as a first column when possible; do not include it in a model formula."),
            tags$li(tags$b("Species traits. "), "Use one row per response column in the same species order. A first species_id column can identify the rows. Enable trait effects only in a workflow that fits them."),
            tags$li(tags$b("Phylogeny. "), "Use a Newick tree or a taxonomy-rank table. Match tree tips or taxonomy names to the response columns. Studio can derive a relatedness matrix; a taxonomy surrogate is not an evolutionary tree with measured branch lengths."),
            tags$li(tags$b("Spatial inputs. "), "Use one coordinate pair per spatial unit and document the CRS and units. Use projected coordinates for spOccupancy. In Hmsc select longitude/latitude handling explicitly when using geographic degrees."),
            tags$li(tags$b("Repeated visits. "), "An observed nondetection is 0; an unobserved visit is not 0. Site-level incidence alone does not identify imperfect detection.")
          )
        ),
        accordion_panel("Interpret effects, associations and predictions",
          tags$ul(
            tags$li(tags$b("Hmsc Beta and Gamma. "), "Beta is a species-specific coefficient on the selected model scale, not a correlation bounded by -1 and 1. Gamma describes how traits modify those coefficients."),
            tags$li(tags$b("Hmsc rho. "), "The posterior describes phylogenetic dependence in species responses under the fitted model and supplied tree or relatedness structure; it is not a rate of speciation."),
            tags$li(tags$b("Associations. "), "Residual covariance or correlation describes model-dependent co-variation after fitted predictors. It does not demonstrate interactions. Raw latent loadings are not pairwise species correlations."),
            tags$li(tags$b("Prediction. "), "Training fits, predictions at new sites and validation against independently observed outcomes are different. Compare scores only for the same data split, response definition and prediction target."),
            tags$li(tags$b("Occupancy. "), "Latent occupancy, detection given occupancy and observed incidence across visits are different quantities. They should not share an accuracy comparison unless the observation process is explicitly matched."),
            tags$li(tags$b("Universal Benchmark. "), "Simulation modes use documented Studio generating structures, not native simulators for every package. Observed-data mode has no known truth. It skips spOccupancy without genuine repeat-visit data.")
          )
        ),
        accordion_panel("Diagnostics and replay",
          tags$ul(
            tags$li(tags$b("fitted. "), "The fit finished; check convergence, uncertainty and model adequacy before interpreting results."),
            tags$li(tags$b("model_defined. "), "A model description or script was saved without fitted estimates."),
            tags$li(tags$b("check_failed / fit_failed. "), "Use diagnostics to distinguish a data/settings failure from a fitting or post-processing failure."),
            tags$li(tags$b("MCMC and optimization. "), "Check trace plots, ESS and chain agreement for Bayesian models. Thinning reduces stored draws and does not fix poor mixing. For sjSDM, check the optimization and prediction diagnostics instead."),
            tags$li(tags$b("Reproduce. "), "Start with used_config.yml, input snapshots, diagnostics/engine_status.json and the exported replay script. Engine-specific files can be absent or contain NA when an output is unavailable. A ZIP is not proof that the model converged.")
          )
        ),
        accordion_panel("Parameter reference",
          p(class = "small-muted", "The index uses the explanations beside the settings in this version. Conditional controls apply only to their active model choices."),
          downloadButton("parameter_dictionary_download", "Export parameter reference CSV", icon = NULL),
          br(), br(), DTOutput("parameter_dictionary_table")
        ),
        accordion_panel("Package documentation",
          tags$ul(
            tags$li(tags$a(href = "https://cran.r-project.org/package=Hmsc", "Hmsc"), " and ",
                    tags$a(href = "https://github.com/hmsc-r/Hmsc-HPC", "Hmsc-HPC")),
            tags$li(tags$a(href = "https://cran.r-project.org/package=jSDM", "jSDM")),
            tags$li(tags$a(href = "https://cran.r-project.org/package=gjam", "GJAM")),
            tags$li(tags$a(href = "https://doserlab.com/files/spoccupancy-web/", "spOccupancy")),
            tags$li(tags$a(href = "https://cran.r-project.org/package=sjSDM", "sjSDM")),
            tags$li(tags$a(href = "https://cran.r-project.org/package=boral", "boral"))
          )
        )
      )
    )
  )
)

# ------------------------------------------------------------
# Server
# ------------------------------------------------------------

server <- function(input, output, session) {

  studio_bind_preferences(input, session)

  observeEvent(input$home_open_project, {
    updateNavbarPage(session, "main_nav", selected = "Project")
  })
  observeEvent(input$home_open_engine_guide, {
    updateNavbarPage(session, "main_nav", selected = "Engine Guide")
  })

  rv <- reactiveValues(
    hmsc = list(Y=NULL, X=NULL, Tr=NULL, study=NULL, coord=NULL, newdata=NULL, newstudy=NULL, newcoord=NULL, check=NULL, checked=FALSE, log=character(), outdir=NULL, zip=NULL, status="Waiting"),
    hmschpc = list(Y=NULL, X=NULL, Tr=NULL, study=NULL, coord=NULL, taxonomy=NULL, phylo_tree=NULL, newdata=NULL, check=NULL, checked=FALSE, log=character(), outdir=NULL, zip=NULL, status="Waiting"),
    jsdm = list(Y=NULL, X=NULL, Tr=NULL, long=NULL, trials=NULL, check=NULL, checked=FALSE, log=character(), outdir=NULL, zip=NULL, status="Waiting"),
    gjam = list(Y=NULL, X=NULL, types=NULL, censor=NULL, effort=NULL, newdata=NULL, check=NULL, checked=FALSE, log=character(), outdir=NULL, zip=NULL, status="Waiting"),
    spocc = list(y=NULL, occ=NULL, det=NULL, coords=NULL, species=NULL, integrated=NULL, check=NULL, checked=FALSE, log=character(), outdir=NULL, zip=NULL, status="Waiting"),
    sjsdm = list(Y=NULL, env=NULL, spatial=NULL, traits=NULL, newdata=NULL, check=NULL, checked=FALSE, log=character(), outdir=NULL, zip=NULL, status="Waiting"),
    boral = list(Y=NULL, X=NULL, traits=NULL, rowids=NULL, ranefids=NULL, distmat=NULL, offset=NULL, newdata=NULL, check=NULL, checked=FALSE, log=character(), outdir=NULL, zip=NULL, status="Waiting"),
    univ = list(log=character(), data_dir=file.path(app_dir, "examples", "universal_benchmark"), data_source="synthetic", master_dir=NULL, zip=NULL, status="Waiting", preflight=NULL),
    compare = NULL
  )

  add_log <- function(engine, ...) {
    line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", paste(..., collapse = " "))
    if (engine == "hmsc") rv$hmsc$log <- c(rv$hmsc$log, line)
    if (engine == "hmschpc") rv$hmschpc$log <- c(rv$hmschpc$log, line)
    if (engine == "jsdm") rv$jsdm$log <- c(rv$jsdm$log, line)
    if (engine == "gjam") rv$gjam$log <- c(rv$gjam$log, line)
    if (engine == "spocc") rv$spocc$log <- c(rv$spocc$log, line)
    if (engine == "sjsdm") rv$sjsdm$log <- c(rv$sjsdm$log, line)
    if (engine == "boral") rv$boral$log <- c(rv$boral$log, line)
    if (engine == "univ") rv$univ$log <- c(rv$univ$log, line)
  }

  observeEvent(input$hmsc_preset, {
    preset <- as.character(input$hmsc_preset %||% "Custom")
    vals <- switch(
      preset,
      "Quick test" = list(samples = 50, transient = 25, thin = 1, nChains = 2, nParallel = 1, verbose = 10),
      "Standard" = list(samples = 1000, transient = 500, thin = 5, nChains = 4, nParallel = 1, verbose = 100),
      "Publication" = list(samples = 10000, transient = 5000, thin = 10, nChains = 4, nParallel = 1, verbose = 1000),
      NULL
    )
    if (!is.null(vals)) {
      updateNumericInput(session, "hmsc_samples", value = vals$samples)
      updateNumericInput(session, "hmsc_transient", value = vals$transient)
      updateNumericInput(session, "hmsc_thin", value = vals$thin)
      updateNumericInput(session, "hmsc_nChains", value = vals$nChains)
      updateNumericInput(session, "hmsc_nParallel", value = vals$nParallel)
      updateNumericInput(session, "hmsc_verbose", value = vals$verbose)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$jsdm_preset, {
    preset <- switch(input$jsdm_preset,
      "Quick test" = c(burnin = 100L, mcmc = 100L, thin = 1L),
      "Standard" = c(burnin = 2000L, mcmc = 4000L, thin = 10L),
      "Publication" = c(burnin = 10000L, mcmc = 20000L, thin = 10L),
      NULL)
    for (nm in names(preset)) updateNumericInput(session, paste0("jsdm_", nm), value = preset[[nm]])
  }, ignoreInit = TRUE)

  file_upload_specs <- data.frame(
    Engine = c(
      rep("Hmsc", 9),
      rep("Hmsc-HPC", 8),
      rep("jSDM", 7),
      rep("GJAM", 10),
      rep("spOccupancy", 8),
      rep("sjSDM", 9),
      rep("boral", 9)
    ),
    Id = c(
      "hmsc_Y_file", "hmsc_X_file", "hmsc_newdata_file", "hmsc_Tr_file", "hmsc_study_file", "hmsc_newstudy_file", "hmsc_coord_file", "hmsc_newcoords_file", "hmsc_phylo_file",
      "hmschpc_Y_file", "hmschpc_X_file", "hmschpc_traits_file", "hmschpc_study_file", "hmschpc_coord_file", "hmschpc_taxonomy_file", "hmschpc_phylo_tree_file", "hmschpc_newdata_file",
      "jsdm_Y_file", "jsdm_X_file", "jsdm_trait_file", "jsdm_long_file", "jsdm_trials_file", "jsdm_newdata_file", "jsdm_prediction_ids_file",
      "gjam_Y_file", "gjam_X_file", "gjam_type_file", "gjam_censor_file", "gjam_effort_file", "gjam_newdata_file", "gjam_trait_spec_file", "gjam_trait_types_file", "gjam_holdout_file", "gjam_prior_file",
      "spocc_y_file", "spocc_occ_file", "spocc_det_file", "spocc_coords_file", "spocc_species_file", "spocc_integrated_file", "spocc_newdata_file", "spocc_newcoords_file",
      "sjsdm_Y_file", "sjsdm_env_file", "sjsdm_spatial_file", "sjsdm_traits_file", "sjsdm_newdata_file", "sjsdm_group_file", "sjsdm_spatial_new_file", "sjsdm_weights_file", "sjsdm_cv_file",
      "boral_Y_file", "boral_X_file", "boral_traits_file", "boral_rowids_file", "boral_ranefids_file", "boral_distmat_file", "boral_offset_file", "boral_newdata_file", "boral_trials_file"
    ),
    Label = c(
      "Y.csv response matrix", "XData.csv predictors", "newdata.csv prediction covariates", "TrData.csv traits", "studyDesign.csv random effects", "prediction_studyDesign.csv prediction random effects", "coordinates.csv spatial coordinates", "newcoords.csv prediction spatial coordinates", "phylogeny or taxonomy file",
      "Y.csv response matrix", "XData.csv predictors", "traits.csv species traits", "studyDesign.csv random effects", "coordinates.csv spatial coordinates", "taxonomy.csv rank table", "phylo_tree.nwk Newick tree", "newdata.csv prediction covariates",
      "Y.csv response matrix", "XData.csv site data", "trait_data.csv traits", "long_format.csv long data", "trials.csv binomial trials", "newdata.csv prediction covariates", "prediction_ids.csv site/species IDs",
      "Y.csv response matrix", "XData.csv predictors", "typeNames.csv response types", "censor.csv censoring intervals", "effort.csv sampling effort", "newdata.csv prediction covariates", "specByTrait.csv traits", "traitTypes.csv trait response types", "holdoutIndex.csv holdouts", "priorTemplate.csv prior file",
      "y.csv detection data", "occ.covs.csv occurrence covariates", "det.covs.csv detection covariates", "coords.csv spatial coordinates", "species.csv metadata", "integrated_sources.csv data sources", "newdata.csv prediction covariates", "newcoords.csv prediction coordinates",
      "Y.csv response matrix", "env.csv environmental predictors", "spatial.csv spatial predictors", "traits.csv species traits", "newdata.csv prediction covariates", "species_groups.csv groups", "new_spatial.csv prediction spatial data", "weights.rds pretrained weights", "folds.csv validation folds",
      "Y.csv response matrix", "XData.csv covariates", "traits.csv species traits", "row.ids.csv row effects", "ranef.ids.csv random intercept IDs", "distmat.csv distance matrix", "offset.csv offset matrix", "newdata.csv prediction covariates", "trial.size.csv binomial trials"
    ),
    stringsAsFactors = FALSE
  )

  uploaded_files <- reactive({
    rows <- lapply(seq_len(nrow(file_upload_specs)), function(i) {
      spec <- file_upload_specs[i, , drop = FALSE]
      f <- input[[spec$Id]]
      if (is.null(f) || !is.data.frame(f) || nrow(f) < 1) return(NULL)
      data.frame(
        Engine = spec$Engine,
        Id = spec$Id,
        Label = spec$Label,
        Name = as.character(f$name[1]),
        Size = as.numeric(f$size[1]),
        Type = as.character(f$type[1] %||% ""),
        Datapath = as.character(f$datapath[1]),
        stringsAsFactors = FALSE
      )
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) {
      return(data.frame(Engine=character(), Id=character(), Label=character(), Name=character(), Size=numeric(), Type=character(), Datapath=character(), stringsAsFactors = FALSE))
    }
    do.call(rbind, rows)
  })

  preview_uploaded_file_data <- function(row) {
    if (is.null(row) || !nrow(row)) return(data.frame(Message = "No uploaded file available for preview.", stringsAsFactors = FALSE))
    path <- row$Datapath
    ext <- tolower(tools::file_ext(row$Name))
    dat <- tryCatch({
      if (!file.exists(path)) {
        data.frame(Message = "Uploaded temporary file is no longer available.", stringsAsFactors = FALSE)
      } else if (ext %in% c("csv")) {
        read.csv(path, nrows = 100, check.names = FALSE, stringsAsFactors = FALSE)
      } else if (ext %in% c("tsv", "tab")) {
        read.delim(path, nrows = 100, check.names = FALSE, stringsAsFactors = FALSE)
      } else if (ext %in% c("txt", "tre", "tree", "nwk", "log", "yml", "yaml", "json")) {
        lines <- readLines(path, n = 100, warn = FALSE, encoding = "UTF-8")
        data.frame(Line = seq_along(lines), Text = lines, stringsAsFactors = FALSE)
      } else {
        data.frame(Message = paste0("Preview is not parsed for .", ext, " files. The file is still available to the workflow."), stringsAsFactors = FALSE)
      }
    }, error = function(e) {
      data.frame(Message = paste("Preview failed:", conditionMessage(e)), stringsAsFactors = FALSE)
    })
    if (!nrow(dat)) dat <- data.frame(Message = "The file was readable but has no rows to preview.", stringsAsFactors = FALSE)
    dat
  }

  register_engine_file_preview <- function(prefix, engine_name) {
    local({
      p <- prefix
      eng <- engine_name
      choice_id <- paste0(p, "_file_preview_choice")
      selector_id <- paste0(p, "_file_preview_selector")
      meta_id <- paste0(p, "_file_preview_meta")
      table_id <- paste0(p, "_file_preview_table")

      engine_files <- reactive({
        dat <- uploaded_files()
        if (!nrow(dat)) return(dat)
        dat[tolower(dat$Engine) == tolower(eng), , drop = FALSE]
      })

      selected_engine_file <- reactive({
        dat <- engine_files()
        if (!nrow(dat)) return(NULL)
        idx <- suppressWarnings(as.integer(input[[choice_id]] %||% 1))
        if (!is.finite(idx) || idx < 1 || idx > nrow(dat)) idx <- 1
        dat[idx, , drop = FALSE]
      })

      output[[selector_id]] <- renderUI({
        dat <- engine_files()
        if (!nrow(dat)) {
          return(div(class="note", tags$b("No uploaded files yet."), paste("Upload", eng, "files in step 1, then preview them here before checking or fitting.")))
        }
        labels <- paste(dat$Label, dat$Name, sep = " - ")
        selectInput(choice_id, "Uploaded file", choices = stats::setNames(seq_len(nrow(dat)), labels), selected = 1)
      })

      output[[meta_id]] <- renderText({
        row <- selected_engine_file()
        if (is.null(row)) return(paste("No", eng, "uploaded file selected."))
        paste(
          paste("Engine:", row$Engine),
          paste("Input:", row$Label),
          paste("File:", row$Name),
          paste("Size:", format(row$Size, big.mark = ","), "bytes"),
          paste("Temporary path:", row$Datapath),
          sep = "\n"
        )
      })

      output[[table_id]] <- renderDT({
        row <- selected_engine_file()
        datatable(preview_uploaded_file_data(row), options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
      })
    })
  }

  register_engine_file_preview("hmsc", "Hmsc")
  register_engine_file_preview("hmschpc", "Hmsc-HPC")
  register_engine_file_preview("jsdm", "jSDM")
  register_engine_file_preview("gjam", "GJAM")
  register_engine_file_preview("spocc", "spOccupancy")
  register_engine_file_preview("sjsdm", "sjSDM")
  register_engine_file_preview("boral", "boral")

  output$comparison_meaning_table <- renderDT({
    dat <- data.frame(
      Item = c("Run status", "Status meaning", "Runtime", "Prediction performance", "Main response directions", "Residual associations", "Output completeness", "Raw parameters", "Failed runs", "Scientific adequacy"),
      Compare = c("Yes", "Interpret before comparing", "Yes, as practical cost", "Only under matched validation", "Cautiously", "Qualitatively only", "Yes", "No direct equality", "Diagnostics only", "Never inferred from status alone"),
      Reason = c(
        "fitted, model_defined, check_failed and fit_failed have different meanings and are reported as workflow outcomes.",
        "fitted means fitted outputs were produced; model_defined means a model boundary exists without fitted posterior evidence; check_failed and fit_failed are failure states with diagnostic value.",
        "Runtime is engine-specific but comparable as practical cost.",
        "Comparable only when the same training/testing split, response definition, response scale and metric definition were used.",
        "Signs and broad patterns can be compared only after checking link functions, scaling and response family.",
        "Hmsc Omega, Hmsc-HPC Eta/Lambda draws, jSDM residual correlations, GJAM corMu/sigMu, spOccupancy latent factors, sjSDM covariance and boral residual correlations are related but not numerically interchangeable.",
        "The comparison panel checks used_config.yml, standard tables, workflow maps, reports, table counts, plot counts and diagnostics counts.",
        "Different parameterizations, priors, links and latent structures.",
        "check_failed and fit_failed are valid audit folders. They are not fitted analyses and should point to diagnostics files.",
        "Convergence, identifiability, data design, family choice and ecological interpretation must still be reviewed after a fitted status."
      )
    )
    datatable(dat, options = list(dom="tip", pageLength=10), rownames=FALSE)
  })


  project_profile <- reactive({
    structure <- input$project_response_structure %||% "unknown"
    list(
      n_sites = input$project_n_sites, n_responses = input$project_n_responses,
      n_predictors = input$project_n_predictors, response_structure = structure,
      response_family = if (identical(structure, "single response type")) input$project_response_family %||% "unknown",
      has_traits = isTRUE(input$project_has_traits), has_phylogeny = isTRUE(input$project_has_phylogeny),
      has_spatial = isTRUE(input$project_has_spatial), many_zeros = isTRUE(input$project_many_zeros),
      need_prediction = isTRUE(input$project_need_prediction),
      has_integrated_sources = identical(input$question_template, "Integrated occupancy with multiple data sources") &&
        isTRUE(input$project_has_integrated_sources)
    )
  })

  project_recommendation <- reactive({
    studio_project_recommend(project_profile(), input$question_template)
  })

  output$project_recommendation_summary <- renderUI({
    result <- project_recommendation()
    leading <- result$top
    div(class = "studio-match-summary", role = "status",
      `data-leading-engines` = paste(leading$Engine, collapse = ","),
      if (nrow(leading)) tagList(
        tags$strong(if (nrow(leading) == 1L) "Recommended starting point" else "Tied starting points"),
        p(sprintf("%s (%d pts)", paste(leading$Engine, collapse = " / "), leading$Score[[1L]])),
        if (nrow(leading) > 1L) p(class = "studio-match-note", "The scores tie; compare the model structures below.")
      ) else tagList(tags$strong("No primary match yet"), p(result$message)),
      if (length(result$notices)) lapply(result$notices, function(note) p(class = "studio-match-note", note)),
      tags$a(href = "#project_matches", class = "studio-match-jump", "View all match scores")
    )
  })

  output$project_recommendation_cards <- renderUI({
    result <- project_recommendation()
    rows <- result$rows
    div(class = "studio-project-recs", lapply(seq_len(nrow(rows)), function(i) {
      studio_project_match_row(rows[i, , drop = FALSE], leading = rows$Engine[i] %in% result$top$Engine)
    }))
  })

  output$project_metadata_table <- renderDT({
    fields <- c(list(Project = input$project_name, "Question template" = input$question_template), project_profile())
    values <- vapply(fields, function(value) as.character(value %||% "Not specified")[[1L]], character(1))
    datatable(data.frame(Field = gsub("_", " ", names(fields)), Value = unname(values)),
              options = list(dom = "tip", pageLength = 15), rownames = FALSE)
  })

  # Universal Benchmark
  univ_runner <- file.path(app_dir, "workflow_scripts", "universal_benchmark_runner.R")
  univ_rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

  univ_selected_engines <- reactive({
    input$univ_engines %||% c("Hmsc", "Hmsc-HPC", "jSDM", "GJAM", "spOccupancy", "sjSDM", "boral")
  })

  univ_generated_dir <- tempfile("studio_simulation_")

  univ_required_files <- function() {
    c("Y_occurrence.csv", "Y_count.csv", "Y_normal.csv", "XData.csv", "traits.csv",
      "studyDesign.csv", "coordinates.csv", "phylogeny.nwk", "taxonomy.csv",
      "newdata.csv", "newcoords.csv", "folds.csv", "trial.size.csv", "offset.csv",
      "row.ids.csv", "ranef.ids.csv", "distmat.csv", "spOccupancy_y_detection.csv",
      "occ.covs.csv", "det.covs.csv", file.path("truth", "latent_occurrence.csv"),
      file.path("truth", "occurrence_probability.csv"),
      file.path("truth", "true_environment_effects.csv"),
      file.path("truth", "true_species_associations.csv"),
      file.path("truth", "true_predictions.csv"),
      file.path("truth", "data_generation_config.yml"),
      file.path("truth", "data_generation_script.R"))
  }

  univ_validate_benchmark_dir <- function(dir) {
    if (is.null(dir) || !dir.exists(dir)) {
      return(list(ok = FALSE, missing = univ_required_files(), root = dir %||% ""))
    }
    req <- univ_required_files()
    missing <- req[!file.exists(file.path(dir, req))]
    list(ok = length(missing) == 0L, missing = missing,
         root = normalizePath(dir, winslash = "/", mustWork = FALSE))
  }

  univ_find_benchmark_root <- function(dir) {
    direct <- univ_validate_benchmark_dir(dir)
    if (isTRUE(direct$ok)) return(direct$root)
    kids <- list.dirs(dir, recursive = FALSE, full.names = TRUE)
    hits <- vapply(kids, function(k) isTRUE(univ_validate_benchmark_dir(k)$ok), logical(1))
    if (any(hits)) return(normalizePath(kids[which(hits)[1]], winslash = "/", mustWork = FALSE))
    stop("Uploaded ZIP does not contain a complete Universal Benchmark input contract. Missing examples include: ",
         paste(utils::head(direct$missing, 12), collapse = ", "), call. = FALSE)
  }

  univ_builtin_case_root <- function() file.path(app_dir, "examples", "universal_benchmark_cases")
  univ_builtin_cases <- function() {
    root <- univ_builtin_case_root()
    if (!dir.exists(root)) return(character())
    dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
    dirs <- dirs[vapply(dirs, function(d) isTRUE(univ_validate_benchmark_dir(d)$ok), logical(1))]
    stats::setNames(normalizePath(dirs, winslash = "/", mustWork = FALSE), basename(dirs))
  }

  run_univ_runner <- function(action, extra = character(), progress_message = NULL) {
    if (!file.exists(univ_runner)) stop("Universal benchmark runner is missing: ", univ_runner, call. = FALSE)
    if (action %in% c("all", "empirical")) {
      rv$univ$master_dir <- rv$univ$zip <- NULL
    }
    args <- if (identical(action, "empirical")) {
      c(shQuote(univ_runner), "--action=empirical", extra)
    } else c(
      shQuote(univ_runner),
      paste0("--action=", action),
      paste0("--n_sites=", as.integer(input$univ_n_sites %||% 120)),
      paste0("--n_species=", as.integer(input$univ_n_species %||% 12)),
      paste0("--n_visits=", as.integer(input$univ_n_visits %||% 3)),
      paste0("--seed=", as.integer(input$univ_seed %||% 20260601)),
      paste0("--engines=", paste(univ_selected_engines(), collapse = ",")),
      paste0("--allow_large=", tolower(as.character(isTRUE(input$univ_allow_large)))),
      paste0("--generating_structure=", input$univ_generating_structure %||% "probit_latent"),
      paste0("--truth_profile=", input$univ_truth_profile %||% "recovery"),
      paste0("--effect_scale=", as.numeric(input$univ_effect_scale %||% 1.6)),
      paste0("--association_strength=", as.numeric(input$univ_association_strength %||% 1.65)),
      paste0("--residual_sd=", as.numeric(input$univ_residual_sd %||% 1.05)),
      paste0("--observation_noise_sd=", as.numeric(input$univ_observation_noise_sd %||% 0.15)),
      paste0("--hmsc_samples=", as.integer(input$univ_hmsc_samples %||% 500)),
      paste0("--hmsc_transient=", as.integer(input$univ_hmsc_transient %||% 500)),
      paste0("--hmsc_thin=", as.integer(input$univ_hmsc_thin %||% 1)),
      paste0("--hmsc_nChains=", as.integer(input$univ_hmsc_nChains %||% 2)),
      paste0("--hmschpc_samples=", as.integer(input$univ_hmschpc_samples %||% 500)),
      paste0("--hmschpc_transient=", as.integer(input$univ_hmschpc_transient %||% 500)),
      paste0("--hmschpc_thin=", as.integer(input$univ_hmschpc_thin %||% 1)),
      paste0("--hmschpc_chains=", as.integer(input$univ_hmschpc_chains %||% 1)),
      paste0("--jsdm_burnin=", as.integer(input$univ_jsdm_burnin %||% 200)),
      paste0("--jsdm_mcmc=", as.integer(input$univ_jsdm_mcmc %||% 300)),
      paste0("--jsdm_thin=", as.integer(input$univ_jsdm_thin %||% 5)),
      paste0("--gjam_ng=", as.integer(input$univ_gjam_ng %||% 700)),
      paste0("--gjam_burnin=", as.integer(input$univ_gjam_burnin %||% 150)),
      paste0("--spocc_n_batch=", as.integer(input$univ_spocc_n_batch %||% 30)),
      paste0("--spocc_batch_length=", as.integer(input$univ_spocc_batch_length %||% 10)),
      paste0("--spocc_n_burn=", as.integer(input$univ_spocc_n_burn %||% 100)),
      paste0("--spocc_n_thin=", as.integer(input$univ_spocc_n_thin %||% 5)),
      paste0("--spocc_n_chains=", as.integer(input$univ_spocc_n_chains %||% 1)),
      paste0("--sjsdm_iter=", as.integer(input$univ_sjsdm_iter %||% 100)),
      paste0("--sjsdm_sampling=", as.integer(input$univ_sjsdm_sampling %||% 1000)),
      paste0("--sjsdm_step_size=", as.integer(input$univ_sjsdm_step_size %||% 32)),
      paste0("--sjsdm_learning_rate=", as.numeric(input$univ_sjsdm_learning_rate %||% 0.010)),
      paste0("--sjsdm_weight_decay=", as.numeric(input$univ_sjsdm_weight_decay %||% 0)),
      paste0("--sjsdm_biotic_lambda=", as.numeric(input$univ_sjsdm_biotic_lambda %||% 0.001)),
      paste0("--boral_n_burnin=", as.integer(input$univ_boral_n_burnin %||% 250)),
      paste0("--boral_n_iteration=", as.integer(input$univ_boral_n_iteration %||% 750)),
      paste0("--boral_n_thin=", as.integer(input$univ_boral_n_thin %||% 10)),
      extra
    )
    stdout <- tempfile("universal_benchmark_stdout_", fileext = ".log")
    stderr <- tempfile("universal_benchmark_stderr_", fileext = ".log")
    add_log("univ", "Running universal benchmark action:", action)
    if (!is.null(progress_message)) incProgress(0.1, detail = progress_message)
    neutralize_windows_r_locale_env()
    code <- tryCatch(system2(univ_rscript, args = args, stdout = stdout, stderr = stderr), error = function(e) {
      writeLines(conditionMessage(e), stderr)
      1L
    })
    out_lines <- if (file.exists(stdout)) readLines(stdout, warn = FALSE) else character()
    err_lines <- if (file.exists(stderr)) readLines(stderr, warn = FALSE) else character()
    matrix_line <- grep("^DESIGN_MATRIX_DIR=", out_lines, value = TRUE)
    if (length(matrix_line)) rv$univ$matrix_dir <- sub("^DESIGN_MATRIX_DIR=", "", tail(matrix_line, 1))
    if (length(out_lines)) add_log("univ", paste(tail(out_lines, 25), collapse = "\n"))
    if (length(err_lines)) add_log("univ", paste(tail(err_lines, 25), collapse = "\n"))
    if (!identical(as.integer(code), 0L)) {
      rv$univ$status <- "fit_failed"
      add_log("univ", "Universal benchmark action failed with exit code:", code)
    } else {
      reported_path <- function(key) {
        hit <- grep(paste0("^", key, "="), out_lines, value = TRUE)
        if (length(hit)) sub(paste0("^", key, "="), "", tail(hit, 1))
      }
      master <- reported_path("MASTER_DIR")
      archive <- reported_path("MASTER_ZIP")
      if (!is.null(master) && dir.exists(master)) {
        rv$univ$master_dir <- master
        rv$univ$result_kind <- if (action == "empirical") "empirical" else "simulation"
        rv$univ$status <- "Run complete. Review per-engine status and diagnostics."
      }
      if (!is.null(archive) && file.exists(archive)) rv$univ$zip <- archive
      add_log("univ", "Universal benchmark action completed:", action)
    }
    invisible(code)
  }

  output$univ_builtin_case_ui <- renderUI({
    cases <- univ_builtin_cases()
    if (!length(cases)) {
      return(div(class = "warn", "No built-in real-style benchmark cases have been generated yet."))
    }
    selectInput("univ_builtin_case", "Built-in real-style example case",
                choices = cases, selected = cases[[1]])
  })

  univ_preflight_df <- function() {
    pkgs <- c("shiny", "DT", "yaml", "jsonlite", "zip", "Hmsc", "coda", "ape", "jSDM",
              "gjam", "spOccupancy", "sjSDM", "reticulate", "boral", "rjags", "R2jags")
    rows <- lapply(pkgs, function(pkg) {
      ok <- requireNamespace(pkg, quietly = TRUE)
      data.frame(Component = pkg, Kind = "R package", Available = ok,
                 Version = if (ok) tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "") else "",
                 Note = if (ok) "available" else "not installed or not loadable",
                 stringsAsFactors = FALSE)
    })
    py <- hmschpc_python_bin(list())
    if (!nzchar(py) || !file.exists(py)) py <- Sys.which("python")
    rows[[length(rows) + 1L]] <- data.frame(Component = "JSDM Studio shared Python", Kind = "runtime",
                                            Available = nzchar(py), Version = py,
                                            Note = if (nzchar(py)) "used by Hmsc-HPC and sjSDM" else "not found; run Full Offline Setup",
                                            stringsAsFactors = FALSE)
    python_module_row <- function(module, import_name = module, engine = "Hmsc-HPC", needs_source = FALSE) {
      if (!nzchar(py)) {
        return(data.frame(Component = paste0(engine, " Python:", module),
                          Kind = "Python package", Available = FALSE, Version = "",
                          Note = "Python runtime not found", stringsAsFactors = FALSE))
      }
      source_dir <- normalizePath(file.path(app_dir, "external_packages", "hmsc-hpc-main"),
                                  winslash = "/", mustWork = FALSE)
      source_stmt <- if (isTRUE(needs_source)) sprintf("sys.path.insert(0, %s);", shQuote(source_dir)) else ""
      code <- sprintf(
        "import os, sys, importlib; os.environ['TF_CPP_MIN_LOG_LEVEL']='2'; %s m=importlib.import_module(%s); print(getattr(m, '__version__', 'available'))",
        source_stmt, shQuote(import_name)
      )
      out <- tryCatch(system2(py, c("-c", shQuote(code)), stdout = TRUE, stderr = FALSE),
                      error = function(e) structure(conditionMessage(e), status = 1L))
      ok <- is.null(attr(out, "status")) || identical(attr(out, "status"), 0L)
      data.frame(Component = paste0(engine, " Python:", module),
                 Kind = "Python package", Available = ok,
                 Version = if (ok && length(out)) trimws(out[[1]]) else "",
                 Note = if (ok) "available" else "run Full Offline Setup or inspect dependency diagnostics",
                 stringsAsFactors = FALSE)
    }
    for (item in list(c("numpy", "numpy"), c("pandas", "pandas"), c("patsy", "patsy"),
                      c("PyYAML", "yaml"), c("h5py", "h5py"), c("scipy", "scipy"),
                      c("tensorflow", "tensorflow"), c("tensorflow_probability", "tensorflow_probability"),
                      c("ujson", "ujson"), c("pyhmsc", "pyhmsc"))) {
      rows[[length(rows) + 1L]] <- python_module_row(item[[1]], item[[2]], needs_source = TRUE)
    }
    for (item in list(c("torch", "torch"), c("torch_optimizer", "torch_optimizer"),
                      c("pyro", "pyro"), c("madgrad", "madgrad"))) {
      rows[[length(rows) + 1L]] <- python_module_row(item[[1]], item[[2]], engine = "sjSDM")
    }
    rows[[length(rows) + 1L]] <- data.frame(Component = "JAGS", Kind = "system dependency",
                                            Available = requireNamespace("rjags", quietly = TRUE),
                                            Version = if (requireNamespace("rjags", quietly = TRUE)) "rjags loadable" else "",
                                            Note = "boral requires system JAGS plus rjags/R2jags",
                                            stringsAsFactors = FALSE)
    do.call(rbind, rows)
  }

  observeEvent(input$univ_generate, {
    code <- withProgress(message = "Generating simulation preview", value = 0, {
      status <- run_univ_runner("generate", extra = paste0("--target_dir=", shQuote(univ_generated_dir)),
                                progress_message = "Writing a session-local preview")
      incProgress(1)
      status
    })
    if (code != 0L) return()
    rv$univ$data_dir <- univ_generated_dir
    rv$univ$data_source <- "synthetic"
    rv$univ$status <- "Data generated"
    updateRadioButtons(session, "univ_data_source", selected = "synthetic")
  })

  observeEvent(input$univ_use_builtin, {
    cases <- univ_builtin_cases()
    selected <- input$univ_builtin_case %||% ""
    if (!selected %in% unname(cases) || !dir.exists(selected)) {
      showNotification("No valid built-in benchmark case is selected.", type = "error")
      return()
    }
    chk <- univ_validate_benchmark_dir(selected)
    if (!isTRUE(chk$ok)) {
      showNotification(paste("Selected case is incomplete:", paste(chk$missing, collapse = ", ")), type = "error")
      return()
    }
    rv$univ$data_dir <- chk$root
    rv$univ$data_source <- "builtin"
    rv$univ$status <- paste("Using built-in case", basename(chk$root))
    updateRadioButtons(session, "univ_data_source", selected = "builtin")
    add_log("univ", "Using built-in benchmark case:", chk$root)
  })

  observeEvent(input$univ_import_zip, {
    if (is.null(input$univ_data_zip) || !file.exists(input$univ_data_zip$datapath)) {
      showNotification("Please upload a benchmark ZIP first.", type = "error")
      return()
    }
    import_root <- tempfile("studio_benchmark_zip_")
    dir.create(import_root, recursive = TRUE, showWarnings = FALSE)
    ok <- tryCatch({
      archive <- input$univ_data_zip$datapath
      contents <- zip::zip_list(archive)
      paths <- contents$filename
      if (!length(paths) || length(paths) > 5000L ||
          any(grepl("(^/|\\\\|:|(^|/)\\.\\.?(/|$))", paths)) ||
          anyDuplicated(tolower(paths)) ||
          any(grepl("^l", contents$permissions)) ||
          sum(contents$uncompressed_size) > 256 * 1024^2) {
        stop("Rejecting an unsafe or oversized benchmark ZIP.", call. = FALSE)
      }
      zip::unzip(archive, exdir = import_root)
      root <- univ_find_benchmark_root(import_root)
      rv$univ$data_dir <- root
      rv$univ$data_source <- "upload"
      rv$univ$status <- "Uploaded benchmark ZIP imported"
      updateRadioButtons(session, "univ_data_source", selected = "upload")
      add_log("univ", "Imported benchmark ZIP:", input$univ_data_zip$name)
      add_log("univ", "Using uploaded benchmark input directory:", root)
      TRUE
    }, error = function(e) {
      add_log("univ", "Benchmark ZIP import failed:", conditionMessage(e))
      showNotification(conditionMessage(e), type = "error", duration = 10)
      FALSE
    })
    invisible(ok)
  })

  observeEvent(input$univ_preflight, {
    rv$univ$preflight <- univ_preflight_df()
    rv$univ$status <- "Preflight complete"
    add_log("univ", "Dependency preflight table refreshed.")
  })

  stage_empirical <- function() {
    dir <- tempfile("studio_observed_data_")
    dir.create(dir)
    uploads <- list(y = input$univ_empirical_y, x = input$univ_empirical_x,
      visits = input$univ_empirical_visits, new_x = input$univ_empirical_new_x, test_y = input$univ_empirical_test_y)
    for (nm in names(empirical_input_names)) {
      upload <- uploads[[nm]]
      if (!is.null(upload) && !file.copy(upload$datapath, file.path(dir, empirical_input_names[[nm]]))) {
        stop("Could not stage ", empirical_input_names[[nm]], ".", call. = FALSE)
      }
    }
    dat <- load_empirical_benchmark(dir)
    rv$univ$data_dir <- dir
    rv$univ$data_source <- "empirical"
    rv$univ$empirical_summary <- paste0(
      nrow(dat$Y), " training sites | ", ncol(dat$Y), " species | ", ncol(dat$X), " design columns\n",
      NROW(dat$newdata), " new sites | ", if (is.null(dat$test_y)) "No independent outcomes" else "Independent observed outcomes supplied", "\n",
      if (is.null(dat$det)) "spOccupancy: not applicable without observed repeat visits."
      else paste0("spOccupancy: ", min(dat$visit_count), "-", max(dat$visit_count), " observed visits per site."))
    dir
  }

  check_empirical <- function() {
    tryCatch(stage_empirical(), error = function(e) {
      rv$univ$status <- "check_failed"
      rv$univ$empirical_summary <- conditionMessage(e)
      showNotification(conditionMessage(e), type = "error", duration = 10)
      NULL
    })
  }

  observeEvent(input$univ_empirical_check, {
    if (identical(input$univ_data_source, "empirical") && !is.null(check_empirical())) {
      rv$univ$status <- "Observed inputs passed checks"
    }
  })

  observeEvent(input$univ_run_all, {
    rv$univ$status <- "Running"
    mode <- input$univ_data_source %||% "empirical"
    if (identical(mode, "empirical")) {
      dir <- check_empirical()
      if (is.null(dir)) return()
      profile <- input$univ_empirical_profile %||% "standard"
      empirical_profile(profile)
      withProgress(message = "Fitting observed data", value = 0, {
        run_univ_runner("empirical", extra = c(paste0("--profile=", profile),
                         paste0("--input_dir=", shQuote(dir))),
                         progress_message = "Running eligible workflows in sequence")
        incProgress(1)
      })
      return()
    }
    extra <- paste0("--target_dir=", shQuote(univ_generated_dir))
    if (!identical(mode, "synthetic")) {
      if (!identical(mode, rv$univ$data_source)) {
        showNotification("Load a demonstration case or import a valid ZIP first.", type = "error")
        rv$univ$status <- "check_failed"
        return()
      }
      chk <- univ_validate_benchmark_dir(rv$univ$data_dir)
      if (!isTRUE(chk$ok)) {
        showNotification(paste("Current benchmark input is incomplete:", paste(chk$missing, collapse = ", ")), type = "error", duration = 10)
        rv$univ$status <- "check_failed"
        return()
      }
      extra <- paste0("--input_dir=", shQuote(chk$root))
    } else {
      rv$univ$data_dir <- univ_generated_dir
      rv$univ$data_source <- "synthetic"
    }
    withProgress(message = "Running universal benchmark", value = 0, {
      run_univ_runner("all", extra = extra, progress_message = "Sequential all-engine run started")
      incProgress(1)
    })
  })

  observeEvent(input$univ_run_matrix, {
    if (!identical(input$univ_data_source %||% "synthetic", "synthetic")) {
      showNotification("The generating-structure matrix requires synthetic mode. Select synthetic data first.", type = "error", duration = 10)
      return()
    }
    rv$univ$status <- "Running design matrix"
    withProgress(message = "Running both generating structures", value = 0, {
      run_univ_runner("matrix", progress_message = "Sequential structure-by-engine fits started")
      incProgress(1)
    })
  })

  output$univ_matrix_summary <- renderText({
    d <- rv$univ$matrix_dir %||% ""
    if (!nzchar(d)) return("No two-structure matrix has been run in this session.")
    p <- file.path(d, "design_matrix_status.csv")
    if (!file.exists(p)) return(paste("Matrix directory:", d, "(status table missing)"))
    x <- tryCatch(read.csv(p, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(x)) return(paste("Matrix directory:", d, "(status table unreadable)"))
    paste("Matrix directory:", d, "\nEngine-by-structure rows:", nrow(x),
          "\nFitted:", sum(x$status == "fitted", na.rm = TRUE),
          "\nOther statuses:", sum(x$status != "fitted", na.rm = TRUE))
  })

  observeEvent(input$univ_compare, {
    if (identical(rv$univ$result_kind, "empirical")) {
      showNotification("The current run used observed data; its comparison is already saved.", type = "message")
      return()
    }
    if (!is.null(rv$univ$master_dir) && dir.exists(rv$univ$master_dir)) {
      withProgress(message = "Refreshing benchmark comparison", value = 0, {
        run_univ_runner("compare", extra = paste0("--master_dir=", shQuote(rv$univ$master_dir)),
                        progress_message = "Rebuilding master comparison tables")
        incProgress(1)
      })
    } else {
      add_log("univ", "No universal benchmark master folder exists yet.")
    }
  })

  output$univ_log <- renderText(paste(rv$univ$log, collapse = "\n"))

  output$univ_data_summary <- renderText({
    if (identical(input$univ_data_source, "empirical")) {
      return(rv$univ$empirical_summary %||% "No observed data checked yet.")
    }
    d <- rv$univ$data_dir %||% file.path(app_dir, "examples", "universal_benchmark")
    y <- file.path(d, "Y_occurrence.csv")
    x <- file.path(d, "XData.csv")
    tr <- file.path(d, "traits.csv")
    if (!file.exists(y)) return("Benchmark data have not been generated yet.")
    yy <- tryCatch(read.csv(y, row.names = 1, check.names = FALSE), error = function(e) NULL)
    xx <- tryCatch(read.csv(x, row.names = 1, check.names = FALSE), error = function(e) NULL)
    tt <- tryCatch(read.csv(tr, row.names = 1, check.names = FALSE), error = function(e) NULL)
    generation_cfg <- tryCatch(yaml::read_yaml(file.path(d, "truth", "data_generation_config.yml")), error = function(e) list())
    paste0(
      "Data folder: ", normalizePath(d, winslash = "/", mustWork = FALSE), "\n",
      "Data source: ", rv$univ$data_source %||% input$univ_data_source %||% "synthetic", "\n",
      "Generating structure: ", generation_cfg$generating_structure %||% "not recorded", "\n",
      "Native package simulation: no\n",
      "Sites: ", if (is.null(yy)) NA else nrow(yy), "\n",
      "Species: ", if (is.null(yy)) NA else ncol(yy), "\n",
      "Predictor columns: ", if (is.null(xx)) NA else ncol(xx), "\n",
      "Trait columns: ", if (is.null(tt)) NA else ncol(tt), "\n",
      "Truth files: ", length(list.files(file.path(d, "truth"), recursive = TRUE))
    )
  })

  output$univ_required_files <- renderDT({
    dat <- data.frame(
      File = univ_required_files(),
      Role = c(
        "Presence/absence community matrix", "Count response matrix", "Continuous response matrix",
        "Site environmental predictors", "Species traits", "Sampling design and grouping IDs",
        "Site coordinates", "Newick phylogeny", "Taxonomy rank table; correlation derived internally",
        "Prediction covariates", "Prediction coordinates", "Cross-validation folds",
        "Binomial trial sizes", "Offsets", "boral row effect IDs", "boral random-effect IDs",
        "Distance matrix", "Replicated detection-nondetection table", "Occupancy covariates",
        "Detection covariates", "Latent occurrence states", "Truth occurrence probabilities",
        "Truth environment effects", "Truth species associations", "Truth predictions",
        "Generation/config metadata", "Reproducible data-generation script"
      ),
      stringsAsFactors = FALSE
    )
    datatable(dat, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$univ_input_files <- renderDT({
    d <- rv$univ$data_dir %||% file.path(app_dir, "examples", "universal_benchmark")
    if (!dir.exists(d)) return(datatable(data.frame(Message = "Benchmark data folder does not exist yet."), options = list(dom = "tip"), rownames = FALSE))
    files <- list.files(d, recursive = TRUE, full.names = TRUE)
    dat <- data.frame(File = gsub(paste0("^", normalizePath(d, winslash="/"), "/?"), "", normalizePath(files, winslash="/", mustWork=FALSE)),
                      Size_KB = round(file.info(files)$size / 1024, 2),
                      stringsAsFactors = FALSE)
    datatable(dat, options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
  })

  output$univ_preflight_table <- renderDT({
    dat <- rv$univ$preflight
    if (is.null(dat)) dat <- data.frame(Status = "Dependency preflight has not been run.")
    datatable(dat, options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
  })

  output$univ_run_summary <- renderText({
    paste0(
      "Status: ", rv$univ$status %||% "Waiting", "\n",
      "Result mode: ", rv$univ$result_kind %||% "None yet", "\n",
      "Benchmark data: ", rv$univ$data_dir %||% "None yet", "\n",
      "Master folder: ", rv$univ$master_dir %||% "None yet", "\n",
      "Master ZIP: ", rv$univ$zip %||% "None yet"
    )
  })

  univ_read_master_table <- function(rel, fallback = "No universal benchmark output has been generated yet.") {
    if (is.null(rv$univ$master_dir)) return(data.frame(Message = fallback, stringsAsFactors = FALSE))
    p <- file.path(rv$univ$master_dir, rel)
    if (!file.exists(p)) return(data.frame(Message = paste("Missing:", rel), stringsAsFactors = FALSE))
    tryCatch(read.csv(p, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) data.frame(Message = conditionMessage(e), stringsAsFactors = FALSE))
  }

  output$univ_status_table <- renderDT({
    name <- if (identical(rv$univ$result_kind, "empirical")) "engine_status_matrix.csv" else "model_status_matrix.csv"
    datatable(univ_read_master_table(file.path("04_unified_standard_results", name)),
              options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$univ_effects_table <- renderDT({
    dat <- univ_read_master_table(file.path("04_unified_standard_results", "effects_all_models_long.csv"))
    datatable(dat, filter = "top", options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
  })

  output$univ_predictions_table <- renderDT({
    dat <- univ_read_master_table(file.path("04_unified_standard_results", "predictions_all_models_long.csv"))
    datatable(dat, filter = "top", options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
  })

  output$univ_associations_table <- renderDT({
    dat <- univ_read_master_table(file.path("04_unified_standard_results", "associations_all_models_long.csv"))
    datatable(dat, filter = "top", options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
  })

  output$univ_detection_table <- renderDT({
    if (!identical(rv$univ$result_kind, "empirical")) {
      return(datatable(data.frame(Status = "No observed-data detection model has been fitted."), rownames = FALSE))
    }
    datatable(univ_read_master_table(file.path("04_unified_standard_results", "detection_effects.csv")),
              options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
  })

  output$univ_compare_table <- renderDT({
    name <- if (identical(rv$univ$result_kind, "empirical")) "universal_benchmark_comparison.csv" else "compare_models_summary.csv"
    datatable(univ_read_master_table(file.path("06_compare_models", name)),
              options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$univ_comparable_table <- renderDT({
    name <- if (identical(rv$univ$result_kind, "empirical")) "comparability_rules.csv" else "comparable_results_matrix.csv"
    datatable(univ_read_master_table(file.path("06_compare_models", name)),
              options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$univ_files <- renderDT({
    if (is.null(rv$univ$master_dir) || !dir.exists(rv$univ$master_dir)) {
      return(datatable(data.frame(File = character(), Size_KB = numeric()), options = list(dom = "tip"), rownames = FALSE))
    }
    files <- list.files(rv$univ$master_dir, recursive = TRUE, full.names = TRUE)
    dat <- data.frame(File = gsub(paste0("^", normalizePath(rv$univ$master_dir, winslash="/"), "/?"), "", normalizePath(files, winslash="/", mustWork=FALSE)),
                      Size_KB = round(file.info(files)$size / 1024, 1),
                      stringsAsFactors = FALSE)
    datatable(dat, options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
  })

  output$univ_download <- downloadHandler(
    filename = function() if (is.null(rv$univ$zip)) paste0("JSDMStudio_universal_benchmark_", timestamp_id(), ".zip") else basename(rv$univ$zip),
    content = function(file) {
      copy_zip_to_download(rv$univ$zip, rv$univ$master_dir, "Universal_Benchmark", file, input$project_name)
    }
  )


  # Hmsc check
  hmsc_group_rows <- reactiveVal(integer())
  hmsc_next_group <- 0L
  hmsc_study_choices <- reactive({
    dat <- tryCatch(hmsc_read_design_csv(input$hmsc_study_file$datapath), error = function(e) NULL)
    c("Choose a study-design column" = "", names(dat))
  })
  hmsc_group_selection <- reactive({
    vapply(hmsc_group_rows(), function(i) input[[paste0("hmsc_group_col_", i)]] %||% "", character(1))
  })
  observeEvent(input$hmsc_add_group, {
    hmsc_next_group <<- hmsc_next_group + 1L
    i <- hmsc_next_group
    row_id <- paste0("hmsc-group-row-", i)
    column_id <- paste0("hmsc_group_col_", i)
    remove_id <- paste0("hmsc_group_remove_", i)
    insertUI("#hmsc_group_levels", where = "beforeEnd",
      ui = div(id = row_id, class = "studio-level-item",
        selectInput(column_id, "Study-design column", choices = hmsc_study_choices(), selected = ""),
        help_text("Rows sharing this label share a random effect. Example: plot_07 within the plot column."),
        actionButton(remove_id, label = NULL, class = "studio-level-remove",
                     title = "Remove grouping level", `aria-label` = "Remove grouping level")))
    hmsc_group_rows(c(hmsc_group_rows(), i))
    observeEvent(input[[remove_id]], {
      hmsc_group_rows(setdiff(hmsc_group_rows(), i))
      removeUI(paste0("#", row_id))
    }, once = TRUE, ignoreInit = TRUE)
  })
  observeEvent(hmsc_study_choices(), {
    choices <- hmsc_study_choices()
    for (i in hmsc_group_rows()) {
      id <- paste0("hmsc_group_col_", i)
      value <- input[[id]] %||% ""
      updateSelectInput(session, id, choices = choices, selected = if (value %in% choices) value else "")
    }
  }, ignoreInit = TRUE)

  refresh_hmsc_check <- function() {
    rv$hmsc$Y <- read_csv_safe(input$hmsc_Y_file$datapath)
    rv$hmsc$X <- read_csv_safe(input$hmsc_X_file$datapath)
    rv$hmsc$Tr <- read_csv_safe(input$hmsc_Tr_file$datapath)
    rv$hmsc$newdata <- read_csv_safe(input$hmsc_newdata_file$datapath)
    rv$hmsc$check <- tryCatch({
      rv$hmsc$study <- hmsc_read_design_csv(input$hmsc_study_file$datapath)
      rv$hmsc$coord <- hmsc_read_design_csv(input$hmsc_coord_file$datapath, c("", "row_id", "site_id", "sample_id", "id"))
      rv$hmsc$newstudy <- hmsc_read_design_csv(input$hmsc_newstudy_file$datapath)
      rv$hmsc$newcoord <- hmsc_read_design_csv(input$hmsc_newcoords_file$datapath, c("", "row_id", "site_id", "sample_id", "id"))
      validate_hmsc(rv$hmsc$Y, rv$hmsc$X, rv$hmsc$Tr, rv$hmsc$study, rv$hmsc$coord,
      rv$hmsc$newdata, rv$hmsc$newstudy, rv$hmsc$newcoord,
      input$hmsc_distr, input$hmsc_XFormula, input$hmsc_TrFormula, input$hmsc_use_traits, input$hmsc_use_phylogeny,
      input$hmsc_random_mode, input$hmsc_random_effect_column, input$hmsc_spatial_method, input$hmsc_nNeighbours,
      input$hmsc_lon_col, input$hmsc_lat_col, input$hmsc_samples, input$hmsc_transient, input$hmsc_thin,
      input$hmsc_nChains, input$hmsc_nParallel, input$hmsc_nfMin, input$hmsc_nfMax,
      input$hmsc_nfolds, input$hmsc_partition_column,
      additional_random_levels = hmsc_group_selection())
    }, error = function(e) list(ok = FALSE, messages = conditionMessage(e)))
    rv$hmsc$checked <- TRUE
  }

  observeEvent(input$hmsc_check, {
    refresh_hmsc_check()
    add_log("hmsc", "Hmsc data check completed.")
  })

  output$hmsc_check_messages <- renderText({
    if (!isTRUE(rv$hmsc$checked)) return("Hmsc data have not been checked yet.")
    paste(rv$hmsc$check$messages, collapse="\n")
  })

  output$hmsc_data_table <- renderDT({
    extra <- list(
      data.frame(File="newdata.csv", Required="Optional", Role="Prediction environmental covariates", Rows=if(is.null(rv$hmsc$newdata)) NA else nrow(rv$hmsc$newdata), Columns=if(is.null(rv$hmsc$newdata)) NA else ncol(rv$hmsc$newdata), Status=if(is.null(rv$hmsc$newdata)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="prediction_studyDesign.csv", Required="Conditional", Role="Prediction random-effect design", Rows=if(is.null(rv$hmsc$newstudy)) NA else nrow(rv$hmsc$newstudy), Columns=if(is.null(rv$hmsc$newstudy)) NA else ncol(rv$hmsc$newstudy), Status=if(is.null(rv$hmsc$newstudy)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="newcoords.csv", Required="Conditional", Role="Prediction spatial coordinates", Rows=if(is.null(rv$hmsc$newcoord)) NA else nrow(rv$hmsc$newcoord), Columns=if(is.null(rv$hmsc$newcoord)) NA else ncol(rv$hmsc$newcoord), Status=if(is.null(rv$hmsc$newcoord)) "Not uploaded or unreadable" else "Readable")
    )
    datatable(data_summary(rv$hmsc$Y, rv$hmsc$X, rv$hmsc$Tr, rv$hmsc$study, rv$hmsc$coord, extra),
              options = list(dom="tip", pageLength=10), rownames=FALSE)
  })

  # Hmsc-HPC check
  refresh_hmschpc_check <- function() {
    rv$hmschpc$Y <- read_csv_safe(input$hmschpc_Y_file$datapath)
    rv$hmschpc$X <- read_csv_safe(input$hmschpc_X_file$datapath)
    rv$hmschpc$Tr <- read_csv_safe(input$hmschpc_traits_file$datapath)
    rv$hmschpc$study <- read_csv_safe(input$hmschpc_study_file$datapath, id_columns = c("", "row.names", "row_id", "site_id", "sample_id"))
    rv$hmschpc$coord <- read_csv_safe(input$hmschpc_coord_file$datapath)
    rv$hmschpc$taxonomy <- read_csv_safe(input$hmschpc_taxonomy_file$datapath)
    rv$hmschpc$newdata <- read_csv_safe(input$hmschpc_newdata_file$datapath)
    rv$hmschpc$phylo_tree <- if (!is.null(input$hmschpc_phylo_tree_file$datapath) && file.exists(input$hmschpc_phylo_tree_file$datapath)) input$hmschpc_phylo_tree_file$datapath else NULL
    rv$hmschpc$check <- validate_hmschpc_full(
      rv$hmschpc$Y, rv$hmschpc$X, rv$hmschpc$Tr, rv$hmschpc$study, rv$hmschpc$coord,
      phylo_tree_file = rv$hmschpc$phylo_tree,
      taxonomy = rv$hmschpc$taxonomy,
      newdata = rv$hmschpc$newdata,
      distribution = input$hmschpc_distribution,
      XFormula = input$hmschpc_XFormula,
      use_traits = input$hmschpc_use_traits,
      trait_formula = input$hmschpc_trait_formula,
      phylogeny_mode = input$hmschpc_phylogeny_mode,
      random_mode = input$hmschpc_random_mode,
      random_column = input$hmschpc_random_column,
      coord_x = input$hmschpc_coord_x,
      coord_y = input$hmschpc_coord_y,
      nf = input$hmschpc_nf,
      nfMin = input$hmschpc_nfMin,
      nfMax = input$hmschpc_nfMax,
      samples = input$hmschpc_samples,
      transient = input$hmschpc_transient,
      thin = input$hmschpc_thin,
      chains = input$hmschpc_chains,
      verbose = input$hmschpc_verbose,
      run_sampler = input$hmschpc_run_sampler,
      random_slope_formula = input$hmschpc_random_slope_formula,
      random_name = input$hmschpc_random_name,
      alpha = input$hmschpc_alpha,
      chains_to_run = hmschpc_parse_chain_ids(input$hmschpc_chains_to_run, input$hmschpc_chains),
      seed = input$hmschpc_seed,
      precision = as.integer(input$hmschpc_precision),
      truncated_normal_library = input$hmschpc_tnlib,
      hmcleapfrog = input$hmschpc_hmcleapfrog,
      hmcthin = input$hmschpc_hmcthin,
      python = input$hmschpc_python,
      python_source = input$hmschpc_python_source,
      predictions = input$hmschpc_out_predictions,
      diagnostics = input$hmschpc_out_diagnostics,
      plots = input$hmschpc_out_plots
    )
    rv$hmschpc$checked <- TRUE
  }

  observeEvent(input$hmschpc_check, {
    refresh_hmschpc_check()
    add_log("hmschpc", "Hmsc-HPC data check completed.")
  })

  output$hmschpc_check_messages <- renderText({
    if (!isTRUE(rv$hmschpc$checked)) return("Hmsc-HPC data have not been checked yet.")
    paste(rv$hmschpc$check$messages, collapse="\n")
  })

  output$hmschpc_data_table <- renderDT({
    extra <- list(
      data.frame(File="taxonomy.csv", Required="Optional", Role="Taxonomic rank table used to derive phylogenetic correlation", Rows=if(is.null(rv$hmschpc$taxonomy)) NA else nrow(rv$hmschpc$taxonomy), Columns=if(is.null(rv$hmschpc$taxonomy)) NA else ncol(rv$hmschpc$taxonomy), Status=if(is.null(rv$hmschpc$taxonomy)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="phylo_tree.nwk", Required="Optional", Role="Newick phylogeny", Rows=NA, Columns=NA, Status=if(is.null(rv$hmschpc$phylo_tree)) "Not uploaded or unreadable" else "Readable text file"),
      data.frame(File="newdata.csv", Required="Optional", Role="Prediction covariates", Rows=if(is.null(rv$hmschpc$newdata)) NA else nrow(rv$hmschpc$newdata), Columns=if(is.null(rv$hmschpc$newdata)) NA else ncol(rv$hmschpc$newdata), Status=if(is.null(rv$hmschpc$newdata)) "Not uploaded or unreadable" else "Readable")
    )
    datatable(data_summary(rv$hmschpc$Y, rv$hmschpc$X, rv$hmschpc$Tr, rv$hmschpc$study, rv$hmschpc$coord, extra),
              options = list(dom="tip", pageLength=10), rownames=FALSE)
  })

  # jSDM check
  refresh_jsdm_check <- function() {
    rv$jsdm$Y <- read_csv_safe(input$jsdm_Y_file$datapath)
    rv$jsdm$X <- read_csv_safe(input$jsdm_X_file$datapath)
    rv$jsdm$Tr <- read_csv_safe(input$jsdm_trait_file$datapath)
    rv$jsdm$long <- read_csv_safe(input$jsdm_long_file$datapath, id_columns = c("", "row.names", "...1"),
                                  text_columns = c(input$jsdm_long_site_col, input$jsdm_long_species_col))
    rv$jsdm$trials <- read_csv_safe(input$jsdm_trials_file$datapath)
    rv$jsdm$newdata <- read_csv_safe(input$jsdm_newdata_file$datapath)
    rv$jsdm$prediction_ids <- read_csv_safe(input$jsdm_prediction_ids_file$datapath)
    rv$jsdm$check <- validate_jsdm_full(
      rv$jsdm$Y, rv$jsdm$X, rv$jsdm$Tr, rv$jsdm$long, rv$jsdm$trials,
      rv$jsdm$newdata, rv$jsdm$prediction_ids,
      input$jsdm_model_type, input$jsdm_site_formula, input$jsdm_trait_formula,
      input$jsdm_n_latent, input$jsdm_site_effect,
      input$jsdm_burnin, input$jsdm_mcmc, input$jsdm_thin, input$jsdm_trials,
      input$jsdm_allow_traits, input$jsdm_do_predict,
      input$jsdm_long_site_col, input$jsdm_long_species_col, input$jsdm_long_response_col,
      input$jsdm_constrained_nchains
    )
    rv$jsdm$checked <- TRUE
  }

  observeEvent(input$jsdm_check, {
    refresh_jsdm_check()
    add_log("jsdm", "jSDM data check completed.")
  })

  output$jsdm_check_messages <- renderText({
    if (!isTRUE(rv$jsdm$checked)) return("jSDM data have not been checked yet.")
    paste(rv$jsdm$check$messages, collapse="\n")
  })

  output$jsdm_data_table <- renderDT({
    extra <- list(
      data.frame(File="long_format.csv", Required="Only for long format", Role="Long-format observations", Rows=if(is.null(rv$jsdm$long)) NA else nrow(rv$jsdm$long), Columns=if(is.null(rv$jsdm$long)) NA else ncol(rv$jsdm$long), Status=if(is.null(rv$jsdm$long)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="trials.csv", Required="Optional", Role="Binomial trial counts", Rows=if(is.null(rv$jsdm$trials)) NA else nrow(rv$jsdm$trials), Columns=if(is.null(rv$jsdm$trials)) NA else ncol(rv$jsdm$trials), Status=if(is.null(rv$jsdm$trials)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="newdata.csv", Required="Optional", Role="Prediction covariates", Rows=if(is.null(rv$jsdm$newdata)) NA else nrow(rv$jsdm$newdata), Columns=if(is.null(rv$jsdm$newdata)) NA else ncol(rv$jsdm$newdata), Status=if(is.null(rv$jsdm$newdata)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="prediction_ids.csv", Required="Optional", Role="Prediction site/species IDs", Rows=if(is.null(rv$jsdm$prediction_ids)) NA else nrow(rv$jsdm$prediction_ids), Columns=if(is.null(rv$jsdm$prediction_ids)) NA else ncol(rv$jsdm$prediction_ids), Status=if(is.null(rv$jsdm$prediction_ids)) "Not uploaded or unreadable" else "Readable")
    )
    datatable(data_summary(rv$jsdm$Y, rv$jsdm$X, rv$jsdm$Tr, NULL, NULL, extra),
              options = list(dom="tip", pageLength=10), rownames=FALSE)
  })

  hmsc_config <- reactive({
    list(
      project_name = input$project_name,
      project_question = input$project_question,
      question_template = input$question_template,
      project_design = project_profile(),
      engine = "Hmsc",
      data = list(
        Y = if(!is.null(input$hmsc_Y_file$name)) input$hmsc_Y_file$name else NULL,
        XData = if(!is.null(input$hmsc_X_file$name)) input$hmsc_X_file$name else NULL,
        newdata = if(!is.null(input$hmsc_newdata_file$name)) input$hmsc_newdata_file$name else NULL,
        TrData = if(!is.null(input$hmsc_Tr_file$name)) input$hmsc_Tr_file$name else NULL,
        studyDesign = if(!is.null(input$hmsc_study_file$name)) input$hmsc_study_file$name else NULL,
        prediction_studyDesign = if(!is.null(input$hmsc_newstudy_file$name)) input$hmsc_newstudy_file$name else NULL,
        coordinates = if(!is.null(input$hmsc_coord_file$name)) input$hmsc_coord_file$name else NULL,
        newcoords = if(!is.null(input$hmsc_newcoords_file$name)) input$hmsc_newcoords_file$name else NULL,
        phylogeny = if(!is.null(input$hmsc_phylo_file$name)) input$hmsc_phylo_file$name else NULL
      ),
      model = list(
        distr = input$hmsc_distr,
        XFormula = input$hmsc_XFormula,
        TrFormula = input$hmsc_TrFormula,
        use_traits = input$hmsc_use_traits,
        use_phylogeny = input$hmsc_use_phylogeny,
        random_mode = input$hmsc_random_mode,
        random_effect_column = input$hmsc_random_effect_column,
        additional_random_levels = hmsc_group_selection(),
        spatial_method = input$hmsc_spatial_method,
        nNeighbours = input$hmsc_nNeighbours,
        lon_col = input$hmsc_lon_col,
        lat_col = input$hmsc_lat_col,
        seed = input$hmsc_seed,
        XScale = input$hmsc_XScale,
        TrScale = input$hmsc_TrScale,
        YScale = input$hmsc_YScale,
        truncateNumberOfFactors = input$hmsc_truncateNumberOfFactors,
        Loff_file = input$hmsc_Loff_file,
        ranLevelsUsed = input$hmsc_ranLevelsUsed,
        C_file = "",
        use_XRRR = input$hmsc_use_XRRR,
        ncRRR = input$hmsc_ncRRR,
        XRRRFormula = input$hmsc_XRRRFormula,
        XRRRScale = input$hmsc_XRRRScale,
        XRRR_file = input$hmsc_XRRR_file,
        random_level_type = input$hmsc_random_level_type,
        sMethod = input$hmsc_sMethod,
        random_N = input$hmsc_random_N,
        longlat = input$hmsc_longlat,
        units_column = input$hmsc_units_column,
        distMat_file = input$hmsc_distMat_file,
        xData_file = input$hmsc_xData_file,
        sKnot_file = input$hmsc_sKnot_file,
        nfMin = input$hmsc_nfMin,
        nfMax = input$hmsc_nfMax,
        priors = list(
          setDefault = input$hmsc_use_default_priors,
          a1 = input$hmsc_a1,
          b1 = input$hmsc_b1,
          a2 = input$hmsc_a2,
          b2 = input$hmsc_b2,
          alphapw = input$hmsc_alphapw
        )
      ),
      mcmc = list(
        samples = input$hmsc_samples,
        transient = input$hmsc_transient,
        thin = input$hmsc_thin,
        nChains = input$hmsc_nChains,
        nParallel = input$hmsc_nParallel,
        verbose = input$hmsc_verbose,
        preset = input$hmsc_preset,
        initPar = input$hmsc_initPar,
        alignPost = input$hmsc_alignPost,
        updater = list(
          GammaEta = input$hmsc_updater_GammaEta,
          Beta = input$hmsc_updater_Beta,
          Gamma = input$hmsc_updater_Gamma,
          Omega = input$hmsc_updater_Omega
        ),
        sample_prior = input$hmsc_sample_prior,
        pool_chains = input$hmsc_pool_chains
      ),
      outputs = list(
        save_model = input$hmsc_save_model,
        predicted = input$hmsc_out_predicted,
        fit = input$hmsc_out_fit,
        cv = input$hmsc_out_cv,
        waic = input$hmsc_out_waic,
        diagnostics = input$hmsc_out_diag,
        parameters = input$hmsc_out_params,
        variance_partitioning = input$hmsc_out_vp,
        omega = input$hmsc_out_omega,
        gradients = input$hmsc_out_gradients,
        beta_support = input$hmsc_beta_support,
        gamma_support = input$hmsc_gamma_support,
        omega_support = input$hmsc_omega_support,
        convergence = list(
          showBeta = input$hmsc_conv_beta,
          showGamma = input$hmsc_conv_gamma,
          showOmega = input$hmsc_conv_omega,
          maxOmega = input$hmsc_maxOmega,
          showRho = input$hmsc_conv_rho,
          showAlpha = input$hmsc_conv_alpha,
          effectiveSize = input$hmsc_effective_size,
          gelmanPSRF = input$hmsc_gelman_psrf
        ),
        plotting = list(
          var.part.order.explained = input$hmsc_varpart_order_explained,
          var.part.order.raw = input$hmsc_varpart_order_raw,
          show.sp.names.beta = input$hmsc_show_sp_names_beta,
          plotTree = input$hmsc_plotTree,
          omega.order = input$hmsc_omega_order,
          show.sp.names.omega = input$hmsc_show_sp_names_omega,
          plotBeta = input$hmsc_plot_beta,
          plotGamma = input$hmsc_plot_gamma
        ),
        predictions = list(
          species.list = input$hmsc_species_list,
          trait.list = input$hmsc_trait_list,
          env.list = input$hmsc_env_list,
          nfolds = input$hmsc_nfolds,
          partition_column = input$hmsc_partition_column,
          computeSAIR = input$hmsc_compute_sair
        )
      )
    )
  })

  hmschpc_config <- reactive({
    chain_ids <- trimws(input$hmschpc_chains_to_run %||% "")
    chain_ids <- hmschpc_parse_chain_ids(chain_ids, input$hmschpc_chains)
    list(
      project_name = input$project_name,
      project_question = input$project_question,
      question_template = input$question_template,
      project_design = project_profile(),
      engine = "Hmsc-HPC",
      data = list(
        Y = if(!is.null(input$hmschpc_Y_file$name)) input$hmschpc_Y_file$name else NULL,
        XData = if(!is.null(input$hmschpc_X_file$name)) input$hmschpc_X_file$name else NULL,
        traits = if(!is.null(input$hmschpc_traits_file$name)) input$hmschpc_traits_file$name else NULL,
        studyDesign = if(!is.null(input$hmschpc_study_file$name)) input$hmschpc_study_file$name else NULL,
        coordinates = if(!is.null(input$hmschpc_coord_file$name)) input$hmschpc_coord_file$name else NULL,
        taxonomy = if(!is.null(input$hmschpc_taxonomy_file$name)) input$hmschpc_taxonomy_file$name else NULL,
        phylo_tree = if(!is.null(input$hmschpc_phylo_tree_file$name)) input$hmschpc_phylo_tree_file$name else NULL,
        newdata = if(!is.null(input$hmschpc_newdata_file$name)) input$hmschpc_newdata_file$name else NULL
      ),
        model = list(
        distribution = hmschpc_distribution_key(input$hmschpc_distribution),
        XFormula = input$hmschpc_XFormula,
        use_traits = input$hmschpc_use_traits,
        trait_formula = input$hmschpc_trait_formula,
        phylogeny_mode = input$hmschpc_phylogeny_mode
      ),
      random_effects = list(
        mode = input$hmschpc_random_mode,
        name = input$hmschpc_random_name,
        column = input$hmschpc_random_column,
        coord_x = input$hmschpc_coord_x,
        coord_y = input$hmschpc_coord_y,
        x_formula = input$hmschpc_random_slope_formula,
        nf = input$hmschpc_nf,
        nfMin = input$hmschpc_nfMin,
        nfMax = input$hmschpc_nfMax,
        alpha = input$hmschpc_alpha
      ),
      sampler = list(
        run_sampler = input$hmschpc_run_sampler,
        samples = input$hmschpc_samples,
        transient = input$hmschpc_transient,
        thin = input$hmschpc_thin,
        chains = input$hmschpc_chains,
        chains_to_run = as.list(chain_ids),
        verbose = input$hmschpc_verbose,
        seed = input$hmschpc_seed,
        precision = as.integer(input$hmschpc_precision),
        truncated_normal_library = input$hmschpc_tnlib,
        hmcleapfrog = input$hmschpc_hmcleapfrog,
        hmcthin = input$hmschpc_hmcthin,
        update_beta_eta = input$hmschpc_update_beta_eta,
        save_eta = input$hmschpc_save_eta,
        eager = input$hmschpc_eager,
        profile = input$hmschpc_profile
      ),
      outputs = list(
        predictions = input$hmschpc_out_predictions,
        diagnostics = input$hmschpc_out_diagnostics,
        plots = input$hmschpc_out_plots,
        zip = input$hmschpc_out_zip
      ),
      runtime = list(
        python = input$hmschpc_python,
        python_source = input$hmschpc_python_source
      )
    )
  })

  jsdm_config <- reactive({
    list(
      project_name = input$project_name,
      project_question = input$project_question,
      question_template = input$question_template,
      project_design = project_profile(),
      engine = "jSDM",
      data = list(
        Y = if(!is.null(input$jsdm_Y_file$name)) input$jsdm_Y_file$name else NULL,
        XData = if(!is.null(input$jsdm_X_file$name)) input$jsdm_X_file$name else NULL,
        trait_data = if(!is.null(input$jsdm_trait_file$name)) input$jsdm_trait_file$name else NULL,
        long_format = if(!is.null(input$jsdm_long_file$name)) input$jsdm_long_file$name else NULL,
        trials = if(!is.null(input$jsdm_trials_file$name)) input$jsdm_trials_file$name else NULL,
        newdata = if(!is.null(input$jsdm_newdata_file$name)) input$jsdm_newdata_file$name else NULL,
        prediction_ids = if(!is.null(input$jsdm_prediction_ids_file$name)) input$jsdm_prediction_ids_file$name else NULL
      ),
      model = list(
        model_type = input$jsdm_model_type,
        response_argument = input$jsdm_response_argument,
        site_formula = input$jsdm_site_formula,
        trait_formula = input$jsdm_trait_formula,
        n_latent = input$jsdm_n_latent,
        site_effect = input$jsdm_site_effect,
        trials = input$jsdm_trials,
        constrained_latent = input$jsdm_constrained_latent,
        constrained_nchains = input$jsdm_constrained_nchains,
        long_site_col = input$jsdm_long_site_col,
        long_species_col = input$jsdm_long_species_col,
        long_response_col = input$jsdm_long_response_col,
        scale_site_data = input$jsdm_scale_site_data,
        include_intercept = input$jsdm_include_intercept,
        allow_traits = input$jsdm_allow_traits
      ),
      prediction = list(
        do_predict = input$jsdm_do_predict,
        predict_type = input$jsdm_predict_type,
        predict_probs = input$jsdm_predict_probs,
        max_prediction_sites = input$jsdm_predict_max_sites,
        Id_sites = input$jsdm_Id_sites,
        Id_species = input$jsdm_Id_species,
        prediction_histograms = input$jsdm_prediction_histograms
      ),
      diagnostics = list(
        cor_prob = input$jsdm_cor_prob,
        cor_type = input$jsdm_cor_type,
        plot_residual_cor = input$jsdm_plot_residual_cor,
        plot_associations = input$jsdm_plot_associations,
        diag_beta = input$jsdm_diag_beta,
        diag_lambda = input$jsdm_diag_lambda,
        diag_W = input$jsdm_diag_W,
        diag_alpha = input$jsdm_diag_alpha,
        diag_Valpha = input$jsdm_diag_Valpha,
        diag_V = input$jsdm_diag_V,
        diag_deviance = input$jsdm_diag_deviance,
        coda_summary = input$jsdm_coda_summary
      ),
      mcmc = list(
        burnin = input$jsdm_burnin,
        mcmc = input$jsdm_mcmc,
        thin = input$jsdm_thin,
        seed = input$jsdm_seed,
        verbose = input$jsdm_verbose,
        preset = input$jsdm_preset,
        ropt = input$jsdm_ropt
      ),
      starts = list(
        beta_start = input$jsdm_beta_start,
        gamma_start = input$jsdm_gamma_start,
        lambda_start = input$jsdm_lambda_start,
        W_start = input$jsdm_W_start,
        alpha_start = input$jsdm_alpha_start,
        V_alpha = input$jsdm_V_alpha,
        V_start = input$jsdm_V_start
      ),
      priors = list(
        shape_Valpha = input$jsdm_shape_Valpha,
        rate_Valpha = input$jsdm_rate_Valpha,
        shape_V = input$jsdm_shape_V,
        rate_V = input$jsdm_rate_V,
        mu_beta = input$jsdm_mu_beta,
        V_beta = input$jsdm_V_beta,
        mu_gamma = input$jsdm_mu_gamma,
        V_gamma = input$jsdm_V_gamma,
        mu_lambda = input$jsdm_mu_lambda,
        V_lambda = input$jsdm_V_lambda
      ),
      outputs = list(
        real_fit = input$jsdm_real_fit,
        residual_cor = input$jsdm_out_resid,
        enviro_cor = input$jsdm_out_env,
        traceplots = input$jsdm_out_trace,
        predictions = input$jsdm_out_pred,
        save_model = input$jsdm_out_model,
        report = input$jsdm_out_report,
        save_mcmc_rds = input$jsdm_out_mcmc_rds,
        csv_tables = input$jsdm_out_csv_tables,
        figures = input$jsdm_out_figures,
        table_model_spec = input$jsdm_table_model_spec,
        table_beta = input$jsdm_table_beta,
        table_lambda = input$jsdm_table_lambda,
        table_gamma = input$jsdm_table_gamma,
        table_alpha = input$jsdm_table_alpha,
        table_predictions = input$jsdm_table_predictions,
        copy_inputs = input$jsdm_out_copy_inputs,
        save_config = input$jsdm_out_config,
        zip = input$jsdm_out_zip
      )
    )
  })

  run_safe_engine <- function(engine) {
    is_hmsc <- engine == "Hmsc"
    state <- if (is_hmsc) rv$hmsc else rv$jsdm
    cfg <- if (is_hmsc) hmsc_config() else jsdm_config()
    pkg <- if (is_hmsc) "Hmsc" else "jSDM"
    outdir <- make_engine_run_dir(engine, input$project_name)
    if (is_hmsc) { rv$hmsc$outdir <- outdir; rv$hmsc$status <- "Running" } else { rv$jsdm$outdir <- outdir; rv$jsdm$status <- "Running" }
    add_log(if(is_hmsc) "hmsc" else "jsdm", "Created output folder:", outdir)

    # Copy inputs
    if (is_hmsc && isTRUE(input$hmsc_out_copy_inputs)) {
      copy_upload(input$hmsc_Y_file, outdir, "Y.csv")
      copy_upload(input$hmsc_X_file, outdir, "XData.csv")
      copy_upload(input$hmsc_newdata_file, outdir, "newdata.csv")
      copy_upload(input$hmsc_Tr_file, outdir, "TrData.csv")
      copy_upload(input$hmsc_study_file, outdir, "studyDesign.csv")
      copy_upload(input$hmsc_newstudy_file, outdir, "prediction_studyDesign.csv")
      copy_upload(input$hmsc_coord_file, outdir, "coordinates.csv")
      copy_upload(input$hmsc_newcoords_file, outdir, "newcoords.csv")
      if (!is.null(input$hmsc_phylo_file) && file.exists(input$hmsc_phylo_file$datapath)) {
        file.copy(input$hmsc_phylo_file$datapath, file.path(outdir, "inputs", input$hmsc_phylo_file$name), overwrite=TRUE)
        file.copy(input$hmsc_phylo_file$datapath, file.path(outdir, "data", input$hmsc_phylo_file$name), overwrite=TRUE)
      }
    }
    if (!is_hmsc && (isTRUE(input$jsdm_out_copy_inputs) || isTRUE(cfg$outputs$real_fit))) {
      copy_upload(input$jsdm_Y_file, outdir, "Y.csv")
      copy_upload(input$jsdm_X_file, outdir, "XData.csv")
      copy_upload(input$jsdm_trait_file, outdir, "trait_data.csv")
      copy_upload(input$jsdm_long_file, outdir, "long_format.csv")
      copy_upload(input$jsdm_trials_file, outdir, "trials.csv")
      copy_upload(input$jsdm_newdata_file, outdir, "newdata.csv")
      copy_upload(input$jsdm_prediction_ids_file, outdir, "prediction_ids.csv")
    }

    if (is_hmsc || isTRUE(input$jsdm_out_config) || isTRUE(cfg$outputs$real_fit)) {
      yaml::write_yaml(cfg, file.path(outdir, "used_config.yml"))
    }

    # data summary and checks
    if (is_hmsc) {
      write.csv(data_summary(rv$hmsc$Y, rv$hmsc$X, rv$hmsc$Tr, rv$hmsc$study, rv$hmsc$coord),
                file.path(outdir, "tables", "data_summary.csv"), row.names=FALSE)
      write_data_check_messages(outdir, if (!is.null(rv$hmsc$check)) rv$hmsc$check$messages else character())
    } else {
      write.csv(data_summary(rv$jsdm$Y, rv$jsdm$X, rv$jsdm$Tr),
                file.path(outdir, "tables", "data_summary.csv"), row.names=FALSE)
      write_data_check_messages(outdir, if (!is.null(rv$jsdm$check)) rv$jsdm$check$messages else character())
    }

    # Safe engine status
    status <- list(engine=engine, status="ready", runtime_seconds=0, warnings=character(), errors=character())
    if (engine == "jSDM") {
      if (((cfg$mcmc$burnin + cfg$mcmc$mcmc) %% 10) != 0 || (cfg$mcmc$burnin + cfg$mcmc$mcmc) < 100) {
        status$warnings <- c(status$warnings, "jSDM warning: burnin + mcmc should be divisible by 10 and at least 100 for the progress bar / sampler requirements.")
      }
      if ((cfg$mcmc$mcmc %% cfg$mcmc$thin) != 0) {
        status$warnings <- c(status$warnings, "jSDM warning: mcmc must be divisible by thin.")
      }
      if (isTRUE(cfg$model$constrained_latent) && cfg$model$n_latent <= 0) {
        status$warnings <- c(status$warnings, "jSDM warning: constrained probit model requires n_latent > 0.")
      }
      if (cfg$model$model_type == "binomial_logit" && cfg$model$trials < 1) {
        status$warnings <- c(status$warnings, "jSDM warning: trials should be >= 1 for binomial logit.")
      }
    }
    real_jsdm_requested <- identical(engine, "jSDM") && isTRUE(cfg$outputs$real_fit %||% FALSE)
    if (!safe_require(pkg)) {
      status$status <- "check_failed"
      status$warnings <- paste0(pkg, " package is not available. A diagnostic output scaffold was saved; no model fit was attempted.")
    } else if (is_hmsc && isTRUE(input$hmsc_real_fit %||% FALSE)) {
      status$status <- "model_defined"
      status$warnings <- character()
    } else if (real_jsdm_requested) {
      status$status <- "model_defined"
      status$warnings <- character()
    } else {
      status$status <- "model_defined"
      status$warnings <- paste0(pkg, " package appears available, but this workflow is currently configured as scaffold-only unless a production adapter is connected.")
    }

    status <- normalize_engine_status(status)
    write_engine_status(outdir, status)

    if (engine == "Hmsc") {
      # Create a transparent Hmsc workflow map so users know what the downloaded result ZIP contains.
      workflow_files <- data.frame(
        Step = c("S1 Define models", "S2 Fit models", "S3 Evaluate convergence", "S4 Compute model fit", "S5 Show model fit", "S6 Parameter estimates", "S7 Predictions"),
        Expected_file = c(
          "models/unfitted_models.RData",
          "models/models_thin_[thin]_samples_[samples]_chains_[chains].Rdata",
          "results/MCMC_convergence.pdf and results/MCMC_convergence.txt",
          "models/MF_thin_[thin]_samples_[samples]_chains_[chains]_nfolds_[nfolds].Rdata",
          "results/model_fit.pdf",
          "results/parameter_estimates.pdf, parameter_estimates.txt, parameter_estimates_[parameter].csv",
          "results/predictions.pdf"
        ),
        Meaning = c(
          "Model object(s) constructed but not fitted.",
          "Fitted Hmsc model list after MCMC sampling.",
          "MCMC diagnostics for selected Beta/Gamma/Omega/rho/alpha parameters.",
          "Model fit / cross-validation object.",
          "Visual summary of model explanatory and predictive fit.",
          "Parameter support, plots and CSV tables for Beta/Gamma/Omega/variance partitioning.",
          "Predictions along environmental gradients for selected species/traits."
        )
      )
      write.csv(workflow_files, file.path(outdir, "tables", "Hmsc_result_workflow_map.csv"), row.names = FALSE)
      writeLines(c(
        "Hmsc result workflow",
        "====================",
        "This SAFE build creates the output scaffold and configuration.",
        "When the production Hmsc adapter is connected, it should write the following standard files:",
        "",
        paste(workflow_files$Step, workflow_files$Expected_file, sep = " -> ")
      ), file.path(outdir, "results", "README_Hmsc_results.txt"))

      # Reproducible script export required for GUI-based reproducibility.
      if (isTRUE(input$hmsc_export_scripts %||% TRUE)) {
        script_lines <- c(
          "# Reproducible HMSC script generated by JSDMWorkbench",
          "library(Hmsc)",
          "library(yaml)",
          "cfg <- yaml::read_yaml('../used_config.yml')",
          "Y <- read.csv('../inputs/Y.csv', row.names = 1, check.names = FALSE)",
          "XData <- read.csv('../inputs/XData.csv', row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)",
          "for (nm in names(XData)) if (is.character(XData[[nm]])) { x <- trimws(XData[[nm]]); x[x == ''] <- NA; nx <- suppressWarnings(as.numeric(x)); XData[[nm]] <- if (all(is.na(x) | !is.na(nx))) nx else as.factor(x) }",
          "XFormula <- as.formula(cfg$model$XFormula)",
          "m <- Hmsc(Y = as.matrix(Y), XData = XData, XFormula = XFormula, distr = cfg$model$distr)",
          "m <- sampleMcmc(m, thin = cfg$mcmc$thin, samples = cfg$mcmc$samples, transient = cfg$mcmc$transient, nChains = cfg$mcmc$nChains, nParallel = cfg$mcmc$nParallel, verbose = cfg$mcmc$verbose)",
          "saveRDS(m, '../models/hmsc_model.rds')",
          "preds <- computePredictedValues(m)",
          "fit <- evaluateModelFit(hM = m, predY = preds)",
          "saveRDS(preds, '../predictions/hmsc_predicted_values.rds')",
          "saveRDS(fit, '../results/hmsc_model_fit.rds')"
        )
        writeLines(script_lines, file.path(outdir, "reproducible_script", "run_this_HMSC_analysis.R"))
        for (i in seq_len(nrow(workflow_files))) {
          writeLines(c(
            paste0("# ", workflow_files$Step[i]),
            paste0("# Expected output: ", workflow_files$Expected_file[i]),
            paste0("# Meaning: ", workflow_files$Meaning[i]),
            "# This script is a transparent placeholder in the SAFE build.",
            "# Production fitting can replace this file with executable HMSC code."
          ), file.path(outdir, "workflow_scripts", paste0("S", i, "_", gsub("[^A-Za-z0-9]+", "_", workflow_files$Step[i]), ".R")))
        }
      }

      # Save the original uploaded S1-S7 scripts for traceability and a master runner.
      original_dir <- file.path("workflow_templates", "HMSC_user_S1S7_original")
      if (dir.exists(original_dir)) {
        for (ff in list.files(original_dir, pattern="\\.R$", full.names=TRUE)) {
          file.copy(ff, file.path(outdir, "workflow_scripts", basename(ff)), overwrite=TRUE)
        }
      }
      writeLines(c(
        "# Master HMSC S1-S7 runner generated by JSDMWorkbench",
        "# This script is designed to run inside the downloaded output folder.",
        "source('reproducible_script/run_this_HMSC_analysis.R')",
        "# The GUI also saves original user-provided S1-S7 scripts in workflow_scripts/ for comparison."
      ), file.path(outdir, "reproducible_script", "run_all_HMSC_S1_to_S7.R"))

      if (isTRUE(input$hmsc_create_standard %||% TRUE)) {
        write.csv(data.frame(
          run_id = basename(outdir),
          engine = "Hmsc",
          status = status$status,
          n_sites = if (is.null(rv$hmsc$Y)) NA else nrow(rv$hmsc$Y),
          n_responses = if (is.null(rv$hmsc$Y)) NA else ncol(rv$hmsc$Y),
          n_predictors = if (is.null(rv$hmsc$X)) NA else ncol(rv$hmsc$X),
          distr = cfg$model$distr,
          random_mode = cfg$model$random_mode,
          has_traits = isTRUE(cfg$model$use_traits),
          has_phylogeny = isTRUE(cfg$model$use_phylogeny),
          stringsAsFactors = FALSE
        ), file.path(outdir, "standard", "run_summary.csv"), row.names = FALSE)
        write.csv(data.frame(engine=character(), response_id=character(), predictor=character(), direction=character(), estimate=numeric(), lower=numeric(), upper=numeric(), notes=character()), file.path(outdir, "standard", "effects_long.csv"), row.names = FALSE)
        write.csv(data.frame(engine=character(), site_id=character(), response_id=character(), observed=numeric(), predicted_mean=numeric(), predicted_lower=numeric(), predicted_upper=numeric()), file.path(outdir, "standard", "predictions_long.csv"), row.names = FALSE)
        write.csv(data.frame(engine=character(), response_1=character(), response_2=character(), association_type=character(), estimate=numeric(), comparable_level=character()), file.path(outdir, "standard", "associations_long.csv"), row.names = FALSE)
        write.csv(data.frame(file=list.files(outdir, recursive=TRUE), stringsAsFactors=FALSE), file.path(outdir, "standard", "output_manifest.csv"), row.names = FALSE)
      }
      if (isTRUE(input$hmsc_real_fit %||% FALSE)) {
        add_log("hmsc", "Starting real HMSC S1-S7 results workflow.")
        fit_result <- run_hmsc_s1s7_pipeline(outdir, cfg, rv$hmsc$Y, rv$hmsc$X,
                                              TrData = if (isTRUE(input$hmsc_use_traits)) rv$hmsc$Tr else NULL,
                                              studyDesign = rv$hmsc$study,
                                              coord = rv$hmsc$coord,
                                              newdata = rv$hmsc$newdata,
                                              newStudyDesign = rv$hmsc$newstudy,
                                              newcoord = rv$hmsc$newcoord,
                                              log_fun = function(txt) add_log("hmsc", txt))
        status$status <- fit_result$status
        status$warnings <- c(status$warnings, fit_result$warnings)
        status$errors <- c(status$errors, fit_result$errors)
        status <- normalize_engine_status(status)
        write_engine_status(outdir, status)
        add_log("hmsc", "HMSC S1-S7 results workflow status:", fit_result$status)
      }
    }

    if (engine == "jSDM") {
      jsdm_files <- data.frame(
        Section = c("Configuration", "Inputs", "Model object", "MCMC objects", "Parameter tables", "Correlation tables", "Predictions", "Plots", "Diagnostics", "Report"),
        Expected_file_or_folder = c(
          "used_config.yml",
          "inputs/Y.csv, inputs/XData.csv, inputs/trait_data.csv, inputs/trials.csv, inputs/newdata.csv",
          "models/jsdm_model.rds",
          "mcmc/mcmc_sp.rds, mcmc/mcmc_gamma.rds, mcmc/mcmc_latent.rds, mcmc/mcmc_alpha.rds, mcmc/mcmc_V_alpha.rds, mcmc/mcmc_V.rds, mcmc/mcmc_Deviance.rds",
          "tables/beta_summary.csv, lambda_summary.csv, gamma_summary.csv, alpha_summary.csv, V_alpha_summary.csv, V_summary.csv, Deviance_summary.csv",
          "tables/residual_cor_*.csv, tables/enviro_cor_*.csv",
          "tables/predictions.csv or results/predictions.rds",
          "plots/jSDM_traceplots.pdf, plots/jSDM_densityplots.pdf, plots/residual_cor.pdf, plots/enviro_cor.pdf, plots/prediction_histograms.pdf",
          "diagnostics/engine_status.json, diagnostics/data_check_messages.csv",
          "report/jSDM_report.html"
        ),
        Meaning = c(
          "Exact settings used for the jSDM run.",
          "Copied data used by the jSDM workflow.",
          "Serialized fitted jSDM model object when fitting is enabled.",
          "Raw coda-compatible posterior samples returned by jSDM.",
          "Posterior summaries for species effects, latent loadings, traits, site effects, variance and deviance.",
          "Correlation matrices derived from latent variables and shared environmental responses.",
          "Predicted posterior mean/quantiles/posterior from predict.jSDM depending on selected type.",
          "Visual diagnostics and result figures.",
          "Machine-readable status and warnings/errors.",
          "Human-readable browser report."
        )
      )
      write.csv(jsdm_files, file.path(outdir, "tables", "jSDM_result_workflow_map.csv"), row.names = FALSE)
      writeLines(c(
        "jSDM result workflow",
        "====================",
        "This workflow exports a standalone R script and, when enabled, runs the installed jSDM package.",
        "The reproducible script writes the standard files listed in tables/jSDM_result_workflow_map.csv.",
        "",
        paste(jsdm_files$Section, jsdm_files$Expected_file_or_folder, sep = " -> ")
      ), file.path(outdir, "results", "README_jSDM_results.txt"))

      if (exists("write_jsdm_reproducible_script", mode = "function")) {
        write_jsdm_reproducible_script(outdir)
      } else {
        status$status <- "fit_failed"
        status$errors <- c(status$errors, "Internal jSDM adapter function write_jsdm_reproducible_script() is missing.")
        writeLines(status$errors, file.path(outdir, "diagnostics", "jSDM_adapter_error.txt"))
      }

      if (isTRUE(cfg$outputs$real_fit %||% FALSE)) {
        if (!safe_require("jSDM")) {
          status$status <- "fit_failed"
          status$errors <- c(status$errors, "Package jSDM is not installed, so real fitting cannot run.")
          writeLines(status$errors, file.path(outdir, "diagnostics", "jSDM_reproducible_error.txt"))
          write_standard_outputs(outdir, "jSDM", status, rv$jsdm$Y, rv$jsdm$X)
          write_engine_status(outdir, status)
        } else if (file.exists(file.path(outdir, "reproducible_script", "run_this_jSDM_analysis.R"))) {
          add_log("jsdm", "Starting real jSDM package workflow.")
          script_file <- file.path(outdir, "reproducible_script", "run_this_jSDM_analysis.R")
          log_file <- file.path(outdir, "diagnostics", "jSDM_real_fit_stdout_stderr.txt")
          rscript_bin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
          started <- Sys.time()
          exit_code <- tryCatch({
            neutralize_windows_r_locale_env()
            system2(rscript_bin, shQuote(script_file), stdout = log_file, stderr = log_file)
          },
            error = function(e) {
              writeLines(e$message, file.path(outdir, "diagnostics", "jSDM_system2_error.txt"))
              1L
            }
          )
          status$runtime_seconds <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2)
          engine_json <- file.path(outdir, "diagnostics", "engine_status.json")
          if (file.exists(engine_json) && requireNamespace("jsonlite", quietly = TRUE)) {
            script_status <- tryCatch(jsonlite::fromJSON(engine_json), error = function(e) NULL)
            if (!is.null(script_status)) {
              if (!is.null(script_status$status)) status$status <- as.character(script_status$status)
              if (!is.null(script_status$warnings)) status$warnings <- unique(c(status$warnings, as.character(script_status$warnings)))
              if (!is.null(script_status$errors)) status$errors <- unique(c(status$errors, as.character(script_status$errors)))
            }
          }
          if (!identical(as.integer(exit_code), 0L) && !identical(status$status, "fitted")) {
            status$status <- "fit_failed"
            status$errors <- unique(c(status$errors, paste0("Rscript exited with status ", as.integer(exit_code), ". See diagnostics/jSDM_real_fit_stdout_stderr.txt.")))
          }
          write_engine_status(outdir, status)
          add_log("jsdm", "Real jSDM package workflow status:", status$status)
        } else {
          status$status <- "fit_failed"
          status$errors <- c(status$errors, "The jSDM reproducible script was not created.")
          writeLines(status$errors, file.path(outdir, "diagnostics", "jSDM_reproducible_error.txt"))
          write_standard_outputs(outdir, "jSDM", status, rv$jsdm$Y, rv$jsdm$X)
          write_engine_status(outdir, status)
        }
      }
    }

    status <- normalize_engine_status(status)
    if (!(engine == "Hmsc" && identical(status$status, "fitted")) &&
        !(engine == "jSDM" && isTRUE(cfg$outputs$real_fit %||% FALSE))) {
      write_standard_outputs(outdir, engine, status,
                             Y = if (is_hmsc) rv$hmsc$Y else rv$jsdm$Y,
                             X = if (is_hmsc) rv$hmsc$X else rv$jsdm$X)
      if (!is_hmsc) write_reproducible_stub(outdir, engine)
      if (!is_hmsc) write_engine_scaffold_outputs(outdir, engine, cfg, status, rv$jsdm$Y, rv$jsdm$X)
    }
    ensure_output_contract(outdir, engine, status,
                           Y = if (is_hmsc) rv$hmsc$Y else rv$jsdm$Y,
                           X = if (is_hmsc) rv$hmsc$X else rv$jsdm$X)
    make_html_report(outdir, engine, cfg, paste(status$status, paste(status$warnings, collapse="; ")))
    if (identical(status$status, "fitted")) {
      writeLines("RUN COMPLETE", file.path(outdir, "RUN_COMPLETE.txt"))
    } else if (identical(status$status, "fit_failed")) {
      writeLines("RUN FAILED", file.path(outdir, "RUN_FAILED.txt"))
    } else if (identical(status$status, "model_defined")) {
      writeLines("MODEL DEFINED - NO FITTED POSTERIOR", file.path(outdir, "RUN_MODEL_DEFINED.txt"))
    } else if (identical(status$status, "check_failed")) {
      writeLines("CHECK FAILED - NO FIT WAS ATTEMPTED", file.path(outdir, "RUN_CHECK_FAILED.txt"))
    } else {
      writeLines("FIT FAILED - UNRECOGNIZED TERMINAL STATUS", file.path(outdir, "RUN_FAILED.txt"))
    }
    zipfile <- if (is_hmsc || isTRUE(input$jsdm_out_zip)) make_zip(outdir) else NULL

    if (is_hmsc) {
      rv$hmsc$zip <- zipfile
      rv$hmsc$status <- status$status
      add_log("hmsc", "Hmsc workflow status:", status$status, outdir)
    } else {
      rv$jsdm$zip <- zipfile
      rv$jsdm$status <- status$status
      add_log("jsdm", "jSDM workflow status:", status$status, outdir)
    }
  }


  create_hmsc_check_failed_zip <- function() {
    outdir <- make_engine_run_dir("Hmsc_check_failed", input$project_name)
    rv$hmsc$outdir <- outdir
    cfg <- tryCatch(hmsc_config(), error = function(e) list(engine="Hmsc", error=paste("Could not create config:", e$message)))
    tryCatch(yaml::write_yaml(cfg, file.path(outdir, "used_config.yml")), error=function(e) NULL)
    copy_upload(input$hmsc_Y_file, outdir, "Y.csv")
    copy_upload(input$hmsc_X_file, outdir, "XData.csv")
    copy_upload(input$hmsc_newdata_file, outdir, "newdata.csv")
    copy_upload(input$hmsc_Tr_file, outdir, "TrData.csv")
    copy_upload(input$hmsc_study_file, outdir, "studyDesign.csv")
    copy_upload(input$hmsc_newstudy_file, outdir, "prediction_studyDesign.csv")
    copy_upload(input$hmsc_coord_file, outdir, "coordinates.csv")
    copy_upload(input$hmsc_newcoords_file, outdir, "newcoords.csv")
    if (!is.null(input$hmsc_phylo_file) && file.exists(input$hmsc_phylo_file$datapath)) {
      file.copy(input$hmsc_phylo_file$datapath, file.path(outdir, "inputs", input$hmsc_phylo_file$name), overwrite=TRUE)
      file.copy(input$hmsc_phylo_file$datapath, file.path(outdir, "data", input$hmsc_phylo_file$name), overwrite=TRUE)
    }
    extra <- list(
      data.frame(File="newdata.csv", Required="Optional", Role="Prediction environmental covariates", Rows=if(is.null(rv$hmsc$newdata)) NA else nrow(rv$hmsc$newdata), Columns=if(is.null(rv$hmsc$newdata)) NA else ncol(rv$hmsc$newdata), Status=if(is.null(rv$hmsc$newdata)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="prediction_studyDesign.csv", Required="Conditional", Role="Prediction random-effect design", Rows=if(is.null(rv$hmsc$newstudy)) NA else nrow(rv$hmsc$newstudy), Columns=if(is.null(rv$hmsc$newstudy)) NA else ncol(rv$hmsc$newstudy), Status=if(is.null(rv$hmsc$newstudy)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="newcoords.csv", Required="Conditional", Role="Prediction spatial coordinates", Rows=if(is.null(rv$hmsc$newcoord)) NA else nrow(rv$hmsc$newcoord), Columns=if(is.null(rv$hmsc$newcoord)) NA else ncol(rv$hmsc$newcoord), Status=if(is.null(rv$hmsc$newcoord)) "Not uploaded or unreadable" else "Readable")
    )
    write.csv(data_summary(rv$hmsc$Y, rv$hmsc$X, rv$hmsc$Tr, rv$hmsc$study, rv$hmsc$coord, extra),
              file.path(outdir, "tables", "data_summary.csv"), row.names=FALSE)
    msgs <- if (!is.null(rv$hmsc$check)) rv$hmsc$check$messages else "HMSC check failed before messages were created."
    write_data_check_messages(outdir, msgs)
    status <- list(engine="Hmsc", status="check_failed", runtime_seconds=0,
                   warnings=msgs, errors="HMSC workflow stopped because input data or settings did not pass validation.")
    write_engine_status(outdir, status)
    writeLines(c(
      "HMSC check failed",
      "=================",
      "The application created this ZIP so the problem can be diagnosed instead of returning an empty download.",
      "",
      "Open diagnostics/data_check_messages.csv first.",
      "Common reasons:",
      "- Y.csv was not uploaded.",
      "- XData.csv was not uploaded.",
      "- Y and XData row numbers do not match.",
      "- distr = probit but Y is not 0/1.",
      "- XFormula contains variables not present in XData.csv.",
      "- random effect mode requires studyDesign.csv or coordinates.csv."
    ), file.path(outdir, "results", "README_HMSC_CHECK_FAILED.txt"))
    write_standard_outputs(outdir, "Hmsc", status, rv$hmsc$Y, rv$hmsc$X)
    rv$hmsc$zip <- tryCatch(make_zip(outdir), error = function(e) {
      writeLines(e$message, file.path(outdir, "diagnostics", "ZIP_creation_error.txt"))
      NULL
    })
    rv$hmsc$status <- "check_failed"
    add_log("hmsc", "HMSC check failed; diagnostic ZIP created:", outdir)
    invisible(outdir)
  }


  observeEvent(input$hmsc_run, {
    refresh_hmsc_check()
    if (!isTRUE(rv$hmsc$check$ok)) {
      rv$hmsc$status <- "Stopped: Hmsc data check failed"
      add_log("hmsc", "Stopped: Hmsc data check failed. Creating diagnostic ZIP.")
      create_hmsc_check_failed_zip()
      return(NULL)
    }
    run_safe_engine("Hmsc")
  })

  observeEvent(input$hmschpc_run, {
    refresh_hmschpc_check()
    if (!isTRUE(rv$hmschpc$check$ok)) {
      diag <- create_check_failed_zip("Hmsc-HPC", input$project_name, hmschpc_config(), rv$hmschpc$check,
        files = list("Y.csv" = input$hmschpc_Y_file, "XData.csv" = input$hmschpc_X_file, "traits.csv" = input$hmschpc_traits_file,
                     "studyDesign.csv" = input$hmschpc_study_file, "coordinates.csv" = input$hmschpc_coord_file,
                     "taxonomy.csv" = input$hmschpc_taxonomy_file, "phylo_tree.nwk" = input$hmschpc_phylo_tree_file,
                     "newdata.csv" = input$hmschpc_newdata_file),
        Y = rv$hmschpc$Y, X = rv$hmschpc$X)
      rv$hmschpc$outdir <- diag$outdir
      rv$hmschpc$zip <- diag$zip
      rv$hmschpc$status <- "check_failed"
      add_log("hmschpc", "Hmsc-HPC check failed; diagnostic ZIP created:", diag$outdir)
      return(NULL)
    }
    outdir <- make_engine_run_dir("Hmsc-HPC", input$project_name)
    rv$hmschpc$outdir <- outdir
    rv$hmschpc$status <- "Running"
    add_log("hmschpc", "Created Hmsc-HPC output folder:", outdir)
    status <- run_hmschpc_workflow(
      outdir, hmschpc_config(), rv$hmschpc$Y, rv$hmschpc$X,
      Tr = if (isTRUE(input$hmschpc_use_traits)) rv$hmschpc$Tr else NULL,
      study = rv$hmschpc$study,
      coord = rv$hmschpc$coord,
      taxonomy = if (identical(input$hmschpc_phylogeny_mode, "taxonomy")) rv$hmschpc$taxonomy else NULL,
      phylo_tree_file = if (identical(input$hmschpc_phylogeny_mode, "newick")) rv$hmschpc$phylo_tree else NULL,
      newdata = rv$hmschpc$newdata,
      log_fun = function(txt) add_log("hmschpc", txt)
    )
    status <- normalize_engine_status(status)
    rv$hmschpc$status <- status$status %||% "unknown"
    if (!file.exists(file.path(outdir, "standard", "run_summary.csv"))) {
      write_standard_outputs(outdir, "Hmsc-HPC", status, rv$hmschpc$Y, rv$hmschpc$X)
      write_engine_scaffold_outputs(outdir, "Hmsc-HPC", hmschpc_config(), status, rv$hmschpc$Y, rv$hmschpc$X)
    }
    ensure_output_contract(outdir, "Hmsc-HPC", status, rv$hmschpc$Y, rv$hmschpc$X)
    if (identical(rv$hmschpc$status, "fitted")) writeLines("RUN COMPLETE", file.path(outdir, "RUN_COMPLETE.txt"))
    if (identical(rv$hmschpc$status, "model_defined")) writeLines("MODEL DEFINED - SAMPLER SKIPPED", file.path(outdir, "RUN_MODEL_DEFINED.txt"))
    if (identical(rv$hmschpc$status, "check_failed")) writeLines("CHECK FAILED - NO FIT WAS ATTEMPTED", file.path(outdir, "RUN_CHECK_FAILED.txt"))
    if (identical(rv$hmschpc$status, "fit_failed")) writeLines("RUN FAILED", file.path(outdir, "RUN_FAILED.txt"))
    rv$hmschpc$zip <- if (isTRUE(input$hmschpc_out_zip)) make_zip(outdir) else NULL
    add_log("hmschpc", "Hmsc-HPC workflow status:", rv$hmschpc$status, outdir)
  })

  observeEvent(input$jsdm_run, {
    refresh_jsdm_check()
    if (!isTRUE(rv$jsdm$check$ok)) {
      diag <- create_check_failed_zip("jSDM", input$project_name, jsdm_config(), rv$jsdm$check,
        files = list("Y.csv" = input$jsdm_Y_file, "XData.csv" = input$jsdm_X_file, "trait_data.csv" = input$jsdm_trait_file,
                     "long_format.csv" = input$jsdm_long_file, "trials.csv" = input$jsdm_trials_file,
                     "newdata.csv" = input$jsdm_newdata_file, "prediction_ids.csv" = input$jsdm_prediction_ids_file),
        Y = rv$jsdm$Y, X = rv$jsdm$X)
      rv$jsdm$outdir <- diag$outdir
      rv$jsdm$zip <- diag$zip
      rv$jsdm$status <- "check_failed"
      add_log("jsdm", "jSDM check failed; diagnostic ZIP created:", diag$outdir)
      return(NULL)
    }
    run_safe_engine("jSDM")
  })

  output$hmsc_log <- renderText(paste(rv$hmsc$log, collapse="\n"))
  output$hmschpc_log <- renderText(paste(rv$hmschpc$log, collapse="\n"))
  output$jsdm_log <- renderText(paste(rv$jsdm$log, collapse="\n"))

  output$hmsc_run_summary <- renderText({
    check_msg <- if (!is.null(rv$hmsc$check) && !isTRUE(rv$hmsc$check$ok)) paste0("\n\nCheck messages:\n", paste(rv$hmsc$check$messages, collapse="\n")) else ""
    paste0("Status: ", rv$hmsc$status, "\nOutput folder: ", rv$hmsc$outdir %||% "None yet", "\nZIP: ", rv$hmsc$zip %||% "None yet", check_msg)
  })
  output$jsdm_run_summary <- renderText({
    paste0("Status: ", rv$jsdm$status, "\nOutput folder: ", rv$jsdm$outdir %||% "None yet", "\nZIP: ", rv$jsdm$zip %||% "None yet")
  })
  output$hmschpc_run_summary <- renderText({
    check_msg <- if (!is.null(rv$hmschpc$check) && !isTRUE(rv$hmschpc$check$ok)) paste0("\n\nCheck messages:\n", paste(rv$hmschpc$check$messages, collapse="\n")) else ""
    paste0("Status: ", rv$hmschpc$status, "\nOutput folder: ", rv$hmschpc$outdir %||% "None yet", "\nZIP: ", rv$hmschpc$zip %||% "None yet", check_msg)
  })

  output$hmsc_files <- renderDT({
    if (is.null(rv$hmsc$outdir) || !dir.exists(rv$hmsc$outdir)) return(datatable(data.frame(File=character(), Size_KB=numeric()), options=list(dom="tip"), rownames=FALSE))
    files <- list.files(rv$hmsc$outdir, recursive=TRUE, full.names=TRUE)
    dat <- data.frame(File=gsub(paste0("^", normalizePath(rv$hmsc$outdir, winslash="/"), "/?"), "", normalizePath(files, winslash="/", mustWork=FALSE)),
                      Size_KB=round(file.info(files)$size/1024,1))
    datatable(dat, options=list(pageLength=12), rownames=FALSE)
  })

  output$jsdm_files <- renderDT({
    if (is.null(rv$jsdm$outdir) || !dir.exists(rv$jsdm$outdir)) return(datatable(data.frame(File=character(), Size_KB=numeric()), options=list(dom="tip"), rownames=FALSE))
    files <- list.files(rv$jsdm$outdir, recursive=TRUE, full.names=TRUE)
    dat <- data.frame(File=gsub(paste0("^", normalizePath(rv$jsdm$outdir, winslash="/"), "/?"), "", normalizePath(files, winslash="/", mustWork=FALSE)),
                      Size_KB=round(file.info(files)$size/1024,1))
    datatable(dat, options=list(pageLength=12), rownames=FALSE)
  })

  output$hmschpc_files <- renderDT({
    if (is.null(rv$hmschpc$outdir) || !dir.exists(rv$hmschpc$outdir)) return(datatable(data.frame(File=character(), Size_KB=numeric()), options=list(dom="tip"), rownames=FALSE))
    files <- list.files(rv$hmschpc$outdir, recursive=TRUE, full.names=TRUE)
    dat <- data.frame(File=gsub(paste0("^", normalizePath(rv$hmschpc$outdir, winslash="/"), "/?"), "", normalizePath(files, winslash="/", mustWork=FALSE)),
                      Size_KB=round(file.info(files)$size/1024,1))
    datatable(dat, options=list(pageLength=12), rownames=FALSE)
  })

  output$hmsc_download <- downloadHandler(
    filename = function() if (is.null(rv$hmsc$zip)) paste0("Hmsc_results_", timestamp_id(), ".zip") else basename(rv$hmsc$zip),
    content = function(file) {
      copy_zip_to_download(rv$hmsc$zip, rv$hmsc$outdir, "Hmsc", file, input$project_name)
    }
  )
  output$jsdm_download <- downloadHandler(
    filename = function() if (is.null(rv$jsdm$zip)) paste0("jSDM_results_", timestamp_id(), ".zip") else basename(rv$jsdm$zip),
    content = function(file) {
      copy_zip_to_download(rv$jsdm$zip, rv$jsdm$outdir, "jSDM", file, input$project_name)
    }
  )
  output$hmschpc_download <- downloadHandler(
    filename = function() if (is.null(rv$hmschpc$zip)) paste0("Hmsc-HPC_results_", timestamp_id(), ".zip") else basename(rv$hmschpc$zip),
    content = function(file) {
      copy_zip_to_download(rv$hmschpc$zip, rv$hmschpc$outdir, "Hmsc-HPC", file, input$project_name)
    }
  )


  # GJAM check
  refresh_gjam_check <- function() {
    rv$gjam$Y <- read_csv_safe(input$gjam_Y_file$datapath)
    rv$gjam$X <- read_csv_safe(input$gjam_X_file$datapath)
    rv$gjam$types <- read_csv_safe(input$gjam_type_file$datapath)
    rv$gjam$censor <- read_csv_safe(input$gjam_censor_file$datapath)
    rv$gjam$effort <- read_csv_safe(input$gjam_effort_file$datapath)
    rv$gjam$newdata <- read_csv_safe(input$gjam_newdata_file$datapath)
    rv$gjam$specByTrait <- read_csv_safe(input$gjam_trait_spec_file$datapath)
    rv$gjam$traitTypes <- read_csv_safe(input$gjam_trait_types_file$datapath)
    rv$gjam$holdoutIndex <- read_csv_safe(input$gjam_holdout_file$datapath)
    rv$gjam$prior <- read_csv_safe(input$gjam_prior_file$datapath)
    rv$gjam$check <- validate_gjam_full(
      rv$gjam$Y, rv$gjam$X, rv$gjam$types, input$gjam_typeNames_text,
      input$gjam_type_single, input$gjam_FCgroups, input$gjam_CCgroups,
      input$gjam_ng, input$gjam_burnin, input$gjam_holdoutN,
      input$gjam_random, input$gjam_notStandard, input$gjam_formula,
      rv$gjam$censor, rv$gjam$effort, rv$gjam$newdata, rv$gjam$specByTrait,
      rv$gjam$traitTypes, rv$gjam$holdoutIndex, rv$gjam$prior,
      input$gjam_USE_CENSOR, input$gjam_USE_EFFORT, input$gjam_do_traits,
      input$gjam_trimY, input$gjam_trim_minObs, input$gjam_REDUCT,
      input$gjam_reduct_N, input$gjam_reduct_r
    )
    rv$gjam$checked <- TRUE
  }

  observeEvent(input$gjam_check, {
    refresh_gjam_check()
    add_log("gjam", "GJAM data and settings check completed.")
  })

  output$gjam_check_messages <- renderText({
    if (!isTRUE(rv$gjam$checked)) return("GJAM data have not been checked yet.")
    paste(rv$gjam$check$messages, collapse="\n")
  })

  output$gjam_data_table <- renderDT({
    extra <- list(
      data.frame(File="typeNames.csv", Required="Recommended", Role="Response type per Y column", Rows=if(is.null(rv$gjam$types)) NA else nrow(rv$gjam$types), Columns=if(is.null(rv$gjam$types)) NA else ncol(rv$gjam$types), Status=if(is.null(rv$gjam$types)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="censor.csv", Required="Optional", Role="Censoring values and intervals", Rows=if(is.null(rv$gjam$censor)) NA else nrow(rv$gjam$censor), Columns=if(is.null(rv$gjam$censor)) NA else ncol(rv$gjam$censor), Status=if(is.null(rv$gjam$censor)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="effort.csv", Required="Optional", Role="Sampling effort / plot area / sequencing depth", Rows=if(is.null(rv$gjam$effort)) NA else nrow(rv$gjam$effort), Columns=if(is.null(rv$gjam$effort)) NA else ncol(rv$gjam$effort), Status=if(is.null(rv$gjam$effort)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="newdata.csv", Required="Optional", Role="Prediction covariates", Rows=if(is.null(rv$gjam$newdata)) NA else nrow(rv$gjam$newdata), Columns=if(is.null(rv$gjam$newdata)) NA else ncol(rv$gjam$newdata), Status=if(is.null(rv$gjam$newdata)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="specByTrait.csv", Required="Traits optional", Role="Species-by-trait table", Rows=if(is.null(rv$gjam$specByTrait)) NA else nrow(rv$gjam$specByTrait), Columns=if(is.null(rv$gjam$specByTrait)) NA else ncol(rv$gjam$specByTrait), Status=if(is.null(rv$gjam$specByTrait)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="traitTypes.csv", Required="Traits optional", Role="Trait response type table", Rows=if(is.null(rv$gjam$traitTypes)) NA else nrow(rv$gjam$traitTypes), Columns=if(is.null(rv$gjam$traitTypes)) NA else ncol(rv$gjam$traitTypes), Status=if(is.null(rv$gjam$traitTypes)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="holdoutIndex.csv", Required="Optional", Role="Explicit holdout row indices", Rows=if(is.null(rv$gjam$holdoutIndex)) NA else nrow(rv$gjam$holdoutIndex), Columns=if(is.null(rv$gjam$holdoutIndex)) NA else ncol(rv$gjam$holdoutIndex), Status=if(is.null(rv$gjam$holdoutIndex)) "Not uploaded or unreadable" else "Readable")
    )
    datatable(data_summary(rv$gjam$Y, rv$gjam$X, NULL, NULL, NULL, extra),
              options = list(dom="tip", pageLength=10), rownames=FALSE)
  })

  gjam_config <- reactive({
    list(
      project_name = input$project_name,
      project_question = input$project_question,
      question_template = input$question_template,
      project_design = project_profile(),
      engine = "GJAM",
      data = list(
        Y = if(!is.null(input$gjam_Y_file$name)) input$gjam_Y_file$name else NULL,
        XData = if(!is.null(input$gjam_X_file$name)) input$gjam_X_file$name else NULL,
        typeNames = if(!is.null(input$gjam_type_file$name)) input$gjam_type_file$name else NULL,
        censor = if(!is.null(input$gjam_censor_file$name)) input$gjam_censor_file$name else NULL,
        effort = if(!is.null(input$gjam_effort_file$name)) input$gjam_effort_file$name else NULL,
        newdata = if(!is.null(input$gjam_newdata_file$name)) input$gjam_newdata_file$name else NULL,
        specByTrait = if(!is.null(input$gjam_trait_spec_file$name)) input$gjam_trait_spec_file$name else NULL,
        traitTypes = if(!is.null(input$gjam_trait_types_file$name)) input$gjam_trait_types_file$name else NULL,
        holdoutIndex = if(!is.null(input$gjam_holdout_file$name)) input$gjam_holdout_file$name else NULL
      ),
      model = list(
        formula = input$gjam_formula,
        type_single = input$gjam_type_single,
        typeNames_text = input$gjam_typeNames_text,
        notStandard = input$gjam_notStandard,
        ng = input$gjam_ng,
        burnin = input$gjam_burnin,
        holdoutN = input$gjam_holdoutN,
        seed = input$gjam_seed,
        random = input$gjam_random,
        FULL = input$gjam_FULL,
        PREDICTX = input$gjam_PREDICTX,
        REDUCT = input$gjam_REDUCT,
        reductList = list(N = input$gjam_reduct_N, r = input$gjam_reduct_r),
        ematAlpha = input$gjam_ematAlpha,
        use_censor = input$gjam_USE_CENSOR,
        use_effort = input$gjam_USE_EFFORT
      ),
      response_types = list(
        FCgroups = input$gjam_FCgroups,
        CCgroups = input$gjam_CCgroups,
        composition_reference = input$gjam_composition_reference,
        trimY = input$gjam_trimY,
        trim_minObs = input$gjam_trim_minObs
      ),
      priors_censoring = list(
        use_prior_template = input$gjam_use_prior_template,
        prior_file = if(!is.null(input$gjam_prior_file$name)) input$gjam_prior_file$name else NULL,
        prior_mode = input$gjam_prior_mode,
        censor_columns = input$gjam_censor_columns,
        censor_values = input$gjam_censor_values,
        censor_intervals = input$gjam_censor_intervals
      ),
      analysis = list(
        do_predict = input$gjam_do_predict,
        do_sensitivity = input$gjam_do_sensitivity,
        do_ordination = input$gjam_do_ordination,
        do_conditional = input$gjam_do_conditional,
        do_iie = input$gjam_do_iie,
        do_traits = input$gjam_do_traits,
        missingX = input$gjam_missingX,
        missingY = input$gjam_missingY,
        inverse_prediction = input$gjam_inverse_prediction
      ),
      outputs = list(
        real_fit = input$gjam_real_fit,
        save_model = input$gjam_out_model,
        save_chains = input$gjam_out_chains,
        save_parameters = input$gjam_out_parameters,
        save_fit = input$gjam_out_fit,
        save_prediction = input$gjam_out_prediction,
        save_missing = input$gjam_out_missing,
        save_plots = input$gjam_out_plots,
        report = input$gjam_out_report,
        zip = input$gjam_out_zip,
        copy_inputs = input$gjam_out_copy_inputs,
        save_config = input$gjam_out_config,
        csv_tables = TRUE
      )
    )
  })

  run_gjam_safe <- function() {
    cfg <- gjam_config()
    outdir <- make_engine_run_dir("GJAM", input$project_name)
    rv$gjam$outdir <- outdir
    rv$gjam$status <- "Running"
    add_log("gjam", "Created output folder:", outdir)

    if (isTRUE(input$gjam_out_copy_inputs) || isTRUE(cfg$outputs$real_fit %||% FALSE)) {
      copy_upload(input$gjam_Y_file, outdir, "Y.csv")
      copy_upload(input$gjam_X_file, outdir, "XData.csv")
      copy_upload(input$gjam_type_file, outdir, "typeNames.csv")
      copy_upload(input$gjam_censor_file, outdir, "censor.csv")
      copy_upload(input$gjam_effort_file, outdir, "effort.csv")
      copy_upload(input$gjam_newdata_file, outdir, "newdata.csv")
      copy_upload(input$gjam_trait_spec_file, outdir, "specByTrait.csv")
      copy_upload(input$gjam_trait_types_file, outdir, "traitTypes.csv")
      copy_upload(input$gjam_holdout_file, outdir, "holdoutIndex.csv")
      copy_upload(input$gjam_prior_file, outdir, "priorTemplate.csv")
    }
    if (isTRUE(input$gjam_out_config) || isTRUE(cfg$outputs$real_fit)) yaml::write_yaml(cfg, file.path(outdir, "used_config.yml"))

    write.csv(data_summary(rv$gjam$Y, rv$gjam$X), file.path(outdir, "tables", "data_summary.csv"), row.names=FALSE)
    write_data_check_messages(outdir, if (!is.null(rv$gjam$check)) rv$gjam$check$messages else character())

    status <- list(engine="GJAM", status="ready", runtime_seconds=0, warnings=character(), errors=character())
    if (cfg$model$burnin >= cfg$model$ng) {
      status$warnings <- c(status$warnings, "GJAM warning: burnin must be less than ng.")
    }
    if (cfg$model$holdoutN > 0 && !is.null(rv$gjam$Y) && cfg$model$holdoutN >= nrow(rv$gjam$Y)) {
      status$warnings <- c(status$warnings, "GJAM warning: holdoutN must be smaller than the number of rows in Y.")
    }
    if (!safe_require("gjam")) {
      status$status <- "check_failed"
      status$warnings <- c(status$warnings, "gjam package is not available. A diagnostic output scaffold was saved; no model fit was attempted.")
    } else if (isTRUE(cfg$outputs$real_fit %||% FALSE)) {
      status$status <- "model_defined"
    } else {
      status$status <- "model_defined"
      status$warnings <- c(status$warnings, "gjam package appears available, but production GJAM fitting is not connected in this build.")
    }
    status <- normalize_engine_status(status)

    gjam_files <- data.frame(
      Section = c("Configuration", "Inputs", "Model object", "Chains", "Parameters", "Fit diagnostics", "Prediction", "Sensitivity", "Ordination", "Missing data", "Plots", "Diagnostics", "Report"),
      Expected_file_or_folder = c(
        "used_config.yml",
        "inputs/Y.csv, inputs/XData.csv, inputs/typeNames.csv, inputs/censor.csv, inputs/effort.csv, inputs/newdata.csv",
        "models/gjam_model.rds",
        "chains/bgibbs.rds, bgibbsUn.rds, fgibbs.rds, fbgibbs.rds, sgibbs.rds, ygibbs.rds if FULL=TRUE",
        "tables/betaMu.csv, betaSe.csv, betaMuUn.csv, betaSeUn.csv, fBetaMu.csv, fBetaSd.csv, corMu.csv, corSe.csv, sigMu.csv, sigSe.csv, fmatrix.csv, fMu.csv, fSe.csv, ematrix.csv",
        "tables/fit_DIC_rmspe_xscore_yscore.csv",
        "predictions/richness.csv, ypred.csv, xpred.csv, prediction_uncertainty.csv",
        "tables/sensitivity_*.csv and plots/sensitivity.pdf",
        "tables/ordination_scores.csv and plots/ordination.pdf",
        "tables/missing_x.csv, missing_y.csv, xmissMu.csv, ymissMu.csv",
        "plots/gjamPlot.pdf, sensitivity.pdf, ordination.pdf, IIEplot.pdf",
        "diagnostics/engine_status.json, diagnostics/data_check_messages.csv",
        "report/GJAM_report.html"
      ),
      Meaning = c(
        "Exact GJAM settings used for the run.",
        "Copied GJAM input files.",
        "Serialized fitted GJAM object when production fitting is connected.",
        "Raw Gibbs-sampling chains from the GJAM object.",
        "Parameter estimates and covariance/correlation/sensitivity matrices on the observation scale.",
        "Model fit diagnostics such as DIC and RMSPE.",
        "Predicted responses, richness and inverse predictions.",
        "Sensitivity of joint responses to predictors.",
        "Ordination outputs for response structure.",
        "Predictions and uncertainty for missing X/Y values.",
        "Human-readable figures from GJAM plotting functions.",
        "Machine-readable status, warnings and data checks.",
        "HTML report for browser viewing."
      )
    )
    write.csv(gjam_files, file.path(outdir, "tables", "GJAM_result_workflow_map.csv"), row.names = FALSE)
    writeLines(c(
      "GJAM result workflow",
      "====================",
      "This workflow exports a standalone R script and, when enabled, runs the installed gjam package.",
      "The reproducible script writes the files listed in tables/GJAM_result_workflow_map.csv.",
      "",
      paste(gjam_files$Section, gjam_files$Expected_file_or_folder, sep = " -> ")
    ), file.path(outdir, "results", "README_GJAM_results.txt"))

    if (exists("write_gjam_reproducible_script", mode = "function")) {
      write_gjam_reproducible_script(outdir)
    } else {
      status$status <- "fit_failed"
      status$errors <- c(status$errors, "Internal GJAM adapter function write_gjam_reproducible_script() is missing.")
      writeLines(status$errors, file.path(outdir, "diagnostics", "GJAM_adapter_error.txt"))
    }

    if (isTRUE(cfg$outputs$real_fit %||% FALSE)) {
      if (!safe_require("gjam")) {
        status$status <- "fit_failed"
        status$errors <- c(status$errors, "Package gjam is not installed, so real fitting cannot run.")
        writeLines(status$errors, file.path(outdir, "diagnostics", "GJAM_reproducible_error.txt"))
        write_standard_outputs(outdir, "GJAM", status, rv$gjam$Y, rv$gjam$X)
        write_engine_status(outdir, status)
      } else if (file.exists(file.path(outdir, "reproducible_script", "run_this_GJAM_analysis.R"))) {
        add_log("gjam", "Starting real GJAM package workflow.")
        script_file <- file.path(outdir, "reproducible_script", "run_this_GJAM_analysis.R")
        log_file <- file.path(outdir, "diagnostics", "GJAM_real_fit_stdout_stderr.txt")
        rscript_bin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
        started <- Sys.time()
        exit_code <- tryCatch({
          neutralize_windows_r_locale_env()
          system2(rscript_bin, shQuote(script_file), stdout = log_file, stderr = log_file)
        },
          error = function(e) {
            writeLines(e$message, file.path(outdir, "diagnostics", "GJAM_system2_error.txt"))
            1L
          }
        )
        status$runtime_seconds <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2)
        engine_json <- file.path(outdir, "diagnostics", "engine_status.json")
        if (file.exists(engine_json) && requireNamespace("jsonlite", quietly = TRUE)) {
          script_status <- tryCatch(jsonlite::fromJSON(engine_json), error = function(e) NULL)
          if (!is.null(script_status)) {
            if (!is.null(script_status$status)) status$status <- as.character(script_status$status)
            if (!is.null(script_status$warnings)) status$warnings <- unique(c(status$warnings, as.character(script_status$warnings)))
            if (!is.null(script_status$errors)) status$errors <- unique(c(status$errors, as.character(script_status$errors)))
          }
        }
        if (!identical(as.integer(exit_code), 0L) && !identical(status$status, "fitted")) {
          status$status <- "fit_failed"
          status$errors <- unique(c(status$errors, paste0("Rscript exited with status ", as.integer(exit_code), ". See diagnostics/GJAM_real_fit_stdout_stderr.txt.")))
        }
        write_engine_status(outdir, status)
        add_log("gjam", "Real GJAM package workflow status:", status$status)
      } else {
        status$status <- "fit_failed"
        status$errors <- c(status$errors, "The GJAM reproducible script was not created.")
        writeLines(status$errors, file.path(outdir, "diagnostics", "GJAM_reproducible_error.txt"))
        write_standard_outputs(outdir, "GJAM", status, rv$gjam$Y, rv$gjam$X)
        write_engine_status(outdir, status)
      }
    }

    status <- normalize_engine_status(status)
    if (!isTRUE(cfg$outputs$real_fit %||% FALSE)) {
      write_engine_status(outdir, status)
      write_standard_outputs(outdir, "GJAM", status, rv$gjam$Y, rv$gjam$X)
      write_reproducible_stub(outdir, "GJAM")
      write_engine_scaffold_outputs(outdir, "GJAM", cfg, status, rv$gjam$Y, rv$gjam$X)
    }

    ensure_output_contract(outdir, "GJAM", status, rv$gjam$Y, rv$gjam$X)
    make_html_report(outdir, "GJAM", cfg, paste(status$status, paste(status$warnings, collapse="; ")))
    if (identical(status$status, "fitted")) {
      writeLines("RUN COMPLETE", file.path(outdir, "RUN_COMPLETE.txt"))
    } else if (identical(status$status, "fit_failed")) {
      writeLines("RUN FAILED", file.path(outdir, "RUN_FAILED.txt"))
    } else if (identical(status$status, "model_defined")) {
      writeLines("MODEL DEFINED - NO FITTED POSTERIOR", file.path(outdir, "RUN_MODEL_DEFINED.txt"))
    } else if (identical(status$status, "check_failed")) {
      writeLines("CHECK FAILED - NO FIT WAS ATTEMPTED", file.path(outdir, "RUN_CHECK_FAILED.txt"))
    } else {
      writeLines("FIT FAILED - UNRECOGNIZED TERMINAL STATUS", file.path(outdir, "RUN_FAILED.txt"))
    }
    rv$gjam$zip <- if (isTRUE(input$gjam_out_zip)) make_zip(outdir) else NULL
    rv$gjam$status <- status$status
    add_log("gjam", "GJAM workflow status:", status$status, outdir)
  }

  observeEvent(input$gjam_run, {
    refresh_gjam_check()
    if (!isTRUE(rv$gjam$check$ok)) {
      diag <- create_check_failed_zip("GJAM", input$project_name, gjam_config(), rv$gjam$check,
        files = list("Y.csv" = input$gjam_Y_file, "XData.csv" = input$gjam_X_file, "typeNames.csv" = input$gjam_type_file,
                     "censor.csv" = input$gjam_censor_file, "effort.csv" = input$gjam_effort_file,
                     "newdata.csv" = input$gjam_newdata_file, "specByTrait.csv" = input$gjam_trait_spec_file,
                     "traitTypes.csv" = input$gjam_trait_types_file, "holdoutIndex.csv" = input$gjam_holdout_file),
        Y = rv$gjam$Y, X = rv$gjam$X)
      rv$gjam$outdir <- diag$outdir
      rv$gjam$zip <- diag$zip
      rv$gjam$status <- "check_failed"
      add_log("gjam", "GJAM check failed; diagnostic ZIP created:", diag$outdir)
      return(NULL)
    }
    run_gjam_safe()
  })

  output$gjam_log <- renderText(paste(rv$gjam$log, collapse="\n"))
  output$gjam_run_summary <- renderText({
    paste0("Status: ", rv$gjam$status, "\nOutput folder: ", rv$gjam$outdir %||% "None yet", "\nZIP: ", rv$gjam$zip %||% "None yet")
  })
  output$gjam_files <- renderDT({
    if (is.null(rv$gjam$outdir) || !dir.exists(rv$gjam$outdir)) return(datatable(data.frame(File=character(), Size_KB=numeric()), options=list(dom="tip"), rownames=FALSE))
    files <- list.files(rv$gjam$outdir, recursive=TRUE, full.names=TRUE)
    dat <- data.frame(File=gsub(paste0("^", normalizePath(rv$gjam$outdir, winslash="/"), "/?"), "", normalizePath(files, winslash="/", mustWork=FALSE)),
                      Size_KB=round(file.info(files)$size/1024,1))
    datatable(dat, options=list(pageLength=12), rownames=FALSE)
  })
  output$gjam_download <- downloadHandler(
    filename = function() if (is.null(rv$gjam$zip)) paste0("GJAM_results_", timestamp_id(), ".zip") else basename(rv$gjam$zip),
    content = function(file) {
      copy_zip_to_download(rv$gjam$zip, rv$gjam$outdir, "GJAM", file, input$project_name)
    }
  )



  # spOccupancy check
  refresh_spocc_check <- function() {
    rv$spocc$y <- read_csv_safe(input$spocc_y_file$datapath)
    rv$spocc$occ <- read_csv_safe(input$spocc_occ_file$datapath)
    rv$spocc$det <- read_csv_safe(input$spocc_det_file$datapath)
    rv$spocc$coords <- read_csv_safe(input$spocc_coords_file$datapath)
    rv$spocc$species <- read_csv_safe(input$spocc_species_file$datapath)
    rv$spocc$integrated <- read_csv_safe(input$spocc_integrated_file$datapath)
    rv$spocc$newdata <- read_csv_safe(input$spocc_newdata_file$datapath)
    rv$spocc$newcoords <- read_csv_safe(input$spocc_newcoords_file$datapath)
    rv$spocc$check <- validate_spoccupancy_full(
      rv$spocc$y, rv$spocc$occ, rv$spocc$det, rv$spocc$coords, rv$spocc$species, rv$spocc$integrated,
      if (isTRUE(input$spocc_do_predict)) rv$spocc$newdata else NULL, rv$spocc$newcoords, NULL,
      input$spocc_model_type, input$spocc_n_batch, input$spocc_batch_length, input$spocc_n_burn,
      input$spocc_n_thin, input$spocc_n_chains, input$spocc_n_factors, input$spocc_NNGP,
      input$spocc_n_neighbors, 0, input$spocc_svc_cols, input$spocc_occ_formula,
      input$spocc_det_formula
    )
    rv$spocc$checked <- TRUE
  }

  observeEvent(input$spocc_check, {
    refresh_spocc_check()
    add_log("spocc", "spOccupancy data and settings check completed.")
  })

  output$spocc_check_messages <- renderText({
    if (!isTRUE(rv$spocc$checked)) return("spOccupancy data have not been checked yet.")
    paste(rv$spocc$check$messages, collapse="\n")
  })

  output$spocc_data_table <- renderDT({
    extra <- list(
      data.frame(File="det.covs.csv", Required="Recommended", Role="Detection covariates / replicate-level covariates", Rows=if(is.null(rv$spocc$det)) NA else nrow(rv$spocc$det), Columns=if(is.null(rv$spocc$det)) NA else ncol(rv$spocc$det), Status=if(is.null(rv$spocc$det)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="coords.csv", Required="Spatial models", Role="Spatial coordinates", Rows=if(is.null(rv$spocc$coords)) NA else nrow(rv$spocc$coords), Columns=if(is.null(rv$spocc$coords)) NA else ncol(rv$spocc$coords), Status=if(is.null(rv$spocc$coords)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="species.csv", Required="Multi-species optional", Role="Species names / metadata", Rows=if(is.null(rv$spocc$species)) NA else nrow(rv$spocc$species), Columns=if(is.null(rv$spocc$species)) NA else ncol(rv$spocc$species), Status=if(is.null(rv$spocc$species)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="integrated_sources.csv", Required="Integrated models", Role="Data-source metadata", Rows=if(is.null(rv$spocc$integrated)) NA else nrow(rv$spocc$integrated), Columns=if(is.null(rv$spocc$integrated)) NA else ncol(rv$spocc$integrated), Status=if(is.null(rv$spocc$integrated)) "Not uploaded or unreadable" else "Readable")
    )
    datatable(data_summary(rv$spocc$y, rv$spocc$occ, NULL, NULL, NULL, extra),
              options = list(dom="tip", pageLength=10), rownames=FALSE)
  })

  spocc_config <- reactive({
    list(
      project_name = input$project_name,
      project_question = input$project_question,
      question_template = input$question_template,
      project_design = project_profile(),
      engine = "spOccupancy",
      data = list(
        y = if(!is.null(input$spocc_y_file$name)) input$spocc_y_file$name else NULL,
        occ.covs = if(!is.null(input$spocc_occ_file$name)) input$spocc_occ_file$name else NULL,
        det.covs = if(!is.null(input$spocc_det_file$name)) input$spocc_det_file$name else NULL,
        coords = if(!is.null(input$spocc_coords_file$name)) input$spocc_coords_file$name else NULL,
        species = if(!is.null(input$spocc_species_file$name)) input$spocc_species_file$name else NULL,
        integrated_sources = if(!is.null(input$spocc_integrated_file$name)) input$spocc_integrated_file$name else NULL,
        newdata = if(!is.null(input$spocc_newdata_file$name)) input$spocc_newdata_file$name else NULL,
        newcoords = if(!is.null(input$spocc_newcoords_file$name)) input$spocc_newcoords_file$name else NULL
      ),
      model = list(
        model_type = input$spocc_model_type,
        occ.formula = input$spocc_occ_formula,
        det.formula = input$spocc_det_formula
      ),
      spatial_latent_svc = list(
        cov.model = input$spocc_cov_model,
        NNGP = input$spocc_NNGP,
        n.neighbors = input$spocc_n_neighbors,
        search.type = input$spocc_search_type,
        n.factors = input$spocc_n_factors,
        svc.cols = input$spocc_svc_cols
      ),
      mcmc = list(
        n.batch = input$spocc_n_batch,
        batch.length = input$spocc_batch_length,
        n.burn = input$spocc_n_burn,
        n.thin = input$spocc_n_thin,
        n.chains = input$spocc_n_chains,
        accept.rate = input$spocc_accept_rate,
        n.report = input$spocc_n_report,
        n.omp.threads = input$spocc_n_omp_threads,
        verbose = input$spocc_verbose,
        seed = input$spocc_seed
      ),
      validation_prediction_outputs = list(
        real_fit = input$spocc_real_fit,
        ppcOcc = input$spocc_do_ppc,
        waicOcc = input$spocc_do_waic,
        predict = input$spocc_do_predict,
        save_model = input$spocc_save_model,
        save_samples = input$spocc_save_samples,
        zip = input$spocc_out_zip,
        copy_inputs = input$spocc_out_copy_inputs,
        save_config = input$spocc_out_config,
        report = input$spocc_out_report
      )
    )
  })

  run_spocc_safe <- function() {
    cfg <- spocc_config()
    outdir <- make_engine_run_dir("spOccupancy", input$project_name)
    rv$spocc$outdir <- outdir
    rv$spocc$status <- "Running"
    add_log("spocc", "Created output folder:", outdir)

    if (isTRUE(input$spocc_out_copy_inputs) || isTRUE(cfg$validation_prediction_outputs$real_fit %||% FALSE)) {
      copy_upload(input$spocc_y_file, outdir, "y.csv")
      copy_upload(input$spocc_occ_file, outdir, "occ.covs.csv")
      copy_upload(input$spocc_det_file, outdir, "det.covs.csv")
      copy_upload(input$spocc_coords_file, outdir, "coords.csv")
      copy_upload(input$spocc_species_file, outdir, "species.csv")
      copy_upload(input$spocc_integrated_file, outdir, "integrated_sources.csv")
      copy_upload(input$spocc_newdata_file, outdir, "newdata.csv")
      copy_upload(input$spocc_newcoords_file, outdir, "newcoords.csv")
    }
    if (isTRUE(input$spocc_out_config) || isTRUE(cfg$validation_prediction_outputs$real_fit)) yaml::write_yaml(cfg, file.path(outdir, "used_config.yml"))
    write.csv(data_summary(rv$spocc$y, rv$spocc$occ), file.path(outdir, "tables", "data_summary.csv"), row.names=FALSE)
    write_data_check_messages(outdir, if (!is.null(rv$spocc$check)) rv$spocc$check$messages else character())

    status <- list(engine="spOccupancy", status="ready", runtime_seconds=0, warnings=character(), errors=character())
    total <- cfg$mcmc$n.batch * cfg$mcmc$batch.length
    if (cfg$mcmc$n.burn >= total) status$warnings <- c(status$warnings, "spOccupancy warning: n.burn must be less than n.batch * batch.length.")
    if (cfg$model$model_type %in% c("spPGOcc","spMsPGOcc","sfJSDM","sfMsPGOcc","spIntPGOcc","stPGOcc","stMsPGOcc","svcPGBinom","svcPGOcc","svcMsPGOcc","svcTPGBinom","svcTPGOcc","svcTMsPGOcc") && is.null(rv$spocc$coords)) {
      status$warnings <- c(status$warnings, "spOccupancy warning: selected spatial model requires coords.csv.")
    }
    if (!safe_require("spOccupancy")) {
      status$status <- "check_failed"
      status$warnings <- c(status$warnings, "spOccupancy package is not available. A diagnostic output scaffold was saved; no model fit was attempted.")
    } else if (isTRUE(cfg$validation_prediction_outputs$real_fit %||% FALSE)) {
      status$status <- "model_defined"
    } else {
      status$status <- "model_defined"
      status$warnings <- c(status$warnings, "spOccupancy package appears available, but real fitting is disabled for this run.")
    }
    status <- normalize_engine_status(status)

    spocc_files <- data.frame(
      Section = c("Configuration", "Inputs", "Model object", "Posterior samples", "Summary tables", "Fitted values", "Prediction", "Model assessment", "Spatial outputs", "Diagnostics", "Report"),
      Expected_file_or_folder = c(
        "used_config.yml",
        "inputs/y.csv, occ.covs.csv, det.covs.csv, coords.csv, integrated_sources.csv, newdata.csv",
        "models/spOccupancy_model.rds",
        "samples/beta_samples.rds, alpha_samples.rds, z_samples.rds, psi_samples.rds, p_samples.rds, w_samples.rds, lambda_samples.rds, svc_samples.rds",
        "tables/summary_beta.csv, summary_alpha.csv, community_summary.csv, species_summary.csv, random_effects.csv, latent_factor_loadings.csv",
        "tables/fitted_values.csv or results/fitted_values.rds",
        "predictions/occupancy_predictions.csv, detection_predictions.csv, prediction_uncertainty.csv, maps-ready prediction tables",
        "tables/ppcOcc_results.csv, waicOcc_results.csv, kfold_results.csv",
        "tables/spatial_parameters.csv, spatial_random_effects.csv, SVC_surfaces.csv",
        "diagnostics/engine_status.json, data_check_messages.csv, convergence_Rhat_ESS.csv",
        "report/spOccupancy_report.html"
      ),
      Meaning = c(
        "Exact settings used for the spOccupancy workflow.",
        "Copied detection, occupancy, detection-covariate, spatial and prediction inputs.",
        "Serialized fitted spOccupancy object when production fitting is connected.",
        "Posterior samples from occupancy, detection, latent occurrence, spatial, latent-factor and SVC components depending on model type.",
        "Posterior means, intervals, Rhat and effective sample sizes from summary methods.",
        "Fitted detection/occupancy values for posterior predictive checks.",
        "Out-of-sample occurrence and detection predictions with uncertainty.",
        "Posterior predictive checks, WAIC and k-fold validation outputs.",
        "Spatial random effects, spatial parameters, latent factors and spatially varying coefficient outputs.",
        "Machine-readable status, warnings and data checks.",
        "HTML report for browser viewing."
      )
    )
    write.csv(spocc_files, file.path(outdir, "tables", "spOccupancy_result_workflow_map.csv"), row.names = FALSE)
    writeLines(c(
      "spOccupancy result workflow",
      "===========================",
      "This workflow exports a standalone R script and, when enabled, runs the installed spOccupancy package.",
      "The reproducible script writes the files listed in tables/spOccupancy_result_workflow_map.csv.",
      "",
      paste(spocc_files$Section, spocc_files$Expected_file_or_folder, sep = " -> ")
    ), file.path(outdir, "results", "README_spOccupancy_results.txt"))

    if (exists("write_spoccupancy_reproducible_script", mode = "function")) {
      write_spoccupancy_reproducible_script(outdir)
    } else {
      status$status <- "fit_failed"
      status$errors <- c(status$errors, "Internal spOccupancy adapter function write_spoccupancy_reproducible_script() is missing.")
      writeLines(status$errors, file.path(outdir, "diagnostics", "spOccupancy_adapter_error.txt"))
    }

    if (isTRUE(cfg$validation_prediction_outputs$real_fit %||% FALSE)) {
      if (!safe_require("spOccupancy")) {
        status$status <- "fit_failed"
        status$errors <- c(status$errors, "Package spOccupancy is not installed, so real fitting cannot run.")
        writeLines(status$errors, file.path(outdir, "diagnostics", "spOccupancy_reproducible_error.txt"))
        write_standard_outputs(outdir, "spOccupancy", status, rv$spocc$y, rv$spocc$occ)
        write_engine_status(outdir, status)
      } else if (file.exists(file.path(outdir, "reproducible_script", "run_this_spOccupancy_analysis.R"))) {
        add_log("spocc", "Starting real spOccupancy package workflow.")
        script_file <- file.path(outdir, "reproducible_script", "run_this_spOccupancy_analysis.R")
        log_file <- file.path(outdir, "diagnostics", "spOccupancy_real_fit_stdout_stderr.txt")
        rscript_bin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
        started <- Sys.time()
        exit_code <- tryCatch({
          neutralize_windows_r_locale_env()
          system2(rscript_bin, shQuote(script_file), stdout = log_file, stderr = log_file)
        },
          error = function(e) {
            writeLines(e$message, file.path(outdir, "diagnostics", "spOccupancy_system2_error.txt"))
            1L
          }
        )
        status$runtime_seconds <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2)
        engine_json <- file.path(outdir, "diagnostics", "engine_status.json")
        if (file.exists(engine_json) && requireNamespace("jsonlite", quietly = TRUE)) {
          script_status <- tryCatch(jsonlite::fromJSON(engine_json), error = function(e) NULL)
          if (!is.null(script_status)) {
            if (!is.null(script_status$status)) status$status <- as.character(script_status$status)
            if (!is.null(script_status$warnings)) status$warnings <- unique(c(status$warnings, as.character(script_status$warnings)))
            if (!is.null(script_status$errors)) status$errors <- unique(c(status$errors, as.character(script_status$errors)))
          }
        }
        if (!identical(as.integer(exit_code), 0L) && !identical(status$status, "fitted")) {
          status$status <- "fit_failed"
          status$errors <- unique(c(status$errors, paste0("Rscript exited with status ", as.integer(exit_code), ". See diagnostics/spOccupancy_real_fit_stdout_stderr.txt.")))
        }
        write_engine_status(outdir, status)
        add_log("spocc", "Real spOccupancy package workflow status:", status$status)
      } else {
        status$status <- "fit_failed"
        status$errors <- c(status$errors, "The spOccupancy reproducible script was not created.")
        writeLines(status$errors, file.path(outdir, "diagnostics", "spOccupancy_reproducible_error.txt"))
        write_standard_outputs(outdir, "spOccupancy", status, rv$spocc$y, rv$spocc$occ)
        write_engine_status(outdir, status)
      }
    }

    status <- normalize_engine_status(status)
    if (!isTRUE(cfg$validation_prediction_outputs$real_fit %||% FALSE)) {
      write_engine_status(outdir, status)
      write_standard_outputs(outdir, "spOccupancy", status, rv$spocc$y, rv$spocc$occ)
      write_reproducible_stub(outdir, "spOccupancy")
      write_engine_scaffold_outputs(outdir, "spOccupancy", cfg, status, rv$spocc$y, rv$spocc$occ)
    }

    ensure_output_contract(outdir, "spOccupancy", status, rv$spocc$y, rv$spocc$occ)
    if (isTRUE(cfg$validation_prediction_outputs$report)) {
      make_html_report(outdir, "spOccupancy", cfg, paste(status$status, paste(status$warnings, collapse="; ")))
    }
    if (identical(status$status, "fitted")) {
      writeLines("RUN COMPLETE", file.path(outdir, "RUN_COMPLETE.txt"))
    } else if (identical(status$status, "fit_failed")) {
      writeLines("RUN FAILED", file.path(outdir, "RUN_FAILED.txt"))
    } else if (identical(status$status, "model_defined")) {
      writeLines("MODEL DEFINED - NO FITTED POSTERIOR", file.path(outdir, "RUN_MODEL_DEFINED.txt"))
    } else if (identical(status$status, "check_failed")) {
      writeLines("CHECK FAILED - NO FIT WAS ATTEMPTED", file.path(outdir, "RUN_CHECK_FAILED.txt"))
    } else {
      writeLines("FIT FAILED - UNRECOGNIZED TERMINAL STATUS", file.path(outdir, "RUN_FAILED.txt"))
    }
    rv$spocc$zip <- if (isTRUE(input$spocc_out_zip)) make_zip(outdir) else NULL
    rv$spocc$status <- status$status
    add_log("spocc", "spOccupancy workflow status:", status$status, outdir)
  }

  observeEvent(input$spocc_run, {
    refresh_spocc_check()
    if (!isTRUE(rv$spocc$check$ok)) {
      diag <- create_check_failed_zip("spOccupancy", input$project_name, spocc_config(), rv$spocc$check,
        files = list("y.csv" = input$spocc_y_file, "occ.covs.csv" = input$spocc_occ_file, "det.covs.csv" = input$spocc_det_file,
                     "coords.csv" = input$spocc_coords_file, "species.csv" = input$spocc_species_file,
                     "integrated_sources.csv" = input$spocc_integrated_file, "newdata.csv" = input$spocc_newdata_file,
                     "newcoords.csv" = input$spocc_newcoords_file),
        Y = rv$spocc$y, X = rv$spocc$occ)
      rv$spocc$outdir <- diag$outdir
      rv$spocc$zip <- diag$zip
      rv$spocc$status <- "check_failed"
      add_log("spocc", "spOccupancy check failed; diagnostic ZIP created:", diag$outdir)
      return(NULL)
    }
    run_spocc_safe()
  })

  output$spocc_log <- renderText(paste(rv$spocc$log, collapse="\n"))
  output$spocc_run_summary <- renderText({
    paste0("Status: ", rv$spocc$status, "\nOutput folder: ", rv$spocc$outdir %||% "None yet", "\nZIP: ", rv$spocc$zip %||% "None yet")
  })
  output$spocc_files <- renderDT({
    if (is.null(rv$spocc$outdir) || !dir.exists(rv$spocc$outdir)) return(datatable(data.frame(File=character(), Size_KB=numeric()), options=list(dom="tip"), rownames=FALSE))
    files <- list.files(rv$spocc$outdir, recursive=TRUE, full.names=TRUE)
    dat <- data.frame(File=gsub(paste0("^", normalizePath(rv$spocc$outdir, winslash="/"), "/?"), "", normalizePath(files, winslash="/", mustWork=FALSE)),
                      Size_KB=round(file.info(files)$size/1024,1))
    datatable(dat, options=list(pageLength=12), rownames=FALSE)
  })
  output$spocc_download <- downloadHandler(
    filename = function() if (is.null(rv$spocc$zip)) paste0("spOccupancy_results_", timestamp_id(), ".zip") else basename(rv$spocc$zip),
    content = function(file) {
      copy_zip_to_download(rv$spocc$zip, rv$spocc$outdir, "spOccupancy", file, input$project_name)
    }
  )



  # sjSDM check
  refresh_sjsdm_check <- function() {
    rv$sjsdm$Y <- read_csv_safe(input$sjsdm_Y_file$datapath)
    rv$sjsdm$env <- read_csv_safe(input$sjsdm_env_file$datapath)
    rv$sjsdm$spatial <- read_csv_safe(input$sjsdm_spatial_file$datapath)
    rv$sjsdm$traits <- read_csv_safe(input$sjsdm_traits_file$datapath)
    rv$sjsdm$newdata <- read_csv_safe(input$sjsdm_newdata_file$datapath)
    rv$sjsdm$check <- validate_sjsdm(rv$sjsdm$Y, rv$sjsdm$env, rv$sjsdm$spatial, rv$sjsdm$traits, rv$sjsdm$newdata,
      input$sjsdm_family, input$sjsdm_env_model, input$sjsdm_spatial_model, input$sjsdm_biotic_lambda, input$sjsdm_biotic_alpha,
      env_formula = input$sjsdm_env_formula, spatial_formula = input$sjsdm_spatial_formula,
      iter = input$sjsdm_iter, sampling = input$sjsdm_sampling, learning_rate = input$sjsdm_learning_rate,
      step_size = input$sjsdm_step_size, parallel = input$sjsdm_parallel, cv_k = input$sjsdm_cv_k,
      dnn_hidden = input$sjsdm_dnn_hidden, dropout = input$sjsdm_dropout, device = input$sjsdm_device,
      dtype = input$sjsdm_dtype, verbose = input$sjsdm_verbose, optimizer = input$sjsdm_optimizer,
      weight_decay = input$sjsdm_weight_decay, scheduler = input$sjsdm_scheduler,
      lr_reduce_factor = input$sjsdm_lr_reduce_factor, early_stopping_training = input$sjsdm_early_stopping,
      mixed = input$sjsdm_mixed, generate_spatial_ev = input$sjsdm_generate_spatial_ev,
      spatial_ev_threshold = input$sjsdm_spatial_ev_threshold,
      anova_samples = input$sjsdm_anova_samples, tune_steps = input$sjsdm_tune_steps)
    rv$sjsdm$checked <- TRUE
  }

  observeEvent(input$sjsdm_check, {
    refresh_sjsdm_check()
    add_log("sjsdm", "sjSDM data and settings check completed.")
  })
  output$sjsdm_check_messages <- renderText({ if (!isTRUE(rv$sjsdm$checked)) return("sjSDM data have not been checked yet."); paste(rv$sjsdm$check$messages, collapse="\n") })
  output$sjsdm_data_table <- renderDT({
    extra <- list(
      data.frame(File="spatial.csv", Required="Optional / spatial models", Role="Spatial coordinates, eigenvectors or predictors", Rows=if(is.null(rv$sjsdm$spatial)) NA else nrow(rv$sjsdm$spatial), Columns=if(is.null(rv$sjsdm$spatial)) NA else ncol(rv$sjsdm$spatial), Status=if(is.null(rv$sjsdm$spatial)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="traits.csv", Required="Optional", Role="Species traits / metadata for post-hoc interpretation", Rows=if(is.null(rv$sjsdm$traits)) NA else nrow(rv$sjsdm$traits), Columns=if(is.null(rv$sjsdm$traits)) NA else ncol(rv$sjsdm$traits), Status=if(is.null(rv$sjsdm$traits)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="newdata.csv", Required="Optional", Role="Prediction environmental covariates", Rows=if(is.null(rv$sjsdm$newdata)) NA else nrow(rv$sjsdm$newdata), Columns=if(is.null(rv$sjsdm$newdata)) NA else ncol(rv$sjsdm$newdata), Status=if(is.null(rv$sjsdm$newdata)) "Not uploaded or unreadable" else "Readable")
    )
    datatable(data_summary(rv$sjsdm$Y, rv$sjsdm$env, NULL, NULL, NULL, extra), options=list(dom="tip", pageLength=10), rownames=FALSE)
  })
  sjsdm_config <- reactive({
    list(project_name=input$project_name, project_question=input$project_question, question_template=input$question_template,
      project_design=project_profile(),
      engine="sjSDM",
      data=list(Y=if(!is.null(input$sjsdm_Y_file$name)) input$sjsdm_Y_file$name else NULL, env=if(!is.null(input$sjsdm_env_file$name)) input$sjsdm_env_file$name else NULL, spatial=if(!is.null(input$sjsdm_spatial_file$name)) input$sjsdm_spatial_file$name else NULL, traits=if(!is.null(input$sjsdm_traits_file$name)) input$sjsdm_traits_file$name else NULL, newdata=if(!is.null(input$sjsdm_newdata_file$name)) input$sjsdm_newdata_file$name else NULL, new_spatial=if(!is.null(input$sjsdm_spatial_new_file$name)) input$sjsdm_spatial_new_file$name else NULL, species_groups=if(!is.null(input$sjsdm_group_file$name)) input$sjsdm_group_file$name else NULL, pretrained_weights=if(!is.null(input$sjsdm_weights_file$name)) input$sjsdm_weights_file$name else NULL, folds=if(!is.null(input$sjsdm_cv_file$name)) input$sjsdm_cv_file$name else NULL),
      model=list(family=input$sjsdm_family, env_model=input$sjsdm_env_model, env_formula=input$sjsdm_env_formula, spatial_model=input$sjsdm_spatial_model, spatial_formula=input$sjsdm_spatial_formula, se=input$sjsdm_se, iter=input$sjsdm_iter, step_size=input$sjsdm_step_size, sampling=input$sjsdm_sampling, parallel=input$sjsdm_parallel, dtype=input$sjsdm_dtype, verbose=input$sjsdm_verbose, seed=input$sjsdm_seed),
      regularization_biotic=list(env_lambda=input$sjsdm_env_lambda, env_alpha=input$sjsdm_env_alpha, spatial_lambda=input$sjsdm_spatial_lambda, spatial_alpha=input$sjsdm_spatial_alpha, biotic_lambda=input$sjsdm_biotic_lambda, biotic_alpha=input$sjsdm_biotic_alpha, biotic_df=input$sjsdm_biotic_df, on_diag=input$sjsdm_on_diag, reg_on_Cov=input$sjsdm_reg_on_Cov, inverse=input$sjsdm_inverse, tune_regularization=input$sjsdm_tune_regularization, tune_steps=input$sjsdm_tune_steps, cv_k=input$sjsdm_cv_k),
      dnn_optimizer=list(hidden=input$sjsdm_dnn_hidden, activation=input$sjsdm_activation, dropout=input$sjsdm_dropout, bias=input$sjsdm_bias, optimizer=input$sjsdm_optimizer, learning_rate=input$sjsdm_learning_rate, weight_decay=input$sjsdm_weight_decay, device=input$sjsdm_device),
      control=list(scheduler=input$sjsdm_scheduler, lr_reduce_factor=input$sjsdm_lr_reduce_factor, early_stopping_training=input$sjsdm_early_stopping, mixed=input$sjsdm_mixed),
      spatial_anova_metacommunity=list(generateSpatialEV=input$sjsdm_generate_spatial_ev, spatial_ev_k=input$sjsdm_spatial_ev_k, spatial_ev_threshold=input$sjsdm_spatial_ev_threshold, include_space_in_anova=input$sjsdm_include_space_in_anova, do_anova=input$sjsdm_do_anova, anova_samples=input$sjsdm_anova_samples, do_internal=input$sjsdm_do_internal, internal_fractions=input$sjsdm_internal_fractions, do_assembly=input$sjsdm_do_assembly, assembly_predictor=input$sjsdm_assembly_predictor),
      outputs=list(real_fit=input$sjsdm_real_fit, predict=input$sjsdm_do_predict, Rsquared=input$sjsdm_do_rsquared, importance=input$sjsdm_do_importance, weights=input$sjsdm_do_weights, coef=input$sjsdm_do_coef, covariance_correlation=input$sjsdm_do_cov_cor, residuals=input$sjsdm_do_residuals, plots=input$sjsdm_do_plots, zip=input$sjsdm_out_zip, copy_inputs=input$sjsdm_out_copy_inputs, save_config=input$sjsdm_out_config, report=input$sjsdm_out_report))
  })
  run_sjsdm_safe <- function() {
    cfg <- sjsdm_config()
    outdir <- make_engine_run_dir("sjSDM", input$project_name)
    rv$sjsdm$outdir <- outdir
    rv$sjsdm$status <- "Running"
    add_log("sjsdm", "Created output folder:", outdir)
    if (isTRUE(input$sjsdm_out_copy_inputs) || isTRUE(cfg$outputs$real_fit)) {
      copy_upload(input$sjsdm_Y_file,outdir,"Y.csv"); copy_upload(input$sjsdm_env_file,outdir,"env.csv"); copy_upload(input$sjsdm_spatial_file,outdir,"spatial.csv"); copy_upload(input$sjsdm_traits_file,outdir,"traits.csv"); copy_upload(input$sjsdm_newdata_file,outdir,"newdata.csv"); copy_upload(input$sjsdm_spatial_new_file,outdir,"new_spatial.csv"); copy_upload(input$sjsdm_group_file,outdir,"species_groups.csv"); copy_upload(input$sjsdm_weights_file,outdir,"pretrained_weights.rds"); copy_upload(input$sjsdm_cv_file,outdir,"folds.csv")
    }
    if (isTRUE(input$sjsdm_out_config) || isTRUE(cfg$outputs$real_fit)) yaml::write_yaml(cfg, file.path(outdir,"used_config.yml"))
    write.csv(data_summary(rv$sjsdm$Y, rv$sjsdm$env), file.path(outdir,"tables","data_summary.csv"), row.names=FALSE)
    write_data_check_messages(outdir, if (!is.null(rv$sjsdm$check)) rv$sjsdm$check$messages else character())
    status <- list(engine="sjSDM", status="ready", runtime_seconds=0, warnings=character(), errors=character())
    if (cfg$dnn_optimizer$device %in% c("gpu")) status$warnings <- c(status$warnings, "sjSDM warning: GPU selected. sjSDM 1.0.7 expects device='gpu' or a numeric GPU id, and PyTorch CUDA must be available before fitting.")
    if (cfg$regularization_biotic$tune_regularization && cfg$regularization_biotic$cv_k < 2) status$warnings <- c(status$warnings, "sjSDM warning: regularization tuning needs CV folds >= 2.")
    real_sjsdm_fit <- isTRUE(input$sjsdm_real_fit)
    if (!safe_require("sjSDM")) {
      status$status <- "check_failed"
      status$warnings <- c(status$warnings, "sjSDM package is not available. A diagnostic output scaffold was saved; no model fit was attempted.")
    } else if (isTRUE(real_sjsdm_fit)) {
      status$status <- "model_defined"
      status$warnings <- c(status$warnings, "sjSDM package is available. The GUI will run reproducible_script/run_this_sjSDM_analysis.R for real fitting.")
    } else {
      status$status <- "model_defined"
      status$warnings <- c(status$warnings, "sjSDM package is available, but real fitting was disabled by the user. A reproducible scaffold and executable script were saved.")
    }
    status <- normalize_engine_status(status)
    sjsdm_files <- data.frame(Section=c("Configuration","Inputs","Model object","Coefficients","Associations","Standard errors","Predictions","R-squared","ANOVA / variation partitioning","Internal structure","Importance","Weights","Residuals","Plots","Diagnostics","Report"), Expected_file_or_folder=c("used_config.yml","inputs/Y.csv, env.csv, spatial.csv, traits.csv, newdata.csv","models/sjSDM_model.rds","tables/coef_environment.csv, coef_spatial.csv","tables/covariance_matrix.csv, correlation_matrix.csv","tables/standard_errors.csv, p_values.csv","predictions/predictions.csv or results/predictions.rds","tables/Rsquared_total.csv, Rsquared_species.csv, Rsquared_sites.csv","anova/sjSDM_anova_results.csv, anova_species.csv, anova_sites.csv","internal_structure/internal_structure_species.csv, internal_structure_sites.csv, assembly_effects.csv","importance/importance_summary.csv","weights/env_weights.rds, spatial_weights.rds, model_weights.rds","tables/residuals.csv","plots/sjSDM_plot.pdf, anova_plot.pdf, importance_plot.pdf, internal_structure_plot.pdf, assembly_effects.pdf","diagnostics/engine_status.json, data_check_messages.csv, torch_diagnostic.txt","report/sjSDM_report.html"), Meaning=c("Exact sjSDM settings used for the run.","Copied response, environmental, spatial, trait, prediction and grouping inputs.","Serialized fitted sjSDM object when production fitting is connected.","Environmental and spatial coefficient matrices from coef.sjSDM.","Species covariance/correlation matrices from getCov/getCor.","Post-hoc standard errors and p-values from getSe or se=TRUE workflows.","Predicted responses from predict.sjSDM.","Total, species-level and site-level R-squared / pseudo-R2 outputs.","Variation partitioning into environment, space and associations using anova.sjSDM.","Internal metacommunity structure and assembly-effect summaries.","Predictor importance from getImportance / importance.","DNN/model weights from getWeights for reproducibility.","Residual diagnostics from residuals.sjSDM.","Human-readable figures for model, ANOVA, importance and internal structure.","Machine-readable status, warnings and PyTorch/reticulate diagnostics.","HTML report for browser viewing."))
    write.csv(sjsdm_files, file.path(outdir,"tables","sjSDM_result_workflow_map.csv"), row.names=FALSE)
    sjsdm_api_map <- data.frame(
      GUI_setting = c("family", "env module", "env formula", "spatial module", "spatial formula", "biotic lambda/alpha/df", "iter", "step_size", "sampling", "parallel", "learning_rate", "optimizer + weight_decay", "scheduler / early stopping / mixed", "device", "dtype", "se", "generateSpatialEV", "ANOVA / internalStructure", "weights"),
      sjSDM_1_0_7_API = c("family = binomial('probit'/'logit'), poisson('log'), 'nbinom', or gaussian('identity')", "env = linear(...) or DNN(...); intercept-only creates linear(data, ~1)", "formula argument inside linear/DNN", "spatial = NULL, linear(...), or DNN(...)", "formula argument inside spatial linear/DNN; usually ~ 0 + .", "biotic = bioticStruct(df, lambda, alpha, on_diag, reg_on_Cov, inverse)", "sjSDM(iter = ...)", "sjSDM(step_size = ...), internally used as batch_size", "sjSDM(sampling = ...)", "sjSDM(parallel = ...)", "sjSDM(learning_rate = ...)", "sjSDMControl(optimizer = Adamax/RMSprop/SGD/AccSGD/AdaBound/madgrad(...))", "sjSDMControl(scheduler, lr_reduce_factor, early_stopping_training, mixed)", "sjSDM(device = 'cpu' or 'gpu')", "sjSDM(dtype = 'float32' or 'float64')", "sjSDM(se = TRUE/FALSE), getSe for post-hoc SE", "generateSpatialEV(coords, threshold) then spatial = linear(SPV[,1:k], ~0+.)", "anova(model), internalStructure(anova)", "getWeights(model), setWeights(model, weights)"),
      Audit_note = c("negative_binomial_log was corrected to nbinom.", "The previous 'none' option was corrected to intercept-only because sjSDM requires env.", "Character predictors are converted to factor and model.matrix encodes them.", "The value 'eigenvectors' is treated as a linear spatial module.", "Spatial intercept should usually be removed.", "alpha=0 is lasso, alpha=1 is ridge in sjSDM.", "True sjSDM argument.", "Replaces the previous misleading batch_size control.", "True sjSDM argument.", "Newly exposed; 0 is Windows-safe.", "True sjSDM argument.", "DiffGrad removed because it is not exported in sjSDM 1.0.7.", "Newly exposed.", "cuda label removed; sjSDM uses gpu.", "Newly exposed.", "True sjSDM argument.", "Threshold control added.", "Uses exported sjSDM functions.", "Uses exported sjSDM functions."),
      stringsAsFactors = FALSE
    )
    write.csv(sjsdm_api_map, file.path(outdir, "tables", "sjSDM_api_mapping.csv"), row.names = FALSE)
    write_engine_status(outdir, status)
    write_standard_outputs(outdir, "sjSDM", status, rv$sjsdm$Y, rv$sjsdm$env)
    write_reproducible_stub(outdir, "sjSDM")
    write_sjsdm_reproducible_script(outdir)
    write_engine_scaffold_outputs(outdir, "sjSDM", cfg, status, rv$sjsdm$Y, rv$sjsdm$env)
    writeLines(c("sjSDM result workflow","=====================","This SAFE build creates the output scaffold and configuration.","When the production sjSDM adapter is connected, it should write the files listed in tables/sjSDM_result_workflow_map.csv.","",paste(sjsdm_files$Section, sjsdm_files$Expected_file_or_folder, sep=" -> ")), file.path(outdir,"results","README_sjSDM_results.txt"))
    if (isTRUE(cfg$outputs$report)) make_html_report(outdir, "sjSDM", cfg, paste(status$status, paste(status$warnings, collapse="; ")))
    if (identical(status$status, "model_defined")) {
      writeLines("MODEL DEFINED - NO FITTED POSTERIOR", file.path(outdir, "RUN_MODEL_DEFINED.txt"))
    } else if (identical(status$status, "check_failed")) {
      writeLines("CHECK FAILED - NO FIT WAS ATTEMPTED", file.path(outdir, "RUN_CHECK_FAILED.txt"))
    } else if (identical(status$status, "fit_failed")) {
      writeLines("RUN FAILED", file.path(outdir, "RUN_FAILED.txt"))
    }
    if (isTRUE(real_sjsdm_fit) && safe_require("sjSDM")) {
      script_file <- file.path(outdir, "reproducible_script", "run_this_sjSDM_analysis.R")
      rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
      add_log("sjsdm", "Running real sjSDM reproducible script:", script_file)
      run_log <- tryCatch({
        neutralize_windows_r_locale_env()
        system2(rscript, shQuote(script_file), stdout = TRUE, stderr = TRUE)
      },
        error = function(e) e
      )
      if (inherits(run_log, "error")) {
        status$status <- "fit_failed"
        status$errors <- c(status$errors, conditionMessage(run_log))
        writeLines(conditionMessage(run_log), file.path(outdir, "diagnostics", "sjSDM_real_fit_error.txt"))
        writeLines("REAL FIT FAILED", file.path(outdir, "RUN_REAL_FIT_FAILED.txt"))
      } else {
        writeLines(as.character(run_log), file.path(outdir, "diagnostics", "sjSDM_real_fit_stdout_stderr.txt"))
        exit_status <- attr(run_log, "status") %||% 0
        status_file <- file.path(outdir, "diagnostics", "engine_status.json")
        script_payload <- tryCatch(jsonlite::fromJSON(status_file), error = function(e) NULL)
        script_status <- script_payload$status %||% NA_character_
        if (!is.null(script_payload$warnings) && length(script_payload$warnings) > 0) {
          status$warnings <- c(status$warnings, as.character(script_payload$warnings))
        }
        if (!is.null(script_payload$errors) && length(script_payload$errors) > 0) {
          status$errors <- c(status$errors, as.character(script_payload$errors))
        }
        if (identical(as.integer(exit_status), 0L) && file.exists(file.path(outdir, "models", "sjSDM_model.rds")) && identical(script_status, "fitted")) {
          status$status <- "fitted"
          status$warnings <- c(status$warnings, paste("Real sjSDM fit completed. See diagnostics/sjSDM_real_fit_stdout_stderr.txt."))
          if (file.exists(file.path(outdir, "RUN_MODEL_DEFINED.txt"))) unlink(file.path(outdir, "RUN_MODEL_DEFINED.txt"))
          writeLines("REAL FIT COMPLETE", file.path(outdir, "RUN_REAL_FIT_COMPLETE.txt"))
        } else {
          status$status <- "fit_failed"
          status$errors <- c(status$errors, paste("Real sjSDM script failed or did not save models/sjSDM_model.rds. Exit status:", exit_status))
          writeLines("REAL FIT FAILED", file.path(outdir, "RUN_REAL_FIT_FAILED.txt"))
        }
      }
      status <- normalize_engine_status(status)
      write_engine_status(outdir, status)
      if (isTRUE(cfg$outputs$report)) make_html_report(outdir, "sjSDM", cfg, paste(status$status, paste(c(status$warnings, status$errors), collapse="; ")))
    }
    status <- normalize_engine_status(status)
    ensure_output_contract(outdir, "sjSDM", status, rv$sjsdm$Y, rv$sjsdm$env)
    rv$sjsdm$zip <- if (isTRUE(input$sjsdm_out_zip)) make_zip(outdir) else NULL
    rv$sjsdm$status <- status$status
    add_log("sjsdm", "sjSDM workflow status:", status$status, outdir)
  }
  observeEvent(input$sjsdm_run, {
    refresh_sjsdm_check()
    if (!isTRUE(rv$sjsdm$check$ok)) {
      diag <- create_check_failed_zip("sjSDM", input$project_name, sjsdm_config(), rv$sjsdm$check,
        files = list("Y.csv" = input$sjsdm_Y_file, "env.csv" = input$sjsdm_env_file, "spatial.csv" = input$sjsdm_spatial_file,
                     "traits.csv" = input$sjsdm_traits_file, "newdata.csv" = input$sjsdm_newdata_file,
                     "new_spatial.csv" = input$sjsdm_spatial_new_file, "species_groups.csv" = input$sjsdm_group_file,
                     "folds.csv" = input$sjsdm_cv_file),
        Y = rv$sjsdm$Y, X = rv$sjsdm$env)
      rv$sjsdm$outdir <- diag$outdir
      rv$sjsdm$zip <- diag$zip
      rv$sjsdm$status <- "check_failed"
      add_log("sjsdm", "sjSDM check failed; diagnostic ZIP created:", diag$outdir)
      return(NULL)
    }
    run_sjsdm_safe()
  })
  output$sjsdm_log <- renderText(paste(rv$sjsdm$log, collapse="\n"))
  output$sjsdm_run_summary <- renderText({ paste0("Status: ", rv$sjsdm$status, "\nOutput folder: ", rv$sjsdm$outdir %||% "None yet", "\nZIP: ", rv$sjsdm$zip %||% "None yet") })
  output$sjsdm_files <- renderDT({ if (is.null(rv$sjsdm$outdir) || !dir.exists(rv$sjsdm$outdir)) return(datatable(data.frame(File=character(), Size_KB=numeric()), options=list(dom="tip"), rownames=FALSE)); files <- list.files(rv$sjsdm$outdir, recursive=TRUE, full.names=TRUE); dat <- data.frame(File=gsub(paste0("^", normalizePath(rv$sjsdm$outdir, winslash="/"), "/?"), "", normalizePath(files, winslash="/", mustWork=FALSE)), Size_KB=round(file.info(files)$size/1024,1)); datatable(dat, options=list(pageLength=12), rownames=FALSE) })
  output$sjsdm_download <- downloadHandler(
    filename = function() if (is.null(rv$sjsdm$zip)) paste0("sjSDM_results_", timestamp_id(), ".zip") else basename(rv$sjsdm$zip),
    content = function(file) copy_zip_to_download(rv$sjsdm$zip, rv$sjsdm$outdir, "sjSDM", file, input$project_name)
  )



  # boral check
  refresh_boral_check <- function() {
    rv$boral$Y <- read_csv_safe(input$boral_Y_file$datapath)
    rv$boral$X <- read_csv_safe(input$boral_X_file$datapath)
    rv$boral$traits <- read_csv_safe(input$boral_traits_file$datapath)
    rv$boral$rowids <- read_csv_safe(input$boral_rowids_file$datapath, id_columns = c("", "row.names", "row_id", "site_id"))
    rv$boral$ranefids <- read_csv_safe(input$boral_ranefids_file$datapath, id_columns = c("", "row.names", "row_id", "site_id"))
    rv$boral$distmat <- read_csv_safe(input$boral_distmat_file$datapath)
    rv$boral$offset <- read_csv_safe(input$boral_offset_file$datapath)
    rv$boral$newdata <- read_csv_safe(input$boral_newdata_file$datapath)
    rv$boral$check <- validate_boral(
      rv$boral$Y, rv$boral$X, rv$boral$traits, rv$boral$rowids, rv$boral$ranefids,
      rv$boral$distmat, rv$boral$offset, input$boral_family, input$boral_family_text,
      input$boral_num_lv, input$boral_lv_type, input$boral_row_eff,
      input$boral_n_burnin, input$boral_n_iteration, input$boral_n_thin,
      input$boral_trial_size, input$boral_calc_ics, input$boral_use_traits,
      input$boral_which_traits, input$boral_use_ssvs, input$boral_ssvs_index,
      input$boral_save_model, input$boral_do_fit, input$boral_formula_X
    )
    rv$boral$checked <- TRUE
  }

  observeEvent(input$boral_check, {
    refresh_boral_check()
    add_log("boral", "boral data and settings check completed.")
  })

  output$boral_check_messages <- renderText({
    if (!isTRUE(rv$boral$checked)) return("boral data have not been checked yet.")
    paste(rv$boral$check$messages, collapse="\n")
  })

  output$boral_data_table <- renderDT({
    extra <- list(
      data.frame(File="traits.csv", Required="Optional", Role="Species traits for fourth-corner style model", Rows=if(is.null(rv$boral$traits)) NA else nrow(rv$boral$traits), Columns=if(is.null(rv$boral$traits)) NA else ncol(rv$boral$traits), Status=if(is.null(rv$boral$traits)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="row.ids.csv", Required="Optional", Role="Row-effect grouping IDs", Rows=if(is.null(rv$boral$rowids)) NA else nrow(rv$boral$rowids), Columns=if(is.null(rv$boral$rowids)) NA else ncol(rv$boral$rowids), Status=if(is.null(rv$boral$rowids)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="ranef.ids.csv", Required="Optional", Role="Response-specific random intercept IDs", Rows=if(is.null(rv$boral$ranefids)) NA else nrow(rv$boral$ranefids), Columns=if(is.null(rv$boral$ranefids)) NA else ncol(rv$boral$ranefids), Status=if(is.null(rv$boral$ranefids)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="distmat.csv", Required="Spatial latent variables", Role="Distance matrix for structured latent variables", Rows=if(is.null(rv$boral$distmat)) NA else nrow(rv$boral$distmat), Columns=if(is.null(rv$boral$distmat)) NA else ncol(rv$boral$distmat), Status=if(is.null(rv$boral$distmat)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="offset.csv", Required="Optional", Role="Offset matrix", Rows=if(is.null(rv$boral$offset)) NA else nrow(rv$boral$offset), Columns=if(is.null(rv$boral$offset)) NA else ncol(rv$boral$offset), Status=if(is.null(rv$boral$offset)) "Not uploaded or unreadable" else "Readable"),
      data.frame(File="newdata.csv", Required="Optional", Role="Prediction covariates", Rows=if(is.null(rv$boral$newdata)) NA else nrow(rv$boral$newdata), Columns=if(is.null(rv$boral$newdata)) NA else ncol(rv$boral$newdata), Status=if(is.null(rv$boral$newdata)) "Not uploaded or unreadable" else "Readable")
    )
    datatable(data_summary(rv$boral$Y, rv$boral$X, NULL, NULL, NULL, extra),
              options = list(dom="tip", pageLength=10), rownames=FALSE)
  })

  boral_config <- reactive({
    list(
      project_name = input$project_name,
      project_question = input$project_question,
      question_template = input$question_template,
      project_design = project_profile(),
      engine = "boral",
      data = list(
        Y = if(!is.null(input$boral_Y_file$name)) input$boral_Y_file$name else NULL,
        XData = if(!is.null(input$boral_X_file$name)) input$boral_X_file$name else NULL,
        traits = if(!is.null(input$boral_traits_file$name)) input$boral_traits_file$name else NULL,
        row.ids = if(!is.null(input$boral_rowids_file$name)) input$boral_rowids_file$name else NULL,
        ranef.ids = if(!is.null(input$boral_ranefids_file$name)) input$boral_ranefids_file$name else NULL,
        distmat = if(!is.null(input$boral_distmat_file$name)) input$boral_distmat_file$name else NULL,
        offset = if(!is.null(input$boral_offset_file$name)) input$boral_offset_file$name else NULL,
        newdata = if(!is.null(input$boral_newdata_file$name)) input$boral_newdata_file$name else NULL,
        trial.size.file = if(!is.null(input$boral_trials_file$name)) input$boral_trials_file$name else NULL
      ),
      model = list(
        model_mode = input$boral_model_mode,
        family = input$boral_family,
        family_vector_override = input$boral_family_text,
        lv.control = list(num.lv = input$boral_num_lv, type = input$boral_lv_type),
        model.name = input$boral_model_name,
        formula.X = input$boral_formula_X,
        X.ind = input$boral_X_ind,
        trial.size = input$boral_trial_size,
        row.eff = input$boral_row_eff,
        use_offset = input$boral_use_offset,
        do.fit = input$boral_do_fit
      ),
      traits_random_ssvs = list(
        use_traits = input$boral_use_traits,
        which.traits = input$boral_which_traits,
        traits_no_intercept = input$boral_traits_intercept_warning,
        use_ranef = input$boral_use_ranef,
        use_ssvs = input$boral_use_ssvs,
        ssvs.index = input$boral_ssvs_index,
        ssvs.traitsindex = input$boral_ssvs_traitsindex,
        ssvs.g = input$boral_ssvs_g,
        save.model = input$boral_save_model
      ),
      mcmc_prior = list(
        n.burnin = input$boral_n_burnin,
        n.iteration = input$boral_n_iteration,
        n.thin = input$boral_n_thin,
        seed = input$boral_seed,
        prior.type = input$boral_prior_type,
        hypparams = input$boral_hypparams,
        calc.ics = input$boral_calc_ics,
        save_mcmc_samples = input$boral_save_mcmc_samples,
        save_hpd = input$boral_save_hpd,
        save_dic = input$boral_save_dic
      ),
      diagnostics_outputs = list(
        summary = input$boral_do_summary,
        residual_plot = input$boral_do_residual_plot,
        lvsplot = input$boral_do_lvsplot,
        ind.spp = input$boral_ind_spp,
        ranefsplot = input$boral_do_ranefsplot,
        coefsplot = input$boral_do_coefsplot,
        enviro_cor = input$boral_do_env_cor,
        residual_cor = input$boral_do_resid_cor,
        varpart = input$boral_do_varpart,
        predict = input$boral_do_predict,
        fitted = input$boral_do_fitted,
        tidyboral = input$boral_do_tidy
      ),
      outputs = list(
        save_model = input$boral_out_model,
        save_jags = input$boral_out_jags,
        save_tables = input$boral_out_tables,
        save_plots = input$boral_out_plots,
        report = input$boral_out_report,
        zip = input$boral_out_zip,
        copy_inputs = TRUE,
        save_config = TRUE
      )
    )
  })

  run_boral_safe <- function() {
    cfg <- boral_config()
    outdir <- make_engine_run_dir("boral", input$project_name)
    rv$boral$outdir <- outdir
    rv$boral$status <- "Running"
    add_log("boral", "Created output folder:", outdir)

    {
      copy_upload(input$boral_Y_file, outdir, "Y.csv")
      copy_upload(input$boral_X_file, outdir, "XData.csv")
      copy_upload(input$boral_traits_file, outdir, "traits.csv")
      copy_upload(input$boral_rowids_file, outdir, "row.ids.csv")
      copy_upload(input$boral_ranefids_file, outdir, "ranef.ids.csv")
      copy_upload(input$boral_distmat_file, outdir, "distmat.csv")
      copy_upload(input$boral_offset_file, outdir, "offset.csv")
      copy_upload(input$boral_newdata_file, outdir, "newdata.csv")
      copy_upload(input$boral_trials_file, outdir, "trial.size.csv")
    }
    yaml::write_yaml(cfg, file.path(outdir, "used_config.yml"))
    write.csv(data_summary(rv$boral$Y, rv$boral$X), file.path(outdir, "tables", "data_summary.csv"), row.names=FALSE)
    write_data_check_messages(outdir, if (!is.null(rv$boral$check)) rv$boral$check$messages else character())

    if (exists("write_boral_result_workflow_map", mode = "function")) write_boral_result_workflow_map(outdir)
    if (exists("write_boral_reproducible_script", mode = "function")) {
      script_file <- write_boral_reproducible_script(outdir)
    } else {
      script_file <- file.path(outdir, "reproducible_script", "run_this_boral_analysis.R")
      write_reproducible_stub(outdir, "boral")
      status <- list(engine="boral", status="fit_failed", runtime_seconds=0, warnings=character(),
                     errors="Internal error: R/boral_adapter.R was not loaded, so real boral fitting could not be attempted.")
      write_engine_status(outdir, status)
      write_standard_outputs(outdir, "boral", status, rv$boral$Y, rv$boral$X)
    }

    stdout_file <- file.path(outdir, "diagnostics", "boral_Rscript_stdout.log")
    stderr_file <- file.path(outdir, "diagnostics", "boral_Rscript_stderr.log")
    rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    exit_code <- tryCatch({
      neutralize_windows_r_locale_env()
      system2(rscript, args = shQuote(script_file), stdout = stdout_file, stderr = stderr_file)
    },
      error = function(e) {
        writeLines(conditionMessage(e), stderr_file)
        1L
      }
    )
    if (file.exists(stdout_file)) {
      sout <- readLines(stdout_file, warn = FALSE)
      if (length(sout) > 0) add_log("boral", paste(tail(sout, 20), collapse = "\n"))
    }
    if (file.exists(stderr_file)) {
      serr <- readLines(stderr_file, warn = FALSE)
      if (length(serr) > 0) add_log("boral", paste(tail(serr, 20), collapse = "\n"))
    }

    status_file <- file.path(outdir, "diagnostics", "engine_status.json")
    status <- NULL
    if (file.exists(status_file) && requireNamespace("jsonlite", quietly = TRUE)) {
      status <- tryCatch(jsonlite::fromJSON(status_file, simplifyVector = FALSE), error = function(e) NULL)
    }
    if (is.null(status)) {
      status <- list(engine="boral", status="fit_failed",
                     runtime_seconds=0, warnings=character(),
                     errors=paste0("boral Rscript exited with status ", exit_code, " without writing diagnostics/engine_status.json."))
      write_engine_status(outdir, status)
      write_standard_outputs(outdir, "boral", status, rv$boral$Y, rv$boral$X)
      write_engine_scaffold_outputs(outdir, "boral", cfg, status, rv$boral$Y, rv$boral$X)
    } else if (!identical(exit_code, 0L) && !(status$status %in% c("fit_failed", "check_failed"))) {
      status$status <- "fit_failed"
      status$errors <- c(status$errors %||% character(), paste0("boral Rscript exited with status ", exit_code, "."))
      write_engine_status(outdir, status)
      write_standard_outputs(outdir, "boral", status, rv$boral$Y, rv$boral$X)
    }
    status <- normalize_engine_status(status)
    write_engine_status(outdir, status)
    ensure_output_contract(outdir, "boral", status, rv$boral$Y, rv$boral$X)
    writeLines(c(
      "boral result workflow",
      "=====================",
      paste0("Status: ", status$status %||% "unknown"),
      "This folder contains copied inputs, the exact used_config.yml, executable reproducible_script/run_this_boral_analysis.R, diagnostics, standard comparison tables and any boral outputs created by the run.",
      "If status is fit_failed, start with diagnostics/boral_dependency_error.txt, diagnostics/boral_fit_error.txt, diagnostics/JAGS_status.txt and diagnostics/engine_status.json."
    ), file.path(outdir, "results", "README_boral_results.txt"))
    if (identical(status$status, "fitted")) {
      writeLines("RUN COMPLETE", file.path(outdir, "RUN_COMPLETE.txt"))
    } else if (identical(status$status, "model_defined")) {
      writeLines("MODEL DEFINED - NO FITTED POSTERIOR", file.path(outdir, "RUN_MODEL_DEFINED.txt"))
    } else if (identical(status$status, "check_failed")) {
      writeLines("CHECK FAILED - NO FIT WAS ATTEMPTED", file.path(outdir, "RUN_CHECK_FAILED.txt"))
    } else {
      writeLines("RUN FAILED", file.path(outdir, "RUN_FAILED.txt"))
    }
    rv$boral$zip <- if (isTRUE(input$boral_out_zip)) make_zip(outdir) else NULL
    rv$boral$status <- status$status
    add_log("boral", "boral workflow status:", status$status, outdir)
  }

  observeEvent(input$boral_run, {
    refresh_boral_check()
    if (!isTRUE(rv$boral$check$ok)) {
      diag <- create_check_failed_zip("boral", input$project_name, boral_config(), rv$boral$check,
        files = list("Y.csv" = input$boral_Y_file, "XData.csv" = input$boral_X_file, "traits.csv" = input$boral_traits_file,
                     "row.ids.csv" = input$boral_rowids_file, "ranef.ids.csv" = input$boral_ranefids_file,
                     "distmat.csv" = input$boral_distmat_file, "offset.csv" = input$boral_offset_file,
                     "newdata.csv" = input$boral_newdata_file, "trial.size.csv" = input$boral_trials_file),
        Y = rv$boral$Y, X = rv$boral$X)
      rv$boral$outdir <- diag$outdir
      rv$boral$zip <- diag$zip
      rv$boral$status <- "check_failed"
      add_log("boral", "boral check failed; diagnostic ZIP created:", diag$outdir)
      return(NULL)
    }
    run_boral_safe()
  })

  output$boral_log <- renderText(paste(rv$boral$log, collapse="\n"))
  output$boral_run_summary <- renderText({
    paste0("Status: ", rv$boral$status, "\nOutput folder: ", rv$boral$outdir %||% "None yet", "\nZIP: ", rv$boral$zip %||% "None yet")
  })
  output$boral_files <- renderDT({
    if (is.null(rv$boral$outdir) || !dir.exists(rv$boral$outdir)) return(datatable(data.frame(File=character(), Size_KB=numeric()), options=list(dom="tip"), rownames=FALSE))
    files <- list.files(rv$boral$outdir, recursive=TRUE, full.names=TRUE)
    dat <- data.frame(File=gsub(paste0("^", normalizePath(rv$boral$outdir, winslash="/"), "/?"), "", normalizePath(files, winslash="/", mustWork=FALSE)),
                      Size_KB=round(file.info(files)$size/1024,1))
    datatable(dat, options=list(pageLength=12), rownames=FALSE)
  })
  output$boral_download <- downloadHandler(
    filename = function() if (is.null(rv$boral$zip)) paste0("boral_results_", timestamp_id(), ".zip") else basename(rv$boral$zip),
    content = function(file) {
      copy_zip_to_download(rv$boral$zip, rv$boral$outdir, "boral", file, input$project_name)
    }
  )


  observeEvent(input$compare_run, {
    selected <- input$compare_engines
    if (!length(selected)) {
      rv$compare <- NULL
      showNotification("Select at least one engine.", type = "warning")
      return(invisible(NULL))
    }
    hdir <- trimws(input$compare_hmsc_dir %||% "")
    hpdir <- trimws(input$compare_hmschpc_dir %||% "")
    jdir <- trimws(input$compare_jsdm_dir %||% "")
    gdir <- trimws(input$compare_gjam_dir %||% "")
    sdir <- trimws(input$compare_spocc_dir %||% "")
    sjdir <- trimws(input$compare_sjsdm_dir %||% "")
    bdir <- trimws(input$compare_boral_dir %||% "")
    if (!nzchar(hdir) && !is.null(rv$hmsc$outdir)) hdir <- rv$hmsc$outdir
    if (!nzchar(hpdir) && !is.null(rv$hmschpc$outdir)) hpdir <- rv$hmschpc$outdir
    if (!nzchar(jdir) && !is.null(rv$jsdm$outdir)) jdir <- rv$jsdm$outdir
    if (!nzchar(gdir) && !is.null(rv$gjam$outdir)) gdir <- rv$gjam$outdir
    if (!nzchar(sdir) && !is.null(rv$spocc$outdir)) sdir <- rv$spocc$outdir
    if (!nzchar(sjdir) && !is.null(rv$sjsdm$outdir)) sjdir <- rv$sjsdm$outdir
    if (!nzchar(bdir) && !is.null(rv$boral$outdir)) bdir <- rv$boral$outdir

    add <- function(engine, dir) {
      exists <- nzchar(dir) && dir.exists(dir)
      folder <- if (exists) normalizePath(dir, winslash = "/", mustWork = FALSE) else dir
      safe_file <- function(path) {
        if (!exists || !length(path) || !nzchar(path[1])) return(FALSE)
        info <- file.info(path[1])
        isTRUE(!info$isdir && info$size > 0)
      }
      safe_count <- function(subdir) {
        path <- file.path(dir, subdir)
        if (isTRUE(exists) && dir.exists(path)) length(list.files(path)) else 0
      }
      status_file <- file.path(dir, "diagnostics", "engine_status.json")
      config_file <- file.path(dir, "used_config.yml")
      standard_run_file <- file.path(dir, "standard", "run_summary.csv")
      standard_fit_file <- file.path(dir, "standard", "fit_metrics.csv")
      standard_effects_file <- file.path(dir, "standard", "effects_long.csv")
      standard_predictions_file <- file.path(dir, "standard", "predictions_long.csv")
      standard_associations_file <- file.path(dir, "standard", "associations_long.csv")
      standard_diagnostics_file <- file.path(dir, "standard", "diagnostics_long.csv")
      standard_effects_species_file <- file.path(dir, "standard", "effects_species_environment.csv")
      standard_predictions_site_file <- file.path(dir, "standard", "predictions_site_species.csv")
      standard_associations_species_file <- file.path(dir, "standard", "associations_species_species.csv")
      standard_manifest_file <- file.path(dir, "standard", "output_manifest.csv")
      map_file <- switch(engine,
        Hmsc = file.path(dir, "tables", "Hmsc_result_workflow_map.csv"),
        "Hmsc-HPC" = file.path(dir, "tables", "Hmsc-HPC_result_workflow_map.csv"),
        jSDM = file.path(dir, "tables", "jSDM_result_workflow_map.csv"),
        GJAM = file.path(dir, "tables", "GJAM_result_workflow_map.csv"),
        spOccupancy = file.path(dir, "tables", "spOccupancy_result_workflow_map.csv"),
        sjSDM = file.path(dir, "tables", "sjSDM_result_workflow_map.csv"),
        boral = file.path(dir, "tables", "boral_result_workflow_map.csv")
      )
      workflow_index_file <- switch(engine,
        Hmsc = file.path(dir, "workflow_scripts", "S1_define_models.R"),
        "Hmsc-HPC" = file.path(dir, "tables", "HmscHPC_S1S7_result_index.csv"),
        jSDM = file.path(dir, "tables", "jSDM_result_workflow_map.csv"),
        GJAM = file.path(dir, "tables", "GJAM_result_workflow_map.csv"),
        spOccupancy = file.path(dir, "tables", "spOccupancy_result_workflow_map.csv"),
        sjSDM = file.path(dir, "tables", "sjSDM_result_workflow_map.csv"),
        boral = file.path(dir, "tables", "boral_result_workflow_map.csv")
      )
      report_file <- switch(engine,
        Hmsc = file.path(dir, "report", "Hmsc_report.html"),
        "Hmsc-HPC" = file.path(dir, "report", "Hmsc-HPC_report.html"),
        jSDM = file.path(dir, "report", "jSDM_report.html"),
        GJAM = file.path(dir, "report", "GJAM_report.html"),
        spOccupancy = file.path(dir, "report", "spOccupancy_report.html"),
        sjSDM = file.path(dir, "report", "sjSDM_report.html"),
        boral = file.path(dir, "report", "boral_report.html")
      )
      standard_run <- if (safe_file(standard_run_file)) tryCatch(read.csv(standard_run_file, check.names = FALSE), error = function(e) NULL) else NULL
      status_json <- if (safe_file(status_file)) {
        tryCatch(jsonlite::fromJSON(status_file, simplifyVector = FALSE), error = function(e) NULL)
      } else NULL
      if (!is.list(status_json)) status_json <- NULL
      json_status <- status_json$status
      if (!is.character(json_status) || length(json_status) != 1L || !nzchar(json_status)) json_status <- NULL
      csv_status <- standard_run$status[1]
      status_text_raw <- as.character(json_status %||% csv_status %||% "")
      recorded_engines <- c(status_json$engine, standard_run$engine[1])
      recorded_engines <- recorded_engines[!is.na(recorded_engines) & nzchar(recorded_engines)]
      conflict <- length(recorded_engines) && any(tolower(recorded_engines) != tolower(engine))
      disagree <- length(json_status) == 1L && length(csv_status) == 1L &&
        !is.na(csv_status) && normalize_engine_status(json_status)$status != normalize_engine_status(csv_status)$status
      status_text <- if (!exists) {
        if (nzchar(dir)) "missing_folder" else "not_run"
      } else if (conflict) {
        "engine_mismatch"
      } else if (disagree || (file.exists(status_file) && !length(json_status))) {
        "status_conflict"
      } else if (!nzchar(status_text_raw)) {
        "missing_status"
      } else {
        normalize_engine_status(status_text_raw)$status
      }
      status_class <- if (identical(status_text, "fitted")) {
        "fitted_result"
      } else if (identical(status_text, "model_defined")) {
        "model_defined_no_posterior"
      } else if (identical(status_text, "check_failed")) {
        "check_failed_before_fit"
      } else if (identical(status_text, "fit_failed")) {
        "fit_failed_after_start"
      } else {
        status_text
      }
      standard_complete <- all(c(
        safe_file(status_file),
        safe_file(config_file),
        safe_file(standard_run_file),
        safe_file(standard_fit_file),
        safe_file(standard_effects_file),
        safe_file(standard_predictions_file),
        safe_file(standard_associations_file),
        safe_file(standard_diagnostics_file),
        safe_file(standard_effects_species_file),
        safe_file(standard_predictions_site_file),
        safe_file(standard_associations_species_file),
        safe_file(standard_manifest_file)
      ))
      data_check_file <- file.path(dir, "diagnostics", "data_check_messages.csv")
      session_file <- file.path(dir, "diagnostics", "session_info.txt")
      error_files <- if (isTRUE(exists) && dir.exists(file.path(dir, "diagnostics"))) {
        list.files(file.path(dir, "diagnostics"), pattern = "(_error|error)\\.txt$", full.names = FALSE)
      } else character()
      fitted_status <- identical(status_text, "fitted")
      has_values <- function(path, field) {
        if (!fitted_status || !safe_file(path)) return(FALSE)
        dat <- standard_read_csv(path)
        any(is.finite(suppressWarnings(as.numeric(dat[[field]]))))
      }
      has_metrics <- has_values(standard_fit_file, "value")
      has_effects <- has_values(standard_effects_file, "estimate")
      has_predictions <- has_values(standard_predictions_file, "predicted_mean")
      has_associations <- has_values(standard_associations_file, "estimate")
      metric_note <- if (has_metrics) {
        "Metrics present; validation split, target, scale and metric definition were not verified."
      } else "No finite fit metrics from a verified-status run."
      association_note <- switch(engine,
        Hmsc = "Omega/random-level associations; compare only broad patterns with other engines.",
        "Hmsc-HPC" = "Random-level correlations are derived from loadings where available; raw Eta and Lambda are not association matrices.",
        jSDM = "Latent-variable residual correlations; not the same parameterization as Hmsc or boral.",
        GJAM = "corMu/sigMu on the GJAM observation/model scale; typeNames matter.",
        spOccupancy = "Selected functions estimate residual occupancy correlations; loadings alone are not species correlations.",
        sjSDM = "Covariance or correlation from the fitted biotic structure; check the selected parameterization.",
        boral = "Latent-variable residual correlations from JAGS model; check family and LV structure.",
        "Engine-specific association output."
      )
      failure_diag <- if (!exists) {
        ""
      } else if (identical(status_text, "check_failed")) {
        "diagnostics/data_check_messages.csv; used_config.yml"
      } else if (identical(status_text, "fit_failed")) {
        paste(c("diagnostics/engine_status.json", "diagnostics/session_info.txt", paste0("diagnostics/", head(error_files, 3))), collapse = "; ")
      } else {
        ""
      }
      effects_note <- if (has_effects) "Estimates present; confirm response, term, contrast and link before comparing directions." else "No finite fitted effect estimates."
      prediction_note <- if (has_predictions) "Predictions present; a prediction table does not establish independent validation." else "No finite fitted predictions."
      association_numeric_note <- paste(
        "Not directly comparable numerically across engines:",
        "Check the response scale, random level and statistical quantity before any comparison. Do not compare latent-factor loadings to species correlations."
      )
      association_pattern_note <- if (has_associations) {
        paste(association_note, "Inspect within fit; cross-model pattern agreement has not been checked.")
      } else {
        "No fitted standard/associations_long.csv available; compare diagnostics only."
      }
      collapse_nonempty <- function(x) paste(x[nzchar(x)], collapse = "; ")
      comparable_outputs <- collapse_nonempty(c(
        "status and output-contract completeness",
        if (safe_file(config_file)) "used_config.yml settings" else "",
        if (safe_file(standard_manifest_file)) "output manifest and file counts" else ""
      ))
      non_comparable_outputs <- collapse_nonempty(c(
        if (!fitted_status) "performance/effects are not comparable because status is not fitted" else "",
        "raw coefficients without checking scaling/link/family",
        "raw residual association parameters across engines",
        "WAIC/DIC/AUC/RMSE-like values across mismatched response scales or validation units",
        if (identical(engine, "spOccupancy")) "detection effects versus occurrence/environment effects" else "",
        if (identical(engine, "GJAM")) "GJAM observation-scale typeNames outputs versus link-scale JSDM coefficients" else ""
      ))
      primary_diagnostics <- if (identical(status_text, "not_run")) {
        "No run folder selected"
      } else if (!exists) {
        "Folder does not exist"
      } else if (status_text %in% c("engine_mismatch", "status_conflict", "missing_status")) {
        "Check the selected engine and folder metadata; this folder was not modified."
      } else if (identical(status_text, "check_failed")) {
        "diagnostics/data_check_messages.csv; used_config.yml; diagnostics/session_info.txt"
      } else if (identical(status_text, "fit_failed")) {
        paste(c("diagnostics/engine_status.json", "diagnostics/session_info.txt", paste0("diagnostics/", head(error_files, 3))), collapse = "; ")
      } else {
        "diagnostics/engine_status.json; standard/run_summary.csv; standard/output_manifest.csv; diagnostics/session_info.txt"
      }
      data.frame(
        Engine = engine,
        Folder = folder,
        Exists = exists,
        Engine_status = status_text,
        Engine_status_raw = status_text_raw,
        Status_class = status_class,
        Fitted_for_metric_comparison = FALSE,
        Metric_comparison_allowed = FALSE,
        Has_finite_fit_metrics = has_metrics,
        Has_finite_effects = has_effects,
        Has_finite_predictions = has_predictions,
        Has_finite_associations = has_associations,
        Metric_comparison_note = metric_note,
        Effects_direction_comparison_allowed = FALSE,
        Effects_direction_comparison_note = effects_note,
        Prediction_comparison_allowed = FALSE,
        Prediction_comparison_note = prediction_note,
        Association_numeric_comparison_allowed = FALSE,
        Association_pattern_comparison_allowed = FALSE,
        Association_comparison_note = association_note,
        Association_pattern_note = association_pattern_note,
        Association_numeric_note = association_numeric_note,
        Comparable_outputs = comparable_outputs,
        Non_comparable_outputs = non_comparable_outputs,
        Primary_diagnostics = primary_diagnostics,
        Failure_diagnostic_file = failure_diag,
        Standard_contract_complete = standard_complete,
        Status_file = safe_file(status_file),
        Config = safe_file(config_file),
        Data_check_messages = safe_file(data_check_file),
        Session_info = safe_file(session_file),
        Standard_run_summary = safe_file(standard_run_file),
        Standard_fit_metrics = safe_file(standard_fit_file),
        Standard_effects = safe_file(standard_effects_file),
        Standard_predictions = safe_file(standard_predictions_file),
        Standard_associations = safe_file(standard_associations_file),
        Standard_diagnostics = safe_file(standard_diagnostics_file),
        Standard_effects_species_environment = safe_file(standard_effects_species_file),
        Standard_predictions_site_species = safe_file(standard_predictions_site_file),
        Standard_associations_species_species = safe_file(standard_associations_species_file),
        Output_manifest = safe_file(standard_manifest_file),
        Workflow_map = safe_file(map_file),
        Workflow_step_index = safe_file(workflow_index_file),
        Report = safe_file(report_file),
        Zip_file = if (exists) file.exists(paste0(normalizePath(dir, winslash="/", mustWork=FALSE), ".zip")) else FALSE,
        N_sites = if (!is.null(standard_run) && "n_sites" %in% names(standard_run)) standard_run$n_sites[1] else NA,
        N_responses = if (!is.null(standard_run) && "n_responses" %in% names(standard_run)) standard_run$n_responses[1] else NA,
        N_predictors = if (!is.null(standard_run) && "n_predictors" %in% names(standard_run)) standard_run$n_predictors[1] else NA,
        Tables_count = safe_count("tables"),
        Results_count = safe_count("results"),
        Plots_count = safe_count("plots"),
        Predictions_count = safe_count("predictions"),
        Diagnostics_count = safe_count("diagnostics"),
        stringsAsFactors = FALSE
      )
    }

    rows <- list()
    if ("Hmsc" %in% selected) rows[[length(rows)+1]] <- add("Hmsc", hdir)
    if ("Hmsc-HPC" %in% selected) rows[[length(rows)+1]] <- add("Hmsc-HPC", hpdir)
    if ("jSDM" %in% selected) rows[[length(rows)+1]] <- add("jSDM", jdir)
    if ("GJAM" %in% selected) rows[[length(rows)+1]] <- add("GJAM", gdir)
    if ("spOccupancy" %in% selected) rows[[length(rows)+1]] <- add("spOccupancy", sdir)
    if ("sjSDM" %in% selected) rows[[length(rows)+1]] <- add("sjSDM", sjdir)
    if ("boral" %in% selected) rows[[length(rows)+1]] <- add("boral", bdir)
    rv$compare <- do.call(rbind, rows)
  })

  output$compare_table <- renderDT({
    runs <- rv$compare
    if (is.null(runs)) return(datatable(data.frame(Message = "No run folders inspected yet."),
                                        options = list(dom = "tip"), rownames = FALSE))
    present <- function(x) ifelse(x, "Present", "Unavailable")
    summary <- data.frame(
      Engine = runs$Engine, Status = runs$Engine_status,
      Effects = present(runs$Has_finite_effects),
      Predictions = present(runs$Has_finite_predictions),
      Associations = present(runs$Has_finite_associations),
      Metrics = present(runs$Has_finite_fit_metrics),
      Diagnostics = runs$Primary_diagnostics, check.names = FALSE
    )
    datatable(summary, options = list(pageLength = 7, dom = "tip", scrollX = TRUE), rownames = FALSE)
  })

  output$compare_outputs_matrix <- renderDT({
    dat <- do.call(rbind, lapply(names(studio_engine_catalog), function(engine) {
      entry <- studio_engine_catalog[[engine]]
      data.frame(Engine = engine, Scope = entry[[2]], Inputs = entry[[3]], Outputs = entry[[4]])
    }))
    datatable(dat, options = list(pageLength = 7, dom = "tip", scrollX = TRUE), rownames = FALSE)
  })

  output$compare_notes <- renderText({
    if (is.null(rv$compare)) return("")
    paste("Availability only: statistical comparability has not been verified.",
          "Use the full audit CSV for paths, status and file checks.",
          sep = "\n")
  })

  output$parameter_dictionary_table <- renderDT({
    dat <- parameter_dictionary_data()
    datatable(
      dat,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 15, scrollX = TRUE, autoWidth = TRUE)
    )
  })

  output$parameter_dictionary_download <- downloadHandler(
    filename = function() paste0("JSDMStudio_parameter_dictionary_", timestamp_id(), ".csv"),
    content = function(file) {
      write.csv(parameter_dictionary_data(), file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  output$compare_download <- downloadHandler(
    filename = function() paste0("JSDMStudio_model_comparison_", timestamp_id(), ".csv"),
    content = function(file) {
      if (is.null(rv$compare)) {
        write.csv(data.frame(Message="No comparison created yet."), file, row.names=FALSE)
      } else {
        write.csv(rv$compare, file, row.names=FALSE)
      }
    }
  )

}

shinyApp(ui, server)
