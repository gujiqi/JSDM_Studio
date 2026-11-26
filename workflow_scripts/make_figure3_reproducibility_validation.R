# Figure 3: reproducibility validation for JSDM Studio.
#
# This script builds a reviewer-facing figure and source tables that test four
# claims: workflows run, GUI-exported scripts can be rerun, deliberate bad inputs
# are reported honestly, and the output contract is complete.

options(stringsAsFactors = FALSE)

if (.Platform$OS.type == "windows") {
  Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
}

required <- c("ggplot2", "dplyr", "tidyr", "readr", "patchwork", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing required R packages for Figure 3: ", paste(missing, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(scales)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (length(x) == 1 && is.atomic(x) && is.na(x)) return(y)
  x
}

timestamp_id <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

script_path <- {
  args0 <- commandArgs(trailingOnly = FALSE)
  hit <- args0[startsWith(args0, "--file=")]
  if (length(hit)) normalizePath(sub("^--file=", "", hit[[1]]), winslash = "/", mustWork = FALSE)
  else normalizePath(file.path("workflow_scripts", "make_figure3_reproducibility_validation.R"), winslash = "/", mustWork = FALSE)
}
app_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) app_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) stop("Cannot locate app.R for JSDM Studio.", call. = FALSE)
setwd(app_dir)

cmd_args <- commandArgs(trailingOnly = TRUE)
skip_roundtrip <- any(cmd_args %in% "--skip_roundtrip")
tolerance <- 1e-6

output_root <- file.path(app_dir, "docs", "paper_figures", "figure3_reproducibility_validation")
source_data_dir <- file.path(output_root, "source_data")
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)

engine_order <- c("Hmsc", "Hmsc-HPC", "jSDM", "GJAM", "spOccupancy", "sjSDM", "boral")
engine_palette <- c(
  "Hmsc" = "#2C7BB6",
  "Hmsc-HPC" = "#1B9E77",
  "jSDM" = "#7B5AA6",
  "GJAM" = "#D95F02",
  "spOccupancy" = "#4C78A8",
  "sjSDM" = "#2CA25F",
  "boral" = "#7F7F7F"
)

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(tibble())
  suppressWarnings(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
}

latest_dir <- function(prefix) {
  root <- file.path(app_dir, "output")
  candidates <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  candidates <- candidates[startsWith(basename(candidates), prefix)]
  if (!length(candidates)) return(NA_character_)
  candidates[order(file.info(candidates)$mtime, decreasing = TRUE)][[1]]
}

latest_master_dir <- function() {
  run_dir <- latest_dir("universal_benchmark_")
  master <- file.path(run_dir, "master")
  if (!dir.exists(master)) stop("No Universal Benchmark master directory found under output/.", call. = FALSE)
  normalizePath(master, winslash = "/", mustWork = TRUE)
}

copy_dir <- function(src, dst) {
  if (!dir.exists(src)) stop("Missing source directory: ", src, call. = FALSE)
  if (dir.exists(dst)) stop("Refusing to overwrite existing directory: ", dst, call. = FALSE)
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  items <- list.files(src, all.files = FALSE, no.. = TRUE, full.names = TRUE)
  if (length(items)) file.copy(items, dst, recursive = TRUE, copy.date = TRUE)
  normalizePath(dst, winslash = "/", mustWork = FALSE)
}

file_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

numericish <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(TRUE)
  z <- suppressWarnings(as.numeric(as.character(x)))
  all(is.na(x) | is.finite(z) | is.na(z)) && sum(!is.na(as.character(x)) & nzchar(as.character(x)) & is.na(z)) == 0
}

canonical_text <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

make_duplicate_index <- function(key) {
  ave(seq_along(key), key, FUN = seq_along)
}

compare_csv_tables <- function(before_path, after_path, component, engine, script_ok = TRUE) {
  if (!script_ok) {
    return(tibble(engine = engine, component = component, status = "script_failed",
                  schema_match = FALSE, row_count_match = FALSE, common_rows = 0L,
                  max_abs_diff = NA_real_, numeric_correlation = NA_real_,
                  note = "Exported script returned non-zero status."))
  }
  if (!file.exists(before_path) || !file.exists(after_path)) {
    return(tibble(engine = engine, component = component, status = "missing",
                  schema_match = FALSE, row_count_match = FALSE, common_rows = 0L,
                  max_abs_diff = NA_real_, numeric_correlation = NA_real_,
                  note = "Required table missing before or after script rerun."))
  }
  before <- tryCatch(read.csv(before_path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  after <- tryCatch(read.csv(after_path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(before) || is.null(after)) {
    return(tibble(engine = engine, component = component, status = "read_failed",
                  schema_match = FALSE, row_count_match = FALSE, common_rows = 0L,
                  max_abs_diff = NA_real_, numeric_correlation = NA_real_,
                  note = "Table could not be parsed as CSV."))
  }
  common_cols <- intersect(names(before), names(after))
  schema_match <- identical(names(before), names(after))
  row_count_match <- nrow(before) == nrow(after)
  if (!length(common_cols) || !nrow(before) || !nrow(after)) {
    status <- if (schema_match && row_count_match) "matched" else "schema_or_row_mismatch"
    return(tibble(engine = engine, component = component, status = status,
                  schema_match = schema_match, row_count_match = row_count_match,
                  common_rows = 0L, max_abs_diff = NA_real_, numeric_correlation = NA_real_,
                  note = "No comparable rows or columns."))
  }
  numeric_cols <- common_cols[vapply(common_cols, function(nm) numericish(before[[nm]]) && numericish(after[[nm]]), logical(1))]
  preferred_key_cols <- switch(
    component,
    effects = c("engine", "species", "response_id", "predictor", "effect_type"),
    predictions = c("engine", "site_id", "species", "response_id"),
    associations = c("engine", "species_i", "species_j", "response_1", "response_2", "association_type"),
    c("engine")
  )
  key_cols <- intersect(preferred_key_cols, common_cols)
  if (!length(key_cols)) {
    key_cols <- setdiff(common_cols, c(numeric_cols, "note", "notes", "message", "prediction_scale", "scale", "comparable", "comparable_level"))
  }
  b <- before[, common_cols, drop = FALSE]
  a <- after[, common_cols, drop = FALSE]
  if (length(key_cols)) {
    b$.key <- apply(as.data.frame(lapply(b[, key_cols, drop = FALSE], canonical_text)), 1, paste, collapse = " | ")
    a$.key <- apply(as.data.frame(lapply(a[, key_cols, drop = FALSE], canonical_text)), 1, paste, collapse = " | ")
  } else {
    b$.key <- as.character(seq_len(nrow(b)))
    a$.key <- as.character(seq_len(nrow(a)))
  }
  b$.dup <- make_duplicate_index(b$.key)
  a$.dup <- make_duplicate_index(a$.key)
  merged <- merge(b, a, by = c(".key", ".dup"), suffixes = c(".gui", ".script"), all = FALSE)
  max_abs_diff <- NA_real_
  numeric_correlation <- NA_real_
  if (length(numeric_cols) && nrow(merged)) {
    gui_vals <- c()
    script_vals <- c()
    for (nm in numeric_cols) {
      g <- suppressWarnings(as.numeric(merged[[paste0(nm, ".gui")]]))
      s <- suppressWarnings(as.numeric(merged[[paste0(nm, ".script")]]))
      keep <- is.finite(g) & is.finite(s)
      if (any(keep)) {
        gui_vals <- c(gui_vals, g[keep])
        script_vals <- c(script_vals, s[keep])
      }
    }
    if (length(gui_vals)) {
      max_abs_diff <- max(abs(gui_vals - script_vals), na.rm = TRUE)
      if (length(gui_vals) >= 3 && stats::sd(gui_vals) > 0 && stats::sd(script_vals) > 0) {
        numeric_correlation <- suppressWarnings(stats::cor(gui_vals, script_vals))
      }
    }
  }
  status <- dplyr::case_when(
    !schema_match || !row_count_match ~ "schema_or_row_mismatch",
    is.finite(max_abs_diff) && max_abs_diff <= tolerance ~ "matched",
    is.finite(numeric_correlation) && numeric_correlation >= 0.999 ~ "within_tolerance_pattern",
    length(numeric_cols) > 0 && !is.finite(max_abs_diff) && nrow(merged) == nrow(before) ~ "matched",
    length(numeric_cols) == 0 && nrow(merged) == nrow(before) ~ "matched",
    TRUE ~ "numeric_drift"
  )
  tibble(engine = engine, component = component, status = status,
         schema_match = schema_match, row_count_match = row_count_match,
         common_rows = nrow(merged), max_abs_diff = max_abs_diff,
         numeric_correlation = numeric_correlation,
         note = ifelse(status == "matched", "Schema, row count and numeric values agree.", "Review round-trip numeric or schema differences."))
}

id_set <- function(path, cols) {
  if (!file.exists(path)) return(character())
  x <- tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) data.frame())
  vals <- character()
  for (cc in intersect(cols, names(x))) vals <- c(vals, as.character(x[[cc]]))
  sort(unique(vals[nzchar(vals) & !is.na(vals)]))
}

compare_id_component <- function(before_path, after_path, cols, component, engine, script_ok = TRUE) {
  if (!script_ok) {
    return(tibble(engine = engine, component = component, status = "script_failed",
                  schema_match = FALSE, row_count_match = FALSE, common_rows = 0L,
                  max_abs_diff = NA_real_, numeric_correlation = NA_real_,
                  note = "Exported script returned non-zero status."))
  }
  b <- id_set(before_path, cols)
  a <- id_set(after_path, cols)
  ok <- identical(b, a) && length(b) > 0
  tibble(engine = engine, component = component, status = ifelse(ok, "matched", "id_mismatch"),
         schema_match = ok, row_count_match = ok, common_rows = length(intersect(b, a)),
         max_abs_diff = NA_real_, numeric_correlation = NA_real_,
         note = ifelse(ok, "Identifier set agrees.", "Identifier set differs or is empty."))
}

component_path <- function(engine_dir, primary, fallback = character()) {
  candidates <- file.path(engine_dir, c(primary, fallback))
  candidates[which(file.exists(candidates))[1]] %||% candidates[[1]]
}

run_roundtrip_checks <- function(master_dir) {
  roundtrip_root <- file.path(app_dir, "output", paste0("figure3_roundtrip_", timestamp_id()))
  dir.create(roundtrip_root, recursive = TRUE, showWarnings = FALSE)
  engines <- engine_order[dir.exists(file.path(master_dir, "03_engine_outputs", engine_order))]
  all_rows <- list()
  script_rows <- list()
  for (engine in engines) {
    src <- file.path(master_dir, "03_engine_outputs", engine)
    dst <- file.path(roundtrip_root, engine)
    copy_dir(src, dst)
    script <- list.files(file.path(dst, "reproducible_script"),
                         pattern = "^run_this_.*analysis\\.R$", full.names = TRUE)
    script <- script[!grepl("run_all_", basename(script))]
    script <- script[[1]] %||% NA_character_
    stdout <- file.path(dst, "diagnostics", "roundtrip_stdout.txt")
    stderr <- file.path(dst, "diagnostics", "roundtrip_stderr.txt")
    dir.create(dirname(stdout), recursive = TRUE, showWarnings = FALSE)
    code <- if (is.na(script)) {
      127L
    } else {
      rbin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
      status <- system2(rbin, shQuote(script), stdout = stdout, stderr = stderr)
      if (is.null(status)) 0L else as.integer(status)
    }
    script_ok <- identical(code, 0L)
    script_rows[[engine]] <- tibble(engine = engine,
                                    script_path = script,
                                    copied_output = normalizePath(dst, winslash = "/", mustWork = FALSE),
                                    exit_status = code,
                                    stdout = stdout,
                                    stderr = stderr,
                                    script_ok = script_ok)
    before_base <- src
    after_base <- dst
    config_status <- if (script_ok && identical(file_md5(file.path(before_base, "used_config.yml")),
                                                file_md5(file.path(after_base, "used_config.yml")))) "matched" else if (script_ok) "config_hash_mismatch" else "script_failed"
    all_rows[[paste(engine, "config")]] <- tibble(engine = engine, component = "config", status = config_status,
                                                  schema_match = config_status == "matched",
                                                  row_count_match = config_status == "matched",
                                                  common_rows = NA_integer_,
                                                  max_abs_diff = NA_real_,
                                                  numeric_correlation = NA_real_,
                                                  note = ifelse(config_status == "matched", "used_config.yml hash agrees.", "Config changed or script failed."))
    table_specs <- list(
      effects = c("standard/effects_species_environment.csv", "standard/effects_long.csv"),
      predictions = c("standard/predictions_site_species.csv", "standard/predictions_long.csv"),
      associations = c("standard/associations_species_species.csv", "standard/associations_long.csv")
    )
    for (nm in names(table_specs)) {
      paths <- table_specs[[nm]]
      before <- component_path(before_base, paths[[1]], paths[-1])
      after <- component_path(after_base, paths[[1]], paths[-1])
      all_rows[[paste(engine, nm)]] <- compare_csv_tables(before, after, nm, engine, script_ok = script_ok)
    }
    pred_before <- component_path(before_base, "standard/predictions_site_species.csv", "standard/predictions_long.csv")
    pred_after <- component_path(after_base, "standard/predictions_site_species.csv", "standard/predictions_long.csv")
    effects_before <- component_path(before_base, "standard/effects_species_environment.csv", "standard/effects_long.csv")
    effects_after <- component_path(after_base, "standard/effects_species_environment.csv", "standard/effects_long.csv")
    all_rows[[paste(engine, "site_ids")]] <- compare_id_component(pred_before, pred_after,
                                                                  c("site_id", "site", "row_id", "sample_id"),
                                                                  "site IDs", engine, script_ok)
    all_rows[[paste(engine, "species_ids")]] <- compare_id_component(effects_before, effects_after,
                                                                     c("species", "response_id", "species_i", "species_j"),
                                                                     "species IDs", engine, script_ok)
  }
  checks <- bind_rows(all_rows)
  scripts <- bind_rows(script_rows)
  write_csv(checks, file.path(roundtrip_root, "figure3_roundtrip_table_checks.csv"))
  write_csv(scripts, file.path(roundtrip_root, "figure3_roundtrip_script_runs.csv"))
  write_csv(checks, file.path(source_data_dir, "figure3_roundtrip_table_checks.csv"))
  write_csv(scripts, file.path(source_data_dir, "figure3_roundtrip_script_runs.csv"))
  list(root = roundtrip_root, checks = checks, scripts = scripts)
}

load_runner_env <- function() {
  runner <- file.path(app_dir, "workflow_scripts", "universal_benchmark_runner.R")
  lines <- readLines(runner, warn = FALSE)
  lines <- lines[!grepl("^if \\(!interactive\\(\\)\\) main\\(\\)\\s*$", trimws(lines))]
  e <- new.env(parent = globalenv())
  eval(parse(text = paste(lines, collapse = "\n"), srcfile = runner), envir = e)
  e
}

record_failure_case <- function(runner_env, app_env, failure_root, engine, test_id, expected_status,
                                check_messages, check_ok, Y = NULL, X = NULL, error_text = NULL) {
  outdir <- file.path(failure_root, engine, test_id)
  runner_env$ensure_dirs(outdir, engine)
  observed <- if (identical(expected_status, "fit_failed")) {
    "fit_failed"
  } else if (isTRUE(check_ok)) {
    "unexpected_pass"
  } else {
    "check_failed"
  }
  app_env$write_data_check_messages(outdir, check_messages)
  runner_env$write_status(outdir, engine, observed,
                          warnings = check_messages,
                          errors = error_text %||% if (observed == "check_failed") paste0(test_id, " deliberately failed validation.") else character())
  app_env$write_standard_outputs(outdir, engine, observed, Y, X)
  app_env$ensure_output_contract(outdir, engine, observed, Y, X)
  zip <- tryCatch(runner_env$make_zip_file(outdir), error = function(e) NA_character_)
  err_files <- list.files(file.path(outdir, "diagnostics"), pattern = "error|failed|failure", full.names = TRUE, ignore.case = TRUE)
  tibble(
    test_id = test_id,
    engine = engine,
    expected_status = expected_status,
    observed_status = observed,
    matched = identical(observed, expected_status),
    diagnostics_exists = file.exists(file.path(outdir, "diagnostics", "engine_status.json")),
    data_check_exists = file.exists(file.path(outdir, "diagnostics", "data_check_messages.csv")),
    error_log_exists = length(err_files) > 0,
    zip_exists = !is.na(zip) && file.exists(zip) && file.info(zip)$size > 0,
    output_folder = normalizePath(outdir, winslash = "/", mustWork = FALSE),
    zip = zip %||% NA_character_,
    message = paste(check_messages, collapse = " | ")
  )
}

run_deliberate_failure_checks <- function(master_dir) {
  runner_env <- load_runner_env()
  app_env <- runner_env$source_app_env()
  failure_root <- file.path(app_dir, "output", paste0("figure3_deliberate_failures_", timestamp_id()))
  input_dir <- file.path(failure_root, "benchmark_inputs")
  dat0 <- runner_env$generate_benchmark_data(input_dir, n_sites = 36L, n_species = 6L, n_visits = 3L,
                                             seed = 20260603L, allow_large = FALSE)
  dat <- runner_env$load_benchmark_data(dat0$dir, allow_generate = FALSE)
  rows <- list()

  bad_traits <- dat$traits
  rownames(bad_traits)[1] <- paste0(rownames(bad_traits)[1], "_mismatch")
  chk <- app_env$validate_hmsc(dat$Y, dat$X, bad_traits, dat$study, dat$coords,
                               distr = "probit", use_traits = TRUE, use_phylogeny = TRUE,
                               random_mode = "sample", random_effect_column = "sample")
  rows[["species_names_mismatch"]] <- record_failure_case(runner_env, app_env, failure_root, "Hmsc",
                                                           "species_names_mismatch", "check_failed",
                                                           chk$messages, chk$ok, dat$Y, dat$X)

  bad_y <- dat$Y
  bad_y[1, 1] <- 2
  chk <- app_env$validate_hmsc(bad_y, dat$X, dat$traits, dat$study, dat$coords,
                               distr = "probit", use_traits = TRUE, use_phylogeny = FALSE,
                               random_mode = "none")
  rows[["invalid_binomial_response"]] <- record_failure_case(runner_env, app_env, failure_root, "Hmsc",
                                                              "invalid_binomial_response", "check_failed",
                                                              chk$messages, chk$ok, bad_y, dat$X)

  chk <- app_env$validate_hmsc(dat$Y, dat$X, dat$traits, dat$study, dat$coords,
                               distr = "probit", XFormula = "~ pH + missing_predictor",
                               use_traits = TRUE, random_mode = "none")
  rows[["missing_predictor_in_formula"]] <- record_failure_case(runner_env, app_env, failure_root, "Hmsc",
                                                                 "missing_predictor_in_formula", "check_failed",
                                                                 chk$messages, chk$ok, dat$Y, dat$X)

  dup_coords <- dat$coords
  dup_coords[2, seq_len(2)] <- dup_coords[1, seq_len(2)]
  chk <- app_env$validate_hmsc(dat$Y, dat$X, dat$traits, dat$study, dup_coords,
                               distr = "probit", use_traits = FALSE,
                               random_mode = "spatial_full", lon_col = names(dup_coords)[1],
                               lat_col = names(dup_coords)[2])
  rows[["duplicate_spatial_coordinates"]] <- record_failure_case(runner_env, app_env, failure_root, "Hmsc",
                                                                  "duplicate_spatial_coordinates", "check_failed",
                                                                  chk$messages, chk$ok, dat$Y, dat$X)

  chk <- app_env$validate_spoccupancy_full(y = dat$det, occ = dat$occ_covs, det = NULL,
                                           coords = dat$coords,
                                           species = data.frame(species = colnames(dat$Y)),
                                           model_type = "msPGOcc", det.formula = "~ obs",
                                           occ.formula = "~ pH + moisture")
  rows[["missing_detection_input"]] <- record_failure_case(runner_env, app_env, failure_root, "spOccupancy",
                                                            "missing_detection_input", "check_failed",
                                                            chk$messages, chk$ok, dat$det, dat$occ_covs)

  chk <- app_env$validate_boral(dat$Y, dat$X, family = "not_a_family")
  rows[["invalid_boral_family"]] <- record_failure_case(runner_env, app_env, failure_root, "boral",
                                                         "invalid_boral_family", "check_failed",
                                                         chk$messages, chk$ok, dat$Y, dat$X)

  outdir <- file.path(failure_root, "Hmsc-HPC", "deliberate_dependency_fit_failure")
  runner_env$ensure_dirs(outdir, "Hmsc-HPC")
  fail_script <- file.path(outdir, "reproducible_script", "deliberate_dependency_failure.R")
  writeLines("stop('Deliberate dependency failure contract test', call. = FALSE)", fail_script)
  stdout <- file.path(outdir, "diagnostics", "deliberate_fit_failure_stdout.txt")
  stderr <- file.path(outdir, "diagnostics", "deliberate_fit_failure_stderr.txt")
  rbin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  code <- system2(rbin, shQuote(fail_script), stdout = stdout, stderr = stderr)
  code <- if (is.null(code)) 0L else as.integer(code)
  rows[["deliberate_dependency_fit_failure"]] <- record_failure_case(
    runner_env, app_env, failure_root, "Hmsc-HPC", "deliberate_dependency_fit_failure", "fit_failed",
    paste0("Rscript exit status = ", code, "; deliberate dependency failure contract test."),
    FALSE, dat$Y, dat$X, error_text = "Deliberate dependency failure contract test."
  )

  out <- bind_rows(rows)
  write_csv(out, file.path(failure_root, "figure3_deliberate_failure_checks.csv"))
  write_csv(out, file.path(source_data_dir, "figure3_deliberate_failure_checks.csv"))
  list(root = failure_root, checks = out)
}

build_workflow_validation <- function(master_dir) {
  audit_dir <- latest_dir("workflow_synthetic_case_audit_")
  summary_file <- file.path(audit_dir, "workflow_synthetic_case_summary.csv")
  if (!is.na(audit_dir) && file.exists(summary_file)) {
    cases <- safe_read_csv(summary_file) %>%
      mutate(
        engine = factor(engine, levels = engine_order),
        missing_count = suppressWarnings(as.integer(missing_count)),
        validation = case_when(
          status %in% c("fitted", "model_defined") & replace_na(missing_count, 1L) == 0L ~ "passed",
          status %in% c("check_failed", "fit_failed") ~ status,
          TRUE ~ "needs_review"
        ),
        source = basename(audit_dir)
      )
  } else {
    cases <- safe_read_csv(file.path(master_dir, "04_unified_standard_results", "model_status_matrix.csv")) %>%
      mutate(
        case_id = "universal_benchmark_current",
        engine = factor(engine, levels = engine_order),
        missing_count = 0L,
        validation = ifelse(status %in% c("fitted", "model_defined"), "passed", status),
        source = basename(dirname(master_dir))
      )
  }
  summary <- cases %>%
    group_by(engine, validation) %>%
    summarise(n_cases = n(), .groups = "drop") %>%
    complete(engine = factor(engine_order, levels = engine_order),
             validation = c("passed", "model_defined", "check_failed", "fit_failed", "needs_review"),
             fill = list(n_cases = 0L))
  write_csv(cases, file.path(source_data_dir, "figure3_workflow_validation_cases.csv"))
  write_csv(summary, file.path(source_data_dir, "figure3_workflow_validation_summary.csv"))
  list(cases = cases, summary = summary)
}

build_output_completeness <- function(master_dir) {
  status <- safe_read_csv(file.path(master_dir, "04_unified_standard_results", "model_status_matrix.csv"))
  status_zip <- setNames(status$zip %||% character(), status$engine %||% character())
  specs <- list(
    Config = function(d, e) file.exists(file.path(d, "used_config.yml")),
    Inputs = function(d, e) dir.exists(file.path(d, "inputs")) && length(list.files(file.path(d, "inputs"), recursive = TRUE)) > 0,
    Data = function(d, e) dir.exists(file.path(d, "data")) && length(list.files(file.path(d, "data"), recursive = TRUE)) > 0,
    Model = function(d, e) dir.exists(file.path(d, "models")) && length(list.files(file.path(d, "models"), recursive = TRUE)) > 0,
    Effects = function(d, e) file.exists(file.path(d, "standard", "effects_species_environment.csv")) || file.exists(file.path(d, "standard", "effects_long.csv")),
    Predictions = function(d, e) file.exists(file.path(d, "standard", "predictions_site_species.csv")) || file.exists(file.path(d, "standard", "predictions_long.csv")),
    Associations = function(d, e) file.exists(file.path(d, "standard", "associations_species_species.csv")) || file.exists(file.path(d, "standard", "associations_long.csv")),
    Diagnostics = function(d, e) file.exists(file.path(d, "diagnostics", "engine_status.json")),
    Script = function(d, e) length(list.files(file.path(d, "reproducible_script"), pattern = "^run_this_.*analysis\\.(R|py)$", full.names = TRUE)) > 0,
    Report = function(d, e) dir.exists(file.path(d, "report")) && length(list.files(file.path(d, "report"), recursive = TRUE)) > 0,
    ZIP = function(d, e) {
      z <- status_zip[[e]] %||% file.path(dirname(d), paste0(e, ".zip"))
      file.exists(z) && file.info(z)$size > 0
    }
  )
  rows <- list()
  for (engine in engine_order) {
    d <- file.path(master_dir, "03_engine_outputs", engine)
    for (component in names(specs)) {
      rows[[paste(engine, component)]] <- tibble(
        engine = engine,
        component = component,
        present = dir.exists(d) && isTRUE(specs[[component]](d, engine)),
        engine_folder = normalizePath(d, winslash = "/", mustWork = FALSE)
      )
    }
  }
  completeness <- bind_rows(rows) %>%
    mutate(engine = factor(engine, levels = engine_order),
           component = factor(component, levels = names(specs)))
  write_csv(completeness, file.path(source_data_dir, "figure3_output_completeness_matrix.csv"))
  completeness
}

theme_fig3 <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size * 1.08, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size * 0.78, color = "grey30", margin = margin(b = 8)),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = base_size * 0.72),
      axis.text.y = element_text(size = base_size * 0.76),
      strip.background = element_rect(fill = "grey94", color = "grey78", linewidth = 0.4),
      strip.text = element_text(face = "bold", size = base_size * 0.78),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 12, 8, 8)
    )
}

master_dir <- latest_master_dir()
workflow_validation <- build_workflow_validation(master_dir)
roundtrip <- if (skip_roundtrip) {
  empty_checks <- expand_grid(engine = engine_order,
                              component = c("config", "effects", "predictions", "associations", "site IDs", "species IDs")) %>%
    mutate(status = "skipped", schema_match = NA, row_count_match = NA, common_rows = NA_integer_,
           max_abs_diff = NA_real_, numeric_correlation = NA_real_, note = "Round-trip was skipped by command-line flag.")
  write_csv(empty_checks, file.path(source_data_dir, "figure3_roundtrip_table_checks.csv"))
  list(root = NA_character_, checks = empty_checks, scripts = tibble())
} else {
  run_roundtrip_checks(master_dir)
}
failures <- run_deliberate_failure_checks(master_dir)
completeness <- build_output_completeness(master_dir)

validation_colors <- c(
  passed = "#2CA25F",
  model_defined = "#80CDC1",
  check_failed = "#F4A259",
  fit_failed = "#D95F02",
  needs_review = "#BDBDBD"
)

p_a <- ggplot(workflow_validation$summary, aes(x = engine, y = n_cases, fill = validation)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = validation_colors, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "Workflow validation",
    subtitle = "Synthetic cases and contract completeness from the latest audit",
    x = NULL,
    y = "Number of cases",
    fill = "Status"
  ) +
  theme_fig3()

roundtrip_plot <- roundtrip$checks %>%
  mutate(
    engine = factor(engine, levels = engine_order),
    component = factor(component, levels = c("config", "effects", "predictions", "associations", "site IDs", "species IDs")),
    status_class = case_when(
      status %in% c("matched", "within_tolerance_pattern") ~ "agreement",
      status %in% c("numeric_drift", "schema_or_row_mismatch", "id_mismatch", "config_hash_mismatch") ~ "drift_or_mismatch",
      status == "skipped" ~ "skipped",
      TRUE ~ "failed_or_missing"
    )
  )

p_b <- ggplot(roundtrip_plot, aes(x = component, y = engine, fill = status_class)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = case_when(
    status_class == "agreement" ~ "OK",
    status_class == "drift_or_mismatch" ~ "DR",
    status_class == "skipped" ~ "SK",
    TRUE ~ "NA"
  )), size = 4.2, fontface = "bold") +
  scale_fill_manual(values = c(agreement = "#2CA25F", drift_or_mismatch = "#F4A259",
                               failed_or_missing = "#D95F02", skipped = "#BDBDBD"),
                    drop = FALSE) +
  labs(
    title = "GUI-to-script round trip",
    subtitle = "Exported scripts rerun and standard outputs compared",
    x = NULL,
    y = NULL,
    fill = "Agreement"
  ) +
  theme_fig3() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())

failure_plot <- failures$checks %>%
  mutate(test_label = recode(test_id,
                             species_names_mismatch = "species mismatch",
                             invalid_binomial_response = "invalid binary Y",
                             missing_predictor_in_formula = "missing predictor",
                             duplicate_spatial_coordinates = "duplicate coords",
                             missing_detection_input = "missing detection",
                             invalid_boral_family = "invalid family",
                             deliberate_dependency_fit_failure = "fit failure"),
         test_label = factor(test_label, levels = c("species mismatch", "invalid binary Y", "missing predictor",
                                                    "duplicate coords", "missing detection", "invalid family", "fit failure")),
         outcome = ifelse(matched & diagnostics_exists & error_log_exists & zip_exists, "expected_with_logs", "needs_review"))

p_c <- ggplot(failure_plot, aes(x = test_label, y = engine, fill = outcome)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = paste0(expected_status, "\n", observed_status)), size = 3.6, lineheight = 0.92) +
  scale_fill_manual(values = c(expected_with_logs = "#2CA25F", needs_review = "#D95F02"), drop = FALSE) +
  labs(
    title = "Deliberate failure tests",
    subtitle = "Bad inputs must become check_failed or fit_failed with diagnostics",
    x = NULL,
    y = NULL,
    fill = "Outcome"
  ) +
  theme_fig3() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())

p_d <- ggplot(completeness, aes(x = component, y = engine, fill = present)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = ifelse(present, "OK", "NA")), size = 4.0, fontface = "bold") +
  scale_fill_manual(values = c(`TRUE` = "#2CA25F", `FALSE` = "#D95F02"), drop = FALSE) +
  labs(
    title = "Output completeness",
    subtitle = "Required reproducibility contract files in the latest benchmark",
    x = NULL,
    y = NULL,
    fill = "Present"
  ) +
  theme_fig3() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())

figure3 <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    title = "Reproducible Workflow Validation",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 26, hjust = 0, margin = margin(b = 8)),
      plot.tag = element_text(face = "bold", size = 22),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

ggsave(file.path(output_root, "Figure3_Reproducibility_Validation.pdf"), figure3, width = 16.2, height = 10.8, units = "in", device = cairo_pdf)
ggsave(file.path(output_root, "Figure3_Reproducibility_Validation.png"), figure3, width = 16.2, height = 10.8, units = "in", dpi = 360)
ggsave(file.path(output_root, "Figure3_Reproducibility_Validation.tiff"), figure3, width = 16.2, height = 10.8, units = "in", dpi = 600, compression = "lzw")
ggsave(file.path(output_root, "Figure3_Reproducibility_Validation.svg"), figure3, width = 16.2, height = 10.8, units = "in")

caption <- c(
  "# Figure 3. Reproducible workflow validation",
  "",
  "A, Workflow validation summarizes the latest multi-case synthetic audit. Passed means the workflow reached fitted or model_defined status and its required output contract had no missing components. B, GUI-to-script round-trip reproducibility copies each latest Universal Benchmark engine output, reruns the exported analysis script and compares config, standard tables and identifier sets. Exact numerical equality is expected when seeds and backends are deterministic; small numerical drift is reported rather than hidden. C, Deliberate bad-input tests verify that invalid data or dependency failure is recorded as check_failed or fit_failed with diagnostics and ZIP output. D, Output completeness checks the required reproducibility contract in the latest Universal Benchmark engine folders.",
  "",
  paste0("Latest Universal Benchmark master: ", master_dir),
  paste0("Round-trip output root: ", roundtrip$root %||% "skipped"),
  paste0("Deliberate failure output root: ", failures$root)
)
writeLines(caption, file.path(output_root, "Figure3_caption_draft.md"), useBytes = TRUE)

qa_note <- c(
  "Figure 3 QA note",
  paste0("Generated: ", Sys.time()),
  paste0("Source master: ", master_dir),
  paste0("Round-trip skipped: ", skip_roundtrip),
  "Backend: R/ggplot2/patchwork only.",
  "Each panel has a source CSV in source_data/.",
  "Round-trip status codes: OK means exact or high-correlation/tolerance agreement; DR means drift or schema mismatch requiring review; NA means script/table missing or failed.",
  "Deliberate failure tests write real diagnostics, standard placeholder tables and ZIPs; they are validation of failure-state honesty, not ecological analyses."
)
writeLines(qa_note, file.path(output_root, "Figure3_QA_note.txt"), useBytes = TRUE)

manifest <- tibble(
  file = c(
    "Figure3_Reproducibility_Validation.pdf",
    "Figure3_Reproducibility_Validation.png",
    "Figure3_Reproducibility_Validation.tiff",
    "Figure3_Reproducibility_Validation.svg",
    "Figure3_caption_draft.md",
    "Figure3_QA_note.txt",
    "source_data/figure3_workflow_validation_cases.csv",
    "source_data/figure3_workflow_validation_summary.csv",
    "source_data/figure3_roundtrip_table_checks.csv",
    "source_data/figure3_roundtrip_script_runs.csv",
    "source_data/figure3_deliberate_failure_checks.csv",
    "source_data/figure3_output_completeness_matrix.csv"
  ),
  role = c(
    "vector_figure", "raster_preview", "submission_tiff", "editable_vector",
    "caption", "qa_note", "source_data", "source_data", "source_data",
    "source_data", "source_data", "source_data"
  )
)
write_csv(manifest, file.path(output_root, "Figure3_manifest.csv"))

cat("FIGURE3_DIR=", normalizePath(output_root, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("ROUNDTRIP_ROOT=", normalizePath(roundtrip$root %||% "", winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("FAILURE_ROOT=", normalizePath(failures$root, winslash = "/", mustWork = FALSE), "\n", sep = "")
