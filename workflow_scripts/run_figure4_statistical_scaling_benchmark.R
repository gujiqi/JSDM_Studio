# Runs the Universal Benchmark at multiple data sizes, with longer sampler
# settings, and rebuilds Figure 4 from statistical recovery trends.

options(stringsAsFactors = FALSE)

if (.Platform$OS.type == "windows") {
  Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (length(x) == 1 && is.atomic(x) && is.na(x)) return(y)
  x
}

timestamp_id <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(
    project_root = "",
    collect_log_root = "",
    sizes = "Small:30:10,Medium:100:20,Large:200:50,Ultra-large:300:100",
    reps = "1",
    n_visits = "3",
    seed = "20260830",
    engines = "Hmsc,Hmsc-HPC,jSDM,GJAM,spOccupancy,sjSDM,boral",
    figure_master = "",
    truth_profile = "recovery",
    effect_scale = "1.60",
    association_strength = "1.65",
    residual_sd = "1.05",
    observation_noise_sd = "0.15",
    hmsc_samples = "500",
    hmsc_transient = "500",
    hmsc_thin = "1",
    hmsc_nChains = "2",
    hmschpc_samples = "500",
    hmschpc_transient = "500",
    hmschpc_thin = "1",
    hmschpc_chains = "1",
    jsdm_burnin = "200",
    jsdm_mcmc = "300",
    jsdm_thin = "5",
    gjam_ng = "700",
    gjam_burnin = "150",
    spocc_n_batch = "30",
    spocc_batch_length = "10",
    spocc_n_burn = "100",
    spocc_n_thin = "5",
    spocc_n_chains = "1",
    sjsdm_iter = "100",
    sjsdm_sampling = "1000",
    sjsdm_step_size = "32",
    sjsdm_learning_rate = "0.010",
    sjsdm_weight_decay = "0",
    sjsdm_biotic_lambda = "0.001",
    boral_n_burnin = "250",
    boral_n_iteration = "750",
    boral_n_thin = "10"
  )
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "true"
    out[[key]] <- val
  }
  out$reps <- max(1L, suppressWarnings(as.integer(out$reps)))
  out$n_visits <- max(1L, suppressWarnings(as.integer(out$n_visits)))
  out$seed <- suppressWarnings(as.integer(out$seed))
  if (!is.finite(out$seed)) out$seed <- 20260830L
  out
}

cmd <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", cmd, value = TRUE))
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg[1], winslash = "/", mustWork = TRUE))
} else {
  getwd()
}

args <- parse_args()
project_root <- if (nzchar(args$project_root)) {
  normalizePath(args$project_root, winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
}

runner_script <- file.path(project_root, "workflow_scripts", "universal_benchmark_runner.R")
figure_script <- file.path(project_root, "workflow_scripts", "make_figure4_statistical_universal_benchmark.R")
if (!file.exists(runner_script)) stop("Missing Universal Benchmark runner: ", runner_script)
if (!file.exists(figure_script)) stop("Missing Figure 4 script: ", figure_script)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
output_root <- file.path(project_root, "docs", "paper_figures", "figure4_statistical_universal_benchmark")
source_data_dir <- file.path(output_root, "source_data")
collect_only <- nzchar(trimws(args$collect_log_root %||% ""))
log_root <- if (collect_only) {
  normalizePath(args$collect_log_root, winslash = "/", mustWork = TRUE)
} else {
  file.path(project_root, "output", paste0("figure4_statistical_scaling_logs_", timestamp_id()))
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)

parse_sizes <- function(sizes) {
  parts <- trimws(unlist(strsplit(sizes, ",")))
  parts <- parts[nzchar(parts)]
  rows <- lapply(seq_along(parts), function(i) {
    bits <- strsplit(parts[[i]], ":", fixed = TRUE)[[1]]
    if (length(bits) != 3) stop("Bad --sizes entry: ", parts[[i]], ". Use Label:n_sites:n_species.")
    data.frame(
      size_order = i,
      size_label = bits[[1]],
      n_sites = as.integer(bits[[2]]),
      n_species = as.integer(bits[[3]]),
      stringsAsFactors = FALSE
    )
  })
  out <- bind_rows(rows)
  if (any(!is.finite(out$n_sites)) || any(!is.finite(out$n_species))) {
    stop("All --sizes n_sites and n_species values must be integers.")
  }
  out
}

sizes <- parse_sizes(args$sizes)
sampler_keys <- c(
  "hmsc_samples", "hmsc_transient", "hmsc_thin", "hmsc_nChains",
  "hmschpc_samples", "hmschpc_transient", "hmschpc_thin", "hmschpc_chains",
  "jsdm_burnin", "jsdm_mcmc", "jsdm_thin",
  "gjam_ng", "gjam_burnin",
  "spocc_n_batch", "spocc_batch_length", "spocc_n_burn", "spocc_n_thin", "spocc_n_chains",
  "sjsdm_iter", "sjsdm_sampling", "sjsdm_step_size",
  "sjsdm_learning_rate", "sjsdm_weight_decay", "sjsdm_biotic_lambda",
  "boral_n_burnin", "boral_n_iteration", "boral_n_thin"
)
truth_keys <- c("truth_profile", "effect_scale", "association_strength", "residual_sd", "observation_noise_sd")

stale_files <- c(
  "figure4_statistical_scaling_effect_metrics.csv",
  "figure4_statistical_scaling_prediction_metrics.csv",
  "figure4_statistical_scaling_association_metrics.csv",
  "figure4_statistical_scaling_status.csv",
  "figure4_statistical_scaling_run_log.csv",
  "figure4_runtime_scaling_metrics.csv",
  "figure4_runtime_scaling_run_log.csv",
  "figure4_runtime_scaling_summary.csv"
)
unlink(file.path(source_data_dir, stale_files), force = TRUE)
unlink(file.path(output_root, stale_files), force = TRUE)

write_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(enc2utf8(as.character(x)), path, useBytes = TRUE)
}

extract_master_dir <- function(lines) {
  hit <- grep("^MASTER_DIR=", lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub("^MASTER_DIR=", "", hit[[length(hit)]]))
}

run_system <- function(command, command_args, log_path) {
  lines <- tryCatch(
    system2(command, command_args, stdout = TRUE, stderr = TRUE),
    error = function(e) structure(conditionMessage(e), status = 1L)
  )
  status <- attr(lines, "status")
  if (is.null(status)) status <- 0L
  write_lines(lines, log_path)
  list(status = as.integer(status), lines = as.character(lines))
}

read_optional_csv <- function(path) {
  if (!file.exists(path)) return(tibble())
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

add_context <- function(tbl, scenario, rep_id, master_dir) {
  if (!nrow(tbl)) return(tbl)
  tbl %>%
    mutate(
      size_order = scenario$size_order,
      size_label = scenario$size_label,
      n_sites = scenario$n_sites,
      n_species = scenario$n_species,
      replicate = rep_id,
      benchmark_run = basename(dirname(master_dir)),
      master_dir = normalizePath(master_dir, winslash = "/", mustWork = FALSE),
      .before = 1
    )
}

make_figure_for_master <- function(master_dir, scenario, rep_id) {
  log_path <- file.path(log_root, sprintf("make_figure4_%s_rep%02d.log", scenario$size_label, rep_id))
  make_args <- c(
    normalizePath(figure_script, winslash = "/", mustWork = TRUE),
    paste0("--project_root=", project_root),
    paste0("--master_dir=", normalizePath(master_dir, winslash = "/", mustWork = FALSE))
  )
  res <- run_system(rscript, make_args, log_path)
  if (!identical(res$status, 0L)) {
    warning("Figure metric extraction failed for ", master_dir, "; see ", log_path)
  }
  list(
    effect = add_context(read_optional_csv(file.path(output_root, "figure4_environment_effect_recovery_metrics.csv")), scenario, rep_id, master_dir),
    prediction = add_context(read_optional_csv(file.path(output_root, "figure4_prediction_recovery_metrics.csv")), scenario, rep_id, master_dir),
    association = add_context(read_optional_csv(file.path(output_root, "figure4_association_recovery_metrics.csv")), scenario, rep_id, master_dir),
    status = add_context(read_optional_csv(file.path(master_dir, "04_unified_standard_results", "model_status_matrix.csv")), scenario, rep_id, master_dir),
    figure_status = res$status,
    figure_log = log_path
  )
}

run_one_benchmark <- function(scenario, rep_id) {
  run_seed <- args$seed + scenario$size_order * 1000L + rep_id
  target_dir <- file.path(log_root, "generated_inputs", sprintf("%02d_%s_rep%02d", scenario$size_order, scenario$size_label, rep_id))
  runner_args <- c(
    normalizePath(runner_script, winslash = "/", mustWork = TRUE),
    "--action=all",
    paste0("--n_sites=", scenario$n_sites),
    paste0("--n_species=", scenario$n_species),
    paste0("--n_visits=", args$n_visits),
    paste0("--seed=", run_seed),
    "--allow_large=true",
    paste0("--engines=", args$engines),
    paste0("--target_dir=", normalizePath(target_dir, winslash = "/", mustWork = FALSE))
  )
  runner_args <- c(
    runner_args,
    unlist(lapply(truth_keys, function(key) {
      val <- args[[key]] %||% ""
      if (!nzchar(as.character(val))) character() else paste0("--", key, "=", val)
    }), use.names = FALSE),
    unlist(lapply(sampler_keys, function(key) {
      val <- args[[key]] %||% ""
      if (!nzchar(as.character(val))) character() else paste0("--", key, "=", val)
    }), use.names = FALSE)
  )
  log_path <- file.path(log_root, sprintf("universal_%s_rep%02d.log", scenario$size_label, rep_id))
  message("[", format(Sys.time(), "%H:%M:%S"), "] Running ", scenario$size_label,
          " benchmark: ", scenario$n_sites, " sites x ", scenario$n_species,
          " species, rep ", rep_id)
  res <- run_system(rscript, runner_args, log_path)
  master_dir <- extract_master_dir(res$lines)
  if (is.na(master_dir) || !nzchar(master_dir) || !dir.exists(master_dir)) {
    warning("No master directory recovered for ", scenario$size_label, " rep ", rep_id, "; see ", log_path)
    return(list(
      effect = tibble(),
      prediction = tibble(),
      association = tibble(),
      status = tibble(),
      log = data.frame(
        size_order = scenario$size_order,
        size_label = scenario$size_label,
        n_sites = scenario$n_sites,
        n_species = scenario$n_species,
        replicate = rep_id,
        seed = run_seed,
        runner_status = res$status,
        master_dir = NA_character_,
        runner_log = log_path,
        figure_status = NA_integer_,
        figure_log = NA_character_,
        stringsAsFactors = FALSE
      )
    ))
  }
  metrics <- make_figure_for_master(master_dir, scenario, rep_id)
  metrics$log <- data.frame(
    size_order = scenario$size_order,
    size_label = scenario$size_label,
    n_sites = scenario$n_sites,
    n_species = scenario$n_species,
    replicate = rep_id,
    seed = run_seed,
    runner_status = res$status,
    master_dir = normalizePath(master_dir, winslash = "/", mustWork = FALSE),
    runner_log = log_path,
    figure_status = metrics$figure_status,
    figure_log = metrics$figure_log,
    stringsAsFactors = FALSE
  )
  metrics
}

collect_one_existing <- function(log_path) {
  base <- basename(log_path)
  m <- regexec("^universal_(.*)_rep([0-9]+)\\.log$", base)
  hit <- regmatches(base, m)[[1]]
  if (length(hit) != 3) return(NULL)
  size_label <- hit[[2]]
  rep_id <- as.integer(hit[[3]])
  scenario <- sizes[sizes$size_label == size_label, , drop = FALSE]
  if (nrow(scenario)) scenario <- scenario[1, , drop = FALSE]
  if (!nrow(scenario)) {
    lines0 <- readLines(log_path, warn = FALSE)
    run_line <- grep("Running .* benchmark:", lines0, value = TRUE)
    if (length(run_line)) {
      parsed <- regexec("Running (.*) benchmark: ([0-9]+) sites x ([0-9]+) species", run_line[[1]])
      fields <- regmatches(run_line[[1]], parsed)[[1]]
      if (length(fields) == 4) {
        scenario <- data.frame(
          size_order = length(unique(sizes$size_order)) + 1L,
          size_label = fields[[2]],
          n_sites = as.integer(fields[[3]]),
          n_species = as.integer(fields[[4]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!nrow(scenario)) {
    warning("Could not match existing scaling log to a scenario: ", log_path)
    return(NULL)
  }
  lines <- readLines(log_path, warn = FALSE)
  master_dir <- extract_master_dir(lines)
  if (is.na(master_dir) || !nzchar(master_dir) || !dir.exists(master_dir)) {
    warning("Existing scaling log has no usable master directory: ", log_path)
    return(NULL)
  }
  metrics <- make_figure_for_master(master_dir, scenario, rep_id)
  metrics$log <- data.frame(
    size_order = scenario$size_order,
    size_label = scenario$size_label,
    n_sites = scenario$n_sites,
    n_species = scenario$n_species,
    replicate = rep_id,
    seed = NA_integer_,
    runner_status = 0L,
    master_dir = normalizePath(master_dir, winslash = "/", mustWork = FALSE),
    runner_log = log_path,
    figure_status = metrics$figure_status,
    figure_log = metrics$figure_log,
    stringsAsFactors = FALSE
  )
  metrics
}

all_runs <- list()
if (collect_only) {
  message("Collecting existing Figure 4 scaling outputs from: ", log_root)
  logs <- list.files(log_root, pattern = "^universal_.*_rep[0-9]+\\.log$", full.names = TRUE)
  if (!length(logs)) stop("No universal_*_rep*.log files found in collect_log_root: ", log_root)
  for (log_path in logs) {
    one <- collect_one_existing(log_path)
    if (!is.null(one)) all_runs[[length(all_runs) + 1L]] <- one
  }
} else {
  for (i in seq_len(nrow(sizes))) {
    scenario <- sizes[i, , drop = FALSE]
    for (rep_id in seq_len(args$reps)) {
      all_runs[[length(all_runs) + 1L]] <- run_one_benchmark(scenario, rep_id)
    }
  }
}

if (!length(all_runs)) {
  stop("No scaling benchmark runs were available for aggregation.")
}

bind_component <- function(name) {
  rows <- lapply(all_runs, `[[`, name)
  rows <- rows[vapply(rows, nrow, integer(1)) > 0]
  if (!length(rows)) tibble() else bind_rows(rows)
}

effect_all <- bind_component("effect")
prediction_all <- bind_component("prediction")
association_all <- bind_component("association")
status_all <- bind_component("status")
run_log <- bind_rows(lapply(all_runs, `[[`, "log"))

readr::write_csv(effect_all, file.path(source_data_dir, "figure4_statistical_scaling_effect_metrics.csv"))
readr::write_csv(prediction_all, file.path(source_data_dir, "figure4_statistical_scaling_prediction_metrics.csv"))
readr::write_csv(association_all, file.path(source_data_dir, "figure4_statistical_scaling_association_metrics.csv"))
readr::write_csv(status_all, file.path(source_data_dir, "figure4_statistical_scaling_status.csv"))
readr::write_csv(run_log, file.path(source_data_dir, "figure4_statistical_scaling_run_log.csv"))
readr::write_csv(run_log, file.path(source_data_dir, "figure4_runtime_scaling_run_log.csv"))
readr::write_csv(effect_all, file.path(output_root, "figure4_statistical_scaling_effect_metrics.csv"))
readr::write_csv(prediction_all, file.path(output_root, "figure4_statistical_scaling_prediction_metrics.csv"))
readr::write_csv(association_all, file.path(output_root, "figure4_statistical_scaling_association_metrics.csv"))
readr::write_csv(status_all, file.path(output_root, "figure4_statistical_scaling_status.csv"))
readr::write_csv(run_log, file.path(output_root, "figure4_statistical_scaling_run_log.csv"))
readr::write_csv(run_log, file.path(output_root, "figure4_runtime_scaling_run_log.csv"))

runtime_scaling <- status_all %>%
  mutate(missing_count = if ("missing_count" %in% names(.)) missing_count else NA_real_) %>%
  transmute(
    size_order,
    size_label,
    n_sites,
    n_species,
    replicate,
    engine,
    status,
    runtime_seconds = suppressWarnings(as.numeric(runtime_seconds)),
    file_count = suppressWarnings(as.numeric(file_count)),
    zip_size_kb = suppressWarnings(as.numeric(zip_size_kb)),
    missing_count = suppressWarnings(as.numeric(missing_count)),
    master_dir
  )
readr::write_csv(runtime_scaling, file.path(source_data_dir, "figure4_runtime_scaling_metrics.csv"))
readr::write_csv(runtime_scaling, file.path(output_root, "figure4_runtime_scaling_metrics.csv"))

figure_master <- trimws(args$figure_master %||% "")
if (!nzchar(figure_master)) {
  successful <- run_log %>% filter(!is.na(master_dir), nzchar(master_dir))
  largest <- successful %>% arrange(desc(size_order), desc(replicate)) %>% slice(1)
  figure_master <- if (nrow(largest)) largest$master_dir[[1]] else NA_character_
}
if (!is.na(figure_master) && nzchar(figure_master) && dir.exists(figure_master)) {
  final_log <- file.path(log_root, "make_figure4_final.log")
  final_res <- run_system(
    rscript,
    c(
      normalizePath(figure_script, winslash = "/", mustWork = TRUE),
      paste0("--project_root=", project_root),
      paste0("--master_dir=", normalizePath(figure_master, winslash = "/", mustWork = FALSE))
    ),
    final_log
  )
} else {
  final_log <- NA_character_
  final_res <- list(status = 1L)
}

report <- c(
  "# Figure 4 Statistical Scaling Benchmark",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Project root: ", project_root),
  paste0("Scenarios: ", paste(sprintf("%s=%sx%s", sizes$size_label, sizes$n_sites, sizes$n_species), collapse = "; ")),
  paste0("Replicates per scenario: ", args$reps),
  paste0("Engines: ", args$engines),
  paste0("Truth profile: ", args$truth_profile),
  paste0("Truth parameters: effect_scale=", args$effect_scale,
         ", association_strength=", args$association_strength,
         ", residual_sd=", args$residual_sd,
         ", observation_noise_sd=", args$observation_noise_sd),
  "",
  "## Sampler Settings",
  "",
  paste(sprintf("- `%s`: %s", sampler_keys, unlist(args[sampler_keys], use.names = FALSE)), collapse = "\n"),
  "",
  "## Outputs",
  "",
  paste0("- Source data: ", source_data_dir),
  paste0("- Logs: ", log_root),
  paste0("- Figure master: ", figure_master),
  paste0("- Final Figure 4 status: ", final_res$status),
  "",
  "## Interpretation Note",
  "",
  "Effect and association panels now test recovery across increasing benchmark sizes. Panel A uses continuous predictors for the main cross-engine recovery score; categorical substrate contrasts remain exported in the standard tables and truth files, but are not mixed into the plotted recovery metric because factor baselines can differ across engines. A monotonically improving trend is expected only when the engine-specific estimand, coefficient scale, data-generating truth and exported standard table are compatible; prediction can improve even when parameter recovery remains weak."
)
write_lines(report, file.path(output_root, "Figure4_statistical_scaling_benchmark_report.md"))

message("FIGURE4_STATISTICAL_SCALING_DIR=", output_root)
message("SCALING_LOG_ROOT=", log_root)
message("FIGURE_MASTER=", figure_master)
message("MAKE_FIGURE4_EXIT_STATUS=", final_res$status)
