# Run small/medium/large Universal Benchmark timings for Figure 4D.
#
# The runs are intentionally sequential to avoid Windows conflicts among JAGS,
# reticulate/Python, TensorFlow/PyTorch and package-level temporary files.

options(stringsAsFactors = FALSE)

if (.Platform$OS.type == "windows") {
  Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
}

required <- c("readr", "dplyr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (length(x) == 1 && is.atomic(x) && is.na(x)) return(y)
  x
}

timestamp_id <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

parse_args <- function(args) {
  out <- list(
    project_root = "",
    figure_master = "",
    reps = 1L,
    engines = "Hmsc,Hmsc-HPC,jSDM,GJAM,spOccupancy,sjSDM,boral",
    hmsc_samples = 30L,
    hmsc_transient = 30L,
    hmsc_nChains = 2L,
    hmschpc_samples = 30L,
    hmschpc_transient = 30L,
    hmschpc_chains = 1L
  )
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "true"
    out[[key]] <- val
  }
  ints <- c("reps", "hmsc_samples", "hmsc_transient", "hmsc_nChains",
            "hmschpc_samples", "hmschpc_transient", "hmschpc_chains")
  for (nm in ints) out[[nm]] <- as.integer(out[[nm]])
  out
}

script_path <- {
  args0 <- commandArgs(trailingOnly = FALSE)
  hit <- args0[startsWith(args0, "--file=")]
  if (length(hit)) normalizePath(sub("^--file=", "", hit[[1]]), winslash = "/", mustWork = FALSE)
  else normalizePath(file.path("workflow_scripts", "run_figure4_runtime_scaling_benchmark.R"), winslash = "/", mustWork = FALSE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
app_dir <- if (nzchar(args$project_root %||% "")) {
  normalizePath(args$project_root, winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
}
if (!file.exists(file.path(app_dir, "app.R"))) stop("Cannot locate app.R for JSDM Studio.", call. = FALSE)
setwd(app_dir)

latest_master <- function() {
  dirs <- list.dirs(file.path(app_dir, "output"), recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^universal_benchmark_\\d", basename(dirs))]
  masters <- file.path(dirs, "master")
  keep <- file.exists(file.path(masters, "04_unified_standard_results", "model_status_matrix.csv"))
  masters <- masters[keep]
  if (!length(masters)) stop("No Universal Benchmark master directory found.", call. = FALSE)
  normalizePath(masters[which.max(file.info(masters)$mtime)], winslash = "/", mustWork = TRUE)
}

figure_master <- if (nzchar(args$figure_master %||% "")) {
  normalizePath(args$figure_master, winslash = "/", mustWork = TRUE)
} else {
  latest_master()
}

output_root <- file.path(app_dir, "docs", "paper_figures", "figure4_statistical_universal_benchmark")
source_data_dir <- file.path(output_root, "source_data")
log_root <- file.path(app_dir, "output", paste0("figure4_runtime_scaling_logs_", timestamp_id()))
input_root <- file.path(log_root, "inputs")
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
dir.create(input_root, recursive = TRUE, showWarnings = FALSE)

scenarios <- data.frame(
  size_order = c(1, 2, 3),
  size_label = c("Small 36x6", "Medium 60x10", "Large 90x14"),
  n_sites = c(36L, 60L, 90L),
  n_species = c(6L, 10L, 14L),
  n_visits = c(3L, 3L, 3L),
  stringsAsFactors = FALSE
)

rbin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
runner <- file.path(app_dir, "workflow_scripts", "universal_benchmark_runner.R")
fig4_script <- file.path(app_dir, "workflow_scripts", "make_figure4_statistical_universal_benchmark.R")
engines <- paste(trimws(unlist(strsplit(args$engines, ","))), collapse = ",")

parse_master_from_stdout <- function(stdout) {
  if (!file.exists(stdout)) return(NA_character_)
  lines <- readLines(stdout, warn = FALSE)
  hit <- grep("^MASTER_DIR=", lines, value = TRUE)
  if (length(hit)) sub("^MASTER_DIR=", "", hit[[length(hit)]]) else NA_character_
}

rows <- list()
run_rows <- list()
for (replicate_id in seq_len(max(1L, args$reps))) {
  for (i in seq_len(nrow(scenarios))) {
    sc <- scenarios[i, ]
    tag <- paste0("size", sc$size_order, "_rep", replicate_id)
    stdout <- file.path(log_root, paste0(tag, "_stdout.txt"))
    stderr <- file.path(log_root, paste0(tag, "_stderr.txt"))
    target_dir <- file.path(input_root, tag)
    seed <- 20260700L + sc$size_order * 100L + replicate_id
    flags <- c(
      runner,
      "--action=all",
      "--allow_large=true",
      paste0("--target_dir=", target_dir),
      paste0("--n_sites=", sc$n_sites),
      paste0("--n_species=", sc$n_species),
      paste0("--n_visits=", sc$n_visits),
      paste0("--seed=", seed),
      paste0("--engines=", engines),
      paste0("--hmsc_samples=", args$hmsc_samples),
      paste0("--hmsc_transient=", args$hmsc_transient),
      "--hmsc_thin=1",
      paste0("--hmsc_nChains=", args$hmsc_nChains),
      paste0("--hmschpc_samples=", args$hmschpc_samples),
      paste0("--hmschpc_transient=", args$hmschpc_transient),
      "--hmschpc_thin=1",
      paste0("--hmschpc_chains=", args$hmschpc_chains)
    )
    start <- Sys.time()
    code <- system2(rbin, shQuote(flags), stdout = stdout, stderr = stderr)
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    code <- if (is.null(code)) 0L else as.integer(code)
    master_dir <- parse_master_from_stdout(stdout)
    if (!is.na(master_dir) && dir.exists(master_dir)) {
      status_file <- file.path(master_dir, "04_unified_standard_results", "model_status_matrix.csv")
      status <- if (file.exists(status_file)) read.csv(status_file, check.names = FALSE, stringsAsFactors = FALSE) else data.frame()
      if (nrow(status)) {
        rows[[length(rows) + 1L]] <- status %>%
          mutate(size_order = sc$size_order, size_label = sc$size_label,
                 n_sites = sc$n_sites, n_species = sc$n_species,
                 n_visits = sc$n_visits, replicate = replicate_id,
                 benchmark_run = basename(dirname(master_dir)),
                 master_dir = normalizePath(master_dir, winslash = "/", mustWork = FALSE),
                 run_exit_status = code)
      }
    }
    run_rows[[length(run_rows) + 1L]] <- data.frame(
      size_order = sc$size_order,
      size_label = sc$size_label,
      n_sites = sc$n_sites,
      n_species = sc$n_species,
      replicate = replicate_id,
      seed = seed,
      exit_status = code,
      elapsed_seconds_total = elapsed,
      master_dir = master_dir %||% NA_character_,
      stdout = stdout,
      stderr = stderr,
      stringsAsFactors = FALSE
    )
  }
}

scaling <- if (length(rows)) bind_rows(rows) else data.frame()
run_log <- bind_rows(run_rows)
write_csv(scaling, file.path(source_data_dir, "figure4_runtime_scaling_metrics.csv"))
write_csv(run_log, file.path(source_data_dir, "figure4_runtime_scaling_run_log.csv"))
write_csv(scaling, file.path(output_root, "figure4_runtime_scaling_metrics.csv"))
write_csv(run_log, file.path(output_root, "figure4_runtime_scaling_run_log.csv"))

fig_status <- system2(rbin, c(shQuote(fig4_script), paste0("--project_root=", app_dir),
                              paste0("--master_dir=", figure_master)),
                      stdout = file.path(log_root, "make_figure4_stdout.txt"),
                      stderr = file.path(log_root, "make_figure4_stderr.txt"))
fig_status <- if (is.null(fig_status)) 0L else as.integer(fig_status)

cat("FIGURE4_SCALING_DIR=", normalizePath(output_root, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("SCALING_LOG_ROOT=", normalizePath(log_root, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("FIGURE_MASTER=", normalizePath(figure_master, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("MAKE_FIGURE4_EXIT_STATUS=", fig_status, "\n", sep = "")
