# Hmsc-HPC adapter for JSDM Studio.
# This adapter uses the Python-native pyhmsc CPU workflow bundled in external_packages/hmsc-hpc-main.

if (!exists("%||%")) {
  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0) return(y)
    if (length(x) == 1 && is.atomic(x) && is.na(x)) return(y)
    x
  }
}

hmschpc_blank <- function(x) {
  is.null(x) || length(x) == 0 || all(is.na(x)) || !nzchar(trimws(as.character(x)[1]))
}

hmschpc_clean_table <- function(x, role = "table") {
  if (is.null(x)) return(NULL)
  x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  for (nm in names(x)) {
    if (is.character(x[[nm]])) {
      z <- trimws(x[[nm]])
      z[z == ""] <- NA
      num <- suppressWarnings(as.numeric(z))
      if (!grepl("study|design|group|random", role, ignore.case = TRUE) &&
          all(is.na(z) | !is.na(num))) {
        x[[nm]] <- num
      } else {
        x[[nm]] <- as.character(z)
      }
    }
  }
  x
}

hmschpc_backtick_names <- function(x) {
  x <- as.character(x)
  ok <- grepl("^[A-Za-z_][A-Za-z0-9_]*$", x)
  out <- x
  out[!ok] <- paste0("Q(\"", gsub("\"", "\\\\\"", out[!ok]), "\")")
  out
}

hmschpc_expand_formula <- function(formula_text, data) {
  formula_text <- trimws(as.character(formula_text %||% "~ ."))
  if (!nzchar(formula_text)) formula_text <- "~ ."
  if (!grepl("^~", formula_text)) formula_text <- paste("~", formula_text)
  if (!is.null(data) && grepl("^~\\s*\\.\\s*$", formula_text)) {
    terms <- hmschpc_backtick_names(names(data))
    formula_text <- paste("~", paste(terms, collapse = " + "))
  }
  formula_text
}

hmschpc_formula_vars <- function(formula_text) {
  formula_text <- trimws(as.character(formula_text %||% "~ ."))
  if (!nzchar(formula_text) || grepl("^~\\s*\\.\\s*$", formula_text)) return(character())
  vars <- tryCatch(all.vars(stats::as.formula(formula_text)), error = function(e) character())
  unique(vars)
}

hmschpc_parse_chain_ids <- function(text = "", chains = NULL) {
  text <- trimws(as.character(text %||% ""))
  if (!nzchar(text)) return(integer())
  ids <- suppressWarnings(as.integer(trimws(unlist(strsplit(text, ",")))))
  ids <- ids[is.finite(ids)]
  if (!is.null(chains) && is.finite(as.numeric(chains))) {
    ids <- ids[ids >= 0 & ids < as.integer(chains)]
  } else {
    ids <- ids[ids >= 0]
  }
  unique(ids)
}

hmschpc_distribution_key <- function(distribution = "poisson") {
  key <- tolower(trimws(as.character(distribution %||% "poisson")))
  if (identical(key, "bernoulli")) key <- "probit"
  if (identical(key, "gaussian")) key <- "normal"
  key
}

hmschpc_python_source_dir <- function() {
  root <- get0("app_dir", ifnotfound = normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  src <- file.path(root, "external_packages", "hmsc-hpc-main")
  normalizePath(src, winslash = "/", mustWork = FALSE)
}

hmschpc_python_bin <- function(cfg = list()) {
  runtime <- cfg$runtime %||% list()
  py <- runtime$python %||% ""
  if (!hmschpc_blank(py)) return(unname(as.character(py)[1]))
  found <- unname(Sys.which("python"))
  if (hmschpc_blank(found)) "python" else found
}

validate_hmschpc_full <- function(Y = NULL, X = NULL, Tr = NULL, study = NULL, coord = NULL,
                                  phylo_cov = NULL, phylo_tree_file = NULL, newdata = NULL,
                                  distribution = "poisson", XFormula = "~ .",
                                  use_traits = FALSE, trait_formula = "~ .",
                                  phylogeny_mode = "none",
                                  random_mode = "none", random_column = "plot",
                                  coord_x = "x", coord_y = "y",
                                  nf = 1, nfMin = 1, nfMax = 4,
                                  samples = 10, transient = 10, thin = 1,
                                  chains = 2, verbose = 5,
                                  run_sampler = TRUE, random_slope_formula = "",
                                  random_name = "plot", alpha = 1,
                                  chains_to_run = integer(), seed = 1234,
                                  precision = 64, truncated_normal_library = "tf",
                                  hmcleapfrog = 10, hmcthin = 0,
                                  python = "python", python_source = hmschpc_python_source_dir(),
                                  predictions = TRUE, diagnostics = TRUE, plots = TRUE) {
  check <- list(ok = TRUE, messages = character())
  add <- function(ok, msg) {
    check$ok <<- isTRUE(check$ok) && isTRUE(ok)
    check$messages <<- c(check$messages, msg)
  }
  if (is.null(Y)) add(FALSE, "Y.csv is required.")
  if (is.null(X)) add(FALSE, "XData.csv is required.")
  if (!is.null(Y) && !is.null(X) && nrow(Y) != nrow(X)) {
    add(FALSE, paste0("Y rows (", nrow(Y), ") must match XData rows (", nrow(X), ")."))
  }
  if (!is.null(Y)) {
    ym <- suppressWarnings(as.matrix(Y))
    storage.mode(ym) <- "numeric"
    if (any(!is.finite(ym), na.rm = TRUE)) add(FALSE, "Y contains non-numeric, Inf or NaN values.")
    if (anyDuplicated(colnames(Y))) add(FALSE, "Y species/response column names must be unique.")
    key <- hmschpc_distribution_key(distribution)
    if (key %in% c("probit", "bernoulli")) {
      add(all(ym %in% c(0, 1), na.rm = TRUE), "probit response check: Y must be 0/1.")
    } else if (key == "poisson") {
      add(all(ym >= 0 & abs(ym - round(ym)) < 1e-8, na.rm = TRUE), "poisson response check: Y must be non-negative integers.")
    } else if (key %in% c("normal", "gaussian")) {
      add(TRUE, "normal/gaussian response check: continuous numeric Y is allowed.")
    } else {
      add(FALSE, paste0("Unsupported Hmsc-HPC distribution: ", distribution, ". Use poisson, probit or normal."))
    }
  }
  if (!is.null(X)) {
    x_clean <- hmschpc_clean_table(X, "XData")
    if (anyDuplicated(names(x_clean))) add(FALSE, "XData column names must be unique.")
    vars <- hmschpc_formula_vars(XFormula)
    missing <- setdiff(vars, names(x_clean))
    add(length(missing) == 0, if (length(missing)) paste0("XFormula variables missing from XData: ", paste(missing, collapse = ", ")) else "XFormula variables are present.")
    char_cols <- names(x_clean)[vapply(x_clean, is.character, logical(1))]
    if (length(char_cols) > 0) check$messages <- c(check$messages, paste0("Categorical XData predictors will be encoded by patsy/pyhmsc: ", paste(char_cols, collapse = ", ")))
  }
  if (!is.null(newdata) && !is.null(X)) {
    vars <- hmschpc_formula_vars(XFormula)
    if (length(vars) == 0) vars <- names(hmschpc_clean_table(X, "XData"))
    missing_new <- setdiff(vars, names(newdata))
    add(length(missing_new) == 0, if (length(missing_new)) paste0("newdata.csv is missing XFormula columns: ", paste(missing_new, collapse = ", ")) else "newdata.csv contains the XFormula columns needed for prediction.")
  }
  if (isTRUE(use_traits)) {
    if (is.null(Tr)) {
      add(FALSE, "traits.csv is required when Use traits is enabled.")
    } else if (!is.null(Y) && nrow(Tr) != ncol(Y)) {
      add(FALSE, paste0("traits.csv rows (", nrow(Tr), ") must match Y species columns (", ncol(Y), ")."))
    } else {
      if (anyDuplicated(names(Tr))) add(FALSE, "traits.csv column names must be unique.")
      vars <- hmschpc_formula_vars(trait_formula)
      missing <- setdiff(vars, names(Tr))
      add(length(missing) == 0, if (length(missing)) paste0("trait_formula variables missing from traits.csv: ", paste(missing, collapse = ", ")) else "trait formula variables are present.")
    }
  }
  if (identical(phylogeny_mode, "covariance")) {
    if (is.null(phylo_cov)) {
      add(FALSE, "phylo_cov.csv is required when phylogeny mode is covariance.")
    } else if (!is.null(Y)) {
      add(nrow(phylo_cov) == ncol(phylo_cov), "phylo_cov.csv must be a square species-by-species matrix.")
      missing <- setdiff(colnames(Y), intersect(rownames(phylo_cov), colnames(phylo_cov)))
      add(length(missing) == 0, if (length(missing)) paste0("phylo_cov.csv missing species: ", paste(missing, collapse = ", ")) else "phylogenetic covariance names match Y species.")
      cov_mat <- suppressWarnings(as.matrix(phylo_cov))
      storage.mode(cov_mat) <- "numeric"
      add(all(is.finite(cov_mat), na.rm = TRUE), "phylo_cov.csv values must be numeric and finite.")
      if (nrow(cov_mat) == ncol(cov_mat)) {
        add(max(abs(cov_mat - t(cov_mat)), na.rm = TRUE) < 1e-6, "phylo_cov.csv should be symmetric.")
      }
    }
  }
  if (identical(phylogeny_mode, "newick")) {
    add(!is.null(phylo_tree_file) && file.exists(phylo_tree_file), "Newick tree file must exist when phylogeny mode is Newick.")
  }
  if (!identical(random_mode, "none")) {
    add(!hmschpc_blank(random_name), "Random level name must not be empty when a random level is selected.")
    if (is.null(study)) {
      add(FALSE, "studyDesign.csv is required for iid, spatial_full or random-slope Hmsc-HPC random levels.")
    } else if (!random_column %in% names(study)) {
      add(FALSE, paste0("studyDesign.csv is missing random column: ", random_column))
    } else if (!is.null(Y) && nrow(study) != nrow(Y)) {
      add(FALSE, "studyDesign.csv rows must match Y rows.")
    } else {
      add(TRUE, paste0("Random level column detected: ", random_column))
    }
  }
  if (identical(random_mode, "spatial_full")) {
    has_study_coords <- !is.null(study) && all(c(coord_x, coord_y) %in% names(study))
    has_coord_file <- !is.null(coord) && all(c(coord_x, coord_y) %in% names(coord)) && (!is.null(Y) && nrow(coord) == nrow(Y))
    add(has_study_coords || has_coord_file, paste0("spatial_full requires coordinate columns ", coord_x, ", ", coord_y, " in studyDesign.csv or coordinates.csv."))
    coord_source <- if (has_study_coords) study[, c(coord_x, coord_y), drop = FALSE] else if (has_coord_file) coord[, c(coord_x, coord_y), drop = FALSE] else NULL
    if (!is.null(coord_source)) {
      cx <- suppressWarnings(as.numeric(coord_source[[coord_x]]))
      cy <- suppressWarnings(as.numeric(coord_source[[coord_y]]))
      add(all(is.finite(cx)) && all(is.finite(cy)), "spatial_full coordinate columns must be numeric and finite.")
    }
    add(is.finite(as.numeric(alpha)) && as.numeric(alpha) >= 0, "spatial alpha must be a finite non-negative number.")
  }
  if (identical(random_mode, "random_slope_iid") && isTRUE(run_sampler)) {
    add(FALSE, "Hmsc-HPC random-slope native sampling is currently guarded. Set Run sampler = FALSE to compile this branch only.")
  }
  if (identical(random_mode, "random_slope_iid") && !hmschpc_blank(random_slope_formula) && !is.null(study)) {
    vars <- hmschpc_formula_vars(random_slope_formula)
    missing <- setdiff(vars, names(study))
    add(length(missing) == 0, if (length(missing)) paste0("random slope x_formula variables missing from studyDesign.csv: ", paste(missing, collapse = ", ")) else "random slope x_formula variables are present.")
  }
  add(is.finite(as.numeric(nf)) && nf >= 1, "nf must be >= 1.")
  add(is.finite(as.numeric(nfMin)) && is.finite(as.numeric(nfMax)) && nfMin <= nfMax && nf >= nfMin && nf <= nfMax, "nfMin <= nf <= nfMax must hold.")
  add(samples >= 1 && transient >= 0 && thin >= 1 && chains >= 1, "Sampler counts are valid.")
  if (length(chains_to_run) > 0) {
    add(all(chains_to_run >= 0 & chains_to_run < chains), "chains to run must be zero-based ids smaller than chains.")
  }
  add(verbose >= 1, "verbose must be >= 1 for the current Hmsc-HPC sampler; verbose=0 causes a TensorFlow modulo error.")
  add(is.finite(as.numeric(seed)) && seed >= 0, "rngseed must be a finite non-negative integer.")
  add(as.integer(precision) %in% c(32L, 64L), "floating precision must be 32 or 64.")
  add(as.character(truncated_normal_library) %in% c("tf", "tfd", "scipy"), "truncated normal library must be tf, tfd or scipy.")
  add(hmcleapfrog >= 1 && hmcthin >= 0, "hmcleapfrog must be >= 1 and hmcthin must be >= 0.")
  add(nzchar(trimws(as.character(python %||% ""))), "Python executable must not be empty.")
  add(nzchar(trimws(as.character(python_source %||% ""))) && dir.exists(python_source), "Hmsc-HPC source directory must exist.")
  add(is.logical(predictions) || predictions %in% c(TRUE, FALSE), "Predictions output flag is valid.")
  add(is.logical(diagnostics) || diagnostics %in% c(TRUE, FALSE), "Diagnostics output flag is valid.")
  add(is.logical(plots) || plots %in% c(TRUE, FALSE), "Plots output flag is valid.")
  if (length(check$messages) == 0) check$messages <- "Hmsc-HPC data and settings check passed."
  check
}

hmschpc_write_pair <- function(x, outdir, filename, row.names = TRUE) {
  if (is.null(x)) return(invisible(FALSE))
  dir.create(file.path(outdir, "data"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "inputs"), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file.path(outdir, "data", filename), row.names = row.names)
  write.csv(x, file.path(outdir, "inputs", filename), row.names = row.names)
  TRUE
}

hmschpc_write_text_pair <- function(path, outdir, filename) {
  if (is.null(path) || !file.exists(path)) return(invisible(FALSE))
  dir.create(file.path(outdir, "data"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "inputs"), recursive = TRUE, showWarnings = FALSE)
  file.copy(path, file.path(outdir, "data", filename), overwrite = TRUE)
  file.copy(path, file.path(outdir, "inputs", filename), overwrite = TRUE)
  TRUE
}

hmschpc_write_model_files <- function(outdir, cfg, Y, X, Tr = NULL, study = NULL, coord = NULL,
                                      phylo_cov = NULL, phylo_tree_file = NULL, newdata = NULL) {
  for (d in c("inputs", "data", "models", "results", "tables", "plots", "diagnostics",
              "report", "predictions", "reproducible_script", "standard", "workflow_scripts",
              "samples", "compile", "hdf5", "python_logs")) {
    dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
  }
  Y <- as.data.frame(Y, check.names = FALSE)
  X <- hmschpc_clean_table(X, "XData")
  Tr <- hmschpc_clean_table(Tr, "traits")
  study <- hmschpc_clean_table(study, "studyDesign")
  coord <- hmschpc_clean_table(coord, "coordinates")
  newdata <- hmschpc_clean_table(newdata, "newdata")
  if (is.null(rownames(Y))) rownames(Y) <- paste0("site_", seq_len(nrow(Y)))
  if (is.null(rownames(X)) || any(!nzchar(rownames(X)))) rownames(X) <- rownames(Y)
  if (!is.null(Tr) && (is.null(rownames(Tr)) || any(!nzchar(rownames(Tr))))) rownames(Tr) <- colnames(Y)
  random <- cfg$random_effects %||% list()
  random_mode <- random$mode %||% "none"
  random_column <- random$column %||% "plot"
  coord_x <- random$coord_x %||% "x"
  coord_y <- random$coord_y %||% "y"
  if (!identical(random_mode, "none")) {
    if (is.null(study)) study <- data.frame(plot = rownames(Y), stringsAsFactors = FALSE)
    if (!random_column %in% names(study)) study[[random_column]] <- rownames(Y)
    if (identical(random_mode, "spatial_full") && !all(c(coord_x, coord_y) %in% names(study)) && !is.null(coord)) {
      study[[coord_x]] <- coord[[coord_x]]
      study[[coord_y]] <- coord[[coord_y]]
    }
    rownames(study) <- rownames(Y)
  }
  hmschpc_write_pair(Y, outdir, "Y.csv", row.names = TRUE)
  hmschpc_write_pair(X, outdir, "X.csv", row.names = TRUE)
  if (!is.null(Tr)) hmschpc_write_pair(Tr, outdir, "traits.csv", row.names = TRUE)
  if (!is.null(study)) hmschpc_write_pair(study, outdir, "study_design.csv", row.names = TRUE)
  if (!is.null(coord)) hmschpc_write_pair(coord, outdir, "coordinates.csv", row.names = TRUE)
  if (!is.null(phylo_cov)) hmschpc_write_pair(phylo_cov, outdir, "phylo_cov.csv", row.names = TRUE)
  if (!is.null(newdata)) hmschpc_write_pair(newdata, outdir, "newdata.csv", row.names = TRUE)
  model <- cfg$model %||% list()
  x_formula <- hmschpc_expand_formula(model$XFormula %||% "~ .", X)
  yaml_cfg <- list(
    response = "../data/Y.csv",
    covariates = "../data/X.csv",
    formula = list(X = x_formula),
    distribution = hmschpc_distribution_key(model$distribution %||% "poisson"),
    chains = as.integer((cfg$sampler %||% list())$chains %||% 2L)
  )
  if (isTRUE(model$use_traits %||% FALSE)) {
    yaml_cfg$traits <- "../data/traits.csv"
    yaml_cfg$trait_formula <- hmschpc_expand_formula(model$trait_formula %||% "~ .", Tr)
  }
  if (identical(model$phylogeny_mode %||% "none", "covariance")) yaml_cfg$phylo_cov <- "../data/phylo_cov.csv"
  if (identical(model$phylogeny_mode %||% "none", "newick")) {
    hmschpc_write_text_pair(phylo_tree_file, outdir, "tree.nwk")
    yaml_cfg$phylo_tree <- "../data/tree.nwk"
  }
  if (!identical(random_mode, "none")) {
    yaml_cfg$study_design <- "../data/study_design.csv"
    level_name <- random$name %||% "plot"
    level <- list(
      column = random_column,
      type = if (identical(random_mode, "spatial_full")) "spatial_full" else "iid",
      nf = as.integer(random$nf %||% 1L),
      nfMin = as.integer(random$nfMin %||% random$nf %||% 1L),
      nfMax = as.integer(random$nfMax %||% max(4L, as.integer(random$nf %||% 1L)))
    )
    if (identical(random_mode, "spatial_full")) {
      level$coords <- c(coord_x, coord_y)
      level$alpha <- as.numeric(random$alpha %||% 1)
    }
    if (identical(random_mode, "random_slope_iid") && !hmschpc_blank(random$x_formula)) {
      level$x_formula <- hmschpc_expand_formula(random$x_formula, study)
    }
    yaml_cfg$random_levels <- stats::setNames(list(level), level_name)
  }
  yaml::write_yaml(yaml_cfg, file.path(outdir, "workflow_scripts", "hmschpc_model.yaml"))
  write.csv(data.frame(
    file = c("workflow_scripts/hmschpc_model.yaml", "models/compiled_model/init.json", "samples/posterior.h5"),
    meaning = c("pyhmsc model config", "compiled JSON/HDF5 model metadata", "Hmsc-HPC posterior samples"),
    stringsAsFactors = FALSE
  ), file.path(outdir, "tables", "Hmsc-HPC_result_workflow_map.csv"), row.names = FALSE)
  invisible(yaml_cfg)
}

write_hmschpc_reproducible_script <- function(outdir, cfg = list()) {
  runner_src <- file.path(get0("app_dir", ifnotfound = getwd()), "workflow_scripts", "hmschpc_runner.py")
  runner_dst <- file.path(outdir, "reproducible_script", "run_this_HmscHPC_analysis.py")
  dir.create(dirname(runner_dst), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(runner_src)) file.copy(runner_src, runner_dst, overwrite = TRUE)
  source_dir <- cfg$runtime$python_source %||% hmschpc_python_source_dir()
  py <- unname(hmschpc_python_bin(cfg))
  r_script <- c(
    "# Reproducible Hmsc-HPC CPU workflow generated by JSDM Studio",
    "args <- commandArgs(trailingOnly = FALSE)",
    "file_arg <- '--file='",
    "script_arg <- args[startsWith(args, file_arg)]",
    "if (length(script_arg) > 0) {",
    "  script_path <- normalizePath(sub(file_arg, '', script_arg[[1]]), winslash = '/', mustWork = FALSE)",
    "  setwd(dirname(dirname(script_path)))",
    "}",
    sprintf("python <- %s", deparse(py)),
    sprintf("source_dir <- %s", deparse(source_dir)),
    "runner <- file.path('reproducible_script', 'run_this_HmscHPC_analysis.py')",
    "status <- system2(python, c(runner, '--workdir', normalizePath('.', winslash='/', mustWork=TRUE), '--python-source', source_dir))",
    "if (!identical(as.integer(status), 0L)) stop('Hmsc-HPC reproducible Python workflow failed with status ', status, call. = FALSE)"
  )
  writeLines(r_script, file.path(outdir, "reproducible_script", "run_this_HmscHPC_analysis.R"))
  writeLines(c(
    "# Workflow wrapper for Hmsc-HPC",
    "source(file.path('reproducible_script', 'run_this_HmscHPC_analysis.R'))"
  ), file.path(outdir, "workflow_scripts", "run_HmscHPC_workflow.R"))
  invisible(file.path(outdir, "reproducible_script", "run_this_HmscHPC_analysis.R"))
}

run_hmschpc_workflow <- function(outdir, cfg, Y, X, Tr = NULL, study = NULL, coord = NULL,
                                 phylo_cov = NULL, phylo_tree_file = NULL, newdata = NULL,
                                 log_fun = message) {
  started <- Sys.time()
  status <- list(engine = "Hmsc-HPC", status = "fit_failed", runtime_seconds = 0,
                 warnings = character(), errors = character())
  tryCatch({
    yaml::write_yaml(cfg, file.path(outdir, "used_config.yml"))
    hmschpc_write_model_files(outdir, cfg, Y, X, Tr, study, coord, phylo_cov, phylo_tree_file, newdata)
    script <- write_hmschpc_reproducible_script(outdir, cfg)
    log_fun("Hmsc-HPC compile/sample workflow starting")
    rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    stdout <- file.path(outdir, "diagnostics", "HmscHPC_Rscript_stdout.log")
    stderr <- file.path(outdir, "diagnostics", "HmscHPC_Rscript_stderr.log")
    exit_code <- tryCatch(system2(rscript, shQuote(script), stdout = stdout, stderr = stderr), error = function(e) {
      writeLines(conditionMessage(e), stderr)
      1L
    })
    status_file <- file.path(outdir, "diagnostics", "engine_status.json")
    if (file.exists(status_file) && requireNamespace("jsonlite", quietly = TRUE)) {
      parsed <- tryCatch(jsonlite::fromJSON(status_file, simplifyVector = FALSE), error = function(e) NULL)
      if (!is.null(parsed$status)) status$status <- as.character(parsed$status)
      if (!is.null(parsed$warnings)) status$warnings <- as.character(parsed$warnings)
      if (!is.null(parsed$errors)) status$errors <- as.character(parsed$errors)
    }
    if (!identical(as.integer(exit_code), 0L) && !status$status %in% c("fitted", "model_defined")) {
      status$status <- "fit_failed"
      status$errors <- unique(c(status$errors, paste0("Rscript exited with status ", as.integer(exit_code), ".")))
    }
    status$runtime_seconds <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2)
    if (exists("write_engine_status", mode = "function")) write_engine_status(outdir, status)
    log_fun(paste("Hmsc-HPC workflow status:", status$status))
    status
  }, error = function(e) {
    status$status <- "fit_failed"
    status$runtime_seconds <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2)
    status$errors <- c(status$errors, conditionMessage(e))
    dir.create(file.path(outdir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
    writeLines(status$errors, file.path(outdir, "diagnostics", "HmscHPC_R_adapter_error.txt"))
    if (exists("write_engine_status", mode = "function")) write_engine_status(outdir, status)
    status
  })
}
