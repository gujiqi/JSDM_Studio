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

hmschpc_align_newdata <- function(newdata, X, outdir = NULL) {
  if (is.null(newdata)) return(NULL)
  out <- as.data.frame(newdata, check.names = FALSE, stringsAsFactors = FALSE)
  template <- as.data.frame(X, check.names = FALSE, stringsAsFactors = FALSE)
  if (is.null(rownames(out)) || any(!nzchar(rownames(out)))) rownames(out) <- paste0("new_", seq_len(nrow(out)))
  adjustments <- list()
  add_adj <- function(column, row_id, original, replacement, reason) {
    adjustments[[length(adjustments) + 1L]] <<- data.frame(
      column = column,
      row_id = row_id,
      original = as.character(original),
      replacement = as.character(replacement),
      reason = reason,
      stringsAsFactors = FALSE
    )
  }
  for (nm in names(template)) {
    train <- template[[nm]]
    if (!nm %in% names(out)) {
      if (is.numeric(train)) {
        repl <- suppressWarnings(stats::median(train, na.rm = TRUE))
        if (!is.finite(repl)) repl <- 0
      } else {
        vals <- unique(trimws(as.character(train)))
        vals <- vals[!is.na(vals) & nzchar(vals)]
        repl <- if (length(vals)) vals[[1]] else ""
      }
      out[[nm]] <- rep(repl, nrow(out))
      add_adj(nm, rownames(out), NA_character_, repl, "column_missing_in_newdata_filled_from_training_reference")
    }
    if (is.numeric(train)) {
      vals <- suppressWarnings(as.numeric(out[[nm]]))
      bad <- !is.finite(vals)
      if (any(bad)) {
        repl <- suppressWarnings(stats::median(train, na.rm = TRUE))
        if (!is.finite(repl)) repl <- 0
        for (rid in rownames(out)[bad]) add_adj(nm, rid, out[[nm]][rownames(out) == rid][1], repl, "non_numeric_or_missing_value_imputed_from_training_median")
        vals[bad] <- repl
      }
      out[[nm]] <- vals
    } else {
      levels_train <- unique(trimws(as.character(train)))
      levels_train <- levels_train[!is.na(levels_train) & nzchar(levels_train)]
      if (length(levels_train)) {
        vals <- trimws(as.character(out[[nm]]))
        vals[vals == ""] <- NA
        bad <- is.na(vals) | !(vals %in% levels_train)
        if (any(bad)) {
          repl <- levels_train[[1]]
          for (rid in rownames(out)[bad]) add_adj(nm, rid, out[[nm]][rownames(out) == rid][1], repl, "unseen_or_missing_factor_level_replaced_by_training_reference")
          vals[bad] <- repl
        }
        out[[nm]] <- vals
      }
    }
  }
  if (!is.null(outdir)) {
    tab <- if (length(adjustments)) do.call(rbind, adjustments) else data.frame(
      column = NA_character_, row_id = NA_character_, original = NA_character_,
      replacement = NA_character_, reason = "no_newdata_adjustments_needed",
      stringsAsFactors = FALSE
    )
    dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
    write.csv(tab, file.path(outdir, "tables", "HmscHPC_newdata_adjustments.csv"), row.names = FALSE)
  }
  out
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

hmschpc_make_positive_definite_correlation <- function(C, eps = 1e-6) {
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

hmschpc_taxonomy_to_correlation <- function(taxonomy, species_names) {
  taxonomy <- as.data.frame(taxonomy, stringsAsFactors = FALSE, check.names = FALSE)
  nm_low <- tolower(names(taxonomy))
  species_col <- names(taxonomy)[nm_low %in% c("species", "species_name", "species_id", "taxon", "taxon_id", "taxon_name", "response", "response_id")]
  rownames_match <- !is.null(rownames(taxonomy)) && length(setdiff(species_names, rownames(taxonomy))) == 0
  if (!rownames_match && length(species_col) > 0) {
    rownames(taxonomy) <- as.character(taxonomy[[species_col[1]]])
  }
  if (is.null(rownames(taxonomy)) || any(!nzchar(rownames(taxonomy)))) {
    stop("taxonomy.csv must have species names either as row names or in a species/taxon column.", call. = FALSE)
  }
  missing_species <- setdiff(species_names, rownames(taxonomy))
  if (length(missing_species) > 0) {
    stop(paste0("taxonomy.csv is missing species: ", paste(missing_species, collapse = ", ")), call. = FALSE)
  }
  taxonomy <- taxonomy[species_names, , drop = FALSE]
  rank_order <- c("kingdom", "phylum", "division", "class", "order", "family", "genus", "species")
  rank_cols <- names(taxonomy)[tolower(names(taxonomy)) %in% rank_order]
  rank_cols <- rank_cols[order(match(tolower(rank_cols), rank_order))]
  rank_cols <- setdiff(rank_cols, species_col[1] %||% character())
  if (length(rank_cols) == 0) {
    stop("taxonomy.csv must contain taxonomic-rank columns such as kingdom, phylum, class, order, family, genus or species.", call. = FALSE)
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
      C[i, j] <- if (!any(comparable)) as.numeric(i == j) else mean(unlist(ranks[i, comparable, drop = TRUE]) == unlist(ranks[j, comparable, drop = TRUE]))
    }
  }
  hmschpc_make_positive_definite_correlation(C)
}

hmschpc_newick_to_correlation <- function(phylo_tree_file, species_names) {
  if (is.null(phylo_tree_file) || !file.exists(phylo_tree_file)) return(NULL)
  if (!requireNamespace("ape", quietly = TRUE)) return(NULL)
  tree <- tryCatch(ape::read.tree(phylo_tree_file), error = function(e) NULL)
  if (is.null(tree) || is.null(tree$tip.label)) return(NULL)
  missing_species <- setdiff(species_names, tree$tip.label)
  if (length(missing_species) > 0) {
    stop(paste0("phylo_tree.nwk is missing species: ", paste(missing_species, collapse = ", ")), call. = FALSE)
  }
  extra <- setdiff(tree$tip.label, species_names)
  if (length(extra) > 0) tree <- ape::keep.tip(tree, species_names)
  C <- tryCatch(ape::vcv.phylo(tree, corr = TRUE), error = function(e) NULL)
  if (is.null(C)) return(NULL)
  C <- C[species_names, species_names, drop = FALSE]
  hmschpc_make_positive_definite_correlation(C)
}

hmschpc_write_derived_phylogeny <- function(outdir, C = NULL, source = "none", input_file = NA_character_) {
  if (is.null(C)) return(invisible(FALSE))
  C <- hmschpc_make_positive_definite_correlation(C)
  for (subdir in c("data", "tables", "standard")) {
    dir.create(file.path(outdir, subdir), recursive = TRUE, showWarnings = FALSE)
    write.csv(C, file.path(outdir, subdir, "derived_phylogeny_correlation.csv"))
  }
  write.csv(
    data.frame(
      source = source,
      input_file = input_file,
      derived_file = "derived_phylogeny_correlation.csv",
      user_upload_contract = "Upload a Newick tree or taxonomy rank table; do not upload a species-by-species covariance/correlation matrix.",
      note = "Hmsc-HPC uses this internally derived matrix where the CPU pyhmsc path requires a species relatedness matrix.",
      stringsAsFactors = FALSE
    ),
    file.path(outdir, "tables", "phylogeny_input_contract.csv"),
    row.names = FALSE
  )
  invisible(TRUE)
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
  windows <- .Platform$OS.type == "windows"
  app_python <- if (windows) file.path(get0("app_dir", ifnotfound = getwd()), "runtimes", "python310", "python.exe") else character()
  local_appdata <- Sys.getenv("LOCALAPPDATA")
  bundled <- if (windows && nzchar(local_appdata)) {
    file.path(local_appdata, "JSDMStudio", "hmschpc-python310", "Scripts", "python.exe")
  } else character()
  paths <- c(Sys.getenv("JSDM_HMSCHPC_PYTHON"), app_python, bundled, Sys.getenv("RETICULATE_PYTHON"))
  paths <- paths[nzchar(paths) & file.exists(paths)]
  if (length(paths)) {
    path <- paths[[1]]
    # Resolving a venv's Python symlink itself would select system Python.
    return(file.path(normalizePath(dirname(path), winslash = "/", mustWork = TRUE), basename(path)))
  }
  commands <- if (windows) c("python", "python3") else c("python3", "python")
  found <- unname(Sys.which(commands))
  found <- found[nzchar(found)][1]
  if (hmschpc_blank(found)) "python" else found
}

validate_hmschpc_full <- function(Y = NULL, X = NULL, Tr = NULL, study = NULL, coord = NULL,
                                  phylo_cov = NULL, phylo_tree_file = NULL, taxonomy = NULL,
                                  newdata = NULL,
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
  pairing <- c(studio_row_order_issue(Y, X, "XData.csv"),
               if (use_traits) studio_trait_order_issue(Y, Tr),
               if (random_mode != "none") studio_row_order_issue(Y, study, "studyDesign.csv"),
               if (random_mode == "spatial_full") studio_row_order_issue(Y, coord, "coordinates.csv"))
  check <- list(ok = !length(pairing), messages = pairing)
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
    ym <- as.matrix(Y)
    observed <- !is.na(ym)
    suppressWarnings(storage.mode(ym) <- "numeric")
    if (any(!is.finite(ym[observed]))) add(FALSE, "Y contains non-numeric or infinite observations.")
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
  if (identical(phylogeny_mode, "taxonomy")) {
    if (is.null(taxonomy)) {
      add(FALSE, "taxonomy.csv is required when phylogeny mode is taxonomy.")
    } else if (!is.null(Y)) {
      tax_ok <- tryCatch({
        hmschpc_taxonomy_to_correlation(taxonomy, colnames(Y))
        TRUE
      }, error = function(e) {
        add(FALSE, conditionMessage(e))
        FALSE
      })
      if (isTRUE(tax_ok)) add(TRUE, "taxonomy.csv species and taxonomic-rank columns match Y species; a derived correlation matrix can be generated internally.")
    }
  }
  if (identical(phylogeny_mode, "covariance")) {
    if (is.null(phylo_cov)) {
      add(FALSE, "Deprecated internal covariance mode requested, but no derived matrix is available. Use taxonomy.csv or phylo_tree.nwk as the user input.")
    } else if (!is.null(Y)) {
      add(nrow(phylo_cov) == ncol(phylo_cov), "Deprecated internal derived correlation matrix must be square.")
      missing <- setdiff(colnames(Y), intersect(rownames(phylo_cov), colnames(phylo_cov)))
      add(length(missing) == 0, if (length(missing)) paste0("Deprecated internal derived correlation matrix missing species: ", paste(missing, collapse = ", ")) else "Derived phylogenetic/taxonomic correlation names match Y species.")
      cov_mat <- suppressWarnings(as.matrix(phylo_cov))
      storage.mode(cov_mat) <- "numeric"
      add(all(is.finite(cov_mat), na.rm = TRUE), "Deprecated internal derived correlation values must be numeric and finite.")
      if (nrow(cov_mat) == ncol(cov_mat)) {
        add(max(abs(cov_mat - t(cov_mat)), na.rm = TRUE) < 1e-6, "Deprecated internal derived correlation matrix should be symmetric.")
      }
    }
  }
  if (identical(phylogeny_mode, "newick")) {
    add(!is.null(phylo_tree_file) && file.exists(phylo_tree_file), "Newick tree file must exist when phylogeny mode is Newick.")
    if (!is.null(phylo_tree_file) && file.exists(phylo_tree_file) && !is.null(Y)) {
      ok_tree <- tryCatch({
        hmschpc_newick_to_correlation(phylo_tree_file, colnames(Y))
        TRUE
      }, error = function(e) {
        add(FALSE, conditionMessage(e))
        FALSE
      })
      if (isTRUE(ok_tree)) add(TRUE, "Newick tree tips match Y species; a reviewer-facing derived correlation matrix can be exported.")
    }
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
      if (all(is.finite(cx)) && all(is.finite(cy))) {
        add(!any(duplicated(data.frame(x = cx, y = cy))), "spatial_full coordinate rows must not contain exact duplicates.")
      }
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
                                      phylo_cov = NULL, phylo_tree_file = NULL, taxonomy = NULL,
                                      newdata = NULL) {
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
  taxonomy <- hmschpc_clean_table(taxonomy, "taxonomy")
  newdata <- hmschpc_clean_table(newdata, "newdata")
  newdata <- hmschpc_align_newdata(newdata, X, outdir)
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
  if (!is.null(taxonomy)) hmschpc_write_pair(taxonomy, outdir, "taxonomy.csv", row.names = TRUE)
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
  phylogeny_mode <- model$phylogeny_mode %||% "none"
  if (identical(phylogeny_mode, "taxonomy")) {
    C_derived <- hmschpc_taxonomy_to_correlation(taxonomy, colnames(Y))
    hmschpc_write_derived_phylogeny(outdir, C_derived, source = "taxonomy_table", input_file = "taxonomy.csv")
    yaml_cfg$phylo_cov <- "../data/derived_phylogeny_correlation.csv"
    yaml_cfg$phylogeny_source <- "taxonomy.csv"
  }
  if (identical(phylogeny_mode, "covariance")) {
    C_derived <- hmschpc_make_positive_definite_correlation(phylo_cov)
    rownames(C_derived) <- colnames(Y)
    colnames(C_derived) <- colnames(Y)
    hmschpc_write_derived_phylogeny(outdir, C_derived, source = "deprecated_internal_matrix", input_file = "not_user_upload")
    yaml_cfg$phylo_cov <- "../data/derived_phylogeny_correlation.csv"
    yaml_cfg$phylogeny_source <- "deprecated_internal_matrix"
  }
  if (identical(model$phylogeny_mode %||% "none", "newick")) {
    hmschpc_write_text_pair(phylo_tree_file, outdir, "tree.nwk")
    C_derived <- hmschpc_newick_to_correlation(phylo_tree_file, colnames(Y))
    hmschpc_write_derived_phylogeny(outdir, C_derived, source = "newick_tree", input_file = "tree.nwk")
    yaml_cfg$phylo_cov <- "../data/derived_phylogeny_correlation.csv"
    yaml_cfg$phylogeny_source <- "tree.nwk"
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
    file = c("workflow_scripts/hmschpc_model.yaml", "data/taxonomy.csv or data/tree.nwk", "data/derived_phylogeny_correlation.csv", "models/compiled_model/init.json", "samples/posterior.h5"),
    meaning = c("pyhmsc model config", "raw phylogeny input when enabled", "internally derived species relatedness matrix", "compiled JSON/HDF5 model metadata", "Hmsc-HPC posterior samples"),
    stringsAsFactors = FALSE
  ), file.path(outdir, "tables", "Hmsc-HPC_result_workflow_map.csv"), row.names = FALSE)
  invisible(yaml_cfg)
}

hmschpc_rscript_clean_env <- function() {
  if (.Platform$OS.type == "windows") {
    Sys.unsetenv(c("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME", "LANG", "LANGUAGE"))
  }
  invisible(TRUE)
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
    "if (.Platform$OS.type == 'windows') Sys.unsetenv(c('LC_ALL','LC_COLLATE','LC_CTYPE','LC_MONETARY','LC_TIME','LANG','LANGUAGE'))",
    "args <- commandArgs(trailingOnly = FALSE)",
    "file_arg <- '--file='",
    "script_arg <- args[startsWith(args, file_arg)]",
    "if (length(script_arg) > 0) {",
    "  script_path <- sub(file_arg, '', script_arg[[1]])",
    "  if (!file.exists(script_path)) script_path <- gsub('~+~', ' ', script_path, fixed = TRUE)",
    "  setwd(dirname(dirname(normalizePath(script_path, winslash = '/', mustWork = TRUE))))",
    "}",
    sprintf("python <- %s", deparse(py)),
    sprintf("source_dir <- %s", deparse(source_dir)),
    "runner <- file.path('reproducible_script', 'run_this_HmscHPC_analysis.py')",
    "status <- system2(python, shQuote(c(runner, '--workdir', normalizePath('.', winslash='/', mustWork=TRUE), '--python-source', source_dir)))",
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
                                 phylo_cov = NULL, phylo_tree_file = NULL, taxonomy = NULL, newdata = NULL,
                                 log_fun = message) {
  started <- Sys.time()
  status <- list(engine = "Hmsc-HPC", status = "fit_failed", runtime_seconds = 0,
                 warnings = character(), errors = character())
  tryCatch({
    yaml::write_yaml(cfg, file.path(outdir, "used_config.yml"))
    hmschpc_write_model_files(outdir, cfg, Y, X, Tr, study, coord, phylo_cov, phylo_tree_file, taxonomy, newdata)
    script <- write_hmschpc_reproducible_script(outdir, cfg)
    log_fun("Hmsc-HPC compile/sample workflow starting")
    rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    stdout <- file.path(outdir, "diagnostics", "HmscHPC_Rscript_stdout.log")
    stderr <- file.path(outdir, "diagnostics", "HmscHPC_Rscript_stderr.log")
    hmschpc_rscript_clean_env()
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
