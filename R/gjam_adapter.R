# Real GJAM adapter for JSDM Studio.
# This file is sourced by app.R and used by automated GJAM examples.

parse_gjam_type_vector <- function(Y = NULL, type_table = NULL, type_text = "", single_type = "DA") {
  if (nzchar(trimws(type_text %||% ""))) {
    type_vec <- toupper(trimws(unlist(strsplit(type_text, ","))))
  } else if (!is.null(type_table) && ncol(type_table) >= 1) {
    cn <- tolower(colnames(type_table))
    idx <- which(cn %in% c("typename", "typenames", "type", "types"))[1]
    if (is.na(idx)) idx <- ncol(type_table)
    type_vec <- toupper(trimws(as.character(type_table[[idx]])))
  } else {
    type_vec <- toupper(trimws(single_type %||% "DA"))
  }
  type_vec <- type_vec[nzchar(type_vec)]
  if (length(type_vec) == 1 && !is.null(Y)) type_vec <- rep(type_vec, ncol(Y))
  type_vec
}

parse_gjam_group_vector <- function(text = "", S = 0, type_vec = NULL, target = "FC") {
  txt <- trimws(text %||% "")
  if (nzchar(txt)) {
    out <- suppressWarnings(as.integer(trimws(unlist(strsplit(txt, ",")))))
  } else {
    out <- rep(0L, S)
    if (!is.null(type_vec) && any(type_vec == target)) out[type_vec == target] <- 1L
  }
  out
}

validate_gjam_full <- function(Y, X, type_table = NULL, type_text = "", single_type = "DA",
                               fcgroups = "", ccgroups = "", ng = 2000, burnin = 500,
                               holdoutN = 0, random = "", notStandard = "", formula_text = "~ .",
                               censor = NULL, effort = NULL, newdata = NULL, specByTrait = NULL,
                               traitTypes = NULL, holdoutIndex = NULL, prior_file = NULL,
                               use_censor = FALSE, use_effort = FALSE, do_traits = FALSE,
                               trimY = FALSE, trim_minObs = 5, REDUCT = FALSE,
                               reduct_N = 20, reduct_r = 3) {
  ok <- TRUE
  msg <- character()
  allowed <- c("PA", "CON", "CA", "DA", "FC", "CC", "OC", "CAT")
  if (is.null(Y)) { ok <- FALSE; msg <- c(msg, "Y.csv is required and could not be read.") }
  if (is.null(X)) { ok <- FALSE; msg <- c(msg, "XData.csv is required and could not be read.") }
  if (!is.null(Y) && !is.null(X) && nrow(Y) != nrow(X)) {
    ok <- FALSE
    msg <- c(msg, paste0("Y and XData row numbers differ: Y=", nrow(Y), ", XData=", nrow(X), "."))
  }
  if (!is.null(X)) {
    X_clean <- clean_predictor_types(X, "GJAM XData")
    factor_cols <- names(X_clean)[vapply(X_clean, is.factor, logical(1))]
    if (length(factor_cols) > 0) msg <- c(msg, paste0("GJAM categorical predictors will be converted to factors before fitting: ", paste(factor_cols, collapse = ", ")))
    ftxt <- trimws(formula_text %||% "~ .")
    if (identical(ftxt, "~ .")) {
      xvars <- setdiff(names(X_clean), "intercept")
      msg <- c(msg, paste0("GJAM formula '~ .' will be expanded before fitting to avoid gjam() formula-dot failures: ~ ", paste(xvars, collapse = " + ")))
    } else {
      fcheck <- validate_one_sided_formula(ftxt, X_clean, "GJAM formula")
      ok <- ok && isTRUE(fcheck$ok)
      msg <- c(msg, fcheck$messages)
    }
    random <- trimws(random %||% "")
    if (nzchar(random)) {
      if (!(random %in% names(X_clean))) {
        ok <- FALSE
        msg <- c(msg, paste0("random column '", random, "' is not present in XData."))
      } else {
        tab <- table(as.factor(X_clean[[random]]))
        if (length(tab) < 2) msg <- c(msg, "GJAM random effects should have at least two groups.")
        if (any(tab < 2)) msg <- c(msg, "One or more GJAM random-effect groups has fewer than two observations; fitting can be unstable.")
      }
    }
    if (nzchar(notStandard %||% "")) {
      ns <- trimws(unlist(strsplit(notStandard, ",")))
      ns <- ns[nzchar(ns)]
      bad <- setdiff(ns, names(X_clean))
      if (length(bad) > 0) {
        ok <- FALSE
        msg <- c(msg, paste0("notStandard contains columns not found in XData: ", paste(bad, collapse = ", ")))
      }
      factor_bad <- intersect(ns, factor_cols)
      if (length(factor_bad) > 0) msg <- c(msg, paste0("notStandard should usually contain continuous columns, not factors: ", paste(factor_bad, collapse = ", ")))
    }
    if (!is.null(newdata)) {
      nd <- clean_predictor_types(newdata, "GJAM newdata")
      need_cols <- setdiff(names(X_clean), "intercept")
      missing_nd <- setdiff(need_cols, names(nd))
      if (length(missing_nd) > 0) {
        ok <- FALSE
        msg <- c(msg, paste0("newdata.csv is missing predictor columns used by XData/formula: ", paste(missing_nd, collapse = ", ")))
      }
    }
  }
  type_vec <- parse_gjam_type_vector(Y, type_table, type_text, single_type)
  if (!is.null(Y) && length(type_vec) != ncol(Y)) {
    ok <- FALSE
    msg <- c(msg, paste0("typeNames length must be 1 or equal to number of Y columns. Current length=", length(type_vec), ", Y columns=", ncol(Y), "."))
  }
  bad_types <- setdiff(type_vec, allowed)
  if (length(bad_types) > 0) {
    ok <- FALSE
    msg <- c(msg, paste0("Invalid GJAM typeNames: ", paste(unique(bad_types), collapse = ", "), ". Allowed: ", paste(allowed, collapse = ", ")))
  }
  if (!is.null(Y)) {
    if (length(type_vec) == ncol(Y)) {
      non_cat <- type_vec != "CAT"
      Ym <- NULL
      if (any(non_cat)) {
        ycheck <- numeric_matrix_check(Y[, non_cat, drop = FALSE], "GJAM non-CAT Y")
        ok <- ok && isTRUE(ycheck$ok)
        msg <- c(msg, ycheck$messages)
        Ym <- matrix(NA_real_, nrow = nrow(Y), ncol = ncol(Y), dimnames = list(rownames(Y), colnames(Y)))
        Ym[, non_cat] <- ycheck$matrix
      }
      for (typ in unique(type_vec)) {
        cols <- type_vec == typ
        if (identical(typ, "CA") && any(Ym[, cols, drop = FALSE] < 0, na.rm = TRUE)) {
          ok <- FALSE
          msg <- c(msg, "GJAM CA continuous-abundance responses must be non-negative.")
        } else if (identical(typ, "CON")) {
          if (any(!is.finite(Ym[, cols, drop = FALSE]) & !is.na(Ym[, cols, drop = FALSE]))) {
            ok <- FALSE
            msg <- c(msg, "GJAM CON responses must be finite numeric values or NA.")
          }
        } else if (identical(typ, "FC")) {
          vals <- Ym[, cols, drop = FALSE]
          if (any(vals < 0 | vals > 1, na.rm = TRUE)) {
            ok <- FALSE
            msg <- c(msg, "GJAM FC fractional-composition responses must be between 0 and 1.")
          }
        } else if (identical(typ, "OC")) {
          vals <- Ym[, cols, drop = FALSE]
          if (any(vals < 0, na.rm = TRUE) || any(abs(vals - round(vals)) > 1e-8, na.rm = TRUE)) {
            ok <- FALSE
            msg <- c(msg, "GJAM OC ordinal-count responses must be non-negative integer levels.")
          }
        } else if (identical(typ, "CAT")) {
          vals <- as.data.frame(Y[, cols, drop = FALSE], stringsAsFactors = FALSE)
          bad_empty <- vapply(vals, function(z) any(!is.na(z) & !nzchar(trimws(as.character(z)))), logical(1))
          if (any(bad_empty)) {
            ok <- FALSE
            msg <- c(msg, "GJAM CAT categorical responses contain empty category labels.")
          }
          low_levels <- vapply(vals, function(z) length(unique(na.omit(as.character(z)))) < 2, logical(1))
          if (any(low_levels)) {
            msg <- c(msg, paste0("GJAM CAT columns should have at least two observed categories: ", paste(names(vals)[low_levels], collapse = ", ")))
          }
        } else {
          fam_check <- response_family_messages(Ym[, cols, drop = FALSE], typ, paste("GJAM", typ))
          ok <- ok && isTRUE(fam_check$ok)
          msg <- c(msg, fam_check$messages)
        }
      }
      fc <- parse_gjam_group_vector(fcgroups, ncol(Y), type_vec, "FC")
      cc <- parse_gjam_group_vector(ccgroups, ncol(Y), type_vec, "CC")
      if ("FC" %in% type_vec && length(fc) != ncol(Y)) { ok <- FALSE; msg <- c(msg, "FCgroups must have one integer value per Y column.") }
      if ("CC" %in% type_vec && length(cc) != ncol(Y)) { ok <- FALSE; msg <- c(msg, "CCgroups must have one integer value per Y column.") }
      if ("FC" %in% type_vec && all(fc[type_vec == "FC"] == 0)) { ok <- FALSE; msg <- c(msg, "FC response columns need positive FCgroups values.") }
      if ("CC" %in% type_vec && all(cc[type_vec == "CC"] == 0)) { ok <- FALSE; msg <- c(msg, "CC response columns need positive CCgroups values.") }
      if (any(type_vec == "CAT") && sum(type_vec == "CAT") < 2) {
        ok <- FALSE
        msg <- c(msg, "GJAM CAT categorical response models require at least two CAT columns.")
      }
    }
  }
  ng <- as.integer(ng %||% 2000)
  burnin <- as.integer(burnin %||% 500)
  holdoutN <- as.integer(holdoutN %||% 0)
  if (!is.finite(ng) || ng < 50) { ok <- FALSE; msg <- c(msg, "GJAM ng should be at least 50 for a real sampler run.") }
  if (!is.finite(burnin) || burnin < 0 || burnin >= ng) { ok <- FALSE; msg <- c(msg, "GJAM burnin must be >= 0 and less than ng.") }
  if (!is.null(Y) && holdoutN >= nrow(Y)) { ok <- FALSE; msg <- c(msg, "holdoutN must be smaller than the number of observations.") }
  if (!is.null(holdoutIndex)) {
    idx <- suppressWarnings(as.integer(unlist(holdoutIndex)))
    idx <- idx[is.finite(idx)]
    if (!is.null(Y) && any(idx < 1 | idx > nrow(Y))) {
      ok <- FALSE
      msg <- c(msg, "holdoutIndex.csv contains row indices outside 1:nrow(Y).")
    }
  }
  if (isTRUE(use_effort) && is.null(effort)) {
    ok <- FALSE
    msg <- c(msg, "Use effort is checked, but effort.csv was not uploaded.")
  }
  if (isTRUE(use_censor) && is.null(censor)) {
    msg <- c(msg, "Use censor list is checked, but censor.csv was not uploaded. Text censor settings will be used if supplied.")
  }
  if (isTRUE(do_traits)) {
    if (is.null(specByTrait) || is.null(traitTypes)) {
      ok <- FALSE
      msg <- c(msg, "GJAM trait analysis requires specByTrait.csv and traitTypes.csv.")
    } else if (!is.null(Y) && nrow(specByTrait) != ncol(Y)) {
      ok <- FALSE
      msg <- c(msg, "specByTrait.csv rows must match Y columns/species.")
    }
  }
  if (isTRUE(trimY) && !is.null(Y) && (length(unique(type_vec)) > 1 || !identical(unique(type_vec), "CC"))) {
    msg <- c(msg, "gjamTrimY will only be applied automatically for single-type CC composition-count data. Other response types are fitted without trimming and a warning is written.")
  }
  if (isTRUE(REDUCT)) {
    if (!is.null(Y) && ncol(Y) < 5) msg <- c(msg, "REDUCT is usually useful for larger response matrices; small S can be unstable in gjam 2.7.")
    if (any(type_vec == "CAT")) msg <- c(msg, "REDUCT with CAT responses can be unstable in gjam 2.7; test without REDUCT first.")
    if (reduct_r < 1 || reduct_N < 1) { ok <- FALSE; msg <- c(msg, "reductList$N and reductList$r must be positive.") }
  }
  if (length(msg) == 0) msg <- "GJAM data and settings check passed."
  list(ok = ok, messages = msg, type_vec = type_vec)
}

write_gjam_reproducible_script <- function(outdir) {
  dir.create(file.path(outdir, "reproducible_script"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "workflow_scripts"), recursive = TRUE, showWarnings = FALSE)
  script <- c(
    "# Reproducible GJAM 2.7 analysis script generated by JSDM Studio",
    "# Run from the output folder with: Rscript reproducible_script/run_this_GJAM_analysis.R",
    "args <- commandArgs(trailingOnly = FALSE)",
    "file_arg <- '--file='",
    "script_arg <- args[startsWith(args, file_arg)]",
    "if (length(script_arg) > 0) {",
    "  script_path <- normalizePath(sub(file_arg, '', script_arg[[1]]), winslash='/', mustWork=FALSE)",
    "  setwd(dirname(dirname(script_path)))",
    "}",
    "for (d in c('models','chains','tables','predictions','plots','diagnostics','standard','report','results','sensitivity','ordination','missing_data')) dir.create(d, recursive=TRUE, showWarnings=FALSE)",
    "for (d in c('plots','results','sensitivity','ordination','missing_data')) writeLines(paste('GJAM', d, 'outputs are written here when available.'), file.path(d, paste0('README_GJAM_', d, '.txt')))",
    "`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x",
    "as_bool <- function(x, default=FALSE) { if (is.null(x) || length(x)==0 || is.na(x)) return(default); if (is.logical(x)) return(isTRUE(x)); tolower(as.character(x)) %in% c('true','t','1','yes','y') }",
    "as_num <- function(x, default) { z <- suppressWarnings(as.numeric(x)); if (!length(z) || !is.finite(z[1])) default else z[1] }",
    "as_int <- function(x, default) as.integer(round(as_num(x, default)))",
    "step_warnings <- character()",
    "append_warning <- function(msg) step_warnings <<- unique(c(step_warnings, as.character(msg)))",
    "write_status <- function(status, errors=character()) { payload <- list(engine='GJAM', status=status, warnings=as.character(step_warnings), errors=as.character(errors), time=as.character(Sys.time())); if (requireNamespace('jsonlite', quietly=TRUE)) writeLines(jsonlite::toJSON(payload, pretty=TRUE, auto_unbox=TRUE), file.path('diagnostics','engine_status.json')); write.csv(data.frame(engine='GJAM', status=status, message=paste(c(step_warnings, errors), collapse='; ')), file.path('tables','engine_status.csv'), row.names=FALSE); try(writeLines(capture.output(sessionInfo()), file.path('diagnostics','session_info.txt')), silent=TRUE) }",
    "safe_step <- function(name, expr, required=FALSE) { tryCatch(withCallingHandlers(force(expr), warning=function(w){ append_warning(paste(name, conditionMessage(w), sep=': ')); invokeRestart('muffleWarning') }), error=function(e){ msg <- conditionMessage(e); append_warning(paste(name, msg, sep=': ')); writeLines(msg, file.path('diagnostics', paste0('GJAM_', gsub('[^A-Za-z0-9_]+','_',name), '_error.txt'))); if (required) stop(msg, call.=FALSE); NULL }) }",
    "read_csv_safe <- function(path) { if (!file.exists(path)) return(NULL); dat <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE); if (ncol(dat) > 1) { first <- dat[[1]]; first_name <- names(dat)[1] %||% ''; first_chr <- as.character(first); first_num <- suppressWarnings(as.numeric(first_chr)); sequence_index <- all(!is.na(first_num)) && identical(as.integer(first_num), seq_len(length(first_num))); row_id_name <- first_name %in% c('', 'X', 'X.1', '...1', 'row.names', 'rowname', 'row_id', 'id', 'site_id', 'sample_id'); row_id_text <- !all(!is.na(first_num)); if (!anyDuplicated(first_chr) && (row_id_name || sequence_index || row_id_text)) { dat <- dat[-1]; rownames(dat) <- make.unique(first_chr) } }; dat }",
    "clean_df <- function(df) { if (is.null(df)) return(NULL); df <- as.data.frame(df, check.names=FALSE, stringsAsFactors=FALSE); for (nm in names(df)) { if (is.character(df[[nm]])) { x <- trimws(df[[nm]]); x[x==''] <- NA; nx <- suppressWarnings(as.numeric(x)); df[[nm]] <- if (all(is.na(x) | !is.na(nx))) nx else factor(x, levels=unique(x[!is.na(x)])) } else if (is.logical(df[[nm]])) df[[nm]] <- factor(df[[nm]]) }; df }",
    "response_df_typed <- function(dat, type_vec) { if (is.null(dat)) stop('Y.csv is missing.', call.=FALSE); out <- as.data.frame(dat, check.names=FALSE, stringsAsFactors=FALSE); if (is.null(colnames(out))) colnames(out) <- paste0('y_', seq_len(ncol(out))); if (is.null(rownames(out))) rownames(out) <- paste0('obs_', seq_len(nrow(out))); for (j in seq_along(out)) { if (!identical(unname(type_vec[j]), 'CAT')) { z <- suppressWarnings(as.numeric(out[[j]])); if (any(is.na(z) & !is.na(out[[j]]) & nzchar(trimws(as.character(out[[j]]))))) stop('Non-CAT GJAM response column is not numeric: ', names(out)[j], call.=FALSE); out[[j]] <- z } else { out[[j]] <- as.character(out[[j]]) } }; out }",
    "expand_formula <- function(ftxt, xdata) { ftxt <- trimws(ftxt %||% '~ .'); if (identical(ftxt, '~ .')) { vars <- setdiff(names(xdata), 'intercept'); if (!length(vars)) return(as.formula('~ 1')); return(as.formula(paste('~', paste(vars, collapse=' + ')))) }; as.formula(ftxt) }",
    "parse_type_vec <- function(Y, type_table, type_text, single_type) { if (nzchar(trimws(type_text %||% ''))) tv <- toupper(trimws(unlist(strsplit(type_text, ',')))) else if (!is.null(type_table) && ncol(type_table) >= 1) { cn <- tolower(names(type_table)); idx <- which(cn %in% c('typename','typenames','type','types'))[1]; if (is.na(idx)) idx <- ncol(type_table); tv <- toupper(trimws(as.character(type_table[[idx]]))) } else tv <- toupper(single_type %||% 'DA'); tv <- tv[nzchar(tv)]; if (length(tv) == 1) tv <- rep(tv, ncol(Y)); tv }",
    "parse_group_vec <- function(text, S, type_vec, target) { txt <- trimws(text %||% ''); if (nzchar(txt)) out <- suppressWarnings(as.integer(trimws(unlist(strsplit(txt, ','))))) else { out <- rep(0L, S); out[type_vec == target] <- 1L }; out }",
    "write_any_csv <- function(x, path) { if (is.null(x)) return(FALSE); if (is.atomic(x) && is.null(dim(x))) x <- data.frame(name=names(x) %||% paste0('v', seq_along(x)), value=as.vector(x)); if (is.matrix(x) || is.data.frame(x)) write.csv(as.data.frame(x, check.names=FALSE), path) else write.csv(data.frame(field=names(x), value=vapply(x, function(z) paste(capture.output(str(z, max.level=1)), collapse=' '), character(1))), path, row.names=FALSE); TRUE }",
    "long_matrix <- function(mat, value_name='estimate') { if (is.null(mat) || !is.matrix(mat)) return(data.frame()); if (is.null(rownames(mat))) rownames(mat) <- paste0('row_', seq_len(nrow(mat))); if (is.null(colnames(mat))) colnames(mat) <- paste0('col_', seq_len(ncol(mat))); grid <- expand.grid(row_id=rownames(mat), column_id=colnames(mat), KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE); data.frame(grid, value=as.vector(mat), stringsAsFactors=FALSE) }",
    "standard_predictions <- function(obs, pred) { if (is.null(pred) || !is.matrix(pred)) return(data.frame(engine='GJAM', site_id=NA_character_, response_id=NA_character_, observed=NA_real_, predicted_mean=NA_real_, predicted_lower=NA_real_, predicted_upper=NA_real_)); if (is.null(rownames(pred))) rownames(pred) <- paste0('obs_', seq_len(nrow(pred))); if (is.null(colnames(pred))) colnames(pred) <- paste0('y_', seq_len(ncol(pred))); grid <- expand.grid(site_id=rownames(pred), response_id=colnames(pred), KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE); o <- if (!is.null(obs) && all(dim(obs) == dim(pred))) as.vector(as.matrix(obs)) else NA_real_; data.frame(engine='GJAM', site_id=grid$site_id, response_id=grid$response_id, observed=o, predicted_mean=as.vector(pred), predicted_lower=NA_real_, predicted_upper=NA_real_, stringsAsFactors=FALSE) }",
    "fit_metrics <- function(fit) { if (is.null(fit)) return(data.frame(engine='GJAM', metric='fit_available', response_id=NA_character_, value=0, notes='missing fit')); rows <- list(); for (nm in names(fit)) { z <- fit[[nm]]; if (is.atomic(z)) { vals <- as.numeric(z); ids <- names(z) %||% rep(NA_character_, length(vals)); rows[[length(rows)+1]] <- data.frame(engine='GJAM', metric=nm, response_id=ids, value=vals, notes='gjam fit', stringsAsFactors=FALSE) } }; if (!length(rows)) data.frame(engine='GJAM', metric='fit_available', response_id=NA_character_, value=1, notes='non-scalar fit') else do.call(rbind, rows) }",
    "assoc_long <- function(mat, type='corMu') { if (is.null(mat) || !is.matrix(mat) || nrow(mat) < 2) return(data.frame(engine='GJAM', response_1=NA_character_, response_2=NA_character_, association_type=type, estimate=NA_real_, comparable_level='not_available')); ids <- which(upper.tri(mat), arr.ind=TRUE); data.frame(engine='GJAM', response_1=rownames(mat)[ids[,1]], response_2=colnames(mat)[ids[,2]], association_type=type, estimate=mat[ids], comparable_level='observation_scale', stringsAsFactors=FALSE) }",
    "make_censor <- function(Y, cfg, censor_raw, type_vec) { if (!as_bool(cfg$model$use_censor, FALSE)) return(list(Y=Y, censor=NULL)); cols_txt <- trimws(cfg$priors_censoring$censor_columns %||% ''); vals_txt <- trimws(cfg$priors_censoring$censor_values %||% ''); ints_txt <- trimws(cfg$priors_censoring$censor_intervals %||% ''); cols <- if (nzchar(cols_txt)) trimws(unlist(strsplit(cols_txt, ','))) else character(); values <- if (nzchar(vals_txt)) suppressWarnings(as.numeric(trimws(unlist(strsplit(vals_txt, ','))))) else numeric(); bounds <- if (nzchar(ints_txt)) suppressWarnings(as.numeric(trimws(unlist(strsplit(gsub(';', ',', ints_txt), ','))))) else numeric(); if (!is.null(censor_raw)) { cn <- tolower(names(censor_raw)); if (!length(cols) && any(cn %in% c('column','columns','response'))) cols <- as.character(censor_raw[[which(cn %in% c('column','columns','response'))[1]]]); if (!length(values) && any(cn %in% c('value','values'))) values <- as.numeric(censor_raw[[which(cn %in% c('value','values'))[1]]]); if (!length(bounds) && all(c('lower','upper') %in% cn)) bounds <- as.numeric(rbind(censor_raw[[which(cn=='lower')[1]]], censor_raw[[which(cn=='upper')[1]]])) }; if (!length(cols)) cols <- colnames(Y)[type_vec %in% c('CA','CON')][1]; whichcol <- match(cols, colnames(Y)); whichcol <- whichcol[is.finite(whichcol) & !is.na(whichcol)]; if (!length(whichcol) || !length(values)) { append_warning('censor requested but columns/values were incomplete; censor skipped'); return(list(Y=Y, censor=NULL)) }; if (!length(bounds)) bounds <- c(-Inf, values[1]); intervals <- matrix(bounds, nrow=2); if (ncol(intervals) != length(values)) intervals <- matrix(rep(intervals[,1], length(values)), nrow=2); typ <- unname(type_vec[whichcol[1]]); g <- gjam::gjamCensorY(values=values, intervals=intervals, y=as.matrix(Y[, whichcol, drop=FALSE]), type=typ, whichcol=seq_along(whichcol)); Y2 <- Y; Y2[, whichcol] <- as.data.frame(g$y, check.names=FALSE); list(Y=Y2, censor=g$censor) }",
    "make_effort <- function(Y, type_vec, effort_raw) { if (is.null(effort_raw)) return(NULL); vals <- suppressWarnings(as.numeric(as.matrix(effort_raw))); if (!length(vals) || all(is.na(vals))) return(NULL); cols <- colnames(Y)[type_vec == 'DA']; if (!length(cols)) cols <- colnames(Y); if (length(vals) == nrow(Y)) list(columns=cols, values=vals) else if (length(vals) == length(Y)) list(columns=colnames(Y), values=matrix(vals, nrow=nrow(Y), ncol=ncol(Y), dimnames=dimnames(as.matrix(Y)))) else list(columns=cols, values=rep(vals[1], nrow(Y))) }",
    "try({ if (!requireNamespace('yaml', quietly=TRUE)) stop('yaml is required.'); cfg <- yaml::read_yaml('used_config.yml') }, silent=FALSE)",
    "tryCatch({",
    "  if (!requireNamespace('gjam', quietly=TRUE)) stop('Package gjam is required for real fitting.', call.=FALSE)",
    "  suppressPackageStartupMessages(library(gjam))",
    "  set.seed(as_int(cfg$model$seed, 1234))",
    "  Y_raw <- read_csv_safe(file.path('data','Y.csv'))",
    "  X <- clean_df(read_csv_safe(file.path('data','XData.csv')))",
    "  if (is.null(Y_raw)) stop('Y.csv is missing.', call.=FALSE)",
    "  if (is.null(X)) stop('XData.csv is missing.', call.=FALSE)",
    "  if (nzchar(cfg$model$random %||% '') && cfg$model$random %in% names(X)) X[[cfg$model$random]] <- factor(X[[cfg$model$random]])",
    "  type_table <- read_csv_safe(file.path('data','typeNames.csv')); type_vec <- parse_type_vec(Y_raw, type_table, cfg$response_types$typeNames_text %||% cfg$model$typeNames_text %||% '', cfg$model$type_single %||% 'DA')",
    "  if (length(type_vec) != ncol(Y_raw)) stop('typeNames length does not match Y columns.', call.=FALSE)",
    "  names(type_vec) <- colnames(Y_raw)",
    "  Y <- response_df_typed(Y_raw, type_vec)",
    "  if (nrow(Y) != nrow(X)) stop('Y.csv and XData.csv row counts differ.', call.=FALSE)",
    "  rownames(X) <- rownames(Y)",
    "  if (as_bool(cfg$response_types$trimY, FALSE) && length(unique(type_vec)) == 1 && identical(unname(type_vec[1]), 'CC')) { tr <- safe_step('gjamTrimY', gjam::gjamTrimY(as.matrix(Y), minObs=as_int(cfg$response_types$trim_minObs, 5))); if (!is.null(tr) && !is.null(tr$y)) { Y <- as.data.frame(tr$y, check.names=FALSE); type_vec <- rep(type_vec[1], ncol(Y)); names(type_vec) <- colnames(Y); write.csv(data.frame(original_column=names(tr$nobs), nobs=as.numeric(tr$nobs)), file.path('tables','trimY_nobs.csv'), row.names=FALSE) } } else if (as_bool(cfg$response_types$trimY, FALSE)) append_warning('gjamTrimY skipped: gjam 2.7 trimming is only applied automatically for single-type CC composition-count data')",
    "  censor_raw <- read_csv_safe(file.path('data','censor.csv')); cobj <- make_censor(Y, cfg, censor_raw, type_vec); Y <- cobj$Y",
    "  effort <- if (as_bool(cfg$model$use_effort, FALSE)) make_effort(Y, type_vec, read_csv_safe(file.path('data','effort.csv'))) else NULL",
    "  form <- expand_formula(cfg$model$formula %||% '~ .', X)",
    "  notStandard <- trimws(unlist(strsplit(cfg$model$notStandard %||% '', ','))); notStandard <- notStandard[nzchar(notStandard)]",
    "  modelList <- list(ng=as_int(cfg$model$ng, 2000), burnin=as_int(cfg$model$burnin, 500), typeNames=type_vec, holdoutN=as_int(cfg$model$holdoutN, 0), FULL=as_bool(cfg$model$FULL, FALSE), PREDICTX=as_bool(cfg$model$PREDICTX, TRUE), REDUCT=as_bool(cfg$model$REDUCT, FALSE), ematAlpha=as_num(cfg$model$ematAlpha, 0.5))",
    "  if (length(notStandard)) modelList$notStandard <- notStandard",
    "  if (nzchar(cfg$model$random %||% '')) modelList$random <- cfg$model$random",
    "  if (as_bool(cfg$model$REDUCT, FALSE)) modelList$reductList <- list(N=as_int(cfg$model$reductList$N, 20), r=as_int(cfg$model$reductList$r, 3))",
    "  fc <- parse_group_vec(cfg$response_types$FCgroups %||% '', ncol(Y), type_vec, 'FC'); cc <- parse_group_vec(cfg$response_types$CCgroups %||% '', ncol(Y), type_vec, 'CC')",
    "  if ('FC' %in% type_vec) { modelList$FCgroups <- fc; attr(modelList$typeNames, 'FCgroups') <- fc }",
    "  if ('CC' %in% type_vec) { modelList$CCgroups <- cc; attr(modelList$typeNames, 'CCgroups') <- cc }",
    "  holdout_raw <- read_csv_safe(file.path('data','holdoutIndex.csv')); if (!is.null(holdout_raw)) { hi <- suppressWarnings(as.integer(unlist(holdout_raw))); hi <- hi[is.finite(hi) & hi >= 1 & hi <= nrow(Y)]; if (length(hi)) { modelList$holdoutIndex <- unique(hi); modelList$holdoutN <- length(unique(hi)) } }",
    "  if (!is.null(cobj$censor)) modelList$censor <- cobj$censor",
    "  if (!is.null(effort)) modelList$effort <- effort",
    "  if (as_bool(cfg$priors_censoring$use_prior_template, FALSE)) { pmode <- cfg$priors_censoring$prior_mode %||% 'non-informative/default'; if (identical(pmode, 'sign-constrained')) { if ('CAT' %in% type_vec) append_warning('sign-constrained betaPrior skipped for CAT responses because gjam 2.7 expands categorical columns internally') else { pred_names <- setdiff(names(X), 'intercept'); lo <- as.list(rep(0, length(pred_names))); names(lo) <- pred_names; bp <- safe_step('gjamPriorTemplate', gjam::gjamPriorTemplate(form, X, Y, lo=lo)); if (!is.null(bp)) modelList$betaPrior <- bp } } else append_warning('prior template requested but non-informative/custom prior is recorded without betaPrior injection') }",
    "  write.csv(Y, file.path('tables','Y_used.csv')); write.csv(X, file.path('tables','XData_used.csv')); write.csv(data.frame(response=colnames(Y), typeName=as.character(type_vec), FCgroups=fc, CCgroups=cc), file.path('tables','typeNames_used.csv'), row.names=FALSE); saveRDS(modelList, file.path('models','gjam_modelList.rds'))",
    "  output <- gjam::gjam(form, xdata=X, ydata=Y, modelList=modelList)",
    "  if (as_bool(cfg$outputs$save_model, TRUE)) saveRDS(output, file.path('models','gjam_model.rds'))",
    "  saveRDS(output$modelList, file.path('models','gjam_modelList_fitted.rds')); saveRDS(output$inputs, file.path('models','gjam_inputs.rds'))",
    "  if (as_bool(cfg$outputs$save_chains, TRUE) && !is.null(output$chains)) { for (nm in names(output$chains)) saveRDS(output$chains[[nm]], file.path('chains', paste0(nm, '.rds'))); write.csv(data.frame(chain=names(output$chains), rows=vapply(output$chains, NROW, integer(1)), cols=vapply(output$chains, NCOL, integer(1))), file.path('chains','chain_manifest.csv'), row.names=FALSE) }",
    "  if (as_bool(cfg$outputs$save_parameters, TRUE) && !is.null(output$parameters)) { for (nm in names(output$parameters)) write_any_csv(output$parameters[[nm]], file.path('tables', paste0(nm, '.csv'))) }",
    "  if (as_bool(cfg$outputs$save_fit, TRUE)) { write.csv(fit_metrics(output$fit), file.path('tables','fit_DIC_rmspe_xscore_yscore.csv'), row.names=FALSE); saveRDS(output$fit, file.path('results','gjam_fit.rds')) }",
    "  if (!is.null(output$prediction)) { for (nm in names(output$prediction)) write_any_csv(output$prediction[[nm]], file.path('predictions', paste0(nm, '.csv'))); saveRDS(output$prediction, file.path('predictions','gjam_prediction_from_fit.rds')) }",
    "  pextra <- NULL; new_raw <- clean_df(read_csv_safe(file.path('data','newdata.csv'))); if (as_bool(cfg$analysis$do_predict, TRUE)) { if (!is.null(new_raw)) { pextra <- tryCatch({ rownames(new_raw) <- paste0('new_', seq_len(nrow(new_raw))); gjam::gjamPredict(output, newdata=list(xdata=new_raw, nsim=50), FULL=as_bool(cfg$model$FULL, FALSE)) }, error=function(e) { append_warning(paste('gjamPredict newdata failed; using fitted-data prediction instead', conditionMessage(e), sep=': ')); NULL }) }; if (is.null(pextra)) pextra <- safe_step('gjamPredict_in_sample', gjam::gjamPredict(output, FULL=as_bool(cfg$model$FULL, FALSE))) }; if (!is.null(pextra)) { saveRDS(pextra, file.path('predictions','gjamPredict.rds')); if (!is.null(pextra$sdList)) for (nm in names(pextra$sdList)) write_any_csv(pextra$sdList[[nm]], file.path('predictions', paste0('gjamPredict_sdList_', nm, '.csv'))); if (!is.null(pextra$piList)) for (nm in names(pextra$piList)) write_any_csv(pextra$piList[[nm]], file.path('predictions', paste0('gjamPredict_piList_', nm, '.csv'))); write_any_csv(pextra$ypredMu %||% pextra$sdList$yMu %||% NULL, file.path('predictions','gjamPredict_yMu.csv')) }",
    "  sens <- NULL; if (as_bool(cfg$analysis$do_sensitivity, TRUE)) sens <- safe_step('gjamSensitivity', gjam::gjamSensitivity(output, nsim=min(100, max(10, modelList$ng - modelList$burnin)))); if (!is.null(sens)) { write.csv(as.data.frame(sens), file.path('sensitivity','gjamSensitivity.csv')); write.csv(data.frame(predictor=colnames(sens), mean=colMeans(sens, na.rm=TRUE), sd=apply(sens,2,sd,na.rm=TRUE)), file.path('tables','sensitivity_summary.csv'), row.names=FALSE) }",
    "  ord <- NULL; if (as_bool(cfg$analysis$do_ordination, FALSE)) ord <- safe_step('gjamOrdination', gjam::gjamOrdination(output, PLOT=FALSE)); if (!is.null(ord)) { saveRDS(ord, file.path('ordination','gjamOrdination.rds')); write_any_csv(ord$eVecs, file.path('ordination','ordination_eVecs.csv')); write_any_csv(ord$eValues, file.path('ordination','ordination_eValues.csv')) }",
    "  cond <- NULL; if (as_bool(cfg$analysis$do_conditional, FALSE)) cond <- safe_step('gjamConditionalParameters', gjam::gjamConditionalParameters(output, conditionOn=colnames(Y)[1], nsim=min(200, max(20, modelList$ng - modelList$burnin)))); if (!is.null(cond)) { saveRDS(cond, file.path('tables','gjamConditionalParameters.rds')); for (nm in names(cond)) write_any_csv(cond[[nm]], file.path('tables', paste0('conditional_', nm, '.csv'))) }",
    "  iie <- NULL; if (as_bool(cfg$analysis$do_iie, FALSE)) iie <- safe_step('gjamIIE', { xv <- colMeans(output$inputs$xdata, na.rm=TRUE); gjam::gjamIIE(output, xvector=xv) }); if (!is.null(iie)) saveRDS(iie, file.path('tables','gjamIIE.rds'))",
    "  trait_out <- NULL; if (as_bool(cfg$analysis$do_traits, FALSE)) trait_out <- safe_step('gjamSpec2Trait', { sbyt <- clean_df(read_csv_safe(file.path('data','specByTrait.csv'))); tt_raw <- read_csv_safe(file.path('data','traitTypes.csv')); tt <- rep('CON', ncol(sbyt)); if (!is.null(tt_raw)) { cn <- tolower(names(tt_raw)); idx <- which(cn %in% c('typename','typenames','type','types'))[1]; if (is.na(idx)) idx <- ncol(tt_raw); tt <- as.character(tt_raw[[idx]])[seq_len(ncol(sbyt))] }; gjam::gjamSpec2Trait(pbys=as.matrix(Y), sbyt=sbyt, tTypes=tt) }); if (!is.null(trait_out)) { saveRDS(trait_out, file.path('tables','gjamSpec2Trait.rds')); write_any_csv(trait_out$plotByCWM, file.path('tables','trait_plotByCWM.csv')); write.csv(data.frame(trait=names(trait_out$traitTypes), typeName=as.character(trait_out$traitTypes)), file.path('tables','traitTypes_generated.csv'), row.names=FALSE); write_any_csv(trait_out$specByTrait, file.path('tables','trait_specByTrait_generated.csv')) }",
    "  if (as_bool(cfg$outputs$save_plots, TRUE)) { safe_step('gjamPlot', { pdf(file.path('plots','gjamPlot.pdf'), width=8, height=7); on.exit(dev.off(), add=TRUE); gjam::gjamPlot(output) }); safe_step('summary_plots', { pdf(file.path('plots','GJAM_summary_plots.pdf'), width=8, height=6); on.exit(dev.off(), add=TRUE); if (!is.null(output$parameters$corMu)) image(output$parameters$corMu, main='GJAM corMu'); if (!is.null(output$prediction$ypredMu)) hist(as.vector(output$prediction$ypredMu), main='GJAM predicted y', xlab='ypredMu'); if (!is.null(sens)) boxplot(as.data.frame(sens), main='GJAM sensitivity') }) }",
    "  pred_mat <- output$prediction$ypredMu %||% if (!is.null(pextra$sdList)) pextra$sdList$yMu else NULL",
    "  write.csv(standard_predictions(as.matrix(Y), pred_mat), file.path('standard','predictions_long.csv'), row.names=FALSE)",
    "  beta <- output$parameters$betaMuUn %||% output$parameters$betaMu; effects <- if (!is.null(beta) && is.matrix(beta)) { lg <- long_matrix(beta); data.frame(engine='GJAM', response_id=lg$column_id, predictor=lg$row_id, direction=ifelse(lg$value>0,'positive',ifelse(lg$value<0,'negative','zero')), estimate=lg$value, lower=NA_real_, upper=NA_real_, notes='GJAM beta mean', stringsAsFactors=FALSE) } else data.frame(engine='GJAM', response_id=NA_character_, predictor=NA_character_, direction=NA_character_, estimate=NA_real_, lower=NA_real_, upper=NA_real_, notes='no beta matrix', stringsAsFactors=FALSE); write.csv(effects, file.path('standard','effects_long.csv'), row.names=FALSE)",
    "  write.csv(assoc_long(output$parameters$corMu, 'corMu'), file.path('standard','associations_long.csv'), row.names=FALSE)",
    "  write.csv(fit_metrics(output$fit), file.path('standard','fit_metrics.csv'), row.names=FALSE)",
    "  write.csv(data.frame(run_id=basename(getwd()), engine='GJAM', status='fitted', n_sites=nrow(Y), n_responses=ncol(Y), typeNames=paste(type_vec, collapse=','), ng=modelList$ng, burnin=modelList$burnin, holdoutN=modelList$holdoutN, random=modelList$random %||% '', REDUCT=as_bool(modelList$REDUCT, FALSE), stringsAsFactors=FALSE), file.path('standard','run_summary.csv'), row.names=FALSE)",
    "  for (f in c('run_summary.csv','fit_metrics.csv','effects_long.csv','predictions_long.csv','associations_long.csv')) file.copy(file.path('standard', f), file.path('results', f), overwrite=TRUE)",
    "  writeLines(c('GJAM results', '============', 'Status: fitted', paste0('typeNames: ', paste(type_vec, collapse=', ')), '', 'Important files:', '- models/gjam_model.rds', '- tables/betaMu.csv, corMu.csv, sigMu.csv, fit_DIC_rmspe_xscore_yscore.csv', '- predictions/*.csv', '- sensitivity/gjamSensitivity.csv when enabled', '- diagnostics/engine_status.json', '- reproducible_script/run_this_GJAM_analysis.R'), file.path('results','README_GJAM_results.txt'))",
    "  write.csv(data.frame(file=list.files('.', recursive=TRUE), stringsAsFactors=FALSE), file.path('standard','output_manifest.csv'), row.names=FALSE)",
    "  writeLines(c('<!DOCTYPE html><html><head><meta charset=\"UTF-8\"><title>GJAM Report</title></head><body><h1>GJAM workflow report</h1>', paste0('<p>Status: fitted</p><p>typeNames: ', paste(type_vec, collapse=', '), '</p>'), '<p>Open diagnostics/engine_status.json and standard/*.csv for machine-readable outputs.</p></body></html>'), file.path('report','GJAM_report.html'))",
    "  writeLines('RUN COMPLETE', 'RUN_COMPLETE.txt')",
    "  write_status('fitted')",
    "}, error=function(e) {",
    "  msg <- conditionMessage(e)",
    "  writeLines(msg, file.path('diagnostics','GJAM_reproducible_error.txt'))",
    "  write.csv(data.frame(engine='GJAM', response_id=NA_character_, predictor=NA_character_, direction=NA_character_, estimate=NA_real_, lower=NA_real_, upper=NA_real_, notes='fit_failed'), file.path('standard','effects_long.csv'), row.names=FALSE)",
    "  write.csv(data.frame(engine='GJAM', site_id=NA_character_, response_id=NA_character_, observed=NA_real_, predicted_mean=NA_real_, predicted_lower=NA_real_, predicted_upper=NA_real_), file.path('standard','predictions_long.csv'), row.names=FALSE)",
    "  write.csv(data.frame(engine='GJAM', response_1=NA_character_, response_2=NA_character_, association_type=NA_character_, estimate=NA_real_, comparable_level='fit_failed'), file.path('standard','associations_long.csv'), row.names=FALSE)",
    "  write.csv(data.frame(engine='GJAM', metric=NA_character_, response_id=NA_character_, value=NA_real_, notes='fit_failed'), file.path('standard','fit_metrics.csv'), row.names=FALSE)",
    "  write.csv(data.frame(run_id=basename(getwd()), engine='GJAM', status='fit_failed'), file.path('standard','run_summary.csv'), row.names=FALSE)",
    "  writeLines(c('GJAM results', '============', 'Status: fit_failed', '', 'Open diagnostics/GJAM_reproducible_error.txt and diagnostics/engine_status.json first.'), file.path('results','README_GJAM_results.txt'))",
    "  write_status('fit_failed', msg)",
    "  writeLines('RUN FAILED', 'RUN_FAILED.txt')",
    "  quit(status=1, save='no')",
    "})"
  )
  writeLines(script, file.path(outdir, "reproducible_script", "run_this_GJAM_analysis.R"))
  writeLines(c(
    "# Workflow wrapper for GJAM",
    "source(file.path('reproducible_script', 'run_this_GJAM_analysis.R'))"
  ), file.path(outdir, "workflow_scripts", "run_GJAM_workflow.R"))
}
