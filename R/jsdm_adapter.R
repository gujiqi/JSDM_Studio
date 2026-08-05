# Real jSDM adapter for JSDM Studio.
# This file is sourced by app.R and is also used by automated examples.

validate_jsdm_full <- function(Y = NULL, X = NULL, Tr = NULL, long = NULL, trials = NULL,
                               newdata = NULL, prediction_ids = NULL,
                               model_type = "binomial_probit", site_formula = "~ .",
                               trait_formula = "~ .", n_latent = 0, site_effect = "none",
                               burnin = 50, mcmc = 50, thin = 1, trials_scalar = 1,
                               allow_traits = TRUE, do_predict = TRUE,
                               long_site_col = "site", long_species_col = "species",
                               long_response_col = "presence", constrained_nchains = 2) {
  ok <- TRUE
  msg <- character()
  model_type <- as.character(model_type %||% "binomial_probit")
  n_latent <- as.integer(n_latent %||% 0)
  burnin <- as.integer(burnin %||% 0)
  mcmc <- as.integer(mcmc %||% 0)
  thin <- as.integer(thin %||% 1)
  if ((burnin + mcmc) < 100 || ((burnin + mcmc) %% 10) != 0) {
    ok <- FALSE
    msg <- c(msg, "jSDM requires burnin + mcmc >= 100 and divisible by 10.")
  }
  if (mcmc < 1 || thin < 1 || (mcmc %% thin) != 0) {
    ok <- FALSE
    msg <- c(msg, "jSDM requires mcmc >= 1, thin >= 1, and mcmc divisible by thin.")
  }
  if (!(site_effect %in% c("none", "fixed", "random"))) {
    ok <- FALSE
    msg <- c(msg, "jSDM site_effect must be none, fixed or random.")
  }
  if (!(model_type %in% c("binomial_probit", "binomial_logit", "poisson_log", "gaussian", "binomial_probit_long_format", "binomial_probit_sp_constrained"))) {
    ok <- FALSE
    msg <- c(msg, "Unsupported jSDM model_type.")
  }
  if (identical(model_type, "binomial_logit") && n_latent > 0 && site_effect != "none") {
    ok <- FALSE
    msg <- c(msg, "jSDM 0.2.7 binomial_logit is unstable with both n_latent > 0 and site_effect fixed/random. Use latent variables with site_effect = none, or set n_latent = 0 when using site effects.")
  }
  is_long <- identical(model_type, "binomial_probit_long_format")
  if (is_long) {
    if (is.null(long)) {
      ok <- FALSE
      msg <- c(msg, "long_format.csv is required for jSDM_binomial_probit_long_format.")
    } else {
      missing_cols <- setdiff(c(long_site_col, long_species_col, long_response_col), names(long))
      if (length(missing_cols) > 0) {
        ok <- FALSE
        msg <- c(msg, paste0("long_format.csv is missing required columns: ", paste(missing_cols, collapse = ", ")))
      } else {
        yy <- suppressWarnings(as.numeric(long[[long_response_col]]))
        if (any(is.na(yy)) || any(!(yy %in% c(0, 1)))) {
          ok <- FALSE
          msg <- c(msg, "jSDM long-format response column must be binary 0/1.")
        }
        long_check <- long
        names(long_check)[names(long_check) == long_site_col] <- "site"
        names(long_check)[names(long_check) == long_species_col] <- "species"
        names(long_check)[names(long_check) == long_response_col] <- "Y"
        fcheck <- validate_one_sided_formula(site_formula, clean_predictor_types(long_check, "jSDM long_format"), "jSDM long-format site_formula")
        ok <- ok && isTRUE(fcheck$ok)
        msg <- c(msg, fcheck$messages)
        if (!grepl("species", site_formula, fixed = TRUE)) {
          msg <- c(msg, "jSDM long-format fitting will use species-specific terms internally; for custom formulas use e.g. ~ species + species:x1 + species:x2.")
        }
      }
    }
  } else {
    if (is.null(Y) || is.null(X)) {
      ok <- FALSE
      msg <- c(msg, "Y.csv and XData.csv are required for matrix-format jSDM models.")
    } else {
      ycheck <- numeric_matrix_check(Y, "jSDM Y")
      ok <- ok && isTRUE(ycheck$ok)
      msg <- c(msg, ycheck$messages)
      Ymat <- ycheck$matrix
      if (!is.null(Ymat)) {
        if (nrow(Ymat) != nrow(X)) {
          ok <- FALSE
          msg <- c(msg, "jSDM Y.csv rows must match XData.csv rows.")
        }
        if (n_latent > 0 && ncol(Ymat) <= 1) {
          ok <- FALSE
          msg <- c(msg, "jSDM latent-variable models require at least two response/species columns.")
        }
        if (identical(model_type, "binomial_probit_sp_constrained") && (n_latent < 1 || ncol(Ymat) < (n_latent + 2))) {
          ok <- FALSE
          msg <- c(msg, "jSDM constrained probit needs n_latent > 0 and more species than latent axes; use at least n_latent + 2 species.")
        }
        if (model_type %in% c("binomial_probit", "binomial_probit_sp_constrained")) {
          fam <- response_family_messages(Ymat, "binomial_probit", "jSDM")
          ok <- ok && isTRUE(fam$ok)
          msg <- c(msg, fam$messages)
        } else if (identical(model_type, "binomial_logit")) {
          if (any(Ymat < 0, na.rm = TRUE) || any(abs(Ymat - round(Ymat)) > 1e-8, na.rm = TRUE)) {
            ok <- FALSE
            msg <- c(msg, "jSDM binomial_logit responses must be non-negative integer successes.")
          }
          trial_vec <- NULL
          if (!is.null(trials)) {
            trial_raw <- suppressWarnings(as.numeric(as.matrix(trials)))
            if (length(trial_raw) == 1) {
              trial_vec <- rep(trial_raw, nrow(Ymat))
            } else if (length(trial_raw) == nrow(Ymat)) {
              trial_vec <- trial_raw
            } else if (length(trial_raw) == length(Ymat)) {
              trial_matrix <- matrix(trial_raw, nrow = nrow(Ymat), ncol = ncol(Ymat))
              if (all(apply(trial_matrix, 1, function(z) length(unique(z[is.finite(z)])) <= 1))) {
                trial_vec <- trial_matrix[, 1]
              } else {
                ok <- FALSE
                msg <- c(msg, "jSDM binomial_logit supports one trial count per site, applied to all species. trials.csv has species-specific counts.")
              }
            }
          } else {
            trial_vec <- rep(as.numeric(trials_scalar %||% 1), nrow(Ymat))
          }
          if (is.null(trial_vec) || any(!is.finite(trial_vec)) || any(trial_vec < 1)) {
            ok <- FALSE
            msg <- c(msg, "jSDM binomial_logit trials must provide one positive trial count or one positive count per site.")
          } else if (any(Ymat > matrix(trial_vec, nrow = nrow(Ymat), ncol = ncol(Ymat)), na.rm = TRUE)) {
            ok <- FALSE
            msg <- c(msg, "jSDM binomial_logit successes in Y.csv must be <= trials for each site.")
          }
        } else if (identical(model_type, "poisson_log")) {
          fam <- response_family_messages(Ymat, "poisson_log", "jSDM")
          ok <- ok && isTRUE(fam$ok)
          msg <- c(msg, fam$messages)
        } else if (identical(model_type, "gaussian")) {
          fam <- response_family_messages(Ymat, "gaussian", "jSDM")
          ok <- ok && isTRUE(fam$ok)
          msg <- c(msg, fam$messages)
        }
      }
      x_clean <- clean_predictor_types(X, "jSDM XData")
      factor_cols <- names(x_clean)[vapply(x_clean, is.factor, logical(1))]
      if (length(factor_cols) > 0) msg <- c(msg, paste0("jSDM categorical predictors will be factor/model-matrix encoded: ", paste(factor_cols, collapse = ", ")))
      fcheck <- validate_one_sided_formula(site_formula, x_clean, "jSDM site_formula")
      ok <- ok && isTRUE(fcheck$ok)
      msg <- c(msg, fcheck$messages)
      if (!is.null(Tr) && isTRUE(allow_traits)) {
        if (!is.null(Y) && nrow(Tr) != ncol(Y)) {
          ok <- FALSE
          msg <- c(msg, "jSDM trait_data.csv rows must match Y.csv columns/species.")
        }
        tr_clean <- clean_predictor_types(Tr, "jSDM trait_data")
        trcheck <- validate_one_sided_formula(trait_formula, tr_clean, "jSDM trait_formula")
        ok <- ok && isTRUE(trcheck$ok)
        msg <- c(msg, trcheck$messages)
      }
      if (!is.null(newdata)) {
        nd_clean <- clean_predictor_types(newdata, "jSDM newdata")
        nfcheck <- validate_one_sided_formula(site_formula, nd_clean, "jSDM prediction newdata")
        ok <- ok && isTRUE(nfcheck$ok)
        msg <- c(msg, nfcheck$messages)
      }
    }
  }
  if (isTRUE(do_predict) && n_latent > 0 && identical(site_effect, "random")) {
    msg <- c(msg, "Prediction with latent variables/random site effects uses existing site IDs; newdata should keep row names matching training sites for safest predict.jSDM use.")
  }
  if (length(msg) == 0) msg <- "jSDM data and settings check passed."
  list(ok = ok, messages = msg)
}

write_jsdm_reproducible_script <- function(outdir) {
  dir.create(file.path(outdir, "reproducible_script"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "workflow_scripts"), recursive = TRUE, showWarnings = FALSE)
  script <- c(
    "# Reproducible jSDM 0.2.7 analysis script generated by JSDM Studio",
    "# Run from the output folder with: Rscript reproducible_script/run_this_jSDM_analysis.R",
    "args <- commandArgs(trailingOnly = FALSE)",
    "file_arg <- '--file='",
    "script_arg <- args[startsWith(args, file_arg)]",
    "if (length(script_arg) > 0) {",
    "  script_path <- normalizePath(sub(file_arg, '', script_arg[[1]]), winslash = '/', mustWork = FALSE)",
    "  setwd(dirname(dirname(script_path)))",
    "}",
    "for (d in c('models','mcmc','tables','predictions','plots','diagnostics','standard','report','results')) dir.create(d, showWarnings = FALSE, recursive = TRUE)",
    "writeLines(c('jSDM plot outputs', '=================', 'Diagnostic and result figures are written here when enabled. If figure export is disabled, this README keeps the folder intentionally non-empty.'), file.path('plots','README_jSDM_plots.txt'))",
    "`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x",
    "as_bool <- function(x, default = FALSE) { if (is.null(x) || length(x) == 0 || is.na(x)) return(default); if (is.logical(x)) return(isTRUE(x)); tolower(as.character(x)) %in% c('true','t','1','yes','y') }",
    "as_num <- function(x, default) { z <- suppressWarnings(as.numeric(x)); if (!length(z) || !is.finite(z[1])) default else z[1] }",
    "as_int <- function(x, default) as.integer(round(as_num(x, default)))",
    "append_warning <- function(msg) step_warnings <<- unique(c(step_warnings, as.character(msg)))",
    "write_status <- function(status, errors = character()) {",
    "  payload <- list(engine='jSDM', status=status, warnings=as.character(step_warnings), errors=as.character(errors), time=as.character(Sys.time()))",
    "  if (requireNamespace('jsonlite', quietly = TRUE)) writeLines(jsonlite::toJSON(payload, pretty=TRUE, auto_unbox=TRUE), file.path('diagnostics','engine_status.json'))",
    "  write.csv(data.frame(engine='jSDM', status=status, message=paste(c(step_warnings, errors), collapse='; ')), file.path('tables','engine_status.csv'), row.names=FALSE)",
    "  try(writeLines(capture.output(sessionInfo()), file.path('diagnostics','session_info.txt')), silent=TRUE)",
    "}",
    "safe_step <- function(name, expr, required = FALSE) {",
    "  tryCatch(withCallingHandlers(force(expr), warning = function(w) { append_warning(paste(name, conditionMessage(w), sep=': ')); invokeRestart('muffleWarning') }), error = function(e) { msg <- conditionMessage(e); append_warning(paste(name, msg, sep=': ')); writeLines(msg, file.path('diagnostics', paste0('jSDM_', gsub('[^A-Za-z0-9_]+', '_', name), '_error.txt'))); if (required) stop(msg, call. = FALSE); NULL })",
    "}",
    "read_csv_safe <- function(path) {",
    "  if (!file.exists(path)) return(NULL)",
    "  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)",
    "  if (ncol(dat) > 1) {",
    "    first <- dat[[1]]",
    "    if (!anyDuplicated(first) && !all(suppressWarnings(!is.na(as.numeric(first))))) { dat <- dat[-1]; rownames(dat) <- first }",
    "  }",
    "  dat",
    "}",
    "factor_preserve_order <- function(x) { x <- as.character(x); x[trimws(x) == ''] <- NA; factor(x, levels = unique(x[!is.na(x)])) }",
    "clean_df <- function(df) {",
    "  if (is.null(df)) return(NULL)",
    "  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)",
    "  for (nm in names(df)) {",
    "    if (is.character(df[[nm]])) { x <- trimws(df[[nm]]); x[x == ''] <- NA; nx <- suppressWarnings(as.numeric(x)); df[[nm]] <- if (all(is.na(x) | !is.na(nx))) nx else factor_preserve_order(x) }",
    "    if (is.logical(df[[nm]])) df[[nm]] <- factor_preserve_order(df[[nm]])",
    "  }",
    "  df",
    "}",
    "scale_site_data <- function(df) {",
    "  if (is.null(df)) return(NULL)",
    "  for (nm in names(df)) if (is.numeric(df[[nm]]) && length(unique(df[[nm]][is.finite(df[[nm]])])) > 2 && stats::sd(df[[nm]], na.rm=TRUE) > 0) df[[nm]] <- as.numeric(scale(df[[nm]]))",
    "  df",
    "}",
    "as_response_matrix <- function(dat, label='Y') {",
    "  if (is.null(dat)) stop(label, ' is missing.', call. = FALSE)",
    "  mat <- as.matrix(dat); storage.mode(mat) <- 'numeric'",
    "  if (any(!is.finite(mat))) stop(label, ' must be numeric and finite for jSDM fitting.', call. = FALSE)",
    "  if (is.null(rownames(mat))) rownames(mat) <- paste0('site_', seq_len(nrow(mat)))",
    "  if (is.null(colnames(mat))) colnames(mat) <- paste0('sp_', seq_len(ncol(mat)))",
    "  mat",
    "}",
    "parse_probs <- function(x) { z <- suppressWarnings(as.numeric(strsplit(as.character(x %||% ''), ',')[[1]])); z <- z[is.finite(z) & z >= 0 & z <= 1]; if (length(z) == 0) c(0.025, 0.975) else z }",
    "parse_ids <- function(x, available, max_n = length(available)) {",
    "  if (is.null(available) || length(available) == 0) return(character())",
    "  txt <- trimws(as.character(x %||% 'all'))",
    "  if (!nzchar(txt) || identical(tolower(txt), 'all')) return(available[seq_len(min(length(available), max_n))])",
    "  parts <- trimws(unlist(strsplit(txt, ','))); parts <- parts[nzchar(parts)]",
    "  idx <- suppressWarnings(as.integer(parts)); out <- character()",
    "  if (all(!is.na(idx))) out <- available[idx[idx >= 1 & idx <= length(available)]] else out <- intersect(parts, available)",
    "  if (length(out) == 0) out <- available[seq_len(min(length(available), max_n))]",
    "  out[seq_len(min(length(out), max_n))]",
    "}",
    "prediction_id_value <- function(dat, choices, fallback) {",
    "  if (is.null(dat)) return(fallback)",
    "  hit <- intersect(choices, names(dat))",
    "  if (!length(hit)) return(fallback)",
    "  vals <- unique(trimws(as.character(dat[[hit[1]]]))); vals <- vals[nzchar(vals) & !is.na(vals)]",
    "  if (!length(vals)) fallback else paste(vals, collapse=',')",
    "}",
    "make_trials_vector <- function(trials_raw, Y, scalar = 1) {",
    "  vals <- suppressWarnings(as.numeric(as.matrix(trials_raw)))",
    "  if (length(vals) == 0 || all(is.na(vals))) vals <- scalar",
    "  if (length(vals) == 1) out <- rep(vals, nrow(Y)) else if (length(vals) == nrow(Y)) out <- vals else if (length(vals) == length(Y)) { mat <- matrix(vals, nrow(Y), ncol(Y)); if (!all(apply(mat, 1, function(z) length(unique(z[is.finite(z)])) <= 1))) stop('jSDM_binomial_logit supports one trial count per site, not species-specific trial counts.', call. = FALSE); out <- mat[,1] } else stop('trials.csv must contain one value or one value per site for jSDM_binomial_logit.', call. = FALSE)",
    "  if (any(!is.finite(out)) || any(out < 1)) stop('jSDM binomial_logit trials must be positive finite values.', call. = FALSE)",
    "  out",
    "}",
    "mcmc_summary <- function(obj, component = 'mcmc', group = NA_character_) {",
    "  rows <- list()",
    "  add_mat <- function(x, comp, grp) {",
    "    m <- as.matrix(x); if (is.null(colnames(m))) colnames(m) <- paste0('par_', seq_len(ncol(m)))",
    "    data.frame(component=comp, group=grp, parameter=colnames(m), mean=colMeans(m, na.rm=TRUE), sd=apply(m,2,sd,na.rm=TRUE), median=apply(m,2,median,na.rm=TRUE), lower=apply(m,2,quantile,0.025,na.rm=TRUE), upper=apply(m,2,quantile,0.975,na.rm=TRUE), stringsAsFactors=FALSE)",
    "  }",
    "  if (is.null(obj)) return(data.frame())",
    "  if (is.list(obj) && !inherits(obj, 'mcmc')) { for (nm in names(obj)) rows[[length(rows)+1]] <- mcmc_summary(obj[[nm]], component, nm) } else rows[[length(rows)+1]] <- add_mat(obj, component, group)",
    "  do.call(rbind, rows)",
    "}",
    "write_matrix_if <- function(x, path) { if (!is.null(x)) write.csv(as.data.frame(x), path) }",
    "prediction_slot <- function(model) {",
    "  for (nm in c('theta_latent','Y_pred','probit_theta_latent','logit_theta_latent','log_theta_latent')) if (!is.null(model[[nm]]) && is.matrix(model[[nm]])) return(model[[nm]])",
    "  NULL",
    "}",
    "standard_predictions <- function(pred, observed = NULL) {",
    "  if (is.null(pred) || !is.matrix(pred)) return(data.frame(engine='jSDM', site_id=NA_character_, response_id=NA_character_, observed=NA_real_, predicted_mean=NA_real_, predicted_lower=NA_real_, predicted_upper=NA_real_))",
    "  if (is.null(rownames(pred))) rownames(pred) <- paste0('site_', seq_len(nrow(pred))); if (is.null(colnames(pred))) colnames(pred) <- paste0('sp_', seq_len(ncol(pred)))",
    "  grid <- expand.grid(site_id=rownames(pred), response_id=colnames(pred), KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)",
    "  obs <- if (!is.null(observed) && all(dim(observed) == dim(pred))) as.vector(observed) else NA_real_",
    "  data.frame(engine='jSDM', site_id=grid$site_id, response_id=grid$response_id, observed=obs, predicted_mean=as.vector(pred), predicted_lower=NA_real_, predicted_upper=NA_real_, stringsAsFactors=FALSE)",
    "}",
    "fit_metrics <- function(obs, pred, family) {",
    "  if (is.null(obs) || is.null(pred) || !all(dim(obs) == dim(pred))) return(data.frame(engine='jSDM', metric='prediction_available', response_id=NA_character_, value=as.numeric(!is.null(pred)), notes='dimension mismatch or no observed matrix'))",
    "  rows <- list(); fam <- as.character(family)",
    "  for (j in seq_len(ncol(obs))) { o <- obs[,j]; p <- pred[,j]; nm <- colnames(obs)[j] %||% paste0('sp_',j); ok <- is.finite(o) & is.finite(p); if (!any(ok)) next; o <- o[ok]; p <- p[ok]",
    "    rows[[length(rows)+1]] <- data.frame(engine='jSDM', metric='RMSE', response_id=nm, value=sqrt(mean((o-p)^2)), notes=fam)",
    "    if (stats::var(o) > 0) rows[[length(rows)+1]] <- data.frame(engine='jSDM', metric='R2_like', response_id=nm, value=1 - sum((o-p)^2)/sum((o-mean(o))^2), notes=fam)",
    "    if (all(o %in% c(0,1)) && any(o==1) && any(o==0)) rows[[length(rows)+1]] <- data.frame(engine='jSDM', metric='TjurR2', response_id=nm, value=mean(p[o==1])-mean(p[o==0]), notes=fam)",
    "  }",
    "  if (length(rows) == 0) data.frame(engine='jSDM', metric='prediction_available', response_id=NA_character_, value=as.numeric(!is.null(pred)), notes=fam) else do.call(rbind, rows)",
    "}",
    "assoc_long <- function(mat, type) {",
    "  if (is.null(mat) || !is.matrix(mat) || nrow(mat) < 2) return(data.frame(engine='jSDM', response_1=NA_character_, response_2=NA_character_, association_type=type, estimate=NA_real_, comparable_level='not_available'))",
    "  ids <- which(upper.tri(mat), arr.ind=TRUE)",
    "  data.frame(engine='jSDM', response_1=rownames(mat)[ids[,1]], response_2=colnames(mat)[ids[,2]], association_type=type, estimate=mat[ids], comparable_level='fitted', stringsAsFactors=FALSE)",
    "}",
    "try({ if (!requireNamespace('yaml', quietly=TRUE)) stop('yaml is required.'); cfg <- yaml::read_yaml('used_config.yml') }, silent = FALSE)",
    "step_warnings <- character()",
    "tryCatch({",
    "  if (!requireNamespace('jSDM', quietly=TRUE)) stop('Package jSDM is required for real fitting.', call. = FALSE)",
    "  suppressPackageStartupMessages(library(jSDM))",
    "  suppressPackageStartupMessages(library(coda))",
    "  model_type <- as.character(cfg$model$model_type %||% 'binomial_probit')",
    "  if (isTRUE(as_bool(cfg$model$constrained_latent, FALSE)) && identical(model_type, 'binomial_probit')) model_type <- 'binomial_probit_sp_constrained'",
    "  burnin <- as_int(cfg$mcmc$burnin, 50); mcmc <- as_int(cfg$mcmc$mcmc, 50); thin <- as_int(cfg$mcmc$thin, 1)",
    "  if ((burnin + mcmc) < 100 || ((burnin + mcmc) %% 10) != 0 || (mcmc %% thin) != 0) stop('Invalid jSDM MCMC settings: burnin+mcmc must be >=100 and divisible by 10; mcmc must be divisible by thin.', call. = FALSE)",
    "  Y_raw <- read_csv_safe(file.path('data','Y.csv')); X_raw <- read_csv_safe(file.path('data','XData.csv')); Tr_raw <- read_csv_safe(file.path('data','trait_data.csv')); long_raw <- read_csv_safe(file.path('data','long_format.csv')); trials_raw <- read_csv_safe(file.path('data','trials.csv')); new_raw <- read_csv_safe(file.path('data','newdata.csv')); prediction_ids_raw <- read_csv_safe(file.path('data','prediction_ids.csv'))",
    "  Y <- if (!identical(model_type, 'binomial_probit_long_format')) as_response_matrix(Y_raw, 'Y.csv') else NULL",
    "  site_data <- clean_df(X_raw); if (!is.null(site_data) && nrow(site_data) != nrow(Y %||% site_data)) stop('XData.csv rows must match Y.csv rows.', call. = FALSE)",
    "  if (as_bool(cfg$model$scale_site_data, TRUE)) site_data <- scale_site_data(site_data)",
    "  if (!is.null(Y) && !is.null(site_data)) rownames(site_data) <- rownames(Y)",
    "  write.csv(site_data %||% data.frame(), file.path('tables','site_data_used.csv'))",
    "  trait_data <- if (as_bool(cfg$model$allow_traits, TRUE)) clean_df(Tr_raw) else NULL",
    "  if (!is.null(trait_data) && !is.null(Y)) { if (nrow(trait_data) != ncol(Y)) stop('trait_data.csv rows must match Y species columns.', call. = FALSE); rownames(trait_data) <- colnames(Y) }",
    "  site_formula <- as.formula(cfg$model$site_formula %||% '~ .')",
    "  if (!as_bool(cfg$model$include_intercept, TRUE)) { rhs <- trimws(as.character(site_formula)[2]); if (!grepl('(^|[+ ])0([+ ]|$)|- *1', rhs)) site_formula <- as.formula(paste('~ 0 +', rhs)) }",
    "  trait_formula <- if (!is.null(trait_data)) as.formula(cfg$model$trait_formula %||% '~ .') else NULL",
    "  n_latent <- as_int(cfg$model$n_latent, 0); site_effect <- as.character(cfg$model$site_effect %||% 'none'); seed <- as_int(cfg$mcmc$seed, 1234)",
    "  common <- list(burnin=burnin, mcmc=mcmc, thin=thin, site_formula=site_formula, n_latent=n_latent, site_effect=site_effect, beta_start=as_num(cfg$starts$beta_start,0), gamma_start=as_num(cfg$starts$gamma_start,0), lambda_start=as_num(cfg$starts$lambda_start,0), W_start=as_num(cfg$starts$W_start,0), alpha_start=as_num(cfg$starts$alpha_start,0), V_alpha=as_num(cfg$starts$V_alpha,1), shape_Valpha=as_num(cfg$priors$shape_Valpha,0.5), rate_Valpha=as_num(cfg$priors$rate_Valpha,0.0005), mu_beta=as_num(cfg$priors$mu_beta,0), V_beta=as_num(cfg$priors$V_beta,10), mu_gamma=as_num(cfg$priors$mu_gamma,0), V_gamma=as_num(cfg$priors$V_gamma,10), mu_lambda=as_num(cfg$priors$mu_lambda,0), V_lambda=as_num(cfg$priors$V_lambda,10), seed=seed, verbose=as_int(cfg$mcmc$verbose,0))",
    "  if (!is.null(trait_data) && !identical(model_type, 'binomial_probit_long_format')) { common$trait_data <- trait_data; common$trait_formula <- trait_formula }",
    "  fit_call <- switch(model_type,",
    "    binomial_probit = { common$presence_data <- Y; common$site_data <- site_data; quote(do.call(jSDM::jSDM_binomial_probit, common)) },",
    "    binomial_logit = { common$presence_data <- Y; common$site_data <- site_data; common$trials <- make_trials_vector(trials_raw, Y, as_num(cfg$model$trials,1)); common$ropt <- as_num(cfg$mcmc$ropt,0.44); quote(do.call(jSDM::jSDM_binomial_logit, common)) },",
    "    poisson_log = { common$count_data <- Y; common$site_data <- site_data; common$ropt <- as_num(cfg$mcmc$ropt,0.44); quote(do.call(jSDM::jSDM_poisson_log, common)) },",
    "    gaussian = { common$response_data <- Y; common$site_data <- site_data; common$V_start <- as_num(cfg$starts$V_start,1); common$shape_V <- as_num(cfg$priors$shape_V,0.5); common$rate_V <- as_num(cfg$priors$rate_V,0.0005); quote(do.call(jSDM::jSDM_gaussian, common)) },",
    "    binomial_probit_sp_constrained = { common$presence_data <- Y; common$site_data <- site_data; common$nchains <- as_int(cfg$model$constrained_nchains,2); common$ncores <- 1L; common$seed <- as.integer(seed + seq_len(common$nchains) - 1L); quote(do.call(jSDM::jSDM_binomial_probit_sp_constrained, common)) },",
    "    binomial_probit_long_format = { long_dat <- clean_df(long_raw); names(long_dat)[names(long_dat) == (cfg$model$long_site_col %||% 'site')] <- 'site'; names(long_dat)[names(long_dat) == (cfg$model$long_species_col %||% 'species')] <- 'species'; names(long_dat)[names(long_dat) == (cfg$model$long_response_col %||% 'presence')] <- 'Y'; long_covars <- setdiff(names(long_dat), c('site','species','Y')); sf_txt <- paste(deparse(common$site_formula), collapse=' '); if (!grepl('species', sf_txt, fixed=TRUE)) { rhs <- if (length(long_covars)) paste(c('species', paste0('species:', long_covars)), collapse=' + ') else 'species'; common$site_formula <- as.formula(paste('~', rhs)); append_warning(paste('long-format site_formula was rewritten for jSDM 0.2.7:', deparse(common$site_formula))) }; common$data <- long_dat; common$trait_data <- NULL; common$trait_formula <- NULL; quote(do.call(jSDM::jSDM_binomial_probit_long_format, common)) },",
    "    stop('Unsupported jSDM model_type: ', model_type, call. = FALSE)",
    "  )",
    "  model <- eval(fit_call)",
    "  model_primary <- model",
    "  if (!inherits(model_primary, 'jSDM') && is.list(model_primary)) { idx <- which(vapply(model_primary, function(z) inherits(z, 'jSDM'), logical(1))); if (length(idx) > 0) model_primary <- model_primary[[idx[1]]] }",
    "  if (!inherits(model_primary, 'jSDM')) stop('jSDM fitting returned no usable jSDM object.', call. = FALSE)",
    "  if (as_bool(cfg$outputs$save_model, TRUE)) saveRDS(model, file.path('models','jsdm_model.rds'))",
    "  saveRDS(model_primary, file.path('models','jsdm_model_primary.rds'))",
    "  saveRDS(model_primary$model_spec, file.path('models','jsdm_model_spec.rds'))",
    "  write.csv(data.frame(field=names(model_primary$model_spec), value=vapply(model_primary$model_spec, function(x) paste(utils::capture.output(str(x, max.level=1)), collapse=' '), character(1))), file.path('tables','model_spec.csv'), row.names=FALSE)",
    "  if (as_bool(cfg$outputs$save_mcmc_rds, TRUE)) { saveRDS(model_primary$mcmc.sp, file.path('mcmc','mcmc_sp.rds')); saveRDS(model_primary$mcmc.gamma %||% NULL, file.path('mcmc','mcmc_gamma.rds')); saveRDS(model_primary$mcmc.latent %||% NULL, file.path('mcmc','mcmc_latent.rds')); saveRDS(model_primary$mcmc.alpha %||% NULL, file.path('mcmc','mcmc_alpha.rds')); saveRDS(model_primary$mcmc.V_alpha %||% NULL, file.path('mcmc','mcmc_V_alpha.rds')); saveRDS(model_primary$mcmc.V %||% NULL, file.path('mcmc','mcmc_V.rds')); saveRDS(model_primary$mcmc.Deviance %||% NULL, file.path('mcmc','mcmc_Deviance.rds')) }",
    "  sp_sum <- mcmc_summary(model_primary$mcmc.sp, 'mcmc.sp'); if (!nrow(sp_sum)) sp_sum <- data.frame(component=character(), group=character(), parameter=character(), mean=numeric(), sd=numeric(), median=numeric(), lower=numeric(), upper=numeric(), stringsAsFactors=FALSE); write.csv(sp_sum, file.path('tables','species_parameter_summary.csv'), row.names=FALSE); write.csv(sp_sum[grepl('^beta_', sp_sum$parameter),,drop=FALSE], file.path('tables','beta_summary.csv'), row.names=FALSE); write.csv(sp_sum[grepl('^lambda_', sp_sum$parameter),,drop=FALSE], file.path('tables','lambda_summary.csv'), row.names=FALSE)",
    "  gamma_sum <- mcmc_summary(model_primary$mcmc.gamma %||% NULL, 'mcmc.gamma'); write.csv(gamma_sum, file.path('tables','gamma_summary.csv'), row.names=FALSE)",
    "  alpha_sum <- mcmc_summary(model_primary$mcmc.alpha %||% NULL, 'mcmc.alpha'); write.csv(alpha_sum, file.path('tables','alpha_summary.csv'), row.names=FALSE)",
    "  write.csv(mcmc_summary(model_primary$mcmc.V_alpha %||% NULL, 'mcmc.V_alpha'), file.path('tables','V_alpha_summary.csv'), row.names=FALSE)",
    "  write.csv(mcmc_summary(model_primary$mcmc.V %||% NULL, 'mcmc.V'), file.path('tables','V_summary.csv'), row.names=FALSE)",
    "  write.csv(mcmc_summary(model_primary$mcmc.Deviance %||% NULL, 'mcmc.Deviance'), file.path('tables','Deviance_summary.csv'), row.names=FALSE)",
    "  pred_mat <- prediction_slot(model_primary); write_matrix_if(pred_mat, file.path('predictions','fitted_values.csv')); saveRDS(pred_mat, file.path('predictions','fitted_values.rds'))",
    "  if (as_bool(cfg$prediction$do_predict, TRUE)) safe_step('predict.jSDM', { sites <- rownames(pred_mat %||% model_primary$model_spec$presence_data %||% model_primary$model_spec$count_data %||% model_primary$model_spec$response_data); species <- colnames(pred_mat %||% model_primary$model_spec$presence_data %||% model_primary$model_spec$count_data %||% model_primary$model_spec$response_data); max_sites <- as_int(cfg$prediction$max_prediction_sites, 500); site_req <- prediction_id_value(prediction_ids_raw, c('Id_sites','site','site_id','sites'), cfg$prediction$Id_sites); sp_req <- prediction_id_value(prediction_ids_raw, c('Id_species','species','species_id','responses'), cfg$prediction$Id_species); Id_sites <- parse_ids(site_req, sites, max_sites); Id_species <- parse_ids(sp_req, species, length(species)); nd <- clean_df(new_raw); if (!is.null(nd)) { if (as_bool(cfg$model$scale_site_data, TRUE)) nd <- scale_site_data(nd); if (!is.null(rownames(nd)) && any(rownames(nd) %in% sites)) Id_sites <- intersect(rownames(nd), sites); nd <- nd[seq_len(min(nrow(nd), length(Id_sites))),,drop=FALSE] }; pobj <- predict(model_primary, newdata=nd, Id_species=Id_species, Id_sites=Id_sites, type=cfg$prediction$predict_type %||% 'mean', probs=parse_probs(cfg$prediction$predict_probs)); saveRDS(pobj, file.path('predictions','predict_jSDM.rds')); if (is.matrix(pobj) || is.data.frame(pobj)) write.csv(pobj, file.path('predictions','predict_jSDM.csv')) else if (is.list(pobj)) write.csv(data.frame(object=names(pobj), rows=vapply(pobj, NROW, integer(1)), cols=vapply(pobj, NCOL, integer(1))), file.path('predictions','predict_jSDM_manifest.csv'), row.names=FALSE) })",
    "  obs <- Y; if (is.null(obs) && !is.null(model_primary$model_spec$presence_data)) obs <- model_primary$model_spec$presence_data; if (is.null(obs) && !is.null(model_primary$model_spec$count_data)) obs <- model_primary$model_spec$count_data; if (is.null(obs) && !is.null(model_primary$model_spec$response_data)) obs <- model_primary$model_spec$response_data",
    "  write.csv(standard_predictions(pred_mat, obs), file.path('standard','predictions_long.csv'), row.names=FALSE)",
    "  write.csv(fit_metrics(obs, pred_mat, paste(model_primary$model_spec$family, model_primary$model_spec$link, sep='_')), file.path('standard','fit_metrics.csv'), row.names=FALSE)",
    "  rc <- NULL; if (as_bool(cfg$outputs$residual_cor, TRUE) && as_int(model_primary$model_spec$n_latent,0) > 1) rc <- safe_step('get_residual_cor', jSDM::get_residual_cor(model_primary, prob=as_num(cfg$diagnostics$cor_prob,0.95), type=cfg$diagnostics$cor_type %||% 'mean')) else append_warning('residual correlations skipped: jSDM requires n_latent > 1')",
    "  if (!is.null(rc)) { for (nm in names(rc)) write_matrix_if(rc[[nm]], file.path('tables', paste0('residual_cor_', nm, '.csv'))); write.csv(assoc_long(rc$cor.mean, 'residual_cor'), file.path('standard','associations_long.csv'), row.names=FALSE) } else write.csv(assoc_long(NULL, 'residual_cor'), file.path('standard','associations_long.csv'), row.names=FALSE)",
    "  ec <- NULL; if (as_bool(cfg$outputs$enviro_cor, TRUE)) ec <- safe_step('get_enviro_cor', jSDM::get_enviro_cor(model_primary, prob=as_num(cfg$diagnostics$cor_prob,0.95), type=cfg$diagnostics$cor_type %||% 'mean'))",
    "  if (!is.null(ec)) for (nm in names(ec)) write_matrix_if(ec[[nm]], file.path('tables', paste0('enviro_cor_', nm, '.csv')))",
    "  beta_effects <- if (exists('sp_sum') && nrow(sp_sum)) sp_sum[grepl('^beta_', sp_sum$parameter),,drop=FALSE] else data.frame(); if (nrow(beta_effects)) { effects <- data.frame(engine='jSDM', response_id=beta_effects$group, predictor=sub('^beta_', '', beta_effects$parameter), direction=ifelse(beta_effects$mean>0,'positive',ifelse(beta_effects$mean<0,'negative','zero')), estimate=beta_effects$mean, lower=beta_effects$lower, upper=beta_effects$upper, notes='mcmc.sp beta summary', stringsAsFactors=FALSE) } else effects <- data.frame(engine='jSDM', response_id=NA_character_, predictor=NA_character_, direction=NA_character_, estimate=NA_real_, lower=NA_real_, upper=NA_real_, notes='no beta summary', stringsAsFactors=FALSE); write.csv(effects, file.path('standard','effects_long.csv'), row.names=FALSE)",
    "  if (as_bool(cfg$outputs$figures, TRUE)) safe_step('plots', { pdf(file.path('plots','jSDM_diagnostics.pdf'), width=8, height=6); on.exit(dev.off(), add=TRUE); if (nrow(sp_sum)) { op <- par(mar=c(8,5,3,1)); barplot(head(sp_sum$mean, 20), names.arg=head(paste(sp_sum$group, sp_sum$parameter, sep=':'),20), las=2, main='jSDM posterior means'); par(op) }; if (!is.null(pred_mat)) hist(as.vector(pred_mat), main='jSDM fitted predictions', xlab='prediction'); if (!is.null(rc) && !is.null(rc$cor.mean)) image(rc$cor.mean, main='Residual correlation') })",
    "  if (as_bool(cfg$diagnostics$plot_residual_cor, TRUE) && !is.null(rc)) safe_step('plot_residual_cor', { pdf(file.path('plots','residual_cor.pdf'), width=7, height=7); on.exit(dev.off(), add=TRUE); jSDM::plot_residual_cor(model_primary, prob=as_num(cfg$diagnostics$cor_prob,0.95)) })",
    "  write.csv(data.frame(run_id=basename(getwd()), engine='jSDM', status='fitted', model_type=model_type, n_sites=if(!is.null(obs)) nrow(obs) else NA_integer_, n_responses=if(!is.null(obs)) ncol(obs) else length(unique(model_primary$model_spec$data$species)), n_predictors=if(!is.null(site_data)) ncol(site_data) else NA_integer_, burnin=burnin, mcmc=mcmc, thin=thin, n_latent=n_latent, site_effect=site_effect, stringsAsFactors=FALSE), file.path('standard','run_summary.csv'), row.names=FALSE)",
    "  file.copy(file.path('standard','run_summary.csv'), file.path('results','jSDM_run_summary.csv'), overwrite=TRUE)",
    "  file.copy(file.path('standard','fit_metrics.csv'), file.path('results','fit_metrics.csv'), overwrite=TRUE)",
    "  file.copy(file.path('standard','effects_long.csv'), file.path('results','effects_long.csv'), overwrite=TRUE)",
    "  file.copy(file.path('standard','predictions_long.csv'), file.path('results','predictions_long.csv'), overwrite=TRUE)",
    "  writeLines(c('jSDM results', '============', paste0('Status: fitted'), paste0('Model type: ', model_type), '', 'Important files:', '- standard/run_summary.csv and results/jSDM_run_summary.csv', '- tables/*_summary.csv for posterior summaries', '- predictions/fitted_values.csv and standard/predictions_long.csv', '- diagnostics/engine_status.json and diagnostics/session_info.txt', '- reproducible_script/run_this_jSDM_analysis.R'), file.path('results','README_jSDM_results.txt'))",
    "  write.csv(data.frame(file=list.files('.', recursive=TRUE), stringsAsFactors=FALSE), file.path('standard','output_manifest.csv'), row.names=FALSE)",
    "  writeLines(c('<!DOCTYPE html><html><head><meta charset=\"UTF-8\"><title>jSDM Report</title></head><body><h1>jSDM workflow report</h1>', paste0('<p>Status: fitted</p><p>Model type: ', model_type, '</p>'), '<p>Open diagnostics/engine_status.json and standard/*.csv for machine-readable outputs.</p></body></html>'), file.path('report','jSDM_report.html'))",
    "  writeLines('RUN COMPLETE', 'RUN_COMPLETE.txt')",
    "  write_status('fitted')",
    "}, error = function(e) {",
    "  msg <- conditionMessage(e)",
    "  writeLines(msg, file.path('diagnostics','jSDM_reproducible_error.txt'))",
    "  write.csv(data.frame(engine='jSDM', response_id=NA_character_, predictor=NA_character_, direction=NA_character_, estimate=NA_real_, lower=NA_real_, upper=NA_real_, notes='fit_failed'), file.path('standard','effects_long.csv'), row.names=FALSE)",
    "  write.csv(data.frame(engine='jSDM', site_id=NA_character_, response_id=NA_character_, observed=NA_real_, predicted_mean=NA_real_, predicted_lower=NA_real_, predicted_upper=NA_real_), file.path('standard','predictions_long.csv'), row.names=FALSE)",
    "  write.csv(data.frame(engine='jSDM', response_1=NA_character_, response_2=NA_character_, association_type=NA_character_, estimate=NA_real_, comparable_level='fit_failed'), file.path('standard','associations_long.csv'), row.names=FALSE)",
    "  write.csv(data.frame(engine='jSDM', metric=NA_character_, response_id=NA_character_, value=NA_real_, notes='fit_failed'), file.path('standard','fit_metrics.csv'), row.names=FALSE)",
    "  write.csv(data.frame(run_id=basename(getwd()), engine='jSDM', status='fit_failed'), file.path('standard','run_summary.csv'), row.names=FALSE)",
    "  writeLines(c('jSDM results', '============', 'Status: fit_failed', '', 'Open diagnostics/jSDM_reproducible_error.txt and diagnostics/engine_status.json first.'), file.path('results','README_jSDM_results.txt'))",
    "  write_status('fit_failed', msg)",
    "  writeLines('RUN FAILED', 'RUN_FAILED.txt')",
    "  quit(status=1, save='no')",
    "})"
  )
  writeLines(script, file.path(outdir, "reproducible_script", "run_this_jSDM_analysis.R"))
  writeLines(c(
    "# Workflow wrapper for jSDM",
    "source(file.path('reproducible_script', 'run_this_jSDM_analysis.R'))"
  ), file.path(outdir, "workflow_scripts", "run_jSDM_workflow.R"))
}
