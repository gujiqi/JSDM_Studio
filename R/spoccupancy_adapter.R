# Real spOccupancy adapter for JSDM Studio.

spocc_spatial_models <- c("spPGOcc", "spMsPGOcc", "sfJSDM", "sfMsPGOcc", "spIntPGOcc",
                          "stPGOcc", "stMsPGOcc", "svcPGBinom", "svcPGOcc",
                          "svcTPGBinom", "svcTPGOcc", "svcMsPGOcc", "svcTMsPGOcc")
spocc_multi_models <- c("msPGOcc", "spMsPGOcc", "lfMsPGOcc", "sfMsPGOcc",
                        "tMsPGOcc", "stMsPGOcc", "svcMsPGOcc", "svcTMsPGOcc")
spocc_integrated_models <- c("intPGOcc", "spIntPGOcc", "stIntPGOcc", "intMsPGOcc")
spocc_lf_models <- c("lfMsPGOcc", "sfMsPGOcc")
spocc_svc_models <- c("svcPGBinom", "svcPGOcc", "svcTPGBinom", "svcTPGOcc",
                      "svcMsPGOcc", "svcTMsPGOcc")
spocc_implemented_models <- c("PGOcc", "spPGOcc", "msPGOcc", "spMsPGOcc",
                              "lfMsPGOcc", "intPGOcc", "svcPGOcc")

parse_spocc_species <- function(species = NULL, y = NULL) {
  if (!is.null(species) && ncol(species) > 0) {
    cn <- tolower(names(species))
    idx <- which(cn %in% c("species", "species_id", "sp", "name", "taxon"))[1]
    if (is.na(idx)) idx <- 1
    out <- as.character(species[[idx]])
    out <- out[nzchar(out)]
    if (length(out)) return(out)
  }
  if (!is.null(y)) {
    prefixes <- sub("(__rep|_rep|\\.rep)[0-9]+$", "", colnames(y))
    prefixes <- prefixes[nzchar(prefixes)]
    if (length(unique(prefixes)) > 1) return(unique(prefixes))
  }
  character()
}

validate_spoccupancy_full <- function(y = NULL, occ = NULL, det = NULL, coords = NULL, species = NULL,
                                      integrated = NULL, newdata = NULL, newcoords = NULL, folds = NULL,
                                      model_type = "PGOcc", n.batch = 1000, batch.length = 25,
                                      n.burn = 0, n.thin = 1, n.chains = 1, n.factors = 3,
                                      NNGP = FALSE, n.neighbors = 15, k.fold = 0, svc.cols = "",
                                      occ.formula = "~ 1", det.formula = "~ 1", formula = "~ 1",
                                      data_structure = "single-species replicated") {
  ok <- TRUE
  msg <- character()
  if (!(model_type %in% c(spocc_implemented_models, "sfMsPGOcc", "lfJSDM", "sfJSDM",
                          "tPGOcc", "stPGOcc", "tMsPGOcc", "stMsPGOcc",
                          "svcPGBinom", "svcTPGBinom", "svcTPGOcc", "svcMsPGOcc", "svcTMsPGOcc",
                          "spIntPGOcc", "intMsPGOcc", "stIntPGOcc"))) {
    ok <- FALSE
    msg <- c(msg, paste0("Unknown spOccupancy model function: ", model_type))
  }
  if (!(model_type %in% spocc_implemented_models)) {
    msg <- c(msg, paste0(model_type, " is exposed in the UI but not yet enabled for automatic real fitting in this safe adapter. The adapter will fail with a real diagnostic instead of pretending success. Implemented real fits: ", paste(spocc_implemented_models, collapse = ", "), "."))
  }
  if (is.null(y)) {
    ok <- FALSE
    msg <- c(msg, "y.csv is required.")
  } else {
    ycheck <- suppressWarnings(as.matrix(as.data.frame(y, check.names = FALSE)))
    suppressWarnings(storage.mode(ycheck) <- "numeric")
    if (any(is.na(ycheck) & !is.na(as.matrix(y)))) {
      ok <- FALSE
      msg <- c(msg, "y.csv must contain numeric detection/nondetection values.")
    }
    if (model_type == "svcPGBinom") {
      if (any(ycheck < 0, na.rm = TRUE)) { ok <- FALSE; msg <- c(msg, "svcPGBinom y values must be non-negative counts.") }
    } else if (any(!(ycheck %in% c(0, 1, NA)))) {
      ok <- FALSE
      msg <- c(msg, "Occupancy y values must be 0/1/NA detection histories.")
    }
    if (model_type %in% spocc_multi_models) {
      sp <- parse_spocc_species(species, y)
      if (!length(sp)) {
        ok <- FALSE
        msg <- c(msg, "Multi-species spOccupancy models need species.csv or y columns named like sp1_rep1, sp1_rep2, ...")
      }
    }
  }
  if (!is.null(occ)) {
    occ_clean <- clean_predictor_types(occ, "spOccupancy occ.covs")
    factor_cols <- names(occ_clean)[vapply(occ_clean, is.factor, logical(1))]
    if (length(factor_cols)) msg <- c(msg, paste0("Occurrence categorical covariates will be converted to factors: ", paste(factor_cols, collapse = ", ")))
    fcheck <- validate_one_sided_formula(occ.formula, occ_clean, "spOccupancy occ.formula")
    ok <- ok && isTRUE(fcheck$ok)
    msg <- c(msg, fcheck$messages)
  } else if (!identical(trimws(occ.formula), "~ 1")) {
    ok <- FALSE
    msg <- c(msg, "occ.covs.csv is required when occ.formula is not ~ 1.")
  }
  if (!is.null(det)) {
    det_clean <- clean_predictor_types(det, "spOccupancy det.covs")
    factor_cols <- names(det_clean)[vapply(det_clean, is.factor, logical(1))]
    if (length(factor_cols)) msg <- c(msg, paste0("Detection categorical covariates will be converted to factors: ", paste(factor_cols, collapse = ", ")))
  } else if (!identical(trimws(det.formula), "~ 1")) {
    ok <- FALSE
    msg <- c(msg, "det.covs.csv is required when det.formula is not ~ 1.")
  }
  if (model_type %in% spocc_spatial_models || model_type %in% spocc_lf_models) {
    if (is.null(coords)) {
      ok <- FALSE
      msg <- c(msg, paste0(model_type, " requires coords.csv."))
    } else if (ncol(coords) < 2) {
      ok <- FALSE
      msg <- c(msg, "coords.csv must contain at least two coordinate columns.")
    }
  }
  if (model_type %in% spocc_integrated_models && is.null(integrated)) {
    ok <- FALSE
    msg <- c(msg, paste0(model_type, " requires integrated_sources.csv with at least source and prefix columns."))
  }
  total <- as.integer(n.batch) * as.integer(batch.length)
  if (!is.finite(total) || total < 10) msg <- c(msg, "Very small MCMC settings are for software checks only; publication analyses need much larger samples.")
  if (as.integer(n.burn) >= total) {
    ok <- FALSE
    msg <- c(msg, paste0("n.burn must be smaller than n.batch * batch.length = ", total, "."))
  }
  if (as.integer(n.thin) < 1 || ((total - as.integer(n.burn)) %% as.integer(n.thin)) != 0) {
    msg <- c(msg, "n.thin should divide n.batch * batch.length - n.burn for clean posterior storage.")
  }
  if (as.integer(n.chains) < 1) { ok <- FALSE; msg <- c(msg, "n.chains must be >= 1.") }
  if (model_type %in% spocc_lf_models && as.integer(n.factors) < 1) {
    ok <- FALSE
    msg <- c(msg, "Latent-factor occupancy models require n.factors >= 1.")
  }
  if (isTRUE(NNGP) && as.integer(n.neighbors) < 1) {
    ok <- FALSE
    msg <- c(msg, "NNGP requires n.neighbors >= 1.")
  }
  if (model_type %in% spocc_svc_models && !nzchar(trimws(svc.cols))) {
    msg <- c(msg, "SVC model selected but svc.cols is empty; the adapter will use the first occurrence design column.")
  }
  if (as.integer(k.fold) < 0) { ok <- FALSE; msg <- c(msg, "k.fold must be >= 0.") }
  if (length(msg) == 0) msg <- "spOccupancy data and settings check passed."
  list(ok = ok, messages = msg)
}

write_spoccupancy_reproducible_script <- function(outdir) {
  dir.create(file.path(outdir, "reproducible_script"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "workflow_scripts"), recursive = TRUE, showWarnings = FALSE)
  script <- c(
    "# Reproducible spOccupancy analysis script generated by JSDM Studio",
    "args <- commandArgs(trailingOnly = FALSE)",
    "file_arg <- '--file='",
    "script_arg <- args[startsWith(args, file_arg)]",
    "if (length(script_arg) > 0) setwd(dirname(dirname(normalizePath(sub(file_arg, '', script_arg[[1]]), winslash='/', mustWork=FALSE))))",
    "for (d in c('models','samples','tables','predictions','plots','diagnostics','standard','report','results','workflow_scripts','spatial','model_assessment')) dir.create(d, recursive=TRUE, showWarnings=FALSE)",
    "`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x",
    "as_bool <- function(x, default=FALSE) { if (is.null(x) || length(x)==0 || is.na(x)) return(default); if (is.logical(x)) return(isTRUE(x)); tolower(as.character(x)) %in% c('true','t','1','yes','y') }",
    "as_num <- function(x, default) { z <- suppressWarnings(as.numeric(x)); if (!length(z) || !is.finite(z[1])) default else z[1] }",
    "as_int <- function(x, default) as.integer(round(as_num(x, default)))",
    "step_warnings <- character()",
    "append_warning <- function(msg) step_warnings <<- unique(c(step_warnings, as.character(msg)))",
    "write_status <- function(status, errors=character()) { payload <- list(engine='spOccupancy', status=status, warnings=as.character(step_warnings), errors=as.character(errors), time=as.character(Sys.time())); if (requireNamespace('jsonlite', quietly=TRUE)) writeLines(jsonlite::toJSON(payload, pretty=TRUE, auto_unbox=TRUE), file.path('diagnostics','engine_status.json')); write.csv(data.frame(engine='spOccupancy', status=status, message=paste(c(step_warnings, errors), collapse='; ')), file.path('tables','engine_status.csv'), row.names=FALSE); try(writeLines(capture.output(sessionInfo()), file.path('diagnostics','session_info.txt')), silent=TRUE) }",
    "safe_step <- function(name, expr, required=FALSE) { tryCatch(withCallingHandlers(force(expr), warning=function(w){ append_warning(paste(name, conditionMessage(w), sep=': ')); invokeRestart('muffleWarning') }), error=function(e){ msg <- conditionMessage(e); append_warning(paste(name, msg, sep=': ')); writeLines(msg, file.path('diagnostics', paste0('spOccupancy_', gsub('[^A-Za-z0-9_]+','_',name), '_error.txt'))); if (required) stop(msg, call.=FALSE); NULL }) }",
    "read_csv_safe <- function(path) { if (!file.exists(path)) return(NULL); dat <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE); if (ncol(dat) > 1) { first <- dat[[1]]; first_name <- names(dat)[1] %||% ''; first_chr <- as.character(first); first_num <- suppressWarnings(as.numeric(first_chr)); sequence_index <- all(!is.na(first_num)) && identical(as.integer(first_num), seq_len(length(first_num))); row_id_name <- first_name %in% c('', 'X', 'X.1', '...1', 'row.names', 'rowname', 'row_id', 'id', 'site_id', 'sample_id'); row_id_text <- !all(!is.na(first_num)); if (!anyDuplicated(first_chr) && (row_id_name || sequence_index || row_id_text)) { dat <- dat[-1]; rownames(dat) <- make.unique(first_chr) } }; dat }",
    "clean_df <- function(df) { if (is.null(df)) return(NULL); df <- as.data.frame(df, check.names=FALSE, stringsAsFactors=FALSE); for (nm in names(df)) { if (is.character(df[[nm]])) { x <- trimws(df[[nm]]); x[x==''] <- NA; nx <- suppressWarnings(as.numeric(x)); df[[nm]] <- if (all(is.na(x) | !is.na(nx))) nx else factor(x, levels=unique(x[!is.na(x)])) } else if (is.logical(df[[nm]])) df[[nm]] <- factor(df[[nm]]) }; df }",
    "expand_formula <- function(ftxt, dat) { ftxt <- trimws(ftxt %||% '~ 1'); if (identical(ftxt, '~ .')) { vars <- names(dat); if (!length(vars)) return(as.formula('~ 1')); return(as.formula(paste('~', paste(vars, collapse=' + ')))) }; as.formula(ftxt) }",
    "numeric_matrix <- function(dat, name) { if (is.null(dat)) return(NULL); m <- as.matrix(dat); suppressWarnings(storage.mode(m) <- 'numeric'); if (any(is.na(m) & !is.na(as.matrix(dat)))) stop(name, ' contains non-numeric values.', call.=FALSE); m }",
    "split_rep_cols <- function(cols) { m <- regexec('^(.+?)(?:__rep|_rep|\\\\.rep)([0-9]+)$', cols, perl=TRUE); p <- regmatches(cols, m); species <- vapply(p, function(z) if (length(z) >= 3) z[2] else NA_character_, character(1)); rep <- suppressWarnings(as.integer(vapply(p, function(z) if (length(z) >= 3) z[3] else NA_character_, character(1)))); list(prefix=species, rep=rep) }",
    "species_from_table <- function(species_raw, y_raw) { if (!is.null(species_raw) && ncol(species_raw) > 0) { cn <- tolower(names(species_raw)); idx <- which(cn %in% c('species','species_id','sp','name','taxon'))[1]; if (is.na(idx)) idx <- 1; sp <- as.character(species_raw[[idx]]); sp <- sp[nzchar(sp)]; if (length(sp)) return(sp) }; spl <- split_rep_cols(colnames(y_raw)); unique(spl$prefix[!is.na(spl$prefix)]) }",
    "make_single_y <- function(y_raw) { y <- numeric_matrix(y_raw, 'y.csv'); storage.mode(y) <- 'numeric'; y }",
    "make_multi_y <- function(y_raw, species_raw=NULL) { ymat <- numeric_matrix(y_raw, 'multi-species y.csv'); sp <- species_from_table(species_raw, y_raw); spl <- split_rep_cols(colnames(y_raw)); if (length(sp) && all(!is.na(spl$prefix))) { reps <- sort(unique(spl$rep)); arr <- array(NA_real_, dim=c(length(sp), nrow(ymat), length(reps)), dimnames=list(sp, rownames(ymat), paste0('rep', reps))); for (s in sp) for (r in reps) { idx <- which(spl$prefix == s & spl$rep == r)[1]; if (!is.na(idx)) arr[s,,paste0('rep', r)] <- ymat[, idx] }; return(arr) }; if (!length(sp)) stop('Cannot infer species for multi-species y. Use species.csv or sp1_rep1 style column names.', call.=FALSE); R <- ncol(ymat) / length(sp); if (R %% 1 != 0) stop('Multi-species y columns must equal species * replicates.', call.=FALSE); arr <- array(ymat, dim=c(nrow(ymat), R, length(sp))); arr <- aperm(arr, c(3,1,2)); dimnames(arr) <- list(sp, rownames(ymat), paste0('rep', seq_len(R))); arr }",
    "make_det <- function(det_raw, J, R, formula_text='~ 1', prefix=NULL) { if (is.null(det_raw) || identical(trimws(formula_text), '~ 1')) return(list()); det <- clean_df(det_raw); if (!is.null(prefix)) { keep <- grepl(paste0('^', prefix, '(_|__|\\\\.)'), names(det)); if (any(keep)) { names(det)[keep] <- sub(paste0('^', prefix, '(_|__|\\\\.)'), '', names(det)[keep]); det <- det[, keep, drop=FALSE] } }; out <- list(); spl <- split_rep_cols(names(det)); if (nrow(det) == J && all(!is.na(spl$prefix))) { for (v in unique(spl$prefix)) { m <- matrix(NA_real_, J, R); for (r in seq_len(R)) { idx <- which(spl$prefix == v & spl$rep == r)[1]; if (!is.na(idx)) m[, r] <- as.numeric(det[[idx]]) }; out[[v]] <- m }; return(out) }; if (nrow(det) == J && ncol(det) == R) { out$obs <- numeric_matrix(det, 'det.covs.csv'); return(out) }; if (nrow(det) == J * R) { for (nm in names(det)) out[[nm]] <- matrix(as.numeric(det[[nm]]), nrow=J, ncol=R); return(out) }; if (nrow(det) == J) { for (nm in names(det)) out[[nm]] <- det[[nm]]; return(out) }; stop('det.covs.csv rows must be sites, sites*replicates, or sites x replicate columns.', call.=FALSE) }",
    "integrated_sources <- function(src_raw, y_raw) { if (!is.null(src_raw) && nrow(src_raw) > 0) { cn <- tolower(names(src_raw)); src_idx <- which(cn %in% c('source','source_id','data_source','dataset'))[1]; pre_idx <- which(cn %in% c('prefix','column_prefix'))[1]; if (is.na(src_idx)) src_idx <- 1; if (is.na(pre_idx)) pre_idx <- src_idx; return(data.frame(source=as.character(src_raw[[src_idx]]), prefix=as.character(src_raw[[pre_idx]]), stringsAsFactors=FALSE)) }; spl <- split_rep_cols(colnames(y_raw)); p <- unique(spl$prefix[!is.na(spl$prefix)]); data.frame(source=p, prefix=p, stringsAsFactors=FALSE) }",
    "make_integrated_data <- function(y_raw, occ, det_raw, src_raw) { ymat <- numeric_matrix(y_raw, 'integrated y.csv'); src <- integrated_sources(src_raw, y_raw); ylist <- list(); detlist <- list(); sites <- list(); for (i in seq_len(nrow(src))) { pre <- src$prefix[i]; spl <- split_rep_cols(colnames(y_raw)); idx <- which(spl$prefix == pre); if (!length(idx)) idx <- grep(paste0('^', pre), colnames(y_raw)); if (!length(idx)) stop('No y columns found for integrated source prefix: ', pre, call.=FALSE); ylist[[i]] <- ymat[, idx, drop=FALSE]; sites[[i]] <- seq_len(nrow(ymat)); detlist[[i]] <- make_det(det_raw, nrow(ymat), ncol(ylist[[i]]), '~ obs', prefix=pre) }; list(y=ylist, occ.covs=occ, det.covs=detlist, sites=sites, sources=src) }",
    "parse_svc_cols <- function(txt, X) { txt <- trimws(txt %||% ''); if (!nzchar(txt)) return(1L); vals <- trimws(unlist(strsplit(txt, ','))); nums <- suppressWarnings(as.integer(vals)); if (all(is.finite(nums))) return(nums); mm <- match(vals, colnames(X)); mm <- mm[is.finite(mm) & !is.na(mm)]; if (!length(mm)) 1L else mm }",
    "coord_mat <- function(coords_raw) { if (is.null(coords_raw)) return(NULL); m <- numeric_matrix(coords_raw[, seq_len(min(2, ncol(coords_raw))), drop=FALSE], 'coords.csv'); colnames(m) <- c('x','y')[seq_len(ncol(m))]; m }",
    "max_dist <- function(coords) { if (is.null(coords) || nrow(coords) < 2) return(3); max(stats::dist(coords), na.rm=TRUE) %||% 3 }",
    "single_priors <- function(spatial=FALSE, svc=FALSE, coords=NULL) { p <- list(beta.normal=list(mean=0, var=2.72), alpha.normal=list(mean=0, var=2.72)); if (spatial && !svc) { md <- max_dist(coords); p$sigma.sq.ig <- c(2,1); p$phi.unif <- c(max(0.1, md), max(1, md/0.1)) }; if (svc) { md <- max_dist(coords); p$sigma.sq.ig <- list(a=2,b=1); p$phi.unif <- list(a=max(0.1, md), b=max(1, md/0.1)) }; p }",
    "multi_priors <- function(spatial=FALSE, coords=NULL) { p <- list(beta.comm.normal=list(mean=0,var=2.72), alpha.comm.normal=list(mean=0,var=2.72), tau.sq.beta.ig=list(a=2,b=1), tau.sq.alpha.ig=list(a=2,b=1)); if (spatial) { md <- max_dist(coords); p$sigma.sq.ig <- list(a=2,b=1); p$phi.unif <- list(a=max(0.1, md), b=max(1, md/0.1)) }; p }",
    "int_priors <- function(n.data, p.occ=2, p.det=2) { list(beta.normal=list(rep(0,p.occ), rep(2.72,p.occ)), alpha.normal=list(rep(list(rep(0,p.det)), n.data), rep(list(rep(2.72,p.det)), n.data))) }",
    "summarize_samples <- function(x) { m <- tryCatch(as.matrix(x), error=function(e) NULL); if (is.null(m)) return(data.frame(parameter=NA_character_, mean=NA_real_, sd=NA_real_, q025=NA_real_, q50=NA_real_, q975=NA_real_)); data.frame(parameter=colnames(m) %||% paste0('p', seq_len(ncol(m))), mean=colMeans(m, na.rm=TRUE), sd=apply(m,2,sd,na.rm=TRUE), q025=apply(m,2,quantile,0.025,na.rm=TRUE), q50=apply(m,2,quantile,0.5,na.rm=TRUE), q975=apply(m,2,quantile,0.975,na.rm=TRUE), row.names=NULL) }",
    "write_any_csv <- function(x, path) { if (is.null(x)) return(FALSE); if (is.atomic(x) && length(dim(x)) > 2) x <- data.frame(index=seq_along(as.vector(x)), value=as.vector(x)); if (is.atomic(x) && is.null(dim(x))) x <- data.frame(name=names(x) %||% paste0('v', seq_along(x)), value=as.vector(x)); if (is.matrix(x) || is.data.frame(x)) write.csv(as.data.frame(x, check.names=FALSE), path) else { nms <- names(x) %||% paste0('item_', seq_along(x)); write.csv(data.frame(field=nms, value=vapply(x, function(z) paste(capture.output(str(z, max.level=1)), collapse=' '), character(1))), path, row.names=FALSE) }; TRUE }",
    "std_pred <- function(y, psi=NULL) { yy <- as.vector(y); n <- length(yy); pred <- if (!is.null(psi)) as.vector(apply(as.matrix(psi), 2, mean, na.rm=TRUE)) else rep(NA_real_, n); if (length(pred) != n) pred <- rep(NA_real_, n); data.frame(engine='spOccupancy', site_id=seq_len(n), response_id='occupancy', observed=yy, predicted_mean=pred, predicted_lower=NA_real_, predicted_upper=NA_real_) }",
    "try({ if (!requireNamespace('yaml', quietly=TRUE)) stop('yaml is required.'); cfg <- yaml::read_yaml('used_config.yml') }, silent=FALSE)",
    "tryCatch({",
    "  if (!requireNamespace('spOccupancy', quietly=TRUE)) stop('Package spOccupancy is required for real fitting.', call.=FALSE)",
    "  suppressPackageStartupMessages(library(spOccupancy))",
    "  set.seed(as_int(cfg$mcmc$seed, 1234))",
    "  y_raw <- read_csv_safe(file.path('data','y.csv')); occ <- clean_df(read_csv_safe(file.path('data','occ.covs.csv'))); det_raw <- read_csv_safe(file.path('data','det.covs.csv')); coords <- coord_mat(read_csv_safe(file.path('data','coords.csv'))); species_raw <- read_csv_safe(file.path('data','species.csv')); src_raw <- read_csv_safe(file.path('data','integrated_sources.csv')); newdata <- clean_df(read_csv_safe(file.path('data','newdata.csv'))); newcoords <- coord_mat(read_csv_safe(file.path('data','newcoords.csv')))",
    "  if (is.null(y_raw)) stop('y.csv is missing.', call.=FALSE); if (is.null(occ)) occ <- data.frame(int=rep(1, nrow(y_raw)))",
    "  model_type <- cfg$model$model_type %||% 'PGOcc'; implemented <- c('PGOcc','spPGOcc','msPGOcc','spMsPGOcc','lfMsPGOcc','intPGOcc','svcPGOcc'); if (!(model_type %in% implemented)) stop(paste0(model_type, ' is not enabled for automatic real fitting in this adapter.'), call.=FALSE)",
    "  occ.form <- expand_formula(cfg$model$occ.formula %||% '~ 1', occ); det.form <- expand_formula(cfg$model$det.formula %||% '~ 1', data.frame(obs=1)); n.samples <- as_int(cfg$mcmc$n.batch, 1000) * as_int(cfg$mcmc$batch.length, 25); n.burn <- as_int(cfg$mcmc$n.burn, max(1, floor(0.1*n.samples))); n.thin <- as_int(cfg$mcmc$n.thin, 1); n.chains <- as_int(cfg$mcmc$n.chains, 1); n.report <- as_int(cfg$mcmc$n.report, 100); n.threads <- as_int(cfg$mcmc$n.omp.threads, 1)",
    "  call_args <- list(); y_for_standard <- NULL",
    "  if (model_type %in% c('PGOcc','spPGOcc','svcPGOcc')) { y <- make_single_y(y_raw); J <- nrow(y); R <- ncol(y); det <- make_det(det_raw, J, R, cfg$model$det.formula %||% '~ 1'); data <- list(y=y, occ.covs=occ, det.covs=det); z <- apply(y,1,max,na.rm=TRUE); z[!is.finite(z)] <- 1; y_for_standard <- y; if (model_type %in% c('spPGOcc','svcPGOcc')) data$coords <- coords; if (model_type == 'PGOcc') { call_args <- list(occ.formula=occ.form, det.formula=det.form, data=data, inits=list(z=z), priors=single_priors(FALSE), n.samples=n.samples, n.burn=n.burn, n.thin=n.thin, n.chains=n.chains, n.omp.threads=n.threads, verbose=as_bool(cfg$mcmc$verbose, FALSE), n.report=n.report) } else if (model_type == 'spPGOcc') { call_args <- list(occ.formula=occ.form, det.formula=det.form, data=data, inits=list(z=z, sigma.sq=1, phi=max_dist(coords)), priors=single_priors(TRUE, FALSE, coords), tuning=list(phi=0.5), cov.model=cfg$spatial_latent_svc$cov.model %||% 'exponential', NNGP=as_bool(cfg$spatial_latent_svc$NNGP, TRUE), n.neighbors=as_int(cfg$spatial_latent_svc$n.neighbors, 15), search.type=cfg$spatial_latent_svc$search.type %||% 'cb', n.batch=as_int(cfg$mcmc$n.batch, 1000), batch.length=as_int(cfg$mcmc$batch.length, 25), accept.rate=as_num(cfg$mcmc$accept.rate, 0.43), n.burn=n.burn, n.thin=n.thin, n.chains=n.chains, n.omp.threads=n.threads, verbose=as_bool(cfg$mcmc$verbose, FALSE), n.report=n.report) } else { Xtmp <- model.matrix(occ.form, occ); svc.cols <- parse_svc_cols(cfg$spatial_latent_svc$svc.cols %||% '', Xtmp); call_args <- list(occ.formula=occ.form, det.formula=det.form, data=data, inits=list(z=z, sigma.sq=rep(1,length(svc.cols)), phi=rep(max_dist(coords), length(svc.cols))), priors=single_priors(TRUE, TRUE, coords), tuning=list(phi=rep(0.5, length(svc.cols))), svc.cols=svc.cols, cov.model=cfg$spatial_latent_svc$cov.model %||% 'exponential', NNGP=as_bool(cfg$spatial_latent_svc$NNGP, TRUE), n.neighbors=as_int(cfg$spatial_latent_svc$n.neighbors, 15), search.type=cfg$spatial_latent_svc$search.type %||% 'cb', n.batch=as_int(cfg$mcmc$n.batch, 1000), batch.length=as_int(cfg$mcmc$batch.length, 25), accept.rate=as_num(cfg$mcmc$accept.rate, 0.43), n.burn=n.burn, n.thin=n.thin, n.chains=n.chains, n.omp.threads=n.threads, verbose=as_bool(cfg$mcmc$verbose, FALSE), n.report=n.report) } }",
    "  if (model_type %in% c('msPGOcc','spMsPGOcc','lfMsPGOcc')) { y <- make_multi_y(y_raw, species_raw); J <- dim(y)[2]; R <- dim(y)[3]; det <- make_det(det_raw, J, R, cfg$model$det.formula %||% '~ 1'); data <- list(y=y, occ.covs=occ, det.covs=det); z <- apply(y,c(1,2),max,na.rm=TRUE); z[!is.finite(z)] <- 1; y_for_standard <- apply(y, c(2,3), max, na.rm=TRUE); if (model_type %in% c('spMsPGOcc','lfMsPGOcc')) data$coords <- coords; if (model_type == 'msPGOcc') { call_args <- list(occ.formula=occ.form, det.formula=det.form, data=data, inits=list(z=z), priors=multi_priors(FALSE), n.samples=n.samples, n.burn=n.burn, n.thin=n.thin, n.chains=n.chains, n.omp.threads=n.threads, verbose=as_bool(cfg$mcmc$verbose, FALSE), n.report=n.report) } else if (model_type == 'spMsPGOcc') { S <- dim(y)[1]; call_args <- list(occ.formula=occ.form, det.formula=det.form, data=data, inits=list(z=z, sigma.sq=rep(1,S), phi=rep(max_dist(coords),S)), priors=multi_priors(TRUE, coords), tuning=list(phi=rep(0.5,S)), cov.model=cfg$spatial_latent_svc$cov.model %||% 'exponential', NNGP=as_bool(cfg$spatial_latent_svc$NNGP, TRUE), n.neighbors=as_int(cfg$spatial_latent_svc$n.neighbors, 15), search.type=cfg$spatial_latent_svc$search.type %||% 'cb', n.batch=as_int(cfg$mcmc$n.batch, 1000), batch.length=as_int(cfg$mcmc$batch.length, 25), accept.rate=as_num(cfg$mcmc$accept.rate, 0.43), n.burn=n.burn, n.thin=n.thin, n.chains=n.chains, n.omp.threads=n.threads, verbose=as_bool(cfg$mcmc$verbose, FALSE), n.report=n.report) } else { call_args <- list(occ.formula=occ.form, det.formula=det.form, data=data, inits=list(z=z), priors=multi_priors(FALSE), n.factors=as_int(cfg$spatial_latent_svc$n.factors, 2), n.samples=n.samples, n.burn=n.burn, n.thin=n.thin, n.chains=n.chains, n.omp.threads=n.threads, verbose=as_bool(cfg$mcmc$verbose, FALSE), n.report=n.report) } }",
    "  if (model_type == 'intPGOcc') { idata <- make_integrated_data(y_raw, occ, det_raw, src_raw); pdet <- 2; det.forms <- rep(list(det.form), length(idata$y)); z <- apply(do.call(cbind, lapply(idata$y, function(m) apply(m,1,max,na.rm=TRUE))), 1, max); z[!is.finite(z)] <- 1; data <- list(y=idata$y, occ.covs=idata$occ.covs, det.covs=idata$det.covs, sites=idata$sites); call_args <- list(occ.formula=occ.form, det.formula=det.forms, data=data, inits=list(z=z), priors=int_priors(length(idata$y), ncol(model.matrix(occ.form, occ)), pdet), n.samples=n.samples, n.burn=n.burn, n.thin=n.thin, n.chains=n.chains, n.omp.threads=n.threads, verbose=as_bool(cfg$mcmc$verbose, FALSE), n.report=n.report); y_for_standard <- do.call(cbind, idata$y) }",
    "  saveRDS(call_args, file.path('models','spOccupancy_call_args.rds')); write.csv(data.frame(model_type=model_type, n.samples=n.samples, n.burn=n.burn, n.thin=n.thin, n.chains=n.chains), file.path('tables','model_settings_used.csv'), row.names=FALSE)",
    "  fun <- get(model_type, asNamespace('spOccupancy')); output <- do.call(fun, call_args)",
    "  if (as_bool(cfg$validation_prediction_outputs$save_model, TRUE)) saveRDS(output, file.path('models','spOccupancy_model.rds'))",
    "  capture.output(summary(output), file=file.path('results','spOccupancy_summary.txt')); write.csv(summarize_samples(output$beta.samples %||% output$beta.comm.samples), file.path('tables','summary_beta.csv'), row.names=FALSE); write.csv(summarize_samples(output$alpha.samples %||% output$alpha.comm.samples), file.path('tables','summary_alpha.csv'), row.names=FALSE)",
    "  sample_names <- names(output)[grepl('samples$', names(output))]; man <- data.frame(sample=character(), file=character(), rows=integer(), cols=integer()); if (as_bool(cfg$validation_prediction_outputs$save_samples, TRUE)) for (nm in sample_names) { f <- file.path('samples', paste0(nm, '.rds')); saveRDS(output[[nm]], f); m <- tryCatch(as.matrix(output[[nm]]), error=function(e) matrix(NA,0,0)); man <- rbind(man, data.frame(sample=nm, file=f, rows=nrow(m), cols=ncol(m))) }; write.csv(man, file.path('samples','posterior_sample_manifest.csv'), row.names=FALSE)",
    "  wa <- NULL; if (as_bool(cfg$validation_prediction_outputs$waicOcc, TRUE)) wa <- safe_step('waicOcc', spOccupancy::waicOcc(output)); if (!is.null(wa)) { write_any_csv(wa, file.path('model_assessment','waicOcc_results.csv')); write_any_csv(wa, file.path('tables','waicOcc_results.csv')) }",
    "  pp <- NULL; if (as_bool(cfg$validation_prediction_outputs$ppcOcc, FALSE)) pp <- safe_step('ppcOcc', spOccupancy::ppcOcc(output, fit.stat='freeman-tukey', group=1)); if (!is.null(pp)) { saveRDS(pp, file.path('model_assessment','ppcOcc.rds')); write_any_csv(pp$fit.y, file.path('model_assessment','ppc_fit_y.csv')); write_any_csv(pp$fit.y.rep, file.path('model_assessment','ppc_fit_y_rep.csv')) }",
    "  fit <- NULL; if (as_bool(cfg$validation_prediction_outputs$fitted, TRUE)) fit <- safe_step('fitted', fitted(output)); if (!is.null(fit)) { saveRDS(fit, file.path('results','fitted_values.rds')); for (nm in names(fit)) write_any_csv(fit[[nm]], file.path('tables', paste0('fitted_', nm, '.csv'))) }",
    "  pred <- NULL; if (as_bool(cfg$validation_prediction_outputs$predict, TRUE)) { if (identical(model_type, 'lfMsPGOcc') && (is.null(newdata) || is.null(newcoords))) { append_warning('predict.lfMsPGOcc skipped because it requires newdata.csv and newcoords.csv for genuinely new locations; fitted psi samples were saved instead') } else pred <- safe_step('predict', { nd <- newdata %||% occ; X0 <- model.matrix(occ.form, nd); if (model_type %in% c('spPGOcc','spMsPGOcc','lfMsPGOcc','svcPGOcc')) { c0 <- newcoords %||% coords; if (model_type == 'svcPGOcc') predict(output, X.0=X0, coords.0=c0, weights.0=rep(1,nrow(X0)), verbose=FALSE, n.report=n.report) else if (model_type == 'lfMsPGOcc') predict(output, X.0=X0, coords.0=c0) else predict(output, X.0=X0, coords.0=c0, verbose=FALSE, n.report=n.report) } else predict(output, X.0=X0) }) }; if (!is.null(pred)) { saveRDS(pred, file.path('predictions','spOccupancy_predict.rds')); for (nm in names(pred)) write_any_csv(pred[[nm]], file.path('predictions', paste0(nm, '.csv'))) }",
    "  if (!is.null(output$theta.samples)) write.csv(summarize_samples(output$theta.samples), file.path('spatial','spatial_theta_summary.csv'), row.names=FALSE); if (!is.null(output$w.samples)) saveRDS(output$w.samples, file.path('spatial','spatial_or_factor_w_samples.rds')); if (!is.null(output$lambda.samples)) write.csv(summarize_samples(output$lambda.samples), file.path('tables','latent_factor_loadings.csv'), row.names=FALSE)",
    "  beta_eff <- read.csv(file.path('tables','summary_beta.csv'), check.names=FALSE); effects <- data.frame(engine='spOccupancy', response_id='occupancy', predictor=beta_eff$parameter, direction=ifelse(beta_eff$mean>0,'positive',ifelse(beta_eff$mean<0,'negative','zero')), estimate=beta_eff$mean, lower=beta_eff$q025, upper=beta_eff$q975, notes=paste(model_type, 'posterior mean'), stringsAsFactors=FALSE); write.csv(effects, file.path('standard','effects_long.csv'), row.names=FALSE)",
    "  psi <- output$psi.samples %||% NULL; write.csv(std_pred(y_for_standard, psi), file.path('standard','predictions_long.csv'), row.names=FALSE); write.csv(data.frame(engine='spOccupancy', response_1=NA_character_, response_2=NA_character_, association_type=ifelse(!is.null(output$lambda.samples),'latent_factor','not_available'), estimate=NA_real_, comparable_level='model_specific'), file.path('standard','associations_long.csv'), row.names=FALSE)",
    "  if (!is.null(wa)) { wam <- as.numeric(unlist(wa)); fitm <- data.frame(engine='spOccupancy', metric=names(wam) %||% paste0('waic_', seq_along(wam)), response_id=NA_character_, value=wam, notes='waicOcc') } else fitm <- data.frame(engine='spOccupancy', metric='fit_available', response_id=NA_character_, value=1, notes='waic skipped'); write.csv(fitm, file.path('standard','fit_metrics.csv'), row.names=FALSE)",
    "  write.csv(data.frame(run_id=basename(getwd()), engine='spOccupancy', status='fitted', model_type=model_type, n_sites=ifelse(is.null(y_for_standard), NA_integer_, nrow(as.matrix(y_for_standard))), n_replicates=ifelse(is.null(y_for_standard), NA_integer_, ncol(as.matrix(y_for_standard))), n.samples=n.samples, n.burn=n.burn, n.thin=n.thin, n.chains=n.chains, stringsAsFactors=FALSE), file.path('standard','run_summary.csv'), row.names=FALSE)",
    "  for (f in c('run_summary.csv','fit_metrics.csv','effects_long.csv','predictions_long.csv','associations_long.csv')) file.copy(file.path('standard', f), file.path('results', f), overwrite=TRUE)",
    "  write.csv(data.frame(file=list.files('.', recursive=TRUE), stringsAsFactors=FALSE), file.path('standard','output_manifest.csv'), row.names=FALSE)",
    "  writeLines(c('spOccupancy results','===================', paste0('Status: fitted'), paste0('Model: ', model_type), '', 'Important files:', '- models/spOccupancy_model.rds', '- samples/posterior_sample_manifest.csv', '- tables/summary_beta.csv and summary_alpha.csv', '- model_assessment/waicOcc_results.csv when enabled', '- diagnostics/engine_status.json', '- reproducible_script/run_this_spOccupancy_analysis.R'), file.path('results','README_spOccupancy_results.txt'))",
    "  writeLines(c('<!DOCTYPE html><html><head><meta charset=\"UTF-8\"><title>spOccupancy Report</title></head><body><h1>spOccupancy workflow report</h1>', paste0('<p>Status: fitted</p><p>Model: ', model_type, '</p>'), '<p>Open diagnostics/engine_status.json and standard/*.csv for machine-readable outputs.</p></body></html>'), file.path('report','spOccupancy_report.html'))",
    "  writeLines('RUN COMPLETE', 'RUN_COMPLETE.txt'); write_status('fitted')",
    "}, error=function(e) {",
    "  msg <- conditionMessage(e); writeLines(msg, file.path('diagnostics','spOccupancy_reproducible_error.txt')); write.csv(data.frame(engine='spOccupancy', response_id=NA_character_, predictor=NA_character_, direction=NA_character_, estimate=NA_real_, lower=NA_real_, upper=NA_real_, notes='fit_failed'), file.path('standard','effects_long.csv'), row.names=FALSE); write.csv(data.frame(engine='spOccupancy', site_id=NA_character_, response_id=NA_character_, observed=NA_real_, predicted_mean=NA_real_, predicted_lower=NA_real_, predicted_upper=NA_real_), file.path('standard','predictions_long.csv'), row.names=FALSE); write.csv(data.frame(engine='spOccupancy', response_1=NA_character_, response_2=NA_character_, association_type=NA_character_, estimate=NA_real_, comparable_level='fit_failed'), file.path('standard','associations_long.csv'), row.names=FALSE); write.csv(data.frame(engine='spOccupancy', metric=NA_character_, response_id=NA_character_, value=NA_real_, notes='fit_failed'), file.path('standard','fit_metrics.csv'), row.names=FALSE); write.csv(data.frame(run_id=basename(getwd()), engine='spOccupancy', status='fit_failed'), file.path('standard','run_summary.csv'), row.names=FALSE); writeLines(c('spOccupancy results','===================','Status: fit_failed','','Open diagnostics/spOccupancy_reproducible_error.txt and diagnostics/engine_status.json first.'), file.path('results','README_spOccupancy_results.txt')); write_status('fit_failed', msg); writeLines('RUN FAILED', 'RUN_FAILED.txt'); quit(status=1, save='no')",
    "})"
  )
  writeLines(script, file.path(outdir, "reproducible_script", "run_this_spOccupancy_analysis.R"))
  writeLines(c("# Workflow wrapper for spOccupancy",
               "source(file.path('reproducible_script', 'run_this_spOccupancy_analysis.R'))"),
             file.path(outdir, "workflow_scripts", "run_spOccupancy_workflow.R"))
}
