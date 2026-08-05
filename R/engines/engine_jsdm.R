
# R/engines/engine_jsdm.R
# jSDM engine adapter for JSDMWorkbench.
# This adapter keeps jSDM separate from Hmsc while using the same input/output workflow.

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
}

jsdm_make_output_dirs <- function(outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  for (d in c("models", "results", "plots", "tables", "inputs", "diagnostics")) {
    dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
  }
}

jsdm_write_lines <- function(x, file, append = FALSE) {
  cat(paste(x, collapse = "\n"), file = file, sep = "\n", append = append)
}

jsdm_capture_to_file <- function(expr, file) {
  con <- file(file, open = "wt", encoding = "UTF-8")
  sink(con)
  on.exit({
    try(sink(), silent = TRUE)
    try(close(con), silent = TRUE)
  }, add = TRUE)
  force(expr)
}

jsdm_extract_beta_table <- function(mod) {
  if (is.null(mod$mcmc.sp)) return(NULL)
  species_names <- names(mod$mcmc.sp)
  rows <- list()
  for (sp in species_names) {
    m <- mod$mcmc.sp[[sp]]
    if (is.null(m)) next
    sm <- try(summary(m)[[1]], silent = TRUE)
    if (inherits(sm, "try-error")) next
    pars <- rownames(sm)
    beta_idx <- grepl("^beta", pars)
    if (!any(beta_idx)) next
    out <- data.frame(
      species = sp,
      parameter = pars[beta_idx],
      mean = sm[beta_idx, "Mean"],
      sd = sm[beta_idx, "SD"],
      row.names = NULL,
      check.names = FALSE
    )
    rows[[length(rows) + 1]] <- out
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

jsdm_extract_gamma_table <- function(mod) {
  if (is.null(mod$mcmc.gamma)) return(NULL)
  rows <- list()
  for (nm in names(mod$mcmc.gamma)) {
    m <- mod$mcmc.gamma[[nm]]
    if (is.null(m)) next
    sm <- try(summary(m)[[1]], silent = TRUE)
    if (inherits(sm, "try-error")) next
    out <- data.frame(
      block = nm,
      parameter = rownames(sm),
      mean = sm[, "Mean"],
      sd = sm[, "SD"],
      row.names = NULL,
      check.names = FALSE
    )
    rows[[length(rows) + 1]] <- out
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

fit_jsdm_workflow <- function(Y, XData, TrData = NULL, cfg = list(), outdir,
                              log_fun = function(...) message(...),
                              progress_fun = function(percent, stage, step = NULL, detail = "") {}) {
  jsdm_make_output_dirs(outdir)
  progress_fun(5, "Preparing jSDM", 1, "Checking package and inputs")
  if (!requireNamespace("jSDM", quietly = TRUE)) {
    stop("The R package 'jSDM' is not installed. Run install_packages.bat or install.packages('jSDM'). On some systems GNU GSL/Rtools may be required.")
  }
  if (!requireNamespace("coda", quietly = TRUE)) stop("Package 'coda' is required for jSDM diagnostics.")

  Y <- as.data.frame(Y, check.names = FALSE)
  XData <- as.data.frame(XData, check.names = FALSE)
  ymat <- as.matrix(Y)
  if (nrow(ymat) != nrow(XData)) stop("Y and XData must have the same number of rows.")

  jcfg <- cfg$jsdm %||% cfg
  model_type <- jcfg$model_type %||% "binomial_probit"
  site_formula <- as.formula(jcfg$site_formula %||% "~ .")
  trait_formula <- NULL
  if (isTRUE(jcfg$use_traits) && !is.null(TrData)) trait_formula <- as.formula(jcfg$trait_formula %||% "~ .")

  common <- list(
    burnin = as.integer(jcfg$burnin %||% 100),
    mcmc = as.integer(jcfg$mcmc %||% 100),
    thin = as.integer(jcfg$thin %||% 1),
    site_formula = site_formula,
    site_data = XData,
    n_latent = as.integer(jcfg$n_latent %||% 2),
    site_effect = jcfg$site_effect %||% "random",
    beta_start = as.numeric(jcfg$beta_start %||% 0),
    gamma_start = as.numeric(jcfg$gamma_start %||% 0),
    lambda_start = as.numeric(jcfg$lambda_start %||% 0),
    W_start = as.numeric(jcfg$W_start %||% 0),
    alpha_start = as.numeric(jcfg$alpha_start %||% 0),
    V_alpha = as.numeric(jcfg$V_alpha %||% 1),
    shape_Valpha = as.numeric(jcfg$shape_Valpha %||% 0.5),
    rate_Valpha = as.numeric(jcfg$rate_Valpha %||% 0.0005),
    mu_beta = as.numeric(jcfg$mu_beta %||% 0),
    V_beta = as.numeric(jcfg$V_beta %||% 10),
    mu_gamma = as.numeric(jcfg$mu_gamma %||% 0),
    V_gamma = as.numeric(jcfg$V_gamma %||% 10),
    mu_lambda = as.numeric(jcfg$mu_lambda %||% 0),
    V_lambda = as.numeric(jcfg$V_lambda %||% 10),
    seed = as.integer(jcfg$seed %||% 1234),
    verbose = as.integer(jcfg$verbose %||% 1)
  )
  if (isTRUE(jcfg$use_traits) && !is.null(TrData)) {
    common$trait_data <- as.data.frame(TrData, check.names = FALSE)
    common$trait_formula <- trait_formula
  }

  progress_fun(15, "Fitting jSDM", 2, paste("Model type:", model_type))
  log_fun(paste("Running jSDM engine:", model_type))

  mod <- switch(model_type,
    binomial_probit = do.call(jSDM::jSDM_binomial_probit, c(list(presence_data = ymat), common)),
    binomial_logit = {
      common$ropt <- as.numeric(jcfg$ropt %||% 0.44)
      trials_col <- jcfg$trials_column %||% ""
      if (nzchar(trials_col) && trials_col %in% colnames(XData)) common$trials <- XData[[trials_col]]
      do.call(jSDM::jSDM_binomial_logit, c(list(presence_data = ymat), common))
    },
    binomial_probit_sp_constrained = do.call(jSDM::jSDM_binomial_probit_sp_constrained, c(list(presence_data = ymat), common)),
    poisson_log = {
      common$ropt <- as.numeric(jcfg$ropt %||% 0.44)
      do.call(jSDM::jSDM_poisson_log, c(list(count_data = ymat), common))
    },
    gaussian = {
      common$V_start <- as.numeric(jcfg$V_start %||% 1)
      common$shape_V <- as.numeric(jcfg$shape_V %||% 0.5)
      common$rate_V <- as.numeric(jcfg$rate_V %||% 0.0005)
      do.call(jSDM::jSDM_gaussian, c(list(response_data = ymat), common))
    },
    stop("Unknown jSDM model_type: ", model_type)
  )

  progress_fun(65, "Saving jSDM model", 3, "Writing model object and summaries")
  saveRDS(mod, file.path(outdir, "models", "jsdm_model.rds"))
  jsdm_capture_to_file(print(mod), file.path(outdir, "results", "model_print.txt"))
  jsdm_capture_to_file(str(mod, max.level = 2), file.path(outdir, "results", "model_structure.txt"))

  beta_tab <- jsdm_extract_beta_table(mod)
  if (!is.null(beta_tab)) write.csv(beta_tab, file.path(outdir, "tables", "jsdm_beta_estimates.csv"), row.names = FALSE)
  gamma_tab <- jsdm_extract_gamma_table(mod)
  if (!is.null(gamma_tab)) write.csv(gamma_tab, file.path(outdir, "tables", "jsdm_gamma_estimates.csv"), row.names = FALSE)

  progress_fun(75, "jSDM correlations", 4, "Residual and environmental correlations")
  if (isTRUE(jcfg$compute_residual_cor) && !is.null(mod$model_spec$n_latent) && mod$model_spec$n_latent > 0) {
    rescor <- try(jSDM::get_residual_cor(mod, prob = as.numeric(jcfg$cor_prob %||% 0.95), type = "mean"), silent = TRUE)
    if (!inherits(rescor, "try-error")) {
      saveRDS(rescor, file.path(outdir, "results", "residual_correlation.rds"))
      if (!is.null(rescor$cor.mean)) write.csv(rescor$cor.mean, file.path(outdir, "tables", "residual_correlation_mean.csv"))
      pdf(file.path(outdir, "plots", "residual_correlation.pdf"), width = 8, height = 8)
      try(jSDM::plot_residual_cor(mod, prob = as.numeric(jcfg$cor_prob %||% 0.95)), silent = TRUE)
      dev.off()
    }
  }
  if (isTRUE(jcfg$compute_enviro_cor)) {
    envcor <- try(jSDM::get_enviro_cor(mod, prob = as.numeric(jcfg$cor_prob %||% 0.95), type = "mean"), silent = TRUE)
    if (!inherits(envcor, "try-error")) {
      saveRDS(envcor, file.path(outdir, "results", "environmental_correlation.rds"))
      if (!is.null(envcor$cor)) write.csv(envcor$cor, file.path(outdir, "tables", "environmental_correlation_mean.csv"))
    }
  }

  progress_fun(85, "jSDM predictions", 5, "Saving fitted predictions if available")
  pred_name <- intersect(c("theta_latent", "probit_theta_latent", "logit_theta_latent", "log_theta_latent"), names(mod))
  for (pn in pred_name) {
    obj <- mod[[pn]]
    if (!is.null(obj)) {
      saveRDS(obj, file.path(outdir, "results", paste0(pn, ".rds")))
      if (is.matrix(obj) || is.data.frame(obj)) write.csv(obj, file.path(outdir, "tables", paste0(pn, ".csv")))
    }
  }

  progress_fun(92, "jSDM diagnostics", 6, "Traceplots for first species")
  if (isTRUE(jcfg$make_traceplots) && !is.null(mod$mcmc.sp)) {
    pdf(file.path(outdir, "diagnostics", "traceplots_first_species.pdf"), width = 10, height = 8)
    try({
      par(mfrow = c(2, 2))
      sp1 <- mod$mcmc.sp[[1]]
      np <- min(4, ncol(sp1))
      for (i in seq_len(np)) coda::traceplot(coda::as.mcmc(sp1[, i]), main = colnames(sp1)[i])
    }, silent = TRUE)
    dev.off()
  }

  jsdm_write_lines(c(
    "jSDM run complete",
    paste("Model type:", model_type),
    paste("Rows:", nrow(ymat)),
    paste("Responses/species:", ncol(ymat)),
    paste("Predictors:", ncol(XData)),
    paste("Latent variables:", common$n_latent),
    paste("Site effect:", common$site_effect)
  ), file.path(outdir, "RUN_COMPLETE.txt"))

  if (isTRUE(jcfg$make_report)) {
    html <- paste0(
      "<html><head><meta charset='utf-8'><title>jSDM report</title>",
      "<style>body{font-family:Arial,sans-serif;margin:32px;line-height:1.55}code{background:#f1f5f9;padding:2px 5px;border-radius:4px}</style>",
      "</head><body><h1>jSDM engine report</h1>",
      "<p>This run used the <b>jSDM</b> engine inside JSDMWorkbench.</p>",
      "<h2>Model</h2><ul>",
      "<li>Model type: ", model_type, "</li>",
      "<li>Rows: ", nrow(ymat), "</li>",
      "<li>Responses/species: ", ncol(ymat), "</li>",
      "<li>Latent variables: ", common$n_latent, "</li>",
      "<li>Site effect: ", common$site_effect, "</li>",
      "</ul><h2>Important output files</h2><ul>",
      "<li><code>models/jsdm_model.rds</code>: fitted jSDM object</li>",
      "<li><code>tables/jsdm_beta_estimates.csv</code>: species-level environmental effects</li>",
      "<li><code>tables/jsdm_gamma_estimates.csv</code>: trait-environment effects, if traits were used</li>",
      "<li><code>tables/residual_correlation_mean.csv</code>: residual correlations from latent variables, if requested</li>",
      "</ul></body></html>"
    )
    writeLines(html, file.path(outdir, "report.html"), useBytes = TRUE)
  }
  progress_fun(100, "jSDM complete", 7, "jSDM output generated")
  invisible(mod)
}
