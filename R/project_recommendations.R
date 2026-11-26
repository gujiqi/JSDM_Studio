# Project matching uses the workflows wired in Studio, not the full upstream APIs.
studio_project_questions <- c(
  "Community response to environmental gradients" = "community",
  "Management / treatment effect" = "community",
  "BACI disturbance impact" = "baci",
  "Restoration trajectory" = "community",
  "Invasion impact" = "community",
  "Multi-stressor response" = "community",
  "Diversity and ecosystem function response" = "community",
  "Indicator / ecosystem health assessment" = "community",
  "Mixed-scale biodiversity / attribute analysis" = "mixed",
  "Spatial community structure" = "spatial",
  "Spatial occupancy with imperfect detection" = "spatial_occ",
  "Multi-species occupancy with replicated surveys" = "multi_occ",
  "Integrated occupancy with multiple data sources" = "integrated",
  "Trait-mediated community response" = "traits",
  "Phylogeny-informed community response" = "phylogeny",
  "Species association / residual correlation" = "association",
  "Model-based ordination" = "ordination",
  "Large eDNA / metabarcoding community matrix" = "large",
  "Microbiome / OTU community response" = "large",
  "Composition or proportional response data" = "composition",
  "Mixed presence, count and continuous responses" = "mixed",
  "Zero-heavy abundance or median-zero response matrix" = "zeros",
  "Prediction at new sites" = "prediction",
  "Inverse prediction / environmental reconstruction" = "inverse",
  "Variable selection / sparse environmental effects" = "sparse",
  "Bayesian fourth-corner trait analysis" = "fourth_corner",
  "Detection probability and false absence" = "occupancy",
  "Temporal / multi-season occupancy" = "multi_season",
  "Spatially varying coefficients" = "svc",
  "Custom hypothesis-driven Bayesian community regression" = "community"
)

studio_project_families <- c(
  "Not specified" = "unknown", "Binary occurrence (0/1)" = "binary",
  "Counts (non-negative integers)" = "count", "Continuous (Gaussian scale)" = "continuous",
  "Ordered categories" = "ordinal", "Proportions strictly between 0 and 1" = "proportion"
)

studio_project_engines <- c("Hmsc", "Hmsc-HPC", "jSDM", "GJAM", "spOccupancy", "sjSDM", "boral")

studio_project_response_options <- function(family) {
  switch(family,
    binary = c(Hmsc = "probit", `Hmsc-HPC` = "probit", jSDM = "binomial_probit or binomial_logit",
               GJAM = "typeName PA", sjSDM = "binomial (probit or logit)", boral = "binomial with trial size 1"),
    count = c(Hmsc = "poisson", `Hmsc-HPC` = "poisson", jSDM = "poisson_log", GJAM = "typeName DA",
              sjSDM = "poisson (log link)", boral = "poisson, or negative binomial if justified"),
    continuous = c(Hmsc = "normal", `Hmsc-HPC` = "normal", jSDM = "gaussian", GJAM = "typeName CON",
                   sjSDM = "gaussian (identity link)", boral = "normal"),
    ordinal = c(GJAM = "typeName OC", boral = "ordinal; supply the ordered integer levels"),
    proportion = c(boral = "beta; each value must be strictly between 0 and 1"),
    mixed = c(GJAM = "Supply one typeName per response column.",
              boral = "Supply a family vector; every column must have a supported boral likelihood."),
    composition = c(GJAM = "Use FC (fractional) or CC (count composition) and define the compositional groups."),
    character()
  )
}

studio_project_goal_options <- function(goal, family) {
  community <- c(
    Hmsc = "Hierarchical species responses with optional traits, tree and random levels.",
    `Hmsc-HPC` = "Python HMSC for a supported subset of the Hmsc model.",
    jSDM = "Bayesian species responses with optional latent factors.",
    GJAM = "Joint regression on the observation scale.",
    sjSDM = "Regularised environmental effects and species covariance.",
    boral = "Bayesian multivariate regression with optional latent variables."
  )
  if (goal == "community") {
    direct <- switch(family, ordinal = c("GJAM", "boral"), proportion = "boral",
                     mixed = "GJAM", composition = "GJAM", c("Hmsc", "jSDM"))
    return(list(direct = direct, why = community))
  }
  if (goal %in% c("mixed", "composition")) {
    return(list(direct = "GJAM", why = if (goal == "mixed") community[c("GJAM", "boral")] else community["GJAM"]))
  }
  why <- switch(goal,
    baci = c(Hmsc = "Fit a treatment-by-period interaction and grouping for repeated sites.",
             `Hmsc-HPC` = "A treatment-by-period contrast with supported random levels only."),
    traits = c(Hmsc = "Trait-by-environment effects on species responses (Gamma).",
               `Hmsc-HPC` = "Trait-by-environment effects for the supported HMSC subset.",
               boral = "Trait effects through a fourth-corner model."),
    fourth_corner = c(boral = "Bayesian fourth-corner modelling with a species trait table.",
                      Hmsc = "Trait-by-environment Gamma addresses a related question with a different model."),
    phylogeny = c(Hmsc = "Models phylogenetic structure in species responses and estimates rho.",
                  `Hmsc-HPC` = "Phylogeny-informed HMSC within the supported Python subset."),
    spatial = c(Hmsc = "Spatial random levels for residual community structure.",
                `Hmsc-HPC` = "Full spatial random intercepts; not NNGP or GPP.",
                boral = "Spatial latent effects using distances among sites.",
                sjSDM = "Spatial predictors or eigenvectors; not a spatial random-effect model."),
    ordination = c(boral = "Unconstrained model-based ordination with latent variables.",
                   GJAM = "Ordination of fitted covariance; not the same as latent site scores."),
    association = c(Hmsc = "Random-level Omega; add a suitable random level.",
                    `Hmsc-HPC` = "Associations derived from loadings for the supported random levels.",
                    jSDM = "Residual correlations from fitted latent factors.",
                    GJAM = "corMu on the GJAM observation-model scale.",
                    sjSDM = "Regularised species association structure.",
                    boral = "Residual correlations from fitted latent variables."),
    prediction = c(Hmsc = "Project to new sites and carry through the fitted spatial structure when specified.",
                   GJAM = "Predict responses at new covariates on the GJAM scale.",
                   sjSDM = "Predict from fitted environmental and, if used, supplied spatial predictors.",
                   `Hmsc-HPC` = "Fixed-effect projections; not conditional spatial predictions."),
    large = if (family %in% c("mixed", "composition")) {
      c(GJAM = "Preserves the declared mixed or compositional response structure.")
    } else {
      c(sjSDM = "Regularised covariance for same-family community matrices.",
        Hmsc = "Hierarchical alternative if trait, tree or random-level effects matter.",
        `Hmsc-HPC` = "Consider for a supported HMSC design; assess runtime on your data.")
    },
    zeros = c(GJAM = "Consider a suitable count or continuous observation model; zeros alone do not choose typeName.",
              boral = "For counts, consider negative binomial if overdispersion is supported."),
    inverse = c(GJAM = "Reconstruct environmental predictors from observed responses within the fitted analysis."),
    sparse = c(boral = "Bayesian stochastic-search variable selection (SSVS).",
               sjSDM = "Penalised effects are an alternative, not posterior inclusion probabilities."),
    character()
  )
  direct <- switch(goal, prediction = c("Hmsc", "GJAM", "sjSDM"),
                   association = c("Hmsc", "jSDM", "sjSDM", "boral"), names(why)[1L])
  list(direct = direct, why = why)
}

studio_project_recommend <- function(profile = list(), question = "") {
  engine_names <- studio_project_engines
  goal <- unname(studio_project_questions[match(question, names(studio_project_questions))])
  if (length(goal) != 1L || is.na(goal)) goal <- "community"
  get_flag <- function(key) isTRUE(profile[[key]])
  response <- profile$response_structure
  if (!is.character(response) || length(response) != 1L || is.na(response)) response <- "unknown"
  structure <- switch(response, "single response type" = "single",
    "detection-nondetection occupancy surveys" = "replicated",
    "mixed response types" = "mixed", "presence/count/continuous mixture" = "mixed",
    "composition data" = "composition", "unknown")
  family <- if (structure == "single") profile$response_family else structure
  if (length(family) != 1L || is.na(family) || !(family %in% c(unname(studio_project_families), "mixed", "composition", "replicated"))) {
    family <- "unknown"
  }
  species <- suppressWarnings(as.numeric(profile$n_responses))
  valid_species <- length(species) == 1L && is.finite(species) && species == floor(species) && species >= 1
  coords <- get_flag("has_spatial")
  occ_goals <- c("occupancy", "spatial_occ", "multi_occ", "integrated", "multi_season", "svc")
  notices <- character()
  blocked <- NULL
  pending <- NULL
  if (goal == "multi_season") {
    blocked <- "The multi-season spOccupancy functions are not connected in Studio."
  } else if (structure == "unknown") {
    pending <- "Specify the observed response structure before choosing a workflow."
  } else if (goal %in% occ_goals && structure != "replicated") {
    pending <- "An occupancy or detection model needs real repeated detection/non-detection surveys. A single 0/1 record per site is not enough."
  } else if (structure == "replicated" && !goal %in% c(occ_goals, "community", "spatial", "prediction", "association", "baci")) {
    blocked <- "No connected workflow addresses this goal jointly with imperfect detection."
  } else if (goal == "mixed" && structure != "mixed") {
    pending <- "This question needs different response types in one matrix; confirm the response structure."
  } else if (goal == "composition" && !family %in% c("composition", "proportion")) {
    pending <- "Distinguish a joint composition (e.g. fractions adding to one) from individual bounded proportions."
  } else if (goal == "zeros" && !family %in% c("count", "continuous", "unknown")) {
    pending <- "This template is for zero-heavy count or continuous abundance, not simple incidence or composition."
  } else if (goal %in% c("traits", "fourth_corner") && !get_flag("has_traits")) {
    pending <- "Declare a matched species trait table for trait-mediated inference."
  } else if (goal == "phylogeny" && !get_flag("has_phylogeny")) {
    pending <- "Declare a species-aligned phylogeny or taxonomy table for this question."
  } else if (goal %in% c("integrated", "svc") && valid_species && species != 1) {
    blocked <- if (goal == "integrated") "Studio connects intPGOcc for one species; multi-species integration is not connected." else
      "Studio connects svcPGOcc for one species; multi-species spatially varying coefficients are not connected."
  } else if (goal %in% c("spatial", "spatial_occ", "svc") && !coords) {
    pending <- "Declare training coordinates for a spatial model."
  } else if (goal == "integrated" && !get_flag("has_integrated_sources")) {
    pending <- "Integrated occupancy needs multiple aligned survey sources and their source/prefix metadata."
  } else if (!valid_species) {
    pending <- "Enter the number of response columns or surveyed species."
  } else if (species < 2 &&
             (structure %in% c("mixed", "composition") || goal %in% c("traits", "fourth_corner", "phylogeny", "ordination", "association", "multi_occ"))) {
    pending <- "This question or response structure needs at least two species or response columns."
  } else if (structure == "replicated" && goal == "association" && !coords) {
    pending <- "Studio's lfMsPGOcc validator requires coordinates, although lfMsPGOcc is a non-spatial latent-factor model."
  }

  if (goal == "baci") {
    notices <- c(notices, "BACI requires before/after observations and control/impact sites. A score does not validate the design.")
  }
  if (goal == "zeros" || get_flag("many_zeros") && family == "count") {
    notices <- c(notices, "Many zeros alone do not establish zero inflation. Check effort, dispersion and the observation process.")
  }
  if (get_flag("need_prediction") || goal == "prediction") {
    notices <- c(notices, "New-site prediction requires matching predictor columns and, for spatial predictions, the appropriate spatial inputs.")
  }
  if (goal == "inverse") {
    notices <- c(notices, "Studio's inverse predictions use the fitted response records, not a separate upload of new responses.")
  }

  confirmed <- family != "unknown"
  response_options <- studio_project_response_options(family)
  task <- if (goal == "composition" && family == "proportion") {
    list(direct = "boral", why = c(boral = "Model individual bounded proportions, not a joint composition."))
  } else studio_project_goal_options(goal, family)
  limitation <- c(Hmsc = "Check family, random levels and MCMC convergence.",
    `Hmsc-HPC` = "CPU TensorFlow; iid and Full spatial levels only.",
    jSDM = "Latent-factor residual correlations are not direct interactions.",
    GJAM = "Match typeNames to observation scales and units.",
    spOccupancy = "Check visit histories, site alignment and detection covariates.",
    sjSDM = "PyTorch required. Regularisation is not an MCMC posterior.",
    boral = "JAGS required. Check the latent structure and convergence.")
  trait_engines <- c("Hmsc", "Hmsc-HPC", "boral")
  phylo_engines <- c("Hmsc", "Hmsc-HPC")
  spatial_engines <- c("Hmsc", "Hmsc-HPC", "boral", "sjSDM")
  new_site_engines <- c("Hmsc", "GJAM", "sjSDM", "spOccupancy")

  assess <- function(engine) {
    row <- function(status, why, score = NA_integer_, route = "", points = "", direct = FALSE) {
      data.frame(Engine = engine, Status = status, Score = score, Why = why, Route = route,
                 Points = points, Direct = direct, stringsAsFactors = FALSE)
    }
    if (!is.null(blocked)) return(row("unavailable", blocked))
    if (structure == "replicated" || goal %in% occ_goals) {
      if (engine != "spOccupancy") return(row("not_applicable", "Does not fit a separate detection process for this question."))
      if (!is.null(pending)) return(row("needs_input", pending))
      model <- if (goal == "integrated") "intPGOcc" else if (goal == "svc") "svcPGOcc" else
        if (goal == "association") "lfMsPGOcc" else if (species == 1) {
          if (coords) "spPGOcc" else "PGOcc"
        } else if (coords) "spMsPGOcc" else "msPGOcc"
      predict_bonus <- if (get_flag("need_prediction") && goal != "prediction" && model %in% c("PGOcc", "spPGOcc", "msPGOcc", "spMsPGOcc")) 1L else 0L
      why <- if (model == "lfMsPGOcc") "This non-spatial latent-factor model separates occupancy and detection. Coordinates are required by Studio's input check." else
        "Models occupancy and detection separately using the declared repeat surveys."
      return(row("direct", why, 8L + predict_bonus, paste("Choose", model),
                 paste0("+2 detection data; +6 objective", if (predict_bonus) "; +1 new-site prediction"), TRUE))
    }
    if (engine == "spOccupancy") {
      return(row("not_applicable", "Requires repeated detection/non-detection histories, not a site-by-species matrix."))
    }
    if (!is.null(pending)) return(row("needs_input", pending))
    if (confirmed && !(engine %in% names(response_options))) {
      return(row("not_applicable", "This Studio workflow does not fit the declared response structure or family."))
    }
    if (!(engine %in% names(task$why))) {
      return(row("not_applicable", "Not a match for the selected objective in Studio."))
    }
    if (goal == "zeros" && engine == "boral" && family == "continuous") {
      return(row("not_applicable", "An ordinary continuous likelihood does not address an extra mass at zero."))
    }
    if (family == "mixed" && engine == "boral") {
      return(row("needs_input", "Each response column must have a supported boral likelihood before scoring a mixed matrix.",
                 route = response_options[[engine]]))
    }
    direct <- engine %in% task$direct
    task_score <- if (direct) 6L else 2L
    parts <- stats::setNames(task_score, if (direct) "direct objective" else "alternative objective")
    if (confirmed) parts <- c("response match" = 2L, parts)
    if (!goal %in% c("traits", "fourth_corner") && get_flag("has_traits") && engine %in% trait_engines) parts <- c(parts, traits = 1L)
    if (goal != "phylogeny" && get_flag("has_phylogeny") && engine %in% phylo_engines) parts <- c(parts, phylogeny = 1L)
    if (goal != "spatial" && coords && engine %in% spatial_engines) parts <- c(parts, "spatial input" = 1L)
    if (goal != "prediction" && get_flag("need_prediction") && engine %in% new_site_engines) parts <- c(parts, "new-site prediction" = 1L)
    status <- if (!confirmed) "provisional" else if (direct) "direct" else "alternative"
    route <- if (confirmed) response_options[[engine]] else "Choose a response family before selecting a likelihood."
    if (confirmed && goal == "inverse") route <- paste0(route, ". Enable PREDICTX; inspect xpredMu and xpredSd.")
    if (confirmed && goal == "zeros" && engine == "GJAM" && family == "continuous") {
      route <- "Check how zeros arise before choosing a continuous typeName."
    }
    row(status, task$why[[engine]], sum(parts), route,
        paste(sprintf("+%d %s", parts, names(parts)), collapse = "; "), direct)
  }

  rows <- do.call(rbind, lapply(engine_names, assess))
  rows$Caution <- unname(limitation[rows$Engine])
  tier <- match(rows$Status, c("direct", "alternative", "provisional", "needs_input", "not_applicable", "unavailable"))
  rows <- rows[order(tier, -rows$Score, seq_len(nrow(rows)), na.last = TRUE), , drop = FALSE]
  leading <- rows$Status == "direct"
  top <- if (any(leading)) rows[leading & rows$Score == max(rows$Score[leading]), , drop = FALSE] else rows[FALSE, , drop = FALSE]
  message <- if (!is.null(blocked)) blocked else if (!is.null(pending)) pending else
    if (!confirmed) "Choose a response family. The displayed points are provisional." else
      if (!nrow(top) && any(rows$Status == "alternative")) "No direct match. Compare the assumptions of the alternatives below." else
        if (!nrow(top)) "No scored workflow matches both the declared response and the objective." else ""
  list(rows = rows, top = top, message = message, notices = notices, goal = goal)
}
