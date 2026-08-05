# Engine registry for JSDM Studio.
# This is a developer-facing summary of the engines exposed by the GUI.

jsdmw_available_engines <- function() {
  list(
    hmsc = list(
      id = "hmsc",
      name = "Hmsc",
      label = "Interpretable Bayesian community model",
      language = "R",
      package = "Hmsc",
      best_for = c("traits", "phylogeny", "random effects", "spatial effects", "variance partitioning", "Omega associations", "S1-S7 reproducible Hmsc workflow"),
      response_types = c("probit/binary", "poisson/count", "normal/continuous"),
      native_statuses = c("fitted", "check_failed", "fit_failed"),
      status = "production workflow when Hmsc is installed"
    ),
    hmschpc = list(
      id = "hmschpc",
      name = "Hmsc-HPC",
      label = "CPU pyhmsc / HDF5 HMSC workflow",
      language = "R + Python + TensorFlow",
      package = "pyhmsc / hmsc-hpc source",
      best_for = c("CPU Python-native HMSC", "HDF5 posterior files", "traits", "phylogenetic covariance or Newick", "iid and spatial_full random intercepts", "Hmsc-style S1-S7 output folder"),
      response_types = c("probit/binary", "poisson/count", "normal/continuous"),
      native_statuses = c("fitted", "model_defined", "check_failed", "fit_failed"),
      status = "CPU workflow; GPU/Slurm not exposed"
    ),
    jsdm = list(
      id = "jsdm",
      name = "jSDM",
      label = "Bayesian latent-variable JSDM",
      language = "R + C++",
      package = "jSDM",
      best_for = c("basic Bayesian JSDM", "latent variables", "residual correlations", "environmental correlations", "presence-absence/count/continuous responses"),
      response_types = c("binomial probit", "binomial logit", "poisson log", "gaussian"),
      native_statuses = c("fitted", "model_defined", "check_failed", "fit_failed"),
      status = "production workflow when jSDM is installed"
    ),
    gjam = list(
      id = "gjam",
      name = "GJAM",
      label = "Generalized joint attribute model",
      language = "R",
      package = "gjam",
      best_for = c("mixed response scales", "median-zero data", "observation-scale inference", "prediction", "inverse prediction", "sensitivity"),
      response_types = c("PA", "CON", "CA", "DA", "FC", "CC", "OC", "CAT"),
      native_statuses = c("fitted", "model_defined", "check_failed", "fit_failed"),
      status = "production workflow when gjam is installed"
    ),
    spoccupancy = list(
      id = "spoccupancy",
      name = "spOccupancy",
      label = "Occupancy / imperfect-detection workflow",
      language = "R + C++",
      package = "spOccupancy",
      best_for = c("detection-nondetection data", "imperfect detection", "replicated surveys", "spatial occupancy", "multi-species occupancy", "integrated data"),
      response_types = c("detection/nondetection occupancy arrays"),
      native_statuses = c("fitted", "model_defined", "check_failed", "fit_failed"),
      status = "production workflow when spOccupancy is installed"
    ),
    sjsdm = list(
      id = "sjsdm",
      name = "sjSDM",
      label = "Scalable PyTorch joint species distribution model",
      language = "R + Python/PyTorch via reticulate",
      package = "sjSDM",
      best_for = c("large community matrices", "regularized associations", "environment/spatial/biotic partitioning", "eDNA/OTU data", "fast covariance modelling"),
      response_types = c("binomial", "poisson", "gaussian", "negative binomial"),
      native_statuses = c("fitted", "model_defined", "check_failed", "fit_failed"),
      status = "production workflow when sjSDM and PyTorch are installed"
    ),
    boral = list(
      id = "boral",
      name = "boral",
      label = "Bayesian ordination and regression analysis",
      language = "R + JAGS",
      package = "boral",
      best_for = c("model-based ordination", "latent variables", "fourth-corner traits", "SSVS variable selection", "residual correlations", "JAGS MCMC"),
      response_types = c("binomial", "poisson", "negative binomial", "normal", "tweedie", "gamma", "lognormal", "beta", "ordinal", "zero-truncated"),
      native_statuses = c("fitted", "model_defined", "check_failed", "fit_failed"),
      status = "production workflow when boral, rjags/R2jags and system JAGS are installed"
    )
  )
}
