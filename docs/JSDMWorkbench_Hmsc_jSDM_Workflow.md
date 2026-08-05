
# JSDMWorkbench workflow design: Hmsc + jSDM

This package reorganizes the original JSDM Studio workflow into an engine-based workbench.

## Main workflow

1. **Home**: overview and recommended workflow.
2. **Project**: define the ecological question and choose single-engine or multi-engine mode.
3. **Data**: upload response matrix `Y.csv`, predictor table `XData.csv`, and optional trait, study design, coordinates and phylogeny files.
4. **Check**: validate dimensions and identifiers before modelling.
5. **Engine**: choose `Hmsc`, `jSDM`, or both.
6. **Hmsc Settings**: all original Hmsc parameters are kept and exposed.
7. **jSDM Settings**: jSDM model type, latent variables, site effect, priors, starting values and output options.
8. **Outputs**: shared standardized output choices.
9. **Run & Results**: run selected engine(s), view progress and download a ZIP archive.
10. **Guides**: detailed Hmsc and jSDM explanations.

## Engine design

The workbench follows this principle:

```text
Unified workflow, engine-specific settings.
```

The GUI standardizes input, output, reports and ZIP packaging. The statistical engines keep their own model-specific parameters.

## Engine folders

Each run creates a project output folder and then engine-specific subfolders:

```text
output/jsdm_workbench_YYYYMMDD_HHMMSS/
├── used_config.yml
├── inputs/
├── engines/
│   ├── Hmsc/
│   │   ├── models/
│   │   ├── results/
│   │   ├── tables/
│   │   ├── plots/
│   │   └── report.html
│   └── jSDM/
│       ├── models/
│       ├── results/
│       ├── tables/
│       ├── plots/
│       ├── diagnostics/
│       └── report.html
└── comparison/
```

## When to use Hmsc

Use Hmsc when the goal is an interpretable Bayesian community model with traits, phylogeny, random effects, variance partitioning, predictions, model fit and residual associations.

## When to use jSDM

Use jSDM when the goal is a basic Bayesian JSDM / latent-variable multivariate regression with residual correlations and environmental correlations. jSDM supports binomial logit/probit, Poisson log and Gaussian response models.

## Important installation note for jSDM

The jSDM package uses C++ code and links to GSL through RcppGSL. On Windows, CRAN binary installation is usually easiest. If source installation is required, users may need Rtools and GNU GSL.
