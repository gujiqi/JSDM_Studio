# JSDMWorkbench Nature UI SAFE Build

This build redesigns the interface so it is closer to a high-quality publication-grade workbench.

## Main design principles

1. Shared workflow, separate engines.
2. Hmsc keeps its detailed parameter settings.
3. jSDM is added as an independent engine with its own model type, latent variables, site effects, priors, starting values and MCMC settings.
4. Shared settings are only conceptual and are translated separately.
5. The UI should open safely even when optional modelling packages are not installed.
6. The output folder is always auditable: configuration, inputs, checks, engine status, report and ZIP.

## Navigation

Home
Project
Data
Check
Engine
Shared Settings
Hmsc Settings
jSDM Settings
Outputs
Run & Results
Guides

## Runtime safety note

This SAFE build creates output folders, saves configurations and writes engine status files without crashing if Hmsc or jSDM is missing. Full production fitting should be connected through engine adapters in R/engines after package installation has been confirmed on the target Windows system.

## Why this is safer

The previous design risked mixing Hmsc and jSDM parameters. This version separates:

- Shared Settings
- Hmsc Settings
- jSDM Settings

This avoids pretending that Hmsc and jSDM parameters are identical.
