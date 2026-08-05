# JSDM Studio Quick Start

## Installation for users

1. Install R for Windows:
   https://cran.r-project.org/bin/windows/base/

2. Install `JSDMStudio_Setup.exe`.

3. Open the desktop icon:

```text
JSDM Studio
```

4. On first launch, wait while required R packages are checked or installed.

5. If the window does not open, use:

```text
JSDM Studio Browser Mode
```

---

## First test settings

Use example data first.

Data:

```text
Y.csv     = examples/Y.csv
XData.csv = examples/XData.csv
traits    = leave empty
studyDesign = leave empty unless the example provides it
coordinates = leave empty
phyloTree = leave empty
```

Model settings:

```text
distr = probit
XFormula = ~ pH + moisture + canopy + elevation
random_mode = sample
random_effect_column = sample
```

MCMC settings:

```text
samples = 20
transient = 10
thin = 1
nChains = 2
nParallel = 1
verbose = 5
seed = 123
```

For the first test, enable only:

```text
Save model
Compute predicted values
Explanatory model fit
MCMC diagnostics
Generate HTML report
```

Disable slow outputs such as WAIC, Omega associations, and environmental gradient predictions until the basic run succeeds.
