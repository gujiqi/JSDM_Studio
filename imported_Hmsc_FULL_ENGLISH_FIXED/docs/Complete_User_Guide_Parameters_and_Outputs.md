# JSDM Studio Complete User Guide

# JSDM Studio Parameter Guide

This guide explains every major input file, model setting, MCMC option, and output option in JSDM Studio.

JSDM Studio is a graphical workflow built on top of the R package **Hmsc**. It does not change the Hmsc statistical model; it helps users run Hmsc-based joint species distribution models through a graphical and reproducible workflow.

---

## 1. Input files

### 1.1 Y.csv: species response matrix

**Meaning:** `Y.csv` is the main community matrix.

```text
Rows    = sampling units / plots / sites / samples
Columns = species
Cells   = species responses
```

Examples of cell values:

| Data type | Recommended distr | Meaning |
|---|---|---|
| 0/1 | probit | absence / presence |
| 0, 1, 2, 3... | poisson | counts |
| continuous values | normal | biomass, cover, continuous abundance index |

Important checks:

```text
1. The rows of Y must match the rows of XData.
2. Species names should be valid column names.
3. If traits are used, Y column names must match TrData row names.
4. If a phylogeny is used, Y column names must match tree tip labels.
5. Remove completely empty species columns before formal analysis.
```

---

### 1.2 XData.csv: environmental variables

**Meaning:** `XData.csv` contains plot-level environmental variables or predictors.

```text
Rows    = sampling units
Columns = environmental variables
```

Examples:

```text
pH
moisture
canopy
elevation
temperature
precipitation
substrate
habitat
land_use
```

Important checks:

```text
1. Variable names must match names used in XFormula.
2. Avoid spaces and special characters in column names.
3. Categorical variables should have consistent spelling.
4. Highly correlated variables may make interpretation difficult.
5. Missing values can cause model failure.
```

---

### 1.3 TrData.csv: species traits, optional

**Meaning:** `TrData.csv` contains species-level traits.

```text
Rows    = species
Columns = traits
```

Examples:

```text
body_size
height
life_form
seed_mass
dispersal_type
reproductive_mode
feeding_guild
```

Use traits when you want to ask:

```text
Do species traits explain variation in environmental responses?
```

Important checks:

```text
1. Row names must match species names in Y.
2. Large amounts of missing trait data should be handled before modelling.
3. Too many traits with too few species can make the model unstable.
```

---

### 1.4 studyDesign.csv: random-effect grouping, optional

**Meaning:** `studyDesign.csv` defines grouping variables for random effects.

```text
Rows    = sampling units
Columns = grouping variables
```

Examples:

```text
site
plot
transect
region
year
observer
sample
```

Use random effects when samples are not fully independent, for example:

```text
multiple plots within the same site
repeated surveys across years
nested sampling design
regional clustering
```

---

### 1.5 coordinates.csv: spatial coordinates, optional

**Meaning:** `coordinates.csv` gives the spatial position of each sampling unit.

```text
Rows    = sampling units
Columns = x/y or longitude/latitude
```

Use coordinates when spatial autocorrelation is expected.

Important checks:

```text
1. Coordinate rows must match Y rows.
2. Use consistent coordinate units.
3. Spatial models are slower and more complex.
```

---

### 1.6 phyloTree: phylogenetic tree, optional

Use a phylogenetic tree if you want to model phylogenetic structure in species responses.

Important checks:

```text
1. Tree tip labels must match species names in Y.
2. Name mismatch is the most common error.
3. For a first test run, leave phylogeny off.
```

---

## 2. Model settings

### 2.1 distr: response distribution

| distr | Suitable data | Example |
|---|---|---|
| probit | presence/absence | 0/1 occurrence |
| poisson | counts | number of individuals |
| normal | continuous values | biomass, cover, abundance index |

Choose the distribution according to the values in `Y.csv`.

---

### 2.2 XFormula: environmental formula

Examples:

```r
~ .
~ pH + moisture + canopy + elevation
~ pH + I(pH^2) + moisture
~ pH * moisture
```

Meaning:

```text
~ . uses all columns in XData.
~ pH + moisture uses only selected predictors.
I(pH^2) adds a quadratic term.
pH * moisture adds main effects and interaction.
```

Caution:

```text
Do not use ~ . if XData contains ID columns that should not be predictors.
```

---

### 2.3 TrFormula: trait formula

Examples:

```r
~ .
~ height + life_form + reproductive_mode
```

`TrFormula` controls which traits explain variation in species environmental responses.

---

### 2.4 random_mode

| random_mode | Meaning |
|---|---|
| none | no random effect |
| sample | sample/group random effect |
| spatial | spatial random effect |

For a first test, use `sample` or `none`. Use `spatial` only when coordinates are available and spatial autocorrelation matters.

---

### 2.5 random_effect_column

This is the column in `studyDesign.csv` used for the random effect.

Examples:

```text
site
plot
region
year
sample
```

The name must exactly match a column in `studyDesign.csv`.

---

### 2.6 spatial_method

| Method | Meaning |
|---|---|
| Full | full spatial random effect; suitable for small datasets |
| NNGP | nearest-neighbour Gaussian process approximation |
| GPP | Gaussian predictive process approximation |

Spatial models can be slow. Do not start with spatial models when only testing the program.

---

### 2.7 nNeighbours

Used for NNGP spatial approximation.

Common values:

```text
10
15
20
```

Larger values are more detailed but slower.

---

## 3. MCMC settings

### 3.1 samples

Number of retained posterior samples.

```text
Small test: 20–100
Formal analysis: often 1000–5000 or more
```

Small values are only for testing the software.

---

### 3.2 transient

Burn-in iterations discarded before posterior sampling.

```text
Small test: 10–100
Formal analysis: often 1000 or more
```

Use traceplots and convergence diagnostics to decide whether burn-in is sufficient.

---

### 3.3 thin

Thinning interval.

```text
thin = 1 saves every sample
thin = 10 saves every 10th sample
```

---

### 3.4 nChains

Number of MCMC chains.

```text
Testing: 2
Formal analysis: at least 2, often 4
```

Multiple chains are needed for convergence diagnostics such as PSRF/Rhat.

---

### 3.5 nParallel

Number of parallel workers.

```text
Windows beginners: 1
Stronger computers: 2 or 4
```

If parallel errors occur, set it back to 1.

---

### 3.6 verbose

Console progress interval for Hmsc MCMC.

Smaller values show progress more frequently.

---

### 3.7 seed

Random seed for reproducibility.

Examples:

```text
123
2024
1
```

Record the seed in formal analysis.

---

## 4. Output options

### Save model

Saves fitted model objects such as:

```text
models/hmsc_model.rds
models/models_thin_*_samples_*_chains_*.RData
```

Always save models for formal analysis.

---

### Compute predicted values

Generates model predictions for observed sampling units.

Main output:

```text
results/predicted_values.rds
```

---

### Explanatory model fit

Evaluates how well the model explains the observed data.

Main outputs:

```text
results/model_fit_explanatory.txt
tables/model_fit_explanatory_*.csv
```

This is not independent predictive performance.

---

### Cross-validation

Evaluates predictive performance on held-out data.

Main outputs:

```text
results/model_fit_cross_validation.txt
tables/model_fit_cross_validation_*.csv
```

More reliable for prediction, but slower.

---

### WAIC

WAIC is useful for comparing alternative models. A single WAIC value alone is not very informative.

---

### MCMC diagnostics

Main outputs:

```text
results/MCMC_convergence.txt
plots/MCMC_traceplots.pdf
```

Formal analysis must check convergence. A model finishing successfully does not automatically mean it is reliable.

---

### Beta, Gamma, and Omega

```text
Beta  = environmental effects on species
Gamma = trait effects on environmental responses
Omega = residual species associations
```

Important: Omega is not direct evidence of competition or facilitation. It is a residual correlation after accounting for measured predictors and random effects.

---

### Variance partitioning

Shows how much variation is attributed to different predictor groups or random effects.

---

### Environmental gradient predictions

Shows predicted species responses along environmental gradients. Avoid extrapolating beyond the observed data range.


---

# JSDM Studio Output File Guide

This guide explains the files generated by JSDM Studio.

Results are usually saved in:

```text
output/hmsc_date_time/
```

---

## 1. Top-level files

### used_config.yml

Complete configuration used for the run.

Contains:

```text
input file names
distr
XFormula
TrFormula
random effect settings
MCMC settings
output options
seed
```

This is essential for reproducibility.

---

### RUN_COMPLETE.txt

A completion marker showing that the workflow finished.

---

### JSDM_Studio_report.html / HmscGUI_report.html

An HTML summary report with run information and output links.

---

## 2. inputs/ folder

Contains copies of input files used for this run.

Interpretation:

```text
Y.csv: rows = sampling units, columns = species
XData.csv: rows = sampling units, columns = environmental variables
TrData.csv: rows = species, columns = traits
studyDesign.csv: rows = sampling units, columns = grouping variables
coordinates.csv: rows = sampling units, columns = coordinates
```

Keeping input copies is important for reproducibility.

---

## 3. models/ folder

### unfitted_model.rds

Hmsc model object before MCMC fitting.

### hmsc_model.rds

Fitted Hmsc model object. This is the main model file.

Read it in R:

```r
m <- readRDS("models/hmsc_model.rds")
```

### models_thin_*_samples_*_chains_*.RData

Fitted model saved in RData format. The file name records MCMC settings.

---

## 4. results/ folder

### model_structure.txt

Text description of the model structure.

### predicted_values.rds

Predicted values object.

Expected interpretation:

```text
rows    = sampling units
columns = species
values  = predicted responses
```

The exact internal structure depends on Hmsc output.

### model_fit_explanatory.txt

Explanatory model fit.

Typical table interpretation:

```text
rows    = species
columns = fit metrics
```

For presence/absence data, metrics may include AUC and Tjur R2. For continuous data, metrics may include RMSE and R2.

### model_fit_cross_validation.txt

Cross-validation results.

Typical interpretation:

```text
rows    = species
columns = predictive performance metrics
```

### WAIC.txt

WAIC model comparison result. Useful when comparing multiple models.

### MCMC_convergence.txt

MCMC convergence diagnostics.

Look for:

```text
PSRF/Rhat close to 1
sufficient effective sample size
reasonable chain behaviour
```

### species_associations.rds

Residual species association object.

If converted to a matrix:

```text
rows    = species A
columns = species B
values  = residual association
```

Positive values indicate residual co-occurrence. Negative values indicate residual segregation. Do not interpret these directly as biotic interactions without additional evidence.

---

## 5. tables/ folder

### data_dimensions.csv

Summary of input dimensions.

### species_summary.csv

Species-level summary.

Typical interpretation:

```text
rows    = species
columns = occurrence count, occurrence frequency, mean abundance, missing values, etc.
```

### environment_summary.csv

Environmental variable summary.

Typical interpretation:

```text
rows    = environmental variables
columns = mean, standard deviation, min, max, missing values, variable type
```

### model_fit_explanatory_*.csv

Explanatory model fit table.

```text
rows    = species
columns = explanatory fit metrics
```

### model_fit_cross_validation_*.csv

Cross-validation table.

```text
rows    = species
columns = predictive performance metrics
```

### parameter_estimates_Beta_*.csv

Beta parameter estimates.

Typical interpretation:

```text
rows    = species-environment combinations or species
columns = environmental predictors, posterior means, support, uncertainty
```

Beta describes how environmental variables affect species responses.

### parameter_estimates_Gamma_*.csv

Gamma parameter estimates.

Typical interpretation:

```text
rows    = trait-response combinations
columns = posterior means, support, uncertainty
```

Gamma describes how traits explain species-specific environmental responses.

### variance_partitioning_*.csv

Variance partitioning table.

Typical interpretation:

```text
rows    = species
columns = explanatory groups
values  = proportions of explained variation
```

### Omega_*_*.csv

Residual association matrix.

```text
rows    = species A
columns = species B
cell    = residual association between species A and B
```

---

## 6. plots/ folder

### MCMC_traceplots.pdf

Traceplots for MCMC chains. Good traceplots show stable mixing and no strong trends.

### model_fit_explanatory_vs_predictive.pdf

Compares explanatory and predictive performance.

### Beta_plot.pdf

Environmental response plot.

### Gamma_plot.pdf

Trait-response plot.

### variance_partitioning.pdf

Variance partitioning plot.

### Omega_associations_*.pdf

Residual species association plot.

### predictions_environmental_gradients.pdf

Predicted species responses along environmental gradients.

---

## 7. output ZIP

The downloaded result ZIP contains the full run folder. Save it for collaboration, archiving, and reproducibility.
