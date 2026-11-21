# ODRCC: Outcome-Dependent Right Censored Covariates

ODRCC provides simulation tools and estimation routines for regression models that include outcome-dependent right-censored covariates (ODRCC). The package accompanies the manuscript:

Vazquez, J. et al. (2025). "Robust Estimation Under Outcome-Dependent Right Censoring in Huntington Disease: Estimators for Low and High Censoring Rates." 

---

## Key Features

The package implements several estimation methods for regression parameters under ODRCC:

| Function                      | Estimator                                                        |
|-------------------------------|------------------------------------------------------------------|
| `estimate_beta_ipw()`         | IPW estimator                                                    |
| `estimate_beta_mle()`         | Joint maximum likelihood estimator                               |
| `estimate_beta_aipw()`        | Unguarded AIPW estimator (`fx_method = "joint"` or `"direct"`)   |
| `estimate_beta_aipw_lambda()` | Guarded AIPW estimator                                           |
| `estimate_beta_cc()`          | Complete-case estimator (for comparison)                         |

Every estimator estimates theta = (beta, sigma) and returns sandwich standard errors.

Additional tools included:

- `ess_weights()`: effective sample size of the IPW weights, to guide the choice of estimator
- `data_aft()`: the data generator of the first simulation study in the manuscript
- Weibull AFT models for X | Z and C | (Y, Z), with any number of covariates

---

## Simulation code

The code for the simulation studies and diagnostics in the manuscript (Figure 1, Table 2,
the heat maps of the second simulation study, the sensitivity analysis, and the diagnostics in the Supplementary
Material) is in [`simulations/`](simulations/). See [`simulations/README.md`](simulations/README.md)
for which folder produces which result and how to run it.

---

## Installation

```
# install.packages("devtools")  # if needed
devtools::install_github(
  "jesusepfvazquez/dep-right-censored-covariate",
  subdir = "ODRCC",
  dependencies = TRUE,
  build_vignettes = TRUE
)
```
---

## Tutorial

For a tutorial please check out the vignette which is available here:

https://jesusepfvazquez.github.io/dep-right-censored-covariate/

It can also be accessed through here:

```
# option 1
browseVignettes("ODRCC") # pick "ODRCC-tutorial"

# Or try 
vignette("ODRCC-tutorial")
```
