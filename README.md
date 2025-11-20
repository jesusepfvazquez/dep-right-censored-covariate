# ODRCC: Outcome-Dependent Right Censored Covariates

ODRCC provides simulation tools and estimation routines for regression models that include outcome-dependent right-censored covariates (ODRCC). The package accompanies the manuscript:

Vazquez, J. et al. (2025). "Robust Estimation Under Outcome-Dependent Right Censoring in Huntington Disease: Estimators for Low and High Censoring Rates." 

---

## Key Features

The package implements several estimation methods for regression parameters under ODRCC:

| Estimator      | Description                                         |
|----------------|-----------------------------------------------------|
| CC             | Complete-case estimator                             |
| IPW            | Inverse probability weighting estimator             |
| MLE            | Maximum likelihood estimator                        |
| AIPW           | Augmented inverse probability weighting             |
| AIPW-Lambda    | Closed-form AIPW estimator                          |

Additional tools included:

- Data generation via `data_aft()`
- Modeling for X|Z and C|(Y,Z)
- Robust sandwich variance estimation
- Example workflows and simulation templates

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
