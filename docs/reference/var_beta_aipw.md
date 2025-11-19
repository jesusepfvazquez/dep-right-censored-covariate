# Sandwich variance estimator for the AIPW estimator

Computes a plug-in sandwich variance estimator for the AIPW regression
estimator under outcome-dependent right-censoring of a covariate. The
outcome model is specified by `model`, typically of the form
`y ~ AW + Z1 + ... + Zp`, where `AW = A - X`.

## Usage

``` r
var_beta_aipw(data_yXZ, theta, lbound = 0, ubound = 50)
```

## Arguments

- data_yXZ:

  Data frame containing at least:

  - `y`: outcome,

  - `A`: auxiliary covariate used to form `AW = A - X`,

  - `W`: observed covariate `W = min(X, C)`,

  - `D`: indicator `I(X <= C)`,

  - all covariates appearing in `model`,

  - all covariates appearing in `model_weights` and `model_xz`.

- theta:

  Numeric vector `c(beta, psi)` from the AIPW estimator, where `beta`
  has length equal to the number of columns in
  `model.matrix(model, data_yXZ)` and `psi` is the residual standard
  deviation.

- lbound, ubound:

  Numeric lower and upper bounds for the numerical integration over `X`
  in the augmentation term (defaults: 0, 50).

## Value

A list with components

- beta_est:

  Estimated regression coefficients \\\beta\\.

- psi_est:

  Estimated residual standard deviation \\\psi\\.

- se_beta:

  Sandwich standard errors for \\\beta\\.

## Details

The AIPW estimator combines:

- an IPW component based on a censoring model for `C | (Y, Z, ...)`
  specified by `model_weights`, and

- an augmentation component based on an AFT model for `X | Z` specified
  by `gamma_x` and `model_xz`.

This function takes the estimated parameter vector
`theta = c(beta, psi)` from `estimate_beta_aipw_est` and computes a
sandwich variance for \\\beta\\, treating the nuisance model parameters
as plug-in.
