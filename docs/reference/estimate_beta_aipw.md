# Augmented IPW (AIPW) estimator for regression parameters

Fits a regression model for `Y` using an augmented inverse probability
weighting (AIPW) estimator in the setting of an outcome-dependent
right-censored covariate. The outcome model is specified by `model`,
typically of the form `y ~ AW + Z1 + ... + Zp`, where `AW = A - X`.

## Usage

``` r
estimate_beta_aipw(
  data_yXZ,
  model,
  model_weights,
  model_xz,
  aw_var = "AW",
  lbound = 0,
  ubound = 50
)
```

## Arguments

- data_yXZ:

  A data frame containing at least:

  - `y`: outcome,

  - `A`: auxiliary covariate used to form `AW = A - X`,

  - `W`: observed covariate `W = min(X, C)`,

  - `D`: indicator `I(X <= C)`,

  - columns for the covariates in `model`,

  - columns for the covariates in `model_weights` and `model_xz`.

- model:

  A [`formula`](https://rdrr.io/r/stats/formula.html) specifying the
  outcome regression model, e.g. `y ~ AW + Z` or `y ~ AW + Z1 + Z2`.

- model_weights:

  A [`formula`](https://rdrr.io/r/stats/formula.html) specifying the
  censoring model for `C | (Y, Z, ...)`. Typically a right-hand-side
  only formula, such as `~ y + Z`, which is internally expanded to
  `Surv(W, 1 - D) ~ y + Z`.

- model_xz:

  A [`formula`](https://rdrr.io/r/stats/formula.html) specifying the
  covariate structure for `X | Z`. Typically right-hand-side only, e.g.
  `~ Z`, meaning `log E[X | Z]` depends on those covariates. Only the
  RHS is used.

- aw_var:

  Character string giving the name of the exposure covariate in `model`
  that equals `A - X` (default `"AW"`). This must appear in `model` and
  in `data_yXZ`.

- lbound, ubound:

  Numeric lower and upper bounds for the numerical integration over `X`
  in the augmentation term (defaults: 0 and 50).

## Value

A list with components

- beta_est:

  A 1 x p matrix of estimated regression coefficients \\\hat \beta\\.

- psi_est:

  Scalar, estimated residual standard deviation \\\hat \psi\\.

## Details

The AIPW estimator combines:

- An IPW component based on a censoring model for `C | (Y, Z, ...)`
  specified by `model_weights`.

- An augmentation component that integrates over the conditional
  distribution `X | Z` specified by `gamma_x` and `model_xz`.

The function solves the AIPW estimating equations for the regression
parameters \\\beta\\, treating the nuisance models as plug-in, and
returns \\\hat \beta\\ along with a plug-in estimate of the residual
standard deviation \\\hat \psi\\.
