# Sandwich variance estimator for the IPW estimator

Computes a sandwich variance estimator for the IPW regression estimator,
allowing for a general outcome model `y ~ AW + Z1 + ... + Zp` and a
possibly different censoring model for `C | (Y, Z, ...)` specified by
`model_weights`. The censoring model is fit via a Weibull AFT model
using `survreg`, and the corresponding Gumbel parameterization for
`log(W)` is used to construct the weights \\\pi(Y, W, Z)\\.

## Usage

``` r
var_beta_ipw(data_yXZ, theta, model, model_weights)
```

## Arguments

- data_yXZ:

  A data frame containing at least the outcome `y`, the covariates in
  the outcome model `model`, the observed covariate `W`, the event
  indicator `D`, and the covariates appearing in the censoring model
  `model_weights`.

- theta:

  Numeric vector of parameter estimates \\(\beta, \psi)\\ from the IPW
  estimator. The length of `theta` must be equal to `p_beta + 1`, where
  `p_beta` is the number of regression coefficients in `model`.

- model:

  A [`formula`](https://rdrr.io/r/stats/formula.html) specifying the
  outcome regression model, e.g. `y ~ AW + Z1 + Z2`.

- model_weights:

  A [`formula`](https://rdrr.io/r/stats/formula.html) specifying the
  censoring model used to estimate the IPW weights. Typically of the
  form `~ y + Z1 + Z2` (right-hand-side only). Internally this is
  expanded to `Surv(W, 1 - D) ~ y + Z1 + Z2`. If a full `Surv` formula
  is provided, only its right-hand side is used.

## Value

A list with components

- beta_est:

  Estimated regression coefficients \\\beta\\.

- psi_est:

  Estimated residual standard deviation \\\psi\\.

- se_beta:

  Sandwich standard errors for \\\beta\\.

- sandwich_var:

  Full sandwich variance matrix for the stacked nuisance parameter
  vector \\\xi = (\beta, \alpha)\\, where \\\alpha\\ are the parameters
  of the censoring model.

## Details

The parameter vector `theta` is assumed to be of the form \\\theta =
(\beta, \psi)\\, where:

- \\\beta\\ are the regression coefficients in the outcome model,

- \\\psi\\ is the residual standard deviation on the original scale.
