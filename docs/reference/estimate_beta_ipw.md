# Inverse probability weighted (IPW) estimator for regression parameters

Fits a regression model using inverse probability weighting, where the
weights are derived from a parametric model for the censoring
distribution `C | (Y, Z, ...)`. The user supplies a right-hand-side
formula (e.g. `~ y + Z`) which is used to model `Surv(W, 1 - D)` via a
Weibull AFT model. The resulting estimated survival probabilities at `W`
are used as weights in the IPW estimating equations.

## Usage

``` r
estimate_beta_ipw(data_yXZ, model, model_weights)
```

## Arguments

- data_yXZ:

  A data frame containing the variables in `model`, as well as `W`
  (observed `min(X, C)`), `D` (indicator `I(X <= C)`), and the
  covariates appearing in `model_weights`.

- model:

  A [`formula`](https://rdrr.io/r/stats/formula.html) specifying the
  outcome regression model (e.g. `y ~ AW + Z`).

- model_weights:

  A right-hand-side formula specifying the variables in the censoring
  model, e.g. `~ y + Z`. This will be expanded to
  `Surv(W, 1 - D) ~ y + Z` internally.

## Value

A list with components

- beta_est:

  A 1 x (p + 1) matrix of regression coefficients followed by the
  residual standard deviation (last element).

- se_est:

  A 1 x p matrix of estimated standard errors for the regression
  coefficients.
