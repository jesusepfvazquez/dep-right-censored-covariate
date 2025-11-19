# Sandwich variance estimator for likelihood-based estimator

Computes a sandwich variance estimator for the parameter vector
\\\theta\\ obtained from `estimate_beta_likelihood_optimx`. The function
allows for a general outcome model `y ~ AW + Z1 + ... + Zp` and a
possibly different AFT model for `X | Z` specified by `model_xz`.

## Usage

``` r
var_beta_mle(data_yXZ, theta, model, model_xz = NULL, aw_var = "AW")
```

## Arguments

- data_yXZ:

  A data frame containing at least the outcome `y`, the auxiliary
  covariate `A`, the observed covariate `W`, the event indicator `D`,
  the exposure covariate `AW` (or another name specified by `aw_var`),
  and any additional covariates appearing in `model` and optionally
  `model_xz`.

- theta:

  Numeric vector of parameter estimates, typically taken from
  `estimate_beta_mle(...)[["beta_est"]]`. Must be ordered as \\(\beta,
  \psi, \gamma_x, \text{shape}\_x)\\.

- model:

  A [`formula`](https://rdrr.io/r/stats/formula.html) specifying the
  outcome regression model, e.g. `y ~ AW + Z1 + Z2`.

- model_xz:

  Optional [`formula`](https://rdrr.io/r/stats/formula.html) specifying
  the AFT model for `X | Z`. If `NULL` (default), the model for `X | Z`
  is taken to be Weibull with log-mean linear in `(1, Z1, ..., Zp)`,
  where `Z1, ..., Zp` are all covariates on the right-hand side of
  `model` except `aw_var`. If provided, `model_xz` can be either a
  right-hand-side formula (e.g. `~ Z1 + Z2`) or a full Surv formula
  (e.g. `Surv(W, D) ~ Z1 + Z2`); in both cases, only the RHS is used
  here to construct the design for `X | Z`.

- aw_var:

  Character string giving the name of the exposure covariate that is
  defined as `A - X` (default is `"AW"`). This variable must appear on
  the right-hand side of `model` and as a column in `data_yXZ`.

## Value

A list with components

- beta_est:

  Estimated regression coefficients \\\beta\\.

- psi_est:

  Estimated residual standard deviation \\\psi\\.

- se_beta:

  Sandwich standard errors for \\\beta\\.

- se_psi:

  Sandwich standard error for \\\psi\\.

- sandwich_var:

  Full sandwich variance matrix for \\\theta\\.

## Details

The parameter vector \\\theta\\ is assumed to be of the form \\\theta =
(\beta, \psi, \gamma_x, \text{shape}\_x)\\, where:

- \\\beta\\ are the regression coefficients in the outcome model,

- \\\psi\\ is the residual standard deviation on the original scale,

- \\\gamma_x\\ indexes the (log-)mean of \\X \| Z\\,

- \\\text{shape}\_x\\ is the Weibull shape parameter for \\X \| Z\\.
