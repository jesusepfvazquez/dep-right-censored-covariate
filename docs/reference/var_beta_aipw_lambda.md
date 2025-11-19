# Sandwich variance for the AIPW-lambda estimator

Computes a sandwich variance estimator for the AIPW-lambda estimator
under outcome-dependent right-censoring of a covariate, using the
closed-form augmentation and accounting for estimation of the censoring
model \\C \| (Y, Z)\\ via a Gumbel/Weibull AFT model.

## Usage

``` r
var_beta_aipw_lambda(data_yXZ, mytheta)
```

## Arguments

- data_yXZ:

  Data frame containing at least:

  - `y`: outcome,

  - `A`: auxiliary covariate,

  - `AW`: `A - X` (used in the outcome model),

  - `W`: observed `W = min(X, C)`,

  - `D`: indicator `I(X <= C)`,

  - `Z`: covariate in the outcome and censoring models,

  - `myp_ywz`: weights \\\pi(Y, W, Z)\\ from AFT model (for IPW).

- mytheta:

  Numeric vector `c(beta0, beta1, beta2, psi)` corresponding to the
  AIPW-lambda point estimates from `estimate_beta_aipw_lambda_close()`.

## Value

A list with components:

- beta_est:

  The input parameter vector `mytheta`.

- se_est:

  Sandwich standard errors for \\\beta_0, \beta_1, \beta_2\\.

## Details

This implementation assumes the outcome model is `y ~ AW + Z`, i.e.,
three regression coefficients \\(\beta_0, \beta_1, \beta_2)\\ plus a
residual standard deviation \\\psi\\.
