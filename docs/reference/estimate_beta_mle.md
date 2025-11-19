# Likelihood-based estimator via optimx for general AW + Z model

Fits a parametric likelihood model for the outcome and the
right-censored covariate using direct maximization of the log-likelihood
with `optimx`. The outcome model is specified by `model`, typically of
the form `y ~ AW + Z1 + ... + Zp`, where `AW = A - X`. The distribution
of `X | Z` is modeled via a Weibull AFT model, whose covariate structure
can be specified by `model_xz` or, by default, derived from the
right-hand side of `model` by excluding `AW`.

## Usage

``` r
estimate_beta_mle(
  data_yXZ,
  model,
  aw_var = "AW",
  model_weights,
  model_xz = NULL,
  trace = 0
)
```

## Arguments

- data_yXZ:

  A data frame containing at least the outcome `y`, the auxiliary
  covariate `A`, the observed covariate `W`, the event indicator `D`,
  the exposure covariate `AW` (or another name specified by `aw_var`),
  and any additional covariates appearing in `model` and optionally
  `model_xz`.

- model:

  A [`formula`](https://rdrr.io/r/stats/formula.html) specifying the
  outcome regression model, e.g. `y ~ AW + Z1 + Z2`.

- aw_var:

  Character string giving the name of the exposure covariate that is
  defined as `A - X` (default is `"AW"`). This variable must appear on
  the right-hand side of `model` and as a column in `data_yXZ`.

- model_weights:

  A right-hand-side formula specifying the variables in the censoring
  model, e.g. `~ y + Z`. This will be expanded to
  `Surv(W, 1 - D) ~ y + Z` internally. This is needed so that we can
  compute an initial “good" estimate for our regression parameters using
  the IPW estimator.

- model_xz:

  Optional [`formula`](https://rdrr.io/r/stats/formula.html) specifying
  the AFT model for `X | Z`. If `NULL` (default), the model for `X | Z`
  is taken to be Weibull with log-mean linear in `(1, Z1, ..., Zp)`,
  where `Z1, ..., Zp` are all covariates on the right-hand side of
  `model` except `aw_var`. If provided, `model_xz` can be either a full
  Surv formula such as `Surv(W, D) ~ Z1 + Z2` or a right-hand-side
  formula such as `~ Z1 + Z2`, in which case the left-hand side
  `Surv(W, D)` is filled in automatically. ' @param trace Optional, set
  equal to 1 if you want the trace of the optimx model. Leave at 0 if
  you do not want a trace

## Value

A list with component

- beta_est:

  A numeric vector containing the parameter vector \\\theta = (\beta,
  \psi, \gamma_x, \text{shape}\_x)\\, where \\\beta\\ are the outcome
  regression coefficients (dimension determined by `model`), \\\psi\\ is
  the residual standard deviation on the original scale, \\\gamma_x\\
  indexes the mean of \\X \| Z\\, and \\\text{shape}\_x\\ is the Weibull
  shape parameter for \\X \| Z\\.
