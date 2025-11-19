# AIPW lambda-close estimator

Fits a regression model for `y` using an augmented IPW estimator with a
closed-form ("lambda-close") augmentation under a normal approximation
for `X`. This implementation assumes the outcome model is `y ~ AW + Z`,
where `AW = A - X`, and uses a single covariate `Z` in the outcome
model.

## Usage

``` r
estimate_beta_aipw_lambda(
  data_yXZ = dat,
  model = model,
  model_weights = model_weights,
  aw_var = "AW"
)
```

## Arguments

- data_yXZ:

  Data frame containing at least `y, A, AW, W, D, Z`, and the weight
  variables `myp_ywz_oracle`, `myp_ywz`, `myp_uniform`, `myp_ywz_logit`.
  Assumes the outcome model is `y ~ AW + Z`.

- model:

  Formula for the outcome, currently assumed to be `y ~ AW + Z`. The
  function checks that the design matrix has exactly three columns:
  intercept, AW, and Z, in that order.

- model_weights:

  A [`formula`](https://rdrr.io/r/stats/formula.html) specifying the
  censoring model for `C | (Y, Z, ...)`. Typically a right-hand-side
  only formula, such as `~ y + Z`, which is internally expanded to
  `Surv(W, 1 - D) ~ y + Z`.

- aw_var:

  Character string giving the name of the exposure covariate in `model`
  that equals `A - X` (default `"AW"`). This must appear in `model` and
  in `data_yXZ`.

## Value

A list with components:

- beta_est:

  1 x 4 matrix: \\(\hat\beta_0, \hat\beta_1, \hat\beta_2, \hat\psi)\\.

- se_est:

  1 x 3 matrix of standard errors for \\\beta_0, \beta_1, \beta_2\\
  based on the geex sandwich estimator.

## Details

The estimator combines:

- an IPW component using weights `D / p`, where `p` is chosen via
  `myweight`, and

- a closed-form augmentation term derived under a normal approximation
  for `X`.
