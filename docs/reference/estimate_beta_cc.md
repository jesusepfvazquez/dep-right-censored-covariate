# Complete-case estimator for regression parameters

Fits a regression model using only complete cases, defined by `D = 1`,
via M-estimation using the `geex` framework. This is the complete-case
analogue of
[`estimate_beta`](https://jesusepfvazquez.github.io/dep-right-censored-covariate/reference/estimate_beta.md),
and is primarily intended for comparison in simulation studies. Not
consistent under outcome dependent censoring.

## Usage

``` r
estimate_beta_cc(data_yXZ, model)
```

## Arguments

- data_yXZ:

  A data frame containing at least the variables in `model` and a binary
  indicator `D`, where `D = 1` denotes a complete case.

- model:

  An object of class [`formula`](https://rdrr.io/r/stats/formula.html)
  specifying the regression model for the outcome.

## Value

A list with components

- beta_est:

  A 1 x (p + 1) matrix of regression coefficients followed by the
  residual standard deviation (last element).

- se_est:

  A 1 x p matrix of estimated standard errors for the regression
  coefficients.
