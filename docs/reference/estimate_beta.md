# Estimate regression coefficients via M-estimation

Estimate regression coefficients via M-estimation

## Usage

``` r
estimate_beta(data_yXZ, model)
```

## Arguments

- data_yXZ:

  A data frame containing the variables in `model`. Must include the
  outcome `y` and covariates.

- model:

  An object of class [`formula`](https://rdrr.io/r/stats/formula.html)
  specifying the regression model for `y` in terms of covariates.

## Value

A list with components

- beta_est:

  A 1 x (p + 1) matrix of estimated regression coefficients and residual
  standard deviation (last element).

- se_est:

  A 1 x p matrix of estimated standard errors for the regression
  coefficients.
