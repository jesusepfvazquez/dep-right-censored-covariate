# Simulate AFT data with an outcome-dependent right-censored covariate

Generates a simulated dataset under an accelerated failure time (AFT)
model with a right-censored covariate subject to outcome-dependent
censoring. The function returns the true covariate `X`, censoring time
`C`, observed covariate `W = min(X, C)`, event indicator
`D = I(X <= C)`, and several derived quantities such as weights and
conditional expectations.

## Usage

``` r
data_aft(nSubjects = 10^3)
```

## Arguments

- nSubjects:

  Integer. Number of subjects to simulate. Defaults to `10^3`.

## Value

A data frame with `nSubjects` rows containing:

- y:

  Continuous outcome.

- Z:

  Standard normal covariate.

- D:

  Event indicator `I(X <= C)`.

- A:

  Observed auxiliary covariate.

- X:

  True covariate subject to censoring.

- AX:

  `A - X`.

- C:

  Censoring time.

- W:

  `min(X, C)`.

- AW:

  `A - W`.

- meanCYZ:

  Mean of `C` given `(Y, Z)` on the original scale.

- meanXZ:

  Mean of `X` given `(Z)` on the original scale.

- b:

  Subject ID.

- myp_uniform:

  Random weights (Uniform(0.1, 0.9)).

- myp_ywz_oracle:

  Oracle survival probability `P(C >= W | Y, Z)`.

- myp_ywz:

  Estimated survival probability using a Weibull AFT model.

- meanXD0YZ:

  \\E\[X \| D = 0, Y, Z\]\\ computed by numerical integration.

- AX_yz:

  `A - X` if `D = 1`, otherwise `A - E[X | D = 0, Y, Z]`.

- AX_z:

  `A - X` if `D = 1`, otherwise `A - E[X | Z]`.
