# ODRCC Tutorial

## 1. Introduction

The `ODRCC` package implements robust and efficient estimators for
regression models when a key covariate is **right-censored in an
outcome-dependent way**.

In many applications, we observe

- a continuous outcome $`Y`$
- fully observed covariates $`Z`$,
- an auxiliary variable $`A`$ (represents age),
- a covariate of interest $`X`$ that is **right-censored** by $`C`$,  
  so we only observe
  ``` math

  W = \min(X, C), \quad D = I(X \le C).
  ```

A common modeling strategy is to work with the transformed covariate
(indicating time to event)
``` math

AW = A - X,
```
and fit a regression of the form
``` math

Y = \beta_0 + \beta_1 AW + \beta_2 Z + \epsilon,
```
where $`\epsilon`$ is mean-zero noise.

The complete-case estimator discards subjects with $`D=0`$, leading to:

- loss of efficiency, and  
- potential **bias** when censoring depends on the outcome
  (outcome-dependent censoring).

`ODRCC` provides several estimators that remain consistent and can
improve efficiency:

- **IPW**: inverse probability weighting using a model for
  $`C \mid (Y, Z)`$
- **MLE**: likelihood-based estimator using a model for $`X \mid Z`$
- **AIPW**: augmented IPW estimator
- **AIPW-$`\Lambda`$**: a closed-form augmentation variant based on a
  robust matrix $`\Lambda`$ adjustment

This vignette walks through a typical analysis workflow with **simulated
data**.

------------------------------------------------------------------------

## 2. Load the package and simulate data

``` r
devtools::load_all(".")
#> ℹ Loading ODRCC
#> Warning: replacing previous import 'geex::coef' by 'stats::coef' when loading
#> 'ODRCC'
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> 
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> 
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(tictoc) # to time the functions
```

We will use the internal data-generating function
[`data_aft()`](https://jesusepfvazquez.github.io/dep-right-censored-covariate/reference/data_aft.md)
which simulates data from an AFT-type model with outcome-dependent
censoring of $`X`$.

``` r
set.seed(2025)

n <- 300
dat <- data_aft(nSubjects = n)

str(dat)
#> 'data.frame':    300 obs. of  18 variables:
#>  $ y             : num  -1.53 1.6 -3.16 4.02 1.65 ...
#>  $ Z             : num  0.6208 0.0356 0.7732 1.2725 0.371 ...
#>  $ D             : num  1 0 0 1 0 0 1 1 0 1 ...
#>  $ A             : num  1.1 1.66 2.1 0.97 1.03 ...
#>  $ X             : num  4.261 1.547 6.573 0.417 0.666 ...
#>  $ AX            : num  -3.156 0.116 -4.478 0.553 0.361 ...
#>  $ C             : num  5.143 0.13 0.152 4.78 0.532 ...
#>  $ W             : num  4.261 0.13 0.152 0.417 0.532 ...
#>  $ AW            : num  -3.156 1.532 1.943 0.553 0.495 ...
#>  $ meanCYZ       : num  7.96 1.25 19.41 0.69 1.44 ...
#>  $ meanXZ        : num  2.89 2.73 2.94 3.09 2.82 ...
#>  $ b             : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ myp_uniform   : num  0.612 0.139 0.489 0.238 0.8 ...
#>  $ myp_ywz_oracle: num [1:300, 1] 0.481 0.724 0.915 0.46 0.544 ...
#>  $ myp_ywz       : num  0.57 0.731 0.94 0.525 0.575 ...
#>  $ meanXD0YZ     : num  4.914 1.296 6.547 0.895 1.336 ...
#>  $ AX_yz         : num  -3.156 0.366 -4.452 0.553 -0.309 ...
#>  $ AX_z          : num  -3.156 -1.429 -1.232 0.553 -2.169 ...
summary(dat$D)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0000  0.0000  0.0000  0.4367  1.0000  1.0000
```

In this simulated dataset:

- `y` = outcome  
- `Z` = covariate  
- `A` = auxiliary variable (age)
- `X` = true covariate (subject to censoring)  
- `C` = censoring time  
- `W` = min(X, C)  
- `D` = I(X ≤ C)  
- `AX` = A - X; representing time to event (e.g., diagnosis)
- `AW` = A - W; coarsened version of AX  
- various weight variables (e.g., `myp_ywz`, `myp_ywz_oracle`, etc.)

For regression we’ll use the working model

``` math

Y \sim AW + Z
```

which in R is

``` r
model <- y ~ AW + Z
```

------------------------------------------------------------------------

## 3. Complete-case (CC) estimator

The **complete-case** estimator uses only subjects with $`D=1`$
(uncensored $`X`$), fitting a standard linear model via M-estimation:

``` r
tic()
est_cc <- estimate_beta_cc(dat, model = model)
est_cc
#> $beta_est
#>                                  psi_updated
#> [1,] 0.7097358 0.9449989 1.05793   0.9196435
#> 
#> $se_est
#>            [,1]       [,2]       [,3]
#> [1,] 0.08382001 0.02902404 0.08175951
toc()
#> 0.232 sec elapsed
```

[`estimate_beta_cc()`](https://jesusepfvazquez.github.io/dep-right-censored-covariate/reference/estimate_beta_cc.md)
returns:

- `beta_est`: a matrix with the regression coefficients and $`\hat\psi`$
- `se_est`: the corresponding standard errors for $`\beta`$

You can inspect them as:

``` r
cc_beta <- est_cc$beta_est
cc_se   <- est_cc$se_est

cc_beta
#>                                  psi_updated
#> [1,] 0.7097358 0.9449989 1.05793   0.9196435
cc_se
#>            [,1]       [,2]       [,3]
#> [1,] 0.08382001 0.02902404 0.08175951
```

------------------------------------------------------------------------

## 4. IPW estimator

The **IPW** estimator uses weights based on the estimated probability
that $`X`$ is observed, typically modeling

``` math

f_{C \mid Y, Z}(w,y,z)
```

with an AFT model. In the package version, the user supplies a formula
for the **weighting model**, e.g. `~ y + Z`, and the function internally
fits a Weibull AFT via `survreg()`.

Here we specify an IPW estimator with an AFT-based weighting model:

``` r
# outcome model
model <- y ~ AW + Z

# weights model: C | Y, Z
model_weights <- ~ y + Z

tic("IPW Estimator:")
est_ipw <- estimate_beta_ipw(
  data_yXZ      = dat,
  model         = model,
  model_weights = model_weights
)
est_ipw
#> $beta_est
#>                                  psi_updated
#> [1,] 0.896712 0.9327407 1.041258    0.940482
#> 
#> $se_est
#>            [,1]       [,2]       [,3]
#> [1,] 0.09242807 0.03129192 0.08593062
toc()
#> IPW Estimator:: 0.112 sec elapsed
```

Extract coefficients and standard errors:

``` r
ipw_beta <- est_ipw$beta_est
ipw_se   <- est_ipw$se_est

ipw_beta
#>                                  psi_updated
#> [1,] 0.896712 0.9327407 1.041258    0.940482
ipw_se
#>            [,1]       [,2]       [,3]
#> [1,] 0.09242807 0.03129192 0.08593062
```

If you want the robust sandwich SEs from the dedicated variance
function:

``` r
# Suppose est_ipw$beta_est is a 1 x p matrix including psi at the end
theta_ipw <- as.numeric(est_ipw$beta_est)

tic("IPW Estimator-Robust Sandwich Estimator:")
var_ipw <- var_beta_ipw(
  data_yXZ = dat,
  theta  = theta_ipw,
  model         = model,
  model_weights = model_weights
)
#> IPW sandwich: finished bread (A matrix).
#> IPW sandwich: finished meat (B matrix).
var_ipw$se_est
#> [1] 0.08722834 0.03090214 0.08429156
toc()
#> IPW Estimator-Robust Sandwich Estimator:: 57.557 sec elapsed
```

------------------------------------------------------------------------

## 5. Maximum Likelihood estimator (MLE)

The **MLE** approach uses a joint model, combining:

- the outcome model $`f_{Y \mid X, Z}`$, and  
- a model for $`f_{X \mid Z}`$ (Weibull AFT),

and integrates over censored $`X`$.

``` r
model_xz = as.formula("~Z")

tic("MLE:")
est_mle = estimate_beta_mle(
    data_yXZ = dat,
    model = model,
    aw_var   = "AW",
    model_weights = model_weights,
    model_xz = model_xz,
    trace = 0 # Set to 0 if no trace wanted, 1 otherwise
)
#> Maximizing -- use negfn and neggr
est_mle
#> $beta_est
#> [1] 0.9023008 0.9555097 0.9955888 0.9686213
#> 
#> $gamma_est
#> [1] 0.96285840 0.08246286 0.81313123
toc()
#> MLE:: 1221.958 sec elapsed
```

The corresponding sandwich variance:

``` r
theta_mle <- as.numeric(c(est_mle$beta_est,est_mle$gamma_est))

tic("MLE-Robust Sandwich Estimator:")
var_mle <- var_beta_mle(
    data_yXZ = dat,
    theta = theta_mle,
    model = model,
    aw_var   = "AW",
    model_xz = model_xz
)
var_mle$se_beta
#> [1] 0.07210658 0.02842251 0.07155912
toc()
#> MLE-Robust Sandwich Estimator:: 429.883 sec elapsed
```

------------------------------------------------------------------------

## 6. Augmented IPW (AIPW) estimator

The **AIPW** estimator combines IPW with an augmentation term providing
improved efficiency.

In the package version, you supply:

- `model`: outcome model, e.g. `y ~ AW + Z`  
- `model_weights`: censoring model for $`f_{C | (Y, Z)}`$
- `model_xz`: covariate structure for $`X | Z`$ (RHS-only formula,
  e.g. `~ Z`)

``` r
model_xz <- ~ Z

tic("AIPW:")
est_aipw <- estimate_beta_aipw(
  data_yXZ      = dat,
  model         = model,
  model_weights = model_weights,
  model_xz      = model_xz,
  aw_var        = "AW",
  lbound        = 0,
  ubound        = 50
)
est_aipw
#> $beta_est
#>           [,1]      [,2]    [,3]      [,4]
#> [1,] 0.8801625 0.9578444 1.12362 0.9381438
#> 
#> $beta_se
#>            [,1]      [,2]      [,3]
#> [1,] 0.09291546 0.0327527 0.1126476
toc()
#> AIPW:: 1611.002 sec elapsed
```

AIPW sandwich SEs:

``` r
# est_aipw$beta_est is 1 x p, est_aipw$psi_est is scalar
theta_aipw <- c(as.numeric(est_aipw$beta_est))
# theta_aipw <- c(as.numeric(est_ipw$beta_est))

tic("AIPW Estimator-Robust Sandwich Estimator:")
var_aipw <- var_beta_aipw(
  data_yXZ = dat,
  theta  = theta_aipw,
  lbound        = 0,
  ubound        = 50
)
var_aipw$se_est
#> [1] 0.07964022 0.02990928 0.07936218
toc()
#> AIPW Estimator-Robust Sandwich Estimator:: 63.627 sec elapsed
```

------------------------------------------------------------------------

## 7. AIPW–Lambda estimator (closed-form augmentation)

The **AIPW–Lambda** estimator uses an approximation to the augmentation
term that admits a closed-form expression under specific modeling
assumptions for $`f_{X | Z}`$. This provides a computationally cheaper
alternative to full numerical AIPW.

``` r
tic("AIPW with Lambda Estimator:")
est_aipw_lambda <- estimate_beta_aipw_lambda(
    data_yXZ = dat,
    model = model,
    model_weights = model_weights,
    aw_var  = "AW"
) 
est_aipw_lambda
#> $beta_est
#> [1] 0.8788797 0.9381677 1.0675614 0.9366018
#> 
#> $se_est
#>            [,1]      [,2]       [,3]
#> [1,] 0.08001977 0.0291584 0.07999267
toc()
#> AIPW with Lambda Estimator:: 0.379 sec elapsed
```

Variance via the dedicated robust sandwiich estimator function:

``` r
theta_aipw_lambda <- as.numeric(est_aipw_lambda$beta_est)

tic("AIPW with Lambda Estimator-Robust Sandwich Estimator:")
var_aipw_lambda <- var_beta_aipw_lambda(
  data_yXZ = dat,
  mytheta  = theta_aipw_lambda
)
#> AIPW with Lambda Estimator-Robust Sandwich Estimator:: 1.385 sec elapsed
var_aipw_lambda$se_est
#> [1] 0.08324157 0.02928765 0.08407404
toc()
```

------------------------------------------------------------------------

## 8. Comparing estimators

Once all estimators have been run, we can collect and compare them. Note
that the oracle value is 1 across all.

``` r
compare_beta <- rbind(
  CC      = est_cc$beta_est,
  IPW     = est_ipw$beta_est,
  MLE     = est_mle$beta_est,
  AIPW    = est_aipw$beta_est,
  AIPW_Lambda  = est_aipw_lambda$beta_est
)

row.names(compare_beta) = c("CC", "IPW", "MLE", "AIPW", "AIPW-Lambda")
colnames(compare_beta) = c("beta0", "beta1", "beta2", "sigma")
compare_beta
#>                 beta0     beta1     beta2     sigma
#> CC          0.7097358 0.9449989 1.0579300 0.9196435
#> IPW         0.8967120 0.9327407 1.0412579 0.9404820
#> MLE         0.9023008 0.9555097 0.9955888 0.9686213
#> AIPW        0.8801625 0.9578444 1.1236203 0.9381438
#> AIPW-Lambda 0.8788797 0.9381677 1.0675614 0.9366018
```

The corresponding nuisance-adjusted **standard errors** are:

``` r
se_mat <- rbind(
  CC      = est_cc$se_est,
  IPW     = var_ipw$se_est,
  MLE     = var_mle$se_beta,
  AIPW    = var_aipw$se_est,
  AIPW_Lambda  = var_aipw_lambda$se_est
)

row.names(se_mat) = c("CC", "IPW", "MLE", "AIPW", "AIPW-Lambda")
colnames(se_mat) = c("beta0", "beta1", "beta2")
se_mat
#>                  beta0      beta1      beta2
#> CC          0.08382001 0.02902404 0.08175951
#> IPW         0.08722834 0.03090214 0.08429156
#> MLE         0.07210658 0.02842251 0.07155912
#> AIPW        0.07964022 0.02990928 0.07936218
#> AIPW-Lambda 0.08324157 0.02928765 0.08407404
```
