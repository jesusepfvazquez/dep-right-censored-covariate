# ODRCC 0.2.1

Documentation only; the estimators are unchanged.

* The vignette and help pages follow the revised manuscript: "the MLE" replaces
  "the joint MLE", the help page of `ess_weights()` gives the manuscript's
  recommendation by effective sample size, and the help page of
  `estimate_beta_aipw()` states what the efficiency gain of the unguarded AIPW
  estimator requires.

# ODRCC 0.2.0

The estimators now match the revised manuscript and the code in `simulations/`.

* All estimators estimate theta = (beta, sigma) jointly and return sandwich
  standard errors directly in `se_est`. The `var_beta_*()` functions are removed.
* `estimate_beta_mle()` is the joint MLE: the Weibull model for X | Z is
  estimated jointly with the outcome model (64-node Gauss-Legendre integration).
* `estimate_beta_aipw()` (unguarded AIPW) takes the working model for
  f_{X|Y,Z} from the joint MLE (`fx_method = "joint"`, the default) or from a
  Weibull AFT model for X | (Y, Z) (`fx_method = "direct"`).
* `estimate_beta_aipw_lambda()` (guarded AIPW) uses a working normal model with
  the moments of the fitted Weibull model for X | Z.
* The standard errors of the IPW and AIPW estimators account for estimating the
  model for C | (Y, Z).
* New `ess_weights()`: effective sample size of the IPW weights.
* `data_aft()` now generates the design of the first simulation study, with
  arguments `c_int` (censoring rate) and `p.c` (spread of log C).
* Any number of covariates in the outcome model, `model_xz`, and `model_weights`.
* Dependencies reduced to survival and numDeriv.
