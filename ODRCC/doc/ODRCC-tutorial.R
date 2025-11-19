## ----include = FALSE, eval = TRUE---------------------------------------------
# devtools::document()   # update documentation + NAMESPACE
# devtools::load_all()   # reload package into environment

## ----setup, include = FALSE---------------------------------------------------
rm(list = ls())
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4
)
set.seed(123)

## ----load-package-------------------------------------------------------------
devtools::load_all(".")
library(dplyr)
library(tictoc) # to time the functions

## ----simulate-data------------------------------------------------------------
set.seed(2025)

n <- 300
dat <- data_aft(nSubjects = n)

str(dat)
summary(dat$D)

## ----model-formula------------------------------------------------------------
model <- y ~ AW + Z

## ----cc-estimator-------------------------------------------------------------
tic()
est_cc <- estimate_beta_cc(dat, model = model)
est_cc
toc()

## ----cc-coefs-----------------------------------------------------------------
cc_beta <- est_cc$beta_est
cc_se   <- est_cc$se_est

cc_beta
cc_se

## ----ipw-estimator------------------------------------------------------------
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
toc()

## ----ipw-coefs----------------------------------------------------------------
ipw_beta <- est_ipw$beta_est
ipw_se   <- est_ipw$se_est

ipw_beta
ipw_se

## ----ipw-var, eval = TRUE-----------------------------------------------------
# Suppose est_ipw$beta_est is a 1 x p matrix including psi at the end
theta_ipw <- as.numeric(est_ipw$beta_est)

tic("IPW Estimator-Robust Sandwich Estimator:")
var_ipw <- var_beta_ipw(
  data_yXZ = dat,
  theta  = theta_ipw,
  model         = model,
  model_weights = model_weights
)
var_ipw$se_est
toc()

## ----mle-estimator, eval = TRUE-----------------------------------------------
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
est_mle
toc()

## ----mle-var, eval = TRUE-----------------------------------------------------
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
toc()

## ----aipw-gamma-setup, eval = TRUE--------------------------------------------
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
toc()

## ----aipw-var, eval = TRUE----------------------------------------------------
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
toc()

## ----lambda-estimator, eval = TRUE--------------------------------------------
tic("AIPW with Lambda Estimator:")
est_aipw_lambda <- estimate_beta_aipw_lambda(
    data_yXZ = dat,
    model = model,
    model_weights = model_weights,
    aw_var  = "AW"
) 
est_aipw_lambda
toc()

## ----lambda-var, eval = TRUE--------------------------------------------------
theta_aipw_lambda <- as.numeric(est_aipw_lambda$beta_est)

tic("AIPW with Lambda Estimator-Robust Sandwich Estimator:")
var_aipw_lambda <- var_beta_aipw_lambda(
  data_yXZ = dat,
  mytheta  = theta_aipw_lambda
)
var_aipw_lambda$se_est
toc()

## ----compare-estimators, eval = TRUE------------------------------------------
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

## ----compare-se, eval = TRUE--------------------------------------------------
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

