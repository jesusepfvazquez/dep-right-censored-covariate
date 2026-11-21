## helper_fit_all.R
## Fits every estimator on one dataset and returns one row per estimator:
##   oracle          complete data, y ~ (A - X) + Z
##   naive           y ~ (A - W) + Z, ignoring censoring
##   cc              complete case
##   ipw-est         IPW; SE accounts for estimating the weights (var_beta_ipw)
##   mle-joint       MLE, theta and f_{X|Z} estimated jointly; sandwich SE (var_beta_mle)
##   aipw-est-joint  unguarded AIPW, theta = (beta, sigma), augmentation from the joint
##                   MLE's f_{X|Z}; corrected sandwich (helper_unguarded_sigma.R)
##   aipw-lambda-est guarded AIPW, theta = (beta, sigma), working normal from the joint
##                   MLE's f_{X|Z} moments; corrected sandwich (helper_guarded_sigma.R)
## do_var_ipw / do_var_mle switch off the slow IPW and MLE sandwiches (the AIPW SEs are cheap).

source("helper_mestimate.R")
source("helper_joint_mle.R")
source("helper_guarded_sigma.R")
source("helper_unguarded_sigma.R")

WM <- "mle"   # working normal for the guarded AIPW: "mle" (joint MLE moments) or "obs" (observed W)

fit_all <- function(dat, do_var_ipw = TRUE, do_var_mle = TRUE, sigma_name = "sigma", tag = "") {
  safe <- function(expr, what) tryCatch(expr, error = function(e) {
    message("FAILED ", tag, " ", what, ": ", conditionMessage(e)); NULL })
  na_or <- function(x, k) if (is.null(x)) rep(NA, k) else as.numeric(x)[1:k]
  row <- function(method, est, se, g = NULL, mle = NULL, wm = NA) {
    r <- data.frame(Method = method,
      beta0 = na_or(est, 4)[1], beta1 = na_or(est, 4)[2], beta2 = na_or(est, 4)[3], sigma = na_or(est, 4)[4],
      sebeta0 = na_or(se, 3)[1], sebeta1 = na_or(se, 3)[2], sebeta2 = na_or(se, 3)[3],
      gamma0 = na_or(g, 3)[1], gammaZ = na_or(g, 3)[2], scale_x = na_or(g, 3)[3],
      convcode    = if (is.null(mle)) NA else mle$convcode,
      guard       = if (is.null(mle)) NA else as.integer(mle$restarted) + as.integer(mle$on_edge),
      start_fixed = if (is.null(mle)) NA else as.integer(mle$start_fixed),
      wm = wm)
    names(r)[names(r) == "sigma"] <- sigma_name
    r
  }
  ## guard: 0 = unconstrained fit kept, 1 = refit inside the box, 2 = refit on the box edge
  ## start_fixed: 1 = survreg start for gamma diverged and was refit

  est00 <- safe(estimate_beta(dat, y ~ AX + Z), "oracle")
  est01 <- safe(estimate_beta(dat, y ~ AW + Z), "naive")
  est02 <- safe(estimate_beta_cc(dat, y ~ AW + Z), "cc")

  est10 <- safe(estimate_beta_ipw(dat, y ~ AW + Z, myweight = "aft"), "ipw")
  est11 <- if (do_var_ipw && !is.null(est10)) safe(var_beta_ipw(dat, est10$beta_est), "ipw var") else NULL
  ipw_se <- if (do_var_ipw) est11$se_est else est10$se_est     # geex SE (weights treated as known) if no var

  mle   <- safe(estimate_beta_likelihood_joint(dat, y ~ AW + Z), "mle")
  mlev  <- if (do_var_mle && !is.null(mle)) safe(var_beta_mle(dat, mle$theta_full), "mle var") else NULL
  gam   <- if (!is.null(mle)) mle$gamma_est else NULL

  ung   <- if (!is.null(gam)) safe(estimate_beta_aipw_unguarded_sigma(dat, gam), "unguarded aipw") else NULL
  wm    <- if (WM == "mle" && is.null(gam)) "obs" else WM
  grd   <- safe(estimate_beta_aipw_lambda_sigma(dat, wm = wm, gamma = gam), "guarded aipw")

  rbind(
    row("oracle",          est00$beta_est, est00$se_est),
    row("naive",           est01$beta_est, est01$se_est),
    row("cc",              est02$beta_est, est02$se_est),
    row("ipw-est",         est10$beta_est, ipw_se),
    row("mle-joint",       mle$beta_est,   mlev$se_est, gam, mle),
    row("aipw-est-joint",  ung$beta_est,   ung$se_est,  gam),
    row("aipw-lambda-est", grd$beta_est,   grd$se_est,  if (wm == "mle") gam, wm = wm))
}
