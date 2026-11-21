## helper_joint_mle.R
## MLE with theta = (beta, sigma) and gamma (f_{X|Z}) estimated jointly from the
## full likelihood of (Y, W, Delta, Z). Same likelihood as
## estimate_beta_likelihood_optimx() in helper_mestimate.R, with four changes:
##   1. sigma and the Weibull scale of X|Z are optimized on the log scale
##      (keeps both positive);
##   2. the returned sigma is exp(log sigma); the original returned
##      log(log sigma) in position 4;
##   3. the log-likelihood is evaluated on pre-built vectors instead of
##      filtering the data row by row (same values, much faster);
##   4. Nelder-Mead is followed by a BFGS polish, and the convergence code is
##      returned so non-converged fits can be flagged;
##   5. guard: a fit outside a wide box around the starting values (beta +/- 5;
##      log sigma, gamma, log scale_x +/- 2) is discarded and refit with nlminb
##      inside the box; `restarted` and `on_edge` flag these fits;
##   6. when survreg's default fit for the gamma starting values diverges (or stops with
##      an error), it is refit from least-squares starting values (`start_fixed` flags this).
## Datasets where neither guard triggers give the same fit as the earlier
## version, so only flagged seeds need to be rerun.
## Output parameterization matches var_beta_mle():
##   c(beta0, beta_{A-X}, beta_Z, sigma, gamma0, gammaZ, scale_x),
## where X | Z ~ Weibull(shape = 1/scale_x, scale = exp(gamma0 + gammaZ Z)).

library(optimx)

estimate_beta_likelihood_joint <- function(data_yXZ, model = y ~ AW + Z,
                                           upperbound = Inf, maxit = 5000) {

  ## starting values: weighted least squares for theta; for gamma, the inverse-probability
  ## weighted Weibull fit of W on Z among the uncensored rows (weights 1 / pi-hat). The
  ## weighted fit is consistent for f_{X|Z}. survreg(Surv(W, D) ~ Z) is not (C depends on
  ## Y, so C is not independent of X given Z), and with a wide censoring distribution
  ## (p.c = 3) it either did not converge or converged far from the truth, which pushed
  ## the bounded fit to the edge of its box (beta0 bias of -20 to -40%).
  uo     <- data_yXZ[data_yXZ$D == 1, , drop = FALSE]
  wr     <- tryCatch(survreg(Surv(W) ~ Z, data = uo, weights = 1 / uo$myp_ywz, dist = "w"),
                     error = function(e) NULL)
  sv_wr  <- if (is.null(wr)) NA else c(coef(wr), log(wr$scale))
  start_fixed <- !(all(is.finite(sv_wr)) && all(abs(sv_wr) < 5))
  if (start_fixed) {
    l0 <- lm(log(W) ~ Z, data = uo, weights = 1 / uo$myp_ywz)
    sv_wr <- c(coef(l0), log(sqrt(sum(resid(l0)^2 / uo$myp_ywz) / sum(1 / uo$myp_ywz)) * sqrt(6) / pi))
  }
  mylmer <- lm(model, data = data_yXZ, weights = D / myp_ywz)
  start  <- c(coef(mylmer), log(summary(mylmer)$sigma), sv_wr)

  Xm <- model.matrix(model, data = data_yXZ)          # 1, A - W, Z
  Y  <- model.response(model.frame(model, data = data_yXZ))
  A  <- data_yXZ$A
  z  <- data_yXZ$Z
  w  <- data_yXZ$W
  obs <- data_yXZ$D == 1

  loglik <- function(p) {
    beta <- p[1:3]; sig <- exp(p[4]); g <- p[5:6]; sx <- exp(p[7])
    meanXZ <- exp(g[1] + g[2] * z)
    ll <- numeric(length(Y))
    ## Delta = 1: f(y | x, z) f(x | z) at x = w
    ll[obs] <- dnorm(Y[obs] - Xm[obs, , drop = FALSE] %*% beta, 0, sig, log = TRUE) +
               dweibull(w[obs], shape = 1 / sx, scale = meanXZ[obs], log = TRUE)
    ## Delta = 0: integral over x > w of f(y | x, z) f(x | z)
    for (i in which(!obs)) {
      f <- function(x) dnorm(Y[i] - (beta[1] + beta[2] * (A[i] - x) + beta[3] * z[i]), 0, sig) *
                       dweibull(x, shape = 1 / sx, scale = meanXZ[i])
      v <- tryCatch(integrate(f, lower = w[i], upper = upperbound)$value,
                    error = function(e) 0)
      ll[i] <- log(max(v, 1e-300))
    }
    s <- sum(ll)
    if (!is.finite(s)) s <- -1e10
    s
  }

  ## plausible region for the fit: wide box around the starting values
  lower <- start - c(5, 5, 5, 2, 2, 2, 2)
  upper <- start + c(5, 5, 5, 2, 2, 2, 2)
  inbox <- function(p) all(is.finite(p)) && all(p >= lower & p <= upper)

  fits <- list()
  ## 1. Nelder-Mead from the starting values, then a BFGS polish
  r1 <- tryCatch(optimx::optimx(par = start, fn = loglik, method = "Nelder-Mead",
                                itnmax = maxit, control = list(maximize = TRUE)),
                 error = function(e) NULL)
  if (!is.null(r1)) {
    fits$nm <- r1
    r2 <- tryCatch(optimx::optimx(par = as.numeric(coef(r1)), fn = loglik, method = "BFGS",
                                  control = list(maximize = TRUE)),
                   error = function(e) NULL)
    if (!is.null(r2)) fits$bfgs <- r2
  }
  ## -1e10 is the value loglik() returns when the log-likelihood is not finite
  ok <- Filter(function(r) is.finite(r$value) && r$value > -1e9 &&
                 inbox(as.numeric(coef(r))), fits)

  ## 2. guard: if the unconstrained fit left the box (or failed), refit with bounds
  restarted <- length(ok) == 0
  if (restarted) {
    r3 <- tryCatch(optimx::optimx(par = start, fn = loglik, method = "nlminb",
                                  lower = lower, upper = upper,
                                  control = list(maximize = TRUE)),
                   error = function(e) NULL)
    if (!is.null(r3) && is.finite(r3$value) && r3$value > -1e9) ok$nlminb <- r3
  }
  if (length(ok) == 0) stop("joint MLE failed: no finite fit inside the plausible region")
  best <- ok[[which.max(sapply(ok, function(r) r$value))]]
  p <- as.numeric(coef(best))
  ## a bounded fit that ends on the box edge is not an interior maximum: flag it
  on_edge <- any(abs(p - lower) < 1e-4 | abs(p - upper) < 1e-4)

  est <- c(p[1:3], exp(p[4]), p[5:6], exp(p[7]))
  names(est) <- c("beta0", "beta1", "beta2", "sigma", "gamma0", "gammaZ", "scale_x")
  list(beta_est = est[1:4], gamma_est = est[5:7], theta_full = est,
       convcode = best$convcode, loglik = best$value,
       restarted = restarted, on_edge = on_edge, start_fixed = start_fixed)
}
