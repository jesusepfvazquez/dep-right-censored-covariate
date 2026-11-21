## Unguarded AIPW estimator with theta = (beta, sigma) estimated together.
##
## Estimating equation:
##   sum_i  D_i/pi_i S(theta) + (1 - D_i/pi_i) Psi(theta) = 0,
## with S the full-data score for (beta0, beta1, beta2, sigma) and, as in
## estimate_beta_aipw_est_joint,
##   Psi_j = int S_j(t) w(t) dt / int w(t) dt,
##   w(t)  = f(y | t, z; theta) f_{X|Z}(t; gamma) {1 - 1/pi(t)},
## where f_{X|Z} is the joint MLE's Weibull (gamma = mle$gamma_est) and
## pi(t) = P(C > t | Y, Z) from the Weibull AFT fit of C | (Y, Z).
## The integrals use a fixed Simpson grid on (0, upper) with 2m + 1 nodes (m = 200 matches
## m = 1000 to 6 decimals in the estimates and SEs), so the estimating function is
## smooth in theta and alpha (needed for the Newton steps and the sandwich).
##
## Sandwich: A^{-1} B A^{-T}, with B built from G_i + (dGbar/dalpha) (alpha-hat - alpha)_i,
## so the estimation of alpha is propagated. gamma is held fixed: E{(1 - D/pi) Psi} = 0
## for any gamma when pi is correct, so its estimation does not change the variance.
## include_sigma = FALSE fixes sigma at the IPW estimate and solves only for beta.

library(survival)

.pi_fun <- function(alpha, W, y, Z) {
  lp <- alpha[1] + alpha[2] * y + alpha[3] * Z
  exp(-exp((log(W) - lp) / alpha[4]))
}
.alpha_ll <- function(alpha, W, y, Z, D) {
  lp <- alpha[1] + alpha[2] * y + alpha[3] * Z
  u  <- (log(W) - lp) / alpha[4]
  S  <- exp(-exp(u)); f <- (1 / alpha[4]) * exp(u) * S
  log((1 - D) * f + D * S)
}

.simpson <- function(upper, m = 1000) {           # 2m + 1 nodes on (0, upper)
  t <- seq(0, upper, length.out = 2 * m + 1); h <- t[2] - t[1]
  w <- rep(c(2, 4), length.out = 2 * m + 1); w[1] <- 1; w[2 * m + 1] <- 1
  list(t = t[-1], w = (h / 3 * w)[-1])             # drop t = 0 (log(0) in pi)
}

.G_unguarded <- function(th, alpha, dat, gamma, grid, include_sigma) {
  b <- th[1:3]; s <- th[4]
  y <- dat$y; A <- dat$A; W <- dat$W; Z <- dat$Z; D <- dat$D
  pi_ <- .pi_fun(alpha, W, y, Z)
  ## full-data score at x = W
  e  <- y - (b[1] + b[2] * (A - W) + b[3] * Z)
  S  <- cbind(e, e * (A - W), e * Z) / s^2
  if (include_sigma) S <- cbind(S, -1 / s + e^2 / s^3)
  S[D == 0, ] <- 0
  ## augmentation: integrals over the grid (n x m matrices)
  tt  <- matrix(grid$t, nrow(dat), length(grid$t), byrow = TRUE)
  et  <- y - (b[1] + b[2] * (A - tt) + b[3] * Z)
  sc  <- exp(gamma[1] + gamma[2] * Z)
  wt  <- dnorm(et, 0, s) * dweibull(tt, shape = 1 / gamma[3], scale = sc) *
         (1 - 1 / .pi_fun(alpha, tt, y, Z))
  wt[!is.finite(wt)] <- 0
  wq  <- sweep(wt, 2, grid$w, `*`)
  den <- rowSums(wq)
  Ps  <- cbind(rowSums(wq * et), rowSums(wq * et * (A - tt)), Z * rowSums(wq * et)) / s^2
  if (include_sigma) Ps <- cbind(Ps, -1 / s * den + rowSums(wq * et^2) / s^3)
  Ps  <- Ps / den
  S * (D / pi_) + Ps * (1 - D / pi_)
}

estimate_beta_aipw_unguarded_sigma <- function(dat, gamma, include_sigma = TRUE, do_var = TRUE,
                                               upper = 20, m = 200, sigma_fix = NULL) {
  wr    <- survreg(Surv(W, 1 - D) ~ y + Z, data = dat, dist = "w")
  alpha <- c(coef(wr), wr$scale)
  grid  <- .simpson(upper, m)
  wts   <- dat$D / .pi_fun(alpha, dat$W, dat$y, dat$Z)
  fit   <- lm(y ~ AW + Z, data = dat, weights = wts)
  th0   <- c(coef(fit), sqrt(weighted.mean(resid(fit)^2, wts)))
  if (!is.null(sigma_fix)) th0[4] <- sigma_fix      # only for checks against the old code
  k     <- if (include_sigma) 4 else 3
  full  <- function(p) if (include_sigma) p else c(p, th0[4])
  Gbar  <- function(p) colSums(.G_unguarded(full(p), alpha, dat, gamma, grid, include_sigma))

  par <- th0[1:k]
  for (it in 1:50) {
    g <- Gbar(par); J <- numDeriv::jacobian(Gbar, par, method = "simple")
    step <- solve(J, g); lam <- 1
    repeat {
      pn <- par - lam * step
      ok <- !include_sigma || pn[4] > 0
      if (ok && sum(Gbar(pn)^2) <= sum(g^2) || lam < 1e-4) break
      lam <- lam / 2
    }
    conv <- max(abs(pn - par)) < 1e-9; par <- pn; if (conv) break
  }
  th  <- full(par)
  out <- list(beta_est = th, se_est = rep(NA, k), iter = it)
  if (!do_var) return(out)

  n   <- nrow(dat)
  Gi  <- .G_unguarded(th, alpha, dat, gamma, grid, include_sigma)
  sa  <- numDeriv::jacobian(function(a) .alpha_ll(a, dat$W, dat$y, dat$Z, dat$D), alpha)
  Ha  <- numDeriv::hessian(function(a) sum(.alpha_ll(a, dat$W, dat$y, dat$Z, dat$D)), alpha) / n
  Aai <- -solve(Ha)
  dGa <- numDeriv::jacobian(function(a) colMeans(.G_unguarded(th, a, dat, gamma, grid, include_sigma)), alpha)
  phi <- Gi + sa %*% t(Aai) %*% t(dGa)
  Am  <- numDeriv::jacobian(Gbar, par)
  Ai  <- solve(Am)
  out$se_est <- sqrt(diag(Ai %*% crossprod(phi) %*% t(Ai)))
  out
}
