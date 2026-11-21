## Guarded AIPW estimator with theta = (beta, sigma) estimated together.
##
## Estimating equation (Lambda held fixed at its first-step value):
##   sum_i  D_i/pi_i S(theta) + (1 - D_i/pi_i) Lambda Psi(theta) = 0,
## with S the full-data score for (beta0, beta1, beta2, sigma) and Psi its closed form
## conditional expectation under the normal working model X | Z ~ N(meanXZ, sdXZ^2)
## (supplement S.1.3). meanXZ and sdXZ are the mean and SD of the observed W, as in
## the current Lambda and variance code.
##
## Lambda = -E(U V^T) E(V V^T)^{-1}, with U and V the IPW and augmentation pieces of the
## influence function, each including the correction for estimating alpha (pi).
## The same Lambda is used in the estimate and in the sandwich, so the SE belongs to
## the estimator that was computed.
##
## include_sigma = FALSE gives the beta-only version: sigma is fixed at the IPW
## estimate sqrt{ sum (D/pi) e^2 / sum (D/pi) } and only the three beta rows are solved.

library(survival)

.pi_fun <- function(alpha, W, y, Z) {
  lp <- alpha[1] + alpha[2] * y + alpha[3] * Z
  exp(-exp((log(W) - lp) / alpha[4]))
}
.alpha_ll <- function(alpha, W, y, Z, D) {        # per-row log-likelihood of C | Y, Z
  lp <- alpha[1] + alpha[2] * y + alpha[3] * Z
  u  <- (log(W) - lp) / alpha[4]
  S  <- exp(-exp(u)); f <- (1 / alpha[4]) * exp(u) * S
  log((1 - D) * f + D * S)
}

## full-data score at x = W (only used where D = 1)
.S_full <- function(th, y, A, W, Z, include_sigma) {
  b <- th[1:3]; s <- th[4]
  e <- y - (b[1] + b[2] * (A - W) + b[3] * Z)
  out <- cbind(e, e * (A - W), e * Z) / s^2
  if (include_sigma) out <- cbind(out, -1 / s + e^2 / s^3)
  out
}

## closed form E_{X | Y, Z}{ S(y, X, z; theta) } under X | Z ~ N(m, sd^2)
.Psi <- function(th, y, A, Z, m, sd, include_sigma) {
  b <- th[1:3]; s <- th[4]
  es  <- as.numeric(y - (b[1] + b[2] * A + b[3] * Z))           # e_*
  a   <- -1 / (2 * sd^2) - b[2]^2 / (2 * s^2)
  bb  <- m / sd^2 - b[2] * es / s^2
  mu  <- -bb / (2 * a)                                            # mu_*
  s2  <- (sd^2 * s^2) / (s^2 + b[2]^2 * sd^2)                     # sigma_*^2
  Ee   <- es + b[2] * mu                                          # E(e)
  EeAX <- -b[2] * (s2 + mu^2) + mu * (A * b[2] - es) + A * es     # E{e (A - X)}
  Ee2  <- es^2 + 2 * b[2] * es * mu + b[2]^2 * (s2 + mu^2)        # E(e^2)
  out <- cbind(Ee, EeAX, Z * Ee) / s^2
  if (include_sigma) out <- cbind(out, -1 / s + Ee2 / s^3)
  out
}

.pieces <- function(th, dat, alpha, m, sd, include_sigma) {
  y <- dat$y; A <- dat$A; W <- dat$W; Z <- dat$Z; D <- dat$D
  n   <- nrow(dat)
  pi_ <- .pi_fun(alpha, W, y, Z)
  S   <- .S_full(th, y, A, W, Z, include_sigma); S[D == 0, ] <- 0
  Ps  <- .Psi(th, y, A, Z, m, sd, include_sigma)
  pa  <- numDeriv::jacobian(function(a) .pi_fun(a, W, y, Z), alpha)          # n x 4, d pi / d alpha
  sa  <- numDeriv::jacobian(function(a) .alpha_ll(a, W, y, Z, D), alpha)     # n x 4, alpha scores
  Ha  <- numDeriv::hessian(function(a) sum(.alpha_ll(a, W, y, Z, D)), alpha) / n
  Aai <- -solve(Ha)                                                          # alpha-hat - alpha ~ Aai mean(sa)
  A_ipw <- crossprod(S * (-D / pi_^2), pa) / n          # E d{(D/pi) S}/d alpha
  A_at  <- crossprod(Ps * ( D / pi_^2), pa) / n          # E d{(1 - D/pi) Psi}/d alpha
  corr  <- sa %*% t(Aai)                                  # n x 4
  list(pi = pi_, S = S, Ps = Ps, D = D, A_ipw = A_ipw, A_at = A_at, corr = corr)
}

## Working normal for X | Z, used only inside Psi (efficiency, not consistency):
##   wm = "obs":     N(mean of observed W, SD of observed W), no Z (the current choice)
##   wm = "survreg": moments of a Weibull fit survreg(Surv(W, D) ~ Z) (two-step, cheap)
##   wm = "mle":     moments of the joint MLE's Weibull f_{X|Z}; pass gamma = mle$gamma_est
.working_moments <- function(dat, wm, gamma) {
  if (wm == "obs") return(list(m = mean(dat$W[dat$D == 1]), sd = sd(dat$W[dat$D == 1])))
  if (wm == "survreg") {
    fx <- survreg(Surv(W, D) ~ Z, data = dat, dist = "w")
    gamma <- c(coef(fx), fx$scale)
  }
  sc <- exp(gamma[1] + gamma[2] * dat$Z); k <- 1 / gamma[3]
  list(m = sc * gamma(1 + 1 / k), sd = sc * sqrt(gamma(1 + 2 / k) - gamma(1 + 1 / k)^2))
}

estimate_beta_aipw_lambda_sigma <- function(dat, include_sigma = TRUE, do_var = TRUE, wm = "obs", gamma = NULL) {
  wr    <- survreg(Surv(W, 1 - D) ~ y + Z, data = dat, dist = "w")
  alpha <- c(coef(wr), wr$scale)
  mm <- .working_moments(dat, wm, gamma); m <- mm$m; sd <- mm$sd
  pi_ <- .pi_fun(alpha, dat$W, dat$y, dat$Z)
  wts <- dat$D / pi_

  ## step 1: IPW start, with the IPW-consistent sigma
  fit <- lm(y ~ AW + Z, data = dat, weights = wts)
  th0 <- c(coef(fit), sqrt(weighted.mean(resid(fit)^2, wts)))
  k   <- if (include_sigma) 4 else 3

  ## step 2: Lambda at the start
  p0 <- .pieces(th0, dat, alpha, m, sd, include_sigma)
  U  <- p0$S * wts + p0$corr %*% t(p0$A_ipw)
  V  <- p0$Ps * (1 - wts) + p0$corr %*% t(p0$A_at)
  Lam <- -crossprod(U, V) %*% solve(crossprod(V, V))

  ## step 3: solve the estimating equation (Lambda fixed)
  G <- function(par) {
    th <- if (include_sigma) par else c(par, th0[4])
    S  <- .S_full(th, dat$y, dat$A, dat$W, dat$Z, include_sigma); S[dat$D == 0, ] <- 0
    Ps <- .Psi(th, dat$y, dat$A, dat$Z, m, sd, include_sigma)
    S * wts + (Ps * (1 - wts)) %*% t(Lam)
  }
  par <- th0[1:k]
  for (it in 1:50) {
    g <- colSums(G(par)); J <- numDeriv::jacobian(function(p) colSums(G(p)), par)
    step <- solve(J, g); lam <- 1
    repeat { pn <- par - lam * step; if (pn[k] > 0 || !include_sigma) if (sum(colSums(G(pn))^2) <= sum(g^2) || lam < 1e-4) break; lam <- lam / 2 }
    conv <- max(abs(pn - par)) < 1e-9; par <- pn; if (conv) break
  }
  th <- if (include_sigma) par else c(par, th0[4])
  out <- list(beta_est = th, se_est = rep(NA, k), Lambda = Lam, iter = it)
  if (!do_var) return(out)

  ## step 4: sandwich at theta-hat, alpha correction included, same Lambda
  p1  <- .pieces(th, dat, alpha, m, sd, include_sigma)
  phi <- G(par) + p1$corr %*% t(p1$A_ipw + Lam %*% p1$A_at)
  Am  <- numDeriv::jacobian(function(p) colSums(G(p)), par)
  Ai  <- solve(Am)
  V   <- Ai %*% crossprod(phi) %*% t(Ai)
  out$se_est <- sqrt(diag(V))
  out
}
