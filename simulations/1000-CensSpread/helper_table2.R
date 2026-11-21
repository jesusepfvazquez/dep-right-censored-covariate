## helper_table2.R
## Estimators for the Table 2 rows that helper_known.R does not cover. All use
## theta = (beta0, beta_{A-X}, beta_Z, sigma).
##
## Specifications (as in the original Table 2):
##   correct pi     pi from the true alpha = (1, -0.5, 0.5, 2)
##   estimated pi   pi from survreg(Surv(W, 1 - D) ~ y + Z), SE corrected for alpha-hat
##   incorrect pi   D / myp_uniform, with myp_uniform ~ Uniform(0.1, 0.9) (weights treated as fixed)
##   correct Psi    Psi_eff: true pi(t) and true f_{X|Z} inside the integrals
##   incorrect Psi  Psi_closed: working normal for X | Z (supplement S.1.3), with the
##                  moments of the true Weibull f_{X|Z}, i.e. the same Psi the guarded
##                  AIPW estimator uses, so the two estimators differ only in Lambda
##                  (the unguarded estimator uses Lambda = I). The rows ending in "-n01"
##                  use the N(0, 1) working model of the original Table 2 run instead.
##   MLE            f_{X|Z} correct (gamma fixed at the truth), incorrect (scale_x = 1.2
##                  instead of 0.8), two-step (gamma from survreg(Surv(W, D) ~ Z), which
##                  ignores that C depends on Y), or joint (helper_joint_mle.R)
##
## MLE: the joint MLE is estimate_beta_likelihood_joint() (integrate(), as in the other
## studies). The fixed-gamma MLEs and every MLE sandwich use .ll_rows_gl(), the same
## likelihood with a Gauss-Legendre rule that is smooth in the parameters (stable
## numerical derivatives) and about 10 times faster; it agrees with integrate() to 1e-5
## per row. Sandwich: A^{-1} B A^{-T}, from per-row scores and the Hessian of the summed
## log-likelihood. For the two-step MLE the survreg score for gamma is stacked under the
## score for theta, so the SE accounts for estimating gamma.

source("helper_known.R")       # .newton, estimate_unguarded_known, estimate_guarded_known, truth
source("helper_joint_mle.R")   # estimate_beta_likelihood_joint

gamma_wrong <- c(1, 0.1, 1.2)  # incorrect f_{X|Z}: Weibull scale parameter 1.2 instead of 0.8

## ------------------------------------------------------------------------------------
## Oracle (X observed) and IPW with fixed weights: least squares with the empirical
## sandwich (the same estimates and SEs as estimate_beta() and estimate_beta_ipw() in the
## ODRCC code, without geex)
## ------------------------------------------------------------------------------------
.wls_sandwich <- function(dat, formula, wts) {
  Xm <- model.matrix(formula, data = dat)
  y  <- model.response(model.frame(formula, data = dat))
  b  <- as.numeric(solve(crossprod(Xm * wts, Xm), crossprod(Xm * wts, y)))
  e  <- as.numeric(y - Xm %*% b)
  Ai <- solve(crossprod(Xm * wts, Xm))
  V  <- Ai %*% crossprod(Xm * (wts * e)) %*% Ai
  list(b = b, e = e, se = sqrt(diag(V)))
}
estimate_oracle <- function(dat) {
  f <- .wls_sandwich(dat, y ~ AX + Z, rep(1, nrow(dat)))
  list(beta_est = c(f$b, sqrt(mean(f$e^2))), se_est = f$se)
}
## myp = NULL uses the true pi; otherwise the weights are D / myp (e.g. myp_uniform)
estimate_ipw_fixed <- function(dat, alpha = alpha_true, myp = NULL) {
  pi_ <- if (is.null(myp)) .pi_fun(alpha, dat$W, dat$y, dat$Z) else myp
  f   <- .wls_sandwich(dat, y ~ AW + Z, dat$D / pi_)
  list(beta_est = c(f$b, sqrt(mean(f$e[dat$D == 1]^2))), se_est = f$se)
}

## ------------------------------------------------------------------------------------
## IPW with estimated pi: fast version of var_beta_ipw() (same stacked estimating
## equations, vectorized; checked against var_beta_ipw)
## ------------------------------------------------------------------------------------
estimate_ipw_est <- function(dat) {
  wr    <- survreg(Surv(W, 1 - D) ~ y + Z, data = dat, dist = "w")
  alpha <- c(coef(wr), wr$scale)
  n     <- nrow(dat)
  pi_   <- .pi_fun(alpha, dat$W, dat$y, dat$Z)
  fit   <- lm(y ~ AW + Z, data = dat, weights = dat$D / pi_)
  b     <- coef(fit)
  Xm    <- model.matrix(~ AW + Z, data = dat)
  e     <- as.numeric(dat$y - Xm %*% b)
  U     <- Xm * (dat$D / pi_ * e)                                              # n x 3
  pa    <- numDeriv::jacobian(function(a) .pi_fun(a, dat$W, dat$y, dat$Z), alpha)
  sa    <- numDeriv::jacobian(function(a) .alpha_ll(a, dat$W, dat$y, dat$Z, dat$D), alpha)
  Ha    <- numDeriv::hessian(function(a) sum(.alpha_ll(a, dat$W, dat$y, dat$Z, dat$D)), alpha) / n
  A_ba  <- crossprod(U * (-1 / pi_), pa) / n                                   # d mean(U) / d alpha
  phi   <- U + sa %*% t(-solve(Ha)) %*% t(A_ba)
  Abb   <- -crossprod(Xm * (dat$D / pi_), Xm) / n
  Ai    <- solve(Abb)
  V     <- Ai %*% crossprod(phi / n) %*% t(Ai)
  sig   <- sqrt(mean((e[dat$D == 1])^2))                                       # as in estimate_beta_ipw
  list(beta_est = c(b, sig), se_est = sqrt(diag(V)))
}

## ------------------------------------------------------------------------------------
## Unguarded AIPW with any combination of weights and Psi; alpha, gamma held at the
## values passed (no correction in the sandwich)
##   wts = NULL uses D / pi(alpha); otherwise D / myp for the given myp
##   psi = "eff" (integrals; pi(t) from alpha, f_{X|Z} from gamma) or "closed"
##   working (psi = "closed" only): "gamma" = normal with the moments of the Weibull
##     f_{X|Z}(gamma), the same Psi_closed as the guarded estimator; "n01" = N(0, 1), the
##     working model of the original Table 2 run (estimate_beta_aipw(eta1 = TRUE))
## ------------------------------------------------------------------------------------
.G_closed <- function(th, dat, wts, m, sd, include_sigma) {
  S  <- .S_full(th, dat$y, dat$A, dat$W, dat$Z, include_sigma); S[dat$D == 0, ] <- 0
  Ps <- .Psi(th, dat$y, dat$A, dat$Z, m, sd, include_sigma)
  S * wts + Ps * (1 - wts)
}
.G_eff_wts <- function(th, alpha, dat, gamma, grid, include_sigma, wts) {
  ## same as .G_unguarded() (helper_unguarded_sigma.R) with the weights D/pi replaced by wts;
  ## pi(t) inside Psi still comes from alpha
  b <- th[1:3]; s <- th[4]
  y <- dat$y; A <- dat$A; W <- dat$W; Z <- dat$Z; D <- dat$D
  e  <- y - (b[1] + b[2] * (A - W) + b[3] * Z)
  S  <- cbind(e, e * (A - W), e * Z) / s^2
  if (include_sigma) S <- cbind(S, -1 / s + e^2 / s^3)
  S[D == 0, ] <- 0
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
  S * wts + Ps * (1 - wts)
}

estimate_unguarded_spec <- function(dat, alpha = alpha_true, gamma = gamma_true, myp = NULL,
                                    psi = c("eff", "closed"), working = c("gamma", "n01"),
                                    include_sigma = TRUE,
                                    do_var = TRUE, upper = 40, m = 400) {
  psi  <- match.arg(psi); working <- match.arg(working)
  wts  <- if (is.null(myp)) dat$D / .pi_fun(alpha, dat$W, dat$y, dat$Z) else dat$D / myp
  fit  <- lm(y ~ AW + Z, data = dat, weights = wts)
  th0  <- c(coef(fit), sqrt(weighted.mean(resid(fit)^2, wts)))
  k    <- if (include_sigma) 4 else 3
  full <- function(p) if (include_sigma) p else c(p, th0[4])
  if (psi == "eff") {
    grid <- .simpson(upper, m)
    Gi   <- function(p) .G_eff_wts(full(p), alpha, dat, gamma, grid, include_sigma, wts)
  } else {
    mm <- if (working == "gamma") .working_moments(dat, "mle", gamma) else list(m = 0, sd = 1)
    Gi <- function(p) .G_closed(full(p), dat, wts, mm$m, mm$sd, include_sigma)
  }
  Gbar <- function(p) colSums(Gi(p))
  nw   <- .newton(Gbar, th0[1:k], include_sigma)
  out  <- list(beta_est = full(nw$par), se_est = rep(NA, k), iter = nw$iter)
  if (!do_var) return(out)
  Ai <- solve(numDeriv::jacobian(Gbar, nw$par))
  out$se_est <- sqrt(diag(Ai %*% crossprod(Gi(nw$par)) %*% t(Ai)))
  out
}

## Guarded AIPW with any weights (true alpha or a fixed myp), working normal from gamma.
## Same as estimate_guarded_known() in helper_known.R, with the weights as an argument.
## type = "acc" and "macc" give the (inconsistent) ACC and MACC estimators of Table 2,
## now with theta = (beta, sigma); the original code solved for beta only, with sigma
## fixed at a plug-in value. All three solve
##   sum_i U_i(theta) + Lambda V_i(theta) = 0,   Lambda = -E(U V^T) E(V V^T)^{-1} at the IPW start,
## with U = (D/pi) S, V = (1 - D/pi) Psi    guarded
##      U = D S,      V = (1 - D/pi) Psi    MACC
##      U = D S,      V = (D - pi) Psi      ACC
## and Psi = Psi_closed under the working normal with the moments of f_{X|Z}(gamma).
estimate_guarded_spec <- function(dat, alpha = alpha_true, gamma = gamma_true, myp = NULL,
                                  type = c("guarded", "macc", "acc"),
                                  include_sigma = TRUE, do_var = TRUE) {
  type <- match.arg(type)
  mm  <- .working_moments(dat, "mle", gamma); m <- mm$m; sd <- mm$sd
  pi_ <- if (is.null(myp)) .pi_fun(alpha, dat$W, dat$y, dat$Z) else myp
  wts <- dat$D / pi_
  fit <- lm(y ~ AW + Z, data = dat, weights = wts)
  th0 <- c(coef(fit), sqrt(weighted.mean(resid(fit)^2, wts)))
  k   <- if (include_sigma) 4 else 3
  a_U <- switch(type, guarded = wts, macc = dat$D, acc = dat$D)
  a_V <- switch(type, guarded = 1 - wts, macc = 1 - wts, acc = dat$D - pi_)
  pieces <- function(th) {
    S  <- .S_full(th, dat$y, dat$A, dat$W, dat$Z, include_sigma); S[dat$D == 0, ] <- 0
    Ps <- .Psi(th, dat$y, dat$A, dat$Z, m, sd, include_sigma)
    list(U = S * a_U, V = Ps * a_V)
  }
  p0  <- pieces(th0)
  Lam <- -crossprod(p0$U, p0$V) %*% solve(crossprod(p0$V, p0$V))
  G <- function(par) {
    th <- if (include_sigma) par else c(par, th0[4])
    p  <- pieces(th); p$U + p$V %*% t(Lam)
  }
  nw  <- .newton(function(p) colSums(G(p)), th0[1:k], include_sigma)
  th  <- if (include_sigma) nw$par else c(nw$par, th0[4])
  out <- list(beta_est = th, se_est = rep(NA, k), Lambda = Lam, iter = nw$iter)
  if (!do_var) return(out)
  Ai <- solve(numDeriv::jacobian(function(p) colSums(G(p)), nw$par))
  out$se_est <- sqrt(diag(Ai %*% crossprod(G(nw$par)) %*% t(Ai)))
  out
}

## ------------------------------------------------------------------------------------
## MLE
## p = (beta0, beta1, beta2, log sigma, gamma0, gammaZ, log scale_x), as in
## estimate_beta_likelihood_joint()
## ------------------------------------------------------------------------------------
.mle_data <- function(dat, model = y ~ AW + Z) {
  list(Xm = model.matrix(model, data = dat), Y = model.response(model.frame(model, data = dat)),
       A = dat$A, z = dat$Z, w = dat$W, obs = dat$D == 1)
}

## per-row log-likelihood with integrate() (the code of estimate_beta_likelihood_joint)
.ll_rows_int <- function(p, s, upperbound = Inf) {
  beta <- p[1:3]; sig <- exp(p[4]); g <- p[5:6]; sx <- exp(p[7])
  meanXZ <- exp(g[1] + g[2] * s$z)
  ll <- numeric(length(s$Y)); o <- s$obs
  ll[o] <- dnorm(s$Y[o] - s$Xm[o, , drop = FALSE] %*% beta, 0, sig, log = TRUE) +
           dweibull(s$w[o], shape = 1 / sx, scale = meanXZ[o], log = TRUE)
  for (i in which(!o)) {
    f <- function(x) dnorm(s$Y[i] - (beta[1] + beta[2] * (s$A[i] - x) + beta[3] * s$z[i]), 0, sig) *
                     dweibull(x, shape = 1 / sx, scale = meanXZ[i])
    v <- tryCatch(integrate(f, lower = s$w[i], upper = upperbound)$value, error = function(e) 0)
    ll[i] <- log(max(v, 1e-300))
  }
  ll
}

## Gauss-Legendre nodes and weights on (-1, 1) (Golub-Welsch)
.gauss_legendre <- function(k) {
  j  <- seq_len(k - 1); b <- j / sqrt(4 * j^2 - 1)
  J  <- matrix(0, k, k); J[cbind(j, j + 1)] <- b; J[cbind(j + 1, j)] <- b
  ev <- eigen(J, symmetric = TRUE)
  o  <- order(ev$values)
  list(x = ev$values[o], w = 2 * ev$vectors[1, o]^2)
}
.GL64 <- .gauss_legendre(64)

## per-row log-likelihood with a 64-node Gauss-Legendre rule. For a censored row the
## integral over x > w is taken over (lo, hi) = (max(w, x* - 12 sd), max(w, x* + 12 sd)),
## where x* and sd are the centre and SD of the normal kernel f(y | x, z) as a function
## of x, with the map x = lo + (hi - lo) q^2 so that nodes cluster at the lower end (the
## Weibull density is singular at 0 when its shape is below 1 and w can be close to 0).
## Against integrate(rel.tol = 1e-11) the error is below 1e-7 per row, smaller than the
## error of integrate() with its default tolerance (about 1e-5).
.ll_rows_gl <- function(p, s, gl = .GL64) {
  beta <- p[1:3]; sig <- exp(p[4]); g <- p[5:6]; sx <- exp(p[7])
  meanXZ <- exp(g[1] + g[2] * s$z)
  ll <- numeric(length(s$Y)); o <- s$obs
  ll[o] <- dnorm(s$Y[o] - s$Xm[o, , drop = FALSE] %*% beta, 0, sig, log = TRUE) +
           dweibull(s$w[o], shape = 1 / sx, scale = meanXZ[o], log = TRUE)
  c_ <- which(!o)
  es <- s$Y[c_] - (beta[1] + beta[2] * s$A[c_] + beta[3] * s$z[c_])   # e = es + beta1 x
  xs <- -es / beta[2]; sdx <- sig / abs(beta[2])
  lo <- pmax(s$w[c_], xs - 12 * sdx); hi <- pmax(s$w[c_], xs + 12 * sdx)
  q  <- (gl$x + 1) / 2; wq <- gl$w * q                                  # (w / 2) dx/dq, dx/dq = 2 q (hi - lo)
  xx <- outer(hi - lo, q^2) + lo                                        # n_c x 64 nodes
  fv <- dnorm(es + beta[2] * xx, 0, sig) * dweibull(xx, shape = 1 / sx, scale = meanXZ[c_])
  v  <- as.numeric(fv %*% wq) * (hi - lo)
  ll[c_] <- log(pmax(v, 1e-300))
  ll
}

## per-row log-likelihood of the survreg fit for gamma (X | Z Weibull, censored by W)
.ll_rows_survreg_x <- function(g, s) {
  sc <- exp(g[1] + g[2] * s$z); k <- 1 / exp(g[3])
  ifelse(s$obs, dweibull(s$w, k, sc, log = TRUE), pweibull(s$w, k, sc, lower.tail = FALSE, log.p = TRUE))
}

## MLE of theta with gamma fixed (correct, incorrect or two-step)
estimate_mle_fixed <- function(dat, gamma, model = y ~ AW + Z, maxit = 5000) {
  s     <- .mle_data(dat, model)
  gp    <- c(gamma[1:2], log(gamma[3]))
  fit0  <- lm(model, data = dat, weights = D / myp_ywz)
  start <- c(coef(fit0), log(summary(fit0)$sigma))
  ll    <- function(p) { v <- sum(.ll_rows_gl(c(p, gp), s)); if (is.finite(v)) v else -1e10 }
  r1 <- optim(start, ll, method = "Nelder-Mead", control = list(fnscale = -1, maxit = maxit))
  r2 <- tryCatch(optim(r1$par, ll, method = "BFGS", control = list(fnscale = -1)),
                 error = function(e) r1)
  best <- if (r2$value >= r1$value) r2 else r1
  p    <- best$par
  list(beta_est = c(p[1:3], exp(p[4])), gamma_est = gamma, p_full = c(p, gp),
       convcode = best$convergence, loglik = best$value)
}

## sandwich for the MLE; type = "fixed" (gamma known or misspecified, 4 parameters),
## "twostep" (gamma from survreg, stacked) or "joint" (7 parameters)
var_mle_fast <- function(dat, p_full, type = c("fixed", "twostep", "joint"), model = y ~ AW + Z) {
  type <- match.arg(type)
  s  <- .mle_data(dat, model)
  Ur <- numDeriv::jacobian(function(p) .ll_rows_gl(p, s), p_full)             # n x 7 scores
  H  <- numDeriv::hessian(function(p) sum(.ll_rows_gl(p, s)), p_full)         # 7 x 7
  if (type == "fixed") {
    U <- Ur[, 1:4]; A <- H[1:4, 1:4]
  } else if (type == "joint") {
    U <- Ur; A <- H
  } else {
    g  <- p_full[5:7]
    Ug <- numDeriv::jacobian(function(gg) .ll_rows_survreg_x(gg, s), g)
    Hg <- numDeriv::hessian(function(gg) sum(.ll_rows_survreg_x(gg, s)), g)
    U  <- cbind(Ur[, 1:4], Ug)
    A  <- rbind(cbind(H[1:4, 1:4], H[1:4, 5:7]), cbind(matrix(0, 3, 4), Hg))
  }
  Ai <- solve(A)
  V  <- Ai %*% crossprod(U) %*% t(Ai)
  se <- sqrt(diag(V))
  list(se_est = se[1:3], se_sigma = exp(p_full[4]) * se[4], vcov = V)
}

## gamma from survreg(Surv(W, D) ~ Z): the two-step fit of f_{X|Z} (treats C as
## independent of X given Z, which is false here because C depends on Y)
gamma_twostep <- function(dat) {
  wr <- survreg(Surv(W, D) ~ Z, data = dat, dist = "w")
  sv <- c(coef(wr), log(wr$scale))
  if (!(all(is.finite(sv)) && all(abs(sv) < 5))) {          # same guard as the joint MLE
    l0 <- lm(log(W) ~ Z, data = dat, subset = D == 1)
    wr <- survreg(Surv(W, D) ~ Z, data = dat, dist = "w",
                  init = c(coef(l0), log(sd(resid(l0)) * sqrt(6) / pi)))
  }
  c(coef(wr), wr$scale)
}

## ------------------------------------------------------------------------------------
## Inconsistent estimators: complete case, naive, conditional mean imputation (CMI)
## ------------------------------------------------------------------------------------
estimate_cc <- function(dat) {                    # complete case: y ~ (A - W) + Z on D = 1
  f <- .wls_sandwich(dat, y ~ AW + Z, dat$D)
  list(beta_est = c(f$b, sqrt(mean(f$e[dat$D == 1]^2))), se_est = f$se)
}
estimate_naive <- function(dat) {                 # treats W as X
  f <- .wls_sandwich(dat, y ~ AW + Z, rep(1, nrow(dat)))
  list(beta_est = c(f$b, sqrt(mean(f$e^2))), se_est = f$se)
}

## E[X | Y, Z, X > W] for every row, under the true model (as in ODRCC::data_aft):
## int_W^inf x f(y | x, z) f_{X|Z}(x) dx / int_W^inf f(y | x, z) f_{X|Z}(x) dx,
## with beta = (1, 1, 1), sigma = 1 and f_{X|Z} Weibull(gamma). Same Gauss-Legendre rule
## as .ll_rows_gl().
cond_mean_x <- function(dat, beta = c(1, 1, 1), sigma = 1, gamma = gamma_true, gl = .GL64) {
  sc <- exp(gamma[1] + gamma[2] * dat$Z)
  es <- dat$y - (beta[1] + beta[2] * dat$A + beta[3] * dat$Z)
  xs <- -es / beta[2]; sdx <- sigma / abs(beta[2])
  lo <- pmax(dat$W, xs - 12 * sdx); hi <- pmax(dat$W, xs + 12 * sdx)
  q  <- (gl$x + 1) / 2; wq <- gl$w * q
  xx <- outer(hi - lo, q^2) + lo
  fv <- dnorm(es + beta[2] * xx, 0, sigma) * dweibull(xx, shape = 1 / gamma[3], scale = sc)
  as.numeric((fv * xx) %*% wq) / as.numeric(fv %*% wq)
}

## CMI: censored X replaced by E[X | Y, Z, X > C] ("yz") or by E[X | Z] ("z"), then least
## squares of y on (A - X-imputed) + Z, SE from the sandwich with the imputations fixed.
## E[X | Z] = exp(gamma0 + gammaZ Z) Gamma(1 + scale_x). Note: ODRCC::data_aft() uses
## Gamma(1 + 1 / p.x) = Gamma(2.25) instead of Gamma(1 + p.x) = Gamma(1.8), which is not
## the mean of its own X | Z (E[X | Z] is overstated by about 22%).
## which = "z_odrcc" reproduces the ODRCC version.
estimate_cmi <- function(dat, which = c("yz", "z", "z_odrcc"), gamma = gamma_true) {
  which <- match.arg(which)
  ex <- switch(which, yz = cond_mean_x(dat, gamma = gamma),
               z       = exp(gamma[1] + gamma[2] * dat$Z) * gamma(1 + gamma[3]),
               z_odrcc = exp(gamma[1] + gamma[2] * dat$Z) * gamma(1 + 1 / gamma[3]))
  d2 <- dat; d2$AXimp <- ifelse(dat$D == 1, dat$AX, dat$A - ex)
  f  <- .wls_sandwich(d2, y ~ AXimp + Z, rep(1, nrow(d2)))
  list(beta_est = c(f$b, sqrt(mean(f$e^2))), se_est = f$se)
}
