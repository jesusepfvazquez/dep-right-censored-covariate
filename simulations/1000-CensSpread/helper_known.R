## helper_known.R
## Both AIPW estimators with theta = (beta, sigma) when the nuisance distributions are
## known: pi(.) uses the true alpha for C | (Y, Z), and the augmentation uses the true
## f_{X|Z} (Weibull, gamma = (0.1, 0.1, 0.5)). Same estimating equations as
## helper_unguarded_sigma.R and helper_guarded_sigma.R; the differences are
##   * alpha is passed in instead of fit with survreg, and
##   * the sandwich has no correction for estimating alpha (phi_i = G_i).
## Lambda in the guarded estimator is still estimated (it is part of the estimator,
## not a nuisance distribution), from U and V without the alpha correction.

source("helper_guarded_sigma.R")      # .pi_fun, .S_full, .Psi, .working_moments
source("helper_unguarded_sigma.R")    # .G_unguarded, .simpson

## truth in data_aft() (helper_data_aft.R)
alpha_true <- c(1, -0.5, 0.5, 2)      # log C = 1 - 0.5 y + 0.5 Z + 2 log Exp(1)
gamma_true <- c(1, 0.1, 0.8)          # (gamma0, gammaZ, scale_x); shape = 1 / scale_x = 1.25

.newton <- function(Gbar, par, include_sigma) {
  for (it in 1:50) {
    g <- Gbar(par); J <- numDeriv::jacobian(Gbar, par, method = "simple")
    step <- solve(J, g); lam <- 1
    repeat {
      pn <- par - lam * step
      ok <- all(is.finite(pn)) && (!include_sigma || pn[length(pn)] > 0)
      gn <- if (ok) Gbar(pn) else NA
      ok <- ok && all(is.finite(gn))
      if (ok && sum(gn^2) <= sum(g^2)) break
      if (lam < 1e-4) { if (!ok) stop("line search left the parameter space"); break }
      lam <- lam / 2
    }
    conv <- max(abs(pn - par)) < 1e-9; par <- pn; if (conv) break
  }
  list(par = par, iter = it)
}

## unguarded AIPW, known pi and known f_{X|Z}
estimate_unguarded_known <- function(dat, alpha, gamma = gamma_true, include_sigma = TRUE,
                                     do_var = TRUE, upper = 40, m = 400) {   # X has a longer tail here than in data_aft_int()
  grid <- .simpson(upper, m)
  wts  <- dat$D / .pi_fun(alpha, dat$W, dat$y, dat$Z)
  fit  <- lm(y ~ AW + Z, data = dat, weights = wts)
  th0  <- c(coef(fit), sqrt(weighted.mean(resid(fit)^2, wts)))
  k    <- if (include_sigma) 4 else 3
  full <- function(p) if (include_sigma) p else c(p, th0[4])
  Gbar <- function(p) colSums(.G_unguarded(full(p), alpha, dat, gamma, grid, include_sigma))
  nw   <- .newton(Gbar, th0[1:k], include_sigma)
  th   <- full(nw$par)
  out  <- list(beta_est = th, se_est = rep(NA, k), iter = nw$iter)
  if (!do_var) return(out)
  Gi <- .G_unguarded(th, alpha, dat, gamma, grid, include_sigma)
  Ai <- solve(numDeriv::jacobian(Gbar, nw$par))
  out$se_est <- sqrt(diag(Ai %*% crossprod(Gi) %*% t(Ai)))
  out
}

## guarded AIPW, known pi; working normal = moments of the true f_{X|Z}
estimate_guarded_known <- function(dat, alpha, gamma = gamma_true, include_sigma = TRUE, do_var = TRUE) {
  mm  <- .working_moments(dat, "mle", gamma); m <- mm$m; sd <- mm$sd
  wts <- dat$D / .pi_fun(alpha, dat$W, dat$y, dat$Z)
  fit <- lm(y ~ AW + Z, data = dat, weights = wts)
  th0 <- c(coef(fit), sqrt(weighted.mean(resid(fit)^2, wts)))
  k   <- if (include_sigma) 4 else 3
  pieces <- function(th) {
    S  <- .S_full(th, dat$y, dat$A, dat$W, dat$Z, include_sigma); S[dat$D == 0, ] <- 0
    Ps <- .Psi(th, dat$y, dat$A, dat$Z, m, sd, include_sigma)
    list(U = S * wts, V = Ps * (1 - wts))
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
