## Internal helpers shared by the estimators.
##
## Models
##   outcome    Y = b0 + b1 (A - X) + Zy' bz + e,  e ~ N(0, sigma^2)
##   f_{X|Z}    Weibull AFT:  log X = Zx' gamma + sx * (min-Gumbel error)
##   f_{C|Y,Z}  Weibull AFT:  log C = Zc' alpha + sc * (min-Gumbel error)
##   pi(t) = P(C > t | Y, Z), under the assumption that C is independent of X given (Y, Z).
## theta = (b0, b1, bz, sigma) for every estimator.

## ---------------------------------------------------------------------------------
## setup: data and formulas in, matrices out
## ---------------------------------------------------------------------------------
.odrcc_setup <- function(data_yXZ, model, aw_var = "AW", a_var = "A",
                         model_xz = NULL, model_weights = NULL) {
  if (!inherits(model, "formula") || length(model) != 3L)
    stop("'model' must be a two-sided formula, e.g. y ~ AW + Z.")
  dat <- as.data.frame(data_yXZ)
  y_var <- all.vars(model[[2]])
  if (length(y_var) != 1L) stop("The left-hand side of 'model' must be a single variable.")
  rhs <- attr(stats::terms(model), "term.labels")
  if (!aw_var %in% rhs)
    stop("'", aw_var, "' (the time-to-event covariate A - W) must appear on the right-hand side of 'model'.")
  z_terms <- setdiff(rhs, aw_var)
  if (is.null(model_xz))
    model_xz <- if (length(z_terms)) stats::reformulate(z_terms) else ~ 1
  if (is.null(model_weights))
    model_weights <- stats::reformulate(c(y_var, z_terms))
  if (!inherits(model_xz, "formula") || !inherits(model_weights, "formula"))
    stop("'model_xz' and 'model_weights' must be formulas, e.g. ~ Z and ~ y + Z.")
  model_xz      <- stats::reformulate(attr(stats::terms(model_xz), "term.labels"),
                                      intercept = TRUE)
  model_weights <- stats::reformulate(attr(stats::terms(model_weights), "term.labels"),
                                      intercept = TRUE)

  need <- unique(c(y_var, a_var, "W", "D", aw_var, all.vars(stats::reformulate(c(z_terms, "1"))),
                   all.vars(model_xz), all.vars(model_weights)))
  miss <- setdiff(need, names(dat))
  if (length(miss)) stop("Columns not found in data_yXZ: ", paste(miss, collapse = ", "))
  keep <- stats::complete.cases(dat[, need, drop = FALSE])
  if (any(!keep)) message("Dropping ", sum(!keep), " rows with missing values.")
  dat <- dat[keep, , drop = FALSE]

  s <- list(y = as.numeric(dat[[y_var]]), A = as.numeric(dat[[a_var]]),
            W = as.numeric(dat$W), D = as.numeric(dat$D))
  if (!all(s$D %in% c(0, 1))) stop("'D' must be 0 (censored) or 1 (observed).")
  if (any(s$W <= 0)) stop("'W' must be positive.")
  if (max(abs(dat[[aw_var]] - (s$A - s$W))) > 1e-8)
    stop("'", aw_var, "' must equal ", a_var, " - W.")
  s$Zy <- if (length(z_terms))
    stats::model.matrix(stats::reformulate(z_terms), dat)[, -1, drop = FALSE] else
    matrix(numeric(0), nrow(dat), 0)
  s$Zx <- stats::model.matrix(model_xz, dat)
  s$Zc <- stats::model.matrix(model_weights, dat)
  s$n  <- length(s$y); s$pz <- ncol(s$Zy); s$kb <- 2 + s$pz
  s$iz <- 2 + seq_len(s$pz)
  s$names <- c("(Intercept)", aw_var, colnames(s$Zy), "sigma")
  s$Xw <- cbind(1, s$A - s$W, s$Zy)                 # outcome design at x = W
  s$model_xz <- model_xz; s$model_weights <- model_weights
  s
}

.name <- function(x, s) { names(x) <- s$names; x }

## ---------------------------------------------------------------------------------
## nuisance models
## ---------------------------------------------------------------------------------
## Weibull AFT on a design matrix (first column the intercept): (coefficients, scale)
.aft <- function(time, event, M, weights = NULL, init = NULL) {
  f <- suppressWarnings(survival::survreg(survival::Surv(time, event) ~ M - 1, weights = weights, init = init,
                         dist = "weibull", control = survival::survreg.control(maxiter = 100)))
  co <- stats::coef(f); names(co) <- colnames(M)
  out <- c(co, scale = f$scale)
  if (!all(is.finite(out)) || f$iter >= 100)
    stop("The Weibull AFT fit did not converge.")
  out
}
## P(C > t | Y, Z); t a vector (length n) or an n x K matrix
.pi_t <- function(alpha, t, s) {
  k <- length(alpha); lp <- drop(s$Zc %*% alpha[-k])
  exp(-exp((log(t) - lp) / alpha[k]))
}
## per-row log-likelihood of C | (Y, Z)
.alpha_ll <- function(alpha, s) {
  k <- length(alpha); lp <- drop(s$Zc %*% alpha[-k])
  u <- (log(s$W) - lp) / alpha[k]
  Sv <- exp(-exp(u)); f <- (1 / alpha[k]) * exp(u) * Sv
  log((1 - s$D) * f + s$D * Sv)
}
## n x k: rows of the influence function of alpha-hat
.alpha_corr <- function(alpha, s) {
  sa <- numDeriv::jacobian(function(a) .alpha_ll(a, s), alpha)
  Ha <- numDeriv::hessian(function(a) sum(.alpha_ll(a, s)), alpha) / s$n
  sa %*% t(-solve(Ha))
}
.fit_alpha <- function(s) .aft(s$W, 1 - s$D, s$Zc)

## ---------------------------------------------------------------------------------
## full-data score of theta at x
## ---------------------------------------------------------------------------------
.ez <- function(th, s) drop(s$y - th[1] - th[2] * s$A - s$Zy %*% th[s$iz])   # e = ez + b1 x
.S_at <- function(th, s, x) {
  sg <- th[s$kb + 1]; e <- .ez(th, s) + th[2] * x
  cbind(cbind(e, e * (s$A - x), e * s$Zy) / sg^2, -1 / sg + e^2 / sg^3)
}

## Newton solve of colSums(G(th)) = 0 with step halving
.newton <- function(G, th0, maxit = 50) {
  par <- th0; k <- length(par); conv <- FALSE
  for (it in seq_len(maxit)) {
    g <- colSums(G(par))
    J <- numDeriv::jacobian(function(p) colSums(G(p)), par, method = "simple")
    step <- solve(J, g); lam <- 1
    repeat {
      pn <- par - lam * step
      ok <- all(is.finite(pn)) && pn[k] > 0
      gn <- if (ok) colSums(G(pn)) else NA
      ok <- ok && all(is.finite(gn))
      if (ok && sum(gn^2) <= sum(g^2)) break
      if (lam < 1e-4) { if (!ok) stop("Newton line search left the parameter space."); break }
      lam <- lam / 2
    }
    conv <- max(abs(pn - par)) < 1e-9; par <- pn
    if (conv) break
  }
  list(par = par, iter = it, converged = conv)
}

## sandwich SEs for estimating functions G (n x k), corrected for estimating alpha
.sandwich_alpha <- function(G, th, alpha, s, Gfun_alpha) {
  Gi  <- G(th)
  dGa <- numDeriv::jacobian(function(a) colMeans(Gfun_alpha(th, a)), alpha)
  phi <- Gi + .alpha_corr(alpha, s) %*% t(dGa)
  Ai  <- solve(numDeriv::jacobian(function(p) colSums(G(p)), th))
  sqrt(diag(Ai %*% crossprod(phi) %*% t(Ai)))
}

## ---------------------------------------------------------------------------------
## IPW (also the start for the other estimators)
## ---------------------------------------------------------------------------------
.ipw <- function(s, alpha) {
  Gf <- function(th, a) {
    w <- s$D / .pi_t(a, s$W, s); e <- drop(s$y - s$Xw %*% th[1:s$kb])
    cbind(s$Xw * (w * e), w * (e^2 - th[s$kb + 1]^2))    # sigma: weighted mean of e^2
  }
  w  <- s$D / .pi_t(alpha, s$W, s)
  b  <- drop(solve(crossprod(s$Xw * w, s$Xw), crossprod(s$Xw * w, s$y)))
  th <- c(b, sqrt(stats::weighted.mean(drop(s$y - s$Xw %*% b)^2, w)))
  se <- .sandwich_alpha(function(p) Gf(p, alpha), th, alpha, s, Gf)
  list(est = th, se = se)
}

## ---------------------------------------------------------------------------------
## joint MLE of (theta, gamma); p = (b0, b1, bz, log sigma, gamma, log sx)
## censored-row integrals: 64-node Gauss-Legendre rule around the peak of f(y | x, z)
## ---------------------------------------------------------------------------------
.gauss_legendre <- function(k) {
  j <- seq_len(k - 1); b <- j / sqrt(4 * j^2 - 1)
  J <- matrix(0, k, k); J[cbind(j, j + 1)] <- b; J[cbind(j + 1, j)] <- b
  ev <- eigen(J, symmetric = TRUE); o <- order(ev$values)
  list(x = ev$values[o], w = 2 * ev$vectors[1, o]^2)
}
.ll_rows_mle <- function(p, s, gl) {
  kb <- s$kb; kg <- ncol(s$Zx)
  th <- c(p[1:kb], exp(p[kb + 1])); g <- p[kb + 1 + seq_len(kg)]; sx <- exp(p[kb + 2 + kg])
  sig <- th[kb + 1]; b1 <- th[2]
  ez <- .ez(th, s); scx <- exp(drop(s$Zx %*% g))
  ll <- numeric(s$n); o <- s$D == 1
  ll[o] <- stats::dnorm(ez[o] + b1 * s$W[o], 0, sig, log = TRUE) +
           stats::dweibull(s$W[o], shape = 1 / sx, scale = scx[o], log = TRUE)
  c_ <- which(!o)
  if (length(c_)) {
    xs <- -ez[c_] / b1; sdx <- sig / abs(b1)
    lo <- pmax(s$W[c_], xs - 12 * sdx); hi <- pmax(s$W[c_], xs + 12 * sdx)
    q  <- (gl$x + 1) / 2; wq <- gl$w * q
    xx <- outer(hi - lo, q^2) + lo
    fv <- stats::dnorm(ez[c_] + b1 * xx, 0, sig) *
          stats::dweibull(xx, shape = 1 / sx, scale = scx[c_])
    ll[c_] <- log(pmax(as.numeric(fv %*% wq) * (hi - lo), 1e-300))
  }
  ll
}
.mle <- function(s, alpha, maxit = 5000) {
  gl <- .gauss_legendre(64)
  ## start: IPW theta; IPW-weighted Weibull fit of the uncensored W on Zx for gamma
  ip <- .ipw(s, alpha)$est
  o  <- s$D == 1; wo <- 1 / .pi_t(alpha, s$W, s)[o]
  gx <- .aft(s$W[o], rep(1, sum(o)), s$Zx[o, , drop = FALSE], weights = wo)
  kg <- ncol(s$Zx); kb <- s$kb
  start <- c(ip[1:kb], log(ip[kb + 1]), gx[seq_len(kg)], log(gx[kg + 1]))
  nll <- function(p) { v <- -sum(.ll_rows_mle(p, s, gl)); if (is.finite(v)) v else 1e10 }
  r1 <- stats::optim(start, nll, method = "Nelder-Mead", control = list(maxit = maxit))
  r2 <- tryCatch(stats::optim(r1$par, nll, method = "BFGS", control = list(maxit = 1000)),
                 error = function(e) r1)
  r  <- if (r2$value <= r1$value) r2 else r1
  p  <- r$par
  U  <- numDeriv::jacobian(function(pp) .ll_rows_mle(pp, s, gl), p)
  H  <- numDeriv::hessian(function(pp) sum(.ll_rows_mle(pp, s, gl)), p)
  Ai <- solve(H); V <- Ai %*% crossprod(U) %*% t(Ai)
  se <- sqrt(diag(V)); sig <- exp(p[kb + 1])
  gam <- c(p[kb + 1 + seq_len(kg)], exp(p[kb + 2 + kg]))
  names(gam) <- c(colnames(s$Zx), "scale")
  gse <- c(se[kb + 1 + seq_len(kg)], exp(p[kb + 2 + kg]) * se[kb + 2 + kg])
  names(gse) <- names(gam)
  list(est = c(p[1:kb], sig), se = c(se[1:kb], sig * se[kb + 1]),
       gamma = gam, gamma_se = gse, convergence = r$convergence, loglik = -r$value)
}

## ---------------------------------------------------------------------------------
## unguarded AIPW: Psi = int S(theta; t) w(t) dt / int w(t) dt on a Simpson grid,
##   w(t) = f(y | t, z; theta) f_{X|Z}(t; gamma) {1 - 1/pi(t)}
## ---------------------------------------------------------------------------------
.simpson <- function(upper, m) {
  t <- seq(0, upper, length.out = 2 * m + 1); h <- t[2] - t[1]
  w <- rep(c(2, 4), length.out = 2 * m + 1); w[1] <- 1; w[2 * m + 1] <- 1
  list(t = t[-1], w = (h / 3 * w)[-1])
}
.weib_sd <- function(scale, shape) scale * sqrt(gamma(1 + 2 / shape) - gamma(1 + 1 / shape)^2)
## grid on (0, upper]: upper covers every row's peak of f(y | x, z) plus 12 SDs; spacing is
## 1/20 of the narrower of that kernel's SD (sigma / |b1|) and the smallest SD of f_{X|Z}
.grid_for <- function(s, th, sdx_min = Inf, tau = NULL) {
  sdk <- th[s$kb + 1] / abs(th[2]); xs <- -.ez(th, s) / th[2]
  upper <- max(1.5 * max(s$W), 1.5 * (max(xs) + 12 * sdk))
  if (!is.null(tau)) upper <- tau
  h <- min(max(min(sdk, sdx_min) / 20, 0.005), 0.5)
  .simpson(upper, max(200, ceiling(upper / (2 * h))))
}
## fx(th, tt): n x K density of X | (Y, Z) on the grid, up to a constant per row
.G_ung <- function(th, alpha, s, grid, fx) {
  sg <- th[s$kb + 1]; ez <- .ez(th, s)
  S  <- .S_at(th, s, s$W); S[s$D == 0, ] <- 0
  tt <- matrix(grid$t, s$n, length(grid$t), byrow = TRUE)
  et <- ez + th[2] * tt
  wt <- fx(th, tt) * (1 - 1 / .pi_t(alpha, tt, s))
  wt[!is.finite(wt)] <- 0
  wq <- sweep(wt, 2, grid$w, `*`); den <- rowSums(wq)
  E1 <- rowSums(wq * et) / den
  Ps <- cbind(cbind(E1, rowSums(wq * et * (s$A - tt)) / den, s$Zy * E1) / sg^2,
              -1 / sg + rowSums(wq * et^2) / den / sg^3)
  wD <- s$D / .pi_t(alpha, s$W, s)
  S * wD + Ps * (1 - wD)
}
## fx_method "joint": f_{X|Y,Z} proportional to f(y | x, z; theta) f_{X|Z}(x; gamma)
## fx_method "direct": f_{X|Y,Z} a Weibull AFT, survreg(Surv(W, D) ~ Y + Zx), fitted to all
##   rows (valid because C is independent of X given (Y, Z)); held fixed
.aipw_ung <- function(s, alpha, gamma = NULL, fx_method = "joint", tau = NULL, start = NULL) {
  th0 <- if (is.null(start)) .ipw(s, alpha)$est else start
  if (fx_method == "direct") {
    M  <- cbind(s$Zx[, 1, drop = FALSE], y = s$y, s$Zx[, -1, drop = FALSE])
    o  <- s$D == 1                                  # start from the fit to the uncensored rows
    d0 <- .aft(s$W[o], rep(1, sum(o)), M[o, , drop = FALSE])
    dl <- tryCatch(.aft(s$W, s$D, M, init = d0[-length(d0)]), error = function(e)
      stop("The direct Weibull fit for X | (Y, Z) did not converge; use fx_method = \"joint\".",
           call. = FALSE))
    k  <- length(dl); scd <- exp(drop(M %*% dl[-k]))
    fx <- function(th, tt) stats::dweibull(tt, shape = 1 / dl[k], scale = scd)
    sdx_min <- min(.weib_sd(scd, 1 / dl[k])); fx_par <- dl
  } else {
    k  <- length(gamma); scx <- exp(drop(s$Zx %*% gamma[-k]))
    fx <- function(th, tt) stats::dnorm(.ez(th, s) + th[2] * tt, 0, th[s$kb + 1]) *
                           stats::dweibull(tt, shape = 1 / gamma[k], scale = scx)
    sdx_min <- min(.weib_sd(scx, 1 / gamma[k])); fx_par <- gamma
  }
  grid <- .grid_for(s, th0, sdx_min, tau)
  G  <- function(p) .G_ung(p, alpha, s, grid, fx)
  nw <- .newton(G, th0)
  se <- .sandwich_alpha(G, nw$par, alpha, s, function(p, a) .G_ung(p, a, s, grid, fx))
  list(est = nw$par, se = se, converged = nw$converged, fx_par = fx_par)
}

## ---------------------------------------------------------------------------------
## guarded AIPW: closed-form Psi under a working normal X | Z ~ N(m, sd^2) with the
## moments of the fitted Weibull f_{X|Z}; Lambda estimated at the IPW start and held fixed
## ---------------------------------------------------------------------------------
.Psi_closed <- function(th, s, m, sd) {
  b1 <- th[2]; sg <- th[s$kb + 1]; es <- .ez(th, s)
  a  <- -1 / (2 * sd^2) - b1^2 / (2 * sg^2)
  bb <- m / sd^2 - b1 * es / sg^2
  mu <- -bb / (2 * a); s2 <- (sd^2 * sg^2) / (sg^2 + b1^2 * sd^2)
  Ee   <- es + b1 * mu
  EeAX <- -b1 * (s2 + mu^2) + mu * (s$A * b1 - es) + s$A * es
  Ee2  <- es^2 + 2 * b1 * es * mu + b1^2 * (s2 + mu^2)
  cbind(cbind(Ee, EeAX, s$Zy * Ee) / sg^2, -1 / sg + Ee2 / sg^3)
}
.aipw_grd <- function(s, alpha, gamma) {
  k  <- length(gamma); sc <- exp(drop(s$Zx %*% gamma[-k])); kk <- 1 / gamma[k]
  m  <- sc * base::gamma(1 + 1 / kk)
  sd <- sc * sqrt(base::gamma(1 + 2 / kk) - base::gamma(1 + 1 / kk)^2)
  th0 <- .ipw(s, alpha)$est
  parts <- function(th, a) {
    wD <- s$D / .pi_t(a, s$W, s); S <- .S_at(th, s, s$W); S[s$D == 0, ] <- 0
    list(U = S * wD, V = .Psi_closed(th, s, m, sd) * (1 - wD))
  }
  ## Lambda = -E(U V') E(V V')^{-1}, U and V with their alpha corrections
  cr <- .alpha_corr(alpha, s); p0 <- parts(th0, alpha)
  dU <- numDeriv::jacobian(function(a) colMeans(parts(th0, a)$U), alpha)
  dV <- numDeriv::jacobian(function(a) colMeans(parts(th0, a)$V), alpha)
  U <- p0$U + cr %*% t(dU); V <- p0$V + cr %*% t(dV)
  Lam <- -crossprod(U, V) %*% solve(crossprod(V, V))
  Gf <- function(th, a) { p <- parts(th, a); p$U + p$V %*% t(Lam) }
  G  <- function(th) Gf(th, alpha)
  nw <- .newton(G, th0)
  se <- .sandwich_alpha(G, nw$par, alpha, s, Gf)
  dimnames(Lam) <- list(s$names, s$names)
  list(est = nw$par, se = se, converged = nw$converged, Lambda = Lam)
}
