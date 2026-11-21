## Shared roxygen blocks -------------------------------------------------------------

#' Shared arguments
#'
#' @param data_yXZ A data frame with one row per subject. It must contain the
#'   outcome, the covariates in \code{model}, the age \code{a_var}, the observed
#'   \code{W = min(X, C)}, and the indicator \code{D = I(X <= C)}.
#' @param model A formula for the outcome model, \code{y ~ AW + Z1 + ...}, with
#'   the time-to-event covariate \code{aw_var} on the right-hand side. The outcome
#'   model is \eqn{Y = \beta_0 + \beta_1 (A - X) + \beta_Z' Z + \epsilon},
#'   \eqn{\epsilon \sim N(0, \sigma^2)}.
#' @param aw_var Name of the column \code{A - W} in \code{model}. Defaults to \code{"AW"}.
#' @param a_var Name of the age column \code{A}. Defaults to \code{"A"}.
#' @name odrcc-args
#' @keywords internal
NULL

#' Shared arguments: censoring model
#'
#' @param model_weights A one-sided formula for the Weibull AFT model of
#'   \eqn{C | (Y, Z)}, for example \code{~ y + Z}. Defaults to the outcome and the
#'   covariates of \code{model} other than \code{aw_var}.
#' @name odrcc-args-weights
#' @keywords internal
NULL

#' Shared arguments: covariate model
#'
#' @param model_xz A one-sided formula for the Weibull AFT model of
#'   \eqn{X | Z}, for example \code{~ Z}. Defaults to the covariates of
#'   \code{model} other than \code{aw_var}.
#' @name odrcc-args-xz
#' @keywords internal
NULL

.out <- function(fit, s, extra = list()) {
  c(list(beta_est = .name(fit$est, s), se_est = .name(fit$se, s)), extra)
}

## Estimators --------------------------------------------------------------------------

#' Least-squares estimator with a sandwich standard error
#'
#' Fits \code{model} by least squares to all rows of \code{data_yXZ}. Use it
#' with the true covariate (the oracle estimator, e.g. \code{y ~ AX + Z}) or with
#' \code{W} in place of \code{X} (the naive estimator, \code{y ~ AW + Z}).
#'
#' @param data_yXZ A data frame containing the variables in \code{model}.
#' @param model A formula, e.g. \code{y ~ AX + Z}.
#'
#' @return A list with \code{beta_est} (coefficients and \code{sigma}) and
#'   \code{se_est} (sandwich standard errors; \code{NA} for \code{sigma}).
#' @examples
#' set.seed(1)
#' dat <- data_aft(500)
#' estimate_beta(dat, y ~ AX + Z)
#' @importFrom stats lm model.matrix resid
#' @export
estimate_beta <- function(data_yXZ, model) {
  f <- stats::lm(model, data = data_yXZ)
  X <- stats::model.matrix(f); e <- stats::resid(f)
  Ai <- solve(crossprod(X)); V <- Ai %*% crossprod(X * e) %*% Ai
  list(beta_est = c(stats::coef(f), sigma = sqrt(mean(e^2))),
       se_est = c(sqrt(diag(V)), sigma = NA))
}

#' Complete-case estimator
#'
#' Least squares on the rows with \code{D = 1}. Under outcome-dependent
#' right-censoring, the complete-case estimator is not consistent; it is included
#' for comparison.
#'
#' @inheritParams odrcc-args
#' @return A list with \code{beta_est} (coefficients and \code{sigma}) and
#'   \code{se_est} (sandwich standard errors; \code{NA} for \code{sigma}).
#' @examples
#' set.seed(1)
#' dat <- data_aft(500)
#' estimate_beta_cc(dat, y ~ AW + Z)
#' @export
estimate_beta_cc <- function(data_yXZ, model, aw_var = "AW", a_var = "A") {
  s <- .odrcc_setup(data_yXZ, model, aw_var, a_var)
  o <- s$D == 1; X <- s$Xw[o, , drop = FALSE]; y <- s$y[o]
  b <- solve(crossprod(X), crossprod(X, y)); e <- drop(y - X %*% b)
  Ai <- solve(crossprod(X)); V <- Ai %*% crossprod(X * e) %*% Ai
  .out(list(est = c(drop(b), sqrt(mean(e^2))), se = c(sqrt(diag(V)), NA)), s)
}

#' IPW estimator
#'
#' Inverse probability weighting estimator of \eqn{\theta = (\beta, \sigma)}.
#' The weights are \eqn{D / \hat\pi(W, Y, Z)}, where
#' \eqn{\hat\pi(t, Y, Z) = \hat P(C > t | Y, Z)} comes from a Weibull AFT model for
#' \eqn{C | (Y, Z)} (\code{model_weights}). The IPW estimator is consistent when
#' the model for \eqn{C | (Y, Z)} and the mean of the outcome model are correctly
#' specified. The sandwich standard errors account for estimating the weights.
#'
#' @inheritParams odrcc-args
#' @inheritParams odrcc-args-weights
#' @return A list with
#' \describe{
#'   \item{beta_est}{Estimates of \eqn{(\beta, \sigma)}.}
#'   \item{se_est}{Sandwich standard errors.}
#'   \item{alpha_est}{Estimated coefficients and scale of the Weibull AFT model for \eqn{C | (Y, Z)}.}
#' }
#' @examples
#' set.seed(1)
#' dat <- data_aft(500)
#' estimate_beta_ipw(dat, y ~ AW + Z, model_weights = ~ y + Z)
#' @export
estimate_beta_ipw <- function(data_yXZ, model, model_weights = NULL, aw_var = "AW", a_var = "A") {
  s <- .odrcc_setup(data_yXZ, model, aw_var, a_var, model_weights = model_weights)
  alpha <- .fit_alpha(s)
  .out(.ipw(s, alpha), s, list(alpha_est = alpha))
}

#' Joint maximum likelihood estimator
#'
#' Estimates \eqn{\theta = (\beta, \sigma)} and the Weibull AFT model for
#' \eqn{X | Z} (\code{model_xz}) jointly, by maximizing the observed-data
#' likelihood. For a censored row, the likelihood integrates
#' \eqn{f_{Y|X,Z} f_{X|Z}} over \eqn{x > W} with a 64-node Gauss-Legendre rule.
#' Because the model for \eqn{X | Z} is estimated jointly with the outcome model,
#' the joint MLE remains consistent under outcome-dependent censoring (a model
#' for \eqn{X | Z} fitted to \eqn{(W, D, Z)} alone would not be). The joint MLE
#' does not use a model for \eqn{C | (Y, Z)}; \code{model_weights} is used only
#' for the starting values. It is consistent when \eqn{f_{Y|X,Z}} and
#' \eqn{f_{X|Z}} are correctly specified.
#'
#' @inheritParams odrcc-args
#' @inheritParams odrcc-args-xz
#' @inheritParams odrcc-args-weights
#' @param maxit Maximum number of Nelder-Mead iterations (followed by BFGS).
#' @return A list with
#' \describe{
#'   \item{beta_est}{Estimates of \eqn{(\beta, \sigma)}.}
#'   \item{se_est}{Sandwich standard errors.}
#'   \item{gamma_est}{Estimated coefficients and scale of the Weibull AFT model for \eqn{X | Z}.}
#'   \item{gamma_se}{Their sandwich standard errors.}
#'   \item{convergence}{Convergence code from \code{\link[stats]{optim}} (0 = converged).}
#'   \item{loglik}{Maximized log-likelihood.}
#' }
#' @examples
#' set.seed(1)
#' dat <- data_aft(500)
#' estimate_beta_mle(dat, y ~ AW + Z, model_xz = ~ Z)
#' @export
estimate_beta_mle <- function(data_yXZ, model, model_xz = NULL, model_weights = NULL,
                              aw_var = "AW", a_var = "A", maxit = 5000) {
  s <- .odrcc_setup(data_yXZ, model, aw_var, a_var, model_xz, model_weights)
  f <- .mle(s, .fit_alpha(s), maxit = maxit)
  .out(f, s, list(gamma_est = f$gamma, gamma_se = f$gamma_se,
                  convergence = f$convergence, loglik = f$loglik))
}

#' Unguarded AIPW estimator
#'
#' Augmented IPW estimator of \eqn{\theta = (\beta, \sigma)} with the efficient
#' augmentation function \eqn{\Psi^{eff}_{AIPW}(Y, Z; \theta)} of the manuscript,
#' computed by numerical integration over a working model for \eqn{f_{X|Y,Z}}:
#' \describe{
#'   \item{\code{fx_method = "joint"}}{(default, the estimator in the manuscript)
#'     \eqn{f_{X|Y,Z} \propto f_{Y|X,Z}(y | x, z; \theta) f_{X|Z}(x; \hat\gamma)}, with
#'     \eqn{\hat\gamma} from the joint MLE (\code{\link{estimate_beta_mle}}) unless
#'     \code{gamma} is supplied.}
#'   \item{\code{fx_method = "direct"}}{\eqn{f_{X|Y,Z}} is a Weibull AFT model,
#'     \code{survreg(Surv(W, D) ~ y + Z)}, fitted to all rows; the fit is valid
#'     because \eqn{C} is independent of \eqn{X} given \eqn{(Y, Z)}.}
#' }
#' The unguarded AIPW estimator is consistent when the model for
#' \eqn{C | (Y, Z)} and the mean of the outcome model are correctly specified,
#' whatever the working model for \eqn{f_{X|Y,Z}}. It is more efficient than the
#' IPW estimator when the working model is correct, but it has no such guarantee
#' otherwise (see \code{\link{estimate_beta_aipw_lambda}}). The sandwich standard
#' errors account for estimating the weights.
#'
#' @inheritParams odrcc-args
#' @inheritParams odrcc-args-xz
#' @inheritParams odrcc-args-weights
#' @param fx_method \code{"joint"} or \code{"direct"}; see Details.
#' @param gamma Optional coefficients and scale of the Weibull model for
#'   \eqn{X | Z} (as in \code{estimate_beta_mle()$gamma_est}). Used only with
#'   \code{fx_method = "joint"}; if \code{NULL}, the joint MLE is fitted first.
#' @param tau Optional upper limit of the integral over \eqn{X}. The estimator
#'   stays consistent for any \code{tau}; a finite \code{tau} only changes the
#'   working model. Useful when the integral is unstable in the right tail.
#' @param start Optional starting value for \eqn{(\beta, \sigma)}. Defaults to the
#'   IPW estimate.
#' @return A list with \code{beta_est}, \code{se_est}, \code{converged}, and
#'   \code{fx_par} (the parameters of the working model for \eqn{f_{X|Y,Z}}).
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data_aft(500)
#' estimate_beta_aipw(dat, y ~ AW + Z, model_xz = ~ Z, model_weights = ~ y + Z)
#' estimate_beta_aipw(dat, y ~ AW + Z, fx_method = "direct")
#' }
#' @export
estimate_beta_aipw <- function(data_yXZ, model, model_xz = NULL, model_weights = NULL,
                               fx_method = c("joint", "direct"), gamma = NULL, tau = NULL,
                               start = NULL, aw_var = "AW", a_var = "A") {
  fx_method <- match.arg(fx_method)
  s <- .odrcc_setup(data_yXZ, model, aw_var, a_var, model_xz, model_weights)
  alpha <- .fit_alpha(s)
  if (fx_method == "joint" && is.null(gamma)) gamma <- .mle(s, alpha)$gamma
  f <- .aipw_ung(s, alpha, gamma, fx_method, tau, start)
  .out(f, s, list(converged = f$converged, fx_par = f$fx_par))
}

#' Guarded AIPW estimator
#'
#' Augmented IPW estimator of \eqn{\theta = (\beta, \sigma)} with a closed-form
#' augmentation \eqn{\Psi_{closed}} (a working normal model for \eqn{X | Z} with
#' the mean and standard deviation of the fitted Weibull model for \eqn{X | Z})
#' and the matrix \eqn{\Lambda} that minimizes the asymptotic variance within the
#' class \eqn{\Lambda \Psi_{closed}}. Because of \eqn{\Lambda}, the guarded AIPW
#' estimator is at least as efficient as the IPW estimator asymptotically,
#' whether or not the working model is correct. It is consistent when the model
#' for \eqn{C | (Y, Z)} and the mean of the outcome model are correctly specified.
#' \eqn{\Lambda} is estimated once, at the IPW estimate, and held fixed. The
#' sandwich standard errors account for estimating the weights.
#'
#' @inheritParams odrcc-args
#' @inheritParams odrcc-args-xz
#' @inheritParams odrcc-args-weights
#' @param gamma Optional coefficients and scale of the Weibull model for
#'   \eqn{X | Z} that supply the working moments (as in
#'   \code{estimate_beta_mle()$gamma_est}). If \code{NULL}, the joint MLE is
#'   fitted first.
#' @return A list with \code{beta_est}, \code{se_est}, \code{converged}, and
#'   \code{Lambda}.
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- data_aft(500)
#' estimate_beta_aipw_lambda(dat, y ~ AW + Z, model_xz = ~ Z, model_weights = ~ y + Z)
#' }
#' @export
estimate_beta_aipw_lambda <- function(data_yXZ, model, model_xz = NULL, model_weights = NULL,
                                      gamma = NULL, aw_var = "AW", a_var = "A") {
  s <- .odrcc_setup(data_yXZ, model, aw_var, a_var, model_xz, model_weights)
  alpha <- .fit_alpha(s)
  if (is.null(gamma)) gamma <- .mle(s, alpha)$gamma
  f <- .aipw_grd(s, alpha, gamma)
  .out(f, s, list(converged = f$converged, Lambda = f$Lambda))
}

#' Effective sample size of the IPW weights
#'
#' Computes the Kish effective sample size of the weights \eqn{1 / \hat\pi(W, Y, Z)}
#' among the uncensored rows, \eqn{(\sum w)^2 / \sum w^2}, and reports it as a
#' fraction of the number of uncensored rows. In the manuscript's simulations, the
#' IPW and unguarded AIPW estimators had small bias and near-nominal coverage when
#' this fraction was at least about 65\%, and the guarded AIPW estimator when it was
#' at least about 80\%; below these values, the joint MLE is preferable when the model
#' for \eqn{X | Z} is credible.
#'
#' @param data_yXZ A data frame with \code{W}, \code{D}, and the variables in
#'   \code{model_weights}.
#' @param model_weights A one-sided formula for the Weibull AFT model of
#'   \eqn{C | (Y, Z)}, e.g. \code{~ y + Z}.
#' @return A list with \code{ess} (effective sample size), \code{ess_frac}
#'   (\code{ess} divided by the number of uncensored rows), \code{n_uncensored},
#'   \code{max_weight}, and \code{max_weight_share} (largest weight over the sum of
#'   the weights).
#' @examples
#' set.seed(1)
#' dat <- data_aft(1000, c_int = -1.5)   # about 85% right-censoring
#' ess_weights(dat, ~ y + Z)
#' @importFrom stats model.matrix
#' @export
ess_weights <- function(data_yXZ, model_weights) {
  dat <- as.data.frame(data_yXZ)
  s <- list(W = dat$W, D = dat$D, Zc = stats::model.matrix(model_weights, dat))
  alpha <- .aft(s$W, 1 - s$D, s$Zc)
  w <- 1 / .pi_t(alpha, s$W, s)[s$D == 1]
  ess <- sum(w)^2 / sum(w^2)
  list(ess = ess, ess_frac = ess / length(w), n_uncensored = length(w),
       max_weight = max(w), max_weight_share = max(w) / sum(w))
}
