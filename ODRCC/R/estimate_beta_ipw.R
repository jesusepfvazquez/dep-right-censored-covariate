#' Inverse probability weighted (IPW) estimator for regression parameters
#'
#' Fits a regression model using inverse probability weighting, where the
#' weights are derived from a parametric model for the censoring distribution
#' \code{C | (Y, Z, ...)}. The user supplies a right-hand-side formula (e.g.
#' \code{~ y + Z}) which is used to model \code{Surv(W, 1 - D)} via a Weibull
#' AFT model. The resulting estimated survival probabilities at \code{W} are
#' used as weights in the IPW estimating equations.
#'
#' @param data_yXZ A data frame containing the variables in \code{model},
#'   as well as \code{W} (observed \code{min(X, C)}), \code{D} (indicator
#'   \code{I(X <= C)}), and the covariates appearing in \code{model_weights}.
#' @param model A \code{\link[stats]{formula}} specifying the outcome regression
#'   model (e.g. \code{y ~ AW + Z}).
#' @param model_weights A right-hand-side formula specifying the variables in
#'   the censoring model, e.g. \code{~ y + Z}. This will be expanded to
#'   \code{Surv(W, 1 - D) ~ y + Z} internally.
#'
#' @return A list with components
#' \describe{
#'   \item{beta_est}{A 1 x (p + 1) matrix of regression coefficients followed by
#'                   the residual standard deviation (last element).}
#'   \item{se_est}{A 1 x p matrix of estimated standard errors for the
#'                 regression coefficients.}
#' }
#'
#' @importFrom stats lm coef model.matrix model.response model.frame predict
#' @importFrom survival Surv survreg
#' @importFrom geex m_estimate setup_root_control setup_deriv_control
#' @export
estimate_beta_ipw <- function(data_yXZ, model, model_weights) {

  #######################################################
  ## 1. Build censoring model formula: Surv(W, 1 - D) ~ ...
  #######################################################
  if (!inherits(model_weights, "formula")) {
    stop("'model_weights' must be a formula, e.g. ~ y + Z.")
  }

  # model_weights is like ~ y + Z; extract RHS
  rhs <- deparse(model_weights[[2]])

  cens_formula <- stats::as.formula(
    paste("Surv(W, 1 - D) ~", rhs)
  )

  if (!all(c("W", "D") %in% names(data_yXZ))) {
    stop("data_yXZ must contain columns 'W' and 'D'.")
  }

  #######################################################
  ## 2. Fit Weibull AFT model for C | (Y, Z, ...)
  #######################################################
  wr <- survival::survreg(
    formula = cens_formula,
    data    = data_yXZ,
    dist    = "weibull"
  )

  # Linear predictor for each subject (scale on original survreg scale)
  lp <- stats::predict(wr, newdata = data_yXZ)

  # Shape parameter of Weibull in survreg parameterization
  shape_gamma <- 1 / wr$scale

  # Estimated survival probability at W: S_C(W | covariates)
  data_yXZ$myp <- exp(-(data_yXZ$W / lp)^shape_gamma)

  # Basic sanity check to avoid dividing by 0
  if (any(data_yXZ$myp <= 0)) {
    stop("Some estimated weights are non-positive. Check the censoring model.")
  }

  #######################################################
  ## 3. Define IPW estimating function for geex
  #######################################################
  gee_estfun <- function(data, formula) {

    X   <- stats::model.matrix(object = formula, data = data)
    Y   <- stats::model.response(stats::model.frame(formula = formula, data = data))
    D   <- data$D
    myp <- data$myp

    function(theta) {

      # unweighted score contribution
      mu      <- X %*% theta
      e       <- Y - mu
      dotbeta <- (e / psi^2) %*% X

      # IPW score: sum_i D_i / p_i * S_i
      (D / myp) * c(dotbeta)
    }
  }

  #######################################################
  ## 4. Initial estimates via weighted least squares
  #######################################################
  D   <- data_yXZ$D
  myp <- data_yXZ$myp

  mylmer <- stats::lm(model, data = data_yXZ, weights = D / myp)
  g      <- stats::coef(mylmer)
  psi    <- summary(mylmer)$sigma

  #######################################################
  ## 5. Run M-estimation
  #######################################################
  results <- geex::m_estimate(
    estFUN        = gee_estfun,
    data          = data_yXZ,
    root_control  = geex::setup_root_control(start = c(g)),
    outer_args    = list(formula = model),
    deriv_control = geex::setup_deriv_control(method = "Richardson")
  )

  #######################################################
  ## 6. Return estimates (recompute psi on D = 1 only)
  #######################################################
  beta_estimates <- matrix(t(results@estimates), nrow = 1)

  data_cc <- data_yXZ[data_yXZ$D == 1, , drop = FALSE]
  X_cc    <- stats::model.matrix(object = model, data = data_cc)

  psi_updated <- sqrt(
    sum((data_cc$y - X_cc %*% beta_estimates)^2 / nrow(data_cc))
  )

  beta_estimates <- rbind(beta_estimates, psi_updated)

  se_estimates <- matrix(sqrt(diag(results@vcov)), nrow = 1)

  list(
    beta_est = t(beta_estimates),
    se_est   = t(se_estimates)
  )
}
