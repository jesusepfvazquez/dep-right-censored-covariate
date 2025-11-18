#' Sandwich variance estimator for the IPW estimator
#'
#' Computes a sandwich variance estimator for the IPW regression estimator,
#' allowing for a general outcome model \code{y ~ AW + Z1 + ... + Zp} and a
#' possibly different censoring model for \code{C | (Y, Z, ...)} specified by
#' \code{model_weights}. The censoring model is fit via a Weibull AFT model
#' using \code{survreg}, and the corresponding Gumbel parameterization for
#' \code{log(W)} is used to construct the weights \eqn{\pi(Y, W, Z)}.
#'
#' The parameter vector \code{theta} is assumed to be of the form
#' \eqn{\theta = (\beta, \psi)}, where:
#' \itemize{
#'   \item \eqn{\beta} are the regression coefficients in the outcome model,
#'   \item \eqn{\psi} is the residual standard deviation on the original scale.
#' }
#'
#' @param data_yXZ A data frame containing at least the outcome \code{y}, the
#'   covariates in the outcome model \code{model}, the observed covariate
#'   \code{W}, the event indicator \code{D}, and the covariates appearing in the
#'   censoring model \code{model_weights}.
#' @param theta Numeric vector of parameter estimates \eqn{(\beta, \psi)} from
#'   the IPW estimator. The length of \code{theta} must be equal to
#'   \code{p_beta + 1}, where \code{p_beta} is the number of regression
#'   coefficients in \code{model}.
#' @param model A \code{\link[stats]{formula}} specifying the outcome regression
#'   model, e.g. \code{y ~ AW + Z1 + Z2}.
#' @param model_weights A \code{\link[stats]{formula}} specifying the censoring
#'   model used to estimate the IPW weights. Typically of the form
#'   \code{~ y + Z1 + Z2} (right-hand-side only). Internally this is expanded to
#'   \code{Surv(W, 1 - D) ~ y + Z1 + Z2}. If a full \code{Surv} formula is
#'   provided, only its right-hand side is used.
#'
#' @return A list with components
#' \describe{
#'   \item{beta_est}{Estimated regression coefficients \eqn{\beta}.}
#'   \item{psi_est}{Estimated residual standard deviation \eqn{\psi}.}
#'   \item{se_beta}{Sandwich standard errors for \eqn{\beta}.}
#'   \item{sandwich_var}{Full sandwich variance matrix for the stacked nuisance
#'                       parameter vector \eqn{\xi = (\beta, \alpha)}, where
#'                       \eqn{\alpha} are the parameters of the censoring model.}
#' }
#'
#' @importFrom stats model.matrix model.response model.frame terms reformulate
#' @importFrom survival Surv survreg
#' @importFrom numDeriv jacobian hessian
#' @export
var_beta_ipw <- function(
    data_yXZ,
    theta,
    model,
    model_weights
) {

  ## basic checks
  if (!inherits(model, "formula")) {
    stop("'model' must be a formula, e.g. y ~ AW + Z1 + Z2.")
  }
  if (!inherits(model_weights, "formula")) {
    stop("'model_weights' must be a formula, e.g. ~ y + Z or Surv(W, 1 - D) ~ y + Z.")
  }

  # Outcome design: to infer p_beta
  X_full <- stats::model.matrix(model, data_yXZ)
  p_beta <- ncol(X_full)

  if (length(theta) != p_beta + 1) {
    stop("Length of 'theta' must be p_beta + 1 (beta + psi). Here p_beta = ",
         p_beta, ", but length(theta) = ", length(theta), ".")
  }

  mybeta <- theta[1:p_beta]
  mypsi  <- theta[p_beta + 1]

  # Add subject ID if useful later
  data_yXZ$b <- seq_len(nrow(data_yXZ))

  #########################################################
  ## 1. Fit censoring model Surv(W, 1 - D) ~ covariates
  #########################################################
  # Build censoring formula
  if (length(model_weights) == 2L) {
    # RHS-only: ~ y + Z1 + Z2
    rhs_w   <- deparse(model_weights[[2]])
    cens_formula <- stats::as.formula(paste("Surv(W, 1 - D) ~", rhs_w))
  } else {
    # Full formula: Surv(W, 1 - D) ~ ...
    rhs_w   <- deparse(model_weights[[3]])
    cens_formula <- stats::as.formula(paste("Surv(W, 1 - D) ~", rhs_w))
  }

  wr <- survival::survreg(
    formula = cens_formula,
    data    = data_yXZ,
    dist    = "weibull"
  )

  # alpha = (gamma, shape), using Gumbel parametrization for log(W)
  myalpha <- c(stats::coef(wr), wr$scale)
  len_alpha <- length(myalpha)

  # extract terms used in censoring model (RHS only)
  weight_terms <- attr(wr$terms, "term.labels")

  # helper to build design matrix for weights model
  get_weight_design <- function(data.) {
    if (length(weight_terms) > 0) {
      stats::model.matrix(
        stats::reformulate(weight_terms),
        data = data.
      )
    } else {
      matrix(1, nrow = nrow(data.), ncol = 1)
    }
  }

  #########################################################
  ## 2. Gumbel helper functions for log(W)
  #########################################################
  dgumbel <- function(x, shape, scale) {
    # density of Gumbel(location=scale, scale=shape)
    (1 / shape) * exp((x - scale) / shape) * exp(-exp((x - scale) / shape))
  }

  pgumbel <- function(x, shape, scale) {
    # CDF of Gumbel(location=scale, scale=shape)
    exp(-exp((x - scale) / shape))
  }

  # pi(Y, W, Z) = P(C >= W | covariates) in Gumbel param
  pi_xz_func_b <- function(data., myalpha.) {

    Wlog  <- log(data.$W)

    gamma_vec <- myalpha.[1:(length(myalpha.) - 1)]
    shape_x   <- myalpha.[length(myalpha.)]

    Xw       <- get_weight_design(data.)
    linPred  <- as.numeric(Xw %*% gamma_vec)

    pgumbel(Wlog, shape = shape_x, scale = linPred)
  }

  # log-likelihood contribution for alpha (censoring model) for one or more rows
  alpha_logLik_b <- function(data., myalpha.) {

    Wlog  <- log(data.$W)

    gamma_vec <- myalpha.[1:(length(myalpha.) - 1)]
    shape_x   <- myalpha.[length(myalpha.)]

    Xw       <- get_weight_design(data.)
    linPred  <- as.numeric(Xw %*% gamma_vec)

    val <- (1 - data.$D) * dgumbel(Wlog, shape_x, linPred) +
      data.$D * pgumbel(Wlog, shape_x, linPred)

    log(val)
  }

  #########################################################
  ## 3. Stacked score s_i(xi) = (score_beta_i, score_alpha_i)
  #########################################################
  # xi = (beta, alpha)
  myxi <- c(mybeta, myalpha)

  ipwscore_b <- function(myxi., data.b) {

    beta_len <- p_beta
    beta.    <- myxi.[1:beta_len]
    alpha.   <- myxi.[(beta_len + 1):length(myxi.)]

    # outcome pieces
    X <- stats::model.matrix(object = model, data = data.b)
    Y <- stats::model.response(
      stats::model.frame(formula = model, data = data.b)
    )
    D <- data.b$D

    # weights
    myp <- pi_xz_func_b(data. = data.b, myalpha. = alpha.)

    # outcome score (beta)
    mu <- as.numeric(X %*% beta.)
    e  <- Y - mu
    dotbeta <- as.numeric((e / mypsi^2) %*% X) * (D / myp)

    # alpha score: gradient of alpha log-likelihood
    dotalpha <- numDeriv::jacobian(
      func   = function(x) alpha_logLik_b(data. = data.b, myalpha. = x),
      x      = alpha.,
      method = "Richardson"
    )
    dotalpha <- as.numeric(dotalpha)

    c(dotbeta, dotalpha)
  }

  #########################################################
  ## 4. Sandwich matrices A (bread) and B (meat)
  #########################################################
  # A = sum_i d s_i / d xi^T
  calculate_A <- function() {
    parts <- lapply(seq_len(nrow(data_yXZ)), function(i) {
      numDeriv::jacobian(
        func   = function(x) ipwscore_b(myxi. = x,
                                        data.b = data_yXZ[i, , drop = FALSE]),
        x      = myxi,
        method = "Richardson"
      )
    })
    Reduce("+", parts)
  }
  A_mat <- calculate_A()
  A_inv <- solve(A_mat)
  message("IPW sandwich: finished bread (A matrix).")

  # B = sum_i s_i s_i^T
  calculate_B <- function() {
    parts <- lapply(seq_len(nrow(data_yXZ)), function(i) {
      s_i <- ipwscore_b(
        myxi.   = myxi,
        data.b  = data_yXZ[i, , drop = FALSE]
      )
      s_i %*% t(s_i)
    })
    Reduce("+", parts)
  }
  B_mat <- calculate_B()
  message("IPW sandwich: finished meat (B matrix).")

  sand_var <- A_inv %*% B_mat %*% t(A_inv)

  # Extract SE for beta (first p_beta components of xi)
  se_all  <- sqrt(diag(sand_var))
  se_beta <- se_all[1:p_beta]

  return(
    list(
      beta_est     = mybeta,
      psi_est      = mypsi,
      se_est      = se_beta
      # sandwich_var = sand_var
    )
  )
}
