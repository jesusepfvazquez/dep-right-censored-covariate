#' Augmented IPW (AIPW) estimator for regression parameters
#'
#' Fits a regression model for \code{Y} using an augmented inverse probability
#' weighting (AIPW) estimator in the setting of an outcome-dependent
#' right-censored covariate. The outcome model is specified by \code{model},
#' typically of the form \code{y ~ AW + Z1 + ... + Zp}, where \code{AW = A - X}.
#'
#' The AIPW estimator combines:
#' \itemize{
#'   \item An IPW component based on a censoring model for
#'         \code{C | (Y, Z, ...)} specified by \code{model_weights}.
#'   \item An augmentation component that integrates over the conditional
#'         distribution \code{X | Z} specified by \code{gamma_x} and
#'         \code{model_xz}.
#' }
#'
#' The function solves the AIPW estimating equations for the regression
#' parameters \eqn{\beta}, treating the nuisance models as plug-in, and returns
#' \eqn{\hat \beta} along with a plug-in estimate of the residual standard
#' deviation \eqn{\hat \psi}.
#'
#' @param data_yXZ A data frame containing at least:
#'   \itemize{
#'     \item \code{y}: outcome,
#'     \item \code{A}: auxiliary covariate used to form \code{AW = A - X},
#'     \item \code{W}: observed covariate \code{W = min(X, C)},
#'     \item \code{D}: indicator \code{I(X <= C)},
#'     \item columns for the covariates in \code{model},
#'     \item columns for the covariates in \code{model_weights} and
#'           \code{model_xz}.
#'   }
#' @param model A \code{\link[stats]{formula}} specifying the outcome
#'   regression model, e.g. \code{y ~ AW + Z} or \code{y ~ AW + Z1 + Z2}.
#' @param model_weights A \code{\link[stats]{formula}} specifying the censoring
#'   model for \code{C | (Y, Z, ...)}. Typically a right-hand-side only formula,
#'   such as \code{~ y + Z}, which is internally expanded to
#'   \code{Surv(W, 1 - D) ~ y + Z}.
#' @param model_xz A \code{\link[stats]{formula}} specifying the covariate
#'   structure for \code{X | Z}. Typically right-hand-side only, e.g.
#'   \code{~ Z}, meaning \code{log E[X | Z]} depends on those covariates. Only
#'   the RHS is used.
#' @param aw_var Character string giving the name of the exposure covariate in
#'   \code{model} that equals \code{A - X} (default \code{"AW"}). This must
#'   appear in \code{model} and in \code{data_yXZ}.
#' @param lbound,ubound Numeric lower and upper bounds for the numerical
#'   integration over \code{X} in the augmentation term (defaults: 0 and 50).
#'
#' @return A list with components
#' \describe{
#'   \item{beta_est}{A 1 x p matrix of estimated regression coefficients
#'                   \eqn{\hat \beta}.}
#'   \item{psi_est}{Scalar, estimated residual standard deviation \eqn{\hat \psi}.}
#' }
#'
#' @importFrom stats lm coef model.matrix model.response model.frame dnorm dweibull pweibull integrate
#' @importFrom stats terms reformulate as.formula
#' @importFrom survival Surv survreg
#' @importFrom geex m_estimate setup_root_control setup_deriv_control
#' @export
#'
estimate_beta_aipw <- function(
    data_yXZ,
    model,
    model_weights,
    model_xz,
    aw_var  = "AW",
    lbound  = 0,
    ubound  = 50
) {

  ## ----------------- basic checks & setup ----------------- ##
  if (!inherits(model, "formula")) {
    stop("'model' must be a formula, e.g. y ~ AW + Z.")
  }
  if (!inherits(model_weights, "formula")) {
    stop("'model_weights' must be a formula, e.g. ~ y + Z.")
  }
  if (!inherits(model_xz, "formula")) {
    stop("'model_xz' must be a formula, e.g. ~ Z.")
  }
  if (!aw_var %in% names(data_yXZ)) {
    stop("The variable given by 'aw_var' (", aw_var,
         ") must be a column of data_yXZ.")
  }

  # Outcome design to get p_beta
  X_full <- stats::model.matrix(model, data_yXZ)
  p_beta <- ncol(X_full)

  ## ----------------- 0. Event model X | Z ----------------- ##
  # Build censoring formula: Surv(W, 1 - D) ~ ...
  if (length(model_xz) == 2L) {
    rhs_w        <- deparse(model_xz[[2]])
    event_formula <- stats::as.formula(paste("Surv(W, D) ~", rhs_w))
  } else {
    rhs_w        <- deparse(model_xz[[3]])
    event_formula <- stats::as.formula(paste("Surv(W, D) ~", rhs_w))
  }

  wr_xz <- survival::survreg(
    formula = event_formula,
    data    = data_yXZ,
    dist    = "weibull"
  )

  gamma_coef  <- as.numeric(coef(wr_xz))
  shape.x     <- 1 / wr_xz$scale

  ## ----------------- 1. Censoring model C | (Y, Z, ...) ----------------- ##
  # Build censoring formula: Surv(W, 1 - D) ~ ...
  if (length(model_weights) == 2L) {
    rhs_w        <- deparse(model_weights[[2]])
    cens_formula <- stats::as.formula(paste("Surv(W, 1 - D) ~", rhs_w))
  } else {
    rhs_w        <- deparse(model_weights[[3]])
    cens_formula <- stats::as.formula(paste("Surv(W, 1 - D) ~", rhs_w))
  }

  wr <- survival::survreg(
    formula = cens_formula,
    data    = data_yXZ,
    dist    = "weibull"
  )

  alpha_coef  <- as.numeric(coef(wr))
  shape.c <- 1 / wr$scale

  # meanCYZ: on original AFT scale, used in augmentation
  data_yXZ$meanCYZ <- exp(wr$linear.predictors)

  # IPW weights: S_C(W | covariates) under Weibull(AFT)
  lp_c  <- exp(wr$linear.predictors)  # scale parameter
  data_yXZ$myp <- exp(-(data_yXZ$W / lp_c)^shape.c)

  ## ----------------- 2. X | Z model: meanXZ ----------------- ##
  # Extract RHS terms for model_xz
  if (length(model_xz) == 2L) {
    xz_rhs   <- deparse(model_xz[[2]])
    xz_terms <- all.vars(stats::reformulate(xz_rhs))
  } else {
    xz_rhs   <- deparse(model_xz[[3]])
    xz_terms <- all.vars(stats::reformulate(xz_rhs))
  }

  get_xz_design <- function(data_row) {
    if (length(xz_terms) > 0) {
      stats::model.matrix(stats::reformulate(xz_terms), data = data_row)
    } else {
      matrix(1, nrow = 1, ncol = 1)
    }
  }

  # check gamma length vs design
  test_design <- get_xz_design(data_yXZ[1, , drop = FALSE])
  if (ncol(test_design) != length(gamma_coef)) {
    stop("Length of gamma_coef (", length(gamma_coef),
         ") does not match design columns for X|Z (", ncol(test_design), ").")
  }

  # precompute meanXZ for each row
  data_yXZ$meanXZ <- NA_real_
  for (i in seq_len(nrow(data_yXZ))) {
    row_i <- data_yXZ[i, , drop = FALSE]
    xz_i  <- get_xz_design(row_i)
    data_yXZ$meanXZ[i] <- as.numeric(exp(xz_i %*% gamma_coef))
  }

  ## ----------------- 3. Initial estimates for beta, psi ----------------- ##
  D   <- data_yXZ$D
  myp <- data_yXZ$myp

  init_fit <- stats::lm(model, data = data_yXZ, weights = D / myp)
  g_init   <- stats::coef(init_fit)
  psi      <- summary(init_fit)$sigma

  ## ----------------- 4. AIPW estimating function for geex ----------------- ##

  # per-observation AIPW score phi_i(beta)
  aipw_score_i <- function(beta., data_row) {

    # outcome design at observed AW
    X_row <- stats::model.matrix(model, data_row)
    Y     <- stats::model.response(
      stats::model.frame(formula = model, data = data_row)
    )
    D_i   <- data_row$D
    myp_i <- data_row$myp

    meanXZ_i  <- data_row$meanXZ
    meanCYZ_i <- data_row$meanCYZ

    # complete-case score (for beta)
    mu <- as.numeric(X_row %*% beta.)
    e  <- Y - mu
    cc_score <- as.numeric((e / psi^2) * X_row)  # length p_beta

    # helper: f_Y|X,Z at t
    likelihood_int <- function(t) {
      sapply(t, function(tt) {
        new_row <- data_row
        # AW = A - X; here X = tt
        new_row[[aw_var]] <- new_row$A - tt
        X_t <- stats::model.matrix(model, new_row)
        mu_t <- as.numeric(X_t %*% beta.)
        e_t  <- Y - mu_t
        stats::dnorm(e_t, mean = 0, sd = psi)
      })
    }

    # helper: score wrt beta_j at t
    score_int <- function(t, j) {
      sapply(t, function(tt) {
        new_row <- data_row
        new_row[[aw_var]] <- new_row$A - tt
        X_t <- stats::model.matrix(model, new_row)
        mu_t <- as.numeric(X_t %*% beta.)
        e_t  <- Y - mu_t
        dotbeta_t <- (e_t / psi^2) * X_t
        dotbeta_t[1, j]
      })
    }

    shape_c_i <- shape.c
    shape_x_i <- shape.x

    # numerator: vector length p_beta
    numerator_j <- function(j) {
      f_num <- function(t) {
        likelihood_int(t) *
          score_int(t, j) *
          stats::dweibull(t, shape = shape_x_i, scale = meanXZ_i) *
          (1 - 1 / stats::pweibull(
            t,
            shape = shape_c_i,
            scale = meanCYZ_i,
            lower.tail = FALSE
          ))
      }
      stats::integrate(f_num, lower = lbound, upper = ubound)$value
    }
    num_vec <- vapply(seq_len(p_beta), numerator_j, numeric(1))

    # denominator: scalar
    denom_fun <- function(t) {
      likelihood_int(t) *
        stats::dweibull(t, shape = shape_x_i, scale = meanXZ_i) *
        (1 - 1 / stats::pweibull(
          t,
          shape = shape_c_i,
          scale = meanCYZ_i,
          lower.tail = FALSE
        ))
    }
    denom_val <- stats::integrate(denom_fun, lower = lbound, upper = ubound)$value

    aug_vec <- num_vec / denom_val

    # AIPW estimating function for this observation
    (D_i / myp_i) * cc_score + (1 - D_i / myp_i) * aug_vec
  }

  # geex estimating function: sum_i phi_i(beta)
  gee_estfun <- function(data, formula) {

    function(theta) {
      # theta is beta vector here
      score_mat <- matrix(0, nrow = nrow(data), ncol = p_beta)
      for (i in seq_len(nrow(data))) {
        score_mat[i, ] <- aipw_score_i(theta, data[i, , drop = FALSE])
      }
      # return vector of estimating equations (sum over i)
      colSums(score_mat)
    }
  }

  ## ----------------- 5. Solve estimating equations via geex ----------------- ##

  results <- geex::m_estimate(
    estFUN        = gee_estfun,
    data          = data_yXZ,
    root_control  = geex::setup_root_control(start = g_init),
    outer_args    = list(formula = model),
    deriv_control = geex::setup_deriv_control(method = "Richardson")
  )

  beta_hat <- as.numeric(results@estimates)
  se_estimates <- matrix(sqrt(diag(results@vcov)), nrow = 1)

  ## ----------------- 6. Re-estimate psi using-beta and D = 1 ----------------- ##
  data_cc <- data_yXZ[data_yXZ$D == 1, , drop = FALSE]
  X_cc    <- stats::model.matrix(model, data_cc)
  psi_hat <- sqrt(
    sum((data_cc$y - as.numeric(X_cc %*% beta_hat))^2 / nrow(data_cc))
  )

  # return
  list(
    beta_est = matrix(c(beta_hat,psi_hat), nrow = 1),
    beta_se = matrix(se_estimates, nrow = 1)
  )
}
