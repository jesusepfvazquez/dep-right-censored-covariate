#' Sandwich variance estimator for the AIPW estimator
#'
#' Computes a plug-in sandwich variance estimator for the AIPW regression
#' estimator under outcome-dependent right-censoring of a covariate. The
#' outcome model is specified by \code{model}, typically of the form
#' \code{y ~ AW + Z1 + ... + Zp}, where \code{AW = A - X}.
#'
#' The AIPW estimator combines:
#' \itemize{
#'   \item an IPW component based on a censoring model for
#'         \code{C | (Y, Z, ...)} specified by \code{model_weights}, and
#'   \item an augmentation component based on an AFT model for
#'         \code{X | Z} specified by \code{gamma_x} and \code{model_xz}.
#' }
#'
#' This function takes the estimated parameter vector \code{theta = c(beta, psi)}
#' from \code{\link{estimate_beta_aipw_est}} and computes a sandwich variance
#' for \eqn{\beta}, treating the nuisance model parameters as plug-in.
#'
#' @param data_yXZ Data frame containing at least:
#'   \itemize{
#'     \item \code{y}: outcome,
#'     \item \code{A}: auxiliary covariate used to form \code{AW = A - X},
#'     \item \code{W}: observed covariate \code{W = min(X, C)},
#'     \item \code{D}: indicator \code{I(X \le C)},
#'     \item all covariates appearing in \code{model},
#'     \item all covariates appearing in \code{model_weights} and \code{model_xz}.
#'   }
#' @param theta Numeric vector \code{c(beta, psi)} from the AIPW estimator,
#'   where \code{beta} has length equal to the number of columns in
#'   \code{model.matrix(model, data_yXZ)} and \code{psi} is the residual
#'   standard deviation.
#' @param model \code{\link[stats]{formula}} for the outcome, e.g.
#'   \code{y ~ AW + Z} or \code{y ~ AW + Z1 + Z2}.
#' @param model_weights \code{\link[stats]{formula}} specifying the censoring
#'   model for \code{C | (Y, Z, ...)}. Typically a RHS-only formula such as
#'   \code{~ y + Z}, which is internally expanded to
#'   \code{Surv(W, 1 - D) ~ y + Z}.
#' @param gamma_x Numeric vector of parameters for the \code{X | Z} AFT model:
#'   \code{gamma_x = c(gamma_coef, shape_par)}, where \code{gamma_coef} indexes
#'   \code{log-mean(X | Z)} and \code{shape_par} is the AFT scale parameter
#'   used via \code{shape.x = 1 / shape_par}.
#' @param model_xz \code{\link[stats]{formula}} for the covariate structure
#'   of \code{X | Z}, typically RHS-only such as \code{~ Z}. Only the RHS is
#'   used.
#' @param aw_var Character, name of the exposure covariate in \code{model} that
#'   equals \code{A - X} (default \code{"AW"}).
#' @param lbound,ubound Numeric lower and upper bounds for the numerical
#'   integration over \code{X} in the augmentation term (defaults: 0, 50).
#'
#' @return A list with components
#' \describe{
#'   \item{beta_est}{Estimated regression coefficients \eqn{\beta}.}
#'   \item{psi_est}{Estimated residual standard deviation \eqn{\psi}.}
#'   \item{se_beta}{Sandwich standard errors for \eqn{\beta}.}
#'   \item{sandwich_var}{Full sandwich variance matrix for \eqn{\beta}.}
#' }
#'
#' @importFrom stats model.matrix model.response model.frame dnorm dweibull pweibull integrate
#' @importFrom stats terms reformulate as.formula
#' @importFrom survival Surv survreg
#' @importFrom numDeriv jacobian
#' @export
var_beta_aipw_est <- function(
    data_yXZ,
    theta,
    model,
    model_weights,
    gamma_x,
    model_xz,
    aw_var  = "AW",
    lbound  = 0,
    ubound  = 50
) {

  ## ----------------- basic checks & unpack theta ----------------- ##
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

  X_full <- stats::model.matrix(model, data_yXZ)
  p_beta <- ncol(X_full)

  if (length(theta) != p_beta + 1) {
    stop("Length(theta) must be p_beta + 1 (beta, psi). Here p_beta = ",
         p_beta, ", but length(theta) = ", length(theta))
  }

  beta <- theta[1:p_beta]
  psi  <- theta[p_beta + 1]

  if (length(gamma_x) < 2) {
    stop("'gamma_x' must contain at least coefficients and a final shape parameter.")
  }
  gamma_coef      <- gamma_x[1:(length(gamma_x) - 1)]
  gamma_shape_par <- gamma_x[length(gamma_x)]
  shape.x         <- 1 / gamma_shape_par

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

  shape.c <- 1 / wr$scale

  # meanCYZ: AFT scale on original scale, used in augmentation + survival
  data_yXZ$meanCYZ <- exp(wr$linear.predictors)

  # IPW weights: S_C(W | covariates) under Weibull(AFT)
  lp_c  <- exp(wr$linear.predictors)  # scale
  data_yXZ$myp <- exp(-(data_yXZ$W / lp_c)^shape.c)

  ## ----------------- 2. X | Z model: meanXZ ----------------- ##
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

  ## ----------------- 3. AIPW score for a single observation ----------------- ##

  aipw_score_i <- function(beta., data_row) {

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

  ## ----------------- 4. Sandwich A and B for β ----------------- ##

  # A = sum_i d phi_i / d beta^T
  calculate_A <- function() {
    parts <- lapply(seq_len(nrow(data_yXZ)), function(i) {
      numDeriv::jacobian(
        func   = function(b) aipw_score_i(b, data_yXZ[i, , drop = FALSE]),
        x      = beta,
        method = "Richardson"
      )
    })
    Reduce("+", parts)
  }

  # B = sum_i phi_i phi_i^T
  calculate_B <- function() {
    parts <- lapply(seq_len(nrow(data_yXZ)), function(i) {
      s_i <- aipw_score_i(beta, data_yXZ[i, , drop = FALSE])
      s_i %*% t(s_i)
    })
    Reduce("+", parts)
  }

  A_mat <- calculate_A()
  B_mat <- calculate_B()
  A_inv <- solve(A_mat)

  sand_var <- A_inv %*% B_mat %*% t(A_inv)

  se_beta <- sqrt(diag(sand_var))

  list(
    beta_est     = beta,
    psi_est      = psi,
    se_beta      = se_beta,
    sandwich_var = sand_var
  )
}
