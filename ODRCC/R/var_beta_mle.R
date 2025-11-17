#' Sandwich variance estimator for likelihood-based estimator
#'
#' Computes a sandwich variance estimator for the parameter vector \eqn{\theta}
#' obtained from \code{\link{estimate_beta_likelihood_optimx}}. The function
#' allows for a general outcome model \code{y ~ AW + Z1 + ... + Zp} and a
#' possibly different AFT model for \code{X | Z} specified by \code{model_xz}.
#'
#' The parameter vector \eqn{\theta} is assumed to be of the form
#' \eqn{\theta = (\beta, \psi, \gamma_x, \text{shape}_x)}, where:
#' \itemize{
#'   \item \eqn{\beta} are the regression coefficients in the outcome model,
#'   \item \eqn{\psi} is the residual standard deviation on the original scale,
#'   \item \eqn{\gamma_x} indexes the (log-)mean of \eqn{X | Z},
#'   \item \eqn{\text{shape}_x} is the Weibull shape parameter for \eqn{X | Z}.
#' }
#'
#' @param data_yXZ A data frame containing at least the outcome \code{y}, the
#'   auxiliary covariate \code{A}, the observed covariate \code{W}, the event
#'   indicator \code{D}, the exposure covariate \code{AW} (or another name
#'   specified by \code{aw_var}), and any additional covariates appearing in
#'   \code{model} and optionally \code{model_xz}.
#' @param theta Numeric vector of parameter estimates, typically taken from
#'   \code{estimate_beta_mle(...)[["beta_est"]]}. Must be ordered
#'   as \eqn{(\beta, \psi, \gamma_x, \text{shape}_x)}.
#' @param model A \code{\link[stats]{formula}} specifying the outcome regression
#'   model, e.g. \code{y ~ AW + Z1 + Z2}.
#' @param aw_var Character string giving the name of the exposure covariate that
#'   is defined as \code{A - X} (default is \code{"AW"}). This variable must
#'   appear on the right-hand side of \code{model} and as a column in
#'   \code{data_yXZ}.
#' @param model_xz Optional \code{\link[stats]{formula}} specifying the AFT model
#'   for \code{X | Z}. If \code{NULL} (default), the model for \code{X | Z} is
#'   taken to be Weibull with log-mean linear in \code{(1, Z1, ..., Zp)}, where
#'   \code{Z1, ..., Zp} are all covariates on the right-hand side of
#'   \code{model} except \code{aw_var}. If provided, \code{model_xz} can be
#'   either a right-hand-side formula (e.g. \code{~ Z1 + Z2}) or a full Surv
#'   formula (e.g. \code{Surv(W, D) ~ Z1 + Z2}); in both cases, only the RHS is
#'   used here to construct the design for \code{X | Z}.
#'
#' @return A list with components
#' \describe{
#'   \item{beta_est}{Estimated regression coefficients \eqn{\beta}.}
#'   \item{psi_est}{Estimated residual standard deviation \eqn{\psi}.}
#'   \item{se_beta}{Sandwich standard errors for \eqn{\beta}.}
#'   \item{se_psi}{Sandwich standard error for \eqn{\psi}.}
#'   \item{sandwich_var}{Full sandwich variance matrix for \eqn{\theta}.}
#' }
#'
#' @importFrom stats model.matrix model.response model.frame dnorm dweibull integrate terms reformulate
#' @importFrom numDeriv jacobian hessian
#' @export
var_beta_mle <- function(
    data_yXZ,
    theta,
    model,
    aw_var   = "AW",
    model_xz = NULL
) {

  # basic checks
  if (!inherits(model, "formula")) {
    stop("'model' must be a formula, e.g. y ~ AW + Z1 + Z2.")
  }
  if (!aw_var %in% names(data_yXZ)) {
    stop("The variable given by 'aw_var' (", aw_var,
         ") must be a column of data_yXZ.")
  }

  # add subject ID
  data_yXZ$b <- seq_len(nrow(data_yXZ))

  # determine dimension of beta from the outcome model
  X_full <- stats::model.matrix(model, data_yXZ)
  p_beta <- ncol(X_full)

  n_par <- length(theta)
  if (n_par <= p_beta + 2) {
    stop("Length of 'theta' is too small relative to the outcome model. ",
         "Expected at least beta, psi, gamma_x, and shape_x.")
  }

  # theta = (beta, psi, gamma_x (len_gamma), shape_x)
  len_gamma <- n_par - (p_beta + 1) - 1
  if (len_gamma < 1) {
    stop("Could not infer length of gamma_x (parameters for X|Z).")
  }

  # Extract outcome RHS and default Z variables (for default X|Z)
  terms_obj <- stats::terms(model)
  rhs_vars  <- attr(terms_obj, "term.labels")
  if (!aw_var %in% rhs_vars) {
    stop("'", aw_var, "' must appear on the right-hand side of 'model'.")
  }
  z_vars_default <- setdiff(rhs_vars, aw_var)

  ####### Build Z design terms for X|Z (xz_terms) ########
  if (!is.null(model_xz)) {

    if (!inherits(model_xz, "formula")) {
      stop("'model_xz' must be a formula, e.g. ~ Z1 + Z2 or Surv(W, D) ~ Z1 + Z2.")
    }

    # if model_xz has RHS only
    if (length(model_xz) == 2L) {
      xz_rhs   <- deparse(model_xz[[2]])
      xz_terms <- all.vars(stats::reformulate(xz_rhs))
    } else {
      # full formula; extract RHS part
      xz_rhs   <- deparse(model_xz[[3]])
      xz_terms <- all.vars(stats::reformulate(xz_rhs))
    }

  } else {
    # default: same covariates as outcome model except AW
    xz_terms <- z_vars_default
  }

  # function to get design matrix for X|Z at one row (may be intercept-only)
  get_xz_design <- function(data_row) {
    if (length(xz_terms) > 0) {
      stats::model.matrix(
        stats::reformulate(xz_terms),
        data = data_row
      )
    } else {
      matrix(1, nrow = 1, ncol = 1)
    }
  }

  # sanity check: first row design length should match len_gamma
  test_design <- get_xz_design(data_yXZ[1, , drop = FALSE])
  if (ncol(test_design) != len_gamma) {
    stop(
      "Mismatch between length(theta) and X|Z design dimension. ",
      "Expected gamma_x length ", len_gamma, " but got design with ",
      ncol(test_design), " columns. Check 'model_xz' and the way 'theta' ",
      "was constructed."
    )
  }

  upperbound <- Inf

  #########################################################
  # define the per-subject log-likelihood and its derivatives
  #########################################################
  b_mle <- function(data_b, d1 = TRUE) {

    # log-likelihood contribution for one subject
    logLik_b <- function(theta. = theta) {

      # unpack theta
      beta    <- theta.[1:p_beta]
      psi     <- theta.[p_beta + 1]
      gamma_x <- theta.[(p_beta + 2):(p_beta + 1 + len_gamma)]
      shape_x <- theta.[p_beta + 1 + len_gamma + 1]

      # design for outcome
      X <- stats::model.matrix(object = model, data = data_b)
      Y <- stats::model.response(
        stats::model.frame(formula = model, data = data_b)
      )

      D <- data_b$D
      w <- data_b$W

      # design for X|Z
      xz_design <- get_xz_design(data_b)
      meanXZ    <- as.numeric(exp(xz_design %*% gamma_x))

      # delta = 1
      mu <- as.numeric(X %*% beta)
      e  <- Y - mu

      delta1 <- function() {
        dens_y <- stats::dnorm(e, mean = 0, sd = psi)
        dens_x <- stats::dweibull(w, shape = 1 / shape_x, scale = meanXZ)
        log(dens_y * dens_x)
      }

      # delta = 0
      delta0 <- function() {

        # f_Y|X,Z as function of latent t
        likelihood_int <- function(t_vec) {
          sapply(t_vec, function(tt) {
            new_row <- data_b
            # AW = A - X; here X = tt
            new_row[[aw_var]] <- new_row$A - tt
            X_design <- stats::model.matrix(object = model, data = new_row)
            mu_t     <- as.numeric(X_design %*% beta)
            e_t      <- Y - mu_t
            stats::dnorm(e_t, mean = 0, sd = psi)
          })
        }

        # integrand: f_Y|X,Z(t) * f_X|Z(t)
        integral_func <- function(t_vec) {
          likelihood_int(t_vec) *
            stats::dweibull(t_vec, shape = 1 / shape_x, scale = meanXZ)
        }

        val <- stats::integrate(
          f     = integral_func,
          lower = w,
          upper = upperbound
        )$value

        log(val)
      }

      if (D == 1) return(delta1())
      if (D == 0) return(delta0())
      stop("D must be 0 or 1.")
    }

    # first derivative (score) or second derivative (Hessian)
    if (d1) {
      # score: gradient of logLik_b w.r.t. theta at current theta
      myderiv <- numDeriv::jacobian(
        func   = function(x) logLik_b(theta. = x),
        x      = theta,
        method = "Richardson"
      )
      myderiv <- matrix(myderiv, ncol = 1)
    } else {
      # Hessian of logLik_b
      myderiv <- numDeriv::hessian(
        func   = function(x) logLik_b(theta. = x),
        x      = theta,
        method = "Richardson"
      )
    }

    myderiv
  }

  #########################################################
  # B matrix (meat): sum_i s_i s_i^T, where s_i is score for obs i
  #########################################################
  calculate_B <- function() {
    parts <- lapply(seq_len(nrow(data_yXZ)), function(i) {
      s_i <- b_mle(data_b = data_yXZ[i, , drop = FALSE], d1 = TRUE)
      s_i %*% t(s_i)
    })
    Reduce("+", parts)
  }
  B_mat <- calculate_B()

  #########################################################
  # A matrix (bread): sum_i H_i, where H_i is Hessian of loglik for obs i
  #########################################################
  calculate_A <- function() {
    parts <- lapply(seq_len(nrow(data_yXZ)), function(i) {
      b_mle(data_b = data_yXZ[i, , drop = FALSE], d1 = FALSE)
    })
    Reduce("+", parts)
  }
  A_mat <- calculate_A()

  A_inv <- solve(A_mat)

  sand_var <- A_inv %*% B_mat %*% t(A_inv)

  # extract components of interest
  beta_hat <- theta[1:p_beta]
  psi_hat  <- theta[p_beta + 1]

  se_all   <- sqrt(diag(sand_var))
  se_beta  <- se_all[1:p_beta]
  se_psi   <- se_all[p_beta + 1]

  list(
    beta_est     = beta_hat,
    psi_est      = psi_hat,
    se_beta      = se_beta,
    se_psi       = se_psi,
    sandwich_var = sand_var
  )
}
