#' Likelihood-based estimator via optimx for general AW + Z model
#'
#' Fits a parametric likelihood model for the outcome and the right-censored
#' covariate using direct maximization of the log-likelihood with
#' \code{optimx}. The outcome model is specified by \code{model}, typically of
#' the form \code{y ~ AW + Z1 + ... + Zp}, where \code{AW = A - X}. The
#' distribution of \code{X | Z} is modeled via a Weibull AFT model, whose
#' covariate structure can be specified by \code{model_xz} or, by default,
#' derived from the right-hand side of \code{model} by excluding \code{AW}.
#'
#' @param data_yXZ A data frame containing at least the outcome \code{y}, the
#'   auxiliary covariate \code{A}, the observed covariate \code{W}, the event
#'   indicator \code{D}, the exposure covariate \code{AW} (or another name
#'   specified by \code{aw_var}), and any additional covariates appearing in
#'   \code{model} and optionally \code{model_xz}.
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
#'   \code{model} except \code{aw_var}. If provided, \code{model_xz} can be either
#'   a full Surv formula such as \code{Surv(W, D) ~ Z1 + Z2} or a right-hand-side
#'   formula such as \code{~ Z1 + Z2}, in which case the left-hand side
#'   \code{Surv(W, D)} is filled in automatically.
#'
#' @return A list with component
#' \describe{
#'   \item{beta_est}{A numeric vector containing the parameter vector
#'   \eqn{\theta = (\beta, \psi, \gamma_x, \text{shape}_x)}, where
#'   \eqn{\beta} are the outcome regression coefficients (dimension determined
#'   by \code{model}), \eqn{\psi} is the residual standard deviation on the
#'   original scale, \eqn{\gamma_x} indexes the mean of \eqn{X | Z}, and
#'   \eqn{\text{shape}_x} is the Weibull shape parameter for \eqn{X | Z}.}
#' }
#'
#' @importFrom stats lm coef model.matrix model.response model.frame dnorm dweibull integrate terms reformulate
#' @importFrom survival Surv survreg
#' @importFrom optimx optimx coef<-
#' @export
estimate_beta_likelihood_optimx <- function(
    data_yXZ,
    model,
    aw_var   = "AW",
    model_xz = NULL
) {

  ## basic checks
  if (!inherits(model, "formula")) {
    stop("'model' must be a formula, e.g. y ~ AW + Z1 + Z2.")
  }
  if (!aw_var %in% names(data_yXZ)) {
    stop("The variable given by 'aw_var' (", aw_var,
         ") must be a column of data_yXZ.")
  }

  # extract RHS labels from model
  terms_obj <- stats::terms(model)
  rhs_vars  <- attr(terms_obj, "term.labels")

  if (!aw_var %in% rhs_vars) {
    stop("'", aw_var, "' must appear on the right-hand side of 'model'.")
  }

  # default Z covariates: all RHS vars except AW
  z_vars_default <- setdiff(rhs_vars, aw_var)

  # add subject ID
  data_yXZ$b <- seq_len(nrow(data_yXZ))

  ######## X | Z via Weibull AFT (for starting values) ########
  # Construct the Surv(W, D) ~ ... formula
  if (!is.null(model_xz)) {

    if (!inherits(model_xz, "formula")) {
      stop("'model_xz' must be a formula, e.g. ~ Z1 + Z2 or Surv(W, D) ~ Z1 + Z2.")
    }

    # if model_xz has RHS only (~ Z1 + Z2)
    if (length(model_xz) == 2L) {
      rhs_xz   <- deparse(model_xz[[2]])
      x_formula <- stats::as.formula(paste("Surv(W, D) ~", rhs_xz))
    } else {
      # full formula; use only RHS but ensure LHS is Surv(W, D)
      rhs_xz   <- deparse(model_xz[[3]])
      x_formula <- stats::as.formula(paste("Surv(W, D) ~", rhs_xz))
    }

  } else {
    # default: Surv(W, D) ~ Z1 + ... + Zp or ~ 1 if no Z
    if (length(z_vars_default) > 0) {
      z_rhs     <- paste(z_vars_default, collapse = " + ")
      x_formula <- stats::as.formula(paste("Surv(W, D) ~", z_rhs))
    } else {
      x_formula <- survival::Surv(W, D) ~ 1
    }
  }

  # Fit Weibull AFT model for X | Z
  wr <- survival::survreg(
    formula = x_formula,
    data    = data_yXZ,
    dist    = "weibull"
  )

  # gamma_x: regression coefficients for X|Z;
  # shape_x: scale parameter (we'll reparametrize via shape_x)
  mygamma <- c(stats::coef(wr), wr$scale)

  # Terms used in X|Z mean; will be used to build design rows
  xz_terms <- attr(wr$terms, "term.labels")  # may be character(0) if intercept-only

  ######## Initial parameter estimates for outcome model ########
  # Use complete-case (D == 1) regression for starting values
  data_cc <- data_yXZ[data_yXZ$D == 1, , drop = FALSE]
  mylmer  <- stats::lm(model, data = data_cc)

  g_beta <- stats::coef(mylmer)                   # length p_beta
  g_psi  <- log(summary(mylmer)$sigma)            # log(psi) as working scale

  p_beta  <- length(g_beta)
  len_gamma <- length(mygamma) - 1  # coefficients for X|Z (excluding scale)

  # Full initial parameter vector:
  # theta = (beta, log(psi), gamma_x (len_gamma), shape_x)
  g <- c(g_beta, g_psi, mygamma)

  ## set up upper bound for integration over X when D = 0
  upperbound <- Inf

  ########  DEFINE LOG-LIKELIHOOD FUNCTION ##########
  loglik_fn <- function(theta = g) {

    mydata <- data_yXZ

    ## unpack theta
    beta    <- theta[1:p_beta]
    psi     <- exp(theta[p_beta + 1])
    gamma_x <- theta[(p_beta + 2):(p_beta + 1 + len_gamma)]
    shape_x <- theta[p_beta + 1 + len_gamma + 1]

    Lik.b <- function(data_row) {

      # data_row is a 1-row data frame
      Xmat <- stats::model.matrix(object = model, data = data_row)
      Y    <- stats::model.response(
        stats::model.frame(formula = model, data = data_row)
      )

      D <- data_row$D
      w <- data_row$W

      # Design for X|Z: intercept + covariates from xz_terms
      if (length(xz_terms) > 0) {
        xz_design <- stats::model.matrix(
          stats::reformulate(xz_terms),
          data = data_row
        )
      } else {
        xz_design <- matrix(1, nrow = 1, ncol = 1)  # intercept-only
      }

      meanXZ <- as.numeric(exp(xz_design %*% gamma_x))

      ## contribution when D = 1
      mu <- as.numeric(Xmat %*% beta)
      e  <- Y - mu

      delta1 <- function() {
        dens_y <- stats::dnorm(e, mean = 0, sd = psi)
        dens_x <- stats::dweibull(w, shape = 1 / shape_x, scale = meanXZ)
        log(dens_y * dens_x)
      }

      ## contribution when D = 0
      delta0 <- function() {

        # f_Y|X,Z as a function of latent t
        likelihood_int <- function(t_vec) {
          sapply(t_vec, function(tt) {
            new_row <- data_row
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

    # sum over subjects
    ids <- unique(mydata$b)
    loglik_v <- vapply(ids, function(id) {
      Lik.b(mydata[mydata$b == id, , drop = FALSE])
    }, numeric(1))

    sum(loglik_v)
  }

  ######## Optimize log-likelihood with optimx ########
  opt_result <- optimx::optimx(
    par     = g,
    fn      = loglik_fn,
    method  = "Nelder-Mead",
    control = list(
      trace    = 0,
      maximize = TRUE,  # maximize log-likelihood
      abstol   = 1e-8,
      reltol   = 1e-8
    ),
    hessian = FALSE
  )

  ######## RETURN ESTIMATES ##########
  theta_hat <- as.numeric(optimx::coef(opt_result)[1, ])

  # transform psi back to original scale (component p_beta + 1 is log(psi))
  theta_hat[p_beta + 1] <- exp(theta_hat[p_beta + 1])

  list(beta_est = theta_hat)
}
