#' AIPW lambda-close estimator
#'
#' Fits a regression model for \code{y} using an augmented IPW estimator
#' with a closed-form ("lambda-close") augmentation under a normal
#' approximation for \code{X}. This implementation assumes the outcome
#' model is \code{y ~ AW + Z}, where \code{AW = A - X}, and uses a
#' single covariate \code{Z} in the outcome model.
#'
#' The estimator combines:
#' \itemize{
#'   \item an IPW component using weights \code{D / p}, where \code{p}
#'         is chosen via \code{myweight}, and
#'   \item a closed-form augmentation term derived under a normal
#'         approximation for \code{X}.
#' }
#'
#' @param data_yXZ Data frame containing at least
#'   \code{y, A, AW, W, D, Z}, and the weight variables
#'   \code{myp_ywz_oracle}, \code{myp_ywz}, \code{myp_uniform},
#'   \code{myp_ywz_logit}. Assumes the outcome model is \code{y ~ AW + Z}.
#' @param model Formula for the outcome, currently assumed to be
#'   \code{y ~ AW + Z}. The function checks that the design matrix has
#'   exactly three columns: intercept, AW, and Z, in that order.
#' @param myweight Character; which weight model to use:
#'   \code{"oracle"}, \code{"aft"}, \code{"uniform"}, or \code{"logit"}.
#'
#' @return A list with components:
#' \describe{
#'   \item{beta_est}{1 x 4 matrix: \eqn{(\hat\beta_0, \hat\beta_1, \hat\beta_2, \hat\psi)}.}
#'   \item{se_est}{1 x 3 matrix of standard errors for \eqn{\beta_0, \beta_1, \beta_2}
#'                 based on the geex sandwich estimator.}
#' }
#'
#' @importFrom stats lm coef model.matrix model.response var
#' @importFrom geex m_estimate setup_root_control setup_deriv_control
#' @export
estimate_beta_aipw_lambda_close <- function(
    data_yXZ,
    model,
    myweight = "oracle"
) {

  ## ----- basic checks -----
  if (!inherits(model, "formula")) {
    stop("'model' must be a formula, e.g. y ~ AW + Z.")
  }

  # outcome design to check structure
  X_check <- stats::model.matrix(model, data_yXZ)
  if (ncol(X_check) != 3L) {
    stop(
      "This lambda-close implementation assumes model with 3 parameters:\n",
      "  y ~ AW + Z (intercept, AW, Z). Currently got ",
      ncol(X_check), " columns in model.matrix."
    )
  }

  ## ------------------------------------------------------
  ## 1. Define parameters for X using complete cases
  ## ------------------------------------------------------
  # approximate distribution of X by N(meanXZ, sdXZ^2)
  idx_cc <- which(data_yXZ$D == 1)
  if (length(idx_cc) < 2L) {
    stop("Not enough complete cases (D == 1) to estimate sdXZ and meanXZ.")
  }
  sdXZ   <- sqrt(stats::var(data_yXZ$W[idx_cc]))
  meanXZ <- mean(data_yXZ$W[idx_cc])

  ## ------------------------------------------------------
  ## 2. Choose weight p_i
  ## ------------------------------------------------------
  if (!("myp_ywz_oracle" %in% names(data_yXZ))) {
    stop("data_yXZ must contain 'myp_ywz_oracle' for myweight = 'oracle'.")
  }
  data_yXZ$myp <- data_yXZ$myp_ywz_oracle
  if (myweight == "aft") {
    if (!("myp_ywz" %in% names(data_yXZ))) {
      stop("data_yXZ must contain 'myp_ywz' for myweight = 'aft'.")
    }
    data_yXZ$myp <- data_yXZ$myp_ywz
  }
  if (myweight == "uniform") {
    if (!("myp_uniform" %in% names(data_yXZ))) {
      stop("data_yXZ must contain 'myp_uniform' for myweight = 'uniform'.")
    }
    data_yXZ$myp <- data_yXZ$myp_uniform
  }
  if (myweight == "logit") {
    if (!("myp_ywz_logit" %in% names(data_yXZ))) {
      stop("data_yXZ must contain 'myp_ywz_logit' for myweight = 'logit'.")
    }
    data_yXZ$myp <- data_yXZ$myp_ywz_logit
  }

  ## ------------------------------------------------------
  ## 3. Initial estimates via IPW
  ## ------------------------------------------------------
  D   <- data_yXZ$D
  myp <- data_yXZ$myp

  ipw_fit <- stats::lm(model, data = data_yXZ, weights = D / myp)
  g       <- stats::coef(ipw_fit)  # (beta0, beta1, beta2)
  psi     <- summary(ipw_fit)$sigma

  ## ------------------------------------------------------
  ## 4. Compute Lambda (A in your notation) at theta = g
  ## ------------------------------------------------------
  calculate_Lambda <- function() {

    theta_loc <- g  # length 3

    # per-observation contribution
    calculate_A_b <- function(data_row, part1 = TRUE) {

      W <- stats::model.matrix(~ A + Z, data = data_row)
      X <- stats::model.matrix(y ~ AW + Z, data = data_row)
      Y <- stats::model.response(
        stats::model.frame(formula = y ~ AW + Z, data = data_row)
      )
      D_i   <- data_row$D
      myp_i <- data_row$myp

      ## complete-case score
      CCscore <- function() {
        mu <- as.numeric(X %*% theta_loc)
        e  <- Y - mu
        as.numeric((e / psi^2) * X)  # length 3
      }

      ## closed-form augmentation psi_hat_i(theta_loc)
      psi_hat_i <- function() {

        e_star   <- Y - as.numeric(W %*% theta_loc)
        a        <- -1 / (2 * sdXZ^2) - theta_loc[2]^2 / (2 * psi^2)
        b        <-  meanXZ / sdXZ^2 - theta_loc[2] * e_star / psi^2
        # c term is not used explicitly in the scores
        mu_star  <- -b / (2 * a)
        sd2_star <- (sdXZ^2 * psi^2) / (psi^2 + theta_loc[2]^2 * sdXZ^2)

        # beta0
        beta0_star <- (theta_loc[2] * mu_star + e_star) / psi^2

        # beta1 (AW coeff)
        beta1_star <- (
          -theta_loc[2] * (sd2_star + mu_star^2) +
            mu_star * (W[2] * theta_loc[2] - e_star) +
            W[2] * e_star
        ) / psi^2

        # beta2 (Z coeff)
        beta2_star <- W[3] * beta0_star

        c(beta0_star, beta1_star, beta2_star)
      }

      if (part1) {
        # part1: (1 - D/p)*(D/p)* CCscore * psi_hat_i^T
        term <- (1 - D_i / myp_i) * (D_i / myp_i)
        outer(CCscore(), psi_hat_i()) * term
      } else {
        # part2: (1 - D/p)^2 * psi_hat_i * psi_hat_i^T
        term <- (1 - D_i / myp_i)^2
        v    <- psi_hat_i()
        outer(v, v) * term
      }
    }

    idx <- seq_len(nrow(data_yXZ))

    part1_list <- lapply(idx, function(i) {
      calculate_A_b(data_yXZ[i, , drop = FALSE], part1 = TRUE)
    })
    part1 <- Reduce("+", part1_list)

    part2_list <- lapply(idx, function(i) {
      calculate_A_b(data_yXZ[i, , drop = FALSE], part1 = FALSE)
    })
    part2 <- Reduce("+", part2_list)

    -part1 %*% solve(part2)
  }

  Lambda <- calculate_Lambda()

  ## ------------------------------------------------------
  ## 5. Define geex estimating function
  ## ------------------------------------------------------
  gee_estfun <- function(data, formula) {

    W <- stats::model.matrix(~ A + Z, data = data)
    X <- stats::model.matrix(formula, data = data)
    Y <- stats::model.response(stats::model.frame(formula = formula, data = data))
    D <- data$D
    myp <- data$myp

    function(theta) {

      # theta is (beta0, beta1, beta2)
      CCscore <- function() {
        mu <- as.numeric(X %*% theta)
        e  <- Y - mu
        as.numeric((e / psi^2) * X)
      }

      psi_hat_i <- function() {

        e_star   <- Y - as.numeric(W %*% theta)
        a        <- -1 / (2 * sdXZ^2) - theta[2]^2 / (2 * psi^2)
        b        <-  meanXZ / sdXZ^2 - theta[2] * e_star / psi^2
        mu_star  <- -b / (2 * a)
        sd2_star <- (sdXZ^2 * psi^2) / (psi^2 + theta[2]^2 * sdXZ^2)

        beta0_star <- (theta[2] * mu_star + e_star) / psi^2

        beta1_star <- (
          -theta[2] * (sd2_star + mu_star^2) +
            mu_star * (W[, 2] * theta[2] - e_star) +
            W[, 2] * e_star
        ) / psi^2

        beta2_star <- W[, 3] * beta0_star

        c(beta0_star, beta1_star, beta2_star)
      }

      term1 <- (D / myp) * CCscore()
      term2 <- (1 - D / myp) * as.numeric(Lambda %*% psi_hat_i())
      term1 + term2
    }
  }

  ## ------------------------------------------------------
  ## 6. Solve estimating equations via geex
  ## ------------------------------------------------------
  results <- geex::m_estimate(
    estFUN        = gee_estfun,
    data          = data_yXZ,
    root_control  = geex::setup_root_control(start = as.numeric(g)),
    outer_args    = list(formula = model),
    deriv_control = geex::setup_deriv_control(method = "Richardson")
  )

  beta_hat <- as.numeric(results@estimates)  # length 3

  ## ------------------------------------------------------
  ## 7. Update psi using complete cases
  ## ------------------------------------------------------
  data_cc <- data_yXZ[data_yXZ$D == 1, , drop = FALSE]
  X_cc    <- stats::model.matrix(model, data_cc)
  resid_cc <- data_cc$y - as.numeric(X_cc %*% beta_hat)
  psi_hat  <- sqrt(sum(resid_cc^2) / nrow(data_cc))

  # geex vcov gives variance only for beta (3x3)
  se_hat <- sqrt(diag(results@vcov))
  se_mat <- matrix(se_hat, nrow = 1)

  beta_mat <- rbind(beta_hat, psi_hat)
  beta_mat <- t(beta_mat)  # 1 x 4 (beta0, beta1, beta2, psi)

  list(
    beta_est = beta_mat,
    se_est   = se_mat
  )
}
