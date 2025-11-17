#' Sandwich variance for the AIPW-lambda estimator
#'
#' Computes a sandwich variance estimator for the AIPW-lambda estimator under
#' outcome-dependent right-censoring of a covariate, using the closed-form
#' augmentation and accounting for estimation of the censoring model
#' \eqn{C | (Y, Z)} via a Gumbel/Weibull AFT model.
#'
#' This implementation assumes the outcome model is \code{y ~ AW + Z}, i.e.,
#' three regression coefficients \eqn{(\beta_0, \beta_1, \beta_2)} plus a
#' residual standard deviation \eqn{\psi}.
#'
#' @param data_yXZ Data frame containing at least:
#'   \itemize{
#'     \item \code{y}: outcome,
#'     \item \code{A}: auxiliary covariate,
#'     \item \code{AW}: \code{A - X} (used in the outcome model),
#'     \item \code{W}: observed \code{W = min(X, C)},
#'     \item \code{D}: indicator \code{I(X \le C)},
#'     \item \code{Z}: covariate in the outcome and censoring models,
#'     \item \code{myp_ywz}: weights \eqn{\pi(Y, W, Z)} from AFT model (for IPW).
#'   }
#' @param mytheta Numeric vector \code{c(beta0, beta1, beta2, psi)} corresponding
#'   to the AIPW-lambda point estimates from \code{estimate_beta_aipw_lambda_close()}.
#'
#' @return A list with components:
#' \describe{
#'   \item{beta_est}{The input parameter vector \code{mytheta}.}
#'   \item{se_est}{Sandwich standard errors for \eqn{\beta_0, \beta_1, \beta_2}.}
#' }
#'
#' @importFrom stats model.matrix model.response model.frame var pweibull dweibull
#' @importFrom survival Surv survreg
#' @importFrom numDeriv jacobian hessian
#' @importFrom parallel mclapply
#' @export
var_beta_aipw_lambda <- function(data_yXZ,
                                 mytheta) {

  ## --- basic checks & unpack theta ---
  if (length(mytheta) != 4L) {
    stop("mytheta must be length 4: c(beta0, beta1, beta2, psi).")
  }
  mybeta <- mytheta[1:3]
  mypsi  <- mytheta[4]

  # outcome model is hard-coded as y ~ AW + Z
  model_outcome <- y ~ AW + Z

  ## --- Step 0: fit censoring model C | (Y, Z) via Weibull AFT and obtain Pr(D=1) ---
  wr <- survival::survreg(
    Surv(W, 1 - D) ~ y + Z,
    data = data_yXZ,
    dist = "weibull"
  )

  # IPW weights: S_C(W | covariates) under Weibull(AFT)
  shape.c <- 1 / wr$scale
  lp_c  <- exp(wr$linear.predictors)  # scale parameter
  data_yXZ$myp <- exp(-(data_yXZ$W / lp_c)^shape.c)

  myalpha <- c(stats::coef(wr), wr$scale)  # length 4
  myxi    <- c(mybeta, myalpha)

  ## --- Step 1: Lambda-related preliminaries ---
  sdXZ   <- sqrt(stats::var(data_yXZ$W[data_yXZ$D == 1]))
  meanXZ <- mean(data_yXZ$W[data_yXZ$D == 1])

  ## --- Step 2: helper functions (Gumbel representation of Weibull AFT) ---

  dgumbel <- function(x, shape, scale) {
    (1 / shape) * exp((x - scale) / shape) * exp(-exp((x - scale) / shape))
  }

  pgumbel <- function(x, shape, scale) {
    exp(-exp((x - scale) / shape))
  }

  pi_xz_func.b <- function(data. = data_yXZ, myalpha. = myalpha) {
    data.$W <- log(data.$W)
    gamma.x.vec <- myalpha.[1:3]
    shape.x     <- myalpha.[4]
    linPred_yz  <- cbind(1, data.$y, data.$Z) %*% gamma.x.vec
    pgumbel(data.$W, shape = shape.x, scale = linPred_yz)
  }

  jacobian_pi.b <- function(data) {
    deriv <- numDeriv::jacobian(
      function(x) pi_xz_func.b(data. = data, myalpha. = x),
      x      = myalpha,
      method = "Richardson"
    )
    matrix(deriv, ncol = 1)
  }

  alpha.logLik.b <- function(data. = data_yXZ, myalpha. = myalpha) {
    data.$W <- log(data.$W)
    gamma.x.vec <- myalpha.[1:3]
    shape.x     <- myalpha.[4]
    linPred_yz  <- cbind(1, data.$y, data.$Z) %*% gamma.x.vec

    val <- (1 - data.$D) * dgumbel(data.$W, shape = shape.x, scale = linPred_yz) +
      data.$D * pgumbel(data.$W, shape = shape.x, scale = linPred_yz)
    log(val)
  }

  alpha.jacobian.b <- function(data) {
    deriv <- numDeriv::jacobian(
      function(x) alpha.logLik.b(data. = data, myalpha. = x),
      x      = myalpha,
      method = "Richardson"
    )
    matrix(deriv, ncol = 1)
  }

  alpha.hessian.b <- function(data) {
    numDeriv::hessian(
      function(x) alpha.logLik.b(data. = data, myalpha. = x),
      x      = myalpha,
      method = "Richardson"
    )
  }

  ## --- 1. A_alpha (Hessian of censoring params) ---
  calculate_A_alpha <- function() {
    parts <- lapply(seq_len(nrow(data_yXZ)), function(i) {
      alpha.hessian.b(data_yXZ[i, , drop = FALSE])
    })
    Reduce("+", parts) / nrow(data_yXZ)
  }
  my_A_alpha_inv <- -solve(calculate_A_alpha())

  ## --- 2. A_ipw for theta (beta) ---
  calculate_A_ipw <- function() {

    theta <- mybeta
    psi   <- mypsi

    A_ipw.b <- function(data) {
      X <- stats::model.matrix(~ AW + Z, data = data)
      Y <- stats::model.response(stats::model.frame(y ~ AW + Z, data = data))
      D <- data$D

      IPWscore <- function(myalpha.) {
        mu <- X %*% theta
        e  <- Y - mu
        dotbeta <- (e / psi^2) %*% X
        myp <- pi_xz_func.b(data. = data, myalpha. = myalpha.)
        matrix(dotbeta*c(D/myp), ncol = 1)
      }

      numDeriv::jacobian(
        function(x) IPWscore(myalpha. = x),
        x      = myalpha,
        method = "Richardson"
      )
    }

    parts <- lapply(seq_len(nrow(data_yXZ)), function(i) {
      A_ipw.b(data_yXZ[i, , drop = FALSE])
    })
    Reduce("+", parts) / nrow(data_yXZ)
  }
  my_A_ipw <- calculate_A_ipw()

  ## --- 3. A_alpha_theta (cross term between alpha and augmentation) ---
  calculate_A_alpha_theta <- function() {

    theta <- mybeta
    psi   <- mypsi

    phi_yz.b <- function(data.) {
      W <- stats::model.matrix(y ~ A + Z, data = data.)
      X <- stats::model.matrix(y ~ AW + Z, data = data.)
      Y <- stats::model.response(stats::model.frame(y ~ AW + Z, data = data.))
      D <- data.$D

      psi_hat_i <- function() {

        e_star <- Y - (W %*% theta)
        a      <- -1 / (2 * sdXZ^2) - theta[2]^2 / (2 * psi^2)
        b      <-  meanXZ / sdXZ^2 - theta[2] * e_star / psi^2
        mu_star  <- -b / (2 * a)
        sd2_star <- (sdXZ^2 * psi^2) / (psi^2 + theta[2]^2 * sdXZ^2)

        beta0_star <- (theta[2] * mu_star + e_star) / psi^2
        beta1_star <- (
          -theta[2] * (sd2_star + mu_star^2) +
            mu_star * (W[2] * theta[2] - e_star) +
            W[2] * e_star
        ) / psi^2
        beta2_star <- W[3] * beta0_star

        c(beta0_star, beta1_star, beta2_star)
      }

      psi_hat_i()
    }

    A_alpha_theta.b <- function(mydata) {
      D   <- mydata$D
      myp <- mydata$myp
      matrix(D * phi_yz.b(mydata) / (myp^2), ncol = 1) %*%
        t(alpha.jacobian.b(mydata))
    }

    parts <- lapply(seq_len(nrow(data_yXZ)), function(i) {
      A_alpha_theta.b(data_yXZ[i, , drop = FALSE])
    })
    Reduce("+", parts) / nrow(data_yXZ)
  }
  my_A_alpha_theta <- calculate_A_alpha_theta()

  ## --- 4. Lambda matrix (myA in your notation) ---
  calculate_Lambda <- function() {

    theta <- mybeta
    psi   <- mypsi

    calculate_A.b <- function(data, part1 = TRUE) {

      W <- stats::model.matrix(~ A + Z, data = data)
      X <- stats::model.matrix(y ~ AW + Z, data = data)
      Y <- stats::model.response(stats::model.frame(y ~ AW + Z, data = data))
      D <- data$D
      myp <- data$myp

      CCscore <- function() {
        mu <- X %*% theta
        e  <- Y - mu
        as.numeric(c(e / psi^2) * X)
      }

      psi_hat_i <- function() {

        e_star <- Y - (W %*% theta)
        a      <- -1 / (2 * sdXZ^2) - theta[2]^2 / (2 * psi^2)
        b      <-  meanXZ / sdXZ^2 - theta[2] * e_star / psi^2
        mu_star  <- -b / (2 * a)
        sd2_star <- (sdXZ^2 * psi^2) / (psi^2 + theta[2]^2 * sdXZ^2)

        beta0_star <- (theta[2] * mu_star + e_star) / psi^2
        beta1_star <- (
          -theta[2] * (sd2_star + mu_star^2) +
            mu_star * (W[2] * theta[2] - e_star) +
            W[2] * e_star
        ) / psi^2
        beta2_star <- W[3] * beta0_star

        c(beta0_star, beta1_star, beta2_star)
      }

      if (part1) {
        # part1: E[ (D/p) * CCscore * { (1 - D/p) psi_hat_i + A_alpha_theta A_alpha_inv alpha' }^T ]
        v2 <- (1 - D / myp) * psi_hat_i() +
          my_A_alpha_theta %*% my_A_alpha_inv %*% alpha.jacobian.b(data)
        v1 <- matrix(D * CCscore() / myp, ncol = 1)
        v1 %*% t(v2)
      } else {
        # part2: E[ { (1 - D/p) psi_hat_i + A_alpha_theta A_alpha_inv alpha' }^{\otimes 2} ]
        v <- (1 - D / myp) * psi_hat_i() +
          my_A_alpha_theta %*% my_A_alpha_inv %*% alpha.jacobian.b(data)
        vv <- matrix(v, ncol = 1)
        vv %*% t(vv)
      }
    }

    part1 <- lapply(seq_len(nrow(data_yXZ)), function(i) {
      calculate_A.b(data_yXZ[i, , drop = FALSE], part1 = TRUE)
    })
    part1 <- Reduce("+", part1) / nrow(data_yXZ)

    part2 <- lapply(seq_len(nrow(data_yXZ)), function(i) {
      calculate_A.b(data_yXZ[i, , drop = FALSE], part1 = FALSE)
    })
    part2 <- Reduce("+", part2) / nrow(data_yXZ)

    -part1 %*% solve(part2)
  }

  Lambda <- calculate_Lambda()

  ## --- 5. Influence function of AIPW-lambda estimator ---
  aipw.lambda.b <- function(theta, data.b) {

    W <- stats::model.matrix(~ A + Z, data = data.b)
    X <- stats::model.matrix(y ~ AW + Z, data = data.b)
    Y <- stats::model.response(stats::model.frame(y ~ AW + Z, data = data.b))
    D <- data.b$D
    myp <- data.b$myp
    psi <- mypsi

    CCscore <- function() {
      mu <- X %*% theta
      e  <- Y - mu
      as.numeric( c(e / psi^2) * X)
    }

    psi_hat_i <- function() {
      e_star <- Y - (W %*% theta)
      a      <- -1 / (2 * sdXZ^2) - theta[2]^2 / (2 * psi^2)
      b      <-  meanXZ / sdXZ^2 - theta[2] * e_star / psi^2
      mu_star  <- -b / (2 * a)
      sd2_star <- (sdXZ^2 * psi^2) / (psi^2 + theta[2]^2 * sdXZ^2)

      beta0_star <- (theta[2] * mu_star + e_star) / psi^2
      beta1_star <- (
        -theta[2] * (sd2_star + mu_star^2) +
          mu_star * (W[2] * theta[2] - e_star) +
          W[2] * e_star
      ) / psi^2
      beta2_star <- W[3] * beta0_star

      c(beta0_star, beta1_star, beta2_star)
    }

    CC_term <- CCscore() * (D / myp)
    aug_core <- Lambda %*% psi_hat_i()
    alpha_adj <- (my_A_ipw + Lambda %*% my_A_alpha_theta) %*%
      (my_A_alpha_inv %*% alpha.jacobian.b(data.b))

    CC_term + (1 - D / myp) * aug_core + alpha_adj
  }

  ## --- 6. Sandwich A and B for beta ---
  calculate.B <- function() {
    parts <- lapply(seq_len(nrow(data_yXZ)), function(i) {
      s_i <- aipw.lambda.b(theta = mybeta, data.b = data_yXZ[i, , drop = FALSE])
      s_i %*% t(s_i)
    })
    Reduce("+", parts)
  }
  myB <- calculate.B()

  calculate.A <- function() {
    parts <- parallel::mclapply(
      X = seq_len(nrow(data_yXZ)),
      FUN = function(i) {
        numDeriv::jacobian(
          function(x) aipw.lambda.b(theta = x, data.b = data_yXZ[i, , drop = FALSE]),
          x      = mybeta,
          method = "Richardson"
        )
      }
    )
    Reduce("+", parts)
  }
  myA <- calculate.A()
  myA.inv <- solve(myA)

  sand.var <- myA.inv %*% myB %*% t(myA.inv)
  se_beta  <- sqrt(diag(sand.var))

  list(
    beta_est = mytheta,
    se_est   = se_beta
  )
}
