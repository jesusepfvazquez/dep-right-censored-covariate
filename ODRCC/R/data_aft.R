#' Simulate AFT data with an outcome-dependent right-censored covariate
#'
#' Generates a simulated dataset under an accelerated failure time (AFT) model
#' with a right-censored covariate subject to outcome-dependent censoring.
#' The function returns the true covariate \code{X}, censoring time \code{C},
#' observed covariate \code{W = min(X, C)}, event indicator \code{D = I(X <= C)},
#' and several derived quantities such as weights and conditional expectations.
#'
#' @param nSubjects Integer. Number of subjects to simulate. Defaults to \code{10^3}.
#'
#' @return A data frame with \code{nSubjects} rows containing:
#' \describe{
#'   \item{y}{Continuous outcome.}
#'   \item{Z}{Standard normal covariate.}
#'   \item{D}{Event indicator \code{I(X <= C)}.}
#'   \item{A}{Observed auxiliary covariate.}
#'   \item{X}{True covariate subject to censoring.}
#'   \item{AX}{\code{A - X}.}
#'   \item{C}{Censoring time.}
#'   \item{W}{\code{min(X, C)}.}
#'   \item{AW}{\code{A - W}.}
#'   \item{meanCYZ}{Mean of \code{C} given \code{(Y, Z)} on the original scale.}
#'   \item{meanXZ}{Mean of \code{X} given \code{(Z)} on the original scale.}
#'   \item{b}{Subject ID.}
#'   \item{myp_uniform}{Random weights (Uniform(0.1, 0.9)).}
#'   \item{myp_ywz_oracle}{Oracle survival probability \code{P(C >= W | Y, Z)}.}
#'   \item{myp_ywz}{Estimated survival probability using a Weibull AFT model.}
#'   \item{meanXD0YZ}{\eqn{E[X | D = 0, Y, Z]} computed by numerical integration.}
#'   \item{AX_yz}{\code{A - X} if \code{D = 1}, otherwise \code{A - E[X | D = 0, Y, Z]}.}
#'   \item{AX_z}{\code{A - X} if \code{D = 1}, otherwise \code{A - E[X | Z]}.}
#' }
#'
#' @importFrom survival Surv survreg
#' @importFrom stats rnorm rweibull runif dweibull dnorm integrate predict
#' @export
data_aft <- function(nSubjects = 10^3) {

  ## fixed parameters
  epsilon_sig <- 1
  beta.vec    <- c(1, 1, 1)

  ## subject IDs and baseline covariates
  b   <- seq_len(nSubjects)
  int <- rep(1, nSubjects)

  Z <- stats::rnorm(nSubjects, mean = 0, sd = 1)
  A <- stats::rnorm(nSubjects, mean = 2, sd = 1)

  ## Generate X | Z
  p.x          <- 0.8
  gamma.x.vec  <- c(1, 0.1)
  linPred_z    <- cbind(int, Z) %*% gamma.x.vec
  errors.x.g   <- log(stats::rweibull(nSubjects, shape = 1, scale = 1))
  X            <- exp(linPred_z + p.x * errors.x.g)
  AX           <- A - X

  ## Generate Y | (A, X, Z)
  er       <- stats::rnorm(nSubjects, mean = 0, sd = epsilon_sig)
  x_mat    <- cbind(int, AX, Z)
  linPred  <- x_mat %*% beta.vec
  y        <- linPred + er

  ## Generate C | (Y, Z)
  p              <- 2
  gamma.vec      <- c(1, -0.5, 0.5)
  linPred_yz     <- cbind(int, y, Z) %*% gamma.vec
  errors.g       <- log(stats::rweibull(nSubjects, shape = 1, scale = 1))
  C              <- exp(linPred_yz + p * errors.g)
  D              <- as.numeric(X <= C)
  W              <- pmin(X, C)
  AW             <- A - W

  ## build base data frame
  dat <- data.frame(
    y        = as.numeric(y),
    Z        = as.numeric(Z),
    D        = as.numeric(D),
    A        = as.numeric(A),
    X        = as.numeric(X),
    AX       = as.numeric(AX),
    C        = as.numeric(C),
    W        = as.numeric(W),
    AW       = as.numeric(AW),
    meanCYZ  = as.numeric(exp(linPred_yz)),
    meanXZ   = as.numeric(exp(linPred_z)),
    b        = as.integer(b)
  )

  #############################################################
  ##################### Calculate weights #####################
  #############################################################

  ## random weights
  dat$myp_uniform <- stats::runif(nSubjects, min = 0.1, max = 0.9)

  ## oracle weights: S_C(W | Y, Z)
  dat$myp_ywz_oracle <- exp(-(dat$W / exp(linPred_yz))^(1 / p))

  ## estimated weights using Weibull AFT model C | (Y, Z)
  wr <- survival::survreg(survival::Surv(W, 1 - D) ~ y + Z, data = dat, dist = "weibull")
  lp <- stats::predict(wr, newdata = dat)
  shape_gamma <- 1 / wr$scale  # shape parameter of Weibull

  dat$myp_ywz <- exp(-(dat$W / lp)^shape_gamma)

  #############################################################
  ########## E[X | D = 0, Y, Z] via numerical integration #####
  #############################################################

  get.meanXD0YZ <- function(i) {

    num <- stats::integrate(
      f = function(x) {
        x *
          stats::dweibull(x, shape = 1 / p.x, scale = dat$meanXZ[i]) *
          stats::dnorm(
            dat$y[i] - (beta.vec[1] +
                          (dat$A[i] - x) * beta.vec[2] +
                          dat$Z[i] * beta.vec[3]),
            mean = 0,
            sd   = epsilon_sig
          )
      },
      lower = dat$W[i],
      upper = Inf
    )$value

    den <- stats::integrate(
      f = function(x) {
        stats::dweibull(x, shape = 1 / p.x, scale = dat$meanXZ[i]) *
          stats::dnorm(
            dat$y[i] - (beta.vec[1] +
                          (dat$A[i] - x) * beta.vec[2] +
                          dat$Z[i] * beta.vec[3]),
            mean = 0,
            sd   = epsilon_sig
          )
      },
      lower = dat$W[i],
      upper = Inf
    )$value

    num / den
  }

  dat$meanXD0YZ <- vapply(seq_len(nSubjects), get.meanXD0YZ, numeric(1))

  ## imputed AX
  dat$AX_yz <- ifelse(dat$D == 1, dat$AX, dat$A - dat$meanXD0YZ)

  ## E[X | Z] = meanXZ * Gamma(1 + 1/p.x)
  dat$AX_z  <- ifelse(
    dat$D == 1,
    dat$AX,
    dat$A - dat$meanXZ * base::gamma(1 + 1 / p.x)
  )

  return(dat)
}
