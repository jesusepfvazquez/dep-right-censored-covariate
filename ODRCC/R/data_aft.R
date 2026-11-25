#' Simulate data with an outcome-dependent right-censored covariate
#'
#' Generates one dataset from the simulation design of the manuscript:
#' \deqn{\log X = 0.1 + 0.1 Z + 0.5 \log E_X,}
#' \deqn{Y = 1 + (A - X) + Z + \epsilon, \quad \epsilon \sim N(0, 1),}
#' \deqn{\log C = c_{int} + 0.1 Y + 0.1 Z + p_c \log E_C,}
#' with \eqn{Z \sim N(0, 1)}, \eqn{A \sim N(2, 1)}, and \eqn{E_X, E_C \sim}
#' Exp(1), so that \eqn{X | Z} and \eqn{C | (Y, Z)} follow Weibull AFT models and
#' \eqn{C} is independent of \eqn{X} given \eqn{(Y, Z)}. The intercept
#' \code{c_int} (\eqn{\tilde\eta_0} in the manuscript) sets the right-censoring
#' rate, and \code{p.c} (\eqn{\tilde\sigma_\eta}) is the scale of \eqn{\log C}.
#' The scale of \eqn{\log X} is \eqn{\tilde\sigma_\gamma = 0.5}, so the scale
#' ratio \eqn{\tilde\sigma_\eta / \tilde\sigma_\gamma} in the manuscript is
#' \code{p.c / 0.5}.
#'
#' With \code{p.c = 1.5}, \code{c_int = 4.5, 4, 1.5, 1, 0.5, -0.5, -1, -1.5} give
#' about 4\%, 6\%, 26\%, 34\%, 44\%, 66\%, 76\%, and 85\% right-censoring.
#' \code{set.seed(s); data_aft(1000, c_int, 1.5)} regenerates dataset \code{s} of
#' the manuscript's simulations that vary the right-censoring rate (the code in
#' \code{simulations/} draws the same random numbers in the same order).
#'
#' @param nSubjects Number of subjects. Defaults to 1,000.
#' @param c_int Intercept of the censoring model. Defaults to -0.5 (about 66\%
#'   right-censoring).
#' @param p.c Scale of \eqn{\log C}, \eqn{\tilde\sigma_\eta}. Defaults to 1.5.
#'
#' @return A data frame with one row per subject:
#' \describe{
#'   \item{y}{Outcome.}
#'   \item{Z}{Fully observed covariate.}
#'   \item{D}{Indicator that \eqn{X} is observed, \eqn{I(X \le C)}.}
#'   \item{A}{Age at baseline.}
#'   \item{X}{Right-censored covariate (age at the event; not observed when \code{D = 0}).}
#'   \item{AX}{\code{A - X}, the time to the event.}
#'   \item{C}{Censoring time.}
#'   \item{W}{Observed \code{min(X, C)}.}
#'   \item{AW}{\code{A - W}.}
#'   \item{b}{Subject ID.}
#' }
#'
#' @examples
#' set.seed(1)
#' dat <- data_aft(500)
#' mean(1 - dat$D)   # right-censoring rate
#'
#' @importFrom stats rnorm rweibull
#' @export
data_aft <- function(nSubjects = 1000, c_int = -0.5, p.c = 1.5) {
  Z <- stats::rnorm(nSubjects, 0, 1)
  A <- stats::rnorm(nSubjects, 2, 1)
  X <- exp(0.1 + 0.1 * Z + 0.5 * log(stats::rweibull(nSubjects, 1, 1)))
  AX <- A - X
  y <- 1 + AX + Z + stats::rnorm(nSubjects, 0, 1)
  C <- exp(c_int + 0.1 * y + 0.1 * Z + p.c * log(stats::rweibull(nSubjects, 1, 1)))
  D <- as.numeric(X <= C)
  W <- pmin(X, C)
  data.frame(y = y, Z = Z, D = D, A = A, X = X, AX = AX, C = C, W = W, AW = A - W,
             b = seq_len(nSubjects))
}
