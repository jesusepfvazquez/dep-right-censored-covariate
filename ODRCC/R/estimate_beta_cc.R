#' Complete-case estimator for regression parameters
#'
#' Fits a regression model using only complete cases, defined by \code{D = 1},
#' via M-estimation using the \code{geex} framework. This is the complete-case
#' analogue of \code{\link{estimate_beta}}, and is primarily intended for
#' comparison in simulation studies. Not consistent under outcome dependent censoring.
#'
#' @param data_yXZ A data frame containing at least the variables in \code{model}
#'   and a binary indicator \code{D}, where \code{D = 1} denotes a complete case.
#' @param model An object of class \code{\link[stats]{formula}} specifying the
#'   regression model for the outcome.
#'
#' @return A list with components
#' \describe{
#'   \item{beta_est}{A 1 x (p + 1) matrix of regression coefficients followed by
#'                   the residual standard deviation (last element).}
#'   \item{se_est}{A 1 x p matrix of estimated standard errors for the
#'                 regression coefficients.}
#' }
#'
#' @importFrom stats lm coef model.matrix model.response model.frame
#' @importFrom geex m_estimate setup_root_control setup_deriv_control
#' @importFrom magrittr %>%
#' @export
estimate_beta_cc <- function(data_yXZ, model) {

  ## M-estimating function for complete-case analysis
  gee_estfun <- function(data, formula) {

    X <- stats::model.matrix(object = formula, data = data)
    Y <- stats::model.response(stats::model.frame(formula = formula, data = data))
    D <- data$D

    function(theta) {

      # linear predictor and residuals
      mu <- X %*% theta
      e  <- Y - mu

      # score equations for betas
      dotbeta <- (e / psi^2) %*% X

      # only keep contributions from complete cases (D = 1)
      D * c(dotbeta)
    }
  }

  ## initial parameter estimates using weighted least squares with D as weights
  mylmer <- stats::lm(model, data = data_yXZ, weights = D)
  g      <- stats::coef(mylmer)
  psi    <- summary(mylmer)$sigma

  ## run M-estimation
  results <- geex::m_estimate(
    estFUN        = gee_estfun,
    data          = data_yXZ,
    root_control  = geex::setup_root_control(start = c(g)),
    outer_args    = list(formula = model),
    deriv_control = geex::setup_deriv_control(method = "Richardson")
  )

  ## extract estimates
  beta_estimates <- results@estimates %>% t() %>% matrix()

  ## recompute residual SD using only D = 1 rows
  data_cc <- subset(data_yXZ, D == 1)
  X_cc    <- stats::model.matrix(object = model, data = data_cc)

  psi_updated <- sum(
    (data_cc$y - X_cc %*% beta_estimates)^2 / nrow(data_cc)
  ) %>%
    sqrt()

  beta_estimates <- rbind(beta_estimates, psi_updated)

  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()

  list(
    beta_est = t(beta_estimates),
    se_est   = t(se_estimates)
  )
}
