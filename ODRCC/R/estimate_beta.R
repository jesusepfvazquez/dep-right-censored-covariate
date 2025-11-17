#' Estimate regression coefficients via M-estimation
#'
#' @param data_yXZ A data frame containing the variables in \code{model}.
#'                 Must include the outcome \code{y} and covariates.
#' @param model An object of class \code{\link[stats]{formula}} specifying
#'   the regression model for \code{y} in terms of covariates.
#'
#' @return A list with components
#' \describe{
#'   \item{beta_est}{A 1 x (p + 1) matrix of estimated regression coefficients and
#'                   residual standard deviation (last element).}
#'   \item{se_est}{A 1 x p matrix of estimated standard errors for the regression
#'                 coefficients.}
#' }
#'
#' @importFrom stats lm model.matrix model.response model.frame
#' @importFrom geex m_estimate setup_root_control setup_deriv_control
#' @importFrom magrittr %>%
#' @export
estimate_beta <- function(data_yXZ, model) {

  ## m-estimating function for geex
  gee_estfun <- function(data, formula) {

    X <- stats::model.matrix(object = formula, data = data)
    Y <- stats::model.response(stats::model.frame(formula = formula, data = data))

    function(theta) {

      # useful quantities
      mu <- X %*% theta
      e  <- Y - mu

      # score equations for betas
      dotbeta <- (e / psi^2) %*% X

      c(dotbeta)
    }
  }

  ## initial parameter estimates from weighted linear model
  mylmer <- stats::lm(model, data = data_yXZ, weights = D)
  g      <- stats::coef(mylmer)
  psi    <- summary(mylmer)$sigma

  ## run M-estimation
  results <- geex::m_estimate(
    estFUN       = gee_estfun,
    data         = data_yXZ,
    root_control = geex::setup_root_control(start = c(g)),
    outer_args   = list(formula = model),
    deriv_control = geex::setup_deriv_control(method = "Richardson")
  )

  ## extract estimates and update psi
  beta_estimates <- results@estimates %>% t() %>% matrix()
  X_full <- stats::model.matrix(object = model, data = data_yXZ)

  psi_updated <- sum((data_yXZ$y - X_full %*% beta_estimates)^2 / nrow(data_yXZ)) %>%
    sqrt()

  beta_estimates <- rbind(beta_estimates, psi_updated)

  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()

  list(
    beta_est = t(beta_estimates),
    se_est   = t(se_estimates)
  )
}
