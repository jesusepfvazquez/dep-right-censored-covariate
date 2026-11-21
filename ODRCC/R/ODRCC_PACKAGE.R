#' ODRCC: estimators for regression with an outcome-dependent right-censored covariate
#'
#' The ODRCC package implements the estimators in the manuscript for a linear
#' model \eqn{Y = \beta_0 + \beta_1 (A - X) + \beta_Z' Z + \epsilon} when the
#' covariate \eqn{X} is right-censored by \eqn{C} and the censoring depends on the
#' outcome (\eqn{C} is independent of \eqn{X} given \eqn{(Y, Z)}):
#' \itemize{
#'   \item \code{\link{estimate_beta_ipw}}: the IPW estimator;
#'   \item \code{\link{estimate_beta_mle}}: the joint MLE;
#'   \item \code{\link{estimate_beta_aipw}}: the unguarded AIPW estimator;
#'   \item \code{\link{estimate_beta_aipw_lambda}}: the guarded AIPW estimator;
#'   \item \code{\link{estimate_beta_cc}} and \code{\link{estimate_beta}}: the
#'     complete-case, oracle, and naive estimators, for comparison;
#'   \item \code{\link{ess_weights}}: the effective sample size of the IPW weights,
#'     which guides the choice between the estimators;
#'   \item \code{\link{data_aft}}: the data generator of the first simulation study.
#' }
#' Every estimator returns estimates of \eqn{\theta = (\beta, \sigma)} and sandwich
#' standard errors.
#'
#' @name ODRCC
#' @keywords internal
#' @author Jesus E. Vazquez
"_PACKAGE"
NULL
