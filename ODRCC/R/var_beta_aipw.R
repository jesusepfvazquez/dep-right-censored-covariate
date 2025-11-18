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
#'     \item \code{D}: indicator \code{I(X <= C)},
#'     \item all covariates appearing in \code{model},
#'     \item all covariates appearing in \code{model_weights} and \code{model_xz}.
#'   }
#' @param theta Numeric vector \code{c(beta, psi)} from the AIPW estimator,
#'   where \code{beta} has length equal to the number of columns in
#'   \code{model.matrix(model, data_yXZ)} and \code{psi} is the residual
#'   standard deviation.
#' @param lbound,ubound Numeric lower and upper bounds for the numerical
#'   integration over \code{X} in the augmentation term (defaults: 0, 50).
#'
#' @return A list with components
#' \describe{
#'   \item{beta_est}{Estimated regression coefficients \eqn{\beta}.}
#'   \item{psi_est}{Estimated residual standard deviation \eqn{\psi}.}
#'   \item{se_beta}{Sandwich standard errors for \eqn{\beta}.}
#' }
#'
#' @importFrom stats model.matrix model.response model.frame dnorm dweibull pweibull integrate
#' @importFrom stats terms reformulate as.formula
#' @importFrom survival Surv survreg
#' @importFrom numDeriv jacobian
#' @export
#'
var_beta_aipw <- function(data_yXZ, theta, lbound = 0, ubound = 50){

  ##############################################################
  # Step 1: define parameters values for integration component #
  ##############################################################
  wr <- survreg(Surv(W, 1-D) ~ y + Z, data = data_yXZ, dist="w")
  myalpha = c(coef(wr), wr$scale)
  wr_z <- survreg(Surv(W, D) ~ Z, data = data_yXZ, dist="w")
  mygamma = c(coef(wr_z), wr_z$scale)

  # theta = est32$beta_est
  mybeta = theta[1:3]
  mypsi = theta[4]
  data_yXZ$meanXZ = exp(mygamma[1] + mygamma[2]*data_yXZ$Z)
  data_yXZ$meanCYZ = exp(wr$linear.predictors)

  shape.c = 1/wr$scale
  shape.x = 1/mygamma[3]
  data_yXZ$myp = data_yXZ$myp_ywz

  ##############################################################
  # Step 2: Define helper functions
  ##############################################################
  dgumbel = function(x, shape, scale){
    (1/shape)*exp((x-scale)/shape)*exp(-exp((x-scale)/shape))
  }

  pgumbel = function(x, shape, scale){
    exp(-exp((x-scale)/shape))
  }

  pi_xz_func.b = function(data. = data_yXZ, myalpha. = myalpha){
    data.$W = log(data.$W)
    gamma.x.vec = myalpha.[1:3]
    shape.x = myalpha.[4]
    linPred_yz = cbind(rep(1,nrow(data.)),data.$y, data.$Z)%*%gamma.x.vec
    return(pgumbel(data.$W, shape=shape.x, scale = linPred_yz))
  }

  jacobian_pi.b = function(data){
    myderiv = numDeriv::jacobian(function(x)
      pi_xz_func.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(matrix(ncol=1, myderiv))
  }

  alpha.logLik.b = function(data.=data_yXZ, myalpha. = myalpha){
    data.$W = log(data.$W)
    gamma.x.vec = myalpha.[1:3]
    shape.x = myalpha.[4]
    linPred_yz = cbind(rep(1,nrow(data.)),data.$y, data.$Z)%*%gamma.x.vec

    myreturn = (1-data.$D)*dgumbel(data.$W, shape=shape.x, scale = linPred_yz)+
      data.$D*pgumbel(data.$W, shape=shape.x, scale = linPred_yz)

    return(log(myreturn))
  }

  alpha.jacobian.b = function(data){
    myderiv = numDeriv::jacobian(function(x)
      alpha.logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(matrix(ncol=1, myderiv))
  }

  alpha.hessian.b = function(data){
    myderiv = numDeriv::hessian(function(x)
      alpha.logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(myderiv)
  }

  ## 1. The A_alpha matrix
  calculate_A_alpha = function(){
    # calculate the hessian matrix and return
    part1= lapply(1:nrow(data_yXZ), function(x) alpha.hessian.b(data_yXZ[x,]))
    part1 = Reduce("+", part1)/nrow(data_yXZ)
    return(part1)
  }
  my_A_alpha_inv = -solve(calculate_A_alpha())

  ## 2. The A_aipw of theta
  calculate_A_aipw = function(){

    theta = mybeta
    psi = mypsi

    A_aipw.b = function(data) {

      W <- model.matrix(object = ~ A + Z, data = data)
      X <- model.matrix(object =  y~ AW+Z, data = data)
      Y <- model.response(model.frame(formula =  y~ AW+Z, data = data))
      D <- data$D
      z <- data$Z

      # need these for integration components
      meanXZ = data$meanXZ
      meanCYZ = data$meanCYZ

      # calculate the IPW score equation
      AIPWscore <- function(myalpha.){

        myp = pi_xz_func.b(data.=data, myalpha. = myalpha.)

        CCscore <- function(){
          # useful quantities
          mu <- X %*% theta
          e = Y - mu
          # score equations for betas
          dotbeta = (e/psi^2)%*%X
          return(c(dotbeta))
        }

        #######################
        # Augmented Component #
        #######################
        # augmentation component
        psi_hat_i <- function(){

          # top integral
          likelihood_int = function(t=1){
            return(dnorm(Y -  cbind(W[,1],W[,2] - t, W[,3]) %*% theta, 0, sd = psi))
          }

          # score of integral
          score_int <- function(t=1,j=1){
            # generate empty array of values to save
            # create temporary X matrix
            Wtemp <- cbind(W[,1], W[,2]- t, W[,3])
            mu <- Wtemp %*% theta
            e = Y - mu
            # score equations for theta
            dotbeta = (e/psi^2)%*%Wtemp
            return(c(dotbeta)[j])
          }

          #top integral#
          ## switch tt instead of t
          integral_func_num <- function(t=1,j=1){
            val = rep(NA, length(t))
            for (i in 1:length(t)){
              val[i] = score_int(t[i],j=j)*likelihood_int(t[i])*
                dweibull(t[i], shape=shape.x, scale = meanXZ) *
                (1-1/pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE))
            }
            return(val)
          }

          ### Evaluate
          j_numerator_integrate <- function(jj=1) {
            integrate(function(y) {integral_func_num(t=y,j=jj)},
                      lower = lbound, upper = ubound)$value}
          v.area_num <- Vectorize(j_numerator_integrate)
          numerator = v.area_num(1:3)

          # bottom integral
          integral_func_denom <- function(t=1){
            val = rep(NA, length(t))
            for (i in 1:length(t)){
              val[i] = likelihood_int(t[i])*
                dweibull(t[i], shape=shape.x, scale = meanXZ)*
                (1-1/pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE)) # (1-1/pi)
            }
            return(val)
          }

          ### Evaluate
          denominator = integrate(integral_func_denom,
                                  lower = lbound, upper = ubound)$value

          # Compute Augmented part
          return(numerator/denominator)
        }


        myreturn = CCscore()*c(D/myp) + c(1-D/myp)*psi_hat_i()
        return(myreturn)
      }

      # calculate the partial derivative
      myderiv = numDeriv::jacobian(function(x)
        AIPWscore(myalpha. = x), myalpha, method = "Richardson")
      return(myderiv)
    }

    # calculate the hessian matrix and return
    part1= lapply(1:nrow(data_yXZ), function(x) A_aipw.b(data_yXZ[x,]))
    part1 = Reduce("+", part1)/nrow(data_yXZ)
    return(part1)
  }
  my_A_aipw = calculate_A_aipw()

  ##############################################################
  # Step 3: Calculate the A and the B matrices
  ##############################################################

  aipw.lambda.b <- function(theta, data.b){

    psi = mypsi

    # W <- model.matrix(object = ~ A + Z + Z2 + Z3, data = data.b)
    # X <- model.matrix(object = y ~ AW + Z + Z2 + Z3, data = data.b)
    # Y <- model.response(model.frame(formula = y ~ AW + Z + Z2 + Z3, data = data.b))
    W <- model.matrix(object = ~ A + Z , data = data.b)
    X <- model.matrix(object = y ~ AW + Z, data = data.b)
    Y <- model.response(model.frame(formula = y ~ AW + Z, data = data.b))
    D <- data.b$D
    myp <- data.b$myp

    meanXZ = data.b$meanXZ
    meanCYZ = data.b$meanCYZ

    #################
    # 1st Component #
    #################
    CCscore <- function(){
      # useful quantities
      mu <- X %*% theta
      e = Y - mu
      # score equations for betas
      dotbeta = (e/psi^2)%*%X
      return(c(dotbeta))
    }

    #######################
    # Augmented Component #
    #######################

    psi_hat_i <- function(){

      # top integral
      likelihood_int = function(t=1){
        return(dnorm(Y -  cbind(W[,1],W[,2] - t, W[,3]) %*% theta, 0, sd = psi))
      }

      # score of integral
      score_int <- function(t=1,j=1){
        # generate empty array of values to save
        # create temporary X matrix
        Wtemp <- cbind(W[,1], W[,2]- t, W[,3])
        mu <- Wtemp %*% theta
        e = Y - mu
        # score equations for theta
        dotbeta = (e/psi^2)%*%Wtemp
        return(c(dotbeta)[j])
      }

      #top integral#
      ## switch tt instead of t
      integral_func_num <- function(t=1,j=1){
        val = rep(NA, length(t))
        for (i in 1:length(t)){
          val[i] = score_int(t[i],j=j)*likelihood_int(t[i])*
            dweibull(t[i], shape=shape.x, scale = meanXZ) *
            (1-1/pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE)) # (1-1/pi)
          # (1-1/(smallnum+pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE))) # (1-1/pi)
        }
        return(val)
      }

      ### Evaluate
      j_numerator_integrate <- function(jj=1) {
        integrate(function(y) {integral_func_num(t=y,j=jj)},
                  lower = lbound, upper = ubound)$value}
      v.area_num <- Vectorize(j_numerator_integrate)
      numerator = v.area_num(1:3)

      # bottom integral
      integral_func_denom <- function(t=1){
        val = rep(NA, length(t))
        for (i in 1:length(t)){
          val[i] = likelihood_int(t[i])*
            # pweibull(C, shape=gamma, scale = lp, lower.tail = FALSE)
            dweibull(t[i], shape=shape.x, scale = meanXZ)*
            (1-1/pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE)) # (1-1/pi)
          # (1-1/(smallnum+pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE))) # (1-1/pi)
        }
        return(val)
      }

      ### Evaluate
      denominator = integrate(integral_func_denom,
                              lower = lbound, upper = ubound)$value

      # Compute Augmented part
      return(numerator/denominator)
    }

    # print(theta)
    AIPW_est = CCscore()*D/myp + (1 - D/myp)*psi_hat_i() +
      my_A_aipw%*%(my_A_alpha_inv%*%alpha.jacobian.b(data.b))
    # print(ACC_est)
    return(AIPW_est)

  }

  #####################################################
  # Step 6: Calculate the A and the B matrices

  calculate.B = function(){
    part1 = lapply(1:nrow(data_yXZ), function(i){
      dotalpha = aipw.lambda.b(theta = mybeta, data.b = data_yXZ[i,])
      return(dotalpha %*% t(dotalpha))}
    )
    part1 = Reduce("+", part1)
    return(part1)
  }
  myB = calculate.B()

  # A matrix
  calculate.A = function(){
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(i)
      dotalpha = numDeriv::jacobian(function(x)
        aipw.lambda.b(theta = x, data.b = data_yXZ[i,]), mybeta, method = "Richardson")
    )
    part1 = Reduce("+", part1) #/nrow(data_yXZ)
  }
  myA = calculate.A()
  myA.inv = solve(myA)

  # calculate sandwich variance estimate
  sand.var = myA.inv%*%myB%*%t(myA.inv)
  sand.var = sqrt(diag(sand.var))
  # return values
  return(list(beta_est = theta, se_est = sand.var))
}
