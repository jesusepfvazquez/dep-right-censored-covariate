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
#'     \item \code{D}: indicator \code{I(X <= C)},
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
var_beta_aipw_lambda <- function(data_yXZ, mytheta){

  # mytheta = est32$beta_est
  wr <- survreg(Surv(W, 1-D) ~ y + Z, data = data_yXZ, dist="w")
  myalpha = c(coef(wr), wr$scale)
  mybeta = mytheta[1:(length(mytheta)-1)]
  mypsi = mytheta[length(mytheta)]
  myxi = c(mybeta,myalpha)

  ######################################################
  # Step 1: Calculate Lambda using the values for theta
  sdXZ = var(data_yXZ$W[data_yXZ$D==1])^0.5
  meanXZ = mean(data_yXZ$W[data_yXZ$D==1])
  data_yXZ$myp = data_yXZ$myp_ywz

  #####################################################
  # Step 2: Define helper functions
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
  # pi_xz_func.b(dat[1,])

  jacobian_pi.b = function(data){
    myderiv = numDeriv::jacobian(function(x)
      pi_xz_func.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(matrix(ncol=1, myderiv))
  }
  # jacobian_pi.b(dat[1,])

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
  # alpha.jacobian.b(data_yXZ[1,])

  alpha.hessian.b = function(data){
    myderiv = numDeriv::hessian(function(x)
      alpha.logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(myderiv)
  }
  # alpha.hessian.b(data_yXZ[1,])

  ## 1. The A_alpha matrix
  calculate_A_alpha = function(){
    # calculate the hessian matrix and return
    part1= lapply(1:nrow(data_yXZ), function(x) alpha.hessian.b(data_yXZ[x,]))
    part1 = Reduce("+", part1)/nrow(data_yXZ)
    return(part1)
  }
  my_A_alpha_inv = -solve(calculate_A_alpha())

  ## 2. The A_ipw of theta
  calculate_A_ipw = function(){

    theta = mybeta
    psi = mypsi

    A_ipw.b = function(data) {

      # X <- model.matrix(object = ~ AW+Z+Z2+Z3, data = data)
      # Y <- model.response(model.frame(formula = y~ AW+Z+Z2+Z3, data = data))
      # D <- data$D
      X <- model.matrix(object = ~ AW+Z, data = data)
      Y <- model.response(model.frame(formula = y~ AW+Z, data = data))
      D <- data$D

      # calculate the IPW score equation
      IPWscore <- function(myalpha.){
        # useful quantities
        mu = X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        myp = pi_xz_func.b(data.=data, myalpha. = myalpha.)
        myreturn = matrix(ncol=1, dotbeta*c(D/myp))
        return(myreturn)
      }

      # calculate the partial derivative
      myderiv = numDeriv::jacobian(function(x)
        IPWscore(myalpha. = x), myalpha, method = "Richardson")
      return(myderiv)
    }

    # calculate the hessian matrix and return
    part1= lapply(1:nrow(data_yXZ), function(x) A_ipw.b(data_yXZ[x,]))
    part1 = Reduce("+", part1)/nrow(data_yXZ)
    return(part1)
  }
  my_A_ipw = calculate_A_ipw()

  ## 3. A_alpha_theta
  calculate_A_alpha_theta = function(){

    theta = mybeta
    psi = mypsi

    # augmentation component
    phi_yz.b = function(data.){

      W <- model.matrix(object = y ~ A + Z, data = data.)
      X <- model.matrix(object = y ~ AW + Z, data = data.)
      Y <- model.response(model.frame(formula = y ~ AW + Z, data = data.))
      D <- data.$D
      z <- data.$Z

      # need these for integration components
      # meanXZ = data.$mymeanXZ

      #######################
      # Augmented Component #
      #######################
      psi_hat_i <- function(){

        # values needed
        e_star = Y - (W %*% theta)
        a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2)
        b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)

        # # denominator
        # mydenom = exp(c-(b^2/(4*a)))/sqrt(2*pi*(psi^2+theta[3]^2*sdXZ^2))

        ### beta0
        beta0_star = (theta[2]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
        ### beta2
        beta2_star = W[3]*beta0_star

        return(c(beta0_star, beta1_star, beta2_star)) #, beta3_star, beta4_star))
      }

      return(psi_hat_i())
    }

    # calculate the value for the ith observations
    A_alpha_theta.b = function(mydata){
      D = mydata$D
      myp = mydata$myp
      return(matrix(ncol=1, D*phi_yz.b(mydata)/(myp^2)) %*% t(alpha.jacobian.b(mydata)))
    }

    part1= lapply(1:nrow(data_yXZ), function(x) A_alpha_theta.b(data_yXZ[x,]))
    part1 = Reduce("+", part1) /nrow(data_yXZ)
    return(part1)
  }
  my_A_alpha_theta = calculate_A_alpha_theta()

  #####################################################
  # Step 4: Calculate Lambda using the values for theta

  # score function for alpha
  calculate_A = function(){

    theta = mybeta
    psi = mypsi

    ## calculate the bth component
    calculate_A.b = function(data, part1 = TRUE){

      # data = data_mvn() %>% subset(b==2)
      # theta = c(1,2,3)
      # W <- model.matrix(object = ~ A + Z + Z2 + Z3, data = data)
      # X <- model.matrix(object = y ~ AW + Z + Z2 + Z3, data = data)
      # Y <- model.response(model.frame(formula = y ~ AW + Z + Z2 + Z3, data = data))
      W <- model.matrix(object = ~ A + Z , data = data)
      X <- model.matrix(object = y ~ AW + Z, data = data)
      Y <- model.response(model.frame(formula = y ~ AW + Z, data = data))
      D <- data$D
      myp <- data$myp

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

        # values needed
        e_star = Y - (W %*% theta)
        a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2)
        b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)

        ### beta0
        beta0_star = (theta[2]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
        ### beta2
        beta2_star = W[3]*beta0_star
        # ### beta3
        # beta3_star = W[4]*beta0_star
        # ### beta4
        # beta4_star = W[5]*beta0_star

        return(c(beta0_star, beta1_star, beta2_star)) #, beta3_star, beta4_star))
      }

      # print(theta)
      if(part1 == TRUE){
        # myreturn = matrix(ncol = 1,D*CCscore()/myp) %*% t( matrix(ncol=1, (1-D/myp)*psi_hat_i())) +
        #   my_A_ipw%*%my_A_alpha_inv%*%jacobian.b(data) %*% t(matrix(ncol=1, (1-D/myp)*psi_hat_i()) + my_A_alpha_theta%*%my_A_alpha_inv%*%jacobian.b(data))
        myreturn = matrix(ncol=1, (1-D/myp)*psi_hat_i()) + my_A_alpha_theta%*%my_A_alpha_inv%*%alpha.jacobian.b(data)
        myreturn = (matrix(ncol = 1,D*CCscore()/myp) + my_A_ipw%*%my_A_alpha_inv%*%alpha.jacobian.b(data)) %*% t(myreturn)
        return(myreturn)
      }
      if(part1 == FALSE){
        myreturn = matrix(ncol=1, (1-D/myp)*psi_hat_i()) + my_A_alpha_theta%*%my_A_alpha_inv%*%alpha.jacobian.b(data)
        myreturn = myreturn%*%t(myreturn)
      }
    }

    # estimate first part of A matrix
    part1 = lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,]))
    part1 = Reduce("+", part1)/nrow(data_yXZ)

    part2= lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,], part1=FALSE))
    part2 = Reduce("+", part2)/nrow(data_yXZ)

    returnA = -part1%*%solve(part2)
    return(returnA)
  }
  myA = calculate_A()
  toc()

  #####################################################
  # Step 6: Calculate the A and the B matrices
  aipw.lambda.b <- function(theta, data.b){

    # W <- model.matrix(object = ~ A + Z + Z2 + Z3, data = data.b)
    # X <- model.matrix(object = y ~ AW + Z + Z2 + Z3, data = data.b)
    # Y <- model.response(model.frame(formula = y ~ AW + Z + Z2 + Z3, data = data.b))
    W <- model.matrix(object = ~ A + Z , data = data.b)
    X <- model.matrix(object = y ~ AW + Z, data = data.b)
    Y <- model.response(model.frame(formula = y ~ AW + Z, data = data.b))
    D <- data.b$D
    myp <- data.b$myp
    psi = mypsi

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

      # values needed
      e_star = Y - (W %*% theta)
      a = - 1/(2*sdXZ^2) -theta[2]^2/(2*psi^2)
      b = (meanXZ/sdXZ^2 -theta[2]*e_star/psi^2)
      c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
      mu_star = -b/(2*a)
      sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[2]^2*sdXZ^2)

      # # denominator
      # mydenom = exp(c-(b^2/(4*a)))/sqrt(2*pi*(psi^2+theta[3]^2*sdXZ^2))

      ### beta0
      beta0_star = (theta[2]*mu_star + e_star)/psi^2
      ### beta1
      beta1_star = (-theta[2]*(sd2_star+mu_star^2) + mu_star*(W[2]*theta[2] - e_star) + W[2]*e_star)/psi^2
      ### beta2
      beta2_star = W[3]*beta0_star
      ### beta3
      beta3_star = W[4]*beta0_star
      ### beta4
      beta4_star = W[5]*beta0_star

      return(c(beta0_star, beta1_star, beta2_star)) #, beta3_star, beta4_star))
    }

    # print(theta)
    AIPW_est = CCscore()*D/myp + (1 - D/myp)* myA %*% psi_hat_i() +
      (my_A_ipw +  myA %*% my_A_alpha_theta)%*%(my_A_alpha_inv%*%alpha.jacobian.b(data.b))
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
  return(list(beta_est = mytheta, se_est = sand.var))
}
