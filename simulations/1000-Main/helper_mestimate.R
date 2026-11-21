# remove(list = ls())
library(dplyr)
library(geex)
library(tictoc)
library(cubature)

############################
## m-estimating equations ##
############################

# oracle or naive
estimate_beta <- function(data_yXZ, model){
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    
    function(theta){
      
      # first component of score
      myscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      return(myscore())
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ, weight = D)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  # print(g); print(psi)
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  # print(beta_estimates)
  # head(model.matrix(object = model, data = data_yXZ))
  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

# complete case
estimate_beta_cc <-  function(data_yXZ, model){
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    
    function(theta){
      
      # first component of score
      myscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(D*c(dotbeta))
      }
      
      return(myscore())
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ, weight = D)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

# ipw
estimate_beta_ipw <- function(data_yXZ, model, myweight = "oracle"){
  
  #######################################################
  # define weight components #
  #######################################################
  if (myweight  == "oracle") data_yXZ$myp = data_yXZ$myp_ywz_oracle
  if (myweight  == "aft") data_yXZ$myp = data_yXZ$myp_ywz
  if (myweight  == "uniform") data_yXZ$myp = data_yXZ$myp_uniform
  if (myweight  == "logit") data_yXZ$myp = data_yXZ$myp_ywz_logit
  
  #######################################################
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    myp <- data$myp
    
    function(theta){
      
      # first component of score
      CCscore <- function(){
        # useful quantities
        mu <- X %*% theta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        return(c(dotbeta))
      }
      
      return(c(D/myp)*CCscore())
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ, weights = D/myp)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

# aipw
estimate_beta_aipw <- function(data_yXZ, model, myweight = "oracle", eta1 = FALSE){
  
  # data_yXZ = dat
  
  ######################################################
  # define parameters values for integration component #
  ######################################################
  mybound = 50
  shape.c = 1/2
  shape.x = 1/0.8
  
  #################### Weight ####################
  
  data_yXZ$myp = data_yXZ$myp_ywz_oracle
  if (myweight  == "aft") data_yXZ$myp = data_yXZ$myp_ywz
  if (myweight  == "uniform") data_yXZ$myp = data_yXZ$myp_uniform
  if (myweight  == "logit") data_yXZ$myp = data_yXZ$myp_ywz_logit
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    # data = data_yXZ[1,]
    # W <- model.matrix(object = ~ A + Z, data = data)
    # X <- model.matrix(object = y ~ AW + Z, data = data)
    # Y <- model.response(model.frame(formula = y ~ AW + Z, data = data))
    # theta = c(-6,1,2); psi =10
    
    W <- model.matrix(object = ~ A + Z, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    z <- data$Z
    myp <- data$myp
    
    # need these for integration components
    meanXZ = data$meanXZ
    meanCYZ = data$meanCYZ
    
    # meanCXZ // needs to be estimated for each case of X
    
    function(theta){
      
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
              (1-1/pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE)) # (1-1/pi)
          }
          return(val)
        }
        
        ### Evaluate
        j_numerator_integrate <- function(jj=1) {
          integrate(function(y) {integral_func_num(t=y,j=jj)}, 
                    lower = 0, upper = mybound, rel.tol = 1e-7, abs.tol = 1e-7)$value}
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
          }
          return(val)
        }
        
        ### Evaluate
        denominator = integrate(integral_func_denom, 
                                lower = 0, upper = mybound, rel.tol = 1e-7, abs.tol = 1e-7)$value
        
        # Compute Augmented part 
        return(numerator/denominator)
      }
      
      if(eta1){
        
        psi_hat_i <- function(){
          
          # values needed
          meanXZ = 0
          sdXZ = 1
          
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
          
          return(c(beta0_star, beta1_star, beta2_star))
        }
      }
      # print(data$b)
      # print(psi_hat_i())
      return(CCscore()*D/myp + (1-D/myp)*psi_hat_i())
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ, weights = D/myp)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

estimate_beta_aipw_lambda_close <- function(data_yXZ, model, myweight = "oracle"){
  
  ######################################################
  # define parameters values for integration component #
  ######################################################
  # mybound = Inf
  # shape.c = 0.15
  # shape.x = 0.15
  # if (eta1) shape.x = 0.2
  sdXZ = var(data_yXZ$W[data_yXZ$D==1]) %>% sqrt()
  meanXZ = mean(data_yXZ$W[data_yXZ$D==1])
  
  #################### Weight ####################
  
  data_yXZ$myp = data_yXZ$myp_ywz_oracle
  if (myweight  == "aft") data_yXZ$myp = data_yXZ$myp_ywz
  if (myweight  == "uniform") data_yXZ$myp = data_yXZ$myp_uniform
  if (myweight  == "logit") data_yXZ$myp = data_yXZ$myp_ywz_logit
  
  ######################################################
  # Step 2: define initial estimates of theta using ipw
  mylmer  = lm(model, data = data_yXZ, weight = D/myp)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  #####################################################
  # Step 3: Calculate Lambda using the values for theta
  calculate_A = function(){
    
    theta = g
    
    ## step a: set data
    calculate_A.b = function(data, part1 = TRUE){
      
      # data = data_yXZ[1,]
      
      W <- model.matrix(object = ~ A + Z, data = data)
      X <- model.matrix(object = y ~ AW + Z, data = data)
      Y <- model.response(model.frame(formula = y ~ AW+Z, data = data))
      D <- data$D
      myp <- data$myp
      z <- data$Z
      
      # need these for integration components
      # sdXZ
      # meanXZ = data$meanXZ
      
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
      # augmentation component
      
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
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      # print(theta)
      if(part1 == TRUE){
        myreturn = psi_hat_i()
        myreturn = matrix(ncol = 1, c((1-D/myp)*(D/myp))*CCscore()) %*% psi_hat_i()
        return(myreturn)
      }
      if(part1 == FALSE){
        myreturn = psi_hat_i() %>% matrix(ncol = 1)
        myreturn = myreturn%*%t(myreturn)
        return(myreturn*c((1-D/myp)^2))
      }
    }
    
    # estimate first part of A matrix
    part1 = lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,], part1=TRUE))
    part1 = Reduce("+", part1) 
    
    part2 = lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,], part1=FALSE))
    part2 = Reduce("+", part2) 
    
    returnA = -part1%*%solve(part2)
    return(returnA)  
  }
  myA = calculate_A()
  # print("done with A"); print(myA)
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    W <- model.matrix(object = ~ A + Z, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    z <- data$Z
    myp <- data$myp
    
    # need these for integration components
    meanXZ = data$meanXZ
    # meanCYZ = data$meanCYZ
    
    # meanCXZ // needs to be estimated for each case of X
    
    function(theta){
      
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
      # augmentation component
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
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      return(c(D/myp)*CCscore() + c(1-D/myp)*myA%*%psi_hat_i())
    }
  }
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

# MLE
estimate_beta_likelihood <- function(data_yXZ, model, shape.x = 1/0.8 ,eta1=FALSE){
  
  # data_yXZ = dat
  # lowerbound= 0 #10^(-2)
  upperbound = Inf
  if (eta1){shape.x = 1/1.2}
  
  ################################
  # define eta1 misspecification #
  ################################
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    # data = data_yXZ[181,]
    # W <- model.matrix(object = ~ A + Z, data = data)
    # X <- model.matrix(object = y ~ AW +Z, data = data)
    # Y <- model.response(model.frame(formula = y ~ AW+Z, data = data))
    # theta = c(1.012072,    1.023458,    1.015727,    1.048713)
    
    W <- model.matrix(object = ~ A + Z, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    
    D <- data$D
    w <- data$W
    z <- data$Z
    
    # need these for integration components
    meanXZ =  data$meanXZ
    
    function(theta){
      
      # set theta
      beta = theta[1:3]
      psi = theta[4]
      
      # first component of score
      CCscore <- function(){
        # useful quantities
        mu <- X %*% beta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        # score equations for sigmas
        dotsigma = -(1/psi) + e^2/psi^3
        return(c(dotbeta,dotsigma))
      }
      
      # second component of score
      e_myscore <- function(){
        
        # top integral
        likelihood_int = function(t=1){
          return(dnorm(Y -  cbind(W[,1],W[,2]-t, W[,3]) %*% beta, 0, sd = psi))
        }
        
        # score of integral
        score_int <- function(t=1,j=1){
          # generate empty array of values to save
          # create temporary X matrix
          Wtemp <- cbind(W[,1],W[,2]-t, W[,3])
          mu <- Wtemp %*% beta
          e = Y - mu
          # score equations for betas
          dotbeta = (e/psi^2)%*%Wtemp
          # score equations for sigmas
          dotsigma = -(1/psi) + e^2/psi^3
          return(c(dotbeta,dotsigma)[j])
        }
        
        #top integral#
        integral_func_psi <- function(t=1,j=1){
          val = rep(NA, length(t))
          for (i in 1:length(t)){
            val[i] = score_int(t[i],j=j)*
              likelihood_int(t[i])*
              dweibull(t[i], shape=shape.x, scale = meanXZ)
          }
          return(val)
        }
        
        ## evaluate
        area_num <- function(j) cubintegrate(integral_func_psi,j=j,
                                             lower = w,
                                             upper = upperbound)$integral
        v.area_num <- Vectorize(area_num)
        numerator = v.area_num(1:4)
        
        # bottom integral
        integral_func_denom <- function(t=1){
          val = rep(NA, length(t))
          for (i in 1:length(t)){
            val[i] = likelihood_int(t[i])*
              dweibull(t[i], shape=shape.x, scale = meanXZ)
          }
          return(val)
        }
        # evaluate
        denominator = cubintegrate(function(x) integral_func_denom(t=x),
                                   lower = w,
                                   upper = upperbound)$integral
        
        # evaluate and return
        # print(numerator)
        return(numerator/denominator)
      }
      
      # print(theta); print(data$b)
      if(D==1) return(CCscore())
      if(D==0) return(e_myscore())
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(y ~ AX+Z, data = data_yXZ)
  # mylmer  = lm(model, data = data_yXZ, weights = D/myp_ywz_oracle)
  g = c(coef(mylmer), summary(mylmer)$sigma)
  print(g)
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model)
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% as.data.frame()
  colnames(beta_estimates) = c("","","","psi_updated")
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% as.data.frame()
  colnames(se_estimates) = c("","","")
  return(list(beta_est = beta_estimates, se_est = se_estimates))
}

# MACC
estimate_beta_macc_lambda_close <- function(data_yXZ, model, myweight = "oracle"){
  
  ######################################################
  sdXZ = var(data_yXZ$W[data_yXZ$D==1]) %>% sqrt()
  meanXZ = mean(data_yXZ$W[data_yXZ$D==1])
  
  # define weight components #
  data_yXZ$myp = data_yXZ$myp_ywz_oracle
  
  ######################################################
  # Step 2: define initial estimates of theta using ipw
  mylmer  = lm(model, data = data_yXZ, weight = D/myp)
  # mylmer  = lm(y ~ AW+Z+Z2+Z3, data = data_yXZ, weight = D)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  #####################################################
  # Step 3: Calculate Lambda using the values for theta
  calculate_A = function(){
    
    theta = g
    
    ## step a:set data
    calculate_A.b = function(data, part1 = TRUE){
      
      # data = data_mvn() %>% subset(b==2)
      # theta = c(1,2,3,1) +0.1
      W <- model.matrix(object = ~ Z + A, data = data)
      X <- model.matrix(object = ~ Z+AW, data = data)
      Y <- model.response(model.frame(formula = y~ Z+AW, data = data))
      D <- data$D
      myp <- data$myp
      z <- data$Z
      
      # need these for integration components
      meanXZ = data$meanXZ
      
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
        a = - 1/(2*sdXZ^2) -theta[3]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[3]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[3]^2*sdXZ^2)
        
        # # denominator
        # mydenom = exp(c-(b^2/(4*a)))/sqrt(2*pi*(psi^2+theta[3]^2*sdXZ^2))
        
        ### beta0
        beta0_star = (theta[3]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = W[2]*beta0_star
        ### beta2 
        beta2_star = (-theta[3]*(sd2_star+mu_star^2) + mu_star*(W[3]*theta[3] - e_star) + W[3]*e_star)/psi^2
        ### sigma
        # sigma_star = (theta[3]^2/psi^3)*(sd2_star+mu_star^2) + mu_star*(2*theta[3]*e_star/psi^3) + (e_star^2/psi^3 - 1/psi)
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      # print(theta)
      if(part1 == TRUE){
        myreturn = matrix(ncol = 1, (1-D/myp)*D*CCscore() ) %*% psi_hat_i()
        return(myreturn)
      }
      if(part1 == FALSE){
        myreturn = psi_hat_i()
        myreturn = (1-D/myp)^2*myreturn%*%t(myreturn)
      }
    }
    
    # estimate first part of A matrix
    part1 = lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,]))
    part1 = Reduce("+", part1) 
    
    part2= lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,], part1=FALSE))
    part2 = Reduce("+", part2) 
    
    returnA = -part1%*%solve(part2)
    return(returnA)  
  }
  myA = calculate_A()
  #print("done with A"); print(myA)
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    W <- model.matrix(object = ~ Z + A, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    z <- data$Z
    myp <- data$myp
    
    # need these for integration components
    meanXZ = data$meanXZ
    
    function(theta){
      
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
        a = - 1/(2*sdXZ^2) -theta[3]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[3]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[3]^2*sdXZ^2)
        
        # # denominator
        # mydenom = exp(c-(b^2/(4*a)))/sqrt(2*pi*(psi^2+theta[3]^2*sdXZ^2))
        
        ### beta0
        beta0_star = (theta[3]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = W[2]*beta0_star
        ### beta2 
        beta2_star = (-theta[3]*(sd2_star+mu_star^2) + mu_star*(W[3]*theta[3] - e_star) + W[3]*e_star)/psi^2
        ### sigma
        # sigma_star = (theta[3]^2/psi^3)*(sd2_star+mu_star^2) + mu_star*(2*theta[3]*e_star/psi^3) + (e_star^2/psi^3 - 1/psi)
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      # print(theta)
      ACC_est = CCscore()*D + (1 - D/myp)* myA%*%psi_hat_i()
      return(ACC_est)
    }
  }
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

# ACC
estimate_beta_acc_lambda_close <- function(data_yXZ, model, myweight = "oracle"){
  
  ######################################################
  sdXZ = var(data_yXZ$W[data_yXZ$D==1]) %>% sqrt()
  meanXZ = mean(data_yXZ$W[data_yXZ$D==1])
  
  # define weight components #
  data_yXZ$myp = data_yXZ$myp_ywz_oracle
  
  ######################################################
  # Step 2: define initial estimates of theta using ipw
  mylmer  = lm(model, data = data_yXZ, weight = D/myp)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  
  #####################################################
  # Step 3: Calculate Lambda using the values for theta
  calculate_A = function(){
    
    theta = g
    
    ## step a:set data
    calculate_A.b = function(data, part1 = TRUE){
      
      # data = data_mvn() %>% subset(b==2)
      # theta = c(1,2,3,1) +0.1
      W <- model.matrix(object = ~ Z + A, data = data)
      X <- model.matrix(object = ~ Z+AW, data = data)
      Y <- model.response(model.frame(formula = y~ Z+AW, data = data))
      D <- data$D
      myp <- data$myp
      z <- data$Z
      
      # need these for integration components
      meanXZ = data$meanXZ
      
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
        a = - 1/(2*sdXZ^2) -theta[3]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[3]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[3]^2*sdXZ^2)
        
        # # denominator
        # mydenom = exp(c-(b^2/(4*a)))/sqrt(2*pi*(psi^2+theta[3]^2*sdXZ^2))
        
        ### beta0
        beta0_star = (theta[3]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = W[2]*beta0_star
        ### beta2 
        beta2_star = (-theta[3]*(sd2_star+mu_star^2) + mu_star*(W[3]*theta[3] - e_star) + W[3]*e_star)/psi^2
        ### sigma
        # sigma_star = (theta[3]^2/psi^3)*(sd2_star+mu_star^2) + mu_star*(2*theta[3]*e_star/psi^3) + (e_star^2/psi^3 - 1/psi)
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      # print(theta)
      if(part1 == TRUE){
        myreturn = matrix(ncol = 1, (D-myp)*D*CCscore() ) %*% psi_hat_i()
        return(myreturn)
      }
      if(part1 == FALSE){
        myreturn = (D-myp)*psi_hat_i()
        myreturn = myreturn%*%t(myreturn)
      }
    }
    
    # estimate first part of A matrix
    part1 = lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,]))
    part1 = Reduce("+", part1) 
    
    part2= lapply(1:nrow(data_yXZ), function(x) calculate_A.b(data_yXZ[x,], part1=FALSE))
    part2 = Reduce("+", part2) 
    
    returnA = -part1%*%solve(part2)
    return(returnA)  
  }
  myA = calculate_A()
  #print("done with A"); print(myA)
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    W <- model.matrix(object = ~ Z + A, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    z <- data$Z
    myp <- data$myp
    
    # need these for integration components
    meanXZ = data$meanXZ
    
    function(theta){
      
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
        a = - 1/(2*sdXZ^2) -theta[3]^2/(2*psi^2) 
        b = (meanXZ/sdXZ^2 -theta[3]*e_star/psi^2)
        c = (-meanXZ^2/(2*sdXZ^2) - e_star^2/(2*psi^2))
        mu_star = -b/(2*a)
        sd2_star = (sdXZ^2*psi^2)/(psi^2 + theta[3]^2*sdXZ^2)
        
        # # denominator
        # mydenom = exp(c-(b^2/(4*a)))/sqrt(2*pi*(psi^2+theta[3]^2*sdXZ^2))
        
        ### beta0
        beta0_star = (theta[3]*mu_star + e_star)/psi^2
        ### beta1
        beta1_star = W[2]*beta0_star
        ### beta2 
        beta2_star = (-theta[3]*(sd2_star+mu_star^2) + mu_star*(W[3]*theta[3] - e_star) + W[3]*e_star)/psi^2
        ### sigma
        # sigma_star = (theta[3]^2/psi^3)*(sd2_star+mu_star^2) + mu_star*(2*theta[3]*e_star/psi^3) + (e_star^2/psi^3 - 1/psi)
        
        return(c(beta0_star, beta1_star, beta2_star))
      }
      
      # print(theta)
      ACC_est = CCscore()*D + (D-myp)*myA%*%psi_hat_i()
      return(ACC_est)
    }
  }
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'Richardson')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}


############################
## adjusting for variance ##
############################

estimate_beta_likelihood_est <- function(data_yXZ, model){
  
  # data_yXZ = dat
  # model = as.formula("y~AW+Z")
  ## estimate X|Z
  wr <- survreg(Surv(W, D) ~ Z, data = data_yXZ, dist="w")
  mygamma = c(coef(wr), wr$scale)
  
  ## set up upper bound
  upperbound = Inf
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    # data = data_yXZ[274,]
    # W <- model.matrix(object = ~ A + Z, data = data)
    # X <- model.matrix(object = y ~ AW +Z, data = data)
    # Y <- model.response(model.frame(formula = y ~ AW+Z, data = data))
    # theta = c(1.012072,    1.023458,    1.015727,    1.048713,
    #           1.2121799,   0.2059623,   1.1844931)
    
    W <- model.matrix(object = ~ A + Z, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    
    D <- data$D
    w <- data$W
    z <- data$Z
    
    # need these for integration components
    # meanXZ =  data$meanXZ
    
    function(theta){
      
      # loglikelihood 
      logLik.b = function(theta. = theta){
        
        # theta. = theta
        
        ###############################
        # set theta
        beta = theta.[1:3]
        psi = exp(theta.[4])
        # psi = theta.[4]
        gamma.x.vec = theta.[5:6]
        shape.x = theta.[7]
        
        # X|Z
        # linPred_z = 
        meanXZ = exp(cbind(1,z)%*%gamma.x.vec)
        
        ################################
        # delta = 1
        mu <- X %*% beta
        e = Y - mu
        
        delta1 = function() {
          # compute likelihood
          myreturn = dnorm(e, 0, sd = psi)*
            dweibull(w, shape=1/shape.x, scale = meanXZ)
          return(log(myreturn))
        }
        
        ################################
        # delta = 0
        delta0 = function() {
          
          # f|X,Z
          likelihood_int = function(t=1){
            return(dnorm(Y -  (cbind(W[,1],W[,2]-t, W[,3]) %*% beta), 0, sd = psi))
          }
          
          # integral
          integral_func_denom <- function(t=1){
            val = rep(NA, length(t))
            for (i in 1:length(t)){
              val[i] = likelihood_int(t[i])*
                dweibull(t[i], shape=1/shape.x, scale = meanXZ)
            }
            return(val)
          }
          
          # evaluate
          myreturn = cubintegrate(function(x) integral_func_denom(t=x),
                                  lower = w,
                                  upper = upperbound)$integral
          return(log(myreturn))
        }
        
        ##############################################
        if(D == 1) return(delta1())
        if(D == 0) return(delta0())
      }
      
      # logLik.b()
      
      # first derivative 
      # tic()
      myderiv = numDeriv::jacobian(function(x) 
        logLik.b(theta. = x), theta, method = "simple")
      # toc()
      # # return my derivative 
      # print(myderiv); print(data$b)
      return(myderiv)
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ, weights = D/myp_ywz)
  g = c(coef(mylmer), log(summary(mylmer)$sigma), mygamma)
  print(g)
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    compute_vcov = FALSE,
    root_control = setup_root_control(start = c(g)),
    deriv_control = setup_deriv_control(method = 'simple'),
    outer_args = list(formula = y ~ AW+Z)
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% as.data.frame()
  beta_estimates[4] = exp(beta_estimates[4])
  # colnames(beta_estimates) = c("","","","psi_updated")
  # se_estimates <- diag(results@vcov)^0.5 %>% t() %>% as.data.frame()
  # colnames(se_estimates) = c("","","")
  return(list(beta_est = beta_estimates))
}

estimate_beta_likelihood_optimx <- function(data_yXZ, model){
  
  # data_yXZ = dat
  # model = as.formula("y~AW+Z")
  
  data_yXZ$b = 1:nrow(data_yXZ)
  ######## X|Z ##########
  wr <- survreg(Surv(W, D) ~ Z, data = data_yXZ, dist="w")
  mygamma = c(coef(wr), wr$scale)
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(model, data = data_yXZ, weights = D/myp_ywz)
  g = c(coef(mylmer), log(summary(mylmer)$sigma), mygamma)
  print(g)
  
  ## set up upper bound
  upperbound = Inf
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(theta=g){
    
    mydata=data_yXZ
    # set theta
    # theta = g
    beta = theta[1:3]
    psi = exp(theta[4])
    gamma.x.vec = theta[5:6]
    shape.x = theta[7]
    
    Lik.b = function(data){
      
      # data = data_yXZ[1,]
      W <- model.matrix(object = ~ A + Z, data = data)
      X <- model.matrix(object = model, data = data)
      Y <- model.response(model.frame(formula = model, data = data))
      
      D <- data$D
      w <- data$W
      z <- data$Z
      
      # X|Z
      # linPred_z = 
      meanXZ = exp(cbind(1,z)%*%gamma.x.vec)
      
      ################################
      # delta = 1
      mu <- X %*% beta
      e = Y - mu
      
      delta1 = function() {
        # compute likelihood
        myreturn = dnorm(e, 0, sd = psi)*
          dweibull(w, shape=1/shape.x, scale = meanXZ)
        return(log(myreturn))
      }
      
      ################################
      # delta = 0
      delta0 = function() {
        
        # f|X,Z
        likelihood_int = function(t=1){
          return(dnorm(Y -  (cbind(W[,1],W[,2]-t, W[,3]) %*% beta), 0, sd = psi))
        }
        
        # integral
        integral_func_denom <- function(t=1){
          val = rep(NA, length(t))
          for (i in 1:length(t)){
            val[i] = likelihood_int(t[i])*
              dweibull(t[i], shape=1/shape.x, scale = meanXZ)
          }
          return(val)
        }
        
        # evaluate
        myreturn = integrate(function(x) integral_func_denom(t=x),
                             lower = w,
                             upper = upperbound)$value
        return(log(myreturn))
      }
      
      ##############################################
      if(D == 1) return(delta1())
      if(D == 0) return(delta0())
    }
    
    
    lapply(unique(mydata$b), function(x) Lik.b( mydata %>% filter(b==x))) %>%
      unlist() %>%
      sum() %>% return()
    
  }
  
  opt_result <- optimx::optimx(
    par = g,  # Initial values
    fn = gee_estfun, 
    method = "Nelder-Mead",  # Or try "BFGS", "L-BFGS-B" if gradient is available or bounds are needed
    control = list(
      trace = 1,       # Set to > 0 for more output
      maximize = TRUE,  # We are minimizing -loglik
      abstol = 1e-8,
      reltol = 1e-8
    ),
    hessian = FALSE 
  )
  
  
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- coef(opt_result) %>% unlist()
  beta_estimates[4] = log(beta_estimates[4])
  # colnames(beta_estimates) = c("","","","psi_updated")
  # se_estimates <- diag(results@vcov)^0.5 %>% t() %>% as.data.frame()
  # colnames(se_estimates) = c("","","")
  return(list(beta_est = beta_estimates))
}

var_beta_mle <- function(data_yXZ, theta){
  
  # data_yXZ = dat
  # theta = as.numeric(est63$beta_est)
  # theta = c(1,1,1,1,1, 0.1, 0.8)
  data_yXZ$b = 1:nrow(data_yXZ)
  upperbound = Inf
  
  #########################################################
  # define the stacked score equation
  b.mle <- function(data.b, d1 = TRUE){
    
    # loglikelihood 
    logLik.b = function(theta. = theta){
      
      # theta. = theta
      # data.b = dat[1,]
      
      W <- model.matrix(object = ~ A + Z, data = data.b)
      X <- model.matrix(object = y ~ AW+Z, data = data.b)
      Y <- model.response(model.frame(formula = y ~ AW+Z, data = data.b))
      
      D <- data.b$D
      w <- data.b$W
      z <- data.b$Z
      
      ###############################
      # set theta
      beta = theta.[1:3]
      psi = theta.[4]
      gamma.x.vec = theta.[5:6]
      shape.x = theta.[7]
      
      # X|Z
      meanXZ = exp(matrix(cbind(1,z),nrow = 1)%*%matrix(gamma.x.vec, ncol=1))
      
      ################################
      # delta = 1
      mu <- X %*% beta
      e = Y - mu
      
      delta1 = function() {
        # compute likelihood
        myreturn = dnorm(e, 0, sd = psi)*
          dweibull(w, shape=1/shape.x, scale = meanXZ)
        return(log(myreturn))
      }
      
      ################################
      # delta = 0
      delta0 = function() {
        
        # f|X,Z
        likelihood_int = function(t=1){
          return(dnorm(Y -  (cbind(W[,1],W[,2]-t, W[,3]) %*% beta), 0, sd = psi))
        }
        
        # integral
        integral_func_denom <- function(t=1){
          val = rep(NA, length(t))
          for (i in 1:length(t)){
            val[i] = likelihood_int(t[i])*
              dweibull(t[i], shape=1/shape.x, scale = meanXZ)
          }
          return(val)
        }
        
        # evaluate
        myreturn = cubintegrate(function(x) integral_func_denom(t=x),
                                lower = w,
                                upper = upperbound)$integral
        return(log(myreturn))
      }
      
      ##############################################
      if(D == 1) return(delta1())
      if(D == 0) return(delta0())
      
    }
    
    
    # first derivative 
    if(d1){
      myderiv = numDeriv::jacobian(function(x) 
        logLik.b(theta. = x), theta, method = "Richardson") %>%
        matrix(ncol=1)
    }
    
    # second derivative 
    if(d1==FALSE){
      myderiv = numDeriv::hessian(function(x) 
        logLik.b(theta. = x), theta, method = "Richardson")
    }
    
    return(myderiv)
  }
  # b.mle(dat[1,], d1=TRUE)
  
  # B matrix
  calculate.B = function(){
    part1 = lapply(1:nrow(data_yXZ), function(i){
      dotalpha = b.mle(data.b = data_yXZ[i,], d1 = TRUE)
      return(dotalpha %*% t(dotalpha))} 
    )
    part1 = Reduce("+", part1)
    return(part1)
  }
  myB = calculate.B()
  print("done Meat")
  # A matrix
  calculate.A = function(){
    part1 = lapply(1:nrow(data_yXZ), function(i) dotalpha = b.mle(data.b = data_yXZ[i,], d1 = FALSE))
    part1 = Reduce("+", part1) 
  }
  myA = calculate.A()
  
  print("done Bread")
  myA.inv = solve(myA)
  
  sand.var = myA.inv%*%myB%*%t(myA.inv)
  # return my sandwich variance estimator 
  
  return(list(beta_est = theta[1:4], se_est = sqrt(diag(sand.var))[1:3]))  
}

var_beta_ipw <- function(data_yXZ, mytheta){
  
  # data_yXZ = dat
  # mytheta = c(1,1,1,1)
  wr <- survreg(Surv(W, 1-D) ~ y + Z, data = data_yXZ, dist="w")
  myalpha = c(coef(wr), wr$scale)
  # myalpha = c(1, -0.5, 0.5, 2)
  mybeta = mytheta[1:(length(mytheta)-1)]
  mypsi = mytheta[length(mytheta)]
  myxi = c(mybeta,myalpha)
  
  #########################################################
  # alpha helper functions
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
  
  #########################################################
  # define the stacked score equation
  ipwscore.b <- function(myxi., data.b){
    
    # myxi. = myxi; data.b = dat[1,]
    mybeta. = myxi.[1:length(mybeta)]
    myalpha. = myxi.[(length(mybeta)+1):length(myxi.)]
    
    ######################
    ## IPW
    X <- model.matrix(object =  y~AW+Z, data = data.b)
    Y <- model.response(model.frame(formula =  y~AW+Z, data = data.b))
    myp = pi_xz_func.b(data.=data.b, myalpha. = myalpha.)
    D <- data.b$D
    
    mu <- X %*% mybeta.
    e = Y - mu
    # score equations for betas
    dotbeta = (e/mypsi^2)%*%X
    dotbeta = dotbeta*c(D/myp)
    
    ########################
    ## alpha 
    dotalpha = numDeriv::jacobian(function(x) 
      alpha.logLik.b(myalpha. = x, data. = data.b), myalpha., method = "Richardson")
    
    return(c(dotbeta,dotalpha))
  }
  
  # A matrix
  calculate.A = function(){
    part1 = lapply(1:nrow(data_yXZ), function(i) 
      dotalpha = numDeriv::jacobian(function(x) 
        ipwscore.b(myxi. = x, data.b = data_yXZ[i,]), myxi, method = "Richardson")
    )
    part1 = Reduce("+", part1) #/nrow(data_yXZ)
  }
  myA = calculate.A()
  myA.inv = solve(myA)
  print("done with bread")
  
  # B matrix
  calculate.B = function(){
    part1 = parallel::mclapply(1:nrow(data_yXZ), function(i){
      dotalpha = ipwscore.b(myxi. = myxi, data.b = data_yXZ[i,])
      return(dotalpha %*% t(dotalpha))} 
    )
    part1 = Reduce("+", part1)
    return(part1)
  }
  myB = calculate.B()
  print("done with meat")
  
  sand.var = myA.inv%*%myB%*%t(myA.inv)
  
  return(list(beta_est = mytheta, se_est = sqrt(diag(sand.var))[1:3]))  
}

estimate_beta_aipw_est <- function(data_yXZ, model, mytheta){
  
  # data_yXZ = dat
  # mytheta = c(unlist(est20$beta_est))
  # model = as.formula("y~AW+Z")
  
  ## estimate C|Y,Z
  wr <- survreg(Surv(W, 1-D) ~ y + Z, data = data_yXZ, dist="w")
  data_yXZ$meanCYZ = exp(wr$linear.predictors)
  data_yXZ$meanXZ = exp(mytheta[5] + mytheta[6]*data_yXZ$Z)
  
  ######################################################
  # define parameters values for integration component #
  ######################################################
  smallnum = 0 #1e-10
  lbound = 0
  mybound = 50
  shape.c = 1/wr$scale
  shape.x = 1/mytheta[length(mytheta)]
  
  #################### Weight ####################
  
  data_yXZ$myp = data_yXZ$myp_ywz
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    # formula = model
    # data = data_yXZ[56,]
    W <- model.matrix(object = ~ A + Z, data = data)
    # X <- model.matrix(object = y ~ AW+Z, data = data)
    # Y <- model.response(model.frame(formula = y ~ AW+Z, data = data))
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    z <- data$Z
    myp <- data$myp
    
    # need these for integration components
    meanXZ = data$meanXZ
    meanCYZ = data$meanCYZ
    
    # meanCXZ // needs to be estimated for each case of X
    
    function(theta){
      
      # theta = g #c(0.61066834, -0.09753506, -0.63524289)
      
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
              (1-1/pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE)) # (1-1/pi)
            # (1-1/(pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE)+smallnum)) # (1-1/pi)
          }
          return(val)
        }
        
        ### Evaluate
        j_numerator_integrate <- function(jj=1) {
          integrate(function(y) {integral_func_num(t=y,j=jj)}, 
                    lower = lbound, upper = mybound)$value}
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
            # (1-1/(pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE)+smallnum)) # (1-1/pi)
          }
          return(val)
        }
        
        ### Evaluate
        denominator = integrate(integral_func_denom, 
                                lower = lbound, upper = mybound)$value
        
        # Compute Augmented part
        # if (denominator==0) denominator = 1e-5
        
        return(numerator/denominator)
      }
      
      myreturn = CCscore()*D/myp + (1-D/myp)*psi_hat_i()
      # print(myreturn); print(data$b)
      return(myreturn)
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  # mylmer  = lm(model, data = data_yXZ, weights = D/myp)
  mylmer  = lm(model, data = data_yXZ, weights = D/myp)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  print(g)
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'simple')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

var_beta_aipw_est <- function(data_yXZ, mytheta, mygamma){
  
  # data_yXZ = dat
  # mytheta = c(1,1,1,1, 1, 0.1, 0.8)
  
  ##############################################################
  # Step 1: define parameters values for integration component #
  ##############################################################
  wr <- survreg(Surv(W, 1-D) ~ y + Z, data = data_yXZ, dist="w")
  myalpha = c(coef(wr), wr$scale)
  
  # mytheta = est32$beta_est
  mybeta = mytheta[1:3]
  mypsi = mytheta[4]
  data_yXZ$meanXZ = exp(mygamma[1] + mygamma[2]*data_yXZ$Z)
  data_yXZ$meanCYZ = exp(wr$linear.predictors)
  
  smallnum = 0 #1e-10
  lbound = 0
  mybound = 50
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
  pi_xz_func.b(data_yXZ[1,])
  
  jacobian_pi.b = function(data){
    myderiv = numDeriv::jacobian(function(x) 
      pi_xz_func.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(matrix(ncol=1, myderiv))
  } 
  jacobian_pi.b(data_yXZ[1,])
  
  alpha.logLik.b = function(data.=data_yXZ, myalpha. = myalpha){
    data.$W = log(data.$W)
    gamma.x.vec = myalpha.[1:3]
    shape.x = myalpha.[4]
    linPred_yz = cbind(rep(1,nrow(data.)),data.$y, data.$Z)%*%gamma.x.vec
    
    myreturn = (1-data.$D)*dgumbel(data.$W, shape=shape.x, scale = linPred_yz)+
      data.$D*pgumbel(data.$W, shape=shape.x, scale = linPred_yz)
    
    return(log(myreturn))
  }
  alpha.logLik.b(data_yXZ[1,])
  
  alpha.jacobian.b = function(data){
    myderiv = numDeriv::jacobian(function(x) 
      alpha.logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(matrix(ncol=1, myderiv))
  }
  alpha.jacobian.b(data_yXZ[1,])
  
  alpha.hessian.b = function(data){
    myderiv = numDeriv::hessian(function(x) 
      alpha.logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(myderiv)
  }
  alpha.hessian.b(data_yXZ[1,])
  
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
      
      # # X <- model.matrix(object = ~ AW+Z+Z2+Z3, data = data)
      # # Y <- model.response(model.frame(formula = y~ AW+Z+Z2+Z3, data = data))
      # # D <- data$D
      # X <- model.matrix(object = ~ AW+Z, data = data)
      # Y <- model.response(model.frame(formula = y~ AW+Z, data = data))
      # D <- data$D
      
      W <- model.matrix(object = ~ A + Z, data = data)
      X <- model.matrix(object =  y~ AW+Z, data = data)
      Y <- model.response(model.frame(formula =  y~ AW+Z, data = data))
      D <- data$D
      z <- data$Z
      # myp <- data$myp
      
      # need these for integration components
      meanXZ = data$meanXZ
      meanCYZ = data$meanCYZ
      
      # calculate the IPW score equation
      AIPWscore <- function(myalpha.){
        # # useful quantities
        # mu = X %*% theta
        # e = Y - mu
        # # score equations for betas
        # dotbeta = (e/psi^2)%*%X
        
        myp = pi_xz_func.b(data.=data, myalpha. = myalpha.)
        # myreturn = matrix(ncol=1, dotbeta*c(D/myp))
        # return(myreturn)
        
        
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
                (1-1/pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE)) # (1-1/pi)
              # (1-1/(smallnum+pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE))) # (1-1/pi)
            }
            return(val)
          }
          
          ### Evaluate
          j_numerator_integrate <- function(jj=1) {
            integrate(function(y) {integral_func_num(t=y,j=jj)}, 
                      lower = lbound, upper = mybound)$value}
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
                                  lower = lbound, upper = mybound)$value
          
          # Compute Augmented part 
          return(numerator/denominator)
        }
        
        
        myreturn = CCscore()*c(D/myp) + c(1-D/myp)*psi_hat_i()
        # myreturn = matrix(ncol=1, dotbeta*c(D/myp))
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
                  lower = lbound, upper = mybound)$value}
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
                              lower = lbound, upper = mybound)$value
      
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
  return(list(beta_est = mytheta, se_est = sand.var))
}

var_beta_aipw_lambda <- function(data_yXZ, mytheta){
  
  # data_yXZ = dat
  # mytheta = c(1,1,1,1)
  
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
      # W <- model.matrix(object = ~ A + Z + Z2 + Z3, data = data.)
      # X <- model.matrix(object = y ~ AW + Z + Z2 + Z3, data = data.)
      # Y <- model.response(model.frame(formula = y ~ AW + Z + Z2 + Z3, data = data.))
      # X <- model.matrix(object = ~ model, data = data.)
      # Y <- model.response(model.frame(formula = model, data = data.))
      D <- data.$D
      # myp <- data.$myp
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
        # ### beta3
        # beta3_star = W[4]*beta0_star
        # ### beta4
        # beta4_star = W[5]*beta0_star
        
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
    
    ## option 2
    # A_alpha_theta.b = function(data){
    #   
    #   # calculate the IPW score equation
    #   augscore <- function(myalpha.){
    #     # useful quantities
    #     myp = pi_xz_func.b(data.=data, myalpha. = myalpha.)
    #     myreturn = matrix(ncol=1, data$D*phi_yz.b(data)/myp)
    #     return(myreturn)
    #   }
    #   
    #   # calculate the partial derivative
    #   myderiv = numDeriv::jacobian(function(x) 
    #     augscore(myalpha. = x), myalpha, method = "Richardson")
    #   return(myderiv)
    # }
    
    
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
  print("A is done"); print(myA)
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

# other for AIPW

estimate_beta_likelihood_est_2 <- function(data_yXZ, model){
  
  # data_yXZ = dat
  # lowerbound= 0 #10^(-2)
  upperbound = 20
  # if (eta1){shape.x = 1/1.2}
  
  # X|Z
  wr2 <- survreg(Surv(W, D) ~ Z, data = data_yXZ, dist="w")
  mygamma = c(coef(wr2), wr2$scale)
  data_yXZ$meanXZ = exp(wr2$linear.predictors)
  shape.x = 1/wr2$scale
  
  ################################
  # define eta1 misspecification #
  ################################
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    W <- model.matrix(object = ~ A + Z, data = data)
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    
    D <- data$D
    w <- data$W
    z <- data$Z
    
    # need these for integration components
    meanXZ =  data$meanXZ
    
    function(theta){
      
      # set theta
      beta = theta[1:3]
      psi = theta[4]
      
      # first component of score
      CCscore <- function(){
        # useful quantities
        mu <- X %*% beta
        e = Y - mu
        # score equations for betas
        dotbeta = (e/psi^2)%*%X
        # score equations for sigmas
        dotsigma = -(1/psi) + e^2/psi^3
        return(c(dotbeta,dotsigma))
      }
      
      # second component of score
      e_myscore <- function(){
        
        # top integral
        likelihood_int = function(t=1){
          return(dnorm(Y -  cbind(W[,1],W[,2]-t, W[,3]) %*% beta, 0, sd = psi))
        }
        
        # score of integral
        score_int <- function(t=1,j=1){
          # generate empty array of values to save
          # create temporary X matrix
          Wtemp <- cbind(W[,1],W[,2]-t, W[,3])
          mu <- Wtemp %*% beta
          e = Y - mu
          # score equations for betas
          dotbeta = (e/psi^2)%*%Wtemp
          # score equations for sigmas
          dotsigma = -(1/psi) + e^2/psi^3
          return(c(dotbeta,dotsigma)[j])
        }
        
        #top integral#
        integral_func_psi <- function(t=1,j=1){
          val = rep(NA, length(t))
          for (i in 1:length(t)){
            val[i] = score_int(t[i],j=j)*
              likelihood_int(t[i])*
              dweibull(t[i], shape=shape.x, scale = meanXZ)
          }
          return(val)
        }
        
        ## evaluate
        area_num <- function(j) cubintegrate(integral_func_psi,j=j,
                                             lower = w,
                                             upper = upperbound)$integral
        v.area_num <- Vectorize(area_num)
        numerator = v.area_num(1:4)
        
        # bottom integral
        integral_func_denom <- function(t=1){
          val = rep(NA, length(t))
          for (i in 1:length(t)){
            val[i] = likelihood_int(t[i])*
              dweibull(t[i], shape=shape.x, scale = meanXZ)
          }
          return(val)
        }
        # evaluate
        denominator = cubintegrate(function(x) integral_func_denom(t=x),
                                   lower = w,
                                   upper = upperbound)$integral
        
        # evaluate and return
        # print(numerator)
        return(numerator/denominator)
      }
      
      # print(theta); print(data$b)
      if(D==1) return(CCscore())
      if(D==0) return(e_myscore())
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  mylmer  = lm(y ~ AX+Z, data = data_yXZ)
  # mylmer  = lm(model, data = data_yXZ, weights = D/myp_ywz_oracle)
  g = c(coef(mylmer), summary(mylmer)$sigma)
  print(g)
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'simple')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% as.data.frame()
  colnames(beta_estimates) = c("","","","psi_updated")
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% as.data.frame()
  colnames(se_estimates) = c("","","")
  return(list(beta_est = beta_estimates, se_est = se_estimates))
}

var_beta_mle_2 <- function(data_yXZ, theta){
  
  # data_yXZ = dat
  # theta = as.numeric(est63$beta_est)
  # theta = c(1,1,1,1,1, 0.1, 0.8)
  wr2 <- survreg(Surv(W, D) ~ Z, data = data_yXZ, dist="w")
  mygamma = c(coef(wr2), wr2$scale)
  theta = c(theta, mygamma)
  data_yXZ$b = 1:nrow(data_yXZ)
  upperbound = 30
  
  #########################################################
  # define the stacked score equation
  b.mle <- function(data.b, d1 = TRUE){
    
    # loglikelihood 
    logLik.b = function(theta. = theta){
      
      # theta. = theta
      # data.b = dat[1,]
      
      W <- model.matrix(object = ~ A + Z, data = data.b)
      X <- model.matrix(object = y ~ AW+Z, data = data.b)
      Y <- model.response(model.frame(formula = y ~ AW+Z, data = data.b))
      
      D <- data.b$D
      w <- data.b$W
      z <- data.b$Z
      
      ###############################
      # set theta
      beta = theta.[1:3]
      psi = theta.[4]
      gamma.x.vec = theta.[5:6]
      shape.x = theta.[7]
      
      # X|Z
      meanXZ = exp(matrix(cbind(1,z),nrow = 1)%*%matrix(gamma.x.vec, ncol=1))
      
      ################################
      # delta = 1
      mu <- X %*% beta
      e = Y - mu
      
      delta1 = function() {
        # compute likelihood
        myreturn = dnorm(e, 0, sd = psi)*
          dweibull(w, shape=1/shape.x, scale = meanXZ)
        return(log(myreturn))
      }
      
      ################################
      # delta = 0
      delta0 = function() {
        
        # f|X,Z
        likelihood_int = function(t=1){
          return(dnorm(Y -  (cbind(W[,1],W[,2]-t, W[,3]) %*% beta), 0, sd = psi))
        }
        
        # integral
        integral_func_denom <- function(t=1){
          val = rep(NA, length(t))
          for (i in 1:length(t)){
            val[i] = likelihood_int(t[i])*
              dweibull(t[i], shape=1/shape.x, scale = meanXZ)
          }
          return(val)
        }
        
        # evaluate
        myreturn = cubintegrate(function(x) integral_func_denom(t=x),
                                lower = w,
                                upper = upperbound)$integral
        return(log(myreturn))
      }
      
      ##############################################
      if(D == 1) return(delta1())
      if(D == 0) return(delta0())
      
    }
    
    
    # first derivative 
    if(d1){
      myderiv = numDeriv::jacobian(function(x) 
        logLik.b(theta. = x), theta, method = "Richardson") %>%
        matrix(ncol=1)
    }
    
    # second derivative 
    if(d1==FALSE){
      myderiv = numDeriv::hessian(function(x) 
        logLik.b(theta. = x), theta, method = "Richardson")
    }
    
    return(myderiv)
  }
  # b.mle(dat[1,], d1=TRUE)
  
  # B matrix
  calculate.B = function(){
    part1 = lapply(1:nrow(data_yXZ), function(i){
      dotalpha = b.mle(data.b = data_yXZ[i,], d1 = TRUE)
      return(dotalpha %*% t(dotalpha))} 
    )
    part1 = Reduce("+", part1)
    return(part1)
  }
  myB = calculate.B()
  print("done Meat")
  # A matrix
  calculate.A = function(){
    part1 = lapply(1:nrow(data_yXZ), function(i) dotalpha = b.mle(data.b = data_yXZ[i,], d1 = FALSE))
    part1 = Reduce("+", part1) 
  }
  myA = calculate.A()
  
  print("done Bread")
  myA.inv = solve(myA)
  
  sand.var = myA.inv%*%myB%*%t(myA.inv)
  # return my sandwich variance estimator 
  
  return(list(beta_est = theta[1:4], se_est = sqrt(diag(sand.var))[1:3]))  
}

estimate_beta_aipw_est_2 <- function(data_yXZ, model){
  
  # data_yXZ = dat
  # mytheta = c(unlist(est10$beta_est))
  # model = as.formula("y~AW+Z")
  
  # C|Y,Z
  wr <- survreg(Surv(W, 1-D) ~ y + Z, data = data_yXZ, dist="w")
  myalpha = c(coef(wr), wr$scale)
  data_yXZ$meanCYZ = exp(wr$linear.predictors)
  
  lp <- predict(wr, data_yXZ)
  gamma <- 1 / wr$scale  # Shape parameter
  data_yXZ$myp_ywz = exp(-(data_yXZ$W /lp)^gamma)  # Survival probability at time t
  
  # X|Z
  wr2 <- survreg(Surv(W, D) ~ Z, data = data_yXZ, dist="w")
  mygamma = c(coef(wr2), wr2$scale)
  data_yXZ$meanXZ = exp(wr2$linear.predictors)
  
  ######################################################
  # define parameters values for integration component #
  ######################################################
  smallnum = 0 # 1e-10
  lbound = 0
  mybound = 20
  shape.c = 1/wr$scale
  shape.x = 1/wr2$scale
  
  #################### Weight ####################
  
  data_yXZ$myp = data_yXZ$myp_ywz
  
  ########  DEFINE ESTIMATING FUNCTION ##########
  gee_estfun <- function(data, formula){
    
    # formula = model
    # data = data_yXZ[56,]
    W <- model.matrix(object = ~ A + Z, data = data)
    # X <- model.matrix(object = y ~ AW+Z, data = data)
    # Y <- model.response(model.frame(formula = y ~ AW+Z, data = data))
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))
    D <- data$D
    z <- data$Z
    myp <- data$myp
    
    # need these for integration components
    meanXZ = data$meanXZ
    meanCYZ = data$meanCYZ
    
    # meanCXZ // needs to be estimated for each case of X
    
    function(theta){
      
      # theta = g #c(0.61066834, -0.09753506, -0.63524289)
      
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
              # (1-1/pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE)) # (1-1/pi)
              (1-1/(pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE)+smallnum)) # (1-1/pi)
          }
          return(val)
        }
        
        ### Evaluate
        j_numerator_integrate <- function(jj=1) {
          integrate(function(y) {integral_func_num(t=y,j=jj)}, 
                    lower = lbound, upper = mybound)$value}
        v.area_num <- Vectorize(j_numerator_integrate)
        numerator = v.area_num(1:3)
        
        # bottom integral
        integral_func_denom <- function(t=1){
          val = rep(NA, length(t))
          for (i in 1:length(t)){
            val[i] = likelihood_int(t[i])*
              # pweibull(C, shape=gamma, scale = lp, lower.tail = FALSE)
              dweibull(t[i], shape=shape.x, scale = meanXZ)*
              # (1-1/pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE)) # (1-1/pi)
              (1-1/(pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE)+smallnum)) # (1-1/pi)
          }
          return(val)
        }
        
        ### Evaluate
        denominator = integrate(integral_func_denom, 
                                lower = lbound, upper = mybound)$value
        
        # Compute Augmented part
        # if (denominator==0) denominator = 1e-5
        
        return(numerator/denominator)
      }
      
      myreturn = CCscore()*D/myp + (1-D/myp)*psi_hat_i()
      # print(myreturn); print(data$b)
      return(myreturn)
    }
  }
  
  ######## EST INITIAL PARAMETER ESTIMATES ##########
  # mylmer  = lm(model, data = data_yXZ, weights = D/myp)
  mylmer  = lm(model, data = data_yXZ, weights = D/myp)
  g = c(coef(mylmer))
  psi = summary(mylmer)$sigma
  print(g)
  
  ######## RUNNING M ESTIMATE FUNCTION ##########
  results <- m_estimate(
    estFUN = gee_estfun,
    data = data_yXZ,
    root_control = setup_root_control(start = c(g)),
    outer_args = list(formula = model),
    deriv_control = setup_deriv_control(method = 'simple')
  )
  
  ######## RETURN ESTIMATES ##########
  beta_estimates <- results@estimates %>% t() %>% matrix()
  data_yXZ = data_yXZ %>% subset(D==1)  
  psi_updated = sum((data_yXZ$y - model.matrix(object = model, data = data_yXZ)%*%beta_estimates)^2/(nrow(data_yXZ))) %>% sqrt
  beta_estimates = rbind(beta_estimates, psi_updated)
  se_estimates <- diag(results@vcov)^0.5 %>% t() %>% matrix()
  return(list(beta_est = t(beta_estimates), se_est = t(se_estimates)))
}

var_beta_aipw_est_2 <- function(data_yXZ, mytheta){
  
  # data_yXZ = dat
  # mytheta = c(unlist(est20$beta_est))
  # model = as.formula("y~AW+Z")
  
  ##############################################################
  # Step 1: define parameters values for integration component #
  ##############################################################
  # wr <- survreg(Surv(W, 1-D) ~ y + Z, data = data_yXZ, dist="w")
  # myalpha = c(coef(wr), wr$scale)
  # 
  # # mytheta = est32$beta_est
  # mybeta = mytheta[1:3]
  # mypsi = mytheta[4]
  # data_yXZ$meanXZ = exp(mygamma[1] + mygamma[2]*data_yXZ$Z)
  # data_yXZ$meanCYZ = exp(wr$linear.predictors)
  # 
  # smallnum = 0 #1e-10
  # lbound = 0
  # mybound = 50
  # shape.c = 1/wr$scale
  # shape.x = 1/mygamma[3]
  # data_yXZ$myp = data_yXZ$myp_ywz
  
  # C|Y,Z
  wr <- survreg(Surv(W, 1-D) ~ y + Z, data = data_yXZ, dist="w")
  data_yXZ$meanCYZ = exp(wr$linear.predictors)
  
  lp <- predict(wr, data_yXZ)
  gamma <- 1 / wr$scale  # Shape parameter
  data_yXZ$myp = exp(-(data_yXZ$W /lp)^gamma)  # Survival probability at time t
  
  # X|Z
  wr2 <- survreg(Surv(W, D) ~ Z, data = data_yXZ, dist="w")
  data_yXZ$meanXZ = exp(wr2$linear.predictors)
  
  # mytheta = est32$beta_est
  mybeta = mytheta[1:3]
  mypsi = mytheta[4]
  myalpha = c(coef(wr), wr$scale)
  mygamma = c(coef(wr2), wr2$scale)
  
  
  smallnum = 0 #1e-10
  lbound = 0
  mybound = 20
  shape.c = 1/wr$scale
  shape.x = 1/wr2$scale
  
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
  pi_xz_func.b(data_yXZ[1,])
  
  jacobian_pi.b = function(data){
    myderiv = numDeriv::jacobian(function(x) 
      pi_xz_func.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(matrix(ncol=1, myderiv))
  } 
  jacobian_pi.b(data_yXZ[1,])
  
  alpha.logLik.b = function(data.=data_yXZ, myalpha. = myalpha){
    data.$W = log(data.$W)
    gamma.x.vec = myalpha.[1:3]
    shape.x = myalpha.[4]
    linPred_yz = cbind(rep(1,nrow(data.)),data.$y, data.$Z)%*%gamma.x.vec
    
    myreturn = (1-data.$D)*dgumbel(data.$W, shape=shape.x, scale = linPred_yz)+
      data.$D*pgumbel(data.$W, shape=shape.x, scale = linPred_yz)
    
    return(log(myreturn))
  }
  alpha.logLik.b(data_yXZ[1,])
  
  alpha.jacobian.b = function(data){
    myderiv = numDeriv::jacobian(function(x) 
      alpha.logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(matrix(ncol=1, myderiv))
  }
  alpha.jacobian.b(data_yXZ[1,])
  
  alpha.hessian.b = function(data){
    myderiv = numDeriv::hessian(function(x) 
      alpha.logLik.b(myalpha. = x, data. = data), myalpha, method = "Richardson")
    return(myderiv)
  }
  alpha.hessian.b(data_yXZ[1,])
  
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
      
      # # X <- model.matrix(object = ~ AW+Z+Z2+Z3, data = data)
      # # Y <- model.response(model.frame(formula = y~ AW+Z+Z2+Z3, data = data))
      # # D <- data$D
      # X <- model.matrix(object = ~ AW+Z, data = data)
      # Y <- model.response(model.frame(formula = y~ AW+Z, data = data))
      # D <- data$D
      
      W <- model.matrix(object = ~ A + Z, data = data)
      X <- model.matrix(object =  y~ AW+Z, data = data)
      Y <- model.response(model.frame(formula =  y~ AW+Z, data = data))
      D <- data$D
      z <- data$Z
      # myp <- data$myp
      
      # need these for integration components
      meanXZ = data$meanXZ
      meanCYZ = data$meanCYZ
      
      # calculate the IPW score equation
      AIPWscore <- function(myalpha.){
        # # useful quantities
        # mu = X %*% theta
        # e = Y - mu
        # # score equations for betas
        # dotbeta = (e/psi^2)%*%X
        
        myp = pi_xz_func.b(data.=data, myalpha. = myalpha.)
        # myreturn = matrix(ncol=1, dotbeta*c(D/myp))
        # return(myreturn)
        
        
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
                (1-1/pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE)) # (1-1/pi)
              # (1-1/(smallnum+pweibull(t[i], shape=shape.c, scale = meanCYZ, lower.tail = FALSE))) # (1-1/pi)
            }
            return(val)
          }
          
          ### Evaluate
          j_numerator_integrate <- function(jj=1) {
            integrate(function(y) {integral_func_num(t=y,j=jj)}, 
                      lower = lbound, upper = mybound)$value}
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
                                  lower = lbound, upper = mybound)$value
          
          # Compute Augmented part 
          return(numerator/denominator)
        }
        
        
        myreturn = CCscore()*c(D/myp) + c(1-D/myp)*psi_hat_i()
        # myreturn = matrix(ncol=1, dotbeta*c(D/myp))
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
                  lower = lbound, upper = mybound)$value}
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
                              lower = lbound, upper = mybound)$value
      
      # Compute Augmented part 
      return(numerator/denominator)
    }
    
    # print(theta)
    AIPW_est = CCscore()*D/myp + (1 - D/myp)*psi_hat_i() + 
      my_A_aipw%*%(my_A_alpha_inv%*%alpha.jacobian.b(data.b))
    # print(ACC_est)
    return(AIPW_est)
    
  }
  aipw.lambda.b(mybeta, data_yXZ[1,])
  
  #####################################################
  # Step 6: Calculate the A and the B matrices
  
  calculate.B = function(){
    part1 = lapply(1:nrow(data_yXZ), function(i){
      # i=1
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


#####################################################

simcalc_theta_efficiency = function(i=1, n, c_int = 0, p.c = 1.5, sim_seed){
  
  # n=100
  # c_int = 0
  # i=1
  
  set.seed(sim_seed)
  
  # empty matrix to save results
  mycoef = matrix(nrow = 10, ncol=8)
  dat <- data_aft_int(nSubjects = n, c_int, p.c)
  print(mean(dat$D))
  
  safe <- function(expr) tryCatch(expr, error = function(e) {
    message("FAILED seed=", sim_seed, " c_int=", c_int, " n=", n, ": ", conditionMessage(e))
    list(beta_est = rep(NA, 4), se_est = rep(NA, 3))
  })
  
  est00 = safe(estimate_beta(dat, y ~ AX + Z))
  est01 = safe(estimate_beta(dat, y ~ AW + Z))
  est02 = safe(estimate_beta_cc(dat, y ~ AW + Z))
  
  # IPW
  print("Start IPW")
  est10 = safe(estimate_beta_ipw(dat, y ~ AW + Z, myweight = "aft"))
  est11 = safe(var_beta_ipw(dat, est10$beta_est))
  
  # MLE
  print("Start MLE")
  est20 = safe(estimate_beta_likelihood_est_2(dat, y~AW+Z))
  est21 = safe(var_beta_mle_2(dat, as.numeric(est20$beta_est)))
  
  # AIPW normal
  print("Start AIPW")
  est30 = safe(estimate_beta_aipw_est_2(dat, y ~ AW + Z))
  est31 = safe(var_beta_aipw_est_2(dat, mytheta = c(est30$beta_est)))
  
  # AIPW Lambda
  est40 = safe(estimate_beta_aipw_lambda_close(dat, y ~ AW + Z, myweight = "aft"))
  est41 = safe(var_beta_aipw_lambda(dat, est40$beta_est))
  
  # save results
  mycoef[,1] = c("oracle", "naive", "cc", "ipw", "ipw-est", "mle-est", "aipw","aipw-est", "aipw-lambda",  "aipw-lambda-est")
  mycoef[,2:5] = rbind(est00$beta_est, est01$beta_est, est02$beta_est, est10$beta_est, est11$beta_est,
                       est21$beta_est, est30$beta_est, est31$beta_est, est40$beta_est, est41$beta_est)
  mycoef[,6:8] = rbind(est00$se_est, est01$se_est, est02$se_est, est10$se_est, est11$se_est,
                       est21$se_est, est30$se_est, est31$se_est, est40$se_est, est41$se_est)
  mycoef = data.frame(mycoef)
  colnames(mycoef) = c("Method","beta0", "beta1", "beta2","sigma", "sebeta0", "sebeta1", "sebeta2")
  
  # compute probability of censoring
  mycoef$pcen = 1-mean(dat$D)
  mycoef$group = c_int
  mycoef$sim = sim_seed
  mycoef$n = n
  
  write.csv(mycoef, paste0("results/est_results_n_",n, "_",c_int,"_seed_", sim_seed, ".csv"), row.names = F)
  
  return(mycoef)
}

simcalc_nuisance = function(n, c_int = 0, p.c = 1.5, sim_seed){
  
  set.seed(sim_seed)
  
  # empty matrix to save results
  dat <- data_aft_int(nSubjects = n, c_int, p.c)
  print(mean(dat$D))
  
  # C|Y,Z
  mycoef = matrix(nrow = 1, ncol=8)
  wr <- survreg(Surv(W, 1-D) ~ y + Z, data = dat, dist="w")
  mycoef[1,1:4] = c(coef(wr), wr$scale)
  mycoef[1,5:8] = summary(wr)$var %>% diag() %>% sqrt()
  mycoef = data.frame(mycoef)
  colnames(mycoef) = c("alpha0", "alpha1", "alpha2", "alpha3", 
                       "alpha0se", "alpha1se", "alpha2se", "alpha3se")            
  # compute probability of censoring
  mycoef$pcen = 1-mean(dat$D)
  mycoef$group = c_int
  mycoef$sim = sim_seed
  mycoef$n = n
  
  
  write.csv(mycoef, paste0("results/est_results_alpha_n_",n, "_",c_int,"_seed_", sim_seed, ".csv"), row.names = F)
  
  # X|Z
  mycoef = matrix(nrow = 1, ncol=6)
  wr <- survreg(Surv(W, D) ~ Z, data = dat, dist="w")
  mycoef[1,1:3] = c(coef(wr), wr$scale)
  mycoef[1,4:6] = summary(wr)$var %>% diag() %>% sqrt()
  mycoef = data.frame(mycoef)
  colnames(mycoef) = c("alpha0", "alpha1", "alpha2", 
                       "alpha0se", "alpha1se", "alpha2se")            
  # compute probability of censoring
  mycoef$pcen = 1-mean(dat$D)
  mycoef$group = c_int
  mycoef$sim = sim_seed
  mycoef$n = n
  
  write.csv(mycoef, paste0("results/est_results_gamma_n_",n, "_",c_int,"_seed_", sim_seed, ".csv"), row.names = F)
  
  return(mycoef)
}