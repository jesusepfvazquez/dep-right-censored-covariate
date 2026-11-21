# winner!!
library(survival)
library(dplyr)

# X ~ Weibull
# C ~ Weiubull
data_aft_int <- function(nSubjects = 10^3, c_int = 1, p.c = 1){

  # # parameter values from function
  # c_int = 1
  # nSubjects = 300

  # defined parameters
  epsilon_sig = 1
  beta.vec = c(1,1,1)

  b             = 1:nSubjects # generate subject IDs;
  int           = cbind(rep(1,nSubjects))
  Z             = rnorm(nSubjects, mean = 0, sd = 1)
  A             = rnorm(nSubjects, mean = 2, sd = 1)

  # Generate X
  p.x = 0.5
  gamma.x.vec     = c(0.1,0.1)
  linPred_z       = cbind(int,Z)%*%gamma.x.vec
  errors.x.g      = log(rweibull(nSubjects, 1, 1))
  X               = exp(linPred_z + p.x*errors.x.g); mean(X); max(X)
  AX              = A-X

  er            = rnorm(nSubjects,mean=0,sd=epsilon_sig)

  ## compute linear predictor space and outcomes
  x             = cbind(int, AX, Z)
  linPred	      = x%*%beta.vec

  ## y is a normal random variable
  y             = linPred + er; mean(y)

  # generate C
  # p.c =1
  gamma.vec     = c(c_int, 0.1, 0.1)
  linPred_yz    = cbind(int, y, Z)%*%gamma.vec
  errors.g      = log(rweibull(nSubjects, 1, 1))
  C             = exp(linPred_yz + p.c*errors.g); mean(C); max(C)
  delta         = as.numeric(X <= C)
  W             = pmin(X, C)
  AW            = A-W
  mean(delta)

  # Create dataset
  dat           = as.data.frame(cbind(y,Z,delta,A,X,AX,C,W,AW,exp(linPred_yz),exp(linPred_z),b))
  colnames(dat) = c("y","Z","D","A","X","AX","C","W","AW","meanCYZ","meanXZ","b")

  #############################################################
  ##################### Calculate weights #####################
  #############################################################

  # random weights
  dat$myp_uniform = runif(nSubjects, min=0.1, max = 0.9)

  # oracle
  dat$myp_ywz_oracle = c(exp(-(dat$W /exp(linPred_yz))^(1/p.c))); summary(dat$myp_ywz_oracle)
  # c(exp(-(dat$W /exp(linPred_yz))^(1/p))) == pweibull(W, shape=1/p, scale = exp(linPred_yz), lower.tail = FALSE)
  # c(exp(-(dat$W /exp(linPred_yz))^(1/p))) == pgumbel(log(W), shape=p, scale = linPred_yz)
  
  # correct but estimated
  wr <- survreg(Surv(W, 1-D) ~ y + Z, data = dat, dist="w")
  lp <- predict(wr, dat)
  gamma <- 1 / wr$scale  # Shape parameter
  dat$myp_ywz = exp(-(dat$W /lp)^gamma)  # Survival probability at time t

  # Version 1
  get.meanXD0YZ = function(i){
    # i=5
    myresult1 = integrate(function(x) 
      x*dweibull(x, shape=1/p.x, scale = dat$meanXZ[i])*
        dnorm(dat$y[i] - (1*beta.vec[1] + (dat$A[i] - x)*beta.vec[2] + dat$Z[i]*beta.vec[3]), mean =0, sd = epsilon_sig),
      lower = dat$W[i], upper = Inf)
    myresult2 = integrate(function(x) 
      dweibull(x, shape=1/p.x, scale = dat$meanXZ[i])*
        dnorm(dat$y[i] - (1*beta.vec[1] + (dat$A[i] - x)*beta.vec[2] + dat$Z[i]*beta.vec[3]), mean =0, sd = epsilon_sig),
      lower = dat$W[i], upper = Inf)
    return(myresult1$value/myresult2$value)
  }
  
  dat$meanXD0YZ = lapply(1:nSubjects, get.meanXD0YZ) %>% unlist()
  dat$AX_yz = ifelse(dat$D == 1, dat$AX, dat$A - dat$meanXD0YZ)
  dat$AX_z = ifelse(dat$D == 1, dat$AX, dat$A - dat$meanXZ*gamma(1+1/p.x))
  
  return(dat)
}