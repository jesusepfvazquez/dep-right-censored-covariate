# winner!!
library(survival)
library(dplyr)

## Solve for the censoring intercept that yields target censoring at a given xi.
## Mirrors the data-generating constants in data_aft_int exactly.
calibrate_c_int <- function(xi, target_pcen, n_cal = 1e5, p.c = 1.5, seed = 12345){
  set.seed(seed)
  epsilon_sig = 1
  beta.vec    = c(1,1,1)
  int         = rep(1, n_cal)
  Z           = rnorm(n_cal, 0, 1)
  A           = rnorm(n_cal, 2, 1)
  
  p.x         = 0.5
  gamma.x.vec = c(0.1, 0.1)
  linPred_z   = cbind(int, Z) %*% gamma.x.vec
  errors.x.g  = log(rweibull(n_cal, 1, 1))
  X           = exp(linPred_z + p.x * errors.x.g)
  AX          = A - X
  er          = rnorm(n_cal, 0, epsilon_sig)
  y           = cbind(int, AX, Z) %*% beta.vec + er
  Xc          = X - mean(X)
  
  errors.g    = log(rweibull(n_cal, 1, 1))   # fixed across the root search
  
  ## censoring proportion as a function of the intercept ci, at this xi
  pcen_at <- function(ci){
    linPred_yzx = cbind(int, y, Z) %*% c(ci, 0.1, 0.1) + xi * Xc
    C  = exp(linPred_yzx + p.c * errors.g)
    mean(X > C)
  }
  
  uniroot(function(ci) pcen_at(ci) - target_pcen,
          interval = c(-10, 10), tol = 1e-4)$root
}

data_aft_int <- function(nSubjects = 10^3, c_int = 1, p.c = 1, xi = 0){
  
  epsilon_sig = 1
  beta.vec = c(1,1,1)
  
  b             = 1:nSubjects
  int           = cbind(rep(1,nSubjects))
  Z             = rnorm(nSubjects, mean = 0, sd = 1)
  A             = rnorm(nSubjects, mean = 2, sd = 1)
  
  p.x = 0.5
  gamma.x.vec     = c(0.1,0.1)
  linPred_z       = cbind(int,Z)%*%gamma.x.vec
  errors.x.g      = log(rweibull(nSubjects, 1, 1))
  X               = exp(linPred_z + p.x*errors.x.g)
  AX              = A-X
  
  er            = rnorm(nSubjects,mean=0,sd=epsilon_sig)
  x             = cbind(int, AX, Z)
  linPred       = x%*%beta.vec
  y             = linPred + er
  
  # generate C  -- xi*(X - mean(X)) breaks C _||_ X | (Y,Z); xi=0 recovers original
  gamma.vec     = c(c_int, 0.1, 0.1)
  Xc            = X - mean(X)
  linPred_yz    = cbind(int, y, Z)%*%gamma.vec          # (Y,Z)-only predictor (used by weights)
  linPred_yzx   = linPred_yz + xi*Xc                    # true predictor incl. X dependence
  errors.g      = log(rweibull(nSubjects, 1, 1))
  C             = exp(linPred_yzx + p.c*errors.g)
  delta         = as.numeric(X <= C)
  W             = pmin(X, C)
  AW            = A-W
  
  dat           = as.data.frame(cbind(y,Z,delta,A,X,AX,C,W,AW,exp(linPred_yz),exp(linPred_z),b))
  colnames(dat) = c("y","Z","D","A","X","AX","C","W","AW","meanCYZ","meanXZ","b")
  
  #############################################################
  ##################### Calculate weights #####################
  #############################################################
  dat$myp_uniform = runif(nSubjects, min=0.1, max = 0.9)
  
  # (Y,Z)-oracle: uses linPred_yz WITHOUT the xi term, as a blind analyst would
  dat$myp_ywz_oracle = c(exp(-(dat$W /exp(linPred_yz))^(1/p.c)))
  
  # correct-but-estimated weights, fit on (Y,Z) only
  wr <- survreg(Surv(W, 1-D) ~ y + Z, data = dat, dist="w")
  lp <- predict(wr, dat)
  gamma <- 1 / wr$scale
  dat$myp_ywz = exp(-(dat$W /lp)^gamma)
  
  # Version 1
  get.meanXD0YZ = function(i){
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
  dat$AX_z  = ifelse(dat$D == 1, dat$AX, dat$A - dat$meanXZ*gamma(1+1/p.x))
  
  return(dat)
}