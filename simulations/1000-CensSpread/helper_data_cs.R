## helper_data_cs.R
## data_aft_int() from 1000-Main/helper_data.R (the Figure 1 design), without the
## E[X | D = 0, Y, Z] integrals (only used by the imputation estimators). Same random draws
## in the same order.
##   log X = 0.1 + 0.1 Z + 0.5 log Exp(1);  Y = 1 + (A - X) + Z + e, e ~ N(0, 1)
##   log C = c_int + 0.1 Y + 0.1 Z + p.c log Exp(1)
## p.c is the scale (spread) of log C; c_int sets the censoring rate.
library(survival)

data_cs <- function(nSubjects = 1000, c_int = 1, p.c = 1) {
  int <- rep(1, nSubjects)
  Z   <- rnorm(nSubjects, 0, 1)
  A   <- rnorm(nSubjects, 2, 1)
  linPred_z <- 0.1 + 0.1 * Z
  X   <- exp(linPred_z + 0.5 * log(rweibull(nSubjects, 1, 1)))
  AX  <- A - X
  y   <- 1 + AX + Z + rnorm(nSubjects, 0, 1)
  linPred_yz <- c_int + 0.1 * y + 0.1 * Z
  C   <- exp(linPred_yz + p.c * log(rweibull(nSubjects, 1, 1)))
  D   <- as.numeric(X <= C); W <- pmin(X, C)
  dat <- data.frame(y = y, Z = Z, D = D, A = A, X = X, AX = AX, C = C, W = W, AW = A - W,
                    meanCYZ = exp(linPred_yz), meanXZ = exp(linPred_z), b = seq_len(nSubjects))
  dat$myp_uniform    <- runif(nSubjects, 0.1, 0.9)
  dat$myp_ywz_oracle <- exp(-(dat$W / exp(linPred_yz))^(1 / p.c))
  wr <- survreg(Surv(W, 1 - D) ~ y + Z, data = dat, dist = "w")
  dat$myp_ywz <- exp(-(dat$W / predict(wr, dat))^(1 / wr$scale))
  dat
}
