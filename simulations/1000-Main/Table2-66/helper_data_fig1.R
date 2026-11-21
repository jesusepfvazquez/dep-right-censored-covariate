## helper_data_fig1.R
## Table 2 on the Figure 1 design at about 66% right-censoring.
## data_aft() here is data_aft_int(n, c_int = -0.5, p.c = 1.5) from 1000-Main/helper_data.R
## (the Figure 1 generator), without the E[X | D = 0, Y, Z] integrals. The random draws are
## the same and in the same order, so set.seed(s); data_aft(1000) gives the same dataset as
## Figure 1's set.seed(s); data_aft_int(1000, -0.5, 1.5) (seeds 0-2999).
##   log X = 0.1 + 0.1 Z + 0.5 log Exp(1)          gamma = (0.1, 0.1), scale_x = 0.5
##   Y = 1 + (A - X) + Z + e, e ~ N(0, 1), A ~ N(2, 1), Z ~ N(0, 1)
##   log C = -0.5 + 0.1 Y + 0.1 Z + 1.5 log Exp(1)  alpha = (-0.5, 0.1, 0.1), scale 1.5
## About 66% right-censoring (65.9% in the Figure 1 run).
library(survival)

C_INT <- -0.5
P_C   <- 1.5

data_aft <- function(nSubjects = 1000, c_int = C_INT, p.c = P_C) {
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
