## 1000-Sensitivity/main.R
## Sensitivity analysis for the assumption C independent of X given (Y, Z).
## Design (Section 5.3): n = 800, about 78% right-censoring (PREDICT-HD), and a
## term xi * (X - mean(X)) added to the log scale of the Weibull model for C.
## The censoring intercept is recalibrated for each xi so the censoring rate stays
## near 78%. Each dataset is set.seed(seed) followed by data_aft_int().
## All estimators are fit (helper_fit_all.R). The figure reports bias only, so the
## slow IPW and MLE sandwiches are off (do_var = FALSE); the AIPW SEs are always computed.
## Usage: Rscript main.R <seed>

source("helper_data_sen.R")
source("helper_fit_all.R")

sim_seed    <- as.integer(commandArgs(TRUE)[1])
n_fixed     <- 800
myp.c       <- 1.5
target_pcen <- 0.78
xi_grid     <- seq(0, 1, by = 0.2)   # range reported in the paper
do_var      <- FALSE                 # TRUE for IPW and MLE sandwich SEs (slow)

dir.create("results", showWarnings = FALSE)

## calibrated censoring intercept per xi (seed-independent)
c_int_by_xi <- sapply(xi_grid, function(xx)
  calibrate_c_int(xx, target_pcen = target_pcen, p.c = myp.c))

for (k in seq_along(xi_grid)) {
  xi <- xi_grid[k]; c_int <- c_int_by_xi[k]
  set.seed(sim_seed)
  dat <- data_aft_int(nSubjects = n_fixed, c_int, myp.c, xi = xi)
  t0  <- Sys.time()
  out <- fit_all(dat, do_var_ipw = do_var, do_var_mle = do_var, tag = paste0("seed=", sim_seed, " xi=", xi))
  out$pcen  <- 1 - mean(dat$D)
  out$group <- c_int
  out$sim   <- sim_seed
  out$n     <- n_fixed
  out$xi    <- xi
  write.csv(out, paste0("results/est_results_xi", xi, "_seed_", sim_seed, ".csv"), row.names = FALSE)
  message("seed ", sim_seed, " xi ", xi, ": ", format(round(Sys.time() - t0, 1)))
}
