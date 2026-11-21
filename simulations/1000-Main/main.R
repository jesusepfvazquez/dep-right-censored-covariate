## 1000-Main/main.R
## First simulation study (Figure 1): n = 1,000 at eight censoring intercepts.
## Each dataset is set.seed(seed) followed by data_aft_int(n, c_int, 1.5), as in the
## original Figure 1 run, so the datasets match seed for seed.
## All estimators are fit on every dataset (helper_fit_all.R).
## Usage: Rscript main.R <seed>

source("helper_data.R")
source("helper_fit_all.R")

sim_seed <- as.integer(commandArgs(TRUE)[1])
n        <- 1000
myp.c    <- 1.5
c_ints   <- c(4.5, 4, 1.5, 1, 0.5, -0.5, -1, -1.5)

dir.create("results", showWarnings = FALSE)

for (c_int in c_ints) {
  set.seed(sim_seed)
  dat <- data_aft_int(nSubjects = n, c_int, myp.c)
  t0  <- Sys.time()
  out <- fit_all(dat, tag = paste0("seed=", sim_seed, " c_int=", c_int))
  out$pcen  <- 1 - mean(dat$D)
  out$group <- c_int
  out$sim   <- sim_seed
  out$n     <- n
  write.csv(out, paste0("results/est_results_", c_int, "_seed_", sim_seed, ".csv"), row.names = FALSE)
  message("seed ", sim_seed, " c_int ", c_int, " (censoring ", round(100 * (1 - mean(dat$D))), "%): ",
          format(round(Sys.time() - t0, 1)))
}
