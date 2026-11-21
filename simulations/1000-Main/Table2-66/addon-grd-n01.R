## addon-grd-n01.R  (run from 1000-Main/Table2/, after or alongside main.R)
## Guarded AIPW estimator with the N(0, 1) working model inside Psi_closed, the same Psi
## as the unguarded rows ung-correct-incorrect-n01 and ung-incorrect-incorrect-n01.
## With the same Psi in both, the guarded and unguarded rows differ only in Lambda, so the
## pair isolates the effect of Lambda under a badly misspecified Psi.
## theta = (beta, sigma); SE from the same sandwich as the other rows (pi known or fixed
## uniform weights, Lambda held at its first-step value).
## Datasets are regenerated from the seeds of main.R (set.seed(s); data_aft(1000)).
## Writes results-addon/est_results_grdn01_seed_<s>.csv (finished seeds are skipped).
## Usage: Rscript addon-grd-n01.R <first seed> <last seed>   (default 0 2999)

library(survival)
source("helper_data_fig1.R")
source("helper_table2.R")

## N(0, 1) working model: estimate_guarded_spec() takes its moments from .working_moments()
.working_moments <- function(dat, wm, gamma) list(m = 0, sd = 1)

args  <- as.integer(commandArgs(TRUE))
seeds <- if (length(args) == 2) args[1]:args[2] else 0:2999
dir.create("results-addon", showWarnings = FALSE)

safe <- function(expr, what, s) tryCatch(expr, error = function(e) {
  message("FAILED seed=", s, " ", what, ": ", conditionMessage(e)); NULL })
row <- function(method, fit) {
  g <- function(x, k) if (is.null(x) || length(x) < k) NA else as.numeric(x)[k]
  data.frame(Method = method, beta0 = g(fit$beta_est, 1), beta1 = g(fit$beta_est, 2),
             beta2 = g(fit$beta_est, 3), sigma = g(fit$beta_est, 4),
             sebeta0 = g(fit$se_est, 1), sebeta1 = g(fit$se_est, 2), sebeta2 = g(fit$se_est, 3))
}

for (s in seeds) {
  f <- paste0("results-addon/est_results_grdn01_seed_", s, ".csv")
  if (file.exists(f)) next
  set.seed(s)
  dat <- data_aft(nSubjects = 1000)
  gc_ <- safe(estimate_guarded_spec(dat), "grd-correct-n01", s)
  gi_ <- safe(estimate_guarded_spec(dat, myp = dat$myp_uniform), "grd-incorrect-n01", s)
  out <- rbind(row("grd-correct-n01", gc_), row("grd-incorrect-n01", gi_))
  out$pcen <- 1 - mean(dat$D); out$sim <- s; out$n <- 1000
  write.csv(out, f, row.names = FALSE)
}
message("done seeds ", min(seeds), "-", max(seeds))
