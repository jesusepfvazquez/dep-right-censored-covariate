## 1000-Main/Fig1-TwoStep/twostep.R
## Add-on for Figure 1: the two-step MLE (f_{X|Z} from survreg(Surv(W, D) ~ Z) first, which
## ignores that C depends on Y) on the same datasets as Figure 1: set.seed(s), then
## data_aft_int(1000, c_int, 1.5) at the eight censoring intercepts, seeds 0-2999.
## data_aft() in helper_data_fig1.R makes the same draws in the same order (checked against
## ../results in Table2-66/compare-fig1.R), without the E[X | D = 0, Y, Z] integrals.
## The fit and its SE are the Table 2 code (helper_table2.R: gamma_twostep,
## estimate_mle_fixed, var_mle_fast(type = "twostep")), so the c_int = -0.5 rows here
## equal the mle-twostep rows of Table2-66.
##
## One task = 50 seeds x 8 intercepts; finished files are skipped, so a task can be resubmitted.
## Output: results/est_results_<c_int>_seed_<s>.csv, one row (Method = "mle-twostep"),
## the same columns that ../Figures-Tables/fig1-bias-coverage.R reads.
## Usage: Rscript twostep.R <task>

library(survival)
source("helper_data_fig1.R")
source("helper_table2.R")

task  <- as.integer(commandArgs(TRUE)[1])
SEEDS_PER_TASK <- as.integer(Sys.getenv("SEEDS_PER_TASK", "50"))
seeds  <- (SEEDS_PER_TASK * (task - 1)):(SEEDS_PER_TASK * task - 1)   # seeds 0-2999
n      <- 1000
c_ints <- c(4.5, 4, 1.5, 1, 0.5, -0.5, -1, -1.5)

dir.create("results", showWarnings = FALSE)
safe <- function(expr, what, s, c_int) tryCatch(expr, error = function(e) {
  message("FAILED seed=", s, " c_int=", c_int, " ", what, ": ", conditionMessage(e)); NULL })
g <- function(x, k) if (is.null(x) || length(x) < k) NA else as.numeric(x)[k]

for (s in seeds) for (c_int in c_ints) {
  f <- paste0("results/est_results_", c_int, "_seed_", s, ".csv")
  if (file.exists(f)) next
  set.seed(s)
  dat <- data_aft(nSubjects = n, c_int = c_int, p.c = 1.5)
  t0  <- Sys.time()
  g2  <- safe(gamma_twostep(dat), "two-step gamma", s, c_int)
  fit <- if (!is.null(g2)) safe(estimate_mle_fixed(dat, g2), "mle-twostep", s, c_int)
  v   <- if (!is.null(fit)) safe(var_mle_fast(dat, fit$p_full, "twostep"), "mle-twostep var", s, c_int)
  out <- data.frame(Method = "mle-twostep",
                    beta0 = g(fit$beta_est, 1), beta1 = g(fit$beta_est, 2), beta2 = g(fit$beta_est, 3),
                    sigma = g(fit$beta_est, 4),
                    sebeta0 = g(v$se_est, 1), sebeta1 = g(v$se_est, 2), sebeta2 = g(v$se_est, 3),
                    gamma0 = g(g2, 1), gammaZ = g(g2, 2), scale_x = g(g2, 3),
                    convcode = if (is.null(fit)) NA else fit$convcode,
                    pcen = 1 - mean(dat$D), group = c_int, sim = s, n = n)
  write.csv(out, f, row.names = FALSE)
  message("seed ", s, " c_int ", c_int, " (censoring ", round(100 * out$pcen), "%): ",
          format(round(Sys.time() - t0, 1)))
}
