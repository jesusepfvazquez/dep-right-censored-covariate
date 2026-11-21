## 1000-CensSpread/main.R
## How the spread of the censoring distribution affects the weighted estimators.
## Design: data_aft_int() (Figure 1 design, helper_data_cs.R), n = 1,000.
## Grid: censoring scale p.c (spread of log C) in {0.75, 1, 1.5, 2, 3} by target censoring
## rate in {50, 65, 75, 85}%: 20 cells, c_int calibrated per cell (calibrate.R).
## p.c = 1.5 is the Figure 1 design. p.c = 0.5 is left out on purpose: log X has scale
## p.x = 0.5, and with p.c <= p.x (and C's scale below X's) 1/pi(x) grows faster than
## f_{X|Z}(x) decays, so E[1/pi(X)] is infinite: the IPW variance is infinite and the
## unguarded Psi integral is undefined (the unguarded fit failed in 5 of 6 test datasets).
## Estimators (theta = (beta, sigma) where it applies), as in Figure 1:
##   oracle; ipw-correct (true pi); ipw-est (SE corrected for alpha-hat);
##   aipw-est-joint (unguarded, f_{X|Z} from the joint MLE); aipw-lambda-est (guarded,
##   working normal from the joint MLE); mle-joint (does not use pi; reference).
## ung_retry = 1: the unguarded fit failed from the IPW start (extreme weights) and was
## refit from the joint MLE start.
## Weight diagnostics per dataset (estimated pi, uncensored rows): largest weight, share
## of the total weight carried by the largest weight, effective sample size / n_obs.
##
## Task t (1-400): cell = ((t - 1) %% 20) + 1, seeds 50 * ((t - 1) %/% 20) + 1:50 (GRID_SET=extra: see below).
## So tasks 1-20 give 50 seeds in every cell; 400 tasks give 1,000 seeds per cell.
## Finished files are skipped. Usage: Rscript main.R <task>

## Byte-code JIT off (kept from the first version; harmless). The joint MLE failures and the
## -20 to -30% MLE bias at p.c = 3 came from the gamma starting value (survreg of W on Z
## ignoring the dependent censoring), fixed in helper_joint_mle.R (Nov 14).
invisible(compiler::enableJIT(0))

source("helper_data_cs.R")
source("helper_table2.R")

C_INT <- c("0.75_0.5" = -0.0231, "1_0.5" = 0.0524, "1.5_0.5" = 0.2155, "2_0.5" = 0.3856, "3_0.5" = 0.7399,
           "0.75_0.65" = -0.4276, "1_0.65" = -0.4361, "1.5_0.65" = -0.4598, "2_0.65" = -0.4849, "3_0.65" = -0.5345,
           "0.75_0.75" = -0.7264, "1_0.75" = -0.7851, "1.5_0.75" = -0.925, "2_0.75" = -1.0777, "3_0.75" = -1.3909,
           "0.75_0.85" = -1.1047, "1_0.85" = -1.2128, "1.5_0.85" = -1.4772, "2_0.85" = -1.7665, "3_0.85" = -2.3714)
GRID <- expand.grid(p.c = c(0.75, 1, 1.5, 2, 3), target = c(0.50, 0.65, 0.75, 0.85))
## Nov 21 extra cells: scale ratios 2.5, 3.5, 5 (p.c = 1.25, 1.75, 2.5), the gap between
## median ESS 67% and 80% in the original grid. c_int from calibrate.R (same seed, 1e6 draws).
## GRID_SET=extra selects these 12 cells: task t (1-240), cell ((t - 1) %% 12) + 1,
## seeds 50 * ((t - 1) %/% 12) + 1:50. Same results/ folder and file names, so the figure
## script picks them up as extra columns. Without GRID_SET the original 20 cells run as before.
if (Sys.getenv("GRID_SET") == "extra") {
  C_INT <- c("1.25_0.5" = 0.1321, "1.75_0.5" = 0.3, "2.5_0.5" = 0.5607,
             "1.25_0.65" = -0.4471, "1.75_0.65" = -0.4728, "2.5_0.65" = -0.5107,
             "1.25_0.75" = -0.8526, "1.75_0.75" = -1.0009, "2.5_0.75" = -1.234,
             "1.25_0.85" = -1.34, "1.75_0.85" = -1.6204, "2.5_0.85" = -2.066)
  GRID <- expand.grid(p.c = c(1.25, 1.75, 2.5), target = c(0.50, 0.65, 0.75, 0.85))
}
SEEDS_PER_TASK <- as.integer(Sys.getenv("SEEDS_PER_TASK", "50"))
n <- 1000

task  <- as.integer(commandArgs(TRUE)[1])
cell  <- ((task - 1) %% nrow(GRID)) + 1
p.c   <- GRID$p.c[cell]; target <- GRID$target[cell]
c_int <- C_INT[[paste0(p.c, "_", target)]]
seeds <- SEEDS_PER_TASK * ((task - 1) %/% nrow(GRID)) + seq_len(SEEDS_PER_TASK)

dir.create("results", showWarnings = FALSE)
safe <- function(expr, what, s) tryCatch(expr, error = function(e) {
  message("FAILED p.c=", p.c, " target=", target, " seed=", s, " ", what, ": ", conditionMessage(e)); NULL })
row <- function(method, fit, gamma = NULL) {
  g <- function(x, k) if (is.null(x) || length(x) < k) NA else as.numeric(x)[k]
  data.frame(Method = method, beta0 = g(fit$beta_est, 1), beta1 = g(fit$beta_est, 2),
             beta2 = g(fit$beta_est, 3), sigma = g(fit$beta_est, 4),
             sebeta0 = g(fit$se_est, 1), sebeta1 = g(fit$se_est, 2), sebeta2 = g(fit$se_est, 3),
             gamma0 = g(gamma, 1), gammaZ = g(gamma, 2), scale_x = g(gamma, 3))
}

for (s in seeds) {
  f <- sprintf("results/est_results_pc%s_t%s_seed_%d.csv", p.c, target, s)
  if (file.exists(f)) next
  set.seed(s)
  dat <- data_cs(n, c_int, p.c)
  t0  <- Sys.time()
  orc   <- safe(estimate_oracle(dat), "oracle", s)
  ipw_c <- safe(estimate_ipw_fixed(dat, myp = dat$myp_ywz_oracle), "ipw-correct", s)
  ipw_e <- safe(estimate_ipw_est(dat), "ipw-est", s)
  mle   <- safe(estimate_beta_likelihood_joint(dat), "mle-joint", s)
  if (!is.null(mle)) {
    th <- mle$theta_full
    v  <- safe(var_mle_fast(dat, c(th[1:3], log(th[4]), th[5:6], log(th[7])), "joint"), "mle-joint var", s)
    mle$se_est <- v$se_est
  }
  gam <- mle$gamma_est
  ## unguarded: IPW start; if that fails (extreme weights), retry from the joint MLE
  ung <- if (!is.null(gam)) tryCatch(estimate_beta_aipw_unguarded_sigma(dat, gam), error = function(e) NULL)
  ung_retry <- is.null(ung) && !is.null(gam)
  if (ung_retry) ung <- safe(estimate_beta_aipw_unguarded_sigma(dat, gam, start = mle$beta_est), "aipw-est-joint", s)
  grd <- if (!is.null(gam)) safe(estimate_beta_aipw_lambda_sigma(dat, wm = "mle", gamma = gam), "aipw-lambda-est", s)

  w   <- 1 / dat$myp_ywz[dat$D == 1]
  out <- rbind(row("oracle", orc), row("ipw-correct", ipw_c), row("ipw-est", ipw_e),
               row("aipw-est-joint", ung, gam), row("aipw-lambda-est", grd, gam),
               row("mle-joint", mle, gam))
  out$p.c <- p.c; out$target <- target; out$c_int <- c_int
  out$pcen <- 1 - mean(dat$D); out$max_w <- max(w); out$max_w_share <- max(w) / sum(w)
  out$ess_frac <- sum(w)^2 / sum(w^2) / length(w)
  out$ung_retry <- ifelse(out$Method == "aipw-est-joint", as.integer(ung_retry), NA)
  out$mle_edge <- ifelse(out$Method == "mle-joint", as.integer(isTRUE(mle$on_edge)), NA)
  out$mle_restarted <- ifelse(out$Method == "mle-joint", as.integer(isTRUE(mle$restarted)), NA)
  out$sim <- s; out$n <- n
  write.csv(out, f, row.names = FALSE)
  message(sprintf("p.c %s target %s seed %d (censoring %.0f%%): %s", p.c, target, s,
                  100 * out$pcen[1], format(round(Sys.time() - t0, 1))))
}
