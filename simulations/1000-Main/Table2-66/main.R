## 1000-Main/Table2-66/main.R
## Table 2 on the Figure 1 design: every estimator under each specification, on the same
## datasets as the 66% point of Figure 1 (set.seed(s); data_aft_int(1000, -0.5, 1.5),
## seeds 0-2999). Data: helper_data_fig1.R. theta = (beta0, beta_{A-X}, beta_Z, sigma).
##
## Method                 Table 2 row
## oracle                 Oracle
## ipw-correct            IPW, correct pi
## ipw-est                IPW, estimated pi (SE corrected for alpha-hat)
## ipw-incorrect          IPW, incorrect pi (1 / Uniform(0.1, 0.9))
## ung-correct-correct    unguarded AIPW, correct pi and correct Psi
## ung-est-est            unguarded AIPW, estimated pi and estimated Psi (f_{X|Z} from the joint MLE)
## ung-incorrect-correct  unguarded AIPW, incorrect pi and correct Psi
## ung-correct-incorrect  unguarded AIPW, correct pi and incorrect Psi (Psi_closed)
## ung-incorrect-incorrect unguarded AIPW, incorrect pi and incorrect Psi
##   Psi_closed uses a working normal with the moments of the true f_{X|Z}: the same Psi as
##   the guarded estimator, so the guarded and unguarded rows differ only in Lambda.
## ung-correct-incorrect-n01, ung-incorrect-incorrect-n01: the same two rows with the
##   N(0, 1) working model of the original Table 2 run (choose in table2.R)
## grd-correct            guarded AIPW, correct pi
## grd-est                guarded AIPW, estimated pi (working normal from the joint MLE)
## grd-incorrect          guarded AIPW, incorrect pi
## mle-correct            MLE, f_{X|Z} known
## mle-twostep            MLE, f_{X|Z} estimated by survreg(Surv(W, D) ~ Z) first (ignores that C depends on Y)
## mle-joint              MLE, f_{X|Z} estimated jointly with theta
## mle-incorrect          MLE, incorrect f_{X|Z} (Weibull scale parameter 1.2 instead of 0.5)
## Inconsistent estimators
## cc                     complete case
## acc                    ACC, correct pi, Psi_closed (theta = (beta, sigma))
## macc                   MACC, correct pi, Psi_closed (theta = (beta, sigma))
## cmi-yz                 CMI, censored X replaced by E[X | Y, Z, X > C] (true model)
## cmi-z                  CMI, censored X replaced by E[X | Z] (true model)
## naive                  W treated as X
##
## One task fits 50 seeds: task k runs seeds 50 (k - 1), ..., 50 k - 1 (60 tasks).
## Finished seeds are skipped, so a task can be resubmitted.
## Usage: Rscript main.R <task>

library(survival)
source("helper_data_fig1.R")
source("helper_table2.R")

task  <- as.integer(commandArgs(TRUE)[1])
SEEDS_PER_TASK <- as.integer(Sys.getenv("SEEDS_PER_TASK", "50"))
seeds <- (SEEDS_PER_TASK * (task - 1)):(SEEDS_PER_TASK * task - 1)   # seeds 0-2999, as in Figure 1
n     <- 1000

dir.create("results", showWarnings = FALSE)

safe <- function(expr, what, s) tryCatch(expr, error = function(e) {
  message("FAILED seed=", s, " ", what, ": ", conditionMessage(e)); NULL })
row <- function(method, fit, gamma = NULL, mle = NULL) {
  g <- function(x, k) if (is.null(x) || length(x) < k) NA else as.numeric(x)[k]
  est <- fit$beta_est; se <- fit$se_est
  data.frame(Method = method, beta0 = g(est, 1), beta1 = g(est, 2), beta2 = g(est, 3),
             sigma = g(est, 4), sebeta0 = g(se, 1), sebeta1 = g(se, 2), sebeta2 = g(se, 3),
             gamma0 = g(gamma, 1), gammaZ = g(gamma, 2), scale_x = g(gamma, 3),
             convcode = if (is.null(mle)) NA else mle$convcode,
             guard = if (is.null(mle$restarted)) NA else as.integer(mle$restarted) + as.integer(mle$on_edge))
}
with_se <- function(fit, v) { if (!is.null(fit)) fit$se_est <- v$se_est; fit }

for (s in seeds) {
  f <- paste0("results/est_results_table2_seed_", s, ".csv")
  if (file.exists(f)) next
  set.seed(s)
  dat <- data_aft(nSubjects = n)
  t0  <- Sys.time()
  unif <- dat$myp_uniform

  orc <- safe(estimate_oracle(dat), "oracle", s)

  ## IPW
  ipw_c <- safe(estimate_ipw_fixed(dat), "ipw-correct", s)
  ipw_e <- safe(estimate_ipw_est(dat), "ipw-est", s)
  ipw_i <- safe(estimate_ipw_fixed(dat, myp = unif), "ipw-incorrect", s)

  ## MLE (the joint fit also supplies f_{X|Z} for the estimated-Psi AIPW rows)
  mle_j  <- safe(estimate_beta_likelihood_joint(dat), "mle-joint", s)
  mle_jv <- if (!is.null(mle_j)) safe(var_mle_fast(dat, c(mle_j$theta_full[1:3], log(mle_j$theta_full[4]),
                                            mle_j$theta_full[5:6], log(mle_j$theta_full[7])), "joint"), "mle-joint var", s)
  mle_c  <- safe(estimate_mle_fixed(dat, gamma_true), "mle-correct", s)
  mle_cv <- if (!is.null(mle_c)) safe(var_mle_fast(dat, mle_c$p_full, "fixed"), "mle-correct var", s)
  g2     <- safe(gamma_twostep(dat), "two-step gamma", s)
  mle_t  <- if (!is.null(g2)) safe(estimate_mle_fixed(dat, g2), "mle-twostep", s)
  mle_tv <- if (!is.null(mle_t)) safe(var_mle_fast(dat, mle_t$p_full, "twostep"), "mle-twostep var", s)
  mle_i  <- safe(estimate_mle_fixed(dat, gamma_wrong), "mle-incorrect", s)
  mle_iv <- if (!is.null(mle_i)) safe(var_mle_fast(dat, mle_i$p_full, "fixed"), "mle-incorrect var", s)
  gam    <- mle_j$gamma_est

  ## unguarded AIPW
  ung_cc <- safe(estimate_unguarded_spec(dat, psi = "eff"), "ung-correct-correct", s)
  ung_ee <- if (!is.null(gam)) safe(estimate_beta_aipw_unguarded_sigma(dat, gam), "ung-est-est", s)   # Figure 1 settings
  ung_ic <- safe(estimate_unguarded_spec(dat, myp = unif, psi = "eff"), "ung-incorrect-correct", s)
  ung_ci <- safe(estimate_unguarded_spec(dat, psi = "closed"), "ung-correct-incorrect", s)
  ung_ii <- safe(estimate_unguarded_spec(dat, myp = unif, psi = "closed"), "ung-incorrect-incorrect", s)
  ung_cn <- safe(estimate_unguarded_spec(dat, psi = "closed", working = "n01"), "ung-correct-incorrect-n01", s)
  ung_in <- safe(estimate_unguarded_spec(dat, myp = unif, psi = "closed", working = "n01"), "ung-incorrect-incorrect-n01", s)

  ## guarded AIPW
  grd_c <- safe(estimate_guarded_spec(dat), "grd-correct", s)
  grd_e <- if (!is.null(gam)) safe(estimate_beta_aipw_lambda_sigma(dat, wm = "mle", gamma = gam), "grd-est", s)
  grd_i <- safe(estimate_guarded_spec(dat, myp = unif), "grd-incorrect", s)

  ## inconsistent estimators
  cc    <- safe(estimate_cc(dat), "cc", s)
  acc   <- safe(estimate_guarded_spec(dat, type = "acc"), "acc", s)
  macc  <- safe(estimate_guarded_spec(dat, type = "macc"), "macc", s)
  cmi1  <- safe(estimate_cmi(dat, "yz"), "cmi-yz", s)
  cmi2  <- safe(estimate_cmi(dat, "z"), "cmi-z", s)
  naive <- safe(estimate_naive(dat), "naive", s)

  out <- rbind(row("oracle", orc),
               row("ipw-correct", ipw_c), row("ipw-est", ipw_e), row("ipw-incorrect", ipw_i),
               row("ung-correct-correct", ung_cc), row("ung-est-est", ung_ee, gam),
               row("ung-incorrect-correct", ung_ic), row("ung-correct-incorrect", ung_ci),
               row("ung-incorrect-incorrect", ung_ii),
               row("ung-correct-incorrect-n01", ung_cn), row("ung-incorrect-incorrect-n01", ung_in),
               row("grd-correct", grd_c), row("grd-est", grd_e, gam), row("grd-incorrect", grd_i),
               row("mle-correct", with_se(mle_c, mle_cv), gamma_true, mle_c),
               row("mle-twostep", with_se(mle_t, mle_tv), g2, mle_t),
               row("mle-joint", with_se(mle_j, mle_jv), gam, mle_j),
               row("mle-incorrect", with_se(mle_i, mle_iv), gamma_wrong, mle_i),
               row("cc", cc), row("acc", acc), row("macc", macc),
               row("cmi-yz", cmi1), row("cmi-z", cmi2),
               row("naive", naive))
  out$pcen <- 1 - mean(dat$D); out$sim <- s; out$n <- n
  write.csv(out, f, row.names = FALSE)
  message("seed ", s, " (censoring ", round(100 * out$pcen[1]), "%): ",
          format(round(Sys.time() - t0, 1)))
}
