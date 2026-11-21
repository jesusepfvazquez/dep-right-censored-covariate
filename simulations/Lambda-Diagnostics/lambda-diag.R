## Lambda-Diagnostics/lambda-diag.R
## Diagnostics for the guarded AIPW estimator's Lambda and for the weights 1/pi-hat, on the
## same datasets as 1000-Main (Figure 1, Table 2) and 1000-CensSpread (scale-ratio heatmaps).
##
## No estimator is refit. For each dataset we
##   1. regenerate the data with the study's seed (data_cs() makes the same draws in the same
##      order as data_aft_int(), so set.seed(s); data_cs(n, c_int, p.c) is the study's dataset;
##      the realized censoring rate is checked against the study's results file),
##   2. read the f_{X|Z} parameters the guarded fit used (columns gamma0, gammaZ, scale_x of the
##      aipw-lambda-est row; working moments "obs" if they are missing, as in fit_all()),
##   3. rebuild Lambda exactly as estimate_beta_aipw_lambda_sigma() does (helper_guarded_sigma.R,
##      steps 1-2: IPW start, alpha correction, Lambda = -E(U V^T) E(V V^T)^{-1}).
## Recorded per dataset:
##   kappa_M    2-norm condition number of M = sum V V^T, the matrix that is inverted
##   rcond_M    reciprocal condition number of M (singular if < 1e-10 or solve() fails)
##   kappa_Lam  2-norm condition number of Lambda-hat
##   max_w, max_w_share, ess_frac   largest weight, its share of the total, ESS / n uncensored
##   estimates of both AIPW estimators from the results file (for failed / diverged counts)
##
## Usage (from this folder; the study folders sit next to it):
##   Rscript lambda-diag.R main <task 1-8>    one censoring intercept, seeds 0-2999
##   Rscript lambda-diag.R cs   <task 1-20>   one (p.c, censoring) cell, seeds 1-1000
##   Rscript lambda-diag.R cs-extra <task 1-12> one extra cell (scale ratios 2.5, 3.5, 5)
## Output: results/diag_main_<c_int>.csv or results/diag_cs_pc<p.c>_t<target>.csv
## Set MAIN_DIR / CS_DIR if the study folders are elsewhere.

suppressMessages(library(survival))
MAIN_DIR <- Sys.getenv("MAIN_DIR", "../1000-Main")
CS_DIR   <- Sys.getenv("CS_DIR",   "../1000-CensSpread")
source(file.path(CS_DIR, "helper_data_cs.R"))
source(file.path(CS_DIR, "helper_guarded_sigma.R"))   # identical in 1000-Main

args  <- commandArgs(TRUE)
study <- args[1]; task <- as.integer(args[2])
n <- 1000

if (study == "main") {
  c_ints <- c(4.5, 4, 1.5, 1, 0.5, -0.5, -1, -1.5)
  c_int  <- c_ints[task]; p.c <- 1.5; target <- NA
  seeds  <- 0:2999
  res_file <- function(s) file.path(MAIN_DIR, "results", paste0("est_results_", c_int, "_seed_", s, ".csv"))
  out_file <- file.path("results", paste0("diag_main_", c_int, ".csv"))
} else if (study %in% c("cs", "cs-extra")) {
  C_INT <- c("0.75_0.5" = -0.0231, "1_0.5" = 0.0524, "1.5_0.5" = 0.2155, "2_0.5" = 0.3856, "3_0.5" = 0.7399,
             "0.75_0.65" = -0.4276, "1_0.65" = -0.4361, "1.5_0.65" = -0.4598, "2_0.65" = -0.4849, "3_0.65" = -0.5345,
             "0.75_0.75" = -0.7264, "1_0.75" = -0.7851, "1.5_0.75" = -0.925, "2_0.75" = -1.0777, "3_0.75" = -1.3909,
             "0.75_0.85" = -1.1047, "1_0.85" = -1.2128, "1.5_0.85" = -1.4772, "2_0.85" = -1.7665, "3_0.85" = -2.3714)
  GRID  <- expand.grid(p.c = c(0.75, 1, 1.5, 2, 3), target = c(0.50, 0.65, 0.75, 0.85))
  if (study == "cs-extra") {   # Nov 21: scale ratios 2.5, 3.5, 5 (tasks 1-12); same c_int as 1000-CensSpread/main.R
    C_INT <- c("1.25_0.5" = 0.1321, "1.75_0.5" = 0.3, "2.5_0.5" = 0.5607,
               "1.25_0.65" = -0.4471, "1.75_0.65" = -0.4728, "2.5_0.65" = -0.5107,
               "1.25_0.75" = -0.8526, "1.75_0.75" = -1.0009, "2.5_0.75" = -1.234,
               "1.25_0.85" = -1.34, "1.75_0.85" = -1.6204, "2.5_0.85" = -2.066)
    GRID  <- expand.grid(p.c = c(1.25, 1.75, 2.5), target = c(0.50, 0.65, 0.75, 0.85))
    study <- "cs"              # rows are summarized with the original cells
  }
  p.c   <- GRID$p.c[task]; target <- GRID$target[task]
  c_int <- C_INT[[paste0(p.c, "_", target)]]
  seeds <- 1:1000
  res_file <- function(s) file.path(CS_DIR, "results", sprintf("est_results_pc%s_t%s_seed_%d.csv", p.c, target, s))
  out_file <- file.path("results", sprintf("diag_cs_pc%s_t%s.csv", p.c, target))
} else stop("study must be 'main', 'cs' or 'cs-extra'")

if (!is.null(seeds_env <- Sys.getenv("SEEDS", unset = NA)) && !is.na(seeds_env))   # e.g. SEEDS=0:9 for a test
  seeds <- eval(parse(text = seeds_env))

diag_one <- function(dat, g) {
  wr    <- survreg(Surv(W, 1 - D) ~ y + Z, data = dat, dist = "w")
  alpha <- c(coef(wr), wr$scale)
  wm    <- if (all(is.finite(g))) "mle" else "obs"
  mm    <- .working_moments(dat, wm, if (wm == "mle") g else NULL)
  pi_   <- .pi_fun(alpha, dat$W, dat$y, dat$Z)
  wts   <- dat$D / pi_
  fit   <- lm(y ~ AW + Z, data = dat, weights = wts)
  th0   <- c(coef(fit), sqrt(weighted.mean(resid(fit)^2, wts)))
  p0    <- .pieces(th0, dat, alpha, mm$m, mm$sd, TRUE)
  U     <- p0$S * wts + p0$corr %*% t(p0$A_ipw)
  V     <- p0$Ps * (1 - wts) + p0$corr %*% t(p0$A_at)
  M     <- crossprod(V, V)
  Mi    <- tryCatch(solve(M), error = function(e) NULL)
  Lam   <- if (is.null(Mi)) NULL else -crossprod(U, V) %*% Mi
  w     <- 1 / pi_[dat$D == 1]
  data.frame(wm = wm,
             kappa_M   = kappa(M, exact = TRUE), rcond_M = rcond(M), solve_failed = as.integer(is.null(Mi)),
             kappa_Lam = if (is.null(Lam)) NA else kappa(Lam, exact = TRUE),
             max_w = max(w), max_w_share = max(w) / sum(w), ess_frac = sum(w)^2 / sum(w^2) / length(w))
}

pick <- function(r, method, cols) {
  x <- r[r$Method == method, cols, drop = FALSE]
  if (nrow(x) == 0) setNames(as.list(rep(NA_real_, length(cols))), cols) else as.list(x[1, ])
}

dir.create("results", showWarnings = FALSE)
t0  <- Sys.time()
out <- do.call(rbind, lapply(seeds, function(s) {
  f <- res_file(s)
  r <- if (file.exists(f)) read.csv(f) else NULL
  g <- if (is.null(r)) rep(NA_real_, 3) else unlist(pick(r, "aipw-lambda-est", c("gamma0", "gammaZ", "scale_x")))
  set.seed(s)
  dat <- data_cs(n, c_int, p.c)
  d <- tryCatch(diag_one(dat, g), error = function(e) {
    message("seed ", s, ": ", conditionMessage(e)); NULL })
  if (is.null(d)) d <- data.frame(wm = NA, kappa_M = NA, rcond_M = NA, solve_failed = NA, kappa_Lam = NA,
                                  max_w = NA, max_w_share = NA, ess_frac = NA)
  grd <- if (is.null(r)) list(beta0 = NA, beta1 = NA, beta2 = NA) else pick(r, "aipw-lambda-est", c("beta0", "beta1", "beta2"))
  ung <- if (is.null(r)) list(beta0 = NA, beta1 = NA, beta2 = NA) else pick(r, "aipw-est-joint",  c("beta0", "beta1", "beta2"))
  cbind(data.frame(study = study, sim = s, c_int = c_int, p.c = p.c, target = target,
                   pcen = 1 - mean(dat$D),
                   pcen_results = if (is.null(r) || !"pcen" %in% names(r)) NA else r$pcen[1],
                   has_results = as.integer(!is.null(r))),
        d,
        grd_b0 = grd$beta0, grd_b1 = grd$beta1, grd_b2 = grd$beta2,
        ung_b0 = ung$beta0, ung_b1 = ung$beta1, ung_b2 = ung$beta2)
}))
write.csv(out, out_file, row.names = FALSE)
mism <- sum(abs(out$pcen - out$pcen_results) > 1e-8, na.rm = TRUE)
message(sprintf("%s task %d: %d datasets, %d without results file, %d censoring-rate mismatches, %s",
                study, task, nrow(out), sum(out$has_results == 0), mism, format(round(Sys.time() - t0, 1))))
if (mism > 0) warning("realized censoring rate differs from the results file: datasets do not match")
