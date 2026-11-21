## Lambda-Diagnostics/lambda-diag-summary.R
## Summarizes results/diag_*.csv into the two supplement tables.
##   tables/lambda-diag-main.csv / .tex   one row per censoring rate (Figure 1 design)
##   tables/lambda-diag-cs.csv   / .tex   one row per (censoring rate, scale ratio) cell
## Failed: no finite estimate. Diverged: a finite estimate more than 10 from the truth (1) on
## any coefficient, the rule in 1000-CensSpread/Figures-Tables/cens-spread.R.
## Singular: rcond(M) < 1e-10 or solve(M) failed. Usage: Rscript lambda-diag-summary.R

suppressMessages(library(dplyr))
dir.create("tables", showWarnings = FALSE)
d <- list.files("results", "^diag_.*\\.csv$", full.names = TRUE) |> lapply(read.csv) |> bind_rows()

fd <- function(b0, b1, b2) {
  fin <- is.finite(b0) & is.finite(b1) & is.finite(b2)
  list(failed = sum(!fin), diverged = sum(fin & (abs(b0 - 1) > 10 | abs(b1 - 1) > 10 | abs(b2 - 1) > 10)))
}
summ <- function(x) x |> summarise(
  mismatch    = sum(abs(pcen - pcen_results) > 1e-8, na.rm = TRUE),
  no_results  = sum(has_results == 0),
  reps        = n(),
  pcen        = 100 * mean(pcen),
  ess         = 100 * median(ess_frac, na.rm = TRUE),
  kappa_M_med = median(kappa_M, na.rm = TRUE),
  kappa_M_p95 = quantile(kappa_M, 0.95, na.rm = TRUE),
  kappa_L_med = median(kappa_Lam, na.rm = TRUE),
  singular    = 100 * mean(rcond_M < 1e-10 | solve_failed == 1, na.rm = TRUE),
  max_w_med   = median(max_w, na.rm = TRUE),
  ## failed / diverged counted only where the study's results file exists
  grd_failed   = with(fd(grd_b0[has_results == 1], grd_b1[has_results == 1], grd_b2[has_results == 1]), failed),
  grd_diverged = with(fd(grd_b0[has_results == 1], grd_b1[has_results == 1], grd_b2[has_results == 1]), diverged),
  ung_failed   = with(fd(ung_b0[has_results == 1], ung_b1[has_results == 1], ung_b2[has_results == 1]), failed),
  ung_diverged = with(fd(ung_b0[has_results == 1], ung_b1[has_results == 1], ung_b2[has_results == 1]), diverged),
  .groups = "drop")

f1 <- function(x) formatC(x, format = "f", digits = 1)
f0 <- function(x) formatC(x, format = "f", digits = 0)

## ---- 1000-Main: one row per censoring rate -------------------------------------------------
m <- d |> filter(study == "main") |> group_by(c_int) |> summ() |> arrange(pcen)
if (nrow(m)) {
  write.csv(m, "tables/lambda-diag-main.csv", row.names = FALSE)
  print(as.data.frame(m), digits = 3)
  body <- sprintf("%s\\%% & %s & %s & %s & %s & %s & %s & %d \\\\",
                  f0(m$pcen), formatC(m$reps, big.mark = ",", format = "d"), f1(m$ess),
                  f1(m$kappa_M_med), f1(m$kappa_M_p95), f1(m$singular), f1(m$max_w_med),
                  m$grd_failed + m$grd_diverged)
  writeLines(c("\\begin{tabular}{rrrrrrrr}", "\\toprule",
               " & & & \\multicolumn{2}{c}{Condition number} & Singular & Max weight & Guarded \\\\",
               "\\cmidrule(lr){4-5}",
               "Right-censoring & Reps & ESS (\\%) & Median & 95th pct & fits (\\%) & Median & failed \\\\",
               "\\midrule", body, "\\bottomrule", "\\end{tabular}"), "tables/lambda-diag-main.tex")
}

## ---- 1000-CensSpread: one row per cell ---------------------------------------------------------
cs <- d |> filter(study == "cs") |> group_by(target, p.c) |> summ() |>
  mutate(ratio = p.c / 0.5) |> arrange(target, ratio)
if (nrow(cs)) {
  write.csv(cs, "tables/lambda-diag-cs.csv", row.names = FALSE)
  print(as.data.frame(cs), digits = 3)
  body <- sprintf("%s & %s & %s & %s & %s & %s & %s & %d / %d & %d / %d \\\\",
                  ifelse(!duplicated(cs$target), paste0(f0(100 * cs$target), "\\%"), ""),
                  formatC(cs$ratio, format = "g"), f1(cs$ess), f1(cs$kappa_M_med), f1(cs$kappa_M_p95),
                  f1(cs$singular), f1(cs$max_w_med),
                  cs$ung_failed, cs$ung_diverged, cs$grd_failed, cs$grd_diverged)
  rule <- which(!duplicated(cs$target))[-1]
  body[rule] <- paste0("\\addlinespace\n", body[rule])
  writeLines(c("\\begin{tabular}{rrrrrrrcc}", "\\toprule",
               " & & & \\multicolumn{2}{c}{Condition number} & Singular & Max weight & \\multicolumn{2}{c}{Failed / diverged} \\\\",
               "\\cmidrule(lr){4-5} \\cmidrule(lr){8-9}",
               "Right-censoring & Scale ratio & ESS (\\%) & Median & 95th pct & fits (\\%) & Median & Unguarded & Guarded \\\\",
               "\\midrule", body, "\\bottomrule", "\\end{tabular}"), "tables/lambda-diag-cs.tex")
}
message("Check: no_results and mismatch should be 0 in every row.")
