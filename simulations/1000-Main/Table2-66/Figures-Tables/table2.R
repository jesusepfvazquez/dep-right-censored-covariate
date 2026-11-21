## table2.R: Table 2 (tab:chp3-sims_1000) from ../results.
## Run from 1000-Main/Table2-66/Figures-Tables/. Values follow Table 2: percent bias, mean
## SE x 100, empirical SD x 100 and 95% coverage, for beta0, beta_{A-X} and beta_Z.
## Writes
##   tables/table2-summary.csv   one line per estimator and coefficient (with failures)
##   tables/table2-rows.tex      the body of Table 2 (Oracle, Consistent and Inconsistent
##                               Estimators), ready to paste between the coefficient
##                               header row and \bottomrule
library(tidyverse)

## incorrect Psi for the unguarded AIPW estimator:
##   "gamma" = Psi_closed with the moments of the true f_{X|Z} (same Psi as the guarded rows)
##   "n01"   = Psi_closed with a N(0, 1) working model (as in the original Table 2)
PSI_WRONG <- "gamma"
## E[X | Z] for the CMI row: scale Gamma(1 + 0.5), the true mean of X | Z (no ODRCC variant here)
CMI_Z <- "correct"

res <- list.files("../results", "^est_results_table2_seed_[0-9]+\\.csv$", full.names = TRUE) %>%
  map_df(read.csv) %>% as_tibble()
## add-on (addon-grd-n01.R): guarded AIPW with the N(0, 1) working model, if it was run
add <- list.files("../results-addon", "^est_results_grdn01_seed_[0-9]+\\.csv$", full.names = TRUE)
if (length(add)) res <- bind_rows(res, map_df(add, read.csv) %>% filter(sim %in% res$sim))
cat("Replicates:", n_distinct(res$sim), "; mean censoring", round(100 * mean(res$pcen), 1), "%\n")

## fits that failed or have no SE are left out of that estimator's summaries
fails <- res %>% group_by(Method) %>%
  summarise(failed = sum(is.na(beta0)), se_NA = sum(!is.na(beta0) & is.na(sebeta0)),
            mle_guard = sum(guard > 0, na.rm = TRUE), mle_nonconv = sum(convcode != 0, na.rm = TRUE),
            .groups = "drop")
print(fails, n = Inf)

z <- qnorm(0.975)
summ <- res %>%
  pivot_longer(c(beta0, beta1, beta2, sebeta0, sebeta1, sebeta2),
               names_to = c(".value", "coef"), names_pattern = "(beta|sebeta)(\\d)") %>%
  rename(est = beta, se = sebeta) %>%
  filter(is.finite(est)) %>%
  group_by(Method, coef) %>%
  summarise(reps = n(), bias = 100 * mean(est - 1),
            cov = 100 * mean(abs(est - 1) < z * se, na.rm = TRUE),   # per-replicate SE (before se is overwritten)
            se = 100 * mean(se, na.rm = TRUE), sd = 100 * sd(est),
            .groups = "drop")

dir.create("tables", showWarnings = FALSE)
write.csv(left_join(summ, fails, by = "Method"), "tables/table2-summary.csv", row.names = FALSE)
print(summ %>% mutate(across(where(is.double), ~ round(.x, 2))), n = Inf)

if (PSI_WRONG == "n01") {
  summ <- summ %>% filter(!Method %in% c("ung-correct-incorrect", "ung-incorrect-incorrect")) %>%
    mutate(Method = sub("-n01$", "", Method))
}
if (CMI_Z == "odrcc") summ <- summ %>% filter(Method != "cmi-z") %>% mutate(Method = sub("^cmi-z-odrcc$", "cmi-z", Method))
cells <- summ %>% arrange(Method, coef) %>% group_by(Method) %>%
  summarise(cells = paste(sprintf("%.2f & %.2f & %.2f & %.2f", bias, se, sd, cov), collapse = " & "),
            .groups = "drop") %>%
  { setNames(.$cells, .$Method) }
r <- function(m) if (is.na(cells[m])) "NA & NA & NA & NA & NA & NA & NA & NA & NA & NA & NA & NA" else cells[[m]]
mr <- function(lab) paste0("\\multirow[t]{-2}{*}{\\raggedright\\arraybackslash \\hspace{1em}", lab, "}")
pi_ <- function(w) paste0(w, " $\\pi_{Y,X,\\bZ}$")
ps_ <- function(w) paste0(w, " $\\bPsi_{\\rm AIPW}(y,\\bz)$")

tex <- c(
  paste0("\\textbf{Oracle} & & ", r("oracle"), " \\\\"),
  "\\addlinespace",
  "\\multicolumn{14}{l}{\\textbf{Consistent Estimators}} \\\\",
  "\\addlinespace",
  "%%% IPW",
  paste0(" & ", pi_("correct"), " & ", r("ipw-correct"), " \\\\"),
  mr("IPW"),
  paste0("& ", pi_("estimated"), " & ", r("ipw-est"), " \\\\"),
  paste0("& ", pi_("incorrect"), " & ", r("ipw-incorrect"), " \\\\"),
  "\\addlinespace",
  "%%% AIPW",
  "&\\multicolumn{5}{l}{\\underline{unguarded}}\\\\",
  mr("AIPW"),
  paste0(" & ", pi_("correct"),   " and ", ps_("correct"),   " & ", r("ung-correct-correct"), " \\\\"),
  paste0(" & ", pi_("estimated"), " and ", ps_("estimated"), " & ", r("ung-est-est"), " \\\\"),
  paste0(" & ", pi_("incorrect"), " and ", ps_("correct"),   " & ", r("ung-incorrect-correct"), " \\\\"),
  paste0(" & ", pi_("correct"),   " and ", ps_("incorrect"), " & ", r("ung-correct-incorrect"), " \\\\"),
  paste0(" & ", pi_("incorrect"), " and ", ps_("incorrect"), " & ", r("ung-incorrect-incorrect"), " \\\\"),
  "&\\multicolumn{5}{l}{\\underline{guarded}}\\\\",
  mr("AIPW"),
  paste0(" & ", pi_("correct"),   " & ", r("grd-correct"), " \\\\"),
  paste0(" & ", pi_("estimated"), " & ", r("grd-est"), " \\\\"),
  paste0(" & ", pi_("incorrect"), " & ", r("grd-incorrect"), " \\\\"),
  ## guarded rows with the N(0, 1) working model (add-on); same Psi as the unguarded
  ## "-n01" rows, so the two differ only in Lambda. Uncommented when PSI_WRONG = "n01".
  paste0(if (PSI_WRONG == "n01") "" else "% ", " & ", pi_("correct"),   " and ", ps_("incorrect"), " & ", r("grd-correct-n01"), " \\\\"),
  paste0(if (PSI_WRONG == "n01") "" else "% ", " & ", pi_("incorrect"), " and ", ps_("incorrect"), " & ", r("grd-incorrect-n01"), " \\\\"),
  "\\addlinespace",
  "%%% MLE",
  paste0(" & correct $f_{X|\\bZ}$ & ", r("mle-correct"), " \\\\"),
  mr("MLE"),
  paste0("& estimated $f_{X|\\bZ}$ (two-step) & ", r("mle-twostep"), " \\\\"),
  paste0("& estimated $f_{X|\\bZ}$ (joint) & ", r("mle-joint"), " \\\\"),
  paste0("& incorrect $f_{X|\\bZ}$ & ", r("mle-incorrect"), " \\\\"),
  "\\addlinespace",
  "\\multicolumn{14}{l}{\\textbf{Inconsistent Estimators}} \\\\",
  "\\addlinespace",
  paste0("\\hspace{1em} CC & & ", r("cc"), " \\\\"),
  "\\addlinespace",
  paste0("\\hspace{1em} ACC & correct $\\pi_{Y,X,\\bZ}$ and correct $\\bPsi_{\\rm ACC}(y,\\bz)$ & ", r("acc"), " \\\\"),
  "\\addlinespace",
  paste0("\\hspace{1em} MACC & correct $\\pi_{Y,X,\\bZ}$ and correct $\\bPsi_{\\rm MACC}(y,\\bz)$ & ", r("macc"), " \\\\"),
  "\\addlinespace",
  "\\multicolumn{14}{l}{\\hspace{1em} Conditional mean imputation} \\\\",
  paste0(" & correct $E_{X|Y,\\bZ, X>C}(X)$ & ", r("cmi-yz"), " \\\\"),
  paste0("& correct $E_{X|\\bZ}(X)$ & ", r("cmi-z"), " \\\\"),
  paste0("\\hspace{1em} Naive & & ", r("naive"), " \\\\"),
  "\\addlinespace")
writeLines(tex, "tables/table2-rows.tex")
cat(tex, sep = "\n")
