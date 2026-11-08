## export-heatmap-cells.R
## Exports the exact cell values behind Figures 2-4 (and the counts figure)
## so the numbers quoted in Section 4.3 can be checked.
## Run from Tables-Revision/:  Rscript export-heatmap-cells.R
## Sources tables-efficiency.R, so the plots are regenerated as well.

source("tables-efficiency.R")

nbin <- function(n) cut(n, breaks = seq(100, 1000, by = 100),
                        include.lowest = TRUE, right = FALSE,
                        labels = c("[100, 200)","[200, 300)","[300, 400)","[400, 500)",
                                   "[500, 600)","[600, 700)","[700, 800)","[800, 900)",
                                   "[900, 1000]"))
keep <- levels(cens_bin10(0.5))[3:9]

## signed bias and empirical SD per cell (same replicates as the bias heatmap)
signed <- est_long %>%
  mutate(cen_bin = cens_bin10(pcen), n_bin = nbin(n)) %>%
  group_by(Method, coef, cen_bin, n_bin) %>%
  summarise(signed_bias = mean(pct_bias, na.rm = TRUE),
            emp_sd = sd(est, na.rm = TRUE), nsim_bias = n(), .groups = "drop") %>%
  filter(nsim_bias >= 100, cen_bin %in% keep)

ipw_sd <- signed %>% filter(Method == "ipw-est") %>%
  select(coef, cen_bin, n_bin, ipw_emp_sd = emp_sd)

signed <- signed %>%
  left_join(ipw_sd, by = c("coef", "cen_bin", "n_bin")) %>%
  mutate(emp_sd_ratio = emp_sd / ipw_emp_sd,
         Method = recode(Method, "ipw-est"="IPW","aipw-est"="AIPW (unguarded)",
                         "aipw-lambda-est"="AIPW (guarded)","mle-est"="MLE"),
         coef = recode(coef, "beta0"="beta[0]",
                       "beta1"="beta[A-X]~(A-X)","beta2"="beta[Z]~(Z)"))

key <- c("Method", "coef", "cen_bin", "n_bin")
fix <- function(d) d %>% mutate(across(all_of(key), as.character))

out <- fix(signed) %>%
  full_join(fix(bias_cells) %>% select(all_of(key), abs_bias = bias), by = key) %>%
  full_join(fix(cov_cells)  %>% select(all_of(key), coverage, nsim_cov = nsim), by = key) %>%
  full_join(fix(eff_cells)  %>% select(all_of(key), se_ratio = eff_ratio), by = key) %>%
  arrange(Method, coef, cen_bin, n_bin)

write.csv(out, "heatmap-cells.csv", row.names = FALSE)
cat("Wrote heatmap-cells.csv with", nrow(out), "rows\n")
