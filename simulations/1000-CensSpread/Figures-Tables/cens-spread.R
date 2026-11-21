## cens-spread.R: figures and tables for 1000-CensSpread. Run from 1000-CensSpread/Figures-Tables/.
## Writes
##   tables/cens-spread-summary.csv     bias, SE, SD, coverage, MSE per method, cell and coefficient
##   tables/cens-spread-weights.csv     weight diagnostics per cell
##   figures/cs-lines-<target>.png      bias and coverage against p.c, one figure per censoring rate
##                                      (same layout as Figure 1, x axis = censoring scale)
##   figures/cs-heat-coverage.png       coverage heat map, p.c by censoring rate, per method and coefficient
##   figures/cs-heat-bias.png           percent bias heat map
##   figures/cs-heat-re-mse.png         relative efficiency, MSE(IPW) / MSE(estimator)
##   figures/cs-heat-re-se.png          relative efficiency, {mean SE(IPW) / mean SE(estimator)}^2
##   tables/cens-spread-re.csv          both relative efficiencies
##   figures/cs-heat-weights.png        effective sample size / n_obs and largest weight share
library(tidyverse)
z <- qnorm(0.975)

res <- list.files("../results", "^est_results_pc.*\\.csv$", full.names = TRUE) %>% map_df(read.csv) %>% as_tibble()

method_levels <- c("oracle", "ipw-correct", "ipw-est", "aipw-est-joint", "aipw-lambda-est", "mle-joint")
method_labels <- c("Oracle", "IPW (true pi)", "IPW", "AIPW (unguarded)", "AIPW (guarded)", "MLE (joint)")
method_colors <- c("Oracle" = "#CC79A7", "IPW (true pi)" = "#F0E442", "IPW" = "#E69F00",
                   "AIPW (unguarded)" = "#56B4E9", "AIPW (guarded)" = "#009E73", "MLE (joint)" = "#000000")
param_labels <- c(beta0 = "beta[0]", beta1 = "beta[A-X]~(A-X)", beta2 = "beta[Z]~(Z)")

## diverged: a finite estimate more than 10 (1,000% of the true value 1) from the truth.
## The guarded estimator gives such values (up to 1e14) when a few weights are extreme
## (p.c <= 1, 65% censoring or more). Those fits are counted per cell (column diverged)
## and left out of bias, SD and coverage, the same way failed fits are.
res <- res %>% mutate(diverged = if_all(c(beta0, beta1, beta2), is.finite) &
                        (abs(beta0 - 1) > 10 | abs(beta1 - 1) > 10 | abs(beta2 - 1) > 10))
P_X <- 0.5   # scale of log X; the x-axes show the ratio p.c / P_X (1.5 to 6)
XLAB <- "Ratio of the scale of log C to the scale of log X"

divg <- res %>% group_by(Method, p.c, target) %>%
  summarise(failed = sum(!is.finite(beta0)), diverged = sum(diverged), .groups = "drop")

long <- res %>%
  filter(if_all(c(beta0, beta1, beta2), is.finite), !diverged) %>%
  rename(est0 = beta0, est1 = beta1, est2 = beta2, se0 = sebeta0, se1 = sebeta1, se2 = sebeta2) %>%
  pivot_longer(c(est0, est1, est2, se0, se1, se2), names_to = c(".value", "Parameter"),
               names_pattern = "(est|se)(\\d)") %>%
  mutate(Parameter = paste0("beta", Parameter))

tab <- long %>% group_by(Method, p.c, target, Parameter) %>%
  summarise(reps = n(), pcen = 100 * mean(pcen), bias = 100 * mean(est - 1),
            se_mean = 100 * mean(se, na.rm = TRUE), sd = 100 * sd(est),
            coverage = 100 * mean(abs(est - 1) < z * se, na.rm = TRUE),   # per-replicate SE
            mse = 1e4 * mean((est - 1)^2), .groups = "drop") %>%
  rename(se = se_mean) %>%
  left_join(divg, by = c("Method", "p.c", "target")) %>%
  mutate(Method = factor(Method, levels = method_levels, labels = method_labels),
         Parameter = factor(Parameter, levels = names(param_labels)))

wts <- res %>% distinct(p.c, target, sim, pcen, max_w, max_w_share, ess_frac) %>%
  group_by(p.c, target) %>%
  summarise(seeds = n(), pcen = 100 * mean(pcen), max_w_median = median(max_w),
            max_w_share_median = median(max_w_share), ess_frac_median = median(ess_frac), .groups = "drop")

dir.create("tables", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
write.csv(tab, "tables/cens-spread-summary.csv", row.names = FALSE)
write.csv(wts, "tables/cens-spread-weights.csv", row.names = FALSE)
print(wts, n = Inf)
tab %>% filter(Parameter == "beta0") %>% select(Method, p.c, target, reps, failed, diverged, bias, sd, se, coverage) %>%
  mutate(across(where(is.double), ~ round(.x, 2))) %>% arrange(target, Method, p.c) %>% print(n = Inf)

## 1. line plots against p.c, one per censoring rate (Figure 1 layout)
show <- c("Oracle", "IPW", "AIPW (unguarded)", "AIPW (guarded)", "MLE (joint)")
for (tg in sort(unique(tab$target))) {
  d <- tab %>% filter(target == tg, Method %in% show) %>%
    select(Method, p.c, Parameter, bias, coverage) %>%
    pivot_longer(c(bias, coverage), names_to = "panel") %>%
    mutate(panel = factor(panel, levels = c("bias", "coverage"), labels = c("Percent bias (%)", "Coverage (%)")))
  ## x-axis labels: scale of log C with the median effective sample size underneath
  ## (ESS depends on the setting only, not on the estimator or the coefficient)
  xlab <- wts %>% filter(target == tg) %>% arrange(p.c) %>%
    mutate(lab = sprintf("%s\n(%.0f%%)", p.c / P_X, 100 * ess_frac_median))
  href <- tibble(panel = factor(c("Percent bias (%)", "Coverage (%)"), levels = c("Percent bias (%)", "Coverage (%)")),
                 yint = c(0, 95))
  p <- ggplot(d, aes(x = p.c, y = value, color = Method, group = Method)) +
    geom_hline(data = href, aes(yintercept = yint), linetype = "dotted", color = "red") +
    geom_point() + geom_line() +
    facet_grid(panel ~ Parameter, scales = "free_y", switch = "y",
               labeller = labeller(Parameter = as_labeller(param_labels, label_parsed), panel = label_value)) +
    scale_color_manual(values = method_colors, name = "Estimator") +
    scale_x_continuous(trans = "log", breaks = xlab$p.c, labels = xlab$lab) +
    labs(x = paste0(XLAB, " (median ESS as % of uncensored n)"), y = NULL,
         title = paste0(100 * tg, "% right-censoring")) +
    theme_minimal() +
    theme(legend.position = "bottom", panel.spacing = unit(0.6, "lines"),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 90))
  ggsave(sprintf("figures/cs-lines-%s.png", 100 * tg), p, height = 6, width = 8, units = "in", bg = "white")
}

## 2. heat maps: p.c by censoring rate, per method (columns) and coefficient (rows)
heat <- function(var, lab, fill_scale, file) {
  d <- tab %>% filter(Method %in% show, Method != "Oracle")
  p <- ggplot(d, aes(x = factor(p.c / P_X), y = factor(100 * target), fill = .data[[var]])) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f", .data[[var]])), size = 2.6) +
    facet_grid(Parameter ~ Method, labeller = labeller(Parameter = as_labeller(param_labels, label_parsed))) +
    fill_scale + labs(x = XLAB, y = "Percent censored (%)", fill = lab) +
    theme_minimal() + theme(panel.grid = element_blank(), legend.position = "bottom")
  ggsave(file, p, height = 6.5, width = 10, units = "in", bg = "white")
}
heat("coverage", "Coverage (%)", scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
     midpoint = 95, limits = c(70, 100), oob = scales::squish), "figures/cs-heat-coverage.png")
heat("bias", "Percent bias (%)", scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
     midpoint = 0, limits = c(-5, 5), oob = scales::squish), "figures/cs-heat-bias.png")

## relative efficiency against the IPW estimator (estimated pi), per cell and coefficient,
## two versions kept separate:
##   RE (MSE) = MSE(IPW) / MSE(estimator)              empirical, bias^2 + variance
##   RE (SE)  = {mean SE(IPW) / mean SE(estimator)}^2  model-based, squared so that the
##                                                     ratio is on the variance scale
## RE > 1: the estimator is more efficient than the IPW estimator (blue); RE < 1: less (red).
## Failed and diverged fits are left out (as in the summaries).
re <- tab %>% filter(Method == "IPW") %>% select(p.c, target, Parameter, mse_ipw = mse, se_ipw = se) %>%
  right_join(tab %>% filter(Method %in% c("AIPW (unguarded)", "AIPW (guarded)", "MLE (joint)")),
             by = c("p.c", "target", "Parameter")) %>%
  mutate(re_mse = mse_ipw / mse, re_se = (se_ipw / se)^2)
write.csv(re %>% select(Method, p.c, target, Parameter, mse, mse_ipw, re_mse, se, se_ipw, re_se),
          "tables/cens-spread-re.csv", row.names = FALSE)
## fill on the log2 scale (1 = white) so that 2x and 0.5x get equal color weight
heat_re <- function(var, lab, file) {
  p <- ggplot(re, aes(x = factor(p.c / P_X), y = factor(100 * target), fill = log2(.data[[var]]))) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", .data[[var]])), size = 2.6) +
    facet_grid(Parameter ~ Method, labeller = labeller(Parameter = as_labeller(param_labels, label_parsed))) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
                         limits = c(-2, 3), breaks = -2:3, labels = c(0.25, 0.5, 1, 2, 4, 8),
                         oob = scales::squish) +
    labs(x = XLAB, y = "Percent censored (%)", fill = lab) +
    guides(fill = guide_colorbar(barwidth = unit(14, "lines"), title.vjust = 0.8)) +
    theme_minimal() + theme(panel.grid = element_blank(), legend.position = "bottom")
  ggsave(file, p, height = 6.5, width = 8, units = "in", bg = "white")
}
heat_re("re_mse", "MSE ratio, MSE(IPW) / MSE", "figures/cs-heat-re-mse.png")
heat_re("re_se", "Relative efficiency, {mean SE(IPW) / mean SE}^2", "figures/cs-heat-re-se.png")

## median ESS of the weights 1/pi-hat, as a percentage of the uncensored n (same measure as Figure 1)
pw <- wts %>% select(p.c, target, ess_frac_median)
p <- ggplot(pw, aes(factor(p.c / P_X), factor(100 * target), fill = 100 * ess_frac_median)) + geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.0f", 100 * ess_frac_median)), size = 3.5) +
  scale_fill_viridis_c(limits = c(0, 100)) +
  labs(x = XLAB, y = "Percent censored (%)", fill = "Median effective sample size (%)") +
  guides(fill = guide_colorbar(barwidth = unit(10, "lines"), title.vjust = 0.8)) +
  theme_minimal() + theme(panel.grid = element_blank(), legend.position = "bottom")
ggsave("figures/cs-heat-weights.png", p, height = 4.5, width = 5.5, units = "in", bg = "white")
