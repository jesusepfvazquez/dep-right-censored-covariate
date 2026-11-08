rm(list = ls())
library(tidyverse)
library(latex2exp)


## Censoring bins: left-closed, right-open, [a, b); the last bin is [90, 100].
## Round first so that floating-point error does not move boundary replicates
## (e.g., in R, 0.7 < 7 * 0.1). No reported replicate falls at exactly 90%.
cens_bin10 <- function(p) {
  pct <- round(100 * p, 6)
  lab <- c("[0, 10)", "[10, 20)", "[20, 30)", "[30, 40)", "[40, 50)",
           "[50, 60)", "[60, 70)", "[70, 80)", "[80, 90)", "[90, 100]")
  idx <- pmin(floor(pct / 10) + 1, 10)
  factor(lab[idx], levels = lab)
}

myn <- "1000-Combined-HV"

list_est <- list.files(
  path = paste0("../Simulations-Revision/", myn, "/results"),
  pattern = "^est_results_n", full.names = TRUE
)
print(length(list_est))

myresults <- list_est %>% map_df(read.csv) %>% as_tibble()

## pivot the three SE columns to long form
se_long <- myresults %>%
  mutate(across(c(sebeta0, sebeta1, sebeta2, pcen, n, group, sim), as.numeric)) %>%
  filter(Method %in% c("ipw-est","mle-est","aipw-est","aipw-lambda-est")) %>%
  select(Method, sim, group, n, pcen, sebeta0, sebeta1, sebeta2) %>%
  pivot_longer(c(sebeta0, sebeta1, sebeta2), names_to = "coef", values_to = "se") %>%
  mutate(coef = recode(coef, "sebeta0"="beta0","sebeta1"="beta1","sebeta2"="beta2")) %>%
  filter(is.finite(se), se > 0)

## paired SE ratio vs IPW within each replicate
ipw_se <- se_long %>% filter(Method == "ipw-est") %>%
  select(sim, group, n, pcen, coef, ipw_se = se)

eff <- se_long %>%
  filter(Method %in% c("mle-est","aipw-est","aipw-lambda-est")) %>%
  inner_join(ipw_se, by = c("sim","group","n","pcen","coef")) %>%
  mutate(se_ratio = se / ipw_se) %>%
  filter(is.finite(se_ratio), se_ratio < 5)

## bin into (censoring rate, n) cells and average the ratio per cell
eff_cells <- eff %>%
  mutate(
    cen_bin = cens_bin10(pcen),
    n_bin   = cut(n, breaks = seq(100, 1000, by = 100),
                  include.lowest = TRUE, right = FALSE,
                  labels = c("[100, 200)","[200, 300)","[300, 400)","[400, 500)",
                             "[500, 600)","[600, 700)","[700, 800)","[800, 900)",
                             "[900, 1000]"))
  ) %>%
  group_by(Method, coef, cen_bin, n_bin) %>%
  summarise(eff_ratio = mean(se_ratio, na.rm = TRUE),
            nsim = n(), .groups = "drop") %>%
  filter(nsim >= 100, cen_bin %in% levels(cen_bin)[3:9]) %>%
  mutate(
    Method  = recode(Method, "aipw-est"="AIPW (unguarded)","aipw-lambda-est"="AIPW (guarded)","mle-est"="MLE"),
    Method  = factor(Method, levels = c("AIPW (unguarded)", "AIPW (guarded)", "MLE")),
    coef    = recode(coef, "beta0"="beta[0]",
                     "beta1"="beta[A-X]~(A-X)","beta2"="beta[Z]~(Z)")
  )

## heatmap: x = censoring rate, y = n, fill = efficiency ratio; estimators in rows, betas in columns
p1 = ggplot(eff_cells, aes(x = cen_bin, y = n_bin, fill = eff_ratio)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", eff_ratio)), size = 2.5,
            color = "grey20") +
  facet_grid(Method ~ coef,
           labeller = labeller(coef = label_parsed, Method = label_value)) +
  scale_fill_gradient2(
    midpoint = 1, low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    name = "Efficiency ratio (vs IPW)"
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "Percent censored (%)", y = "Sample size (n)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.spacing = unit(0.6, "lines"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = guide_colorbar(
    barwidth  = 15,   # length of the bar (default ~ 15)
    barheight = 1   # thickness (default ~ 0.5)
  )) 

ggsave("heatmap-efficiency.png", p1, height = 10, width = 10, units = "in")

## ---- BIAS HEATMAP (addition) ----
theta <- c(1, 1, 1)   # truth for beta0, beta1, beta2
est_long <- myresults %>%
  mutate(across(c(beta0, beta1, beta2, pcen, n, group, sim), as.numeric)) %>%
  filter(Method %in% c("ipw-est","mle-est","aipw-est","aipw-lambda-est")) %>%
  select(Method, sim, group, n, pcen, beta0, beta1, beta2) %>%
  pivot_longer(c(beta0, beta1, beta2), names_to = "coef", values_to = "est") %>%
  filter(is.finite(est)) %>%
  mutate(
    truth = case_when(coef == "beta0" ~ theta[1],
                      coef == "beta1" ~ theta[2],
                      coef == "beta2" ~ theta[3]),
    pct_bias = (est - truth) / truth * 100
  ) %>%
  filter(abs(pct_bias) < 200)

bias_cells <- est_long %>%
  mutate(
    cen_bin = cens_bin10(pcen),
    n_bin   = cut(n, breaks = seq(100, 1000, by = 100),
                  include.lowest = TRUE, right = FALSE,
                  labels = c("[100, 200)","[200, 300)","[300, 400)","[400, 500)",
                             "[500, 600)","[600, 700)","[700, 800)","[800, 900)",
                             "[900, 1000]"))
  ) %>%
  group_by(Method, coef, cen_bin, n_bin) %>%
  summarise(bias = abs(mean(pct_bias, na.rm = TRUE)),
            nsim = n(), .groups = "drop") %>%
  filter(nsim >= 100, cen_bin %in% levels(cen_bin)[3:9]) %>%
  mutate(
    Method  = recode(Method, "ipw-est"="IPW","aipw-est"="AIPW (unguarded)",
                     "aipw-lambda-est"="AIPW (guarded)","mle-est"="MLE"),
    Method  = factor(Method, levels = c("IPW", "AIPW (unguarded)", "AIPW (guarded)", "MLE")),
    coef    = recode(coef, "beta0"="beta[0]",
                     "beta1"="beta[A-X]~(A-X)","beta2"="beta[Z]~(Z)")
  )

p2 = ggplot(bias_cells, aes(x = cen_bin, y = n_bin, fill = bias)) +
  geom_tile() +
  geom_text(aes(label = round(bias, 0)), size = 2.5,
            color = "grey20") +
  facet_grid(Method ~ coef,
           labeller = labeller(coef = label_parsed, Method = label_value)) +
  scale_fill_gradientn(
    colours = c("#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
    limits = c(0, 30), oob = scales::squish,
    name = "Absolute percent bias (%)"
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "Percent censored (%)", y = "Sample size (n)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.spacing = unit(0.6, "lines"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = guide_colorbar(
    barwidth  = 15,   # length of the bar (default ~ 15)
    barheight = 1   # thickness (default ~ 0.5)
  )) 

ggsave("heatmap-bias.png", p2, height = 12, width = 10, units = "in")

## ---- COVERAGE HEATMAP (addition) ----
theta <- c(1, 1, 1)
zval  <- qnorm(1 - (1 - 0.95)/2)   # 1.96

cov_long <- myresults %>%
  mutate(across(c(beta0, beta1, beta2, sebeta0, sebeta1, sebeta2,
                  pcen, n, group, sim), as.numeric)) %>%
  filter(Method %in% c("ipw-est","mle-est","aipw-est","aipw-lambda-est")) %>%
  select(Method, sim, group, n, pcen,
         beta0, beta1, beta2, sebeta0, sebeta1, sebeta2) %>%
  rename(estbeta0 = beta0, estbeta1 = beta1, estbeta2 = beta2) %>%
  pivot_longer(
    cols = c(estbeta0, estbeta1, estbeta2, sebeta0, sebeta1, sebeta2),
    names_to = c(".value", "coef"),
    names_pattern = "(est|se)beta(\\d)"
  ) %>%
  mutate(coef = paste0("beta", coef)) %>%
  filter(is.finite(est), is.finite(se), se > 0) %>%
  mutate(
    truth = case_when(coef == "beta0" ~ theta[1],
                      coef == "beta1" ~ theta[2],
                      coef == "beta2" ~ theta[3]),
    covered = as.numeric(truth >= est - zval*se & truth <= est + zval*se)
  )

cov_cells <- cov_long %>%
  mutate(
   cen_bin = cens_bin10(pcen),
    n_bin   = cut(n, breaks = seq(100, 1000, by = 100),
                  include.lowest = TRUE, right = FALSE,
                  labels = c("[100, 200)","[200, 300)","[300, 400)","[400, 500)",
                             "[500, 600)","[600, 700)","[700, 800)","[800, 900)",
                             "[900, 1000]"))
  ) %>%
  group_by(Method, coef, cen_bin, n_bin) %>%
  summarise(coverage = mean(covered, na.rm = TRUE) * 100,
            nsim = n(), .groups = "drop") %>%
  filter(nsim >= 100, cen_bin %in% levels(cen_bin)[3:9]) %>%
  mutate(
    Method  = recode(Method, "ipw-est"="IPW","aipw-est"="AIPW (unguarded)",
                     "aipw-lambda-est"="AIPW (guarded)","mle-est"="MLE"),
    Method  = factor(Method, levels = c("IPW", "AIPW\ (unguarded)", "AIPW (guarded)", "MLE")),
    coef    = recode(coef, "beta0"="beta[0]",
                     "beta1"="beta[A-X]~(A-X)","beta2"="beta[Z]~(Z)")
  )

p3 = ggplot(cov_cells, aes(x = cen_bin, y = n_bin, fill = coverage)) +
  geom_tile() +
  geom_text(aes(label = round(coverage, 0)), size = 2.5,
            color = "grey20") +
  facet_grid(Method ~ coef,
           labeller = labeller(coef = label_parsed, Method = label_value)) +
  scale_fill_gradientn(
      colours = c("#67001F",  # 30 and below: deep red
                  "#B2182B",  # 50
                  "#F4A582",  # 70
                  "#F7F7F7",  # 90: neutral starts
                  "#F7F7F7"), # 100: still neutral
      values = scales::rescale(c(30, 50, 70, 90, 100), from = c(30, 100)),
      limits = c(30, 100), oob = scales::squish,
      name = "95% CI coverage (%)"
    ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "Percent censored (%)", y = "Sample size (n)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.spacing = unit(0.6, "lines"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = guide_colorbar(
    barwidth  = 15,   # length of the bar (default ~ 15)
    barheight = 1   # thickness (default ~ 0.5)
  )) 

ggsave("heatmap-coverage.png", p3, height = 12, width = 10, units = "in")

## ---- COUNTS HEATMAP (addition) ----
count_cells <- cov_long %>%
  mutate(
   cen_bin = cens_bin10(pcen),
    n_bin   = cut(n, breaks = seq(100, 1000, by = 100),
                  include.lowest = TRUE, right = FALSE,
                  labels = c("[100, 200)","[200, 300)","[300, 400)","[400, 500)",
                             "[500, 600)","[600, 700)","[700, 800)","[800, 900)",
                             "[900, 1000]"))
  ) %>%
  group_by(Method, coef, cen_bin, n_bin) %>%
  summarise(coverage = mean(covered, na.rm = TRUE) * 100,
            nsim = n(), .groups = "drop") %>%
  filter(nsim >= 100, cen_bin %in% levels(cen_bin)[3:9]) %>%
  mutate(
    Method  = recode(Method, "ipw-est"="IPW","aipw-est"="AIPW (unguarded)",
                     "aipw-lambda-est"="AIPW (guarded)","mle-est"="MLE"),
    Method  = factor(Method, levels = c("IPW", "AIPW (unguarded)", "AIPW (guarded)", "MLE")),
    coef    = recode(coef, "beta0"="beta[0]",
                     "beta1"="beta[A-X]~(A-X)","beta2"="beta[Z]~(Z)")
  )

p4 = ggplot(count_cells, aes(x = cen_bin, y = n_bin, fill = nsim)) +
  geom_tile() +
  geom_text(aes(label = nsim), size = 2.5, color = "grey20") +
  facet_grid(Method ~ coef,
           labeller = labeller(coef = label_parsed, Method = label_value)) +
  scale_fill_gradient(
    low = "#F7F7F7", high = "#2166AC",
    limits = c(100, max(count_cells$nsim)),
    breaks = c(100, max(count_cells$nsim)),
    name = "Replicates per cell"
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "Percent censored (%)", y = "Sample size (n)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.spacing = unit(0.6, "lines"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = guide_colorbar(
    barwidth  = 15,   # length of the bar (default ~ 15)
    barheight = 1   # thickness (default ~ 0.5)
  )) 

ggsave("heatmap-counts.png", p4, height = 12, width = 10, units = "in")
