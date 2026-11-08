rm(list = ls())
library(tidyverse)
library(latex2exp)

myn <- "1000-Sensitivity"
list_est <- list.files(
  path = paste0("../Simulations-Revision/", myn, "/results"),
  pattern = "^est_results_xi", full.names = TRUE
)
print(length(list_est))
myresults <- list_est %>% map_df(read.csv) %>% as_tibble()

sens <- myresults %>%
  filter(Method %in% c("ipw-est","aipw-est","aipw-lambda-est","mle-est")) %>%
  mutate(across(c(beta0, beta1, beta2, sebeta0, sebeta1, sebeta2, xi), as.numeric)) %>%
  rename(estbeta0 = beta0, estbeta1 = beta1, estbeta2 = beta2) %>%
  pivot_longer(
    cols = c(estbeta0, estbeta1, estbeta2, sebeta0, sebeta1, sebeta2),
    names_to = c(".value", "coef"),
    names_pattern = "(est|se)beta(\\d)"
  ) %>%
  mutate(coef = paste0("beta", coef)) %>%
  filter(is.finite(est), is.finite(se), se > 0, abs(est) < 5) %>%
  ## truth is 1 for every coefficient; percent bias per replicate
  mutate(truth = 1,
         pct_bias = (est - truth) / truth * 100) %>%
  group_by(Method, coef, xi) %>%
  summarise(
    bias = mean(pct_bias, na.rm = TRUE),
    pcen = mean(pcen, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(xi <= 1) %>%   # drop cells where too many replicates diverged
  mutate(
    Method  = recode(Method, "ipw-est"="IPW","aipw-est"="AIPW (unguarded)",
                     "aipw-lambda-est"="AIPW (guarded)","mle-est"="MLE"),
    Method  = factor(Method, levels = c("IPW", "AIPW (unguarded)", "AIPW (guarded)", "MLE")),
    coef = recode(coef, "beta0"="beta[0]",
                  "beta1"="beta[A-X]~(A-X)","beta2"="beta[Z]~(Z)")
  )

sens %>% filter(coef == "beta[0]") %>% print(n = 100)

p1 = ggplot(sens, aes(x = xi, y = bias, color = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  facet_grid(. ~ coef, labeller = label_parsed) +
  scale_color_manual(
    name = NULL,
    breaks = c("IPW", "AIPW (unguarded)", "AIPW (guarded)", "MLE"),
    labels = c("IPW", "AIPW (unguarded)", "AIPW (guarded)", "MLE"),
    values = c("IPW" = "#E69F00",
               "AIPW (unguarded)" = "#56B4E9",
               "AIPW (guarded)" = "#009E73",
               "MLE" = "#000000")
  ) +
  labs(
    x = expression("Sensitivity parameter " * xi),
    y = "Percent bias (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.spacing = unit(0.8, "lines"))

ggsave("sensitivity.png", p1, height = 4, width = 8, units = "in")
