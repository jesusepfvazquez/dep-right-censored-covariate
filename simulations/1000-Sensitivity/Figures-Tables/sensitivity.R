## sensitivity.R
## Figure sensitivity.png and its table, from ../results (this study).
## Run from 1000-Sensitivity/Figures-Tables/.
library(tidyverse)

res <- list.files("../results", "^est_results_xi[0-9.]+_seed_[0-9]+\\.csv$",
                  full.names = TRUE) %>% map_df(read.csv) %>% as_tibble()

methods <- c("ipw-est" = "IPW", "aipw-est-joint" = "AIPW (unguarded)",
             "aipw-lambda-est" = "AIPW (guarded)", "mle-joint" = "MLE (joint)")
z <- qnorm(0.975)

sens <- res %>%
  filter(Method %in% names(methods)) %>%
  rename(estbeta0 = beta0, estbeta1 = beta1, estbeta2 = beta2) %>%
  pivot_longer(c(estbeta0, estbeta1, estbeta2, sebeta0, sebeta1, sebeta2),
               names_to = c(".value", "coef"), names_pattern = "(est|se)beta(\\d)") %>%
  mutate(coef = paste0("beta", coef)) %>%
  filter(is.finite(est), abs(est) < 5) %>%
  group_by(Method, coef, xi) %>%
  summarise(reps = n(), pcen = mean(pcen), bias = 100 * mean(est - 1), sd = 100 * sd(est),
            coverage = 100 * mean(abs(est - 1) < z * se, na.rm = TRUE),   # NA when do_var = FALSE
            .groups = "drop") %>%
  mutate(Method = factor(recode(Method, !!!methods), levels = unname(methods)))

dir.create("tables", showWarnings = FALSE)
write.csv(sens, "tables/sensitivity.csv", row.names = FALSE)
sens %>% filter(coef != "beta2") %>% arrange(coef, xi, Method) %>%
  mutate(across(where(is.double), ~ round(.x, 2))) %>% print(n = Inf)

p <- sens %>%
  mutate(coef = recode(coef, beta0 = "beta[0]", beta1 = "beta[A-X]~(A-X)", beta2 = "beta[Z]~(Z)")) %>%
  ggplot(aes(x = xi, y = bias, color = Method)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  facet_grid(. ~ coef, labeller = label_parsed) +
  scale_color_manual(name = NULL, values = c("IPW" = "#E69F00", "AIPW (unguarded)" = "#56B4E9",
                                             "AIPW (guarded)" = "#009E73", "MLE (joint)" = "#000000")) +
  labs(x = expression("Sensitivity parameter " * xi), y = "Percent bias (%)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.spacing = unit(0.8, "lines"))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/sensitivity.png", p, height = 4, width = 8, units = "in", bg = "white")
