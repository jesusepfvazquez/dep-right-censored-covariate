## fig1-bias-coverage.R
## Figure 1 (02-bias-coverage-estimators.png) and its table, from the first
## simulation study (n = 1000, eight censoring intercepts, seeds 0-2999).
## All estimators come from ../results (this study's main.R), except the two-step MLE,
## which comes from ../Fig1-TwoStep/results (same seeds and datasets).
## Run from 1000-Main/Figures-Tables/.
library(tidyverse)
library(ggh4x)   # facetted_pos_scales(): fixed y-limits per row

res_dir <- "../results"
twostep_dir <- "../Fig1-TwoStep/results"   # two-step MLE add-on; skipped if absent
ess_dir     <- "../Fig1-TwoStep/results-ess"   # ESS of the weights 1/pi-hat (ess.R); skipped if absent
BIAS_LIM <- c(-30, 5)     # y-limits of the bias row
COV_LIM  <- c(40, 100)    # y-limits of the coverage row; the ESS shading spans this range (0-100%)
est_threshold <- 10
z <- qnorm(0.975)

## censoring intercept taken from the file name: est_results_<c_int>_seed_<seed>.csv
rd <- function(dir) list.files(dir, "^est_results_-?[0-9.]+_seed_[0-9]+\\.csv$", full.names = TRUE) %>%
  map_df(~ read.csv(.x) %>%
           select(Method, beta0, beta1, beta2, sebeta0, sebeta1, sebeta2, pcen, sim) %>%
           mutate(across(-Method, as.numeric),
                  group = as.numeric(sub("^est_results_(-?[0-9.]+)_seed_.*", "\\1", basename(.x)))))

res <- rd(res_dir) %>% filter(Method %in% c("oracle", "ipw-est", "mle-joint", "aipw-est-joint", "aipw-lambda-est"))
if (dir.exists(twostep_dir)) res <- bind_rows(res, rd(twostep_dir) %>% filter(Method == "mle-twostep"))

method_levels <- c("oracle", "ipw-est", "aipw-est-joint", "aipw-lambda-est", "mle-joint", "mle-twostep")
method_labels <- c("Oracle", "IPW", "AIPW (unguarded)", "AIPW (guarded)", "MLE (joint)", "MLE (two-step)")
method_colors <- c("Oracle" = "#CC79A7", "IPW" = "#E69F00", "AIPW (unguarded)" = "#56B4E9",
                   "AIPW (guarded)" = "#009E73", "MLE (joint)" = "#000000", "MLE (two-step)" = "#3A5FCD")   # royal blue; grey is the ESS shading
method_lines  <- c("Oracle" = "solid", "IPW" = "solid", "AIPW (unguarded)" = "solid",
                   "AIPW (guarded)" = "solid", "MLE (joint)" = "solid", "MLE (two-step)" = "dashed")

long <- res %>%
  filter(if_all(c(beta0, beta1, beta2), ~ is.finite(.x) & abs(.x) < est_threshold)) %>%
  rename(estbeta0 = beta0, estbeta1 = beta1, estbeta2 = beta2) %>%
  pivot_longer(c(estbeta0, estbeta1, estbeta2, sebeta0, sebeta1, sebeta2),
               names_to = c(".value", "Parameter"), names_pattern = "(est|se)beta(\\d)") %>%
  mutate(Parameter = paste0("beta", Parameter))

## one rate per intercept: the mean realized censoring rate (not a hard-coded map)
rate <- long %>% group_by(group) %>% summarise(rate = mean(pcen), .groups = "drop")

tab <- long %>%
  group_by(Method, group, Parameter) %>%
  summarise(reps = n(), bias = 100 * mean(est - 1),
            coverage = 100 * mean(abs(est - 1) < z * se, na.rm = TRUE),   # per-replicate SE
            sd = 100 * sd(est), se = 100 * mean(se, na.rm = TRUE),         # x 100, as in Table 2
            mse = 1e4 * mean((est - 1)^2),                                 # x 100^2
            .groups = "drop") %>%
  left_join(rate, by = "group") %>%
  mutate(Method = factor(Method, levels = method_levels, labels = method_labels)) %>%
  arrange(Parameter, desc(rate), Method)

dir.create("tables", showWarnings = FALSE)
write.csv(tab, "tables/fig1-bias-coverage.csv", row.names = FALSE)
tab %>% filter(Parameter %in% c("beta0", "beta1")) %>%
  mutate(across(where(is.double), ~ round(.x, 2))) %>% print(n = Inf)

## results at one censoring rate (default: c_int = -0.5, about 66% censoring),
## laid out like Table 2: bias, SE, SD (x 100) and coverage per coefficient
one_rate <- function(c_int = -0.5) {
  tab %>% filter(group == c_int) %>%
    mutate(rate = round(100 * rate, 1), across(c(bias, se, sd, coverage), ~ round(.x, 2))) %>%
    select(Method, rate, Parameter, bias, se, sd, coverage) %>%
    pivot_wider(names_from = Parameter, values_from = c(bias, se, sd, coverage), names_vary = "slowest") %>%
    arrange(Method)
}
tab66 <- one_rate(-0.5)
print(tab66, width = Inf)
write.csv(tab66, "tables/fig1-cens66.csv", row.names = FALSE)

## points outside the fixed y-limits below are dropped from the plot: list them
off <- tab %>% filter(bias < BIAS_LIM[1] | bias > BIAS_LIM[2] | coverage < COV_LIM[1])
if (nrow(off)) { message("Outside the y-limits (not plotted):"); print(select(off, Method, rate, Parameter, bias, coverage)) }

plot_df <- tab %>%
  select(Method, rate, Parameter, bias, coverage) %>%
  pivot_longer(c(bias, coverage), names_to = "panel") %>%
  mutate(panel = factor(panel, levels = c("bias", "coverage"),
                        labels = c("Percent bias (%)", "Coverage (%)")),
         Parameter = factor(Parameter, levels = c("beta0", "beta1", "beta2")))
href <- tibble(panel = factor(c("Percent bias (%)", "Coverage (%)"),
                              levels = c("Percent bias (%)", "Coverage (%)")), yint = c(0, 95))
param_labels <- c(beta0 = "beta[0]", beta1 = "beta[A-X]~(A-X)", beta2 = "beta[Z]~(Z)")

## median ESS of the weights 1/pi-hat (used by the IPW and AIPW estimators) (as a fraction of the uncensored n) per censoring level,
## drawn as a grey shaded area in both rows and read on the right axis (0-100%): in each row
## ESS = 0% is the bottom of the y-limits and ESS = 100% the top
ESS_NAME <- "Effective sample size (%)"
ESS_GREY <- "grey75"          # shading colour (alpha below)
ESS_LABEL_ALL <- FALSE        # TRUE: "ESS" label in every panel; FALSE: only the first panel
ESS_LABEL_X   <- 19           # centre of the "ESS" label, in percent censored
ESS_LABEL_Y   <- -22.5        # centre of the label in the bias row (same relative height in the coverage row)
ESS_LABEL_SIZE <- 7           # ggplot text size (about 20 pt)
ess_df <- if (dir.exists(ess_dir)) {
  list.files(ess_dir, "^ess_.*\\.csv$", full.names = TRUE) %>% map_df(read.csv) %>%
    group_by(group) %>% summarise(ess = median(ess_frac), .groups = "drop") %>%
    left_join(rate, by = "group") %>%
    crossing(Parameter = factor(c("beta0", "beta1", "beta2"), levels = c("beta0", "beta1", "beta2")),
             panel = factor(c("Percent bias (%)", "Coverage (%)"), levels = c("Percent bias (%)", "Coverage (%)"))) %>%
    mutate(lo = ifelse(panel == "Coverage (%)", COV_LIM[1], BIAS_LIM[1]),
           y  = lo + ifelse(panel == "Coverage (%)", diff(COV_LIM), diff(BIAS_LIM)) * ess)
} else NULL
if (!is.null(ess_df)) print(distinct(ess_df, group, rate, ess) %>% arrange(rate))
ess_lab <- if (!is.null(ess_df)) {
  f <- (ESS_LABEL_Y - BIAS_LIM[1]) / diff(BIAS_LIM)
  l <- distinct(ess_df, Parameter, panel) %>%
    mutate(y = ifelse(panel == "Coverage (%)", COV_LIM[1] + f * diff(COV_LIM), ESS_LABEL_Y))
  if (ESS_LABEL_ALL) l else filter(l, Parameter == "beta0", panel == "Percent bias (%)")
} else NULL

p <- ggplot(plot_df, aes(x = 100 * rate, y = value, color = Method, linetype = Method, group = Method)) +
  { if (!is.null(ess_df)) list(
      geom_ribbon(data = ess_df, aes(x = 100 * rate, ymin = lo, ymax = y),
                  inherit.aes = FALSE, fill = ESS_GREY, alpha = 0.12),
      geom_line(data = ess_df, aes(x = 100 * rate, y = y), inherit.aes = FALSE,
                color = ESS_GREY, linewidth = 0.3),
      geom_point(data = ess_df, aes(x = 100 * rate, y = y), inherit.aes = FALSE,
                 color = ESS_GREY, shape = 15, size = 0.8),
      ## "ESS" label centred at (ESS_LABEL_X, ESS_LABEL_Y), between the -20 and -25 bias lines
      ## (first panel only unless ESS_LABEL_ALL)
      geom_text(data = ess_lab, aes(x = ESS_LABEL_X, y = y, label = "ESS"), inherit.aes = FALSE,
                hjust = 0.5, vjust = 0.5, size = ESS_LABEL_SIZE, fontface = "bold", color = "grey30")) } +
  geom_hline(data = href, aes(yintercept = yint), linetype = "dotted", color = "red") +
  geom_point() + geom_line() +
  facet_grid(panel ~ Parameter, scales = "free_y", switch = "y",
             labeller = labeller(Parameter = as_labeller(param_labels, label_parsed), panel = label_value)) +
  scale_color_manual(values = method_colors, name = "Estimator") +
  scale_linetype_manual(values = method_lines, name = "Estimator") +
  labs(x = "Percent censored (%)", y = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", panel.spacing = unit(0.6, "lines"),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 90)) +
  facetted_pos_scales(y = list(
    ## the right-axis title is drawn from the first panel's scale (the bias row), so both rows
    ## carry the same name
    panel == "Percent bias (%)" ~ scale_y_continuous(limits = BIAS_LIM,
      sec.axis = if (!is.null(ess_df)) sec_axis(~ 100 * (. - BIAS_LIM[1]) / diff(BIAS_LIM), name = ESS_NAME,
                                                 breaks = c(0, 25, 50, 75, 100)) else waiver()),
    panel == "Coverage (%)"     ~ scale_y_continuous(limits = COV_LIM,
      sec.axis = if (!is.null(ess_df)) sec_axis(~ 100 * (. - COV_LIM[1]) / diff(COV_LIM), name = ESS_NAME,
                                                 breaks = c(0, 25, 50, 75, 100)) else waiver())))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/02-bias-coverage-estimators.png", p, height = 6, width = 8, units = "in", bg = "white")


## ------------------------------------------------------------------------------------
## Version with a third row: MSE ratio against the IPW estimator,
## RE = MSE(IPW) / MSE(estimator), from the empirical MSE (bias^2 + variance), so a biased
## estimator such as the two-step MLE is penalized. Ratio > 1: smaller MSE than the IPW
## estimator. (The SE version is re = se[Method == "IPW"] / se.) Log2 scale with fixed limits RE_LIM; the ESS shading is drawn in all three rows
## (in the RE row, ESS = 0% is RE_LIM[1] and 100% is RE_LIM[2], evenly spaced on the log scale).
## Output: figures/02-bias-coverage-re.png (the two-row figure above is unchanged).
## ------------------------------------------------------------------------------------
RE_LIM    <- c(0.25, 32)   # y-limits of the RE row (log2); points outside are listed below
RE_BREAKS <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32)
ess_to_re <- function(e) 2^(log2(RE_LIM[1]) + e * diff(log2(RE_LIM)))
re_to_ess <- function(r) 100 * (log2(r) - log2(RE_LIM[1])) / diff(log2(RE_LIM))
lv3 <- c("Percent bias (%)", "Coverage (%)", "MSE ratio")
re_tab <- tab %>% group_by(group, Parameter) %>%
  mutate(re = mse[Method == "IPW"] / mse) %>% ungroup()
plot_df3 <- re_tab %>%
  select(Method, rate, Parameter, bias, coverage, re) %>%
  pivot_longer(c(bias, coverage, re), names_to = "panel") %>%
  mutate(panel = factor(panel, levels = c("bias", "coverage", "re"), labels = lv3),
         Parameter = factor(Parameter, levels = c("beta0", "beta1", "beta2")))
href3 <- tibble(panel = factor(lv3, levels = lv3), yint = c(0, 95, 1))
ess_df3  <- if (!is.null(ess_df)) {
  bind_rows(ess_df %>% mutate(panel = as.character(panel)),
            ess_df %>% filter(panel == "Coverage (%)") %>%
              mutate(panel = "MSE ratio", lo = RE_LIM[1], y = ess_to_re(ess))) %>%
    mutate(panel = factor(panel, levels = lv3))
} else NULL
ess_lab3 <- if (!is.null(ess_lab)) mutate(ess_lab, panel = factor(as.character(panel), levels = lv3)) else NULL
off3 <- re_tab %>% filter(re < RE_LIM[1] | re > RE_LIM[2])
if (nrow(off3)) { message("RE outside RE_LIM (not plotted):"); print(select(off3, Method, rate, Parameter, re)) }
print(re_tab %>% filter(Parameter == "beta0") %>% select(Method, rate, re) %>%
        mutate(rate = round(100 * rate), re = round(re, 2)) %>%
        pivot_wider(names_from = Method, values_from = re) %>% arrange(rate), width = Inf)

p3 <- ggplot(plot_df3, aes(x = 100 * rate, y = value, color = Method, linetype = Method, group = Method)) +
  { if (!is.null(ess_df3)) list(
      geom_ribbon(data = ess_df3, aes(x = 100 * rate, ymin = lo, ymax = y),
                  inherit.aes = FALSE, fill = ESS_GREY, alpha = 0.12),
      geom_line(data = ess_df3, aes(x = 100 * rate, y = y), inherit.aes = FALSE,
                color = ESS_GREY, linewidth = 0.3),
      geom_point(data = ess_df3, aes(x = 100 * rate, y = y), inherit.aes = FALSE,
                 color = ESS_GREY, shape = 15, size = 0.8),
      geom_text(data = ess_lab3, aes(x = ESS_LABEL_X, y = y, label = "ESS"), inherit.aes = FALSE,
                hjust = 0.5, vjust = 0.5, size = ESS_LABEL_SIZE, fontface = "bold", color = "grey30")) } +
  geom_hline(data = href3, aes(yintercept = yint), linetype = "dotted", color = "red") +
  geom_point() + geom_line() +
  facet_grid(panel ~ Parameter, scales = "free_y", switch = "y",
             labeller = labeller(Parameter = as_labeller(param_labels, label_parsed), panel = label_value)) +
  scale_color_manual(values = method_colors, name = "Estimator") +
  scale_linetype_manual(values = method_lines, name = "Estimator") +
  labs(x = "Percent censored (%)", y = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", panel.spacing = unit(0.6, "lines"),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 90)) +
  facetted_pos_scales(y = list(
    panel == "Percent bias (%)" ~ scale_y_continuous(limits = BIAS_LIM,
      sec.axis = if (!is.null(ess_df)) sec_axis(~ 100 * (. - BIAS_LIM[1]) / diff(BIAS_LIM), name = ESS_NAME,
                                                 breaks = c(0, 25, 50, 75, 100)) else waiver()),
    panel == "Coverage (%)"     ~ scale_y_continuous(limits = COV_LIM,
      sec.axis = if (!is.null(ess_df)) sec_axis(~ 100 * (. - COV_LIM[1]) / diff(COV_LIM), name = ESS_NAME,
                                                 breaks = c(0, 25, 50, 75, 100)) else waiver()),
    panel == "MSE ratio" ~ scale_y_continuous(trans = "log2", limits = RE_LIM, breaks = RE_BREAKS,
                                                        labels = function(x) as.character(x),
      sec.axis = if (!is.null(ess_df)) sec_axis(~ re_to_ess(.), name = ESS_NAME,
                                                 breaks = c(0, 25, 50, 75, 100)) else waiver())))

ggsave("figures/02-bias-coverage-re.png", p3, height = 8.5, width = 8, units = "in", bg = "white")
