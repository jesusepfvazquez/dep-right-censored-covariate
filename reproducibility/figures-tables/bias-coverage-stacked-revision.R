rm(list = ls())
library(tidyverse)
library(xtable)
library(latex2exp)
library(extrafont)
`%!in%` = Negate(`%in%`)

myn = "1000-Combined-HV"
mysig = 95
zval  = 1-(1-mysig/100)/2
est_threshold = 10
#theta = c(0.1,0.1,0.1)
theta = c(1, 1, 1)
myc.p = 1.5

# Estimates for theta
read_mydata = function(x){
  # x=-1
  list_est <- list.files(path = paste0("../Simulations-AFT-Both-Est/",myn, "/results"), 
                          pattern = paste0("est_results_", x, "_seed*"), full.names = T)
  print(length(list_est))
  myresults <- list_est %>% map_df(read.csv)
  myresults$cen_group = x 
  return(myresults)
}

myresults_0 <- read_mydata(-1.5)
myresults_1 <- read_mydata(-1)
myresults_2 <- read_mydata(-0.5)
myresults_3 <- read_mydata(0.5)
myresults_4 <- read_mydata(1)
myresults_5 <- read_mydata(1.5)
myresults_6 <- read_mydata(4)
myresults_7 <- read_mydata(4.5)

myresults = rbind(myresults_0, myresults_1, myresults_2, myresults_3, 
                    myresults_4, myresults_5, myresults_6, myresults_7)

myresults_good = myresults %>%
  bind_rows() %>%
  mutate(across(c(beta0, beta1, beta2, sebeta0, sebeta1, sebeta2), as.numeric)) %>%
  subset(abs(beta0) < est_threshold & abs(beta1) < est_threshold & abs(beta2) < est_threshold) %>%
  subset(abs(sebeta0) < est_threshold & abs(sebeta1) < est_threshold & abs(sebeta2) < est_threshold)
  
  
myresults_first = myresults   %>%
  filter(sim %in% myresults_good$sim) %>%
  mutate(beta0_cov = ifelse(beta0 - qnorm(zval)*sebeta0 < theta[1] &
                          beta0 + qnorm(zval)*sebeta0 > theta[1], 1, 0),
       beta1_cov = ifelse(beta1 - qnorm(zval)*sebeta1 < theta[2] &
                          beta1 + qnorm(zval)*sebeta1 > theta[2], 1, 0),
       beta2_cov = ifelse(beta2 - qnorm(zval)*sebeta2 < theta[3] &
                          beta2 + qnorm(zval)*sebeta2 > theta[3], 1, 0)) %>%
  select(Method, beta0, beta1, beta2, sebeta0, sebeta1, sebeta2, beta0_cov, beta1_cov, beta2_cov, cen_group, pcen) %>%
  group_by(Method, cen_group) %>%
  summarise(
    n_obs = n(),
    # Compute mean and SD for beta columns
    across(
      c(beta0, beta1, beta2),
      list(mean = ~ mean(., na.rm = TRUE), sd = ~ sd(., na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    # Compute mean for sebeta columns
    across(
      c(pcen, sebeta0, sebeta1, sebeta2, beta0_cov, beta1_cov, beta2_cov),
      list(mean = ~ mean(., na.rm = TRUE)),
      .names = "{.col}_mean"
    )
  ) %>%
   select(Method, n_obs, cen_group, pcen_mean,
         #beta0_mean, beta1_mean, beta2_mean,
         #sebeta0_mean, sebeta1_mean, sebeta2_mean,
         #beta0_sd, beta1_sd, beta2_sd,
         #beta0_cov_mean, beta1_cov_mean, beta2_cov_mean)
          beta0_mean, sebeta0_mean, beta0_sd, beta0_cov_mean,
          beta1_mean, sebeta1_mean, beta1_sd, beta1_cov_mean,
          beta2_mean, sebeta2_mean, beta2_sd, beta2_cov_mean)
                              
myresults_first
         
#write.csv(myresults_first, paste0("results-sd-rate_group-",myn,".csv"), row.names = FALSE)

myresults_second = myresults %>%
    bind_rows() %>%
    mutate(across(c(beta0, beta1, beta2, sebeta0, sebeta1, sebeta2), as.numeric)) %>%
    subset(abs(beta0) < est_threshold & abs(beta1) < est_threshold & abs(beta2) < est_threshold) %>%
    subset(abs(sebeta0) < est_threshold & abs(sebeta1) < est_threshold & abs(sebeta2) < est_threshold) %>%
    group_by(Method, cen_group) %>%
    summarise(n_obs = n(),
              across(c(beta0, beta1, beta2), mean, .names = "mean_{.col}")) %>%
    mutate(mean_beta0 = (mean_beta0 - theta[1])/theta[1],
           mean_beta1 = (mean_beta1 - theta[2])/theta[2],
           mean_beta2 = (mean_beta2 - theta[3])/theta[3])

myresults_second

#write.csv(myresults_second, paste0("results-bias-rate_group-",myn,".csv"), row.names = FALSE)


########################################################################
# make a figure about the coverages 

## ---- shared recodes ----
cen_map <- function(g) {
  dplyr::case_when(
    g == -1.5 ~ 0.85, g == -1 ~ 0.75, g == -0.5 ~ 0.65,
    g ==  0.5 ~ 0.45, g == 1.5 ~ 0.25, g == 4 ~ 0.10, TRUE ~ 0.05
  )
}

method_levels <- c("oracle", "ipw-est", "aipw-est", "aipw-lambda-est", "mle-est")
method_labels <- c("Oracle", "IPW", "AIPW (unguarded)", "AIPW (guarded)", "MLE")
method_colors <- c("Oracle" = "#CC79A7", "IPW" = "#E69F00",
                   "AIPW (unguarded)" = "#56B4E9", "AIPW (guarded)" = "#009E73",
                   "MLE" = "#000000")

param_labels <- c(beta0 = "beta[0]",
                  beta1 = "beta[A-X]~(A-X)",
                  beta2 = "beta[Z]~(Z)")

## ---- bias (top panel) ----
bias_long <- myresults_second %>%
  filter(Method %in% method_levels) %>%
  select(Method, cen_group, mean_beta0, mean_beta1, mean_beta2) %>%
  pivot_longer(starts_with("mean"), names_to = "Parameter", values_to = "value") %>%
  mutate(Parameter = recode(Parameter,
                            mean_beta0 = "beta0", mean_beta1 = "beta1", mean_beta2 = "beta2"),
         panel = "Percent bias")

## ---- coverage (bottom panel) ----
cov_long <- myresults_first %>%
  filter(Method %in% method_levels) %>%
  select(Method, cen_group, beta0_cov_mean, beta1_cov_mean, beta2_cov_mean) %>%
  pivot_longer(starts_with("beta"), names_to = "Parameter", values_to = "value") %>%
  mutate(Parameter = recode(Parameter,
                            beta0_cov_mean = "beta0", beta1_cov_mean = "beta1", beta2_cov_mean = "beta2"),
         panel = "Coverage")

## ---- combine ----
plot_df <- bind_rows(bias_long, cov_long) %>%
  mutate(
    cen_group = cen_map(cen_group),
    Method    = factor(Method, levels = method_levels, labels = method_labels),
    Parameter = factor(Parameter, levels = c("beta0", "beta1", "beta2")),
    panel     = factor(panel, levels = c("Percent bias", "Coverage"))
  )

## reference lines: 0 for bias, 0.95 for coverage
## reference lines: 0 for bias, 0.95 for coverage
href <- tibble(panel = factor(c("Percent bias (%)", "Coverage (%)"),
                              levels = c("Percent bias (%)", "Coverage (%)")),
               yint  = c(0, 0.95))

plot_df <- bind_rows(bias_long, cov_long) %>%
  mutate(
    cen_group = cen_map(cen_group),
    Method    = factor(Method, levels = method_levels, labels = method_labels),
    Parameter = factor(Parameter, levels = c("beta0", "beta1", "beta2")),
    panel     = factor(panel, levels = c("Percent bias", "Coverage"),
                       labels = c("Percent bias (%)", "Coverage (%)"))
  )

p_stacked <- ggplot(plot_df, aes(x = cen_group, y = value, color = Method, group = Method)) +
  geom_hline(data = href, aes(yintercept = yint), linetype = "dotted", color = "red") +
  geom_point() +
  geom_line() +
  facet_grid(panel ~ Parameter, scales = "free_y", switch = "y",
             labeller = labeller(Parameter = as_labeller(param_labels, label_parsed),
                                 panel = label_value)) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = method_colors, name = "Estimator") +
  labs(x = "Percent censored (%)", y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.spacing = unit(0.6, "lines"),
    strip.placement = "outside",           # put row strips outside the axis
    strip.text.y.left = element_text(angle = 90)  # rotate left labels to read vertically
  )

ggsave("bias-coverage-estimators.png",
       p_stacked, height = 6, width = 8, units = "in")