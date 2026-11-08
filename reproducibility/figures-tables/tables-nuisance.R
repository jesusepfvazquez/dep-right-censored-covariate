rm(list = ls())
library(tidyverse)
library(xtable)
library(latex2exp)
`%!in%` = Negate(`%in%`)

myn = "500-Nuissance"
mysig = 95
zval  = 1-(1-mysig/100)/2
est_threshold = 10
#theta = c(0.1,0.1,0.1)
# alpha = c(0.1, 0.1, 0.8)
# theta = c(1, 0.1, 0.1, 2)


# Estimates for gamma 
# X|Z

my_alpha_cov = function(x){
  
  theta = c(0.1, 0.1, 0.8)
  #x = -2
  list_est <- list.files(path = paste0("../Simulations-AFT-Both-Est/",myn, "/results"), 
                          pattern = paste0("est_results_gamma_", x, "_*"), full.names = T)
  print(length(list_est))
  myresults <- list_est %>% map_df(read.csv)
  myresults$cen_group = x 
  
  myresults_first = myresults %>%
    bind_rows() %>%
    mutate(across(c(alpha0, alpha1, alpha2, alpha0se, alpha1se, alpha2se), as.numeric)) %>%
    subset(abs(alpha0) < est_threshold & abs(alpha1) < est_threshold & abs(alpha2) < est_threshold) %>%
    subset(abs(alpha0se) < est_threshold & abs(alpha1se) < est_threshold & abs(alpha2se) < est_threshold) %>%
    mutate(alpha0_cov = ifelse(alpha0 - qnorm(zval)*alpha0se < theta[1] &
                            alpha0 + qnorm(zval)*alpha0se > theta[1], 1, 0),
         alpha1_cov = ifelse(alpha1 - qnorm(zval)*alpha1se < theta[2] &
                            alpha1 + qnorm(zval)*alpha1se > theta[2], 1, 0),
         alpha2_cov = ifelse(alpha2 - qnorm(zval)*alpha2se < theta[3] &
                            alpha2 + qnorm(zval)*alpha2se > theta[3], 1, 0)) %>%
    select(alpha0:alpha2, alpha0se:alpha2se, alpha0_cov:alpha2_cov, cen_group) %>%
    group_by(cen_group) %>%
    summarise(
      n_obs = n(),
      # Compute mean and SD for beta columns
      across(
        c(alpha0, alpha1, alpha2),
        list(mean = ~ mean(.), sd = ~ sd(.)),
        .names = "{.col}_{.fn}"
      ),
      # Compute mean for sebeta columns
      across(
        c(alpha0se:alpha2se, alpha0_cov:alpha2_cov),
        mean,
        .names = "{.col}_mean"
      )
    ) 
          
  return(myresults_first)
}

my_gamma_cov_1 <- my_alpha_cov(-2)
my_gamma_cov_2 <- my_alpha_cov(-1)
my_gamma_cov_3 <- my_alpha_cov(0)
my_gamma_cov_4 <- my_alpha_cov(1)
my_gamma_cov_5 <- my_alpha_cov(2)
my_gamma_cov_6 <- my_alpha_cov(3)
my_gamma_cov_7 <- my_alpha_cov(4)

myresults_gamma_cov = rbind(my_gamma_cov_1, my_gamma_cov_2, my_gamma_cov_3, my_gamma_cov_4, my_gamma_cov_5, my_gamma_cov_6, my_gamma_cov_7)
# head(myresults_gamma_cov)

# make a figure about the coverages 


p = myresults_gamma_cov %>%
  select(cen_group, alpha0_cov_mean, alpha1_cov_mean, alpha2_cov_mean) %>%
  pivot_longer(
    cols = starts_with("alpha"),
    names_to = "Parameter",
    values_to = "Coverage"
  ) %>%
  #mutate(Parameter = factor(Parameter, labels = c("alpha0", "beta1", "beta2"))) %>%
  ggplot(aes(x = cen_group, y = Coverage)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "red") +
  facet_wrap(~ Parameter) +
  labs(
    title = "Coverage by Method and Censoring Group",
    x = "Censoring Group",
    y = "Coverage"
    #color = "Method"
  ) +
  ylim(0.5, 1)+
  theme_minimal()
p 
  
# ggsave(paste0("results-efficiency-rate_group-",myn,".png"), p, height = 4, width = 8, units = "in")
  
my_alpha_bias = function(x){
  
  theta = c(0.1, 0.1, 0.8)
  #x = -2
  list_est <- list.files(path = paste0("../Simulations-AFT-Both-Est/",myn, "/results"), 
                          pattern = paste0("est_results_gamma_", x, "_*"), full.names = T)
  print(length(list_est))
  myresults <- list_est %>% map_df(read.csv)
  myresults$cen_group = x 
  
  myresults_first = myresults %>%
    bind_rows() %>%
    mutate(across(c(alpha0, alpha1, alpha2, alpha0se, alpha1se, alpha2se), as.numeric)) %>%
    subset(abs(alpha0) < est_threshold & abs(alpha1) < est_threshold & abs(alpha2) < est_threshold) %>%
    subset(abs(alpha0se) < est_threshold & abs(alpha1se) < est_threshold & abs(alpha2se) < est_threshold) %>%
    group_by(cen_group) %>%
    summarise(n_obs = n(),
              across(c(alpha0, alpha1, alpha2), mean, .names = "mean_{.col}")) %>%
    mutate(mean_alpha0 = 100*(mean_alpha0 - theta[1])/theta[1],
           mean_alpha1 = 100*(mean_alpha1 - theta[2])/theta[2],
           mean_alpha2 = 100*(mean_alpha2 - theta[3])/theta[3])
          
  return(myresults_first)
}

my_gamma_bias_1 <- my_alpha_bias(-2)
my_gamma_bias_2 <- my_alpha_bias(-1)
my_gamma_bias_3 <- my_alpha_bias(0)
my_gamma_bias_4 <- my_alpha_bias(1)
my_gamma_bias_5 <- my_alpha_bias(2)
my_gamma_bias_6 <- my_alpha_bias(3)
my_gamma_bias_7 <- my_alpha_bias(4)

myresults_gamma_bias = rbind(my_gamma_bias_1, my_gamma_bias_2, my_gamma_bias_3, my_gamma_bias_4, my_gamma_bias_5, my_gamma_bias_6, my_gamma_bias_7)
print(myresults_gamma_bias)

#write.csv(myresults_second, paste0("results-bias-rate_group-",myn,".csv"), row.names = FALSE)
  

p1 = myresults_gamma_bias %>%
  select(cen_group, mean_alpha0, mean_alpha1, mean_alpha2) %>%
  pivot_longer(
    cols = starts_with("mean"),
    names_to = "Parameter",
    values_to = "Bias"
  ) %>%
  # mutate(Parameter = factor(Parameter, labels = c("beta0", "beta1", "beta2"))) %>%
  ggplot(aes(x = cen_group, y = Bias)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
  facet_wrap(~ Parameter) +
  labs(
    title = "Bias by Method and Censoring Group",
    x = "Censoring Group",
    y = "Bias (%)"  #color = "Method"
  ) +
  ylim(-130, 10)+
  theme_minimal() #+ theme(legend.position = "bottom")
p1

# ggsave(paste0("results-bias-rate_group-",myn,".png"), p1, height = 4, width = 8, units = "in")


################################################################################################################################
## C|Y,Z


my_gamma_cov = function(x){
  
  theta = c(x, 0.1, 0.1, 2)
  # theta = c(0.1, 0.1, 0.8)
  #x = -2
  list_est <- list.files(path = paste0("../Simulations-AFT-Both-Est/",myn, "/results"), 
                          pattern = paste0("est_results_alpha_", x, "_*"), full.names = T)
  print(length(list_est))
  myresults <- list_est %>% map_df(read.csv)
  myresults$cen_group = x 
  
  myresults_first = myresults %>%
    bind_rows() %>%
    mutate(across(c(alpha0, alpha1, alpha2, alpha3, alpha0se, alpha1se, alpha2se, alpha3se), as.numeric)) %>%
    subset(abs(alpha0) < est_threshold & abs(alpha1) < est_threshold & abs(alpha2) < est_threshold) %>%
    subset(abs(alpha0se) < est_threshold & abs(alpha1se) < est_threshold & abs(alpha2se) < est_threshold) %>%
    mutate(alpha0_cov = ifelse(alpha0 - qnorm(zval)*alpha0se < theta[1] &
                            alpha0 + qnorm(zval)*alpha0se > theta[1], 1, 0),
         alpha1_cov = ifelse(alpha1 - qnorm(zval)*alpha1se < theta[2] &
                            alpha1 + qnorm(zval)*alpha1se > theta[2], 1, 0),
         alpha2_cov = ifelse(alpha2 - qnorm(zval)*alpha2se < theta[3] &
                            alpha2 + qnorm(zval)*alpha2se > theta[3], 1, 0),
         alpha3_cov = ifelse(log(alpha3) - qnorm(zval)*alpha3se < log(theta[4]) &
                            log(alpha3) + qnorm(zval)*alpha3se > log(theta[4]), 1, 0)) %>%
    select(alpha0:alpha3, alpha0se:alpha3se, alpha0_cov:alpha3_cov, cen_group) %>%
    group_by(cen_group) %>%
    summarise(
      n_obs = n(),
      # Compute mean and SD for beta columns
      across(
        c(alpha0, alpha1, alpha2, alpha3),
        list(mean = ~ mean(.), sd = ~ sd(.)),
        .names = "{.col}_{.fn}"
      ),
      # Compute mean for sebeta columns
      across(
        c(alpha0se:alpha3se, alpha0_cov:alpha3_cov),
        mean,
        .names = "{.col}_mean"
      )
    ) 
          
  return(myresults_first)
}

my_alpha_cov_1 <- my_gamma_cov(-2)
my_alpha_cov_2 <- my_gamma_cov(-1)
my_alpha_cov_3 <- my_gamma_cov(0)
my_alpha_cov_4 <- my_gamma_cov(1)
my_alpha_cov_5 <- my_gamma_cov(2)
my_alpha_cov_6 <- my_gamma_cov(3)
my_alpha_cov_7 <- my_gamma_cov(4)

myresults_alpha_cov = rbind(my_alpha_cov_1, my_alpha_cov_2, my_alpha_cov_3, my_alpha_cov_4, my_alpha_cov_5, my_alpha_cov_6, my_alpha_cov_7)
# head(myresults_gamma_cov)

# make a figure about the coverages 


p2 = myresults_alpha_cov %>%
  select(cen_group, alpha0_cov_mean, alpha1_cov_mean, alpha2_cov_mean, alpha3_cov_mean) %>%
  pivot_longer(
    cols = starts_with("alpha"),
    names_to = "Parameter",
    values_to = "Coverage"
  ) %>%
  #mutate(Parameter = factor(Parameter, labels = c("alpha0", "beta1", "beta2"))) %>%
  ggplot(aes(x = cen_group, y = Coverage)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "red") +
  facet_grid( cols = vars(Parameter)) +
  labs(
    title = "Coverage by Method and Censoring Group",
    x = "Censoring Group",
    y = "Coverage"
    #color = "Method"
  ) +
  ylim(0.5, 1)+
  theme_minimal()
p2 
  
# ggsave(paste0("results-efficiency-rate_group-",myn,".png"), p, height = 4, width = 8, units = "in")
  
my_gamma_bias = function(x){
  
  # x = 0
  theta = c(x, 0.1, 0.1, 2)
  list_est <- list.files(path = paste0("../Simulations-AFT-Both-Est/",myn, "/results"), 
                          pattern = paste0("est_results_alpha_", x, "_*"), full.names = T)
  print(length(list_est))
  myresults <- list_est %>% map_df(read.csv)
  myresults$cen_group = x 
  
  myresults_first = myresults %>%
    bind_rows() %>%
    mutate(across(c(alpha0, alpha1, alpha2, alpha3, alpha0se, alpha1se, alpha2se, alpha3se), as.numeric)) %>%
    subset(abs(alpha0) < est_threshold & abs(alpha1) < est_threshold & abs(alpha2) < est_threshold) %>%
    subset(abs(alpha0se) < est_threshold & abs(alpha1se) < est_threshold & abs(alpha2se) < est_threshold) %>%
    group_by(cen_group) %>%
    summarise(n_obs = n(),
              across(c(alpha0, alpha1, alpha2, alpha3), mean, .names = "mean_{.col}")) %>%
    mutate(mean_alpha0 = 100*(mean_alpha0 - theta[1])/theta[1],
           mean_alpha1 = 100*(mean_alpha1 - theta[2])/theta[2],
           mean_alpha2 = 100*(mean_alpha2 - theta[3])/theta[3],
           mean_alpha3 = 100*(mean_alpha3 - theta[4])/theta[4])
          
  return(myresults_first)
}

my_alpha_bias_1 <- my_gamma_bias(-2)
my_alpha_bias_2 <- my_gamma_bias(-1)
my_alpha_bias_3 <- my_gamma_bias(0)
my_alpha_bias_4 <- my_gamma_bias(1)
my_alpha_bias_5 <- my_gamma_bias(2)
my_alpha_bias_6 <- my_gamma_bias(3)
my_alpha_bias_7 <- my_gamma_bias(4)

myresults_alpha_bias = rbind(my_alpha_bias_1, my_alpha_bias_2, my_alpha_bias_3, my_alpha_bias_4, my_alpha_bias_5, my_alpha_bias_6, my_alpha_bias_7)
print(myresults_alpha_bias)

#write.csv(myresults_second, paste0("results-bias-rate_group-",myn,".csv"), row.names = FALSE)
  

p3 = myresults_alpha_bias %>%
  filter(cen_group != 0) %>%
  select(cen_group, mean_alpha0, mean_alpha1, mean_alpha2, mean_alpha3) %>%
  pivot_longer(
    cols = starts_with("mean"),
    names_to = "Parameter",
    values_to = "Bias"
  ) %>%
  # mutate(Parameter = factor(Parameter, labels = c("beta0", "beta1", "beta2"))) %>%
  ggplot(aes(x = cen_group, y = Bias)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
  facet_grid( cols = vars(Parameter)) +
  labs(
    title = "Bias by Method and Censoring Group",
    x = "Censoring Group",
    y = "Bias (%)"  #color = "Method"
  ) +
  ylim(-130, 10)+
  theme_minimal() #+ theme(legend.position = "bottom")
p3



