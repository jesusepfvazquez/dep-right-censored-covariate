## compare-fig1.R: check that Table2-66 reproduces the 66% point of Figure 1.
## Run from 1000-Main/Table2-66/ after the test task. Needs Figure 1's results in ../results
## (files est_results_-0.5_seed_<s>.csv). Same seed -> same dataset -> the shared rows
## should agree (oracle exactly; IPW, joint MLE and both AIPW to numerical tolerance).
library(tidyverse)
t2 <- list.files("results", "^est_results_table2_seed_[0-9]+\\.csv$", full.names = TRUE) %>% map_df(read.csv)
seeds <- unique(t2$sim)
f1 <- map_df(seeds, function(s) {
  f <- paste0("../results/est_results_-0.5_seed_", s, ".csv")
  if (file.exists(f)) read.csv(f) %>% mutate(sim = s) else NULL })
if (!nrow(f1)) stop("no Figure 1 files for these seeds in ../results")
map <- c(oracle = "oracle", "ipw-est" = "ipw-est", "mle-joint" = "mle-joint",
         "ung-est-est" = "aipw-est-joint", "grd-est" = "aipw-lambda-est", cc = "cc", naive = "naive")
cmp <- t2 %>% filter(Method %in% names(map)) %>% mutate(F1 = map[Method]) %>%
  inner_join(f1, by = c("sim", "F1" = "Method"), suffix = c(".t2", ".f1")) %>%
  group_by(Method) %>%
  summarise(seeds = n(), censoring_t2 = mean(pcen.t2), censoring_f1 = mean(pcen.f1),
            max_diff_beta = max(abs(c(beta0.t2 - beta0.f1, beta1.t2 - beta1.f1, beta2.t2 - beta2.f1)), na.rm = TRUE),
            max_diff_se   = max(abs(c(sebeta0.t2 - sebeta0.f1, sebeta1.t2 - sebeta1.f1)), na.rm = TRUE),
            .groups = "drop")
print(cmp, width = Inf)
