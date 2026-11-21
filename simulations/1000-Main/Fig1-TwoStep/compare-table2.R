## compare-table2.R: the c_int = -0.5 rows here should equal the mle-twostep rows of
## ../Table2-66 (same datasets, same code). Run after the test task.
library(tidyverse)
a <- list.files("results", "^est_results_-0.5_seed_", full.names = TRUE) %>% map_df(read.csv)
b <- list.files("../Table2-66/results", "^est_results_table2_seed_", full.names = TRUE) %>%
  map_df(read.csv) %>% filter(Method == "mle-twostep", sim %in% a$sim)
m <- inner_join(a, b, by = "sim", suffix = c("", ".t2"))
cat("seeds compared:", nrow(m), "\n")
print(summarise(m, max_diff_beta = max(abs(c(beta0 - beta0.t2, beta1 - beta1.t2, beta2 - beta2.t2))),
                max_diff_se = max(abs(c(sebeta0 - sebeta0.t2, sebeta1 - sebeta1.t2)), na.rm = TRUE),
                max_diff_pcen = max(abs(pcen - pcen.t2))))
