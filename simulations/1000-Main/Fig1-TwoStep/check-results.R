## check-results.R: failures and unfinished tasks
library(tidyverse)
SEEDS_PER_TASK <- 50
res <- list.files("results", "^est_results_", full.names = TRUE) %>% map_df(read.csv)
res %>% group_by(group) %>%
  summarise(reps = n(), failed = sum(is.na(beta0)), se_NA = sum(!is.na(beta0) & is.na(sebeta0)),
            pcen = round(100 * mean(pcen), 1), .groups = "drop") %>% print(n = Inf)
done <- res %>% count(sim) %>% filter(n == 8) %>% pull(sim)
todo <- sort(unique(ceiling((setdiff(0:2999, done) + 1) / SEEDS_PER_TASK)))
if (length(todo)) cat("Tasks with unfinished seeds: bash submit.sh", paste(todo, collapse = ","), "\n")
