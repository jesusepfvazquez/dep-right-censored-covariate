## check-results.R: run from 1000-Main/Table2-66/ while or after the array runs.
library(tidyverse)
SEEDS_PER_TASK <- 50
res <- list.files("results", "^est_results_table2_seed_[0-9]+\\.csv$", full.names = TRUE) %>%
  map_df(read.csv) %>% as_tibble()
cat("Finished seeds:", n_distinct(res$sim), "of 3000; mean censoring",
    round(mean(res$pcen), 3), "\n")
res %>% group_by(Method) %>%
  summarise(reps = n(), failed = sum(is.na(beta0)), se_NA = sum(!is.na(beta0) & is.na(sebeta0)),
            .groups = "drop") %>% print(n = Inf)
todo <- sort(unique(ceiling((setdiff(0:2999, res$sim) + 1) / SEEDS_PER_TASK)))
if (length(todo)) cat("Tasks with unfinished seeds: bash submit.sh", paste(todo, collapse = ","), "\n")
