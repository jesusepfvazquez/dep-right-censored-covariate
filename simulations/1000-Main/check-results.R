## check-results.R: run from 1000-Main/ while or after the array runs.
## Counts finished files, failures and guard flags, and lists seeds to rerun.
library(tidyverse)
seeds  <- 0:2999
c_ints <- c(4.5, 4, 1.5, 1, 0.5, -0.5, -1, -1.5)

res <- list.files("results", "^est_results_-?[0-9.]+_seed_[0-9]+\\.csv$", full.names = TRUE) %>%
  map_df(read.csv) %>% as_tibble()

cat("\nFinished files per censoring intercept (target:", length(seeds), ")\n")
res %>% distinct(sim, group, pcen) %>% group_by(group) %>%
  summarise(files = n(), pcen = mean(pcen), .groups = "drop") %>% print(n = Inf)

cat("\nFailures and flags by estimator\n")
res %>%
  group_by(Method) %>%
  summarise(reps = n(), failed = sum(is.na(beta0)),
            se_NA = sum(!is.na(beta0) & is.na(sebeta0)),
            conv_not0 = sum(!is.na(convcode) & convcode != 0),
            guard = sum(coalesce(guard, 0) > 0),
            start_fixed = sum(coalesce(start_fixed, 0) > 0),
            .groups = "drop") %>%
  print(n = Inf)

cat("\nJoint MLE f_{X|Z} estimates (truth: gamma0 = 0.1, gammaZ = 0.1, scale_x = 0.5)\n")
res %>% filter(Method == "mle-joint", !is.na(gamma0)) %>%
  group_by(group) %>%
  summarise(reps = n(), gamma0 = mean(gamma0), gamma0_med = median(gamma0),
            gammaZ = mean(gammaZ), scale_x = mean(scale_x), .groups = "drop") %>%
  print(n = Inf)

## seeds to rerun: missing any censoring level, or the joint MLE failed / has no SE
done <- res %>% distinct(sim, group) %>% count(sim) %>% filter(n == length(c_ints)) %>% pull(sim)
bad  <- res %>% filter(Method == "mle-joint", is.na(beta0) | is.na(sebeta0)) %>% pull(sim)
redo <- sort(unique(c(setdiff(seeds, done), bad)))
cat("\nSeeds not finished or with a failed joint MLE:", length(redo), "\n")
if (length(redo)) cat("bash submit.sh", paste(redo, collapse = ","), "\n")
cat("(Seeds still running also show up here; wait for the array to finish first.)\n")
