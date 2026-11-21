## check-results.R: run from 1000-Sensitivity/ while or after the array runs.
## Counts finished files, failures and guard flags, and lists seeds to rerun.
library(tidyverse)
seeds   <- 0:3000
xi_grid <- seq(0, 1, by = 0.2)

res <- list.files("results", "^est_results_xi[0-9.]+_seed_[0-9]+\\.csv$", full.names = TRUE) %>%
  map_df(read.csv) %>% as_tibble()

cat("\nFinished (seed, xi) files per xi\n")
res %>% distinct(sim, xi) %>% count(xi) %>% print(n = Inf)

cat("\nFailures and flags by xi and estimator\n")
res %>%
  group_by(xi, Method) %>%
  summarise(reps = n(), failed = sum(is.na(beta0)),
            se_NA = sum(!is.na(beta0) & is.na(sebeta0)),
            conv_not0 = sum(!is.na(convcode) & convcode != 0),
            guard = sum(coalesce(guard, 0) > 0),
            start_fixed = sum(coalesce(start_fixed, 0) > 0),
            .groups = "drop") %>%
  print(n = Inf)

cat("\nJoint MLE f_{X|Z} estimates (truth: gamma0 = 0.1, gammaZ = 0.1, scale_x = 0.5)\n")
res %>% filter(Method == "mle-joint", !is.na(gamma0)) %>%
  group_by(xi) %>%
  summarise(reps = n(), gamma0 = mean(gamma0), gamma0_med = median(gamma0),
            gammaZ = mean(gammaZ), scale_x = mean(scale_x), .groups = "drop") %>%
  print(n = Inf)

## seeds to rerun: missing any xi, or the joint MLE failed / has no SE
done <- res %>% distinct(sim, xi) %>% count(sim) %>% filter(n == length(xi_grid)) %>% pull(sim)
bad  <- res %>% filter(Method == "mle-joint", is.na(beta0)) %>% pull(sim)
redo <- sort(unique(c(setdiff(seeds, done), bad)))
cat("\nSeeds not finished or with a failed joint MLE:", length(redo), "\n")
if (length(redo)) cat("bash submit.sh", paste(redo, collapse = ","), "\n")
cat("(Seeds still running also show up here; wait for the array to finish first.)\n")
