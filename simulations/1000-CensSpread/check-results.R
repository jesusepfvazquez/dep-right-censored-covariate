## check-results.R: run from 1000-CensSpread/ while or after the array runs.
library(tidyverse)
res <- list.files("results", "^est_results_pc.*\\.csv$", full.names = TRUE) %>% map_df(read.csv) %>% as_tibble()
res %>% distinct(p.c, target, sim, pcen) %>% group_by(p.c, target) %>%
  summarise(seeds = n(), pcen = round(mean(pcen), 3), .groups = "drop") %>% print(n = Inf)
res %>% group_by(Method) %>%
  summarise(fits = n(), failed = sum(is.na(beta0)), se_NA = sum(!is.na(beta0) & is.na(sebeta0)),
            ung_retry = sum(ung_retry, na.rm = TRUE),
            diverged = sum(abs(beta0 - 1) > 10 | abs(beta1 - 1) > 10 | abs(beta2 - 1) > 10, na.rm = TRUE),
            mle_edge = sum(mle_edge, na.rm = TRUE), .groups = "drop") %>% print()
## tasks with unfinished seeds (1,000 per cell)
done <- res %>% distinct(p.c, target, sim)
for (set in c("original", "extra")) {
  grid <- if (set == "original") expand.grid(p.c = c(0.75, 1, 1.5, 2, 3), target = c(0.50, 0.65, 0.75, 0.85)) else
                                 expand.grid(p.c = c(1.25, 1.75, 2.5), target = c(0.50, 0.65, 0.75, 0.85))
  todo <- unlist(lapply(seq_len(nrow(grid)), function(k) {
    miss <- setdiff(1:1000, done$sim[abs(done$p.c - grid$p.c[k]) < 1e-9 & abs(done$target - grid$target[k]) < 1e-9])
    unique((ceiling(miss / 50) - 1) * nrow(grid) + k) }))
  if (length(todo)) cat(set, "grid, tasks with unfinished seeds:", length(todo), "\n  bash submit.sh",
                        if (set == "extra") paste0("extra:", paste(sort(todo), collapse = ",")) else paste(sort(todo), collapse = ","), "\n")
}
