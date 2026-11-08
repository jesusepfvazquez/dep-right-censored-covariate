###############################################################################
## summarize-lambda-diagnostics.R
## Aggregate the per-replicate Lambda diagnostics into per-censoring-bin
## summaries (matching the 10-point bins in Figures 2-4) and print a compact
## table for the response to Reviewer 2.
###############################################################################
suppressMessages(library(dplyr))

## Censoring bins: left-closed, right-open, [a, b); the last bin is [90, 100].
## Round first so that floating-point error does not move boundary replicates
## (e.g., in R, 0.7 < 7 * 0.1). No reported replicate falls at exactly 90%.
cens_bin10 <- function(p) {
  pct <- round(100 * p, 6)
  lab <- c("[0, 10)", "[10, 20)", "[20, 30)", "[30, 40)", "[40, 50)",
           "[50, 60)", "[60, 70)", "[70, 80)", "[80, 90)", "[90, 100]")
  idx <- pmin(floor(pct / 10) + 1, 10)
  factor(lab[idx], levels = lab)
}


inp <- "data/lambda-diagnostics.csv"
d   <- read.csv(inp)

## 10-percentage-point censoring bins, as in the heatmaps
d$cen_bin <- cens_bin10(d$pcen)

q95 <- function(x) quantile(x, 0.95, na.rm = TRUE)

summ <- d %>%
  group_by(cen_bin) %>%
  summarise(
    nrep            = n(),
    singular_pct    = 100 * mean(singular, na.rm = TRUE),
    kappa_L_median  = median(kappa_lambda[is.finite(kappa_lambda)], na.rm = TRUE),
    kappa_L_p95     = q95(kappa_lambda[is.finite(kappa_lambda)]),
    kappa_part2_med = median(kappa_part2[is.finite(kappa_part2)], na.rm = TRUE),
    max_wt_median   = median(max_wt, na.rm = TRUE),
    max_wt_p95      = q95(max_wt),
    ess_frac_median = median(ess_frac, na.rm = TRUE),
    frac_extreme_md = median(frac_extreme, na.rm = TRUE),
    .groups = "drop")

write.csv(summ, "data/lambda-diagnostics-summary.csv",
          row.names = FALSE)

cat("\n=== Lambda-hat diagnostics by censoring bin ===\n")
print(as.data.frame(summ), digits = 3, row.names = FALSE)

cat("\nOverall singular-fit rate: ",
    round(100 * mean(d$singular, na.rm = TRUE), 2), "%\n", sep = "")
cat("Rows: ", nrow(d), "   Bins populated: ", nrow(summ), "\n", sep = "")
