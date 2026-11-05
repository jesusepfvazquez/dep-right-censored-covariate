###############################################################################
## summarize-lambda-diagnostics.R
## Aggregate the per-replicate Lambda diagnostics into per-censoring-bin
## summaries (matching the 10-point bins in Figures 2-4) and print a compact
## table for the response to Reviewer 2.
###############################################################################
suppressMessages(library(dplyr))

inp <- "Diagnostics-Revision/lambda-diagnostics.csv"
d   <- read.csv(inp)

## 10-percentage-point censoring bins, as in the heatmaps
d$cen_bin <- cut(d$pcen, breaks = seq(0, 1, by = 0.1), right = FALSE,
                 include.lowest = TRUE,
                 labels = c("[0,10)","[10,20)","[20,30)","[30,40)","[40,50)",
                            "[50,60)","[60,70)","[70,80)","[80,90)","[90,100]"))

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

write.csv(summ, "Diagnostics-Revision/lambda-diagnostics-summary.csv",
          row.names = FALSE)

cat("\n=== Lambda-hat diagnostics by censoring bin ===\n")
print(as.data.frame(summ), digits = 3, row.names = FALSE)

cat("\nOverall singular-fit rate: ",
    round(100 * mean(d$singular, na.rm = TRUE), 2), "%\n", sep = "")
cat("Rows: ", nrow(d), "   Bins populated: ", nrow(summ), "\n", sep = "")
