## 1000-Main/Fig1-TwoStep/ess.R
## Effective sample size of the weights 1/pi-hat (IPW and AIPW) for every Figure 1 dataset (seeds 0-2999, eight
## censoring intercepts), the same definition as in 1000-CensSpread:
##   w = 1 / pi-hat(W, Y, Z) on the uncensored rows, pi-hat from survreg(Surv(W, 1 - D) ~ y + Z)
##   ess_frac = (sum w)^2 / sum(w^2) / (number uncensored)   (1 = all weights equal)
## Only regenerates the data (no estimators), about 0.03 s per dataset.
## Output: results-ess/ess_<first>_<last>.csv. Usage: Rscript ess.R <first seed> <last seed>
library(survival)
source("helper_data_fig1.R")
a <- as.integer(commandArgs(TRUE)); first <- a[1]; last <- a[2]
c_ints <- c(4.5, 4, 1.5, 1, 0.5, -0.5, -1, -1.5)
dir.create("results-ess", showWarnings = FALSE)
out <- do.call(rbind, lapply(first:last, function(s) do.call(rbind, lapply(c_ints, function(c_int) {
  set.seed(s)
  dat <- data_aft(nSubjects = 1000, c_int = c_int, p.c = 1.5)
  w   <- 1 / dat$myp_ywz[dat$D == 1]
  data.frame(sim = s, group = c_int, pcen = 1 - mean(dat$D),
             ess_frac = sum(w)^2 / sum(w^2) / length(w), max_w_share = max(w) / sum(w))
}))))
write.csv(out, sprintf("results-ess/ess_%d_%d.csv", first, last), row.names = FALSE)
