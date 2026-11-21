1000-CensSpread: does a wider censoring distribution help the weighted estimators?

Design: the Figure 1 design (data_aft_int() from 1000-Main, copied in helper_data_cs.R
without the imputation integrals; same draws), n = 1,000.
  log X = 0.1 + 0.1 Z + 0.5 log Exp(1);  Y = 1 + (A - X) + Z + e
  log C = c_int + 0.1 Y + 0.1 Z + p.c log Exp(1)
Grid: p.c (scale of log C) in {0.75, 1, 1.5, 2, 3} by censoring rate in {50, 65, 75, 85}%,
20 cells; c_int calibrated per cell on 1e6 draws (calibrate.R; values hard-coded in
main.R). p.c = 1.5 is the Figure 1 design. 1,000 seeds per cell.

p.c = 0.5 is left out on purpose: log X has scale 0.5, and with p.c <= 0.5 (and C's
scale below X's) 1/pi(x) grows faster than f_{X|Z}(x) decays, so E[1/pi(X)] is infinite;
the unguarded fit failed in 5 of 6 test datasets. Even at p.c = 0.75 and 85% censoring the
weights are extreme (effective sample size about 5% of the uncensored n) and the
unguarded fit fails in most datasets (also from the joint MLE start); report the failure
rate for that cell (check-results.R).

Nov 14 fix (helper_joint_mle.R): the joint MLE now starts gamma from the IPW-weighted Weibull
fit of W on Z among uncensored rows. The old start (survreg of W on Z, ignoring that C
depends on Y) did not converge or converged far from the truth at p.c = 3; the bounded
refit then stopped on the edge of its box (MLE beta0 bias -20 to -30%, coverage 7%), and
the AIPW rows inherit gamma from the MLE. Results from the first version must be deleted.
New columns: mle_edge, mle_restarted. Estimates more than 10 from the truth (the guarded
fit at p.c <= 1 with heavy censoring) are counted as "diverged" and left out of summaries.

Estimators: oracle, IPW (true pi), IPW (estimated pi, corrected SE), unguarded AIPW
(f_{X|Z} from the joint MLE), guarded AIPW (working normal from the joint MLE), joint MLE
(does not use pi; reference). Per dataset, weight diagnostics with the estimated pi:
largest weight, its share of the total weight, effective sample size / uncensored n.

Run (from this folder; nothing to edit)
  1. bash submit.sh test     task 1 (p.c = 0.75, 50%, seeds 1-50); check myslurm/*.out
  2. bash submit.sh first    tasks 2-20: 50 seeds in every cell (first look)
  3. bash submit.sh rest     tasks 21-400: 1,000 seeds per cell
  4. Rscript check-results.R lists cells, failures and tasks to resubmit
  5. cd Figures-Tables; Rscript cens-spread.R
       figures/cs-lines-<rate>.png   Figure 1 layout with p.c on the x axis, per censoring rate
       figures/cs-heat-coverage.png  coverage heat map (p.c by censoring rate)
       figures/cs-heat-bias.png      percent bias heat map
       figures/cs-heat-weights.png   weight diagnostics heat map
       tables/cens-spread-summary.csv, tables/cens-spread-weights.csv
R packages: survival, optimx, numDeriv; tidyverse for the summary scripts.

Added configurations (Nov 2026): scale ratios 2.5, 3.5 and 5 (p.c = 1.25, 1.75, 2.5) at the
same four censoring rates, 12 cells x 1,000 seeds, which fill the gap between median ESS
67% and 80% in the original grid. c_int from calibrate.R, hard-coded in main.R under
GRID_SET=extra. Same results/ folder and file names, so cens-spread.R adds them as columns.
  bash submit.sh extra-test     task 1 of the added grid
  bash submit.sh extra          tasks 2-240
  Rscript check-results.R       checks both grids and prints resubmit commands
