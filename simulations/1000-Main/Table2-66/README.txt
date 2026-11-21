1000-Main/Table2-66: Table 2 (tab:chp3-sims_1000) on the Figure 1 design, at the 66% point.

Design: the Figure 1 generator, data_aft_int(1000, c_int = -0.5, p.c = 1.5), copied into
helper_data_fig1.R as data_aft() without the imputation integrals (same draws, same order).
Seeds 0-2999, exactly as in Figure 1, so each replicate is the same dataset as Figure 1's
66% point, and the estimated, correctly specified rows reproduce that point.
  log X = 0.1 + 0.1 Z + 0.5 log Exp(1)           gamma = (0.1, 0.1), scale_x = 0.5
  Y = 1 + (A - X) + Z + e, e ~ N(0, 1)
  log C = -0.5 + 0.1 Y + 0.1 Z + 1.5 log Exp(1)   alpha = (-0.5, 0.1, 0.1), scale 1.5
About 66% right-censoring. N = 3,000 replicates.

Changes from 1000-Main/Table2 (which used ODRCC::data_aft(), about 60% censoring):
  helper_data_fig1.R  replaces helper_data_aft.R
  helper_known.R      alpha_true = (-0.5, 0.1, 0.1, 1.5), gamma_true = (0.1, 0.1, 0.5)
  helper_table2.R     gamma_wrong = (0.1, 0.1, 1.2): incorrect Weibull scale 1.2 instead of 0.5
  main.R              seeds 0-2999; ung-est-est uses Figure 1's grid (upper 20, m 200);
                      cmi-z-odrcc dropped (it only reproduced an ODRCC bug in data_aft())
  addon-grd-n01.R, addon.slurm, submit.sh, check-results.R: seeds 0-2999
  compare-fig1.R      new: checks the shared rows against Figure 1's results

Rows: the same as 1000-Main/Table2 (see main.R): oracle; IPW correct / estimated / incorrect
pi; unguarded AIPW under five specifications (plus the two -n01 rows); guarded AIPW
correct / estimated / incorrect pi; MLE correct / two-step / joint / incorrect f_{X|Z};
CC, ACC, MACC, CMI (E[X | Y, Z, X > C] and E[X | Z]), naive.

Run (from this folder; R packages survival, optimx, numDeriv, tidyverse)
  1. bash submit.sh test          task 1 = seeds 0-49
  2. Rscript compare-fig1.R       (if Figure 1's results are in ../results) shared rows
                                  should match Figure 1 to numerical tolerance
  3. bash submit.sh all           tasks 2-60
  4. sbatch --array=1-6 addon.slurm   optional: guarded AIPW with the N(0, 1) Psi
  5. Rscript check-results.R      lists tasks with unfinished seeds
  6. cd Figures-Tables; Rscript table2.R
       tables/table2-summary.csv and tables/table2-rows.tex
