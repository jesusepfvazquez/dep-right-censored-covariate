1000-Sensitivity: sensitivity analysis for C independent of X given (Y, Z)
(Section 5.3 and Figure "sensitivity.png" in the paper). n = 800, about 78% censoring,
xi = 0, 0.2, ..., 1, seeds 0-3000. The IPW and MLE sandwiches are off (the figure
reports bias); set do_var <- TRUE in main.R to compute them. The AIPW SEs are always computed.

Run
  1. bash submit.sh test        then check myslurm/*.out and results/ (6 files)
  2. bash submit.sh all         seeds 0-3000
  3. Rscript check-results.R    rerun any listed seeds with: bash submit.sh <list>
  4. cd Figures-Tables; Rscript sensitivity.R   (figures/sensitivity.png, tables/sensitivity.csv)
