Fig1-TwoStep: two-step MLE add-on for Figure 1 (same datasets: seeds 0-2999, eight
censoring intercepts, n = 1000). Code for the fit and SE is the Table 2 code
(helper_table2.R), so the c_int = -0.5 rows match Table2-66's mle-twostep rows.

Run (from this folder, on the cluster, with Table2-66 next to it in 1000-Main)
  1. bash submit.sh test          task 1 = seeds 0-49 x 8 levels; check myslurm/*.out
  2. Rscript compare-table2.R     max diffs should be ~0 (Table2-66 results must be present)
  3. bash submit.sh all           tasks 2-60
  4. Rscript check-results.R      rerun listed tasks with: bash submit.sh <list>
  5. cd ../Figures-Tables; Rscript fig1-bias-coverage.R   (reads ../Fig1-TwoStep/results)

ESS shading for Figure 1 (optional, separate from the two-step run)
  6. sbatch ess.slurm             writes results-ess/ess_0_2999.csv (about 15 minutes)
  7. cd ../Figures-Tables; Rscript fig1-bias-coverage.R   (shades the median ESS in the coverage row)
