1000-Main: first simulation study (Figure 1). n = 1000, eight censoring intercepts
(about 4% to 85% right-censored), seeds 0-2999. Each dataset is set.seed(seed)
followed by data_aft_int(n, c_int, 1.5). All estimators are fit (helper_fit_all.R).

Run
  1. bash submit.sh test         check myslurm/*.out and results/ (8 files)
  2. bash submit.sh all          seeds 0-2999
  3. Rscript check-results.R     rerun listed seeds with: bash submit.sh <list>
  4. cd Figures-Tables; Rscript fig1-bias-coverage.R
     (figures/02-bias-coverage-estimators.png, tables/fig1-bias-coverage.csv)

Not produced here: Table 2 (60% right-censoring, correct/estimated/incorrect specifications).
