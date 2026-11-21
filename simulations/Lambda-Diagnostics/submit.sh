#!/bin/bash
## Usage (run from this folder; 1000-Main and 1000-CensSpread, with their results/, sit next to it):
##   bash submit.sh test   # 10 seeds of 1000-Main at 66% and of the hardest CensSpread cell, run here
##   bash submit.sh all    # 8 tasks for 1000-Main + 20 tasks for 1000-CensSpread
##   bash submit.sh extra  # 12 tasks for the extra CensSpread cells (after those cells finish)
##   then: Rscript lambda-diag-summary.R
mkdir -p results myslurm
case "$1" in
  test)
    module load conda_R 2>/dev/null || module load r/4.5.0
    SEEDS=0:9 Rscript lambda-diag.R main 6 && mv results/diag_main_-0.5.csv results/test_main.csv
    SEEDS=1:10 Rscript lambda-diag.R cs 16 && mv results/diag_cs_pc0.75_t0.85.csv results/test_cs.csv
    echo "Look at results/test_main.csv and results/test_cs.csv (pcen should equal pcen_results)." ;;
  all)
    sbatch --export=ALL,STUDY=main --array=1-8  my.slurm
    sbatch --export=ALL,STUDY=cs   --array=1-20 my.slurm ;;
  extra)
    sbatch --export=ALL,STUDY=cs-extra --array=1-12 my.slurm ;;
  *)
    echo "bash submit.sh test | all | extra" ;;
esac
