#!/bin/bash
## Usage (run from this folder):
##   bash submit.sh test    # task 1 (p.c = 0.75, 50% censoring, seeds 1-50)
##   bash submit.sh first   # tasks 2-20: 50 seeds in every cell (a quick first look)
##   bash submit.sh rest    # tasks 21-400: up to 1,000 seeds per cell
##   bash submit.sh 4,17    # rerun specific tasks; finished seeds are skipped
## Extra scale ratios 2.5, 3.5, 5 (Nov 21; 12 cells, 240 tasks):
##   bash submit.sh extra-test   # task 1 of the extra grid (p.c = 1.25, 50%, seeds 1-50)
##   bash submit.sh extra        # tasks 2-240: 1,000 seeds in each extra cell
##   bash submit.sh extra:4,17   # rerun specific extra tasks
mkdir -p results myslurm
case "$1" in
  test)  sbatch --array=1 my.slurm ;;
  first) sbatch --array=2-20 my.slurm ;;
  rest)  sbatch --array=21-400 my.slurm ;;
  extra-test) sbatch --export=ALL,GRID_SET=extra --array=1 my.slurm ;;
  extra)      sbatch --export=ALL,GRID_SET=extra --array=2-240 my.slurm ;;
  extra:*)    sbatch --export=ALL,GRID_SET=extra --array="${1#extra:}" my.slurm ;;
  *)     sbatch --array="$1" my.slurm ;;
esac
