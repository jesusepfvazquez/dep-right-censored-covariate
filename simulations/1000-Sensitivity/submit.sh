#!/bin/bash
## Usage:
##   bash submit.sh test         # seed 0 only; check run time and errors first
##   bash submit.sh all          # seeds 0-3000, in chunks of 1000
##   bash submit.sh 24,246,1500  # rerun specific seeds (listed by check-results.R)
mkdir -p results myslurm

case "$1" in
  test)
    sbatch --export=ALL,OFFSET=0 --array=0 my.slurm ;;
  all)
    sbatch --export=ALL,OFFSET=0    --array=0-999 my.slurm
    sbatch --export=ALL,OFFSET=1000 --array=0-999 my.slurm
    sbatch --export=ALL,OFFSET=2000 --array=0-1000 my.slurm ;;
  *)
    sbatch --export=ALL,OFFSET=0 --array="$1" my.slurm ;;
esac
