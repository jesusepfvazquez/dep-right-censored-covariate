#!/bin/bash
## Usage (run from this folder):
##   bash submit.sh test    # task 1 (seeds 0-49); check myslurm/*.out and the run time first
##   bash submit.sh all     # tasks 2-60 (seeds 50-2999)
##   bash submit.sh 4,17    # rerun specific tasks; finished files are skipped
mkdir -p results myslurm
case "$1" in
  test) sbatch --array=1 my.slurm ;;
  all)  sbatch --array=2-60 my.slurm ;;
  *)    sbatch --array="$1" my.slurm ;;
esac
