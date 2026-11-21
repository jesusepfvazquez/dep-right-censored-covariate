#!/bin/bash
## One task = one censoring intercept of 1000-Main (3,000 seeds) or one 1000-CensSpread cell
## (1,000 seeds). No estimator is refit, so a task takes a few minutes. STUDY is set by submit.sh.

#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=03:00:00
#SBATCH --job-name=LambdaDiag
#SBATCH --output=myslurm/%x_%A_%a.out

echo "STUDY: $STUDY  TASK: $SLURM_ARRAY_TASK_ID"
mkdir -p myslurm results
module load conda_R 2>/dev/null || module load r/4.5.0   # Hopkins (JHPCE) first, then Longleaf
Rscript lambda-diag.R $STUDY $SLURM_ARRAY_TASK_ID
