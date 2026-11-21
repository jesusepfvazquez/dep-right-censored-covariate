#!/bin/bash
## One task = 50 seeds x 8 censoring intercepts. 60 tasks cover seeds 0-2999.

#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=10:00:00
#SBATCH --job-name=Fig1-TwoStep
#SBATCH --output=myslurm/%x_%A_%a.out

echo "SLURM_ARRAY_JOB_ID: $SLURM_ARRAY_JOB_ID  TASK: $SLURM_ARRAY_TASK_ID"
mkdir -p myslurm results
module load conda_R 2>/dev/null || module load r/4.5.0   # Hopkins (JHPCE) first, then Longleaf
Rscript twostep.R $SLURM_ARRAY_TASK_ID
