#!/bin/bash
## One task = one (p.c, censoring) cell x 50 seeds (about 1-2 minutes per seed).
## 400 tasks give 1,000 seeds in each of the 20 cells; tasks 1-20 cover every cell once.

#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=06:00:00
#SBATCH --job-name=CensSpread
#SBATCH --output=myslurm/%x_%A_%a.out

echo "SLURM_ARRAY_JOB_ID: $SLURM_ARRAY_JOB_ID  TASK: $SLURM_ARRAY_TASK_ID"
mkdir -p myslurm results
module load conda_R 2>/dev/null || module load r/4.5.0   # Hopkins (JHPCE) first, then Longleaf
Rscript main.R $SLURM_ARRAY_TASK_ID
