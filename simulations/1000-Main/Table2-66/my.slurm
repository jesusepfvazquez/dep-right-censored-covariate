#!/bin/bash
## One task = 50 seeds (about 2-4 minutes per seed, so 2-3.5 hours). 60 tasks give the 3,000 replicates.

#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=10:00:00
#SBATCH --job-name=Table2-66
#SBATCH --output=myslurm/%x_%A_%a.out

echo "SLURM_ARRAY_JOB_ID: $SLURM_ARRAY_JOB_ID  TASK: $SLURM_ARRAY_TASK_ID"
mkdir -p myslurm results
module load conda_R 2>/dev/null || module load r/4.5.0   # Hopkins (JHPCE) first, then Longleaf
Rscript main.R $SLURM_ARRAY_TASK_ID
