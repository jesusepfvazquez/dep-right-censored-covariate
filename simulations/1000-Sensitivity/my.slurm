#!/bin/bash
## One task = one seed (all six xi values). Submit with submit.sh, which sets OFFSET
## so that seed = OFFSET + SLURM_ARRAY_TASK_ID stays under Slurm's MaxArraySize.

#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=Sen
#SBATCH --output=myslurm/%x_%A_%a.out

SEED=$(( ${OFFSET:-0} + SLURM_ARRAY_TASK_ID ))
echo "SLURM_ARRAY_JOB_ID: $SLURM_ARRAY_JOB_ID  TASK: $SLURM_ARRAY_TASK_ID  SEED: $SEED"

mkdir -p myslurm results
module load r/4.5.0        # Hopkins (JHPCE): module load conda_R
Rscript main.R $SEED
