#!/bin/bash
## Guarded AIPW with the N(0, 1) working model (addon-grd-n01.R).
## One task = 500 seeds (about 1 second per seed); 6 tasks cover seeds 0-2999.
##   sbatch --array=1-6 addon.slurm

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --job-name=T2-66addon
#SBATCH --output=myslurm/%x_%A_%a.out

mkdir -p myslurm results-addon
module load conda_R 2>/dev/null || module load r/4.5.0   # Hopkins (JHPCE) first, then Longleaf
FIRST=$(( 500 * (SLURM_ARRAY_TASK_ID - 1) ))
LAST=$(( 500 * SLURM_ARRAY_TASK_ID - 1 ))
Rscript addon-grd-n01.R $FIRST $LAST
