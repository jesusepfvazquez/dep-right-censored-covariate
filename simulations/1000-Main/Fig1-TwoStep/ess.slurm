#!/bin/bash
## ESS of the weights 1/pi-hat for all Figure 1 datasets: one task, seeds 0-2999 (about 15 minutes).
##   sbatch ess.slurm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --job-name=Fig1-ESS
#SBATCH --output=myslurm/%x_%A.out
mkdir -p myslurm results-ess
module load conda_R 2>/dev/null || module load r/4.5.0   # Hopkins (JHPCE) first, then Longleaf
Rscript ess.R 0 2999
