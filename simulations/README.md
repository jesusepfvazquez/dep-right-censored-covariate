# Simulation code

This folder has the code for the two simulation studies and the sensitivity analysis in the
manuscript, and for the diagnostics in the Supplementary Material. Each study folder is self-contained: it sources
only its own helper files and does not need the `ODRCC` package. The per-replicate results
are not included because of their size. Every result is regenerated from its seed.

| Folder | Produces |
|---|---|
| `1000-Main/` | First simulation study (fixed n = 1,000; eight right-censoring rates from about 4% to 85%; 3,000 replicates each): Figure 1 |
| `1000-Main/Fig1-TwoStep/` | The two-step MLE curve and the median ESS shading in Figure 1 (same datasets as `1000-Main/`) |
| `1000-Main/Table2-66/` | Table 2 (the 66% right-censoring point of Figure 1, same 3,000 datasets; correct, estimated, and incorrect specifications) |
| `1000-CensSpread/` | Second simulation study (scale ratio by right-censoring rate; 1,000 replicates per configuration): the heat maps of median ESS, percent bias, coverage, and MSE ratio |
| `1000-Sensitivity/` | Sensitivity analysis for the assumption that $C$ is independent of $X$ given $(Y, Z)$ (n = 800, about 78% right-censoring, 3,001 replicates per value of the sensitivity parameter): the sensitivity figure |
| `Lambda-Diagnostics/` | The diagnostics of the weights and of the guarded AIPW matrix in the Supplementary Material, for both studies |

## Requirements

R (version 4.4 or later) with the packages `survival`, `optimx`, `numDeriv`, `tidyverse`,
`scales`, and `ggh4x`. `1000-Main/helper_mestimate.R` also loads `geex`, `cubature`, and
`tictoc`.

## Running

The studies were run on a SLURM cluster. Each folder has a `README.txt` with the exact
steps, a `my.slurm` job script (one array task = a block of seeds), a `submit.sh` wrapper
(`bash submit.sh test` runs one task first), and a `check-results.R` script that lists
any tasks to resubmit. The figure and table scripts are in each folder's
`Figures-Tables/` and read that folder's `results/`.

Order:

1. `1000-Main/`, then `1000-Main/Fig1-TwoStep/` and `1000-Main/Table2-66/`.
2. `1000-CensSpread/`: the original 20 configurations (`bash submit.sh test`, `first`, `rest`) and the
   12 added configurations at scale ratios 2.5, 3.5, and 5 (`bash submit.sh extra-test`, `extra`).
3. `1000-Sensitivity/` (independent of the other studies).
4. `Lambda-Diagnostics/` (after both studies finish, with the study folders next to it):
   `bash submit.sh all`, `bash submit.sh extra`, then `Rscript lambda-diag-summary.R`.

Every replicate is `set.seed(seed)` followed by the study's data generator, so any single
replicate can be regenerated and refit on a laptop: `Rscript main.R 0` in `1000-Main/` fits
seed 0 at all eight right-censoring rates (and in `1000-Sensitivity/`, at all six values of the
sensitivity parameter), and `SEEDS_PER_TASK=1 Rscript main.R 1` does the
same for one seed of task 1 in `1000-CensSpread/` and `1000-Main/Table2-66/`.
