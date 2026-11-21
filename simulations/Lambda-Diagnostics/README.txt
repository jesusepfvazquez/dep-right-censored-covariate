Lambda-Diagnostics: diagnostics of the guarded AIPW estimator's Lambda and of the weights 1/pi-hat
on the datasets of 1000-Main (Figure 1, 8 censoring rates x 3,000 seeds) and 1000-CensSpread
(20 cells x 1,000 seeds). Replaces the old plane-study diagnostics in Supplement S.6
(Diagnostics-Revision/, 17,862 replicates, before theta = (beta, sigma)).

No estimator is refit. Each dataset is regenerated from its seed, the f_{X|Z} parameters the
guarded fit used are read from the study's results file, and Lambda is rebuilt with the same
code as estimate_beta_aipw_lambda_sigma() (helper_guarded_sigma.R). Tested locally: the rebuilt
Lambda equals the estimator's Lambda to machine precision, and the regenerated datasets match
data_aft_int() seed for seed (same censoring rate).

Place this folder next to 1000-Main and 1000-CensSpread on the cluster (their results/ folders
must be there). Then, from this folder:
  1. bash submit.sh test            10 seeds of each study, run in the login shell (seconds)
                                    check results/test_*.csv: pcen equals pcen_results
  2. bash submit.sh all             8 + 20 array tasks, a few minutes each
     bash submit.sh extra           12 tasks for the added 1000-CensSpread cells (scale ratios 2.5, 3.5, 5)
  3. Rscript lambda-diag-summary.R  tables/lambda-diag-main.{csv,tex}, tables/lambda-diag-cs.{csv,tex}
                                    the columns no_results and mismatch must be 0 everywhere

Columns: ESS = median ESS / n uncensored; condition number = kappa of M = sum V V^T, the matrix
inverted to form Lambda (median, 95th percentile); singular = rcond(M) < 1e-10; max weight =
median largest weight; failed / diverged as in cens-spread.R (no finite estimate / an estimate
more than 10 from the truth). The CensSpread failed / diverged counts should reproduce
cens-spread-summary.csv (unguarded 434 / 16 at scale ratio 1.5 and 85%).
