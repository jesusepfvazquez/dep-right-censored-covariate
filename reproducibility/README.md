# Reproducibility and diagnostic report

This folder documents the computational record behind the simulation study in
the manuscript. It covers three things a reader can verify independently: that
the reported replicates are reproducible from their seeds; how many replicates
completed and how the failures break down; and the full sampling behaviour of
each estimator, including calibration and coverage, rather than bias and
coverage alone. A final section characterizes the stability of the guarded
augmentation matrix.

## Contents

- `odrcc-reproducibility-report.html` — the rendered report; open it directly,
  no code or data required.
- `odrcc-reproducibility-report.Rmd` — the source that produces the report.
- `diagnostics-lambda.R` — recomputes, for each reported replicate, the matrix
  inverted by the guarded AIPW estimator, together with the weight behaviour
  (effective sample size, extreme weights) that drives instability. Each
  replicate is regenerated from its seed and validated against the stored
  output.
- `summarize-lambda-diagnostics.R` — summarizes the per-replicate diagnostics
  by censoring rate.
- `data/lambda-diagnostics.csv` — the per-replicate stability diagnostics
  (one row per reported replicate). The Lambda section of the report runs from
  this file.

## Reproducing the full report

The completeness audit, the sampling-distribution figures, and the coverage
tables read the per-replicate estimator output, named
`est_results_n_<n>_<c_int>_seed_<seed>.csv`, together with the scheduler logs.
That output is about 218 MB and is not committed here; set the `results_dir`
and `slurm_dir` parameters in the report header to its location, then knit with

## Full simulation output

The complete per-replicate output is attached to the repository release, split
into two independently extractable archives:

- [odrcc-sim-main.tar.gz](https://github.com/jesusepfvazquez/dep-right-censored-covariate/releases/download/sim-archive-v1/odrcc-sim-main.tar.gz) (12 MB): main study output plus all driver and estimator code.
- [odrcc-sim-sensitivity.tar.gz](https://github.com/jesusepfvazquez/dep-right-censored-covariate/releases/download/sim-archive-v1/odrcc-sim-sensitivity.tar.gz) (23 MB): sensitivity study output.

Release: https://github.com/jesusepfvazquez/dep-right-censored-covariate/releases/tag/sim-archive-v1

```r
rmarkdown::render("odrcc-reproducibility-report.Rmd")
```

The first knit assembles the estimates and caches them to an `.rds` file, so
later knits are fast. The Lambda section runs from the bundled
`data/lambda-diagnostics.csv` and needs no additional data. The full result
tree is available from the authors on request.

## Seeding

Each replicate is seeded through its `(n, c_int, seed)` triple, which is
encoded in the output file name, so any summary regenerates exactly. The
diagnostics script reseeds and regenerates each dataset and confirms the
regenerated censoring fraction matches the stored value to machine precision.
