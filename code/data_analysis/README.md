# Running the data analyses in the `paramedic` paper

This file describes how to replicate the data analyses in the "A multi-view model for relative and absolute microbial abundances" by Williamson, Hughes, and Willis. While the code in this file assumes that the user is submitting batch jobs to a high-performance computing (HPC) cluster using the Slurm batch scheduing system, minor edits to these commands allow the use of either local or alternative HPC cluster environments.

The main workhorse function for these simulations is `qpcr_analysis.R`, which allows the user to specify all of the arguments that may vary in a data analysis. The shell script `qpcr_data_analysis.sh` is simply tasked with running `qpcr_analysis.R` using `Rscript`, using user-provided command-line arguments. This code may be run locally. The code detailed in the next two sections describes how to run the data analyses from the manuscript on a HPC cluster.

We note here that since the data are not publicly available, we cannot reproduce exactly the results from the data analysis. However, using the `paramedic` R package (version 0.0.0.9000), we can run the data analyses on a dataset inspired by the full data analyzed in the manuscript. The code in the `exact_ms_code` subdirectory is exactly the code used to analyze the data for the manuscript (and supporting information).

## Primary analysis of 17 taxa from 20 samples

The arguments that we pass to `submit_qpcr_data_analysis.sh` (which, in turn, calls `qpcr_data_analysis.sh` and `qpcr_analysis.R`) are:


1. the number of total taxa: 17.
2. whether or not to break the data into folds (1 means no folds)
3. the number of samples: 20.
4. the number of taxa with observed qPCR: 7.
5. the number of iterations per chain: 16000.
6. the number of warmups per chain: 14000.
7. the number of chains: 4.
8. the taxon index to leave out: 999 (a leave-one-out analysis would use 1--7)
9. the number to divide qPCR by to create valid R integers: 1000.
10. whether or not to save the stan model along with output: 1 (0 means no).
11. the maximum treedepth (`max_treedepth` in Stan)
12. the acceptance parameter (`adapt_delta` in Stan)
13. should we use precompiled Stan models from `paramedic`? (0 is no, 1 is yes)
14. should we adjust for covariates? (0 is no, 1 is yes)
15. hyperparameter for efficiency distribution (`alpha_sigma`)
16. hyperparameter for efficiency distribution (`kappa_sigma`)

We then run the data analysis with:
```{sh}
./shell/submit_qpcr_data_analysis.sh 17 1 20 7 16000 14000 4 999 1000 1 20 0.9 1 0 4 3
```

To load the results of the analysis and produce the tables and figures used in the manuscript, run `.R`.

## Leave-one-out analysis of 7 taxa from 20 samples

We then run these analyses with:
```{sh}
for i in {1..7}; do
    ./shell/submit_qpcr_data_analysis.sh 7 1 20 7 16000 14000 4 1 1000 $i 20 0.9 1 0 4 3
done
```

To load the results of the analysis and produce the tables and figures used in the manuscript, run `R/load_qpcr_leave_one_out_analysis.R`.

## Sensitivity analyses

The example data do not have any external covariates (recall that in the main manuscript, HIV-1 status was a relevant covariate that we adjusted for). Thus, we do not describe here how to run an unadjusted vs adjusted analysis.

To investigate the choice of the hyperparameters `alpha_sigma` and `kappa_sigma`, run the data analysis with:
```{sh}
./shell/submit_qpcr_data_analysis.sh 17 1 20 7 16000 14000 4 999 1000 1 20 0.9 1 0 2 3
```
