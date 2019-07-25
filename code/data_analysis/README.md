## Running the data analyses in the `paramedic` paper

This file describes how to reproduce the data analyses in the "Estimation of microbial abundances from 16s data and qPCR abundances" by Williamson, Hughes, and Willis. While the code in this file assumes that the user is submitting batch jobs to a high-performance computing (HPC) cluster using the Slurm batch scheduing system, minor edits to these commands allow the use of either local or alternative HPC cluster environments.

The main workhorse function for these simulations is `qpcr_analysis.R`, which allows the user to specify all of the arguments that may vary in a data analysis. The shell script `qpcr_data_analysis.sh` is simply tasked with running `qpcr_analysis.R` using `Rscript`, using user-provided command-line arguments. This code may be run locally. The code detailed in the next two sections describes how to run the data analyses from the manuscript on a HPC cluster.

## Primary analysis of 17 taxa from 1213 women

The arguments that we pass to `qpcr_analysis.R` are:

1. the estimator: "ve" to run `predict_qpcr_with_varying_efficiency.stan`; "no_ve" to run `predict_qpcr.stan`; and "naive" to run the naive estimator.
2. whether or not to run in parallel: 1 (0 means no).
3. the number of chains: 4.
4. the number of iterations per chain: 16000.
5. the number of warmups per chain: 14000.
6. the number of total taxa: 17.
7. the number of women to sample: 1213.
8. the number of taxa with observed qPCR: 7.
9. the taxon index to leave out: 999 (a leave-one-out analysis would use 1--7)
10. the number to divide qPCR by to create valid R integers: 1000.
11. whether or not to save the stan model along with output: 1 (0 means no).

We then run the data analysis with:
```{sh}
./shell/submit_qpcr_data_analysis.sh 17 7 10500 10000 6 999 1000 1
```

To load the results of the analysis and produce the tables and figures used in the manuscript, run `.R`.

## Leave-one-out analysis of 7 taxa from 1213 women

The arguments that we pass to `qpcr_analysis.R` are:

1. the estimator: "ve" to run `predict_qpcr_with_varying_efficiency.stan`; "no_ve" to run `predict_qpcr.stan`; and "naive" to run the naive estimator.
2. whether or not to run in parallel: 1 (0 means no).
3. the number of chains: 4.
4. the number of iterations per chain: 16000.
5. the number of warmups per chain: 14000.
6. the number of total taxa: 7.
7. the number of women to sample: 1213.
8. the number of taxa with observed qPCR: 7.
9. the taxon index to leave out: 1--7
10. the number to divide qPCR by to create valid R integers: 1000.
11. whether or not to save the stan model along with output: 1 (0 means no).

We then run the data analysis with:
```{sh}
./shell/submit_qpcr_data_analysis.sh 7 7 10500 10000 6 1 1000 1

./shell/submit_qpcr_data_analysis.sh 7 7 10500 10000 6 2 1000 1

./shell/submit_qpcr_data_analysis.sh 7 7 10500 10000 6 3 1000 1

./shell/submit_qpcr_data_analysis.sh 7 7 10500 10000 6 4 1000 1

./shell/submit_qpcr_data_analysis.sh 7 7 10500 10000 6 5 1000 1

./shell/submit_qpcr_data_analysis.sh 7 7 10500 10000 6 6 1000 1

./shell/submit_qpcr_data_analysis.sh 7 7 10500 10000 6 7 1000 1
```

To load the results of the analysis and produce the tables and figures used in the manuscript, run `R/load_qpcr_leave_one_out_analysis.R`.
