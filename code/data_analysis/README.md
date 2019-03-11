## Running the data analyses in the `paramedic` paper

This file describes how to reproduce the data analyses in the "Estimation of microbial abundances from 16s data and qPCR abundances" by Williamson, Willis, and Hughes. While the code in this file assumes that the user is submitting batch jobs to a high-performance computing (HPC) cluster using the Slurm batch scheduing system, minor edits to these commands allow the use of either local or alternative HPC cluster environments.

The main workhorse function for these simulations is `qpcr_analysis.R`, which allows the user to specify all of the arguments that may vary in a data analysis. The shell script `qpcr_data_analysis.sh` is simply tasked with running `qpcr_analysis.R` using `Rscript`, using user-provided command-line arguments. This code may be run locally. The code detailed in the next two sections describes how to run the data analyses from the manuscript on a HPC cluster.

## Primary analysis of 17 taxa from 1213 women

## Leave-one-out analysis of 7 taxa from 1213 women