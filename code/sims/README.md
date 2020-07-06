# Running the numerical experiments for the `paramedic` paper

This file describes how to reproduce the simulations in the "A multi-view model for relative and absolute microbial abundances" by Williamson, Hughes, and Willis. While the code in this file assumes that the user is submitting batch jobs to a high-performance computing (HPC) cluster using the Slurm batch scheduing system, minor edits to these commands allow the use of either local or alternative HPC cluster environments.

The main workhorse function for these simulations is `qpcr_sim.R`, which allows the user to specify all of the arguments that may vary in a simulated example. The shell script `qpcr_sim.sh` is simply tasked with running `qpcr_sim.R` using `Rscript`, using user-provided command-line arguments. This code may be run locally. The code detailed in the next two sections describes how to run the numerical experiments from the manuscript on a HPC cluster.

# The effect of varying $q$ and $q^\text{obs}$ with no varying efficiency

This analysis replicates the results presented in Section 4.1 of the main manuscript, and uses `shell/submit_qpcr_sim_vary_q_all.sh`. All of the necessary command-line arguments are provided in this file; the arguments are:

1. The Stan model file: "../../stan/predict_qpcr_noncentered.stan"
2. The sample size: 100
3. The number of total taxa: varies
4. The covariance of the taxa: 0
5. The variance of efficiency: 0
6. The number of chains: 4 (or the number of chains you wish to run)
7. The number of iterations per chain: 10500
8. The number of warmups per chain: 10000
9. Should jobs be run in parallel: 1 (0 means no)
10. The total number of monte-carlo replicates: 50
11. Should the most abundant taxa be observed: 1 (0 means no)
12. The number of unique values of $q$ in the simulation: 4
13. The number of unique values of $q^\text{obs}$ in the simulation: 6
14. The number of unique values of $n$ in the simulation: 1
15. The number of unique values of $\sigma_e$ in the simulation: 1

Navigate to the `/shell` subdirectory. The simulation may then be executed from the command line as follows:

```{sh}
./submit_qpcr_sim_vary_q_all.sh
```

This code creates 24 job arrays (one for each unique combination of $q$ and $q^\text{obs} \in \{2, \dots, 7\}$), each with 50 jobs. We use the Stan model `../stan/predict_qpcr_noncentered.stan` (which is equivalent to the hierarchichal model proposed in the manuscript) because we observed improved algorithm convergence using this model instead of `../stan/predict_qpcr.stan`. Once you have the results from this simulation, run the code in `R/load_qpcr_sim_vary_q.R` to generate plots and summaries of the results.

# The effect of varying $\sigma$

This analysis replicates the results presented in Section 4.2 of the main manuscript, and uses `shell/submit_qpcr_sim_vary_sigma_e_all.sh`. All of the necessary command-line arguments are provided in this file; the arguments are:

1. The Stan model file: "../../stan/predict_qpcr_noncentered.stan" or "../../stan/predict_qpcr_with_varying_efficiency_noncentered.stan"
2. The sample size: 100
3. The number of total taxa: 40
4. The covariance of the taxa: 0
5. The variance of efficiency: varies
6. The number of chains: 4 (or the number of chains you wish to run)
7. The number of iterations per chain: 10500
8. The number of warmups per chain: 10000
9. Should jobs be run in parallel: 1 (0 means no)
10. The total number of monte-carlo replicates: 50
11. Should the most abundant taxa be observed: 1 (0 means no)
12. The number of unique values of $q$ in the simulation: 1
13. The number of unique values of $q^\text{obs}$ in the simulation: 1
14. The number of unique values of $n$ in the simulation: 1
15. The number of unique values of $\sigma_e$ in the simulation: 6

The simulation may then be executed from the command line as follows:

```{sh}
./submit_qpcr_sim_vary_sigma_e_all.sh
```

This code creates 12 job arrays (one for each unique combination of $\sigma$ and the Stan algorithm), each with 50 jobs. We again use the noncentered versions of the Stan models because we observed improved algorithm convergence using these models instead of the centered versions. Once you have the results from this simulation, run the code in `R/load_qpcr_sim_vary_sigma_e.R` to generate plots and summaries of the results.

# Additional numerical experiments

## Effect of varying priors on $\beta$, $\Sigma$, and $e$

This analysis replicates the results presented in Section 4.3 of the supporting information.

To run the analysis for $\beta$ and $\Sigma$, run `./submit_qpcr_sim_vary_sigma_beta_all.sh`; this submits 12 job arrays, each with 50 jobs. This code, in turn, calls `qpcr_sim.sh`, which in turn calls `R/qpcr_sim.R`. Once you have the results from this simulation, run the code in `R/load_sigma_beta_noncentered.R` to generate plots and summaries of the results.

To run the analysis for $e$, run `./submit_sigma_e_hypers_sim.sh` from the command line; this code submits 6 job arrays, each with 50 jobs. This code, in turn, calls `shell/sigma_e_hyper_sim.sh`, which in turn calls `R/sigma_e_hypers.R` to run the analysis. This code uses `paramedic` version 0.0.3 or higher, which can be installed using

```{r}
remotes::install_github("statdivlab/paramedic", ref = "v0.0.3")
```

Once you have the results from this simulation, run the code in `R/load_sigma_e_hypers_sim.R` to generate plots and summaries of the results.

## Effect of model misspecification

This analysis replicates the results presented in Section 4.5 of the supporting information. To run this analysis, run `./submit_model_misspec_sim.sh` from the command line; this code submits 12 job arrays (each with 50 jobs) corresponding to the unique combinations of data-generating mechanisms considered. This code, in turn, calls `shell/model_misspec_sim.sh` with many arguments at fixed values (e.g., the number of iterations and the number of chains). The `R` code underlying this experiment is in `R/model_misspec_sim.R`. This code uses `paramedic` version 0.0.3 or higher. Once you have the results of this simulation, run the code in `R/load_model_misspec_sim.R` to generate plots and summaries of the results.

## Soft- vs hard-centering for $e$

This analysis replicates the results presented in Section 4.6 of the supporting information. To run this analysis, run `./submit_qpcr_sim_e_constraints.sh` from the command line; this code submits 6 job arrays, each with 50 jobs. This code, in turn, calls `shell/qpcr_sim.R`, which in turn calls `R/qpcr_sim.R` with two different Stan models (one soft-centered on $e$ (the default) and the other hard-centered on $e$). Once you have the results of this simulation, run the code in `R/load_sim_e_constraint.R` to generate plots and summaries of the results.
