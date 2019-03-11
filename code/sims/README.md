## Running the numerical experiments for the `paramedic` paper

This file describes how to reproduce the simulations in the "Estimation of microbial abundances from 16s data and qPCR abundances" by Williamson, Willis, and Hughes. While the code in this file assumes that the user is submitting batch jobs to a high-performance computing (HPC) cluster using the Slurm batch scheduing system, minor edits to these commands allow the use of either local or alternative HPC cluster environments.

The main workhorse function for these simulations is `qpcr_sim.R`, which allows the user to specify all of the arguments that may vary in a simulated example. The shell script `qpcr_sim.sh` is simply tasked with running `qpcr_sim.R` using `Rscript`, using user-provided command-line arguments. This code may be run locally. The code detailed in the next two sections describes how to run the numerical experiments from the manuscript on a HPC cluster.

## The effect of varying $q$ and $q^\text{obs}$ with no varying efficiency

This analysis uses `submit_qpcr_sim_vary_q.sh`. The command-line arguments provided are:

1. The Stan model file: "../../stan/predict_qpcr.stan"
2. The sample size: 50
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

The simulation may then be executed from the command line as follows:

```{sh}
./submit_qpcr_sim_vary_q.sh "../../stan/predict_qpcr.stan" 50 10 0 0 4 10500 10000 1 50 1 4 6 1 1

./submit_qpcr_sim_vary_q.sh "../../stan/predict_qpcr.stan" 50 20 0 0 4 10500 10000 1 50 1 4 6 1 1

./submit_qpcr_sim_vary_q.sh "../../stan/predict_qpcr.stan" 50 40 0 0 4 10500 10000 1 50 1 4 6 1 1

./submit_qpcr_sim_vary_q.sh "../../stan/predict_qpcr.stan" 50 60 0 0 4 10500 10000 1 50 1 4 6 1 1
```

This code creates 24 job arrays (one for each unique combination of $q$ and $q^\text{obs} \in \{2, \dots, 7\}$), each with 50 jobs.  

## The effect of varying $\sigma$

This analysis uses `submit_qpcr_sim_vary_sigma_e.sh`. The command-line arguments provided are:

1. The Stan model file: "../../stan/predict_qpcr.stan" or "../../stan/predict_qpcr_with_varying_efficiency.stan" 
2. The sample size: 50
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
./submit_qpcr_sim_vary_sigma_e.sh "../../stan/predict_qpcr.stan" 50 7 0 4 10500 10000 1 50 1 1 1 1 6

./submit_qpcr_sim_vary_sigma_e.sh "../../stan/predict_qpcr_with_varying_efficiency.stan" 50 7 0 4 10500 10000 1 50 1 1 1 1 6
```

This code creates 12 job arrays (one for each unique combination of $\sigma$ and the Stan algorithm), each with 50 jobs.  