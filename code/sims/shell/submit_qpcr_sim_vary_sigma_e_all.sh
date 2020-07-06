#!/bin/bash

ml R/3.4.3-foss-2016b-fh1
ml jbigkit

## all args specified here
# arg 1 is the stan model
# arg 2 is the sample size
# arg 3 is the q
# arg 4 is the covariance
# arg 5 is the variance of missing efficiency
# arg 6 is number of chains
# arg 7 is number of iterations per chain
# arg 8 is number of warmups per chain
# arg 9 is parallel or not
# arg 10 is the number of jobs
# arg 11 is use most abundant
# arg 12 is number of qs
# arg 13 is number of qobss
# arg 14 is number of ns

## submit batch jobs
## no varying efficiency
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_noncentered.stan" 100 40 7 0 0 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_noncentered.stan" 100 40 7 0 0.1 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_noncentered.stan" 100 40 7 0 0.2 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_noncentered.stan" 100 40 7 0 0.3 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_noncentered.stan" 100 40 7 0 0.4 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_noncentered.stan" 100 40 7 0 0.5 4 10500 10000 1 50 1 1 1 1 6 1

# add on later: sigma_e = 0.6, 0.8, 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_noncentered.stan" 100 40 7 0 0.6 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_noncentered.stan" 100 40 7 0 0.8 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_noncentered.stan" 100 40 7 0 1 4 10500 10000 1 50 1 1 1 1 6 1


## varying efficiency
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_with_varying_efficiency_noncentered.stan" 100 40 7 0 0 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_with_varying_efficiency_noncentered.stan" 100 40 7 0 0.1 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_with_varying_efficiency_noncentered.stan" 100 40 7 0 0.2 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_with_varying_efficiency_noncentered.stan" 100 40 7 0 0.3 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_with_varying_efficiency_noncentered.stan" 100 40 7 0 0.4 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_with_varying_efficiency_noncentered.stan" 100 40 7 0 0.5 4 10500 10000 1 50 1 1 1 1 6 1

# add on later: sigma_e = 0.6, 0.8, 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_with_varying_efficiency_noncentered.stan" 100 40 7 0 0.6 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_with_varying_efficiency_noncentered.stan" 100 40 7 0 0.8 4 10500 10000 1 50 1 1 1 1 6 1
sbatch  -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_e "../stan/predict_qpcr_with_varying_efficiency_noncentered.stan" 100 40 7 0 1 4 10500 10000 1 50 1 1 1 1 6 1
