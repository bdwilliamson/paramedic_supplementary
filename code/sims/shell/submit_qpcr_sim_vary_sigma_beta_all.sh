#!/bin/bash

ml R/3.4.3-foss-2016b-fh1
ml jbigkit

## all args specified here;

## submit batch jobs
## q = 10
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_beta "../stan/predict_qpcr_noncentered.stan" 50 10 7 0 0 4 10500 10000 1 50 1 3 1 1 1 4
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_beta "../stan/predict_qpcr_noncentered_sigma_beta_1.stan" 50 10 7 0 0 4 10500 10000 1 50 1 3 1 1 1 4
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_beta "../stan/predict_qpcr_noncentered_sigma_beta_5.stan" 50 10 7 0 0 4 10500 10000 1 50 1 3 1 1 1 4
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_beta "../stan/predict_qpcr_noncentered_sigma_beta_10.stan" 50 10 7 0 0 4 10500 10000 1 50 1 3 1 1 1 4

## q = 20
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_beta "../stan/predict_qpcr_noncentered.stan" 50 20 7 0 0 4 10500 10000 1 50 1 3 1 1 1 4
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_beta "../stan/predict_qpcr_noncentered_sigma_beta_1.stan" 50 20 7 0 0 4 10500 10000 1 50 1 3 1 1 1 4
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_beta "../stan/predict_qpcr_noncentered_sigma_beta_5.stan" 50 20 7 0 0 4 10500 10000 1 50 1 3 1 1 1 4
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_beta "../stan/predict_qpcr_noncentered_sigma_beta_10.stan" 50 20 7 0 0 4 10500 10000 1 50 1 3 1 1 1 4

## q = 40
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_beta "../stan/predict_qpcr_noncentered.stan" 50 40 7 0 0 4 10500 10000 1 50 1 3 1 1 1 4
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_beta "../stan/predict_qpcr_noncentered_sigma_beta_1.stan" 50 40 7 0 0 4 10500 10000 1 50 1 3 1 1 1 4
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_beta "../stan/predict_qpcr_noncentered_sigma_beta_5.stan" 50 40 7 0 0 4 10500 10000 1 50 1 3 1 1 1 4
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_sigma_beta "../stan/predict_qpcr_noncentered_sigma_beta_10.stan" 50 40 7 0 0 4 10500 10000 1 50 1 3 1 1 1 4