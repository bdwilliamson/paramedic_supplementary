#!/bin/bash

ml R/3.4.3-foss-2016b-fh1
ml jbigkit

## all args specified here

## submit batch jobs
## q = 10
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 10 2 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 10 3 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 10 4 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 10 5 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 10 6 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 10 7 0 0 4 10500 10000 1 50 1 4 6 1 1 1

## q = 20
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 20 2 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 20 3 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 20 4 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 20 5 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 20 6 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 20 7 0 0 4 10500 10000 1 50 1 4 6 1 1 1

## q = 40
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 40 2 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 40 3 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 40 4 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 40 5 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 40 6 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 40 7 0 0 4 10500 10000 1 50 1 4 6 1 1 1

## q = 60
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 60 2 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 60 3 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 60 4 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 60 5 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 60 6 0 0 4 10500 10000 1 50 1 4 6 1 1 1
sbatch -M beagle -c4 --array=1-50 --requeue ./qpcr_sim.sh vary_q "../stan/predict_qpcr_noncentered.stan" 100 60 7 0 0 4 10500 10000 1 50 1 4 6 1 1 1