#!/bin/bash

ml R/3.4.3-foss-2016b-fh1

# arg 1 is the stan model (actually arg 1)
# arg 2 is the sample size (actually arg 2)
# arg 3 is the q (specified in this document)
# arg 4 is the q_obs (actually arg 3)
# arg 5 is the covariance (actually arg 4)
# arg 6 is the variance of missing efficiency (specified in this document)
# arg 7 is number of chains (actually arg 5)
# arg 8 is number of iterations per chain (actually arg 6)
# arg 9 is number of warmups per chain (actually arg 7)
# arg 10 is parallel or not (actually arg 8)
# arg 11 is the number of jobs (actually arg 9)
# arg 12 is use most abundant (actually arg 10)
# arg 13 is number of qs (actually arg 11)
# arg 14 is number of qobss (actually arg 12)
# arg 15 is number of ns (actually arg 13)
# arg 16 is number of sigma e's (actually arg 14)

## submit batch jobs
sbatch -M beagle -c4 --array=1-${9} --requeue ./shell/qpcr_sim.sh vary_sigma_e $1 $2 40 $3 $4 0 $5 $6 ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14}

sbatch -M beagle -c4 --array=1-${9} --requeue ./shell/qpcr_sim.sh vary_sigma_e $1 $2 40 $3 $4 0.1 $5 $6 ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14}

sbatch -M beagle -c4 --array=1-${9} --requeue ./shell/qpcr_sim.sh vary_sigma_e $1 $2 40 $3 $4 0.2 $5 $6 ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14}

sbatch -M beagle -c4 --array=1-${9} --requeue ./shell/qpcr_sim.sh vary_sigma_e $1 $2 40 $3 $4 0.3 $5 $6 ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14}

sbatch -M beagle -c4 --array=1-${9} --requeue ./shell/qpcr_sim.sh vary_sigma_e $1 $2 40 $3 $4 0.4 $5 $6 ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14}

sbatch -M beagle -c4 --array=1-${9} --requeue ./shell/qpcr_sim.sh vary_sigma_e $1 $2 40 $3 $4 0.5 $5 $6 ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14}