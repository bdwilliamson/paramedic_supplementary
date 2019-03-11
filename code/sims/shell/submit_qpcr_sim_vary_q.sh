#!/bin/bash

ml R/3.4.3-foss-2016b-fh1

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

# args specified here:
# number of taxa with qpcr, sixth position

## submit batch jobs
sbatch -M beagle -c4 --array=1-${10} --requeue ./qpcr_sim.sh vary_q $1 $2 $3 2 $4 $5 $6 $7 ${8} ${9} ${10} ${11} ${12} ${13} ${14} ${15}

sbatch -M beagle -c4 --array=1-${10} --requeue ./qpcr_sim.sh vary_q $1 $2 $3 3 $4 $5 $6 $7 ${8} ${9} ${10} ${11} ${12} ${13} ${14} ${15}

sbatch -M beagle -c4 --array=1-${10} --requeue ./qpcr_sim.sh vary_q $1 $2 $3 4 $4 $5 $6 $7 ${8} ${9} ${10} ${11} ${12} ${13} ${14} ${15}

sbatch -M beagle -c4 --array=1-${10} --requeue ./qpcr_sim.sh vary_q $1 $2 $3 5 $4 $5 $6 $7 ${8} ${9} ${10} ${11} ${12} ${13} ${14} ${15}

sbatch -M beagle -c4 --array=1-${10} --requeue ./qpcr_sim.sh vary_q $1 $2 $3 6 $4 $5 $6 $7 ${8} ${9} ${10} ${11} ${12} ${13} ${14} ${15}

sbatch -M beagle -c4 --array=1-${10} --requeue ./qpcr_sim.sh vary_q $1 $2 $3 7 $4 $5 $6 $7 ${8} ${9} ${10} ${11} ${12} ${13} ${14} ${15}