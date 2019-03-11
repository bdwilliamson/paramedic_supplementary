#!/bin/bash

## command line arg 1 is for q (less than 433 selects the q - 7 most abundant taxa outside of 1:7)
## command line arg 2 is for the number of folds to run
## command line arg 3 is how many women to sample

q=${1}
qobs=${2}
niter=${3}
nburn=${4}
nchains=${5}
loo=${6}
div_num=${7}
save_model=${8}


ml R/3.4.3-foss-2016b-fh1

sbatch --time=3-0 ./qpcr_data_analysis.sh "naive" 1 "$nchains" "$niter" "$nburn" "$q" "$sample" "$qobs" "$loo" "$div_num" "$save_model"

sbatch -c6 -p largenode --mem 33G --time=7-0 ./qpcr_data_analysis.sh "no_ve" 1 "$nchains" "$niter" "$nburn" "$q" "$sample" "$qobs" "$loo" "$div_num" "$save_model"

sbatch -c6 -p largenode --mem 33G --time=7-0 ./qpcr_data_analysis.sh "ve" 1 "$nchains" "$niter" "$nburn" "$q" "$sample" "$qobs" "$loo" "$div_num" "$save_model"