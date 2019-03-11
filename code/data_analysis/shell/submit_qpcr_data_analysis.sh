#!/bin/bash

## command line arg 1 is for q (less than 433 selects the q - 7 most abundant taxa outside of 1:7)
## command line arg 2 is for the number of folds to run
## command line arg 3 is how many women to sample

q=${1}
nfolds=${2}
sample=${3}
qobs=${4}
niter=${5}
nburn=${6}
nchains=${7}
loo=${8}
div_num=${9}
save_model=${10}


ml R/3.4.3-foss-2016b-fh1

sbatch --time=3-0 ./qpcr_data_analysis.sh "naive" 0 1 "$nchains" "$niter" "$nburn" "$q" "$nfolds" "$sample" "$qobs" "$loo" "$div_num" "$save_model"

for i in $(seq 1 $nfolds)
do
    sbatch -c6 -p largenode --mem 33G --time=7-0 ./qpcr_data_analysis.sh "no_ve" 1 $i "$nchains" "$niter" "$nburn" "$q" "$nfolds" "$sample" "$qobs" "$loo" "$div_num" "$save_model"

    sbatch -c6 -p largenode --mem 33G --time=7-0 ./qpcr_data_analysis.sh "ve" 1 $i "$nchains" "$niter" "$nburn" "$q" "$nfolds" "$sample" "$qobs" "$loo" "$div_num" "$save_model"
done
