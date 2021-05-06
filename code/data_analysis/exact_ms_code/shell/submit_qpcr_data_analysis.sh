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
max_tree=${11}
adapt_delta=${12}
use_precompiled=${13}
adjust=${14}
alpha_sigma=${15}
kappa_sigma=${16}


ml R/3.5.3-foss-2018b-fh1

sbatch --time=3-0 ./qpcr_data_analysis.sh "naive" 0 1 "$nchains" "$niter" "$nburn" "$q" "$nfolds" "$sample" "$qobs" "$loo" "$div_num" "$save_model" "$max_tree" "$adapt_delta" "$use_precompiled" "$adjust" "$alpha_sigma" "$kappa_sigma"

for i in $(seq 1 $nfolds)
do
    sbatch -c6 --mem 100G --time=7-0 ./qpcr_data_analysis.sh "no_ve" 1 $i "$nchains" "$niter" "$nburn" "$q" "$nfolds" "$sample" "$qobs" "$loo" "$div_num" "$save_model" "$max_tree" "$adapt_delta" "$use_precompiled" "$adjust" "$alpha_sigma" "$kappa_sigma"

    sbatch -c6 --mem 100G --time=7-0 ./qpcr_data_analysis.sh "ve" 1 $i "$nchains" "$niter" "$nburn" "$q" "$nfolds" "$sample" "$qobs" "$loo" "$div_num" "$save_model" "$max_tree" "$adapt_delta" "$use_precompiled" "$adjust" "$alpha_sigma" "$kappa_sigma"
done
