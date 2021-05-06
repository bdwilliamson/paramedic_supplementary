#!/bin/bash

## submit all data analyses for the dataset from McClelland et al.
# command line arg: adjust for case-control status (1) or not (0)

## first, submit the main data analysis
## Args:
## 1: q
## 2: nfolds
## 3: sample
## 4: qobs
## 5: niter
## 6: nburn
## 7: nchains
## 8: leave-one-out
## 9: div_num
## 10: save_model
## 11: max_tree
## 12: adapt_delta
## 13: alpha_sigma
## 14: kappa_sigma
./submit_qpcr_data_analysis.sh 127 1 55 13 20000 18000 4 999 1000 1 18 0.9 1 $1 $2 $3

## next, submit the leave-one-out analyses
for i in {1..13}
do
    ./submit_qpcr_data_analysis.sh 13 1 55 13 16000 14000 4 $i 1000 1 18 0.9 1 $1 $2 $3
done
