#!/bin/bash

Rscript code/R/analyze_data/qpcr_analysis.R --estimator ${1} --do-parallel ${2} --fold-num ${3} --n-chains ${4} --n-iter ${5} --n-burnin ${6} --q ${7} --num-folds ${8} --sample-num ${9} --q-obs ${10} --leave-one-out ${11} --div-num ${12} --save-stan-model ${13}
