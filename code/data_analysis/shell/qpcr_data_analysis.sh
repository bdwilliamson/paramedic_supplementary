#!/bin/bash

Rscript ../R/qpcr_analysis_centered.R --estimator ${1} --do-parallel ${2} --n-chains ${3} --n-iter ${4} --n-burnin ${5} --q ${6} --sample-num ${7} --q-obs ${8} --leave-one-out ${9} --div-num ${10} --save-stan-model ${11}
