#!/bin/sh

Rscript ../R/analyze_data/load_qpcr_analysis.R --q $1 --q-obs $2 --sample-num $3 --div-num 1000 --adjust $4
