#!/bin/bash

# call the misspecification sim a single time

# args:
# 1: sim name (e.g., 'misspec-gamma-gamma')
# 2: sample size (50)
# 3: number of taxa (40)
# 4: number of observed taxa (7)
# 5: number of chains (4)
# 6: number of iterations (10500)
# 7: number of burn-in/chain (10000)
# 8: total number of MC reps (50 * )
# 9: adapt delta (stan control par, .85)
# 10: max treedepth (stan control par, 15)
# 11: min num reads (10000)
# 12: max num reads (100000)

Rscript ../R/sigma_e_hypers.R --sim-name $1 --N $2 --q $3 --q_obs $4 --n-chains $5 --iter $6 --warmup $7 --B $8 --adapt-delta $9 --max-treedepth ${10} --m-min ${11} --m-max ${12}
