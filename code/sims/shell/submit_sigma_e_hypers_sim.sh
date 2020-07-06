#!/bin/bash

ml R/3.5.3-foss-2018b-fh1

# args:
# 1: sim name (e.g., 'sigma-e-1-1')
# 2: sample size (100)
# 3: number of taxa (40)
# 4: number of observed taxa (7)
# 5: number of chains (4)
# 6: number of iterations (10500)
# 7: number of burn-in/chain (10000)
# 8: total number of MC reps (50)
# 9: adapt delta (stan control par, .85)
# 10: max treedepth (stan control par, 15)
# 11: min num reads (10000)
# 12: max num reads (100000)

declare -a alphas=(1 2 3 4 5)
declare -a kappas=(1 0.5)
declare -a sims=("sigma-e-3-0.5" "sigma-e-3-1" "sigma-e-2-0.5" "sigma-e-2-1" "sigma-e-1-0.5" "sigma-e-1-1")

# submit each sim
for sim in "${sims[@]}"; do
  io_file="$sim/slurm-%A_%a.out"
  sbatch -c4 --array=1-50 --requeue -e $io_file -o $io_file ./sigma_e_hyper_sim.sh $sim 100 40 7 4 20000 18000 50 .85 18 10000 100000
done
