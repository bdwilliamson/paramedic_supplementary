#!/bin/bash

sim_name=$1
stan_model=$2
N=$3
q=$4
q_obs=$5
corr=$6
sigma=$7
n_chains=$8
iter=$9
warmup=${10}
parallel=${11}
B=${12}
most_abundant=${13}
qs=${14}
qobss=${15}
ns=${16}
sigmaes=${17}

Rscript ../R/qpcr_sim.R --sim-name "$sim_name" --stan-model "$stan_model" --N "$N" --q "$q" --q_obs "$q_obs" --beta "random" --corr "$corr" --sigma "$sigma" --n-chains "$n_chains" --iter "$iter" --warmup "$warmup" --parallel "$parallel" --B "$B" --num-qs "$qs" --num-qobss "$qobss" --num-ns "$ns" --num-sigmaes "$sigmaes" --use-most-abundant "$most_abundant"