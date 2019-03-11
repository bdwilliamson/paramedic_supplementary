#!/usr/local/bin/Rscript
##################################################################################
## FILE: qpcr_sim.R
##
## CREATED: 30 August 2018 by Brian Williamson
##
## PURPOSE: replicate the qPCR estimation procedure across many simulated datasets
##          for a given Stan model
##
## INPUTS: stan_model -- the Stan model file to use
##         job_id -- which replicate are we?
##         
## OUTPUTS: a dataset and the Stan results
##################################################################################

## load required libraries
library("methods")
library("rstan")
library("MASS")
library("argparse")

## load in required functions
source("data_gen_funcs.R")
source("data_generator.R")

## grab in the job id
job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## set up dynamic simulation arguments
parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "large_q_vary_n", help = "name of the simulation")
parser$add_argument("--stan-model", default = "../stan/predict_qpcr.stan",
    help = "Which Stan model file to use.")
parser$add_argument("--N", type = "double", default = 50, help = "sample size")
parser$add_argument("--q", type = "double", default = 10, help = "how many taxa do we have?")
parser$add_argument("--q_obs", type = "double", default = 7, help = "how many taxa do we have qPCR measured on?")
parser$add_argument("--beta", default = "random", help = "hyperparameter controlling mean of mu distribution")
parser$add_argument("--corr", type = "double", default = 0, help = "hyperparameter controlling the off-diagonal elements of Sigma.")
parser$add_argument("--sigma", type = "double", default = 0, help = "SD of efficiencies (0 = no varying efficiency)")
parser$add_argument("--n-chains", type = "double", default = 1, help = "number of chains for MCMC")
parser$add_argument("--iter", type = "double", default = 10000, help = "number of iterations per chain")
parser$add_argument("--warmup", type = "double", default = 5000, help = "number of warmup iterations per chain")
parser$add_argument("--parallel", type = "integer", default = 0, help = "run in parallel or not")
parser$add_argument("--B", type = "double", default = 50, help = "total number of MC reps per q, q^obs")
parser$add_argument("--num-qs", type = "double", default = 4, help = "total number of qs considered")
parser$add_argument("--num-qobss", type = "double", default = 6, help = "total number of qobs considered")
parser$add_argument("--num-ns", type = "double", default = 1, help = "total number of sample sizes considered")
parser$add_argument("--num-sigmaes", type = "double", default = 1, help = "total number of sigma_es to consider")
parser$add_argument("--use-most-abundant", type = "integer", default = 1, help = "use q^obs most abundant taxa or not")
args <- parser$parse_args()

args$qs <- c(10, 20, 40, 60, 80, 100, 120)[1:args$num_qs]
args$qobss <- (2:7)[1:args$num_qobss]
args$ns <- seq(50, 300, 50)[1:args$num_ns]
args$sigma_es <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)[1:args$num_sigmaes]

print(args)

## set up static simulation arguments
Sigma <- diag(1, nrow = args$q, ncol = args$q) # hyperparameter controlling covariance of mus
adapt_delta <- 0.85 # controls stan HMC
max_treedepth <- 15 # controls stan HMC
m_min <- 10000 # min number of reads
m_max <- 100000 # max number of reads

## set up hyperparameter mean vector
if (args$beta == "random") {
  set.seed(4747)
  beta <- rnorm(args$q, 0, sqrt(50)) # hyperparameter controlling mean of mus
} else {
  beta <- rep(3, args$q)
}
# hyperparameter controls covariance in off-diagonal of mus
Sigma[row(Sigma) != col(Sigma)] <- args$corr 
# if parallel, set up cores
if (args$parallel) {
    options(mc.cores = parallel::detectCores())
}

## get the random number seed
num_q <- length(args$qs)
num_n <- length(args$ns)
num_q_obs <- length(args$qobss)
num_sigma_e <- length(args$sigma_es)
seed_vec <- 1:(args$B*num_q*num_n*num_q_obs*num_sigma_e)
samp_seed_vec <- sample(seed_vec)
seed_mat <- array(samp_seed_vec, dim = c(args$B, num_q_obs, num_q, num_n, num_sigma_e))
logi_q <- ifelse(length(args$q == 1), 1, which(args$q == args$qs))
logi_q_obs <- ifelse(length(args$q_obs == 1), 1, which(args$q_obs == args$q_obss))
current_seed <- seed_mat[job_id, logi_q_obs, logi_q, which(args$N == args$ns), which(args$sigma == args$sigma_es)] + 10000*args$use_most_abundant
print(current_seed)
set.seed(current_seed)

## create the data
dataset_with_truth <- data_generator(args$N, args$q, args$q_obs, current_seed, beta[1:args$q],
    Sigma[1:args$q, 1:args$q], args$sigma, m_min, m_max, args$use_most_abundant)
dataset <- dataset_with_truth[(names(dataset_with_truth) %in% c("V", "W", "N", "q", "q_obs"))]

## create output directory
prefix <- sprintf("%s/ve_%s/cov_%s/n_%s/q_%s/q_obs_%s/", 
    args$sim_name, 
    args$sigma, 
    args$corr,
    args$N,
    args$q,
    args$q_obs)
if (!dir.exists(prefix)) {
    dir.create(prefix, recursive = TRUE)
}

## make a different seed for each column of interest
## job_id is the taxon id * replicate number (i.e., for 7 groups of 2 taxa and 20 reps each, 
## job_id ranges from 1 -- 280)
stan_seed <- current_seed

## set up the parameters to save
params_to_save <- if (grepl("varying_efficiency", args$stan_model)) {
  c("beta", "Sigma", "mu", "e", "sigma")
  } else {
    c("beta", "Sigma", "mu")
}

# run the model
set.seed(stan_seed)
system.time(mod <- stan(file = args$stan_model, data = dataset, iter = args$iter,
    warmup = args$warmup, chains = args$n_chains, seed = stan_seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
        verbose = FALSE, open_progress = FALSE,
    pars = params_to_save))

mod_summ <- summary(mod, probs = c(0.025, 0.975))$summary

# extract samples
samps <- extract(mod)

## save the output
saveRDS(mod_summ, file = sprintf("%s/%s_mod_jobid_%d_ad_%f_mt_%d_ab_%s.rds",
    prefix,
    strsplit(strsplit(args$stan_model, "/")[[1]], ".", fixed = TRUE)[[length(strsplit(strsplit(args$stan_model, "/")[[1]], ".", fixed = TRUE))]][1], 
    job_id,
    adapt_delta,
    max_treedepth,
    as.character(args$use_most_abundant)))
saveRDS(dataset_with_truth, file = sprintf("%s/%s_data_jobid_%d_ad_%f_mt_%d_ab_%s.rds",
    prefix,
    strsplit(strsplit(args$stan_model, "/")[[1]], ".", fixed = TRUE)[[length(strsplit(strsplit(args$stan_model, "/")[[1]], ".", fixed = TRUE))]][1], 
    job_id,
    adapt_delta,
    max_treedepth,
    as.character(args$use_most_abundant)))
saveRDS(samps, file = sprintf("%s/%s_samps_jobid_%d_ad_%f_mt_%d_ab_%s.rds",
    prefix,
    strsplit(strsplit(args$stan_model, "/")[[1]], ".", fixed = TRUE)[[length(strsplit(strsplit(args$stan_model, "/")[[1]], ".", fixed = TRUE))]][1], 
    job_id,
    adapt_delta,
    max_treedepth,
    as.character(args$use_most_abundant)))
# ## make some traceplots (don't do this, because it always errors on beagle)
# trace_plots_dir <- paste0(prefix, "traceplots/")
# options(bitmapType = 'cairo')
# if (!dir.exists(trace_plots_dir)){
#     dir.create(trace_plots_dir, recursive = TRUE)
# }
# trace_plot_nms <- c("mu", "beta", "Sigma")
# if (args$sigma != 0) {
#     trace_plot_nms <- c(trace_plot_nms, "e")
# }
# fig_width <- fig_height <- 2590
# cex <- 1.5
# for (n in trace_plot_nms) {
#     for (i in tail(1:args$q, args$q - args$q_obs)) {
#         logi <- grepl(n, names(mod)) & !grepl("log", names(mod)) & grepl(i, names(mod))
#         if (n == "mu") {
#           logi <- grepl(n, names(mod)) & !grepl("log", names(mod)) & grepl(paste0(",", i, "]"), names(mod), fixed = TRUE)
#         }
#         if (n == "e") {
#           logi <- grepl(n, names(mod)) & !grepl("log", names(mod)) & !grepl("beta", names(mod)) & grepl(paste0("[", i, "]"), names(mod), fixed = TRUE)
#         }
#         png(paste0(trace_plots_dir,strsplit(strsplit(args$stan_model, "/")[[1]], ".", fixed = TRUE)[[length(strsplit(strsplit(args$stan_model, "/")[[1]], ".", fixed = TRUE))]][1],"_ve_", args$sigma, "_id_", job_id, "_", n, "_taxon_", i, ".png"), width = fig_width, height = fig_height, res = 300, units = "px")
#         plot(traceplot(mod, pars = names(mod)[logi]), cex = cex)
#         dev.off()
#   }
# }
# warnings()