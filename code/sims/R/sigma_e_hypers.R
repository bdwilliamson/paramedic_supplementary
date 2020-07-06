#!/usr/local/bin/Rscript

# simulation testing hyperparameters for sigma_e distribution (sigma_e_hypers.R)
# N = 100, q = 40, q_obs = 7
# correctly specified models for mu, e, V, W
# vary alpha_sigma, kappa_sigma

## load required libraries
library("methods")
library("StanHeaders", lib.loc = "/home/bwillia2/R/x86_64-pc-linux-gnu-library/3.5")
library("rstan")
library("MASS")
library("argparse")
# devtools::install("statdivlab/paramedic@v0.0.3.900")
library("paramedic")
library("tibble")
library("dplyr")
library("Cairo")

## grab in the job id
job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## load in required functions
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) { # if running locally
  source("code/R/data_gen_funcs.R")
  source("code/R/data_generator.R")
  source("code/R/naive_qpcr_estimator.R")
  job_id <- 10
} else {
  source("data_gen_funcs.R")
  source("data_generator.R")
  source("naive_qpcr_estimator.R")
}

## set up dynamic simulation arguments
parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "misspec-normal-normal-negbin-mult", help = "name of the simulation")
parser$add_argument("--N", type = "double", default = 100, help = "sample size")
parser$add_argument("--q", type = "double", default = 40, help = "how many taxa do we have?")
parser$add_argument("--q_obs", type = "double", default = 7, help = "how many taxa do we have qPCR measured on?")
parser$add_argument("--n-chains", type = "double", default = 1, help = "number of chains for MCMC")
parser$add_argument("--iter", type = "double", default = 10000, help = "number of iterations per chain")
parser$add_argument("--warmup", type = "double", default = 5000, help = "number of warmup iterations per chain")
parser$add_argument("--B", type = "double", default = 50, help = "total number of MC reps per q, q^obs")
parser$add_argument("--adapt-delta", type = "double", default = 0.85, help = "adapt_delta, for Stan fitting")
parser$add_argument("--max-treedepth", type = "integer", default = 15, help = "max_treedepth, for Stan fitting")
parser$add_argument("--m-min", type = "integer", default = 10000, help = "min M (reads)")
parser$add_argument("--m-max", type = "integer", default = 100000, help = "max M (reads)")
args <- parser$parse_args()

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) { # if running locally
  args$iter <- 100
  args$warmup <- 50
  args$adapt_delta <- 0.85
  args$max_treedepth <- 15
  args$n_chains <- 1
  args$q <- 40
  args$q_obs <- 7
  args$N <- 50
  args$sim_name <- "sigma-e-1-1"
}
print(paste0("Running sim ", args$sim_name, " with q = ", args$q, "; q_obs = ", args$q_obs, "; N = ", args$N, "."))
print(paste0("Stan parameters: adapt_delta = ", args$adapt_delta, "; treedepth = ", args$max_treedepth, "; num. burnin iter = ", args$warmup, "; num. total iter = ", args$iter, "; number of chains = ", args$n_chains))

options(mc.cores = parallel::detectCores())

# ----------------------------------------------------------
# set up static simulation arguments
# ----------------------------------------------------------
distns <- get_mu_e_distn("correct-normal-normal-poisson-mult")
set.seed(4747)
hypers <- specify_beta_sigma(distns)

# ----------------------------------------------------------
# set up dynamic sim arguments, get and set seed
# ----------------------------------------------------------
alpha_kappa_sigma <- expand.grid(alpha_sigma = seq(1, 5, 1),
                                 kappa_sigma = c(1, .5))
sim_names <- paste0("sigma-e-", apply(alpha_kappa_sigma, 1, function(x) paste0(x, collapse = "-")))
# get current alpha/kappa
current_alpha_kappa <- alpha_kappa_sigma[which(args$sim_name == sim_names), ]
# get the random number seed
current_seed <- 10 * args$N + 1000 * which(args$sim_name == sim_names) + job_id
print(paste0("Current seed: ", current_seed))
set.seed(current_seed)

## create the data
dataset_with_truth <- data_generator(sample_size = args$N, num_taxa = args$q, 
                                     num_qpcr = args$q_obs, seed = current_seed, 
                                     hyper_mean_mu = hypers$beta, 
                                     hyper_cov_mu = hypers$Sigma, 
                                     hyper_sigma = hypers$hyper_sigma, 
                                     hyper_m_min = args$m_min, hyper_m_max = args$m_max,
                                     mu_dist = distns$mu, e_dist = distns$e,
                                     v_dist = distns$v, w_dist = distns$w)
dataset <- dataset_with_truth[(names(dataset_with_truth) %in% c("V", "W", "N", "q", "q_obs"))]
# hard-code prior distribution specification for hyperparameters
dataset$sigma_beta <- sqrt(50.0)
dataset$sigma_Sigma <- sqrt(50.0)
dataset$alpha_sigma <- current_alpha_kappa$alpha_sigma
dataset$kappa_sigma <- current_alpha_kappa$kappa_sigma
## create output directory
prefix <- sprintf("%s/n_%s/q_%s/q_obs_%s",
                  args$sim_name,
                  args$N,
                  args$q,
                  args$q_obs)
fast_prefix <- paste0("/fh/fast/huang_y/user/bwillia2/qPCR/", prefix)
if (!dir.exists(fast_prefix)) {
  dir.create(fast_prefix, recursive = TRUE)
}

stan_seed <- current_seed

## set up the parameters to save
params_to_save <- c("mu", "e", "beta_0", "Sigma")

## run the naive estimator for initial values
naive <- cbind(dataset$V, apply(matrix((args$q_obs+1):args$q), 1, naive_estimator, dataset$W, dataset$V, 1:args$q_obs))
log_naive <- ifelse(is.infinite(log(naive)), 0, log(naive))
set.seed(4747)
eps <- 1e-1
error_mu_tilde <- t(replicate(args$N, rnorm(args$q, 0, 10*eps)))
error_beta <- rnorm(args$q, 0, eps)
error_sigma <- rnorm(args$q, 0, eps)
naive_beta <- colMeans(log_naive, na.rm = TRUE) + error_beta
naive_sigma <- diag(var(log_naive, na.rm = TRUE)) + error_sigma
naive_sigma <- ifelse(naive_sigma <= 0, 1e-2, naive_sigma)
## numerator: sweep(log_naive, 2, naive_beta, FUN = "-"); add some random noise to everything, just in case
log_mu_tilde <- sweep(sweep(log_naive, 2, naive_beta, FUN = "-"), 2, naive_sigma, FUN = "/") + error_mu_tilde

# set up inits list
if (args$n_chains > 1) {
  inits_list <- list(list(log_mu_tilde = log_mu_tilde),
                     list(beta_0 = naive_beta),
                     list(Sigma = naive_sigma),
                     list(init = "random"))
} else {
  inits_list <- list(list(log_mu_tilde = log_mu_tilde, beta_0 = naive_beta, Sigma = naive_sigma))
}

# run the model
set.seed(stan_seed)
colnames(dataset$W) <- paste0("X", 1:ncol(dataset$W))
colnames(dataset$V) <- paste0("X", 1:ncol(dataset$V))
W_tib <- tibble::add_column(subj_id = 1:args$N, as_tibble(dataset$W, .name_repair = "minimal"), .before = 1)
V_tib <- tibble::add_column(subj_id = 1:args$N, as_tibble(dataset$V, .name_repair = "minimal"), .before = 1)
system.time(mod <- paramedic::run_paramedic(W = W_tib, V = V_tib,
                                            n_iter = args$iter, n_burnin = args$warmup,
                                            n_chains = args$n_chains, stan_seed = stan_seed,
                                            inits_lst = inits_list, sigma_beta = dataset$sigma_beta,
                                            sigma_Sigma = dataset$sigma_Sigma, alpha_sigma = dataset$alpha_sigma,
                                            kappa_sigma = dataset$kappa_sigma,
                                            control = list(adapt_delta = args$adapt_delta,
                                                           max_treedepth = args$max_treedepth),
                                            verbose = FALSE, open_progress = FALSE))
mod_summ <- summary(mod, probs = c(0.025, 0.975))$summary

# extract samples
samps <- extract(mod)

## save the output
saveRDS(mod_summ, file = sprintf("%s/%s_mod_jobid_%d_ad_%f_mt_%d.rds",
                                 fast_prefix,
                                 "run_paramedic",
                                 job_id,
                                 args$adapt_delta,
                                 args$max_treedepth))
saveRDS(dataset_with_truth, file = sprintf("%s/%s_data_jobid_%d_ad_%f_mt_%d.rds",
                                           fast_prefix,
                                           "run_paramedic",
                                           job_id,
                                           args$adapt_delta,
                                           args$max_treedepth))
saveRDS(samps, file = sprintf("%s/%s_samps_jobid_%d_ad_%f_mt_%d.rds",
                              fast_prefix,
                              "run_paramedic",
                              job_id,
                              args$adapt_delta,
                              args$max_treedepth))

trace_plot_nms <- c("mu", "e", "beta_0", "Sigma")
fig_width <- fig_height <- 2590
cex <- 1.5
## sample 5 participants to check traceplots for
set.seed(4747)
samp <- sample(1:args$N, 2)
taxa_of_interest <- 10
if (args$q == 20) {
  taxa_of_interest <- c(taxa_of_interest, 20)
} else if (args$q == 40) {
  taxa_of_interest <- c(taxa_of_interest, 20, 40)
} else {
  taxa_of_interest <- c(taxa_of_interest, 20, 40, 60)
}

trace_plots_dir <- "/fh/fast/huang_y/user/bwillia2/qPCR/trace_plots/"
if (!dir.exists(trace_plots_dir)) {
  dir.create(trace_plots_dir, recursive = TRUE)
}
for (n in 1:length(trace_plot_nms)) {
  for (i in taxa_of_interest) {
    logi <- grepl(trace_plot_nms[n], names(mod)) & !grepl("log", names(mod)) & grepl(i, names(mod))
    if (trace_plot_nms[n] == "mu") {
      logi <- grepl(trace_plot_nms[n], names(mod)) & !grepl("log", names(mod)) & grepl(paste0(",", i, "]"), names(mod), fixed = TRUE)
    }
    if (trace_plot_nms[n] == "e") {
      logi <- grepl(trace_plot_nms[n], names(mod)) & !grepl("log", names(mod)) & !grepl("beta", names(mod)) & grepl(paste0("[", i, "]"), names(mod), fixed = TRUE)
    }
    if (trace_plot_nms[n] == "mu") {
      CairoPNG(paste0(trace_plots_dir, "/", args$sim_name, "_run_paramedic", "_par_", trace_plot_nms[n], "_taxon_", i, "_mcid_", job_id, "_q_", args$q, "_q_obs_", args$q_obs, ".png"), width = fig_width, height = fig_height, res = 300, units = "px")
      plot(traceplot(mod, pars = names(mod)[logi][samp]), cex = cex)
      dev.off()
    } else {
      CairoPNG(paste0(trace_plots_dir, "/", args$sim_name, "_run_paramedic", "_par_", trace_plot_nms[n], "_taxon_", i, "_mcid_", job_id, "_q_", args$q, "_q_obs_", args$q_obs, ".png"), width = fig_width, height = fig_height, res = 300, units = "px")
      plot(traceplot(mod, pars = names(mod)[logi]), cex = cex)
      dev.off()
    }
  }
}
