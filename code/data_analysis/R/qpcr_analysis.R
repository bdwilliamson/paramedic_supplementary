######################################################################################
##
## FILE: qpcr_analysis.R
##
## CREATED: 15 October 2018 by Brian Williamson
##
## PURPOSE: data analysis of the qPCR + br16s data from Hutch collaborators
##          to showcase qPCR estimation method
##
## INPUTS: ../../data/p2_brngs_qPCR_merge_20180314.csv
##
## OUTPUTS:
######################################################################################
## -----------------------------------------------------------------------------------
## load required functions and libraries
## -----------------------------------------------------------------------------------
code_dir <- "code/R/"
data_dir <- "data/"
stan_dir <- "stan/"
library("methods")
library("StanHeaders")
library("rstan")
library("argparse")
library("Cairo")
remotes::install_github("statdivlab/paramedic", ref = "v0.0.3")
library("paramedic")
source(paste0(code_dir, "analyze_data/qpcr_analysis_helper_functions.R"))
source(paste0(code_dir, "naive_qpcr_estimator.R"))
source(paste0(code_dir, "analyze_data/get_most_abundant_taxa.R"))

parser <- ArgumentParser()
parser$add_argument("--estimator", default = "no_ve", help = "the estimator to calculate")
parser$add_argument("--use-precompiled", type = "integer", default = 1, help = "use precompiled stan code?")
parser$add_argument("--adjust", type = "integer", default = 0, help = "adjust for case/control status?")
parser$add_argument("--do-parallel", type = "integer", default = 1, help = "parallelize?")
parser$add_argument("--fold-num", type = "double", default = 1, help = "data fold to run on")
parser$add_argument("--num-folds", type = "double", default = 1, help = "number of data folds")
parser$add_argument("--n-chains", type = "double", default = 6, help = "number of chains")
parser$add_argument("--n-iter", type = "double", default = 10500, help = "number of iterations")
parser$add_argument("--n-burnin", type = "double", default = 10000, help = "number of burn-in")
parser$add_argument("--q", type = "double", default = 127, help = "number of taxa")
parser$add_argument("--sample-num", type = "double", default = 110, help = "sample sample-num women or not")
parser$add_argument("--q-obs", type = "double", default = 13, help = "number of observed qPCR taxa")
parser$add_argument("--leave-one-out", type = "double", default = 999, help = "index to leave out (99 means do not do leave-one-out)")
parser$add_argument("--div-num", type = "double", default = 1000, help = "number to rescale qPCR by")
parser$add_argument("--save-stan-model", type = "double", default = 0, help = "save stan model (1) or not (0)")
parser$add_argument("--max-treedepth", type = "double", default = 15, help = "max treedepth")
parser$add_argument("--adapt-delta", type = "double", default = 0.85, help = "adapt delta")
parser$add_argument("--alpha-sigma", type = "double", default = 2, help = "alpha_sigma: shape parameter for inverse gamma prior on e variance")
parser$add_argument("--kappa-sigma", type = "double", default = 1, help = "kappa_sigma: scale parameter for inverse gamma prior on e variance")
args <- parser$parse_args()

args$do_parallel <- as.logical(args$do_parallel)
args$save_stan_model <- as.logical(args$save_stan_model)
args$adjust <- as.logical(args$adjust)
args$use_precompiled <- as.logical(args$use_precompiled)
print(args)

## -----------------------------------------------------------------------------------
## load in the data, clean it up
## -----------------------------------------------------------------------------------
data("example_qPCR_data")
data("example_br16S_data")
## calculate the read numbers
br16_mat <- as.matrix(example_br16S_data[, 2:dim(example_br16S_data)[2]])
qPCR_mat <- as.matrix(example_qPCR_data[, 2:dim(example_qPCR_data)[2]])
m <- rowSums(br16_mat)

## -----------------------------------------------------------------------------------
## run the estimators
## -----------------------------------------------------------------------------------

## set q, q_obs
q <- args$q
q_obs <- args$q_obs
observed_taxa <- 1:args$q_obs

## if leave_one_out < 999, then do leave-one-out
if (args$leave_one_out < 999) {
    q_obs <- args$q_obs - 1
    observed_taxa <- (1:args$q_obs)[-args$leave_one_out]
}

cat("\n Observed taxa: \n")
print(observed_taxa)

## if q is smaller than the total number, select the most abundant taxa
if (q < dim(br16_mat)[2]) {
  ## get the order
  ordered_by_abundance <- get_most_abundant_taxa(br16_mat, m)
  ## remove the ones corresponding to observed qPCR
  taxa_to_estimate_init <- ordered_by_abundance[!(ordered_by_abundance %in% observed_taxa)]
  ## if leave-one-out, the taxon to estimate is the left-out one
  if (args$leave_one_out < 999) {
    taxa_to_estimate <- args$leave_one_out
  } else {
    taxa_to_estimate <- taxa_to_estimate_init
  }

  ## select the most abundant taxa (always select the first 7, corresponding to qPCR)
  ## if q == q_obs, then only do observed taxa
  if (q == q_obs) {
    most_abundant_16S <- br16_mat[, observed_taxa]
  } else {
    most_abundant_16S <- br16_mat[, c(observed_taxa, taxa_to_estimate[1:(q - q_obs)])]
  }

  ## rename br16_mat
  br16_mat <- most_abundant_16S
  ## re-normalize
  m <- rowSums(br16_mat)
}

cat("\n Taxa to estimate: \n")
print(taxa_to_estimate)

set.seed(4747)
## break up the sample into even chunks
# folds_init <- rep(seq_len(args$num_folds), length = dim(br16_mat)[1])
# folds <- sample(folds_init)
## sample women
samp <- sample(1:dim(br16_mat)[1], args$sample_num)

## set stan args
stan_seeds <- c(1, 2)
adapt_delta <- args$adapt_delta
max_treedepth <- args$max_treedepth

## set up the data list for stan; V needs to be all counts
stan_v <- qpcr_mat[, observed_taxa]
mode(stan_v) <- "integer"
if (args$use_precompiled) {
  n <- nrow(br16_mat[, , drop = FALSE])
  W <- cbind.data.frame(subj_id = (1:n)[samp], br16_mat[samp, , drop = FALSE])
  V <- cbind.data.frame(subj_id = (1:n)[samp], stan_v[samp, , drop = FALSE])
  names(W)[1:(length(observed_taxa) + 1)] <- names(V)[1:(length(observed_taxa) + 1)]
  X <- cbind.data.frame(subj_id = (1:n)[samp], case_control = case_control[samp])
  p <- 1
} else {
  W <- br16_mat[samp, , drop = FALSE]
  V <- stan_v[samp, , drop = FALSE]
  X <- matrix(case_control, ncol = 1)[samp, , drop = FALSE]
  p <- 1
}
if (args$estimator != "naive") {
  stan_data_lst <- list(W = W, V = V,
                        N = length(samp), q = q, q_obs = q_obs,
                        sigma_beta = 1.62, sigma_Sigma = sqrt(50),
                        alpha_sigma = args$alpha_sigma, kappa_sigma = args$kappa_sigma)
} else {
  stan_data_lst <- list(N = length(samp))
}
if (args$adjust) {
  tmp <- c(stan_data_lst, list(X = X, p = 1))
  stan_data_lst <- tmp
} else {
  tmp <- c(stan_data_lst, list(X = V[, 1, drop = FALSE]))
  stan_data_lst <- tmp
}

## set up parallel
if (args$do_parallel) {
  options(mc.cores = parallel::detectCores())
}

## run the naive estimator for initial values
if (q == q_obs) {
  if (args$leave_one_out < 999) {
    naive <- cbind(stan_v[samp, ], apply(matrix(1:q)[-observed_taxa], 1, naive_estimator, br16_mat[samp, ], stan_v[samp, ], observed_taxa))
  } else {
    naive <- stan_v[samp, ]
  }

} else {
  naive <- cbind(stan_v[samp, ], apply(matrix((q_obs+1):q), 1, naive_estimator, br16_mat[samp, ], stan_v[samp, ], observed_taxa))
}
log_naive <- ifelse(is.infinite(log(naive)), 0, log(naive))
set.seed(4747)
eps <- 1e-1
error_mu_tilde <- t(replicate(args$sample_num, rnorm(args$q, 0, 10*eps)))
error_beta <- rnorm(args$q, 0, eps)
error_sigma <- rnorm(args$q, 0, eps)
naive_beta <- colMeans(log_naive, na.rm = TRUE) + error_beta
naive_sigma <- diag(var(log_naive, na.rm = TRUE)) + error_sigma
naive_sigma <- ifelse(naive_sigma <= 0, 1e-2, naive_sigma)
## numerator: sweep(log_naive, 2, naive_beta, FUN = "-"); add some random noise to everything, just in case
log_mu_tilde <- sweep(sweep(log_naive, 2, naive_beta, FUN = "-"), 2, naive_sigma, FUN = "/") + error_mu_tilde

# set up inits list
if (args$n_chains > 1) {
  if (args$n_chains == 4) {
    inits_list <- c(list(list(log_mu_tilde = log_mu_tilde),
                         list(beta = naive_beta),
                         list(Sigma = naive_sigma),
                         list(init = "random")))
  } else {
    inits_list <- c(list(list(log_mu_tilde = log_mu_tilde),
                         list(beta = naive_beta),
                         list(Sigma = naive_sigma)),
                    replicate(args$n_chains - 3, list(init = "random"), simplify = FALSE))
  }

} else {
  inits_list <- list(list(log_mu_tilde = log_mu_tilde, beta = naive_beta, Sigma = naive_sigma))
}

cat("\n Running estimator", args$estimator, "fold", args$fold_num, "\n")
if (args$estimator == "naive") {
  if (q == q_obs) {
    system.time(mod <- cbind(qpcr_mat[samp, observed_taxa], apply(matrix(1:q), 1, naive_estimator, br16_mat[samp, ], qpcr_mat[samp, ], observed_taxa)))
    colnames(mod) <- br16_mat_nms[observed_taxa]
  } else {
    system.time(mod <- cbind(qpcr_mat[samp, observed_taxa], apply(matrix((q_obs+1):q), 1, naive_estimator, br16_mat[samp, ], qpcr_mat[samp, ], observed_taxa)))
    colnames(mod) <- c(br16_mat_nms[observed_taxa], br16_mat_nms[taxa_to_estimate][1:(q - q_obs)])
  }
  samps <- NA
    mod_summ <- mod
} else if (args$estimator == "no_ve") {
  if (args$use_precompiled) {
      system.time(mod <- paramedic::no_efficiency(W = stan_data_lst$W, V = stan_data_lst$V,
                                                  X = stan_data_lst$X, n_iter = args$n_iter,
                                                  n_burnin = args$n_burnin, n_chains = args$n_chains,
                                                  stan_seed = stan_seeds[1],
                                                  inits_lst = inits_list,
                                                  sigma_beta = stan_data_lst$sigma_beta, sigma_Sigma = stan_data_lst$sigma_Sigma,
                                                  control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
                                                  open_progress = FALSE, verbose = FALSE))
  } else {
    system.time(mod <- stan(file = paste0(stan_dir, "predict_qpcr_noncentered.stan"),
                            data = stan_data_lst,
                            iter = args$n_iter, warmup = args$n_burnin, chains = args$n_chains, seed = stan_seeds[1],
                            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
                            verbose = FALSE, open_progress = FALSE,
                            pars = c("mu", "beta", "Sigma"),
                            init = inits_list))
  }
  mod_summ <- summary(mod, probs = c(0.025, 0.975))$summary
  samps <- rstan::extract(mod)
} else if (args$estimator == "ve") {
  if (args$use_precompiled) {
    system.time(mod <- paramedic::run_paramedic(W = stan_data_lst$W, V = stan_data_lst$V,
                                                X = stan_data_lst$X, n_iter = args$n_iter,
                                                n_burnin = args$n_burnin, n_chains = args$n_chains,
                                                stan_seed = stan_seeds[2],
                                                inits_lst = inits_list,
                                                sigma_beta = stan_data_lst$sigma_beta, sigma_Sigma = stan_data_lst$sigma_Sigma,
                                                alpha_sigma = stan_data_lst$alpha_sigma, kappa_sigma = stan_data_lst$kappa_sigma,
                                                control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
                                                open_progress = FALSE, verbose = FALSE))
  } else {
    system.time(mod <- stan(file = paste0(stan_dir, "predict_qpcr_with_varying_efficiency_noncentered.stan"),
                            data = stan_data_lst,
                            iter = args$n_iter, warmup = args$n_burnin, chains = args$n_chains, seed = stan_seeds[2],
                            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
                            verbose = FALSE, open_progress = FALSE,
                            pars = c("mu", "beta", "Sigma", "e", "sigma"),
                            init = inits_list))
  }
  mod_summ <- summary(mod, probs = c(0.025, 0.975))$summary
  samps <- rstan::extract(mod)
} else {
  stop("the estimator requested isn't currently implemented")
}

## save off stan model objects, naive estimators, data
if (args$save_stan_model) {
    save_lst <- list(data = data_lst, mod = mod_summ, stan_data_lst = stan_data_lst, samps = samps, stan_out = mod)
} else {
    save_lst <- list(data = data_lst, mod = mod_summ, stan_data_lst = stan_data_lst, samps = samps, stan_out = NA)
}

saveRDS(save_lst, paste0("qpcr_data_analysis_est_", args$estimator, "_q_", args$q, "_q_obs_", q_obs, "_sample_", args$sample_num, "_loo_", args$leave_one_out, ".rds"))

trace_plot_nms <- c("mu", "beta", "Sigma")
if (args$estimator == "ve") {
    trace_plot_nms <- c(trace_plot_nms, "e")
}

fig_width <- fig_height <- 2590
cex <- 1.5
if (args$estimator == "naive") {

} else {
  trace_plots_dir <- "trace_plots"
  if (!dir.exists(trace_plots_dir)) {
    dir.create(trace_plots_dir, recursive = TRUE)
  }
  for (n in 1:length(trace_plot_nms)) {
    for (i in 1:args$q) {
      logi <- grepl(trace_plot_nms[n], names(mod)) & !grepl("log", names(mod)) & grepl(i, names(mod))
      if (trace_plot_nms[n] == "mu") {
        logi <- grepl(trace_plot_nms[n], names(mod)) & !grepl("log", names(mod)) & grepl(paste0(",", i, "]"), names(mod), fixed = TRUE)
      }
      if (trace_plot_nms[n] == "e") {
        logi <- grepl(trace_plot_nms[n], names(mod)) & !grepl("log", names(mod)) & !grepl("beta", names(mod)) & grepl(paste0("[", i, "]"), names(mod), fixed = TRUE)
      }
      if (sum(logi) > 10) {
        for (j in 1:ceiling(sum(logi)/10)) {
          CairoPNG(paste0(trace_plots_dir, "/", args$estimator, "_par_", trace_plot_nms[n], "_taxon_", i, "_slice_", j, "_q_", args$q, "_q_obs_", q_obs, "_sample_", as.numeric(args$sample), "_loo_", args$leave_one_out, ".png"), width = fig_width, height = fig_height, res = 300, units = "px")
          plot(traceplot(mod, pars = names(mod)[logi][1:10 + (j-1)*10]), cex = cex)
          dev.off()
        }
      } else {
        CairoPNG(paste0(trace_plots_dir, "/", args$estimator, "_par_", trace_plot_nms[n], "_taxon_", i, "_q_", args$q, "_q_obs_", q_obs, "_sample_", as.numeric(args$sample), "_loo_", args$leave_one_out, ".png"), width = fig_width, height = fig_height, res = 300, units = "px")
        plot(traceplot(mod, pars = names(mod)[logi]), cex = cex)
        dev.off()
      }

    }
  }
}
