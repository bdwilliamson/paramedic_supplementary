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
code_dir <- "code/"
data_dir <- "data/"
stan_dir <- "stan/"
library("methods")
library("rstan")
library("argparse")
library("Cairo")
library("paramedic")

source(paste0(code_dir, "naive_qpcr_estimator.R"))
source(paste0(code_dir, "data_analysis/R/get_most_abundant_taxa.R"))

parser <- ArgumentParser()
parser$add_argument("--estimator", default = "naive", help = "the estimator to calculate")
parser$add_argument("--do-parallel", type = "integer", default = 1, help = "parallelize?")
parser$add_argument("--n-chains", type = "double", default = 6, help = "number of chains")
parser$add_argument("--n-iter", type = "double", default = 10500, help = "number of iterations")
parser$add_argument("--n-burnin", type = "double", default = 10000, help = "number of burn-in")
parser$add_argument("--q", type = "double", default = 60, help = "number of taxa")
parser$add_argument("--sample-num", type = "double", default = 100, help = "sample sample-num women or not")
parser$add_argument("--q-obs", type = "double", default = 7, help = "number of observed qPCR taxa")
parser$add_argument("--leave-one-out", type = "double", default = 999, help = "index to leave out (99 means do not do leave-one-out)")
parser$add_argument("--div-num", type = "double", default = 1000, help = "number to rescale qPCR by")
parser$add_argument("--save-stan-model", type = "double", default = 0, help = "save stan model (1) or not (0)")
args <- parser$parse_args()

args$do_parallel <- as.logical(args$do_parallel)
args$save_stan_model <- as.logical(args$save_stan_model)

print(args)

## -----------------------------------------------------------------------------------
## load in the data, clean it up
## -----------------------------------------------------------------------------------
full_data <- read.csv(file=paste0(data_dir, "p2_brngs_qPCR_merge_20180314.csv"))

## set m_min, llod for processing
m_min <- 1000
llod <- 0
## NEED TO MULTIPLY BY DIV_NUM AT END

## get qpcr, br16s indices in full_data; get the br16s indices that also have qpcr measured
qpcr_inds <- 481:494
br_inds <- 2:476
pcr_plus_br_inds <- c(198, 476, 3, 2, 27, 64, 217) - 1
## process the data, yielding qpcr and br16s matrices
data_lst <- paramedic::process_data(full_data, br_inds, qpcr_inds, pcr_plus_br_inds, llod, m_min, args$div_num)
qpcr <- data_lst$qpcr
br <- data_lst$br

## set up matrices without ids
qpcr_mat <- as.matrix(qpcr[, 2:dim(qpcr)[2]])
br16_mat <- as.matrix(br[, 2:dim(br)[2]])
br16_mat_nms <- colnames(br16_mat)
## calculate the read numbers
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
## sample women
samp <- sample(1:dim(br16_mat)[1], args$sample_num)

## set stan args
stan_seeds <- c(1, 2)
adapt_delta <- 0.85
max_treedepth <- 15

## set up the data list for stan; V needs to be all counts
stan_v <- qpcr_mat[, observed_taxa]
mode(stan_v) <- "integer"
if (args$estimator != "naive") {
  stan_data_lst <- list(W = br16_mat[samp, ], V = stan_v[samp, ], N = args$sample_num, q = q, q_obs = q_obs)    
} else {
  stan_data_lst <- list()
}


## set up parallel
if (args$do_parallel) {
  options(mc.cores = parallel::detectCores())
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
    system.time(mod <- paramedic::paramedic(W = stan_data_lst$W, V = stan_data_lst$V,
                                            N = stan_data_lst$N, q = stan_data_lst$q,
                                            q_obs = stan_data_lst$q_obs,
                                            stan_model = paste0(stan_dir, "predict_qpcr.stan"), 
                                            n_iter = args$n_iter, n_burnin = args$n_burnin, n_chains = args$n_chains, 
                                            stan_seed = stan_seeds[1],
                                            params_to_save = c("mu", "beta", "Sigma"),
                              control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
                              verbose = FALSE, open_progress = FALSE))
    mod_summ <- summary(mod, probs = c(0.025, 0.975))$summary
    samps <- rstan::extract(mod)
} else {
    system.time(mod <- paramedic::paramedic(n_iter = args$n_iter, n_burnin = args$n_burnin, n_chains = args$n_chains, 
                                            stan_seed = stan_seeds[1],
                                            params_to_save = c("mu", "beta", "Sigma", "e", "sigma"),
                              control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
                              verbose = FALSE, open_progress = FALSE))
    mod_summ <- summary(mod, probs = c(0.025, 0.975))$summary
    samps <- rstan::extract(mod)
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
