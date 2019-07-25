#!/usr/local/bin/Rscript
##################################################################################
## FILE: load_sigma_beta.R
##
## CREATED: 02 April 2019 by Brian Williamson
##
## PURPOSE: load in results from a qPCR simulation with varying sigma_beta, sigma_sigma
##
## INPUTS: stan_model -- the Stan model file used
##                 ve -- the sd of the varying efficiency
##                cov -- covariance of the taxa
##                  N -- sample size
##                  q -- total number of taxa
##
## OUTPUTS: a dataset and the Stan results
##################################################################################

## --------------------------------------------------------------------------------------------------
## SETUP
## --------------------------------------------------------------------------------------------------
## load required functions and packages
library("methods")
library("argparse")
library("Cairo")
library("tidyr")
library("dplyr")
library("ggplot2")
library("cowplot")
source("code/R/load_qpcr_sim_helpers.R") # provides get_summaries
source("code/R/naive_qpcr_estimator.R")
source("code/R/gen_wald_interval.R")
source("code/R/gen_quantile_interval.R")
source("code/R/plot_qpcr_sim_helpers.R")
source("code/R/analyze_data/get_most_abundant_taxa.R")
source("code/R/gen_naive_interval.R")

## grab command-line arguments
parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "vary_sigma_beta", help = "name of the simulation")
parser$add_argument("--stan-model", nargs = "+", default = c("predict_qpcr_noncentered_sigma_beta_1",
                                                             "predict_qpcr_noncentered_sigma_beta_5",
                                                             "predict_qpcr_noncentered_sigma_beta_10",
                                                             "predict_qpcr_noncentered",
                                                             "predict_qpcr_noncentered_sigma_beta_100",
                                                             "predict_qpcr_noncentered_sigma_beta_144",
                                                             "predict_qpcr_noncentered_sigma_beta_196"),
                    help = "Which Stan model file to use.")
parser$add_argument("--N", type = "double", default = 50, help = "sample size")
parser$add_argument("--q", nargs = "+", type = "double", default = c(10, 20, 40), help = "how many taxa do we have?")
parser$add_argument("--corr", type = "double", default = 0, help = "hyperparameter controlling the off-diagonal elements of Sigma.")
parser$add_argument("--sigma", type = "double", default = 0, help = "SD of efficiencies (0 = no varying efficiency)")
parser$add_argument("--ad", type = "double", default = 0.850000, help = "adapt delta")
parser$add_argument("--mt", type = "double", default = 15, help = "max treedepth")
parser$add_argument("--num-jobs", type = "double", default = 50, help = "number of jobs run")
parser$add_argument("--taxa", type = "double", default = 10, help = "taxon of interest")
parser$add_argument("--q-obs", nargs = "+", type = "double", default = c(7), help = "number of taxa with observed qPCR")
parser$add_argument("--most-abundant", type = "double", default = 1, help = "whether or not to use most abundant taxa")
parser$add_argument("--read-data", type = "double", default = 0, help = "whether or not to read in raw output, or instead do summaries")
args <- parser$parse_args()

read_func <- function(x) tryCatch(readRDS(x), error = function(e) NA)

plots_dir <- paste0("plots/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/")
output_results_dir <- paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
if (!dir.exists(output_results_dir)) dir.create(output_results_dir, recursive = TRUE)

if (grepl("noncentered", args$stan_model[1])) {
  plot_nm <- "_noncentered"
} else {
  plot_nm <- ""
}

## --------------------------------------------------------------------------------------------------
## CREATE LIST OF OUTPUT FOR EASY SUMMARIES; GET THOSE SUMMARIES
## --------------------------------------------------------------------------------------------------
if (args$read_data) {
  output_performances <- vector(mode = "list", length = length(args$q))
  output_performances_avg_over_taxa_n <- vector(mode = "list", length = length(args$q))
  output_performances_avg_over_n <- vector(mode = "list", length = length(args$q))
  
  for (i in 1:length(args$q)) {
    ## load in results and data
    results_dir <- paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/q_", args$q[i], "/q_obs_", args$q_obs, "/")
    
    dir_mat <- expand.grid(job = 1:args$num_jobs, dir = results_dir, stan_model = args$stan_model)
    mod_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_mod_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
    data_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_data_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
    samp_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_samps_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
    
    mod_lst <- lapply(mod_nms_lst, read_func)
    data_lst <- lapply(data_nms_lst, read_func)
    # samps_lst <- lapply(samp_nms_lst, read_func)
    
    ## (1) pair model summaries with relevant data, for all taxa, all q_obs
    summary_df_1 <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type) get_summaries(w, x, y, z, type), mod_lst[1:args$num_jobs], data_lst[1:args$num_jobs], NA, MoreArgs = list(z = 1:args$q[i], type = "sigma_beta_1"), SIMPLIFY = FALSE))
    summary_df_5 <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type) get_summaries(w, x, y, z, type), mod_lst[(1:args$num_jobs) + args$num_jobs], data_lst[(1:args$num_jobs) + args$num_jobs], NA, MoreArgs = list(z = 1:args$q[i], type = "sigma_beta_5"), SIMPLIFY = FALSE))
    summary_df_10 <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type) get_summaries(w, x, y, z, type), mod_lst[(1:args$num_jobs) + 2*args$num_jobs], data_lst[(1:args$num_jobs) + 2*args$num_jobs], NA, MoreArgs = list(z = 1:args$q[i], type = "sigma_beta_10"), SIMPLIFY = FALSE))
    summary_df_50 <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type) get_summaries(w, x, y, z, type), mod_lst[(1:args$num_jobs) + 3*args$num_jobs], data_lst[(1:args$num_jobs) + 3*args$num_jobs], NA, MoreArgs = list(z = 1:args$q[i], type = "sigma_beta_50"), SIMPLIFY = FALSE))
    summary_df_100 <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type) get_summaries(w, x, y, z, type), mod_lst[(1:args$num_jobs) + 3*args$num_jobs], data_lst[(1:args$num_jobs) + 3*args$num_jobs], NA, MoreArgs = list(z = 1:args$q[i], type = "sigma_beta_100"), SIMPLIFY = FALSE))
    summary_df_144 <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type) get_summaries(w, x, y, z, type), mod_lst[(1:args$num_jobs) + 3*args$num_jobs], data_lst[(1:args$num_jobs) + 3*args$num_jobs], NA, MoreArgs = list(z = 1:args$q[i], type = "sigma_beta_144"), SIMPLIFY = FALSE))
    summary_df_196 <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type) get_summaries(w, x, y, z, type), mod_lst[(1:args$num_jobs) + 3*args$num_jobs], data_lst[(1:args$num_jobs) + 3*args$num_jobs], NA, MoreArgs = list(z = 1:args$q[i], type = "sigma_beta_196"), SIMPLIFY = FALSE))
    ## pair with the monte-carlo id, type of estimator
    summary_df_1$mc_id <- summary_df_5$mc_id <- summary_df_10$mc_id <- rep(NA, dim(summary_df_1)[1])
    summary_df_50$mc_id <- summary_df_100$mc_id <- summary_df_144$mc_id <- summary_df_196$mc_id <- rep(NA, dim(summary_df_50)[1])
    mc_lst_1 <- mc_lst_5 <- mc_lst_10 <- mc_lst_50 <- mc_lst_100 <- mc_lst_144 <- mc_lst_196 <- 1:50
    
    summary_df_1$mc_id[!is.na(summary_df_1$q)] <- rep(rep(mc_lst_1, each = args$N*args$q[i]), 1)
    summary_df_5$mc_id[!is.na(summary_df_5$q)] <- rep(rep(mc_lst_5, each = args$N*args$q[i]), 1)
    summary_df_10$mc_id[!is.na(summary_df_10$q)] <- rep(rep(mc_lst_10, each = args$N*args$q[i]), 1)
    summary_df_50$mc_id[!is.na(summary_df_50$q)] <- rep(rep(mc_lst_50, each = args$N*args$q[i]), 1)
    summary_df_100$mc_id[!is.na(summary_df_100$q)] <- rep(rep(mc_lst_100, each = args$N*args$q[i]), 1)
    summary_df_144$mc_id[!is.na(summary_df_144$q)] <- rep(rep(mc_lst_144, each = args$N*args$q[i]), 1)
    summary_df_196$mc_id[!is.na(summary_df_196$q)] <- rep(rep(mc_lst_196, each = args$N*args$q[i]), 1)
    
    ## remove naive est, change names, merge summaries
    summary_df_1_2 <- summary_df_1 %>%
      mutate(sig_1_est = est, sig_1_sd = sd, sig_1_cil = cil, sig_1_ciu = ciu,
             sig_1_wald_pred_cil = wald_pred_cil, sig_1_wald_pred_ciu = wald_pred_ciu) %>%
      select(q, q_obs, mc_id, subj_id, taxon_id, mu, qpcr, naive_est,
             sig_1_est, sig_1_sd, sig_1_cil, sig_1_ciu, sig_1_wald_pred_cil, sig_1_wald_pred_ciu,
             naive_cil, naive_ciu, naive_wald_pred_cil, naive_wald_pred_ciu)
    summary_df_5_2 <- summary_df_5 %>%
      mutate(sig_5_est = est, sig_5_sd = sd, sig_5_cil = cil, sig_5_ciu = ciu,
             sig_5_wald_pred_cil = wald_pred_cil, sig_5_wald_pred_ciu = wald_pred_ciu) %>%
      select(q, q_obs, mc_id, subj_id, taxon_id, mu, qpcr,
             sig_5_est, sig_5_sd, sig_5_cil, sig_5_ciu, sig_5_wald_pred_cil, sig_5_wald_pred_ciu)
    summary_df_10_2 <- summary_df_10 %>%
      mutate(sig_10_est = est, sig_10_sd = sd, sig_10_cil = cil, sig_10_ciu = ciu,
             sig_10_wald_pred_cil = wald_pred_cil, sig_10_wald_pred_ciu = wald_pred_ciu) %>%
      select(q, q_obs, mc_id, subj_id, taxon_id, mu, qpcr,
             sig_10_est, sig_10_sd, sig_10_cil, sig_10_ciu, sig_10_wald_pred_cil, sig_10_wald_pred_ciu)
    summary_df_50_2 <- summary_df_50 %>%
      mutate(sig_50_est = est, sig_50_sd = sd, sig_50_cil = cil, sig_50_ciu = ciu,
             sig_50_wald_pred_cil = wald_pred_cil, sig_50_wald_pred_ciu = wald_pred_ciu) %>%
      select(q, q_obs, mc_id, subj_id, taxon_id, mu, qpcr,
             sig_50_est, sig_50_sd, sig_50_cil, sig_50_ciu, sig_50_wald_pred_cil, sig_50_wald_pred_ciu)
    summary_df_100_2 <- summary_df_100 %>%
      mutate(sig_100_est = est, sig_100_sd = sd, sig_100_cil = cil, sig_100_ciu = ciu,
             sig_100_wald_pred_cil = wald_pred_cil, sig_100_wald_pred_ciu = wald_pred_ciu) %>%
      select(q, q_obs, mc_id, subj_id, taxon_id, mu, qpcr,
             sig_100_est, sig_100_sd, sig_100_cil, sig_100_ciu, sig_100_wald_pred_cil, sig_100_wald_pred_ciu)
    summary_df_144_2 <- summary_df_144 %>%
      mutate(sig_144_est = est, sig_144_sd = sd, sig_144_cil = cil, sig_144_ciu = ciu,
             sig_144_wald_pred_cil = wald_pred_cil, sig_144_wald_pred_ciu = wald_pred_ciu) %>%
      select(q, q_obs, mc_id, subj_id, taxon_id, mu, qpcr,
             sig_144_est, sig_144_sd, sig_144_cil, sig_144_ciu, sig_144_wald_pred_cil, sig_144_wald_pred_ciu)
    summary_df_196_2 <- summary_df_196 %>%
      mutate(sig_196_est = est, sig_196_sd = sd, sig_196_cil = cil, sig_196_ciu = ciu,
             sig_196_wald_pred_cil = wald_pred_cil, sig_196_wald_pred_ciu = wald_pred_ciu) %>%
      select(q, q_obs, mc_id, subj_id, taxon_id, mu, qpcr,
             sig_196_est, sig_196_sd, sig_196_cil, sig_196_ciu, sig_196_wald_pred_cil, sig_196_wald_pred_ciu)
    
    summary_df_1_5 <- left_join(summary_df_1_2, summary_df_5_2,
                                by = c("q", "q_obs","mc_id", "subj_id", "taxon_id", "mu", "qpcr"))
    summary_df_1_5_10 <- left_join(summary_df_1_5, summary_df_10_2,
                                   by = c("q", "q_obs","mc_id", "subj_id", "taxon_id", "mu", "qpcr"))
    summary_df_1_5_10_50 <- left_join(summary_df_1_5_10, summary_df_50_2,
                                      by = c("q", "q_obs","mc_id", "subj_id", "taxon_id", "mu", "qpcr"))
    summary_df_1_5_10_50_100 <- left_join(summary_df_1_5_10_50, summary_df_100_2,
                                          by = c("q", "q_obs","mc_id", "subj_id", "taxon_id", "mu", "qpcr"))
    summary_df_1_5_10_50_100_144 <- left_join(summary_df_1_5_10_50_100, summary_df_144_2,
                                              by = c("q", "q_obs","mc_id", "subj_id", "taxon_id", "mu", "qpcr"))
    summary_df <- left_join(summary_df_1_5_10_50_100_144, summary_df_196_2,
                            by = c("q", "q_obs","mc_id", "subj_id", "taxon_id", "mu", "qpcr"))
    
    ## (2) compute performance for each row
    performance_df <- summary_df %>%
      mutate(naive_mse = (mu - naive_est)^2,
             sig_1_mse = (mu - sig_1_est)^2,
             sig_5_mse = (mu - sig_5_est)^2,
             sig_10_mse = (mu - sig_10_est)^2,
             sig_50_mse = (mu - sig_50_est)^2,
             sig_100_mse = (mu - sig_100_est)^2,
             sig_144_mse = (mu - sig_144_est)^2,
             sig_196_mse = (mu - sig_196_est)^2,
             naive_mspe = (qpcr - naive_est)^2,
             sig_1_mspe = (qpcr - sig_1_est)^2,
             sig_5_mspe = (qpcr - sig_5_est)^2,
             sig_10_mspe = (qpcr - sig_10_est)^2,
             sig_50_mspe = (qpcr - sig_50_est)^2,
             sig_100_mspe = (qpcr - sig_100_est)^2,
             sig_144_mspe = (qpcr - sig_144_est)^2,
             sig_196_mspe = (qpcr - sig_196_est)^2,
             naive_cover = naive_cil <= mu & naive_ciu >= mu,
             sig_1_cover = sig_1_cil <= mu & sig_1_ciu >= mu,
             sig_5_cover = sig_5_cil <= mu & sig_5_ciu >= mu,
             sig_10_cover = sig_10_cil <= mu & sig_10_ciu >= mu,
             sig_50_cover = sig_50_cil <= mu & sig_50_ciu >= mu,
             sig_100_cover = sig_100_cil <= mu & sig_100_ciu >= mu,
             sig_144_cover = sig_144_cil <= mu & sig_144_ciu >= mu,
             sig_196_cover = sig_196_cil <= mu & sig_196_ciu >= mu,
             naive_wald_pred_cover = naive_wald_pred_cil <= qpcr & naive_wald_pred_ciu >= qpcr,
             naive_wald_width = abs(naive_wald_pred_ciu - naive_wald_pred_cil),
             sig_1_wald_pred_cover = sig_1_wald_pred_cil <= qpcr & sig_1_wald_pred_ciu >= qpcr,
             sig_5_wald_pred_cover = sig_5_wald_pred_cil <= qpcr & sig_5_wald_pred_ciu >= qpcr,
             sig_10_wald_pred_cover = sig_10_wald_pred_cil <= qpcr & sig_10_wald_pred_ciu >= qpcr,
             sig_50_wald_pred_cover = sig_50_wald_pred_cil <= qpcr & sig_50_wald_pred_ciu >= qpcr,
             sig_100_wald_pred_cover = sig_100_wald_pred_cil <= qpcr & sig_100_wald_pred_ciu >= qpcr,
             sig_144_wald_pred_cover = sig_144_wald_pred_cil <= qpcr & sig_144_wald_pred_ciu >= qpcr,
             sig_196_wald_pred_cover = sig_196_wald_pred_cil <= qpcr & sig_196_wald_pred_ciu >= qpcr,
             sig_1_wald_width = abs(sig_1_wald_pred_ciu - sig_1_wald_pred_cil),
             sig_5_wald_width = abs(sig_5_wald_pred_ciu - sig_5_wald_pred_cil),
             sig_10_wald_width = abs(sig_10_wald_pred_ciu - sig_10_wald_pred_cil),
             sig_50_wald_width = abs(sig_50_wald_pred_ciu - sig_50_wald_pred_cil),
             sig_100_wald_width = abs(sig_100_wald_pred_ciu - sig_100_wald_pred_cil),
             sig_144_wald_width = abs(sig_144_wald_pred_ciu - sig_144_wald_pred_cil),
             sig_196_wald_width = abs(sig_196_wald_pred_ciu - sig_196_wald_pred_cil),
             naive_bias = mu - naive_est,
             sig_1_bias = mu - sig_1_est,
             sig_5_bias = mu - sig_5_est,
             sig_10_bias = mu - sig_10_est,
             sig_50_bias = mu - sig_50_est,
             sig_100_bias = mu - sig_100_est,
             sig_144_bias = mu - sig_144_est,
             sig_196_bias = mu - sig_196_est)
    
    performance_df$include_taxa_in_avg <- performance_df$taxon_id >= (args$q_obs + 1)
    
    ## (3) average over MC reps for each taxon
    mc_averaged_performance <- performance_df %>%
      group_by(q, q_obs, subj_id, taxon_id) %>%
      summarize(naive_mse = mean(naive_mse, na.rm = TRUE),
                sig_1_mse = mean(sig_1_mse, na.rm = TRUE),
                sig_5_mse = mean(sig_5_mse, na.rm = TRUE),
                sig_10_mse = mean(sig_10_mse, na.rm = TRUE),
                sig_50_mse = mean(sig_50_mse, na.rm = TRUE),
                sig_100_mse = mean(sig_100_mse, na.rm = TRUE),
                sig_144_mse = mean(sig_144_mse, na.rm = TRUE),
                sig_196_mse = mean(sig_196_mse, na.rm = TRUE),
                naive_mspe = mean(naive_mspe, na.rm = TRUE),
                sig_1_mspe = mean(sig_1_mspe, na.rm = TRUE),
                sig_5_mspe = mean(sig_5_mspe, na.rm = TRUE),
                sig_10_mspe = mean(sig_10_mspe, na.rm = TRUE),
                sig_50_mspe = mean(sig_50_mspe, na.rm = TRUE),
                sig_100_mspe = mean(sig_100_mspe, na.rm = TRUE),
                sig_144_mspe = mean(sig_144_mspe, na.rm = TRUE),
                sig_196_mspe = mean(sig_196_mspe, na.rm = TRUE),
                naive_cover = mean(naive_cover, na.rm = TRUE),
                sig_1_cover = mean(sig_1_cover, na.rm = TRUE),
                sig_5_cover = mean(sig_5_cover, na.rm = TRUE),
                sig_10_cover = mean(sig_10_cover, na.rm = TRUE),
                sig_50_cover = mean(sig_50_cover, na.rm = TRUE),
                sig_100_cover = mean(sig_100_cover, na.rm = TRUE),
                sig_144_cover = mean(sig_144_cover, na.rm = TRUE),
                sig_196_cover = mean(sig_196_cover, na.rm = TRUE),
                naive_wald_pred_cover = mean(naive_wald_pred_cover, na.rm = TRUE),
                naive_wald_width = mean(naive_wald_width, na.rm = TRUE),
                sig_1_wald_pred_cover = mean(sig_1_wald_pred_cover, na.rm = TRUE),
                sig_5_wald_pred_cover = mean(sig_5_wald_pred_cover, na.rm = TRUE),
                sig_10_wald_pred_cover = mean(sig_10_wald_pred_cover, na.rm = TRUE),
                sig_50_wald_pred_cover = mean(sig_50_wald_pred_cover, na.rm = TRUE),
                sig_100_wald_pred_cover = mean(sig_100_wald_pred_cover, na.rm = TRUE),
                sig_144_wald_pred_cover = mean(sig_144_wald_pred_cover, na.rm = TRUE),
                sig_196_wald_pred_cover = mean(sig_196_wald_pred_cover, na.rm = TRUE),
                sig_1_wald_width = mean(sig_1_wald_width, na.rm = TRUE),
                sig_5_wald_width = mean(sig_5_wald_width, na.rm = TRUE),
                sig_10_wald_width = mean(sig_10_wald_width, na.rm = TRUE),
                sig_50_wald_width = mean(sig_50_wald_width, na.rm = TRUE),
                sig_100_wald_width = mean(sig_100_wald_width, na.rm = TRUE),
                sig_144_wald_width = mean(sig_144_wald_width, na.rm = TRUE),
                sig_196_wald_width = mean(sig_196_wald_width, na.rm = TRUE),
                naive_bias = mean(naive_bias, na.rm = TRUE),
                sig_1_bias = mean(sig_1_bias, na.rm = TRUE),
                sig_5_bias = mean(sig_5_bias, na.rm = TRUE),
                sig_10_bias = mean(sig_10_bias, na.rm = TRUE),
                sig_50_bias = mean(sig_50_bias, na.rm = TRUE),
                sig_100_bias = mean(sig_100_bias, na.rm = TRUE),
                sig_144_bias = mean(sig_144_bias, na.rm = TRUE),
                sig_196_bias = mean(sig_196_bias, na.rm = TRUE)) %>%
      ungroup()
    
    ## (4) set flag for taxa of interest
    mc_averaged_performance$include_taxa_in_avg <- apply(mc_averaged_performance, 1,
                                                         function(x) ifelse(!is.na(x[1]), x[4] %in% tail(1:x[1], x[1] - x[2]), NA))
    ## (5) average over taxa of interest for each q_obs
    performance_across_taxa_mu <- mc_averaged_performance %>%
      group_by(q, q_obs, subj_id) %>%
      summarize(naive_mse = mean(naive_mse),
                sig_1_mse = mean(sig_1_mse),
                sig_5_mse = mean(sig_5_mse),
                sig_10_mse = mean(sig_10_mse),
                sig_50_mse = mean(sig_50_mse),
                sig_100_mse = mean(sig_100_mse),
                sig_144_mse = mean(sig_144_mse),
                sig_196_mse = mean(sig_196_mse),
                naive_cover = mean(naive_cover, na.rm = TRUE),
                sig_1_cover = mean(sig_1_cover),
                sig_5_cover = mean(sig_5_cover),
                sig_10_cover = mean(sig_10_cover),
                sig_50_cover = mean(sig_50_cover),
                sig_100_cover = mean(sig_100_cover),
                sig_144_cover = mean(sig_144_cover),
                sig_196_cover = mean(sig_196_cover),
                naive_bias = mean(naive_bias),
                sig_1_bias = mean(sig_1_bias),
                sig_5_bias = mean(sig_5_bias),
                sig_10_bias = mean(sig_10_bias),
                sig_50_bias = mean(sig_50_bias),
                sig_100_bias = mean(sig_100_bias),
                sig_144_bias = mean(sig_144_bias),
                sig_196_bias = mean(sig_196_bias))
    performance_across_taxa_v <- mc_averaged_performance %>%
      group_by(q, q_obs, subj_id) %>%
      filter(include_taxa_in_avg) %>%
      summarize(naive_mspe = mean(naive_mspe),
                sig_1_mspe = mean(sig_1_mspe),
                sig_5_mspe = mean(sig_5_mspe),
                sig_10_mspe = mean(sig_10_mspe),
                sig_50_mspe = mean(sig_50_mspe),
                sig_100_mspe = mean(sig_100_mspe),
                sig_144_mspe = mean(sig_144_mspe),
                sig_196_mspe = mean(sig_196_mspe),
                naive_wald_pred_cover = mean(naive_wald_pred_cover, na.rm = TRUE),
                naive_wald_width = mean(naive_wald_width, na.rm = TRUE),
                sig_1_wald_pred_cover = mean(sig_1_wald_pred_cover),
                sig_5_wald_pred_cover = mean(sig_5_wald_pred_cover),
                sig_10_wald_pred_cover = mean(sig_10_wald_pred_cover),
                sig_50_wald_pred_cover = mean(sig_50_wald_pred_cover),
                sig_100_wald_pred_cover = mean(sig_100_wald_pred_cover),
                sig_144_wald_pred_cover = mean(sig_144_wald_pred_cover),
                sig_196_wald_pred_cover = mean(sig_196_wald_pred_cover),
                sig_1_wald_width = mean(sig_1_wald_width),
                sig_5_wald_width = mean(sig_5_wald_width),
                sig_10_wald_width = mean(sig_10_wald_width),
                sig_50_wald_width = mean(sig_50_wald_width),
                sig_100_wald_width = mean(sig_100_wald_width),
                sig_144_wald_width = mean(sig_144_wald_width),
                sig_196_wald_width = mean(sig_196_wald_width))
    performance_across_taxa <- left_join(performance_across_taxa_mu, performance_across_taxa_v,
                                         by = c("q", "q_obs", "subj_id"))
    ## (6) average over n, for each q_obs
    average_over_n <- performance_across_taxa %>%
      group_by(q_obs) %>%
      summarize(naive_mse = mean(naive_mse),
                sig_1_mse = mean(sig_1_mse),
                sig_5_mse = mean(sig_5_mse),
                sig_10_mse = mean(sig_10_mse),
                sig_50_mse = mean(sig_50_mse),
                sig_100_mse = mean(sig_100_mse),
                sig_144_mse = mean(sig_144_mse),
                sig_196_mse = mean(sig_196_mse),
                naive_cover = mean(naive_cover, na.rm = TRUE),
                sig_1_cover = mean(sig_1_cover),
                sig_5_cover = mean(sig_5_cover),
                sig_10_cover = mean(sig_10_cover),
                sig_50_cover = mean(sig_50_cover),
                sig_100_cover = mean(sig_100_cover),
                sig_144_cover = mean(sig_144_cover),
                sig_196_cover = mean(sig_196_cover),
                naive_mspe = mean(naive_mspe),
                sig_1_mspe = mean(sig_1_mspe),
                sig_5_mspe = mean(sig_5_mspe),
                sig_10_mspe = mean(sig_10_mspe),
                sig_50_mspe = mean(sig_50_mspe),
                sig_100_mspe = mean(sig_100_mspe),
                sig_144_mspe = mean(sig_144_mspe),
                sig_196_mspe = mean(sig_196_mspe),
                naive_wald_pred_cover = mean(naive_wald_pred_cover, na.rm = TRUE),
                naive_wald_width = mean(naive_wald_width, na.rm = TRUE),
                sig_1_wald_pred_cover = mean(sig_1_wald_pred_cover),
                sig_5_wald_pred_cover = mean(sig_5_wald_pred_cover),
                sig_10_wald_pred_cover = mean(sig_10_wald_pred_cover),
                sig_50_wald_pred_cover = mean(sig_50_wald_pred_cover),
                sig_100_wald_pred_cover = mean(sig_100_wald_pred_cover),
                sig_144_wald_pred_cover = mean(sig_144_wald_pred_cover),
                sig_196_wald_pred_cover = mean(sig_196_wald_pred_cover),
                sig_1_wald_width = mean(sig_1_wald_width),
                sig_5_wald_width = mean(sig_5_wald_width),
                sig_10_wald_width = mean(sig_10_wald_width),
                sig_50_wald_width = mean(sig_50_wald_width),
                sig_100_wald_width = mean(sig_100_wald_width),
                sig_144_wald_width = mean(sig_144_wald_width),
                sig_196_wald_width = mean(sig_196_wald_width),
                naive_bias = mean(naive_bias),
                sig_1_bias = mean(sig_1_bias),
                sig_5_bias = mean(sig_5_bias),
                sig_10_bias = mean(sig_10_bias),
                sig_50_bias = mean(sig_50_bias),
                sig_100_bias = mean(sig_100_bias),
                sig_144_bias = mean(sig_144_bias),
                sig_196_bias = mean(sig_196_bias))
    
    ## (5) add rmse, transpose
    performance_matrix <- average_over_n %>%
      mutate(naive_rmse = sqrt(naive_mse),
             sig_1_rmse = sqrt(sig_1_mse),
             sig_5_rmse = sqrt(sig_5_mse),
             sig_10_rmse = sqrt(sig_10_mse),
             sig_50_rmse = sqrt(sig_50_mse),
             sig_100_rmse = sqrt(sig_100_mse),
             sig_144_rmse = sqrt(sig_144_mse),
             sig_196_rmse = sqrt(sig_196_mse),
             naive_rmspe = sqrt(naive_mspe),
             sig_1_rmspe = sqrt(sig_1_mspe),
             sig_5_rmspe = sqrt(sig_5_mspe),
             sig_10_rmspe = sqrt(sig_10_mspe),
             sig_50_rmspe = sqrt(sig_50_mspe),
             sig_100_rmspe = sqrt(sig_100_mspe),
             sig_144_rmspe = sqrt(sig_144_mspe),
             sig_196_rmspe = sqrt(sig_196_mspe))
    q_performance <- performance_matrix %>%
      gather(key, value, -q_obs) %>%
      tidyr::extract(key, c("type", "measure"), regex = "([ns].*.[e15046])_{1}([mrcwb].*)") %>%
      tidyr::extract(type, c("estimator", "sigma_beta"), regex = "([ns].*.[eg])_{0,1}(.*)") %>%
      spread(measure, value)
    q_performance$estimator <- ifelse(q_performance$estimator == "sig", "predict_qpcr", q_performance$estimator)
    
    ## (6) average only over n, taxa
    average_over_n_taxa_mu <- performance_df %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(naive_mse = mean(naive_mse),
                sig_1_mse = mean(sig_1_mse),
                sig_5_mse = mean(sig_5_mse),
                sig_10_mse = mean(sig_10_mse),
                sig_50_mse = mean(sig_50_mse),
                sig_100_mse = mean(sig_100_mse),
                sig_144_mse = mean(sig_144_mse),
                sig_196_mse = mean(sig_196_mse),
                naive_cover = mean(naive_cover, na.rm = TRUE),
                sig_1_cover = mean(sig_1_cover),
                sig_5_cover = mean(sig_5_cover),
                sig_10_cover = mean(sig_10_cover),
                sig_50_cover = mean(sig_50_cover),
                sig_100_cover = mean(sig_100_cover),
                sig_144_cover = mean(sig_144_cover),
                sig_196_cover = mean(sig_196_cover),
                naive_bias = mean(naive_bias),
                sig_1_bias = mean(sig_1_bias),
                sig_5_bias = mean(sig_5_bias),
                sig_10_bias = mean(sig_10_bias),
                sig_50_bias = mean(sig_50_bias),
                sig_100_bias = mean(sig_100_bias),
                sig_144_bias = mean(sig_144_bias),
                sig_196_bias = mean(sig_196_bias)) %>%
      ungroup()
    average_over_n_taxa_v <- performance_df %>%
      filter(include_taxa_in_avg) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(naive_mspe = mean(naive_mspe),
                sig_1_mspe = mean(sig_1_mspe),
                sig_5_mspe = mean(sig_5_mspe),
                sig_10_mspe = mean(sig_10_mspe),
                sig_50_mspe = mean(sig_50_mspe),
                sig_100_mspe = mean(sig_100_mspe),
                sig_144_mspe = mean(sig_144_mspe),
                sig_196_mspe = mean(sig_196_mspe),
                naive_wald_pred_cover = mean(naive_wald_pred_cover, na.rm = TRUE),
                naive_wald_width = mean(naive_wald_width, na.rm = TRUE),
                sig_1_wald_pred_cover = mean(sig_1_wald_pred_cover),
                sig_5_wald_pred_cover = mean(sig_5_wald_pred_cover),
                sig_10_wald_pred_cover = mean(sig_10_wald_pred_cover),
                sig_50_wald_pred_cover = mean(sig_50_wald_pred_cover),
                sig_100_wald_pred_cover = mean(sig_100_wald_pred_cover),
                sig_144_wald_pred_cover = mean(sig_144_wald_pred_cover),
                sig_196_wald_pred_cover = mean(sig_196_wald_pred_cover),
                sig_1_wald_width = mean(sig_1_wald_width),
                sig_5_wald_width = mean(sig_5_wald_width),
                sig_10_wald_width = mean(sig_10_wald_width),
                sig_50_wald_width = mean(sig_50_wald_width),
                sig_100_wald_width = mean(sig_100_wald_width),
                sig_144_wald_width = mean(sig_144_wald_width),
                sig_196_wald_width = mean(sig_196_wald_width)) %>%
      ungroup()
    
    average_over_n_taxa <- left_join(average_over_n_taxa_mu, average_over_n_taxa_v, by = c("q", "q_obs", "mc_id")) %>%
      gather(key, value, -q, -q_obs, -mc_id) %>%
      tidyr::extract(key, c("type", "measure"), regex = "([ns].*.[e15046])_{1}([mrcwb].*)") %>%
      tidyr::extract(type, c("estimator", "sigma_beta"), regex = "([ns].*.[eg])_{0,1}(.*)") %>%
      spread(measure, value)
    average_over_n_taxa$estimator <- ifelse(average_over_n_taxa$estimator == "sig", "predict_qpcr", average_over_n_taxa$estimator)
    
    average_over_n_mu <- performance_df %>%
      group_by(q, q_obs, mc_id, taxon_id) %>%
      summarize(naive_mse = mean(naive_mse),
                sig_1_mse = mean(sig_1_mse),
                sig_5_mse = mean(sig_5_mse),
                sig_10_mse = mean(sig_10_mse),
                sig_50_mse = mean(sig_50_mse),
                sig_100_mse = mean(sig_100_mse),
                sig_144_mse = mean(sig_144_mse),
                sig_196_mse = mean(sig_196_mse),
                naive_cover = mean(naive_cover, na.rm = TRUE),
                sig_1_cover = mean(sig_1_cover),
                sig_5_cover = mean(sig_5_cover),
                sig_10_cover = mean(sig_10_cover),
                sig_50_cover = mean(sig_50_cover),
                sig_100_cover = mean(sig_100_cover),
                sig_144_cover = mean(sig_144_cover),
                sig_196_cover = mean(sig_196_cover),
                naive_bias = mean(naive_bias),
                sig_1_bias = mean(sig_1_bias),
                sig_5_bias = mean(sig_5_bias),
                sig_10_bias = mean(sig_10_bias),
                sig_50_bias = mean(sig_50_bias),
                sig_100_bias = mean(sig_100_bias),
                sig_144_bias = mean(sig_144_bias),
                sig_196_bias = mean(sig_196_bias)) %>%
      ungroup()
    average_over_n_v <- performance_df %>%
      filter(include_taxa_in_avg) %>%
      group_by(q, q_obs, mc_id, taxon_id) %>%
      summarize(naive_mspe = mean(naive_mspe),
                sig_1_mspe = mean(sig_1_mspe),
                sig_5_mspe = mean(sig_5_mspe),
                sig_10_mspe = mean(sig_10_mspe),
                sig_50_mspe = mean(sig_50_mspe),
                sig_100_mspe = mean(sig_100_mspe),
                sig_144_mspe = mean(sig_144_mspe),
                sig_196_mspe = mean(sig_196_mspe),
                naive_wald_pred_cover = mean(naive_wald_pred_cover, na.rm = TRUE),
                naive_wald_width = mean(naive_wald_width, na.rm = TRUE),
                sig_1_wald_pred_cover = mean(sig_1_wald_pred_cover),
                sig_5_wald_pred_cover = mean(sig_5_wald_pred_cover),
                sig_10_wald_pred_cover = mean(sig_10_wald_pred_cover),
                sig_50_wald_pred_cover = mean(sig_50_wald_pred_cover),
                sig_100_wald_pred_cover = mean(sig_100_wald_pred_cover),
                sig_144_wald_pred_cover = mean(sig_144_wald_pred_cover),
                sig_196_wald_pred_cover = mean(sig_196_wald_pred_cover),
                sig_1_wald_width = mean(sig_1_wald_width),
                sig_5_wald_width = mean(sig_5_wald_width),
                sig_10_wald_width = mean(sig_10_wald_width),
                sig_50_wald_width = mean(sig_50_wald_width),
                sig_100_wald_width = mean(sig_100_wald_width),
                sig_144_wald_width = mean(sig_144_wald_width),
                sig_196_wald_width = mean(sig_196_wald_width)) %>%
      ungroup()
    
    average_over_n <- left_join(average_over_n_mu, average_over_n_v, by = c("q", "q_obs", "mc_id", "taxon_id")) %>%
      gather(key, value, -q, -q_obs, -mc_id, -taxon_id) %>%
      tidyr::extract(key, c("type", "measure"), regex = "([ns].*.[e15046])_{1}([mrcwb].*)") %>%
      tidyr::extract(type, c("estimator", "sigma_beta"), regex = "([ns].*.[eg])_{0,1}(.*)") %>%
      spread(measure, value)
    average_over_n$estimator <- ifelse(average_over_n_taxa$estimator == "sig", "predict_qpcr", average_over_n_taxa$estimator)
    
    
    output_performances[[i]] <- q_performance
    output_performances_avg_over_taxa_n[[i]] <- average_over_n_taxa
    output_performances_avg_over_n[[i]] <- average_over_n
    
  }
  saveRDS(output_performances, paste0(output_results_dir, "output_performances_ab_", args$most_abundant, "_q_obs", args$q_obs, "", plot_nm, ".rds"))
  saveRDS(output_performances_avg_over_taxa_n, paste0(output_results_dir, "output_performances_avg_over_taxa_n_ab_", args$most_abundant, "_q_obs", args$q_obs, "", plot_nm, ".rds"))
  saveRDS(output_performances_avg_over_n, paste0(output_results_dir, "output_performances_avg_over_n_ab_", args$most_abundant, "_q_obs", args$q_obs, "", plot_nm, ".rds"))
} else {
  output_performances <- readRDS(paste0(output_results_dir, "output_performances_ab_", args$most_abundant, "_q_obs", args$q_obs, "", plot_nm, ".rds"))
  output_performances_avg_over_taxa_n <- readRDS(paste0(output_results_dir, "output_performances_avg_over_taxa_n_ab_", args$most_abundant, "_q_obs", args$q_obs, "", plot_nm, ".rds"))
  output_performances_avg_over_n <- readRDS(paste0(output_results_dir, "output_performances_avg_over_n_ab_", args$most_abundant, "_q_obs", args$q_obs, "", plot_nm, ".rds"))
}

## transform to lists with a matrix for each q_obs (so that it can be plotted on x-axis)
q_obs_performance <- q_to_q_obs(output_performances)

## transform into a giant matrix for the average over n and taxa
performance_avg_over_taxa_n <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n)
performance_avg_over_taxa_n$grouping <- factor(paste(performance_avg_over_taxa_n$q,
                                                     performance_avg_over_taxa_n$estimator,
                                                     sep = "_"),
                                               levels = c("10_naive", "10_predict_qpcr",
                                                          "20_naive", "20_predict_qpcr",
                                                          "40_naive", "40_predict_qpcr"),
                                               labels = c("naive; 10",
                                                          "predict_qpcr; 10",
                                                          "naive; 20",
                                                          "predict_qpcr; 20",
                                                          "naive; 40",
                                                          "predict_qpcr; 40"))

## get rmse, rmspe
performance_by_type <- performance_avg_over_taxa_n %>%
  filter(!is.na(q)) %>%
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe), numeric_sigma_beta = ifelse(!is.na(as.numeric(sigma_beta)),
                                                                           as.numeric(sigma_beta),
                                                                           0)) %>%
  select(-mse, -mspe)

## do the same thing for the average only over n
performance_avg_over_n <- do.call(rbind.data.frame, output_performances_avg_over_n)
## get grid of possible names
performance_avg_over_n$grouping <- factor(paste(performance_avg_over_n$q,
                                                performance_avg_over_n$taxon_id,
                                                performance_avg_over_taxa_n$estimator,
                                                sep = "_"))

nms_10 <- expand.grid(c("naive; 10; ",
                        "predict_qpcr; 10; "), c(1, 10, 2:9))
nms_20 <- expand.grid(c("naive; 20; ",
                        "predict_qpcr; 20; "), c(1, 10:19, 2, 20, 3:9))
nms_40 <- expand.grid(c("naive; 40; ",
                        "predict_qpcr; 40; "), c(1, 10:19, 2, 20:29, 3, 30:39, 4, 40, 5:9))
all_nms <- rbind(nms_10, nms_20, nms_40)
attr(performance_avg_over_n$grouping, "labels") <- apply(all_nms, 1, function(x) paste0(x[1], x[2]))

performance_by_type_over_n <- performance_avg_over_n %>%
  filter(!is.na(q)) %>%
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe), numeric_sigma_beta = ifelse(!is.na(as.numeric(sigma_beta)),
                                                                           as.numeric(sigma_beta),
                                                                           0)) %>%
  select(-mse, -mspe)

## --------------------------------------------------------------------------------------------------
## SET UP PLOT OPTIONS
## --------------------------------------------------------------------------------------------------
cols <- c("dodgerblue2","#E31A1C", # red
          "#6A3D9A", # purple
          "#FF7F00", # orange
          "black",
          "skyblue2","#FB9A99", # lt pink
          "palegreen2",
          "#CAB2D6", # lt purple
          "#FDBF6F", # lt orange
          "gray70", "khaki2",
          "maroon","orchid1","deeppink1","blue1","steelblue4",
          "darkturquoise","green1","yellow4","yellow3",
          "darkorange4","brown")
pchs <- c(15, 16, 17, 18, 8, 9)
fig_width <- fig_height <- 2590
point_cex <- 3
point_dodge <- 0.5
axis_cex <- 1.75

## --------------------------------------------------------------------------------------------------
## FINAL PLOT: 4 panels (prediction coverage, coverage, rmspe, mspe)
## add widths, non-abundance for supplement
## --------------------------------------------------------------------------------------------------
## create the plots
## first, get prediction interval coverage (only the proposed estimator)
pred_cover <- performance_by_type %>%
  group_by(numeric_sigma_beta, sigma_beta, grouping) %>%
  summarize(mean_cover = mean(wald_pred_cover, na.rm = TRUE)) %>%
  ggplot(aes(x = sqrt(numeric_sigma_beta), y = mean_cover,
             group = paste(as.character(grouping), "_", sigma_beta, sep = ""), shape = grouping)) +
  xlab(expression(sigma[beta])) +
  ylab(expression(Coverage)) +
  labs(shape = expression(paste("Estimator; q"))) +
  ggtitle("Prediction interval coverage") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.7, 1)) +
  guides(shape = FALSE)
## next, get credible interval coverage
cred_cover <- performance_by_type %>%
  group_by(numeric_sigma_beta, sigma_beta, grouping) %>%
  summarize(mean_cover = mean(cover)) %>%
  ggplot(aes(x = sqrt(numeric_sigma_beta), y = mean_cover,
             group = paste(as.character(grouping), "_", sigma_beta, sep = ""), shape = grouping)) +
  xlab(expression(sigma[beta])) +
  ylab(expression(Coverage)) +
  ggtitle("Credible interval coverage") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.7, 1)) +
  theme(legend.position = c(0.5, 0.3),
        legend.box = "horizontal")

## get rmspe
rmspe <- performance_by_type %>%
  group_by(numeric_sigma_beta, sigma_beta, grouping) %>%
  summarize(mean_rme = mean(rmspe[!is.infinite(rmse)])) %>%
  ggplot(aes(x = sqrt(numeric_sigma_beta), y = mean_rme,
             group = paste(as.character(grouping), "_", sigma_beta, sep = ""), shape = grouping)) +
  xlab(expression(sigma[beta])) +
  ylab(expression(RMSPE)) +
  labs(shape = expression(paste("Estimator; q"))) +
  ggtitle("Root mean squared prediction error") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 60)) +
  guides(shape = FALSE)

## get rmse
rmse <- performance_by_type %>%
  group_by(numeric_sigma_beta, sigma_beta, grouping) %>%
  summarize(mean_rme = mean(rmse[!is.infinite(rmse)])) %>%
  ggplot(aes(x = sqrt(numeric_sigma_beta), y = mean_rme,
             group = paste(as.character(grouping), "_", sigma_beta, sep = ""), shape = grouping)) +
  xlab(expression(sigma[beta])) +
  ylab(expression(RMSE)) +
  ggtitle("Root mean squared error") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 60)) +
  guides(shape = FALSE)


## actually plot them
ggsave(paste0(plots_dir, "vary_sigma_beta", plot_nm, ".png"),
       plot = plot_grid(cred_cover, pred_cover,
                        rmse, rmspe),
       width = 30, height = 25, units = "cm", dpi = 300)

## get widths, for supplement
width <- performance_by_type %>%
  filter(estimator == "predict_qpcr") %>%
  group_by(numeric_sigma_beta, sigma_beta, grouping) %>%
  summarize(mean_width = mean(wald_width[!is.infinite(wald_width)])) %>%
  ggplot(aes(x = sqrt(numeric_sigma_beta), y = mean_width,
             group = paste(as.character(grouping), "_", sigma_beta, sep = ""), shape = grouping)) +
  xlab(expression(sigma[beta])) +
  ylab(expression(Width)) +
  ggtitle("Prediction interval width") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  labs(color = "q") +
  scale_fill_manual(values = cols[c(2, 4, 6)]) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 30))

ggsave(paste0(plots_dir, "vary_sigma_beta_width", plot_nm, ".png"), plot = width,
       width = 10, height = 10, units = "cm", dpi = 300)

## plot bias, variance for supplement
bias <- performance_by_type_over_n %>%
  filter(numeric_sigma_beta == 1 | numeric_sigma_beta == 50) %>%
  group_by(q, q_obs, taxon_id, estimator, sigma_beta, numeric_sigma_beta, grouping) %>%
  summarize(bias = mean(bias[!is.infinite(bias)])^2) %>%
  ggplot(aes(x = sqrt(numeric_sigma_beta), y = bias,
             group = paste(as.character(grouping), "_", sigma_beta, sep = ""), color = grouping)) +
  xlab(expression(sigma[beta])) +
  ylab("Estimated squared bias") +
  ggtitle(expression(paste("Estimated squared bias for ", mu))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  labs(color = "q") +
  scale_color_manual(values = rep(cols, 4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 10000)) +
  theme(legend.position = c(0.9, 0.7)) +
  guides(color = FALSE)

ggsave(paste0(plots_dir, "vary_sigma_beta_bias", plot_nm, ".png"), plot = bias,
       width = 10, height = 10, units = "cm", dpi = 300)

## -----------------------------------------------------------------------------------------------------
## plot of coverage vs rank order
## NB: since I'm using the most abundant as qPCR, and the rest are ordered already,
## taxon_id is equivalent to rank order
## -----------------------------------------------------------------------------------------------------

## ------------------------------------------------------------
## q = 40
## ------------------------------------------------------------
cover_vs_rank_order_q_40_sig_beta_1 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 40, estimator == "predict_qpcr", sigma_beta == "1") %>%
  group_by(taxon_id) %>%
  summarize(cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle(expression(paste(sigma[beta], " = ", sqrt(1))))

cover_vs_rank_order_q_40_sig_beta_5 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 40, estimator == "predict_qpcr", sigma_beta == "5") %>%
  group_by(taxon_id) %>%
  summarize(cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle(expression(paste(sigma[beta], " = ", sqrt(5))))

cover_vs_rank_order_q_40_sig_beta_10 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 40, estimator == "predict_qpcr", sigma_beta == "10") %>%
  group_by(taxon_id) %>%
  summarize(cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle(expression(paste(sigma[beta], " = ", sqrt(10))))

cover_vs_rank_order_q_40_sig_beta_50 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 40, estimator == "predict_qpcr", sigma_beta == "50") %>%
  group_by(taxon_id) %>%
  summarize(cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle(expression(paste(sigma[beta], " = ", sqrt(50))))


ggsave(paste0(plots_dir, "vary_sigma_beta_cover_vs_rank_order_40", plot_nm, ".png"),
       plot = plot_grid(cover_vs_rank_order_q_40_sig_beta_1, cover_vs_rank_order_q_40_sig_beta_5,
                        cover_vs_rank_order_q_40_sig_beta_10, cover_vs_rank_order_q_40_sig_beta_50),
       width = 30, height = 25, units = "cm", dpi = 300)

## ------------------------------------------------------------
## q = 20
## ------------------------------------------------------------
cover_vs_rank_order_q_20_sig_beta_1 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 20, estimator == "predict_qpcr", sigma_beta == "1") %>%
  group_by(taxon_id) %>%
  summarize(cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle(expression(paste("Coverage of credible intervals, ", sigma[beta], " = ", sqrt(1))))

cover_vs_rank_order_q_20_sig_beta_5 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 20, estimator == "predict_qpcr", sigma_beta == "5") %>%
  group_by(taxon_id) %>%
  summarize(cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle(expression(paste("Coverage of credible intervals, ", sigma[beta], " = ", sqrt(5))))

cover_vs_rank_order_q_20_sig_beta_10 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 20, estimator == "predict_qpcr", sigma_beta == "10") %>%
  group_by(taxon_id) %>%
  summarize(cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle(expression(paste("Coverage of credible intervals, ", sigma[beta], " = ", sqrt(10))))

cover_vs_rank_order_q_20_sig_beta_50 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 20, estimator == "predict_qpcr", sigma_beta == "50") %>%
  group_by(taxon_id) %>%
  summarize(cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle(expression(paste("Coverage of credible intervals, ", sigma[beta], " = ", sqrt(50))))


ggsave(paste0(plots_dir, "vary_sigma_beta_cover_vs_rank_order_20", plot_nm, ".png"),
       plot = plot_grid(cover_vs_rank_order_q_20_sig_beta_1, cover_vs_rank_order_q_20_sig_beta_5,
                        cover_vs_rank_order_q_20_sig_beta_10, cover_vs_rank_order_q_20_sig_beta_50),
       width = 30, height = 25, units = "cm", dpi = 300)

## ------------------------------------------------------------
## q = 10
## ------------------------------------------------------------
cover_vs_rank_order_q_10_sig_beta_1 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 10, estimator == "predict_qpcr", sigma_beta == "1") %>%
  group_by(taxon_id) %>%
  summarize(cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle(expression(paste("Coverage of credible intervals, ", sigma[beta], " = ", sqrt(1))))

cover_vs_rank_order_q_10_sig_beta_5 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 10, estimator == "predict_qpcr", sigma_beta == "5") %>%
  group_by(taxon_id) %>%
  summarize(cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle(expression(paste("Coverage of credible intervals, ", sigma[beta], " = ", sqrt(5))))

cover_vs_rank_order_q_10_sig_beta_10 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 10, estimator == "predict_qpcr", sigma_beta == "10") %>%
  group_by(taxon_id) %>%
  summarize(cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle(expression(paste("Coverage of credible intervals, ", sigma[beta], " = ", sqrt(10))))

cover_vs_rank_order_q_10_sig_beta_50 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 10, estimator == "predict_qpcr", sigma_beta == "50") %>%
  group_by(taxon_id) %>%
  summarize(cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle(expression(paste("Coverage of credible intervals, ", sigma[beta], " = ", sqrt(50))))


ggsave(paste0(plots_dir, "vary_sigma_beta_cover_vs_rank_order_10", plot_nm, ".png"),
       plot = plot_grid(cover_vs_rank_order_q_10_sig_beta_1, cover_vs_rank_order_q_10_sig_beta_5,
                        cover_vs_rank_order_q_10_sig_beta_10, cover_vs_rank_order_q_10_sig_beta_50),
       width = 30, height = 25, units = "cm", dpi = 300)
