#!/usr/local/bin/Rscript
##################################################################################
## FILE: load_sigma_beta.R
##
## CREATED: 06 August 2019 by Brian Williamson
##
## PURPOSE: load in results from a qPCR simulation with varying renormalization of efficiencies
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
parser$add_argument("--sim-name", default = "e_constraint", help = "name of the simulation")
parser$add_argument("--stan-model", nargs = "+", default = c("predict_qpcr_with_varying_efficiency_noncentered", 
                                                             "predict_qpcr_with_varying_efficiency_noncentered_hard_e_constraint"),
                    help = "Which Stan model file to use.")
parser$add_argument("--N", type = "double", default = 50, help = "sample size")
parser$add_argument("--q", nargs = "+", type = "double", default = 20, help = "how many taxa do we have?")
parser$add_argument("--corr", type = "double", default = 0, help = "hyperparameter controlling the off-diagonal elements of Sigma.")
parser$add_argument("--sigma", type = "double", default = c(0, 0.5, 1), help = "SD of efficiencies (0 = no varying efficiency)")
parser$add_argument("--ad", type = "double", default = 0.850000, help = "adapt delta")
parser$add_argument("--mt", type = "double", default = 15, help = "max treedepth")
parser$add_argument("--num-jobs", type = "double", default = 50, help = "number of jobs run")
parser$add_argument("--taxa", type = "double", default = 10, help = "taxon of interest")
parser$add_argument("--q-obs", nargs = "+", type = "double", default = c(7), help = "number of taxa with observed qPCR")
parser$add_argument("--most-abundant", type = "double", default = 1, help = "whether or not to use most abundant taxa")
parser$add_argument("--read-data", type = "double", default = 0, help = "whether or not to read in raw output, or instead do summaries")
args <- parser$parse_args()

read_func <- function(x) tryCatch(readRDS(x), error = function(e) NA)

plots_dir <- paste0("plots/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/")
output_results_dir <- paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
if (!dir.exists(output_results_dir)) dir.create(output_results_dir, recursive = TRUE)

if (grepl("noncentered", args$stan_model[1])) {
  plot_nm <- "_noncentered"
} else {
  plot_nm <- ""
}

if (args$read_data) {
  output_performances <- vector(mode = "list", length = length(args$sigma))

  output_performances_avg_over_taxa_n <- output_performances_avg_over_taxa_n_mn <- output_performances_avg_over_taxa_n_var <- output_performances_avg_over_n <- vector(mode = "list", length = length(args$sigma))
  output_performances_e <- vector(mode = "list", length = length(args$sigma))
  output_performances_e_avg_over_taxa <- vector(mode = "list", length = length(args$sigma))
  output_log_e_averages <- vector(mod = "list", length = length(args$sigma))
  
  for (i in 1:length(args$sigma)) {
    results_dir <- paste0("results/", args$sim_name, "/ve_", args$sigma[i], "/cov_", args$corr, "/n_", args$N, "/q_", args$q, "/q_obs_", args$q_obs, "/")
    
    dir_mat <- expand.grid(job = 1:args$num_jobs, dir = results_dir, stan_model = args$stan_model)
    mod_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_mod_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
    data_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_data_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
    samp_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_samps_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
    
    mod_lst <- lapply(mod_nms_lst, read_func)
    data_lst <- lapply(data_nms_lst, read_func)
    samp_lst <- lapply(samp_nms_lst, read_func)
    
    ## (0) get mean of the log-transformed efficiency samples
    log_e_averages_lst <- lapply(samp_lst, function(x) colMeans(log(x$e)))
    log_e_averages <- do.call(rbind.data.frame, log_e_averages_lst)
    colnames(log_e_averages) <- paste0("taxon_", 1:args$q)
    log_e_averages$mc_id <- rep(1:50, 2)
    log_e_averages$estimator <- c(rep("soft_constraint", 50), rep("e_normalized", 50))
    output_log_e_averages[[i]] <- log_e_averages
    
    ## (1) pair model summaries with relevant data, for all taxa, all q_obs
    summary_df_no_renorm <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type) get_summaries(w, x, y, z, type), mod_lst[1:args$num_jobs], data_lst[1:args$num_jobs], NA, MoreArgs = list(z = 1:args$q, type = "no_renorm"), SIMPLIFY = FALSE))
    summary_df_renorm <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type) get_summaries(w, x, y, z, type), mod_lst[(args$num_jobs + 1):(2*args$num_jobs)], data_lst[(args$num_jobs + 1):(2*args$num_jobs)], NA, MoreArgs = list(z = 1:args$q, type = "renorm"), SIMPLIFY = FALSE))
    ## pair with the monte-carlo id, type of estimator
    summary_df_no_renorm$mc_id <- rep(NA, dim(summary_df_no_renorm)[1])
    summary_df_renorm$mc_id <- rep(NA, dim(summary_df_renorm)[1])
    no_renorm_mc_lst <- 1:50
    renorm_mc_lst <- 1:50
    summary_df_no_renorm$mc_id[!is.na(summary_df_no_renorm$q)] <- rep(rep(no_renorm_mc_lst, each = args$N*args$q), 1)
    summary_df_renorm$mc_id[!is.na(summary_df_renorm$q)] <- rep(rep(renorm_mc_lst, each = args$N*args$q), 1)
    
    ## (1a) pair model summaries with data for efficiency
    e_summary_df_no_renorm <- do.call(rbind.data.frame, mapply(function(x, y, z) get_summaries_with_ve(x, y, z), mod_lst[1:args$num_jobs], data_lst[1:args$num_jobs], MoreArgs = list(z = 1:args$q), SIMPLIFY = FALSE))
    e_summary_df_renorm <- do.call(rbind.data.frame, mapply(function(x, y, z) get_summaries_with_ve(x, y, z), mod_lst[(args$num_jobs + 1):(2*args$num_jobs)], data_lst[(args$num_jobs + 1):(2*args$num_jobs)], MoreArgs = list(z = 1:args$q), SIMPLIFY = FALSE))
    e_summary_df_no_renorm$mc_id <- rep(NA, dim(e_summary_df_renorm)[1])
    e_summary_df_renorm$mc_id <- rep(NA, dim(e_summary_df_renorm)[1])
    e_summary_df_no_renorm$mc_id[!is.na(e_summary_df_no_renorm$q)] <- rep(1:50, each = args$q)
    e_summary_df_renorm$mc_id[!is.na(e_summary_df_renorm$q)] <- rep(1:50, each = args$q)
    e_summary_df_renorm_2 <- e_summary_df_renorm %>% 
      mutate(ve_est_e = est_e, ve_est_sigma = est_sigma, ve_e_cil = e_cil,
             ve_e_ciu = e_ciu, ve_sigma_cil = sigma_cil, ve_sigma_ciu = sigma_ciu) %>% 
      select(-est_e, -est_sigma, -e_cil, -e_ciu, -sigma_cil, -sigma_ciu)
    e_summary_df <- left_join(e_summary_df_no_renorm, e_summary_df_renorm_2, by = c("q", "q_obs", "taxon_id", "mc_id", "sigma", "e"))
    
    ## remove naive est from ve summaries, change names, and merge with no_renorm summaries
    summary_df_renorm_2 <- summary_df_renorm %>% 
      select(-naive_est) %>% 
      mutate(ve_est = est, ve_sd = sd, ve_cil = cil, ve_ciu = ciu, ve_wald_cil = wald_cred_cil, ve_wald_ciu = wald_cred_ciu,
             ve_wald_pred_cil = wald_pred_cil, ve_wald_pred_ciu = wald_pred_ciu,
             ve_quantile_cred_pred_cil = quantile_cred_pred_cil, ve_quantile_cred_pred_ciu = quantile_cred_pred_ciu,
             ve_quantile_samp_pred_cil = quantile_samp_pred_cil, ve_quantile_samp_pred_ciu = quantile_samp_pred_ciu) %>% 
      select(q, q_obs, subj_id, taxon_id, mc_id, mu, qpcr, ve_est, ve_sd, ve_cil, ve_ciu, ve_wald_cil, ve_wald_ciu,
             ve_wald_pred_cil, ve_wald_pred_ciu, ve_quantile_cred_pred_cil, ve_quantile_cred_pred_ciu,
             ve_quantile_samp_pred_cil, ve_quantile_samp_pred_ciu)
    summary_df_no_renorm_2 <- summary_df_no_renorm %>% 
      mutate(wald_cil = wald_cred_cil, wald_ciu = wald_cred_ciu) %>% 
      select(-type)
    summary_df <- left_join(summary_df_no_renorm_2, summary_df_renorm_2, by = c("q", "q_obs", "subj_id", "taxon_id", "mc_id", "mu", "qpcr"))
    
    ## (2) compute performance for each row
    performance_df <- summary_df %>%
      mutate(mse = (mu - est)^2, ve_mse = (mu - ve_est)^2, naive_mse = (mu - naive_est)^2,
             mspe = (qpcr - est)^2, ve_mspe = (qpcr - ve_est)^2, naive_mspe = (qpcr - naive_est)^2,
             cover = cil <= mu & ciu >= mu,
             wald_cover = wald_cil <= mu & wald_ciu >= mu,
             ve_cover = ve_cil <= mu & ve_ciu >= mu,
             ve_wald_cover = ve_wald_cil <= mu & ve_wald_ciu >= mu,
             naive_cover = naive_cil <= mu & naive_ciu >= mu,
             wald_pred_cover = wald_pred_cil <= qpcr & wald_pred_ciu >= qpcr,
             ve_wald_pred_cover = ve_wald_pred_cil <= qpcr & ve_wald_pred_ciu >= qpcr,
             quantile_cred_pred_cover = quantile_cred_pred_cil <= qpcr & quantile_cred_pred_ciu >= qpcr,
             ve_quantile_cred_pred_cover = ve_quantile_cred_pred_cil <= qpcr & ve_quantile_cred_pred_ciu >= qpcr,
             quantile_samp_pred_cover = quantile_samp_pred_cil <= qpcr & quantile_samp_pred_ciu >= qpcr,
             ve_quantile_samp_pred_cover = ve_quantile_samp_pred_cil <= qpcr & ve_quantile_samp_pred_ciu >= qpcr,
             naive_pred_cover = naive_wald_pred_cil <= qpcr & naive_wald_pred_ciu >= qpcr,
             wald_width = abs(wald_pred_ciu - wald_pred_cil),
             naive_width = abs(naive_wald_pred_ciu - naive_wald_pred_cil),
             quantile_cred_width = abs(quantile_cred_pred_ciu - quantile_cred_pred_cil),
             quantile_samp_width = abs(quantile_samp_pred_ciu - quantile_samp_pred_cil),
             ve_wald_width = abs(ve_wald_pred_ciu - ve_wald_pred_cil),
             ve_quantile_cred_width = abs(ve_quantile_cred_pred_ciu - ve_quantile_cred_pred_cil),
             ve_quantile_samp_width = abs(ve_quantile_samp_pred_ciu - ve_quantile_samp_pred_cil),
             bias = mu - est, ve_bias = mu - ve_est, naive_bias = mu - naive_est,
             sigma = args$sigma[i])
    
    performance_df$include_taxa_in_v_avg <- performance_df$taxon_id >= (args$q_obs + 1)
    
    ## (2a) compute mses, coverages for each row
    e_performance_df <- e_summary_df %>% 
      mutate(mse = (e - est_e)^2, ve_mse = (e - ve_est_e)^2,
             sigma_mse = (sigma - est_sigma)^2, ve_sigma_mse = (sigma - ve_est_sigma)^2,
             cover = e_cil <= e & e_ciu >= e, ve_cover = ve_e_cil <= e & ve_e_ciu >= e,
             sigma_cover = sigma_cil <= sigma & sigma_ciu >= sigma, ve_sigma_cover = ve_sigma_cil <= sigma & ve_sigma_ciu >= sigma,
             bias = e - est_e, sigma_bias = sigma - est_sigma, ve_bias = e - ve_est_e, ve_sigma_bias = sigma - ve_est_sigma,
             sigma = args$sigma[i])
    
    ## (3) average over MC reps for each taxon
    mc_averaged_performance <- performance_df %>%
      filter(include_taxa_in_avg) %>%      
      select(q, q_obs, subj_id, taxon_id, mse, ve_mse, naive_mse, mspe, ve_mspe, naive_mspe, 
             cover, ve_cover, naive_cover, wald_cover, ve_wald_cover,
             wald_pred_cover, ve_wald_pred_cover, naive_pred_cover,
             quantile_cred_pred_cover, ve_quantile_cred_pred_cover,
             quantile_samp_pred_cover, ve_quantile_samp_pred_cover,
             wald_width, naive_width, ve_wald_width, quantile_cred_width, ve_quantile_cred_width,
             quantile_samp_width, ve_quantile_samp_width,
             bias, ve_bias, naive_bias, sigma) %>%
      group_by(q, q_obs, subj_id, taxon_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), ve_mse = mean(ve_mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE), 
                mspe = mean(mspe, na.rm = TRUE), ve_mspe = mean(ve_mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE), 
                cover = mean(cover, na.rm = TRUE), ve_cover = mean(ve_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE), ve_wald_cover = mean(ve_wald_cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE),
                wald_pred_cover = mean(wald_pred_cover, na.rm = TRUE), ve_wald_pred_cover = mean(ve_wald_pred_cover, na.rm = TRUE),
                quantile_cred_pred_cover = mean(quantile_cred_pred_cover, na.rm = TRUE), ve_quantile_cred_pred_cover = mean(ve_quantile_cred_pred_cover, na.rm = TRUE),
                quantile_samp_pred_cover = mean(quantile_samp_pred_cover, na.rm = TRUE), ve_quantile_samp_pred_cover = mean(ve_quantile_samp_pred_cover, na.rm = TRUE),
                naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                naive_width = mean(naive_width, na.rm = TRUE),
                wald_width = mean(wald_width, na.rm = TRUE), ve_wald_width = mean(ve_wald_width, na.rm = TRUE),
                quantile_cred_width = mean(quantile_cred_width, na.rm = TRUE), ve_quantile_cred_width = mean(ve_quantile_cred_width, na.rm = TRUE),
                quantile_samp_width = mean(quantile_samp_width, na.rm = TRUE), ve_quantile_samp_width = mean(ve_quantile_samp_width, na.rm = TRUE),
                bias = mean(bias, na.rm = TRUE), ve_bias = mean(ve_bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE),
                sigma = mean(sigma)) %>%
      ungroup()
    
    ## (3a) average over MC reps for each taxon, for e
    e_mc_averaged_performance <- e_performance_df %>% 
      group_by(q, q_obs, taxon_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), sigma_mse = mean(sigma_mse, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), sigma_cover = mean(sigma_cover, na.rm = TRUE),
                bias = mean(bias, na.rm = TRUE), sigma_bias = mean(sigma_bias, na.rm = TRUE),
                ve_mse = mean(ve_mse, na.rm = TRUE), ve_sigma_mse = mean(ve_sigma_mse, na.rm = TRUE),
                ve_cover = mean(ve_cover, na.rm = TRUE), ve_sigma_cover = mean(ve_sigma_cover, na.rm = TRUE),
                ve_bias = mean(ve_bias, na.rm = TRUE), ve_sigma_bias = mean(ve_sigma_bias, na.rm = TRUE),
                sigma = mean(sigma)) %>%
      ungroup()
    
    ## (4) set flag for taxa of interest
    mc_averaged_performance$include_taxa_in_v_avg <- apply(mc_averaged_performance, 1, 
                                                           function(x) ifelse(!is.na(x[1]), x[4] %in% tail(1:x[1], x[1] - x[2]), NA))
    ## (5) average over taxa of interest for each q_obs
    performance_across_taxa_mu <- mc_averaged_performance %>%
      group_by(q, q_obs, subj_id) %>%
      summarize(mse = mean(mse), ve_mse = mean(ve_mse), naive_mse = mean(naive_mse),
                cover = mean(cover), ve_cover = mean(ve_cover), naive_cover = mean(naive_cover), wald_cover = mean(wald_cover), ve_wald_cover = mean(ve_wald_cover),
                bias = mean(bias), ve_bias = mean(ve_bias), naive_bias = mean(naive_bias),
                sigma = mean(sigma))
    performance_across_taxa_v <- mc_averaged_performance %>% 
      group_by(q, q_obs, subj_id) %>% 
      filter(include_taxa_in_v_avg) %>%
      summarize(mspe = mean(mspe), ve_mspe = mean(ve_mspe), naive_mspe = mean(naive_mspe),
                wald_pred_cover = mean(wald_pred_cover), ve_wald_pred_cover = mean(ve_wald_pred_cover), naive_pred_cover = mean(naive_pred_cover),
                quantile_cred_pred_cover = mean(quantile_cred_pred_cover), ve_quantile_cred_pred_cover = mean(ve_quantile_cred_pred_cover),
                quantile_samp_pred_cover = mean(quantile_samp_pred_cover), ve_quantile_samp_pred_cover = mean(ve_quantile_samp_pred_cover),
                naive_width = mean(naive_width), wald_width = mean(wald_width), ve_wald_width = mean(ve_wald_width),
                quantile_cred_width = mean(quantile_cred_width), ve_quantile_cred_width = mean(ve_quantile_cred_width),
                quantile_samp_width = mean(quantile_samp_width), ve_quantile_samp_width = mean(ve_quantile_samp_width))
    performance_across_taxa <- left_join(performance_across_taxa_mu, performance_across_taxa_v,
                                         by = c("q", "q_obs", "subj_id"))
    ## (6) average over n, for each q_obs
    average_over_n <- performance_across_taxa %>%
      group_by(q_obs) %>%
      summarize(mse = mean(mse), ve_mse = mean(ve_mse), naive_mse = mean(naive_mse),
                mspe = mean(mspe), ve_mspe = mean(ve_mspe), naive_mspe = mean(naive_mspe),
                cover = mean(cover), ve_cover = mean(ve_cover), naive_cover = mean(naive_cover), wald_cover = mean(wald_cover), ve_wald_cover = mean(ve_wald_cover),
                wald_pred_cover = mean(wald_pred_cover), ve_wald_pred_cover = mean(ve_wald_pred_cover),
                quantile_cred_pred_cover = mean(quantile_cred_pred_cover), ve_quantile_cred_pred_cover = mean(ve_quantile_cred_pred_cover),
                quantile_samp_pred_cover = mean(quantile_samp_pred_cover), ve_quantile_samp_pred_cover = mean(ve_quantile_samp_pred_cover),
                naive_width = mean(naive_width), wald_width = mean(wald_width), ve_wald_width = mean(ve_wald_width),
                quantile_cred_width = mean(quantile_cred_width), ve_quantile_cred_width = mean(ve_quantile_cred_width),
                quantile_samp_width = mean(quantile_samp_width), ve_quantile_samp_width = mean(ve_quantile_samp_width),
                bias = mean(bias), ve_bias = mean(ve_bias), naive_bias = mean(naive_bias), sigma = mean(sigma))
    average_over_n_single_taxon <- mc_averaged_performance %>%
      filter(taxon_id == 10) %>%
      group_by(q_obs) %>%
      summarize(mse = mean(mse), ve_mse = mean(ve_mse), naive_mse = mean(naive_mse),
                mspe = mean(mspe), ve_mspe = mean(ve_mspe), naive_mspe = mean(naive_mspe),
                cover = mean(cover), ve_cover = mean(ve_cover), naive_cover = mean(naive_cover), wald_cover = mean(wald_cover), ve_wald_cover = mean(ve_wald_cover),
                wald_pred_cover = mean(wald_pred_cover), ve_wald_pred_cover = mean(ve_wald_pred_cover),
                quantile_cred_pred_cover = mean(quantile_cred_pred_cover), ve_quantile_cred_pred_cover = mean(ve_quantile_cred_pred_cover),
                quantile_samp_pred_cover = mean(quantile_samp_pred_cover), ve_quantile_samp_pred_cover = mean(ve_quantile_samp_pred_cover),
                naive_width = mean(naive_width), wald_width = mean(wald_width), ve_wald_width = mean(ve_wald_width),
                quantile_cred_width = mean(quantile_cred_width), ve_quantile_cred_width = mean(ve_quantile_cred_width),
                quantile_samp_width = mean(quantile_samp_width), ve_quantile_samp_width = mean(ve_quantile_samp_width),
                bias = mean(bias), ve_bias = mean(ve_bias), naive_bias = mean(naive_bias), sigma = mean(sigma))
    
    ## (5) add rmse, transpose
    performance_matrix <- average_over_n %>%
      mutate(rmse = sqrt(mse), ve_mse = sqrt(ve_mse), naive_rmse = sqrt(naive_mse),
             rmspe = sqrt(mspe), ve_rmspe = sqrt(ve_mspe), naive_rmspe = sqrt(naive_mspe))
    performance_matrix_single_taxon <- average_over_n_single_taxon %>%
      mutate(rmse = sqrt(mse), ve_mse = sqrt(ve_mse), naive_rmse = sqrt(naive_mse),
             rmspe = sqrt(mspe), ve_rmspe = sqrt(ve_mspe), naive_rmspe = sqrt(naive_mspe))
    
    ## (6) average only over n, taxa
    ## (a) filter by colSums(W) > 0
    average_over_n_taxa_mu <- performance_df %>% 
      filter(include_taxa_in_avg) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, ve_mse, naive_mse, mspe, ve_mspe, naive_mspe, 
             cover, ve_cover, naive_cover, wald_cover, ve_wald_cover, wald_pred_cover, ve_wald_pred_cover, quantile_cred_pred_cover, ve_quantile_cred_pred_cover, naive_pred_cover,
             quantile_samp_pred_cover, ve_quantile_samp_pred_cover, wald_width, ve_wald_width, quantile_cred_width, ve_quantile_cred_width,
             quantile_samp_width, ve_quantile_samp_width,
             bias, ve_bias, naive_bias, sigma) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), ve_mse = mean(ve_mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), ve_cover = mean(ve_cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE), ve_wald_cover = mean(ve_wald_cover, na.rm = TRUE),
                bias = mean(bias, na.rm = TRUE), ve_bias = mean(ve_bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE), sigma = mean(sigma)) %>% 
      ungroup()
    average_over_n_taxa_v <- performance_df %>%
      filter(include_taxa_in_v_avg) %>% 
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, ve_mse, naive_mse, mspe, ve_mspe, naive_mspe, 
             cover, ve_cover, wald_pred_cover, ve_wald_pred_cover, quantile_cred_pred_cover, ve_quantile_cred_pred_cover, naive_pred_cover,
             quantile_samp_pred_cover, ve_quantile_samp_pred_cover, wald_width, ve_wald_width, quantile_cred_width, ve_quantile_cred_width,
             naive_width, quantile_samp_width, ve_quantile_samp_width,
             bias, ve_bias, naive_bias, sigma) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mspe = mean(mspe, na.rm = TRUE), ve_mspe = mean(ve_mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE), 
                wald_pred_cover = mean(wald_pred_cover, na.rm = TRUE), ve_wald_pred_cover = mean(ve_wald_pred_cover, na.rm = TRUE), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                quantile_cred_pred_cover = mean(quantile_cred_pred_cover, na.rm = TRUE), ve_quantile_cred_pred_cover = mean(ve_quantile_cred_pred_cover, na.rm = TRUE),
                quantile_samp_pred_cover = mean(quantile_samp_pred_cover, na.rm = TRUE), ve_quantile_samp_pred_cover = mean(ve_quantile_samp_pred_cover, na.rm = TRUE),
                naive_width = mean(naive_width, na.rm = TRUE), wald_width = mean(wald_width, na.rm = TRUE), ve_wald_width = mean(ve_wald_width, na.rm = TRUE),
                quantile_cred_width = mean(quantile_cred_width, na.rm = TRUE), ve_quantile_cred_width = mean(ve_quantile_cred_width, na.rm = TRUE),
                quantile_samp_width = mean(quantile_samp_width, na.rm = TRUE), ve_quantile_samp_width = mean(ve_quantile_samp_width, na.rm = TRUE)) %>% 
      ungroup()
    
    average_over_n_taxa <- left_join(average_over_n_taxa_mu, average_over_n_taxa_v, by = c("q", "q_obs", "mc_id"))
    
    ## (b) filter by colMeans(W) > 0.5
    average_over_n_taxa_mu_mn <- performance_df %>% 
      filter(include_taxa_in_avg_w_mn) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, ve_mse, naive_mse, mspe, ve_mspe, naive_mspe, 
             cover, ve_cover, naive_cover, wald_cover, ve_wald_cover, wald_pred_cover, ve_wald_pred_cover, quantile_cred_pred_cover, ve_quantile_cred_pred_cover, naive_pred_cover,
             quantile_samp_pred_cover, ve_quantile_samp_pred_cover, naive_width, wald_width, ve_wald_width, quantile_cred_width, ve_quantile_cred_width,
             quantile_samp_width, ve_quantile_samp_width,
             bias, ve_bias, naive_bias, sigma) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), ve_mse = mean(ve_mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), ve_cover = mean(ve_cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE), ve_wald_cover = mean(ve_wald_cover, na.rm = TRUE),
                bias = mean(bias, na.rm = TRUE), ve_bias = mean(ve_bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE), sigma = mean(sigma)) %>% 
      ungroup()
    average_over_n_taxa_v_mn <- performance_df %>%
      filter(include_taxa_in_v_avg, include_taxa_in_avg_w_mn) %>% 
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, ve_mse, naive_mse, mspe, ve_mspe, naive_mspe, 
             cover, ve_cover, wald_pred_cover, ve_wald_pred_cover, quantile_cred_pred_cover, ve_quantile_cred_pred_cover, naive_pred_cover,
             quantile_samp_pred_cover, ve_quantile_samp_pred_cover, naive_width, wald_width, ve_wald_width, quantile_cred_width, ve_quantile_cred_width,
             quantile_samp_width, ve_quantile_samp_width,
             bias, ve_bias, naive_bias, sigma) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mspe = mean(mspe, na.rm = TRUE), ve_mspe = mean(ve_mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE), 
                wald_pred_cover = mean(wald_pred_cover, na.rm = TRUE), ve_wald_pred_cover = mean(ve_wald_pred_cover, na.rm = TRUE), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                quantile_cred_pred_cover = mean(quantile_cred_pred_cover, na.rm = TRUE), ve_quantile_cred_pred_cover = mean(ve_quantile_cred_pred_cover, na.rm = TRUE),
                quantile_samp_pred_cover = mean(quantile_samp_pred_cover, na.rm = TRUE), ve_quantile_samp_pred_cover = mean(ve_quantile_samp_pred_cover, na.rm = TRUE),
                naive_width = mean(naive_width, na.rm = TRUE), wald_width = mean(wald_width, na.rm = TRUE), ve_wald_width = mean(ve_wald_width, na.rm = TRUE),
                quantile_cred_width = mean(quantile_cred_width, na.rm = TRUE), ve_quantile_cred_width = mean(ve_quantile_cred_width, na.rm = TRUE),
                quantile_samp_width = mean(quantile_samp_width, na.rm = TRUE), ve_quantile_samp_width = mean(ve_quantile_samp_width, na.rm = TRUE)) %>% 
      ungroup()
    
    average_over_n_taxa_mn <- left_join(average_over_n_taxa_mu_mn, average_over_n_taxa_v_mn, by = c("q", "q_obs", "mc_id"))
    
    ## (c) filter by colVars(W) > 1
    average_over_n_taxa_mu_var <- performance_df %>% 
      filter(include_taxa_in_avg_var) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, ve_mse, naive_mse, mspe, ve_mspe, naive_mspe, 
             cover, ve_cover, naive_cover, wald_cover, ve_wald_cover, wald_pred_cover, ve_wald_pred_cover, quantile_cred_pred_cover, ve_quantile_cred_pred_cover, naive_pred_cover,
             quantile_samp_pred_cover, ve_quantile_samp_pred_cover, naive_width, wald_width, ve_wald_width, quantile_cred_width, ve_quantile_cred_width,
             quantile_samp_width, ve_quantile_samp_width,
             bias, ve_bias, naive_bias, sigma) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), ve_mse = mean(ve_mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), ve_cover = mean(ve_cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE), ve_wald_cover = mean(ve_wald_cover, na.rm = TRUE),
                bias = mean(bias, na.rm = TRUE), ve_bias = mean(ve_bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE), sigma = mean(sigma)) %>% 
      ungroup()
    average_over_n_taxa_v_var <- performance_df %>%
      filter(include_taxa_in_v_avg, include_taxa_in_avg_var) %>% 
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, ve_mse, naive_mse, mspe, ve_mspe, naive_mspe, 
             cover, ve_cover, wald_pred_cover, ve_wald_pred_cover, quantile_cred_pred_cover, ve_quantile_cred_pred_cover, naive_pred_cover,
             quantile_samp_pred_cover, ve_quantile_samp_pred_cover, naive_width, wald_width, ve_wald_width, quantile_cred_width, ve_quantile_cred_width,
             quantile_samp_width, ve_quantile_samp_width,
             bias, ve_bias, naive_bias, sigma) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mspe = mean(mspe, na.rm = TRUE), ve_mspe = mean(ve_mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE), 
                wald_pred_cover = mean(wald_pred_cover, na.rm = TRUE), ve_wald_pred_cover = mean(ve_wald_pred_cover, na.rm = TRUE), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                quantile_cred_pred_cover = mean(quantile_cred_pred_cover, na.rm = TRUE), ve_quantile_cred_pred_cover = mean(ve_quantile_cred_pred_cover, na.rm = TRUE),
                quantile_samp_pred_cover = mean(quantile_samp_pred_cover, na.rm = TRUE), ve_quantile_samp_pred_cover = mean(ve_quantile_samp_pred_cover, na.rm = TRUE),
                naive_width = mean(naive_width, na.rm = TRUE), wald_width = mean(wald_width, na.rm = TRUE), ve_wald_width = mean(ve_wald_width, na.rm = TRUE),
                quantile_cred_width = mean(quantile_cred_width, na.rm = TRUE), ve_quantile_cred_width = mean(ve_quantile_cred_width, na.rm = TRUE),
                quantile_samp_width = mean(quantile_samp_width, na.rm = TRUE), ve_quantile_samp_width = mean(ve_quantile_samp_width, na.rm = TRUE)) %>% 
      ungroup()
    
    average_over_n_taxa_var <- left_join(average_over_n_taxa_mu_var, average_over_n_taxa_v_var, by = c("q", "q_obs", "mc_id"))      
    ## (6a) average only over taxa, for e
    e_average_over_taxa <- e_performance_df %>%
      select(q, q_obs, taxon_id, mc_id, mse, sigma_mse, cover, sigma_cover, bias, sigma_bias, sigma) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mse = mean(mse), sigma_mse = mean(sigma_mse), 
                cover = mean(cover), sigma_cover = mean(sigma_cover),
                bias = mean(bias), sigma_bias = mean(sigma_bias), sigma = mean(sigma)) %>% 
      ungroup()
    
    ## (7) average only over n, mc_id
    ## (a) filter by colMeans(mu) > 0
    average_over_n_mu <- mc_averaged_performance %>% 
      select(q, q_obs, subj_id, taxon_id, mse, ve_mse, naive_mse, mspe, ve_mspe, naive_mspe, 
             cover, ve_cover, naive_cover, wald_cover, ve_wald_cover,
             wald_pred_cover, ve_wald_pred_cover, naive_pred_cover,
             quantile_cred_pred_cover, ve_quantile_cred_pred_cover,
             quantile_samp_pred_cover, ve_quantile_samp_pred_cover,
             wald_width, naive_width, ve_wald_width, quantile_cred_width, ve_quantile_cred_width,
             quantile_samp_width, ve_quantile_samp_width,
             bias, ve_bias, naive_bias, sigma) %>%
      group_by(q, q_obs, taxon_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), ve_mse = mean(ve_mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE), 
                cover = mean(cover, na.rm = TRUE), ve_cover = mean(ve_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE), ve_wald_cover = mean(ve_wald_cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE),
                bias = mean(bias, na.rm = TRUE), ve_bias = mean(ve_bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE),
                sigma = mean(sigma)) %>% 
      ungroup()
    average_over_n_v <- mc_averaged_performance %>% 
      select(q, q_obs, subj_id, taxon_id, mse, ve_mse, naive_mse, mspe, ve_mspe, naive_mspe, 
             cover, ve_cover, naive_cover, wald_cover, ve_wald_cover,
             wald_pred_cover, ve_wald_pred_cover, naive_pred_cover,
             quantile_cred_pred_cover, ve_quantile_cred_pred_cover,
             quantile_samp_pred_cover, ve_quantile_samp_pred_cover,
             wald_width, naive_width, ve_wald_width, quantile_cred_width, ve_quantile_cred_width,
             quantile_samp_width, ve_quantile_samp_width,
             bias, ve_bias, naive_bias, sigma) %>%
      group_by(q, q_obs, taxon_id) %>%
      summarize(mspe = mean(mspe, na.rm = TRUE), ve_mspe = mean(ve_mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE), 
                wald_pred_cover = mean(wald_pred_cover, na.rm = TRUE), ve_wald_pred_cover = mean(ve_wald_pred_cover, na.rm = TRUE),
                quantile_cred_pred_cover = mean(quantile_cred_pred_cover, na.rm = TRUE), ve_quantile_cred_pred_cover = mean(ve_quantile_cred_pred_cover, na.rm = TRUE),
                quantile_samp_pred_cover = mean(quantile_samp_pred_cover, na.rm = TRUE), ve_quantile_samp_pred_cover = mean(ve_quantile_samp_pred_cover, na.rm = TRUE),
                naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                naive_width = mean(naive_width, na.rm = TRUE),
                wald_width = mean(wald_width, na.rm = TRUE), ve_wald_width = mean(ve_wald_width, na.rm = TRUE),
                quantile_cred_width = mean(quantile_cred_width, na.rm = TRUE), ve_quantile_cred_width = mean(ve_quantile_cred_width, na.rm = TRUE),
                quantile_samp_width = mean(quantile_samp_width, na.rm = TRUE), ve_quantile_samp_width = mean(ve_quantile_samp_width, na.rm = TRUE)) %>%
      ungroup()
    average_over_n <- left_join(average_over_n_mu, average_over_n_v, by = c("q", "q_obs", "taxon_id"))
    
    
    output_performances[[i]] <- performance_matrix
    output_performances_avg_over_taxa_n[[i]] <- average_over_n_taxa
    output_performances_avg_over_taxa_n_mn[[i]] <- average_over_n_taxa_mn
    output_performances_avg_over_taxa_n_var[[i]] <- average_over_n_taxa_var
    output_performances_avg_over_n[[i]] <- average_over_n
    
    output_performances_e[[i]] <- e_mc_averaged_performance
    output_performances_e_avg_over_taxa[[i]] <- e_average_over_taxa
  }
  saveRDS(output_performances, paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_ab_", args$most_abundant, "_q_obs", args$q_obs, ".rds"))
  saveRDS(output_performances_avg_over_taxa_n, paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_ab_", args$most_abundant, "_q_obs", args$q_obs, ".rds"))
  saveRDS(output_performances_avg_over_taxa_n_mn, paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_mn_ab_", args$most_abundant, "_q_obs", args$q_obs, ".rds"))
  saveRDS(output_performances_avg_over_taxa_n_var, paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_var_ab_", args$most_abundant, "_q_obs", args$q_obs, ".rds"))
  saveRDS(output_performances_avg_over_n, paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_n_ab_", args$most_abundant, "_q_obs", args$q_obs, ".rds"))
  
  saveRDS(output_performances_e, paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_e_n_ab_", args$most_abundant, "_q_obs", args$q_obs, ".rds"))
  saveRDS(output_performances_e_avg_over_taxa, paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_e_avg_over_taxa_ab_", args$most_abundant, "_q_obs", args$q_obs, ".rds"))
  saveRDS(output_log_e_averages, paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_log_e_averages_ab_", args$most_abundant, "_q_obs", args$q_obs, ".rds"))
} else {
  
}
out_perf_lst <- vector("list", length = length(args$q_obs))
out_perf_avg_lst <- out_perf_avg_lst_mn <- out_perf_avg_lst_var <- out_perf_lst
out_perf_e <- out_perf_lst
out_perf_e_avg <- out_perf_lst
out_perf_avg_n_lst <- out_perf_lst
out_log_e_avg <- out_perf_lst
for (j in 1:length(args$q_obs)) {
  out_perf_lst[[j]] <- readRDS(paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
  out_perf_avg_lst[[j]] <- readRDS(paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
  out_perf_avg_lst_mn[[j]] <- readRDS(paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_mn_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
  out_perf_avg_lst_var[[j]] <- readRDS(paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_var_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
  out_perf_avg_n_lst[[j]] <- readRDS(paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_n_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
  
  out_perf_e[[j]] <- readRDS(paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_e_n_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
  out_perf_e_avg[[j]] <- readRDS(paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_e_avg_over_taxa_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
  out_log_e_avg[[j]] <- readRDS(paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_log_e_averages_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
}
output_performances <- do.call(c, out_perf_lst)
output_performances_avg_over_taxa_n <- do.call(c, out_perf_avg_lst)
output_performances_avg_over_taxa_n_mn <- do.call(c, out_perf_avg_lst_mn)
output_performances_avg_over_taxa_n_var <- do.call(c, out_perf_avg_lst_var)
output_performances_e <- do.call(c, out_perf_e)
output_performances_e_avg_over_taxa <- do.call(c, out_perf_e_avg)
output_performances_avg_over_n <- do.call(c, out_perf_avg_n_lst)
output_log_e_averages <- do.call(c, out_log_e_avg)

## transform into matrices
performance_df <- do.call(rbind.data.frame, output_performances)
performance_avg_over_taxa_n <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n)
performance_avg_over_taxa_n_mn <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n_mn)
performance_avg_over_taxa_n_var <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n_var)
performance_e <- do.call(rbind.data.frame, output_performances_e)
performance_e_avg_over_taxa <- do.call(rbind.data.frame, output_performances_e_avg_over_taxa)
performance_df_avg_over_n <- do.call(rbind.data.frame, output_performances_avg_over_n)
log_e_averages <- do.call(rbind.data.frame, output_log_e_averages)

## add new stuff
e_avgs <- lapply(output_performances_e[1:length(output_performances_e)], colMeans)
e_avgs_2 <- lapply(e_avgs, function(x) data.frame(t(x)))
## add on sigma
log_e_averages$sigma <- rep(args$sigma, each = 50)
## make it long
log_es <- log_e_averages %>% 
  gather(key = "taxon", value = "log_e", -mc_id, -sigma, -estimator) %>% 
  mutate(taxon_id = as.numeric(lapply(strsplit(taxon, "_", fixed = TRUE), function(x) x[2]))) %>% 
  select(mc_id, sigma, taxon_id, estimator, log_e)

## create new tibble with rmspe, cover, pred_cover and estimator type
performance_by_type <- performance_avg_over_taxa_n %>% 
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_wald_pred_cover = wald_pred_cover,
         no_ve_quantile_cred_pred_cover = quantile_cred_pred_cover, no_ve_quantile_samp_pred_cover = quantile_samp_pred_cover,
         naive_wald_width = naive_width, no_ve_wald_width = wald_width, no_ve_quantile_cred_width = quantile_cred_width, no_ve_quantile_samp_width = quantile_samp_width,
         no_ve_bias = bias,
         naive_wald_pred_cover = naive_pred_cover,
         no_ve_wald_cover = wald_cover) %>%
  select(-mse, -mspe, -cover, -naive_width, -wald_pred_cover, -quantile_cred_pred_cover, -quantile_samp_pred_cover,
         -wald_width, -quantile_cred_width, -quantile_samp_width, -bias, -naive_pred_cover, -wald_cover) %>% 
  gather(key, value, -q, -q_obs, -mc_id, -sigma) %>% 
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>% 
  spread(measure, value) %>% 
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe), sq_bias = bias^2)

performance_by_type_mn <- performance_avg_over_taxa_n_mn %>% 
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_wald_pred_cover = wald_pred_cover,
         no_ve_quantile_cred_pred_cover = quantile_cred_pred_cover, no_ve_quantile_samp_pred_cover = quantile_samp_pred_cover,
         naive_wald_width = naive_width, no_ve_wald_width = wald_width, no_ve_quantile_cred_width = quantile_cred_width, no_ve_quantile_samp_width = quantile_samp_width,
         no_ve_bias = bias,
         naive_wald_pred_cover = naive_pred_cover,
         no_ve_wald_cover = wald_cover) %>%
  select(-mse, -mspe, -cover, -wald_pred_cover, -quantile_cred_pred_cover, -quantile_samp_pred_cover,
         -naive_width, -wald_width, -quantile_cred_width, -quantile_samp_width, -bias, -naive_pred_cover, -wald_cover) %>% 
  gather(key, value, -q, -q_obs, -mc_id, -sigma) %>% 
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>% 
  spread(measure, value) %>% 
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe), sq_bias = bias^2)

performance_by_type_var <- performance_avg_over_taxa_n_var %>% 
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_wald_pred_cover = wald_pred_cover,
         no_ve_quantile_cred_pred_cover = quantile_cred_pred_cover, no_ve_quantile_samp_pred_cover = quantile_samp_pred_cover,
         naive_wald_width = naive_width, no_ve_wald_width = wald_width, no_ve_quantile_cred_width = quantile_cred_width, no_ve_quantile_samp_width = quantile_samp_width,
         no_ve_bias = bias,
         naive_wald_pred_cover = naive_pred_cover,
         no_ve_wald_cover = wald_cover) %>%
  select(-mse, -mspe, -cover, -wald_pred_cover, -quantile_cred_pred_cover, -quantile_samp_pred_cover,
         -naive_width, -wald_width, -quantile_cred_width, -quantile_samp_width, -bias, -naive_pred_cover, -wald_cover) %>% 
  gather(key, value, -q, -q_obs, -mc_id, -sigma) %>% 
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>% 
  spread(measure, value) %>% 
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe), sq_bias = bias^2)

performance_by_type_over_n <- performance_df_avg_over_n %>% 
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_wald_pred_cover = wald_pred_cover,
         no_ve_quantile_cred_pred_cover = quantile_cred_pred_cover, no_ve_quantile_samp_pred_cover = quantile_samp_pred_cover,
         naive_wald_width = naive_width, no_ve_wald_width = wald_width, no_ve_quantile_cred_width = quantile_cred_width, no_ve_quantile_samp_width = quantile_samp_width,
         no_ve_bias = bias,
         naive_wald_pred_cover = naive_pred_cover,
         no_ve_wald_cover = wald_cover) %>%
  select(-mse, -mspe, -cover, -wald_pred_cover, -quantile_cred_pred_cover, -quantile_samp_pred_cover,
         -naive_width, -wald_width, -quantile_cred_width, -quantile_samp_width, -bias, -naive_pred_cover, -wald_cover) %>% 
  gather(key, value, -q, -q_obs, -taxon_id, -sigma) %>% 
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>% 
  spread(measure, value) %>% 
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe), sq_bias = bias^2)

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
point_dodge <- 0.05
axis_cex <- 1.75

## --------------------------------------------------------------------------------------------------
## MAIN PLOT: prediction interval coverage, credible interval coverage, rms(p)e
## for q_obs == 7 only
## --------------------------------------------------------------------------------------------------
pred_cover_plot <- performance_by_type %>%
  filter(q_obs == 7) %>% 
  group_by(q, q_obs, sigma, estimator) %>% 
  summarize(mean_cover = mean(wald_pred_cover)) %>% 
  ggplot(aes(x = sigma, y = mean_cover, group = paste(sigma, estimator, sep = "_"),
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), ordered = FALSE, labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  xlab(expression(sigma[e])) +
  ylab("Coverage") + 
  ggtitle("Prediction interval coverage for V") + 
  labs(shape = "Estimator type") +
  # geom_boxplot(width = 0.08, position = position_dodge(), outlier.shape = NA) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.2, 1)) +
  theme(legend.position = c(0.185, 0.2)) +
  guides(shape = FALSE)
cover_plot <- performance_by_type %>%
  filter(q_obs == 7) %>% 
  group_by(q, q_obs, sigma, estimator) %>% 
  summarize(mean_cover = mean(cover)) %>% 
  ggplot(aes(x = sigma, y = mean_cover, group = paste(sigma, estimator, sep = "_"), 
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), ordered = FALSE, labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  xlab(expression(sigma[e])) +
  ylab("Coverage") + 
  ggtitle(expression(bold(paste("Interval coverage for ", mu, sep = "")))) + 
  labs(shape = "Estimator type") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) + 
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.2, 1)) +
  theme(legend.position = c(0.05, 0.2))

rmspe_plot <- performance_by_type %>%
  filter(q_obs == 7) %>%
  group_by(q, q_obs, sigma, estimator) %>% 
  summarize(mean_rmspe = mean(rmspe[!is.infinite(rmspe)])) %>% 
  ggplot(aes(x = sigma, y = mean_rmspe, group = paste(sigma, estimator, sep = "_"), 
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), ordered = FALSE, labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  xlab(expression(sigma[e])) +
  ylab("Root mean squared prediction error") + 
  ggtitle("RMSPE") + 
  labs(shape = "Estimator type") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) + 
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 3e1)) +
  guides(shape = FALSE)

rmse_plot <- performance_by_type %>%
  filter(q_obs == 7) %>%
  group_by(q, q_obs, sigma, estimator) %>% 
  summarize(mean_rmse = mean(rmse)) %>% 
  ggplot(aes(x = sigma, y = mean_rmse, group = paste(sigma, estimator, sep = "_"), 
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), ordered = FALSE, labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  xlab(expression(sigma[e])) +
  ylab("Root mean squared error") + 
  ggtitle("RMSE") + 
  labs(shape = "Estimator type") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) + 
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 1e3)) +
  guides(shape = FALSE) 

ggsave(paste0(plots_dir, "e_constraint.png"), plot = plot_grid(cover_plot, pred_cover_plot, 
                                                               rmse_plot, rmspe_plot),
       device = "png", width = 30, height = 25, units = "cm", dpi = 300)

## debug: coverage with filtering by mean, variance
cover_plot_mn <- performance_by_type_mn %>%
  filter(q_obs == 7) %>% 
  group_by(q, q_obs, sigma, estimator) %>% 
  summarize(mean_cover = mean(cover)) %>% 
  ggplot(aes(x = sigma, y = mean_cover, group = paste(sigma, estimator, sep = "_"), 
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), ordered = FALSE, labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  xlab(expression(sigma[e])) +
  ylab("Coverage") + 
  ggtitle(expression(paste("Interval coverage for ", mu, "; filter by mean(W) > 0.5", sep = ""))) + 
  labs(shape = "Estimator type") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) + 
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.2, 1)) +
  theme(legend.position = c(0.185, 0.2)) 
cover_plot_var <- performance_by_type_mn %>%
  filter(q_obs == 7) %>% 
  group_by(q, q_obs, sigma, estimator) %>% 
  summarize(mean_cover = mean(cover)) %>% 
  ggplot(aes(x = sigma, y = mean_cover, group = paste(sigma, estimator, sep = "_"), 
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), ordered = FALSE, labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  xlab(expression(sigma[e])) +
  ylab("Coverage") + 
  ggtitle(expression(paste("Interval coverage for ", mu, "; filter by var(W) > 1", sep = ""))) + 
  labs(shape = "Estimator type") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) + 
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.2, 1)) +
  theme(legend.position = c(0.185, 0.2)) +
  guides(shape = FALSE, alpha = FALSE)
ggsave(paste0(plots_dir, "e_constraint_mn_var_cover.png"), plot = plot_grid(cover_plot_mn, cover_plot_var),
       device = "png", width = 30, height = 15, units = "cm", dpi = 300)
## --------------------------------------------------------------------------------------------------
## SUPPLEMENTARY PLOTS:
## (a) same as main plot, but for q_obs == 3
## (b) widths for both q_obs
## (c) squared bias
## (d) efficiencies?
## --------------------------------------------------------------------------------------------------
pred_cover_plot_3 <- performance_by_type %>%
  filter(q_obs == 3) %>% 
  group_by(q, q_obs, sigma, estimator) %>% 
  summarize(mean_cover = mean(wald_pred_cover)) %>% 
  ggplot(aes(x = sigma, y = mean_cover, group = paste(sigma, estimator, sep = "_"),
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), ordered = FALSE, labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  xlab(expression(sigma[e])) +
  ylab("Coverage") + 
  ggtitle("Prediction interval coverage") + 
  labs(shape = "Estimator type") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) + 
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.4, 1)) +
  theme(legend.position = c(0.185, 0.2)) +
  guides(shape = FALSE, alpha = FALSE)
cover_plot_3 <- performance_by_type %>%
  filter(q_obs == 3) %>% 
  group_by(q, q_obs, sigma, estimator) %>% 
  summarize(mean_cover = mean(cover)) %>% 
  ggplot(aes(x = sigma, y = mean_cover, group = paste(sigma, estimator, sep = "_"), 
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), ordered = FALSE, labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  xlab(expression(sigma[e])) +
  ylab("Coverage") + 
  ggtitle("Credible interval coverage") + 
  labs(shape = "Estimator type") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) + 
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.4, 1)) +
  theme(legend.position = c(0.185, 0.2)) +
  guides(shape = FALSE, alpha = FALSE)

rmspe_plot_3 <- performance_by_type %>%
  filter(q_obs == 3) %>%
  group_by(q, q_obs, sigma, estimator) %>% 
  summarize(mean_rmspe = mean(rmspe)) %>% 
  ggplot(aes(x = sigma, y = mean_rmspe, group = paste(sigma, estimator, sep = "_"), 
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), ordered = FALSE, labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  xlab(expression(sigma[e])) +
  ylab("Root mean squared prediction error") + 
  ggtitle("RMSPE") + 
  labs(shape = "Estimator type") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) + 
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 200)) +
  theme(legend.position = c(0.8, 0.3))

rmse_plot_3 <- performance_by_type %>%
  filter(q_obs == 3) %>%
  group_by(q, q_obs, sigma, estimator) %>% 
  summarize(mean_rmse = mean(rmse)) %>% 
  ggplot(aes(x = sigma, y = mean_rmse, group = paste(sigma, estimator, sep = "_"), 
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), ordered = FALSE, labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  xlab(expression(sigma[e])) +
  ylab("Root mean squared error") + 
  ggtitle("RMSE") + 
  labs(shape = "Estimator type") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) + 
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 2500)) +
  guides(shape = FALSE)

png(paste0(plots_dir, "e_constraint_q_obs_3.png"), width = fig_width, height = fig_height, units = "px", res = 300)
gridExtra::grid.arrange(pred_cover_plot_3, cover_plot_3, 
                        rmspe_plot_3, rmse_plot_3)
dev.off()

## width
width_plot <- performance_by_type %>% 
  filter(q_obs == 7) %>% 
  group_by(sigma, estimator) %>% 
  summarize(mean_width = mean(wald_width[!is.infinite(wald_width)])) %>% 
  ggplot(aes(x = sigma, y = mean_width, group = paste(sigma, estimator, sep = "_"),
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), ordered = FALSE, labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  xlab(expression(sigma[e])) +
  ylab("Prediction interval width") + 
  ggtitle("Prediction interval width") + 
  labs(shape = "Estimator type") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 4e2)) +
  guides(alpha = guide_legend(override.aes = list(shape = "blue"))) 
ggsave(paste0(plots_dir, "e_constraint_width.png"), plot = width_plot,
       width = 15, height = 15, units = "cm", dpi = 300)

## mse as bias^2 + variance
bias_plot <- performance_by_type %>%
  filter(q_obs == 7) %>%
  mutate(sq_bias = bias^2) %>% 
  group_by(sigma, estimator) %>% 
  summarize(mn_bias = mean(sq_bias)) %>% 
  ggplot(aes(x = sigma, y = mn_bias, group = paste(sigma, estimator, sep = "_"), 
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), ordered = FALSE, labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  xlab(expression(sigma[e])) +
  ylab("Squared bias") + 
  ggtitle("Squared bias") + 
  labs(shape = "Estimator type") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) + 
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 2000))
variance_plot <- performance_by_type %>%
  filter(q_obs == 7) %>% 
  mutate(variance = mse - sq_bias) %>% 
  group_by(sigma, estimator) %>% 
  summarize(mn_var = mean(variance)) %>% 
  ggplot(aes(x = sigma, y = mn_var, group = paste(sigma, estimator, sep = "_"), 
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), ordered = FALSE, labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  xlab(expression(sigma[e])) +
  ylab("Estimated variance") + 
  ggtitle("Estimated variance") + 
  labs(shape = "Estimator type") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) + 
  scale_y_log10() +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") 

ggsave(paste0(plots_dir, "e_constraint_bias.png"), plot = plot_grid(bias_plot, variance_plot),
       width = 30, height = 15, units = "cm", dpi = 300)

## bias vs rank order in W
bias_vs_rank_order_0 <- performance_by_type_over_n %>% 
  filter(q_obs == 7, sigma == 0) %>% 
  group_by(taxon_id, estimator) %>% 
  summarize(mn_bias = mean(bias)) %>% 
  ggplot(aes(x = taxon_id, y = mn_bias, group = taxon_id, color = factor(taxon_id), 
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(-30, 30)) +
  labs(color = "Taxon") +
  labs(shape = "Estimator") +
  labs(y = "Bias") +
  labs(x = "Rank order of W") +
  ggtitle(expression(bold(paste("Bias vs. rank order, ", sigma[e], " = 0"))))

bias_vs_rank_order_0.5 <- performance_by_type_over_n %>% 
  filter(q_obs == 7, sigma == 0.5) %>% 
  group_by(taxon_id, estimator) %>% 
  summarize(mn_bias = mean(bias)) %>% 
  ggplot(aes(x = taxon_id, y = mn_bias, group = taxon_id, color = factor(taxon_id), 
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(-30, 30)) +
  labs(color = "Taxon") +
  labs(shape = "Estimator") +
  labs(y = "Bias") +
  labs(x = "Rank order of W") +
  ggtitle(expression(bold(paste("Bias vs. rank order, ", sigma[e], " = 0.5"))))

bias_vs_rank_order_1 <- performance_by_type_over_n %>% 
  filter(q_obs == 7, sigma == 1) %>% 
  group_by(taxon_id, estimator) %>% 
  summarize(mn_bias = mean(bias)) %>% 
  ggplot(aes(x = taxon_id, y = mn_bias, group = taxon_id, color = factor(taxon_id), 
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), labels = c("Naive", "Proposed, no renorm", "Proposed, renorm")))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(-30, 30)) +
  labs(color = "Taxon") +
  labs(shape = "Estimator") +
  labs(y = "Bias") +
  labs(x = "Rank order of W") +
  ggtitle(expression(bold(paste("Bias vs. rank order, ", sigma[e], " = 1"))))

ggsave(paste0(plots_dir, "e_constraint_bias_vs_rank_order.png"), 
       plot = plot_grid(bias_vs_rank_order_0, bias_vs_rank_order_0.5,
                        bias_vs_rank_order_1),
       width = 40, height = 40, units = "cm", dpi = 300)
## ---------------------------------------------------------------------
## distribution of estimated efficiencies when efficiency truly is zero!
## ---------------------------------------------------------------------
results_dir <- paste0("results/", args$sim_name, "/ve_", args$sigma[1], "/cov_", args$corr, "/n_", args$N, "/q_", args$q, "/q_obs_", args$q_obs[1], "/")

dir_mat <- expand.grid(job = 1:args$num_jobs, dir = results_dir, stan_model = args$stan_model)
mod_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_mod_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
data_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_data_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))

mod_lst <- lapply(mod_nms_lst, read_func)
data_lst <- lapply(data_nms_lst, read_func)

## pick out the ones where I've estimated varying efficiency; get MC reps for all taxa
get_e_rows <- function(x) {
  if (any(is.na(x))) {
    return(NA)
  } else {
    return(grepl("e", rownames(x)) & !grepl("beta", rownames(x)))
  }
}
get_sigma_row <- function(x) {
  if (any(is.na(x))) {
    return(NA)
  } else {
    return(grepl("sigma", rownames(x)))
  }
}
e_lst <- lapply(mod_lst[(args$num_jobs + 1):(2*args$num_jobs)], function(x) tryCatch(rbind(x[get_e_rows(x), ], x[get_sigma_row(x), ]), error = function(e) NA))

## pick out the means for both e and sigma
e_df <- do.call(rbind.data.frame, lapply(e_lst, function(x) tryCatch(x[, "mean"], error = function(e) NA)))
colnames(e_df) <- c(paste0("taxon_", 1:args$q), "sigma")
e_df$mc_id <- 1:args$num_jobs

e_df %>% 
  ggplot(aes(y = taxon_1))  +
  geom_histogram() 

png(paste0(plots_dir, "e_constraint_sigma.png"), width = fig_width, height = fig_height, units = "px", res = 300)
e_df %>% 
  ggplot(aes(y = sigma)) +
  geom_boxplot() +
  ggtitle(expression(paste("Estimated ", sigma[e], " when true ", sigma[e], " = 0"))) +
  labs(y = expression(hat(sigma)[e])) 
dev.off()

## ---------------------------------------------------------------------------------
## mean of the log-efficiencies! investigate soft thresholding
## ---------------------------------------------------------------------------------
mean_log_e_plot <- log_es %>% 
  group_by(sigma, mc_id, estimator) %>% 
  summarize(mn_log_e = mean(log_e)) %>% 
  ggplot(aes(x = sigma, y = mn_log_e, group = factor(paste0(sigma, estimator)),
         fill = factor(estimator, levels = c("soft_constraint", "e_normalized"), labels = c("Soft constraint", "e normalized"), ordered = FALSE))) +
  geom_boxplot() +
  scale_fill_manual(values = cols) +
  ggtitle("Log efficiency averaged over q = 20 taxa") +
  labs(y = "Mean log efficiency") +
  labs(x = expression(sigma[e])) +
  labs(fill = "Estimator")
mean_log_e_qpcr_observed_plot <- log_es %>% 
  filter(taxon_id <= 7) %>% 
  group_by(sigma, mc_id, estimator) %>% 
  summarize(mn_log_e = mean(log_e)) %>% 
  ggplot(aes(x = sigma, y = mn_log_e, group = factor(paste0(sigma, estimator)),
             fill = factor(estimator, levels = c("soft_constraint", "e_normalized"), labels = c("Soft constraint", "e normalized"), ordered = FALSE))) +
  geom_boxplot() +
  scale_fill_manual(values = cols) +
  ggtitle("Log efficiency averaged over qPCR-observed taxa") +
  labs(y = "Mean log efficiency") +
  labs(x = expression(sigma[e])) +
  labs(fill = "Estimator")

ggsave(paste0(plots_dir, "mean_log_e.png"),
       plot = plot_grid(mean_log_e_plot, mean_log_e_qpcr_observed_plot),
       width = 18, height = 8)

## ---------------------------------------------------------------------------------
## convergence diagnostics
## ---------------------------------------------------------------------------------
for (j in 1:length(args$q_obs)) {
  rhats <- list(length(args$sigma))
  for (i in 1:length(args$sigma)) {
    ## load in results and data
    results_dir <- paste0("results/", args$sim_name, "/ve_", args$sigma[i], "/cov_", args$corr, "/n_", args$N, "/q_", args$q, "/q_obs_", args$q_obs[j], "/")
    
    dir_mat <- expand.grid(job = 1:args$num_jobs, dir = results_dir, stan_model = args$stan_model)
    mod_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_mod_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
    data_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_data_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
    
    mod_lst <- lapply(mod_nms_lst, read_func)
    data_lst <- lapply(data_nms_lst, read_func)
    print(paste0("Showing sigma = ", args$sigma[i]))
    beta_summ_lst <- lapply(mod_lst, function(x) summary(x[grepl("beta", rownames(x)), 7]))
    beta_mat <- do.call(rbind.data.frame, beta_summ_lst)
    beta_mat$est <- c(rep("no renorm", 50), rep("renorm", 50))
    colnames(beta_mat) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "est")
    Sigma_summ_lst <- lapply(mod_lst, function(x) summary(x[grepl("Sigma", rownames(x)), 7]))
    Sigma_mat <- do.call(rbind.data.frame, Sigma_summ_lst)
    Sigma_mat$est <- c(rep("no renorm", 50), rep("renorm", 50))
    colnames(Sigma_mat) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "est")
    sigma_summ_lst <- lapply(mod_lst, function(x) summary(x[grepl("sigma", rownames(x)), 7]))
    sigma_mat <- do.call(rbind.data.frame, sigma_summ_lst)
    sigma_mat$est <- c(rep("no renorm", 50), rep("renorm", 50))
    colnames(sigma_mat) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "est")
    mu_summ_lst <- lapply(mod_lst, function(x) summary(x[grepl("mu", rownames(x)), 7])[1:6])
    mu_mat <- do.call(rbind.data.frame, mu_summ_lst)
    mu_mat$est <- c(rep("no renorm", 50), rep("renorm", 50))
    colnames(mu_mat) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "est")
    
    beta_means <- rbind(colMeans(subset(beta_mat, est == "no renorm")[, 1:6]), colMeans(subset(beta_mat, est == "renorm")[, 1:6]))
    Sigma_means <- rbind(colMeans(subset(Sigma_mat, est == "no renorm")[, 1:6]), colMeans(subset(Sigma_mat, est == "renorm")[, 1:6]))
    mu_means <- rbind(colMeans(subset(mu_mat, est == "no renorm")[, 1:6]), colMeans(subset(mu_mat, est == "renorm")[, 1:6]))
    sigma_means <- rbind(colMeans(subset(sigma_mat, est == "no renorm")[, 1:6]), colMeans(subset(sigma_mat, est == "renorm")[, 1:6]))
    
    print("beta, no renorm:")
    print(colMeans(subset(beta_mat, est == "no renorm")[, 1:6]))
    print("beta, renorm:")
    print(colMeans(subset(beta_mat, est == "renorm")[, 1:6]))
    
    print("Sigma, no renorm:")
    print(colMeans(subset(Sigma_mat, est == "no renorm")[, 1:6]))
    print("Sigma, renorm:")
    print(colMeans(subset(Sigma_mat, est == "renorm")[, 1:6]))
    
    print("mu, no renorm:")
    print(colMeans(subset(mu_mat, est == "no renorm")[, 1:6]))
    print("mu, renorm:")
    print(colMeans(subset(mu_mat, est == "renorm")[, 1:6]))
    
    print("sigma, no renorm:")
    print(colMeans(subset(sigma_mat, est == "no renorm")[, 1:6]))
    print("sigma, renorm:")
    print(colMeans(subset(sigma_mat, est == "renorm")[, 1:6]))
    rhats[[i]] <- tibble(est = rep(c("no renorm", "renorm"), 4), sigma_e = args$sigma[i], parameter = c(rep("beta", 2), rep("Sigma", 2), rep("mu", 2), rep("sigma", 2)))
    rhats[[i]] <- tibble::add_column(rhats[[i]], median = c(beta_means[, 3], Sigma_means[, 3], mu_means[, 3], sigma_means[, 3]))
    rhats[[i]] <- tibble::add_column(rhats[[i]], iqr_low = c(beta_means[, 2], Sigma_means[, 2], mu_means[, 2], sigma_means[, 2]))
    rhats[[i]] <- tibble::add_column(rhats[[i]], iqr_high = c(beta_means[, 5], Sigma_means[, 5], mu_means[, 5], sigma_means[, 5]))
  }
}
rhats_all <- do.call(rbind.data.frame, rhats)
write.csv(rhats_all, "plots/e_constraint/rhats.csv")
