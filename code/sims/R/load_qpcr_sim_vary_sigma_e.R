#!/usr/local/bin/Rscript
##################################################################################
## FILE: load_vary_sigma_e.R
##
## CREATED: 28 December 2018 by Brian Williamson
##
## PURPOSE: load in results from a qPCR simulation varying sigma_e, create summaries
##
## INPUTS: listed below
##         
## OUTPUTS: plots, tables, etc.
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
library("paramedic")
source("code/R/load_qpcr_sim_helpers.R") # provides function get_summaries


## grab command-line arguments
parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "vary_sigma_e", help = "name of the simulation")
parser$add_argument("--stan-model", default = "predict_qpcr",
                    help = "Which Stan model file to use.")
parser$add_argument("--N", type = "double", default = 50, help = "sample size")
parser$add_argument("--q", nargs = "+", type = "double", default = 40, help = "how many taxa do we have?")
parser$add_argument("--corr", type = "double", default = 0, help = "hyperparameter controlling the off-diagonal elements of Sigma.")
parser$add_argument("--sigma", type = "double", default = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), help = "SD of efficiencies (0 = no varying efficiency)")
parser$add_argument("--ad", type = "double", default = 0.850000, help = "adapt delta")
parser$add_argument("--mt", type = "double", default = 15, help = "max treedepth")
parser$add_argument("--num-jobs", type = "double", default = 50, help = "number of jobs run")
parser$add_argument("--taxa", type = "double", default = 10, help = "taxon of interest")
parser$add_argument("--q-obs", nargs = "+", type = "double", default = c(3,7), help = "number of taxa with observed qPCR")
parser$add_argument("--most-abundant", type = "double", default = 1, help = "whether or not to use most abundant taxa")
parser$add_argument("--read-data", type = "double", default = 0, help = "whether or not to read in raw output, or instead do summaries")
args <- parser$parse_args()

## add in with varying efficiency
args$stan_model <- c("predict_qpcr", "predict_qpcr_with_varying_efficiency")

read_func <- function(x) tryCatch(readRDS(x), error = function(e) NA)

plots_dir <- paste0("plots/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
save_res_dir <- paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/")
if (!dir.exists(save_res_dir)) dir.create(save_res_dir, recursive = TRUE)
## --------------------------------------------------------------------------------------------------
## CREATE LIST OF OUTPUT FOR EASY SUMMARIES; GET THOSE SUMMARIES
## --------------------------------------------------------------------------------------------------
if (args$read_data) {
  for (j in 1:length(args$q_obs)) {
    output_performances <- vector(mode = "list", length = length(args$sigma))
    output_performances_single_taxon <- vector(mode = "list", length = length(args$sigma))
    
    output_performances_avg_over_taxa_n <- vector(mode = "list", length = length(args$sigma))
    output_performances_e <- vector(mode = "list", length = length(args$sigma))
    output_performances_e_avg_over_taxa <- vector(mode = "list", length = length(args$sigma))
    
    for (i in 1:length(args$sigma)) {
      ## load in results and data
      results_dir <- paste0("results/", args$sim_name, "/ve_", args$sigma[i], "/cov_", args$corr, "/n_", args$N, "/q_", args$q, "/q_obs_", args$q_obs[j], "/")
      
      dir_mat <- expand.grid(job = 1:args$num_jobs, dir = results_dir, stan_model = args$stan_model)
      mod_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_mod_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
      data_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_data_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
      samp_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_samps_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
      
      mod_lst <- lapply(mod_nms_lst, read_func)
      data_lst <- lapply(data_nms_lst, read_func)
      samp_lst <- lapply(samp_nms_lst, read_func)
      
      ## (1) pair model summaries with relevant data, for all taxa, all q_obs
      summary_df_no_ve <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type) get_summaries(w, x, y, z, type), mod_lst[1:args$num_jobs], data_lst[1:args$num_jobs], samp_lst[1:args$num_jobs], MoreArgs = list(z = 1:args$q, type = "no_ve"), SIMPLIFY = FALSE))
      summary_df_ve <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type) get_summaries(w, x, y, z, type), mod_lst[(args$num_jobs + 1):(2*args$num_jobs)], data_lst[(args$num_jobs + 1):(2*args$num_jobs)], samp_lst[(args$num_jobs + 1):(2*args$num_jobs)], MoreArgs = list(z = 1:args$q, type = "ve"), SIMPLIFY = FALSE))
      ## pair with the monte-carlo id, type of estimator
      summary_df_no_ve$mc_id <- rep(NA, dim(summary_df_no_ve)[1])
      summary_df_ve$mc_id <- rep(NA, dim(summary_df_ve)[1])
      no_ve_mc_lst <- 1:50
      ve_mc_lst <- 1:50
      if (i == 6 & j == 2) no_ve_mc_lst <- no_ve_mc_lst[-49]
      if (i == 4 & j == 1) ve_mc_lst <- ve_mc_lst[-27]
      summary_df_no_ve$mc_id[!is.na(summary_df_no_ve$q)] <- rep(rep(no_ve_mc_lst, each = args$N*args$q), 1)
      summary_df_ve$mc_id[!is.na(summary_df_ve$q)] <- rep(rep(ve_mc_lst, each = args$N*args$q), 1)
      
      ## remove naive est from ve summaries, change names, and merge with no_ve summaries
      summary_df_ve_2 <- summary_df_ve %>% 
        select(-naive_est) %>% 
        mutate(ve_est = est, ve_sd = sd, ve_cil = cil, ve_ciu = ciu, 
               ve_wald_pred_cil = wald_pred_cil, ve_wald_pred_ciu = wald_pred_ciu,
               ve_quantile_cred_pred_cil = quantile_cred_pred_cil, ve_quantile_cred_pred_ciu = quantile_cred_pred_ciu,
               ve_quantile_samp_pred_cil = quantile_samp_pred_cil, ve_quantile_samp_pred_ciu = quantile_samp_pred_ciu) %>% 
        select(q, q_obs, subj_id, taxon_id, mc_id, mu, qpcr, ve_est, ve_sd, ve_cil, ve_ciu, 
               ve_wald_pred_cil, ve_wald_pred_ciu, ve_quantile_cred_pred_cil, ve_quantile_cred_pred_ciu,
               ve_quantile_samp_pred_cil, ve_quantile_samp_pred_ciu)
      summary_df_no_ve_2 <- summary_df_no_ve %>% 
        select(-type)
      summary_df <- left_join(summary_df_no_ve_2, summary_df_ve_2, by = c("q", "q_obs", "subj_id", "taxon_id", "mc_id", "mu", "qpcr"))
      
      ## (1a) pair model summaries with data for efficiency
      e_summary_df <- do.call(rbind.data.frame, mapply(function(x, y, z) get_summaries_with_ve(x, y, z), mod_lst[grepl("with_varying_efficiency", mod_nms_lst)], data_lst[grepl("with_varying_efficiency", mod_nms_lst)], MoreArgs = list(z = 1:args$q), SIMPLIFY = FALSE))
      e_summary_df$mc_id <- rep(NA, dim(e_summary_df)[1])
      e_summary_df$mc_id[!is.na(e_summary_df$q)] <- rep(ve_mc_lst, each = args$q)
      
      ## (2) compute performance for each row
      performance_df <- summary_df %>%
        mutate(mse = (mu - est)^2, ve_mse = (mu - ve_est)^2, naive_mse = (mu - naive_est)^2,
               mspe = (qpcr - est)^2, ve_mspe = (qpcr - ve_est)^2, naive_mspe = (qpcr - naive_est)^2,
               cover = cil <= mu & ciu >= mu,
               ve_cover = ve_cil <= mu & ve_ciu >= mu,
               wald_pred_cover = wald_pred_cil <= qpcr & wald_pred_ciu >= qpcr,
               ve_wald_pred_cover = ve_wald_pred_cil <= qpcr & ve_wald_pred_ciu >= qpcr,
               quantile_cred_pred_cover = quantile_cred_pred_cil <= qpcr & quantile_cred_pred_ciu >= qpcr,
               ve_quantile_cred_pred_cover = ve_quantile_cred_pred_cil <= qpcr & ve_quantile_cred_pred_ciu >= qpcr,
               quantile_samp_pred_cover = quantile_samp_pred_cil <= qpcr & quantile_samp_pred_ciu >= qpcr,
               ve_quantile_samp_pred_cover = ve_quantile_samp_pred_cil <= qpcr & ve_quantile_samp_pred_ciu >= qpcr,
               wald_width = abs(wald_pred_ciu - wald_pred_cil),
               quantile_cred_width = abs(quantile_cred_pred_ciu - quantile_cred_pred_cil),
               quantile_samp_width = abs(quantile_samp_pred_ciu - quantile_samp_pred_cil),
               ve_wald_width = abs(ve_wald_pred_ciu - ve_wald_pred_cil),
               ve_quantile_cred_width = abs(ve_quantile_cred_pred_ciu - ve_quantile_cred_pred_cil),
               ve_quantile_samp_width = abs(ve_quantile_samp_pred_ciu - ve_quantile_samp_pred_cil),
               bias = mu - est, ve_bias = mu - ve_est, naive_bias = mu - naive_est,
               sigma = args$sigma[i])
      
      performance_df$include_taxa_in_avg <- performance_df$taxon_id >= (args$q_obs[j] + 1)
      
      ## (2a) compute mses, coverages for each row
      e_performance_df <- e_summary_df %>% 
        mutate(mse = (e - est_e)^2, sigma_mse = (sigma - est_sigma)^2,
               cover = e_cil <= e & e_ciu >= e,
               sigma_cover = sigma_cil <= sigma & sigma_ciu >= sigma,
               bias = e - est_e, sigma_bias = sigma - est_sigma,
               sigma = args$sigma[i])
      
      ## (3) average over MC reps for each taxon
      mc_averaged_performance <- performance_df %>%
        select(q, q_obs, subj_id, taxon_id, mse, ve_mse, naive_mse, mspe, ve_mspe, naive_mspe, 
               cover, ve_cover, wald_pred_cover, ve_wald_pred_cover,
               quantile_cred_pred_cover, ve_quantile_cred_pred_cover,
               quantile_samp_pred_cover, ve_quantile_samp_pred_cover,
               wald_width, ve_wald_width, quantile_cred_width, ve_quantile_cred_width,
               quantile_samp_width, ve_quantile_samp_width,
               bias, ve_bias, naive_bias, sigma) %>%
        group_by(q, q_obs, subj_id, taxon_id) %>%
        summarize(mse = mean(mse, na.rm = TRUE), ve_mse = mean(ve_mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE), 
                  mspe = mean(mspe, na.rm = TRUE), ve_mspe = mean(ve_mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE), 
                  cover = mean(cover, na.rm = TRUE), ve_cover = mean(ve_cover, na.rm = TRUE),
                  wald_pred_cover = mean(wald_pred_cover, na.rm = TRUE), ve_wald_pred_cover = mean(ve_wald_pred_cover, na.rm = TRUE),
                  quantile_cred_pred_cover = mean(quantile_cred_pred_cover, na.rm = TRUE), ve_quantile_cred_pred_cover = mean(ve_quantile_cred_pred_cover, na.rm = TRUE),
                  quantile_samp_pred_cover = mean(quantile_samp_pred_cover, na.rm = TRUE), ve_quantile_samp_pred_cover = mean(ve_quantile_samp_pred_cover, na.rm = TRUE),
                  wald_width = mean(wald_width, na.rm = TRUE), ve_wald_width = mean(ve_wald_width, na.rm = TRUE),
                  quantile_cred_width = mean(quantile_cred_width, na.rm = TRUE), ve_quantile_cred_width = mean(ve_quantile_cred_width, na.rm = TRUE),
                  quantile_samp_width = mean(quantile_samp_width, na.rm = TRUE), ve_quantile_samp_width = mean(ve_quantile_samp_width, na.rm = TRUE),
                  bias = mean(bias, na.rm = TRUE), ve_bias = mean(ve_bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE),
                  sigma = mean(sigma)) %>%
        ungroup()
      
      ## (3a) average over MC reps for each taxon, for e
      e_mc_averaged_performance <- e_performance_df %>% 
        select(q, q_obs, taxon_id, mse, sigma_mse, cover, sigma_cover,
               bias, sigma_bias, sigma) %>%
        group_by(q, q_obs, taxon_id) %>%
        summarize(mse = mean(mse, na.rm = TRUE), sigma_mse = mean(sigma_mse, na.rm = TRUE),
                  cover = mean(cover, na.rm = TRUE), sigma_cover = mean(sigma_cover, na.rm = TRUE),
                  bias = mean(bias, na.rm = TRUE), sigma_bias = mean(sigma_bias, na.rm = TRUE),
                  sigma = mean(sigma)) %>%
        ungroup()
      
      ## (4) set flag for taxa of interest
      mc_averaged_performance$include_taxa_in_avg <- apply(mc_averaged_performance, 1, 
                                                           function(x) ifelse(!is.na(x[1]), x[4] %in% tail(1:x[1], x[1] - x[2]), NA))
      ## (5) average over taxa of interest for each q_obs
      performance_across_taxa <- mc_averaged_performance %>%
        group_by(q, q_obs, subj_id) %>%
        filter(include_taxa_in_avg) %>%
        summarize(mse = mean(mse), ve_mse = mean(ve_mse), naive_mse = mean(naive_mse),
                  mspe = mean(mspe), ve_mspe = mean(ve_mspe), naive_mspe = mean(naive_mspe),
                  cover = mean(cover), ve_cover = mean(ve_cover),
                  wald_pred_cover = mean(wald_pred_cover), ve_wald_pred_cover = mean(ve_wald_pred_cover),
                  quantile_cred_pred_cover = mean(quantile_cred_pred_cover), ve_quantile_cred_pred_cover = mean(ve_quantile_cred_pred_cover),
                  quantile_samp_pred_cover = mean(quantile_samp_pred_cover), ve_quantile_samp_pred_cover = mean(ve_quantile_samp_pred_cover),
                  wald_width = mean(wald_width), ve_wald_width = mean(ve_wald_width),
                  quantile_cred_width = mean(quantile_cred_width), ve_quantile_cred_width = mean(ve_quantile_cred_width),
                  quantile_samp_width = mean(quantile_samp_width), ve_quantile_samp_width = mean(ve_quantile_samp_width),
                  bias = mean(bias), ve_bias = mean(ve_bias), naive_bias = mean(naive_bias), sigma = mean(sigma))
      
      ## (6) average over n, for each q_obs
      average_over_n <- performance_across_taxa %>%
        group_by(q_obs) %>%
        summarize(mse = mean(mse), ve_mse = mean(ve_mse), naive_mse = mean(naive_mse),
                  mspe = mean(mspe), ve_mspe = mean(ve_mspe), naive_mspe = mean(naive_mspe),
                  cover = mean(cover), ve_cover = mean(ve_cover),
                  wald_pred_cover = mean(wald_pred_cover), ve_wald_pred_cover = mean(ve_wald_pred_cover),
                  quantile_cred_pred_cover = mean(quantile_cred_pred_cover), ve_quantile_cred_pred_cover = mean(ve_quantile_cred_pred_cover),
                  quantile_samp_pred_cover = mean(quantile_samp_pred_cover), ve_quantile_samp_pred_cover = mean(ve_quantile_samp_pred_cover),
                  wald_width = mean(wald_width), ve_wald_width = mean(ve_wald_width),
                  quantile_cred_width = mean(quantile_cred_width), ve_quantile_cred_width = mean(ve_quantile_cred_width),
                  quantile_samp_width = mean(quantile_samp_width), ve_quantile_samp_width = mean(ve_quantile_samp_width),
                  bias = mean(bias), ve_bias = mean(ve_bias), naive_bias = mean(naive_bias), sigma = mean(sigma))
      average_over_n_single_taxon <- mc_averaged_performance %>%
        filter(taxon_id == 10) %>%
        group_by(q_obs) %>%
        summarize(mse = mean(mse), ve_mse = mean(ve_mse), naive_mse = mean(naive_mse),
                  mspe = mean(mspe), ve_mspe = mean(ve_mspe), naive_mspe = mean(naive_mspe),
                  cover = mean(cover), ve_cover = mean(ve_cover),
                  wald_pred_cover = mean(wald_pred_cover), ve_wald_pred_cover = mean(ve_wald_pred_cover),
                  quantile_cred_pred_cover = mean(quantile_cred_pred_cover), ve_quantile_cred_pred_cover = mean(ve_quantile_cred_pred_cover),
                  quantile_samp_pred_cover = mean(quantile_samp_pred_cover), ve_quantile_samp_pred_cover = mean(ve_quantile_samp_pred_cover),
                  wald_width = mean(wald_width), ve_wald_width = mean(ve_wald_width),
                  quantile_cred_width = mean(quantile_cred_width), ve_quantile_cred_width = mean(ve_quantile_cred_width),
                  quantile_samp_width = mean(quantile_samp_width), ve_quantile_samp_width = mean(ve_quantile_samp_width),
                  bias = mean(bias), ve_bias = mean(ve_bias), naive_bias = mean(naive_bias), sigma = mean(sigma))
      
      ## (5) add rmse, transpose
      performance_matrix <- average_over_n %>%
        mutate(rmse = sqrt(mse), ve_rmse = sqrt(ve_mse), naive_rmse = sqrt(naive_mse),
               rmspe = sqrt(mspe), ve_rmspe = sqrt(ve_mspe), naive_rmspe = sqrt(naive_mspe))
      performance_matrix_single_taxon <- average_over_n_single_taxon %>%
        mutate(rmse = sqrt(mse), ve_rmse = sqrt(ve_mse), naive_rmse = sqrt(naive_mse),
               rmspe = sqrt(mspe), ve_rmspe = sqrt(ve_mspe), naive_rmspe = sqrt(naive_mspe))
      
      ## (6) average only over n, taxa
      average_over_n_taxa <- performance_df %>%
        filter(include_taxa_in_avg) %>% 
        select(q, q_obs, subj_id, taxon_id, mc_id, mse, ve_mse, naive_mse, mspe, ve_mspe, naive_mspe, 
               cover, ve_cover, wald_pred_cover, ve_wald_pred_cover, quantile_cred_pred_cover, ve_quantile_cred_pred_cover,
               quantile_samp_pred_cover, ve_quantile_samp_pred_cover, wald_width, ve_wald_width, quantile_cred_width, ve_quantile_cred_width,
               quantile_samp_width, ve_quantile_samp_width,
               bias, ve_bias, naive_bias, sigma) %>%
        group_by(q, q_obs, mc_id) %>%
        summarize(mse = mean(mse, na.rm = TRUE), ve_mse = mean(ve_mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE), 
                  mspe = mean(mspe, na.rm = TRUE), ve_mspe = mean(ve_mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE), 
                  cover = mean(cover, na.rm = TRUE), ve_cover = mean(ve_cover, na.rm = TRUE),
                  wald_pred_cover = mean(wald_pred_cover, na.rm = TRUE), ve_wald_pred_cover = mean(ve_wald_pred_cover, na.rm = TRUE),
                  quantile_cred_pred_cover = mean(quantile_cred_pred_cover, na.rm = TRUE), ve_quantile_cred_pred_cover = mean(ve_quantile_cred_pred_cover, na.rm = TRUE),
                  quantile_samp_pred_cover = mean(quantile_samp_pred_cover, na.rm = TRUE), ve_quantile_samp_pred_cover = mean(ve_quantile_samp_pred_cover, na.rm = TRUE),
                  wald_width = mean(wald_width, na.rm = TRUE), ve_wald_width = mean(ve_wald_width, na.rm = TRUE),
                  quantile_cred_width = mean(quantile_cred_width, na.rm = TRUE), ve_quantile_cred_width = mean(ve_quantile_cred_width, na.rm = TRUE),
                  quantile_samp_width = mean(quantile_samp_width, na.rm = TRUE), ve_quantile_samp_width = mean(ve_quantile_samp_width, na.rm = TRUE),
                  bias = mean(bias, na.rm = TRUE), ve_bias = mean(ve_bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE), sigma = mean(sigma)) %>%
        ungroup()
      
      ## (6a) average only over taxa, for e
      e_average_over_taxa <- e_performance_df %>%
        select(q, q_obs, taxon_id, mc_id, mse, sigma_mse, cover, sigma_cover, bias, sigma_bias, sigma) %>%
        group_by(q, q_obs, mc_id) %>%
        summarize(mse = mean(mse), sigma_mse = mean(sigma_mse), 
                  cover = mean(cover), sigma_cover = mean(sigma_cover),
                  bias = mean(bias), sigma_bias = mean(sigma_bias), sigma = mean(sigma)) %>% 
        ungroup()
      
      output_performances[[i]] <- performance_matrix
      output_performances_single_taxon[[i]] <- performance_matrix_single_taxon
      output_performances_avg_over_taxa_n[[i]] <- average_over_n_taxa
      
      output_performances_e[[i]] <- e_mc_averaged_performance
      output_performances_e_avg_over_taxa[[i]] <- e_average_over_taxa
    }
    saveRDS(output_performances, paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
    saveRDS(output_performances_single_taxon, paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_single_taxon_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
    saveRDS(output_performances_avg_over_taxa_n, paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
    
    saveRDS(output_performances_e, paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_e_n_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
    saveRDS(output_performances_e_avg_over_taxa, paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_e_avg_over_taxa_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
  }
} else {
  
}
out_perf_lst <- vector("list", length = length(args$q_obs))
out_perf_single_lst <- out_perf_lst
out_perf_avg_lst <- out_perf_lst
out_perf_e <- out_perf_lst
out_perf_e_avg <- out_perf_lst
for (j in 1:length(args$q_obs)) {
  out_perf_lst[[j]] <- readRDS(paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
  out_perf_single_lst[[j]] <- readRDS(paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_single_taxon_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
  out_perf_avg_lst[[j]] <- readRDS(paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
  
  out_perf_e[[j]] <- readRDS(paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_e_n_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
  out_perf_e_avg[[j]] <- readRDS(paste0("results/", args$sim_name, "/cov_", args$corr, "/n_", args$N, "/output_performances_e_avg_over_taxa_ab_", args$most_abundant, "_q_obs", args$q_obs[j], ".rds"))
}
output_performances <- do.call(c, out_perf_lst)
output_performances_avg_over_taxa_n <- do.call(c, out_perf_avg_lst)
output_performances_e <- do.call(c, out_perf_e)
output_performances_e_avg_over_taxa <- do.call(c, out_perf_e_avg)

## transform into matrices
performance_df <- do.call(rbind.data.frame, output_performances)
performance_avg_over_taxa_n <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n)
performance_e <- do.call(rbind.data.frame, output_performances_e)
performance_e_avg_over_taxa <- do.call(rbind.data.frame, output_performances_e_avg_over_taxa)

## add new stuff
e_avgs <- lapply(output_performances_e[1:12], colMeans)
e_avgs_2 <- lapply(e_avgs, function(x) data.frame(t(x)))

## create new tibble with rmspe, cover, pred_cover and estimator type
performance_by_type <- performance_avg_over_taxa_n %>% 
  mutate(naive_cover = 999, naive_pred_cover = 999, naive_width = -999,
         no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_wald_pred_cover = wald_pred_cover,
         no_ve_quantile_cred_pred_cover = quantile_cred_pred_cover, no_ve_quantile_samp_pred_cover = quantile_samp_pred_cover,
         no_ve_wald_width = wald_width, no_ve_quantile_cred_width = quantile_cred_width, no_ve_quantile_samp_width = quantile_samp_width,
         no_ve_bias = bias) %>%
  select(-mse, -mspe, -cover, -wald_pred_cover, -quantile_cred_pred_cover, -quantile_samp_pred_cover,
         -wald_width, -quantile_cred_width, -quantile_samp_width, -bias) %>% 
  gather(key, value, -q, -q_obs, -mc_id, -sigma) %>% 
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
point_cex <- 2
axis_cex <- 1.75

## --------------------------------------------------------------------------------------------------
## MAIN PLOT: prediction interval coverage, credible interval coverage, rms(p)e
## for q_obs == 7 only
## --------------------------------------------------------------------------------------------------
pred_cover_plot <- performance_by_type %>%
  filter(q_obs == 7) %>% 
  ggplot(aes(x = sigma, y = wald_pred_cover, group = paste(sigma, estimator, sep = "_"),
             fill = factor(estimator, levels = c("ve", "no_ve", "naive"), ordered = FALSE, labels = c("Proposed, ve", "Proposed, no ve", "Naive")))) +
  xlab(expression(sigma[e])) +
  ylab("Coverage") + 
  ggtitle("Prediction interval coverage") + 
  labs(fill = "Estimator type") +
  geom_boxplot(width = 0.08, position = position_dodge(), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.4, 1)) +
  theme_bw() +
  theme(legend.position = c(0.185, 0.2)) +
  guides(fill = FALSE, alpha = FALSE)
cover_plot <- performance_by_type %>%
  filter(q_obs == 7) %>% 
  ggplot(aes(x = sigma, y = cover, group = paste(sigma, estimator, sep = "_"), 
             fill = factor(estimator, levels = c("ve", "no_ve", "naive"), ordered = FALSE, labels = c("Proposed, ve", "Proposed, no ve", "Naive")))) +
  xlab(expression(sigma[e])) +
  ylab("Coverage") + 
  ggtitle("Credible interval coverage") + 
  labs(fill = "Estimator type") +
  geom_boxplot(width = 0.08, position = position_dodge(), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.4, 1)) +
  theme_bw() +
  theme(legend.position = c(0.185, 0.2)) +
  guides(fill = FALSE, alpha = FALSE)

rmspe_plot <- performance_by_type %>%
  filter(q_obs == 7) %>%
  ggplot(aes(x = sigma, y = rmspe, group = paste(sigma, estimator, sep = "_"), 
             fill = factor(estimator, levels = c("ve", "no_ve", "naive"), ordered = FALSE, labels = c("Proposed, ve", "Proposed, no ve", "Naive")))) +
  xlab(expression(sigma[e])) +
  ylab("Root mean squared prediction error") + 
  ggtitle("RMSPE") + 
  labs(fill = "Estimator type") +
  geom_boxplot(width = 0.08, position = position_dodge(), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 200)) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8))

rmse_plot <- performance_by_type %>%
  filter(q_obs == 7) %>%
  ggplot(aes(x = sigma, y = rmse, group = paste(sigma, estimator, sep = "_"), 
             fill = factor(estimator, levels = c("ve", "no_ve", "naive"), ordered = FALSE, labels = c("Proposed, ve", "Proposed, no ve", "Naive")))) +
  xlab(expression(sigma[e])) +
  ylab("Root mean squared error") + 
  ggtitle("RMSE") + 
  labs(fill = "Estimator type") +
  geom_boxplot(width = 0.08, position = position_dodge(), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 200)) +
  guides(fill = FALSE) +
  theme_bw() 

png(paste0(plots_dir, "vary_sigma_e.png"), width = fig_width, height = fig_height, units = "px", res = 300)
gridExtra::grid.arrange(pred_cover_plot, cover_plot, 
                        rmspe_plot, rmse_plot)
dev.off()

## --------------------------------------------------------------------------------------------------
## SUPPLEMENTARY PLOTS:
## (a) same as main plot, but for q_obs == 3
## (b) widths for both q_obs
## (c) squared bias
## (d) efficiencies?
## --------------------------------------------------------------------------------------------------
pred_cover_plot_3 <- performance_by_type %>%
  filter(q_obs == 3) %>% 
  ggplot(aes(x = sigma, y = wald_pred_cover, group = paste(sigma, estimator, sep = "_"),
             fill = factor(estimator, levels = c("ve", "no_ve", "naive"), ordered = FALSE, labels = c("Proposed, ve", "Proposed, no ve", "Naive")))) +
  xlab(expression(sigma[e])) +
  ylab("Coverage") + 
  ggtitle("Prediction interval coverage") + 
  labs(fill = "Estimator type") +
  geom_boxplot(width = 0.08, position = position_dodge(), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.4, 1)) +
  theme_bw() +
  theme(legend.position = c(0.185, 0.2)) +
  guides(fill = FALSE, alpha = FALSE)
cover_plot_3 <- performance_by_type %>%
  filter(q_obs == 3) %>% 
  ggplot(aes(x = sigma, y = cover, group = paste(sigma, estimator, sep = "_"), 
             fill = factor(estimator, levels = c("ve", "no_ve", "naive"), ordered = FALSE, labels = c("Proposed, ve", "Proposed, no ve", "Naive")))) +
  xlab(expression(sigma[e])) +
  ylab("Coverage") + 
  ggtitle("Credible interval coverage") + 
  labs(fill = "Estimator type") +
  geom_boxplot(width = 0.08, position = position_dodge(), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.4, 1)) +
  theme_bw() +
  theme(legend.position = c(0.185, 0.2)) +
  guides(fill = FALSE, alpha = FALSE)

rmspe_plot_3 <- performance_by_type %>%
  filter(q_obs == 3) %>%
  ggplot(aes(x = sigma, y = rmspe, group = paste(sigma, estimator, sep = "_"), 
             fill = factor(estimator, levels = c("ve", "no_ve", "naive"), ordered = FALSE, labels = c("Proposed, ve", "Proposed, no ve", "Naive")))) +
  xlab(expression(sigma[e])) +
  ylab("Root mean squared prediction error") + 
  ggtitle("RMSPE") + 
  labs(fill = "Estimator type") +
  geom_boxplot(width = 0.08, position = position_dodge(), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 200)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.3))

rmse_plot_3 <- performance_by_type %>%
  filter(q_obs == 3) %>%
  ggplot(aes(x = sigma, y = rmse, group = paste(sigma, estimator, sep = "_"), 
             fill = factor(estimator, levels = c("ve", "no_ve", "naive"), ordered = FALSE, labels = c("Proposed, ve", "Proposed, no ve", "Naive")))) +
  xlab(expression(sigma[e])) +
  ylab("Root mean squared error") + 
  ggtitle("RMSE") + 
  labs(fill = "Estimator type") +
  geom_boxplot(width = 0.08, position = position_dodge(), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 200)) +
  guides(fill = FALSE) +
  theme_bw() 

png(paste0(plots_dir, "vary_sigma_e_q_obs_3.png"), width = fig_width, height = fig_height, units = "px", res = 300)
gridExtra::grid.arrange(pred_cover_plot_3, cover_plot_3, 
                        rmspe_plot_3, rmse_plot_3)
dev.off()

width_plot <- performance_by_type %>%
  gather(key = "width_type", value = "width", width, wald_width, quantile_cred_width, quantile_samp_width) %>% 
  select(q_obs, mc_id, sigma, estimator, width_type, width) %>% 
  ggplot(aes(x = sigma, y = width, group = paste(q_obs, sigma, estimator, width_type, sep = "_"),
             alpha = factor(q_obs),
             fill = factor(estimator, levels = c("ve", "no_ve", "naive"), ordered = FALSE, labels = c("Proposed, ve", "Proposed, no ve", "Naive")),
             color = factor(width_type, levels = c("quantile_cred_width", "wald_width", "quantile_samp_width", "width"), ordered = FALSE, labels = c("Quantile-cred.", "Wald-type", "Quantile-samp.", "Pred.")))) +
  xlab(expression(sigma[e])) +
  ylab("Prediction interval width") + 
  ggtitle("Prediction interval width") + 
  labs(fill = "Estimator type") +
  labs(color = "Interval type") +
  labs(alpha = expression(q^obs)) +
  geom_boxplot(width = 0.08, position = position_dodge(), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols[c(3, 5)]) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 700)) +
  guides(alpha = guide_legend(override.aes = list(fill = "blue"))) +
  theme_bw() 
png(paste0(plots_dir, "vary_sigma_e_width.png"), width = fig_width, height = fig_height, units = "px", res = 300)
width_plot
dev.off()

## mse as bias^2 + variance
bias_plot <- performance_by_type %>%
  ggplot(aes(x = sigma, y = sq_bias, group = paste(q_obs, sigma, estimator, sep = "_"), 
             fill = factor(estimator, levels = c("ve", "no_ve", "naive"), ordered = FALSE, labels = c("Proposed, ve", "Proposed, no ve", "Naive")))) +
  xlab(expression(sigma[e])) +
  ylab("Squared bias") + 
  ggtitle("Squared bias") + 
  labs(fill = "Estimator type") +
  geom_boxplot(width = 0.08, position = position_dodge(), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 200)) +
  theme_bw() 

png(paste0(plots_dir, "vary_sigma_e_bias.png"), width = fig_width, height = fig_height, units = "px", res = 300)
bias_plot
dev.off()

## ---------------------------------------------------------------------
## distribution of estimated efficiencies when efficiency truly is zero!
## ---------------------------------------------------------------------
results_dir <- paste0("results/", args$sim_name, "/ve_", args$sigma[1], "/cov_", args$corr, "/n_", args$N, "/q_", args$q, "/q_obs_", args$q_obs[2], "/")

dir_mat <- expand.grid(job = 1:args$num_jobs, dir = results_dir, stan_model = args$stan_model)
mod_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_mod_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
data_nms_lst <- as.list(paste0(dir_mat$dir, dir_mat$stan_model, "_data_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))

mod_lst <- lapply(mod_nms_lst, read_func)
data_lst <- lapply(data_nms_lst, read_func)

## pick out the ones where I've estimated varying efficiency; get MC reps for all taxa
get_e_rows <- function(x) grepl("e", rownames(x)) & !grepl("beta", rownames(x))
get_sigma_row <- function(x) grepl("sigma", rownames(x))
e_lst <- lapply(mod_lst[(args$num_jobs + 1):(2*args$num_jobs)], function(x) rbind(x[get_e_rows(x), ], x[get_sigma_row(x), ]))

## pick out the means for both e and sigma
e_df <- do.call(rbind.data.frame, lapply(e_lst, function(x) x[, "mean"]))
colnames(e_df) <- c(paste0("taxon_", 1:args$q), "sigma")
e_df$mc_id <- 1:args$num_jobs

e_df %>% 
  ggplot(aes(y = taxon_1))  +
  geom_histogram() 

png(paste0(plots_dir, "vary_sigma_e_sigma.png"), width = fig_width, height = fig_height, units = "px", res = 300)
e_df %>% 
  ggplot(aes(y = sigma)) +
  geom_boxplot() +
  ggtitle(expression(paste("Estimated ", sigma[e], " when true ", sigma[e], " = 0"))) +
  labs(y = expression(hat(sigma)[e])) +
  theme_bw()
dev.off()
