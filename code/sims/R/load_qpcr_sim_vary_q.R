#!/usr/local/bin/Rscript
##################################################################################
## FILE: load_qpcr_sim_vary_q.R
##
## CREATED: 05 November 2018 by Brian Williamson
##
## PURPOSE: load in results from a qPCR simulation varying q, q_obs
##          and create plots, etc.
##
## INPUTS: listed below
##
## OUTPUTS: plots, tables, etc.
##################################

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
library("paramedic")
source("load_qpcr_sim_helpers.R") # provides function get_summaries
source("code/R/naive_qpcr_estimator.R")
source("code/R/gen_naive_interval.R")


## grab command-line arguments
parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "vary_q", help = "name of the simulation")
parser$add_argument("--stan-model", default = "predict_qpcr_noncentered",
                    help = "Which Stan model file to use.")
parser$add_argument("--N", type = "double", default = 100, help = "sample size")
parser$add_argument("--q", nargs = "+", type = "double", default = c(10, 20, 40, 60), help = "how many taxa do we have?")
parser$add_argument("--corr", type = "double", default = 0, help = "hyperparameter controlling the off-diagonal elements of Sigma.")
parser$add_argument("--sigma", type = "double", default = 0, help = "SD of efficiencies (0 = no varying efficiency)")
parser$add_argument("--ad", type = "double", default = 0.850000, help = "adapt delta")
parser$add_argument("--mt", type = "double", default = 15, help = "max treedepth")
parser$add_argument("--num-jobs", type = "double", default = 50, help = "number of jobs run")
parser$add_argument("--taxa", type = "double", default = 10, help = "taxon of interest")
parser$add_argument("--q-obs", nargs = "+", type = "double", default = c(2, 3, 4, 5, 6, 7), help = "number of taxa with observed qPCR")
parser$add_argument("--most-abundant", type = "double", default = 1, help = "whether or not to use most abundant taxa")
parser$add_argument("--read-data", type = "double", default = 0, help = "whether or not to read in raw output, or instead do summaries")
args <- parser$parse_args()

read_func <- function(x) tryCatch(readRDS(x), error = function(e) NA)

plots_dir <- paste0("plots/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

## --------------------------------------------------------------------------------------------------
## CREATE LIST OF OUTPUT FOR EASY SUMMARIES; GET THOSE SUMMARIES
## --------------------------------------------------------------------------------------------------
if (args$read_data) {
  output_performances <- vector(mode = "list", length = length(args$q))
  output_performances_single_taxon <- vector(mode = "list", length = length(args$q))
  
  output_performances_avg_over_taxa_n <- output_performances_avg_over_taxa_n_mn <- output_performances_avg_over_taxa_n_var <- output_performances_avg_over_taxa_n_nofilter <- vector(mode = "list", length = length(args$q))
  output_performances_avg_over_n <- output_performances_avg_over_n_mn <- output_performances_avg_over_n_var <- output_performances_avg_over_n_nofilter <- vector(mode = "list", length = length(args$q))
  for (i in 1:length(args$q)) {
    ## load in results and data
    results_dir <- paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/q_", args$q[i], "/q_obs_", args$q_obs, "/")
    
    dir_mat <- expand.grid(job = 1:args$num_jobs, dir = results_dir)
    mod_nms_lst <- as.list(paste0(dir_mat$dir, args$stan_model, "_mod_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
    data_nms_lst <- as.list(paste0(dir_mat$dir, args$stan_model, "_data_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
    samp_nms_lst <- as.list(paste0(dir_mat$dir, args$stan_model, "_samps_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
    
    mod_lst <- lapply(mod_nms_lst, read_func)
    data_lst <- lapply(data_nms_lst, read_func)
    # samps_lst <- lapply(samp_nms_lst, read_func)
    
    ## (1) pair model summaries with relevant data, for all taxa, all q_obs; replace NA with samps_lst if you want
    summary_df <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type) get_summaries(w, x, y, z, type), mod_lst, data_lst, NA, MoreArgs = list(z = 1:args$q[i], type = "no_ve"), SIMPLIFY = FALSE))
    ## pair with the monte-carlo id
    summary_df$mc_id <- rep(NA, dim(summary_df)[1])
    mc_lst_2 <- mc_lst_3 <- mc_lst_4 <- mc_lst_5 <- mc_lst_6 <- mc_lst_7 <- 1:50
    
    summary_df$mc_id[!is.na(summary_df$q)] <- c(rep(mc_lst_2, each = args$N*args$q[i]), rep(mc_lst_3, each = args$N*args$q[i]),
                                                rep(mc_lst_4, each = args$N*args$q[i]), rep(mc_lst_5, each = args$N*args$q[i]),
                                                rep(mc_lst_6, each = args$N*args$q[i]), rep(mc_lst_7, each = args$N*args$q[i]))
    
    
    ## (2) compute performance for each row
    performance_df <- summary_df %>%
      mutate(mse = (mu - est)^2, naive_mse = (mu - naive_est)^2,
             mspe = (qpcr - est)^2, naive_mspe = (qpcr - naive_est)^2,
             cover = cil <= mu & ciu >= mu,
             wald_cover = wald_cred_cil <= mu & wald_cred_ciu >= mu,
             naive_cover = naive_cil <= mu & naive_ciu >= mu,
             pred_cover = wald_pred_cil <= qpcr & wald_pred_ciu >= qpcr, # choose wald
             naive_pred_cover = naive_wald_pred_cil <= qpcr & naive_wald_pred_ciu >= qpcr,
             width = abs(wald_pred_ciu - wald_pred_cil),
             bias = mu - est, naive_bias = mu - naive_est)
    
    performance_df$include_taxa_in_v_avg <- performance_df$taxon_id > performance_df$q_obs
    
    ## (3) average over MC reps for each taxon
    mc_averaged_performance <- performance_df %>%
      filter(include_taxa_in_avg) %>%
      select(q, q_obs, subj_id, taxon_id, mse, naive_mse, mspe, naive_mspe, cover, naive_cover, wald_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, subj_id, taxon_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE),
                mspe = mean(mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE),
                pred_cover = mean(pred_cover, na.rm = TRUE), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                width = mean(width, na.rm = TRUE), bias = mean(bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE)) %>%
      ungroup()
    mc_averaged_performance_nofilter <- performance_df %>%
      select(q, q_obs, subj_id, taxon_id, mse, naive_mse, mspe, naive_mspe, cover, naive_cover, wald_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, subj_id, taxon_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE),
                mspe = mean(mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE),
                pred_cover = mean(pred_cover, na.rm = TRUE), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                width = mean(width, na.rm = TRUE), bias = mean(bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE)) %>%
      ungroup()
    
    ## (4) set flag for taxa of interest: for V, need taxon_id > q^obs
    mc_averaged_performance$include_taxa_in_v_avg <- apply(mc_averaged_performance, 1,
                                                           function(x) ifelse(!is.na(x[1]), x[4] %in% tail(1:x[1], x[1] - x[2]), NA))
    mc_averaged_performance_nofilter$include_taxa_in_v_avg <- apply(mc_averaged_performance_nofilter, 1,
                                                                    function(x) ifelse(!is.na(x[1]), x[4] %in% tail(1:x[1], x[1] - x[2]), NA))
    ## (5) average over taxa of interest for each q_obs
    performance_across_taxa_mu <- mc_averaged_performance %>%
      group_by(q, q_obs, subj_id) %>%
      summarize(mse = mean(mse), naive_mse = mean(naive_mse),
                cover = mean(cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE),
                bias = mean(bias), naive_bias = mean(naive_bias))
    performance_across_taxa_v <- mc_averaged_performance %>%
      group_by(q, q_obs, subj_id) %>%
      filter(include_taxa_in_v_avg) %>%
      summarize(mspe = mean(mspe), naive_mspe = mean(naive_mspe),
                pred_cover = mean(pred_cover), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                width = mean(width))
    performance_across_taxa <- left_join(performance_across_taxa_mu, performance_across_taxa_v, by = c("q", "q_obs", "subj_id"))
    ## without filter
    performance_across_taxa_mu_nofilter <- mc_averaged_performance_nofilter %>%
      group_by(q, q_obs, subj_id) %>%
      summarize(mse = mean(mse), naive_mse = mean(naive_mse),
                cover = mean(cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE),
                bias = mean(bias), naive_bias = mean(naive_bias))
    performance_across_taxa_v_nofilter <- mc_averaged_performance_nofilter %>%
      group_by(q, q_obs, subj_id) %>%
      filter(include_taxa_in_v_avg) %>%
      summarize(mspe = mean(mspe), naive_mspe = mean(naive_mspe),
                pred_cover = mean(pred_cover), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                width = mean(width))
    performance_across_taxa_nofilter <- left_join(performance_across_taxa_mu_nofilter, performance_across_taxa_v_nofilter, by = c("q", "q_obs", "subj_id"))
    
    ## (6) average over n, for each q_obs
    average_over_n <- performance_across_taxa %>%
      group_by(q_obs) %>%
      summarize(mse = mean(mse), naive_mse = mean(naive_mse),
                mspe = mean(mspe), naive_mspe = mean(naive_mspe),
                cover = mean(cover), naive_cover = mean(naive_cover), wald_cover = mean(wald_cover),
                pred_cover = mean(pred_cover), naive_pred_cover = mean(naive_pred_cover),
                width = mean(width), bias = mean(bias), naive_bias = mean(naive_bias))
    average_over_n_single_taxon <- mc_averaged_performance %>%
      filter(taxon_id == 10) %>%
      group_by(q_obs) %>%
      summarize(mse = mean(mse), naive_mse = mean(naive_mse),
                mspe = mean(mspe), naive_mspe = mean(naive_mspe),
                cover = mean(cover), naive_cover = mean(naive_cover), wald_cover = mean(wald_cover),
                pred_cover = mean(pred_cover), naive_pred_cover = mean(naive_pred_cover),
                width = mean(width), bias = mean(bias), naive_bias = mean(naive_bias))
    
    ## no filter
    average_over_n_nofilter <- performance_across_taxa_nofilter %>%
      group_by(q_obs) %>%
      summarize(mse = mean(mse), naive_mse = mean(naive_mse),
                mspe = mean(mspe), naive_mspe = mean(naive_mspe),
                cover = mean(cover), naive_cover = mean(naive_cover), wald_cover = mean(wald_cover),
                pred_cover = mean(pred_cover), naive_pred_cover = mean(naive_pred_cover),
                width = mean(width), bias = mean(bias), naive_bias = mean(naive_bias))
    average_over_n_single_taxon_nofilter <- mc_averaged_performance_nofilter %>%
      filter(taxon_id == 10) %>%
      group_by(q_obs) %>%
      summarize(mse = mean(mse), naive_mse = mean(naive_mse),
                mspe = mean(mspe), naive_mspe = mean(naive_mspe),
                cover = mean(cover), naive_cover = mean(naive_cover), wald_cover = mean(wald_cover),
                pred_cover = mean(pred_cover), naive_pred_cover = mean(naive_pred_cover),
                width = mean(width), bias = mean(bias), naive_bias = mean(naive_bias))
    
    ## (5) add rmse, transpose
    performance_matrix <- average_over_n %>%
      mutate(rmse = sqrt(mse), naive_rmse = sqrt(naive_mse),
             rmspe = sqrt(mspe), naive_rmspe = sqrt(naive_mspe))
    performance_matrix_single_taxon <- average_over_n_single_taxon %>%
      mutate(rmse = sqrt(mse), naive_rmse = sqrt(naive_mse),
             rmspe = sqrt(mspe), naive_rmspe = sqrt(naive_mspe))
    
    performance_matrix_nofilter <- average_over_n_nofilter %>%
      mutate(rmse = sqrt(mse), naive_rmse = sqrt(naive_mse),
             rmspe = sqrt(mspe), naive_rmspe = sqrt(naive_mspe))
    performance_matrix_single_taxon_nofilter <- average_over_n_single_taxon_nofilter %>%
      mutate(rmse = sqrt(mse), naive_rmse = sqrt(naive_mse),
             rmspe = sqrt(mspe), naive_rmspe = sqrt(naive_mspe))
    
    ## (6) average only over n, taxa
    ## no filter
    average_over_n_taxa_mu_nofilter <- performance_df %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, naive_cover, wald_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE),
                bias = mean(bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE)) %>%
      ungroup()
    
    average_over_n_taxa_v_nofilter <- performance_df %>%
      filter(include_taxa_in_v_avg) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mspe = mean(mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE),
                pred_cover = mean(pred_cover, na.rm = TRUE), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                width = mean(width, na.rm = TRUE)) %>%
      ungroup()
    average_over_n_taxa_nofilter <- left_join(average_over_n_taxa_mu_nofilter, average_over_n_taxa_v_nofilter, by = c("q", "q_obs", "mc_id"))
    
    ## (a) filter by colSums(W) > 0
    average_over_n_taxa_mu <- performance_df %>%
      filter(include_taxa_in_avg) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, naive_cover, wald_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE),
                bias = mean(bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE)) %>%
      ungroup()
    
    average_over_n_taxa_v <- performance_df %>%
      filter(include_taxa_in_v_avg, include_taxa_in_avg) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mspe = mean(mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE),
                pred_cover = mean(pred_cover, na.rm = TRUE), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                width = mean(width, na.rm = TRUE)) %>%
      ungroup()
    average_over_n_taxa <- left_join(average_over_n_taxa_mu, average_over_n_taxa_v, by = c("q", "q_obs", "mc_id"))
    
    ## (b) filter by colMeans(W) > 0.5
    average_over_n_taxa_mu_mn <- performance_df %>%
      filter(include_taxa_in_avg_w_mn) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, naive_cover, wald_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover),
                bias = mean(bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE)) %>%
      ungroup()
    
    average_over_n_taxa_v_mn <- performance_df %>%
      filter(include_taxa_in_v_avg, include_taxa_in_avg_w_mn) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mspe = mean(mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE),
                pred_cover = mean(pred_cover, na.rm = TRUE), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                width = mean(width, na.rm = TRUE)) %>%
      ungroup()
    average_over_n_taxa_mn <- left_join(average_over_n_taxa_mu_mn, average_over_n_taxa_v_mn, by = c("q", "q_obs", "mc_id"))
    
    ## (c) filter by colVars(W) > 1
    average_over_n_taxa_mu_var <- performance_df %>%
      filter(include_taxa_in_avg_var) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, naive_cover, wald_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE),
                bias = mean(bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE)) %>%
      ungroup()
    
    average_over_n_taxa_v_var <- performance_df %>%
      filter(include_taxa_in_v_avg, include_taxa_in_avg_var) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, mc_id) %>%
      summarize(mspe = mean(mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE),
                pred_cover = mean(pred_cover, na.rm = TRUE), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                width = mean(width, na.rm = TRUE)) %>%
      ungroup()
    average_over_n_taxa_var <- left_join(average_over_n_taxa_mu_var, average_over_n_taxa_v_var, by = c("q", "q_obs", "mc_id"))
    
    ## (7) average only over n
    ## (a) filter by colSums(W) > 0
    average_over_n_mu <- performance_df %>%
      filter(include_taxa_in_avg) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, naive_cover, wald_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, taxon_id, mc_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE),
                bias = mean(bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE)) %>%
      ungroup()
    average_over_n_v <- performance_df %>%
      filter(include_taxa_in_v_avg, include_taxa_in_avg) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, naive_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, taxon_id, mc_id) %>%
      summarize(mspe = mean(mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE),
                pred_cover = mean(pred_cover, na.rm = TRUE), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                width = mean(width, na.rm = TRUE)) %>%
      ungroup()
    average_over_n <- left_join(average_over_n_mu, average_over_n_v, by = c("q", "q_obs", "taxon_id", "mc_id"))
    
    ## (b) filter by colMeans(W) > 0
    average_over_n_mu_mn <- performance_df %>%
      filter(include_taxa_in_avg_w_mn) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, naive_cover, wald_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, taxon_id, mc_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE),
                bias = mean(bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE)) %>%
      ungroup()
    average_over_n_v_mn <- performance_df %>%
      filter(include_taxa_in_v_avg, include_taxa_in_avg_w_mn) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, naive_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, taxon_id, mc_id) %>%
      summarize(mspe = mean(mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE),
                pred_cover = mean(pred_cover, na.rm = TRUE), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                width = mean(width, na.rm = TRUE)) %>%
      ungroup()
    average_over_n_mn <- left_join(average_over_n_mu_mn, average_over_n_v_mn, by = c("q", "q_obs", "taxon_id", "mc_id"))
    
    ## (a) filter by colVars(W) > 0
    average_over_n_mu_var <- performance_df %>%
      filter(include_taxa_in_avg_var) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, naive_cover, wald_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, taxon_id, mc_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE),
                bias = mean(bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE)) %>%
      ungroup()
    average_over_n_v_var <- performance_df %>%
      filter(include_taxa_in_v_avg, include_taxa_in_avg_var) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, naive_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, taxon_id, mc_id) %>%
      summarize(mspe = mean(mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE),
                pred_cover = mean(pred_cover, na.rm = TRUE), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                width = mean(width, na.rm = TRUE)) %>%
      ungroup()
    average_over_n_var <- left_join(average_over_n_mu_var, average_over_n_v_var, by = c("q", "q_obs", "taxon_id", "mc_id"))
    ## no filter
    average_over_n_mu_nofilter <- performance_df %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, naive_cover, wald_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, taxon_id, mc_id) %>%
      summarize(mse = mean(mse, na.rm = TRUE), naive_mse = mean(naive_mse, na.rm = TRUE),
                cover = mean(cover, na.rm = TRUE), naive_cover = mean(naive_cover, na.rm = TRUE), wald_cover = mean(wald_cover, na.rm = TRUE),
                bias = mean(bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE)) %>%
      ungroup()
    average_over_n_v_nofilter <- performance_df %>%
      filter(include_taxa_in_v_avg) %>%
      select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
             cover, naive_cover, pred_cover, naive_pred_cover, width, bias, naive_bias) %>%
      group_by(q, q_obs, taxon_id, mc_id) %>%
      summarize(mspe = mean(mspe, na.rm = TRUE), naive_mspe = mean(naive_mspe, na.rm = TRUE),
                pred_cover = mean(pred_cover, na.rm = TRUE), naive_pred_cover = mean(naive_pred_cover, na.rm = TRUE),
                width = mean(width, na.rm = TRUE)) %>%
      ungroup()
    average_over_n_nofilter <- left_join(average_over_n_mu_nofilter, average_over_n_v_nofilter, by = c("q", "q_obs", "taxon_id", "mc_id"))
    
    output_performances[[i]] <- performance_matrix
    output_performances_single_taxon[[i]] <- performance_matrix_single_taxon
    output_performances_avg_over_taxa_n[[i]] <- average_over_n_taxa
    output_performances_avg_over_n[[i]] <- average_over_n
    output_performances_avg_over_taxa_n_mn[[i]] <- average_over_n_taxa_mn
    output_performances_avg_over_n_mn[[i]] <- average_over_n_mn
    output_performances_avg_over_taxa_n_var[[i]] <- average_over_n_taxa_var
    output_performances_avg_over_n_var[[i]] <- average_over_n_var
    output_performances_avg_over_taxa_n_nofilter[[i]] <- average_over_n_taxa_nofilter
    output_performances_avg_over_n_nofilter[[i]] <- average_over_n_nofilter
  }
  saveRDS(output_performances, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_ab_", args$most_abundant, ".rds"))
  saveRDS(output_performances_single_taxon, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_single_taxon_ab_", args$most_abundant, ".rds"))
  saveRDS(output_performances_avg_over_taxa_n, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_ab_", args$most_abundant, ".rds"))
  saveRDS(output_performances_avg_over_n, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_n_ab_", args$most_abundant, ".rds"))
  saveRDS(output_performances_avg_over_taxa_n_mn, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_mn_ab_", args$most_abundant, ".rds"))
  saveRDS(output_performances_avg_over_n_mn, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_n_mn_ab_", args$most_abundant, ".rds"))
  saveRDS(output_performances_avg_over_taxa_n_var, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_var_ab_", args$most_abundant, ".rds"))
  saveRDS(output_performances_avg_over_n_var, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_n_var_ab_", args$most_abundant, ".rds"))
  saveRDS(output_performances_avg_over_taxa_n_nofilter, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_nofilter_ab_", args$most_abundant, ".rds"))
  saveRDS(output_performances_avg_over_n_nofilter, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_n_nofilter_ab_", args$most_abundant, ".rds"))
} else {
  output_performances <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_ab_", args$most_abundant, ".rds"))
  output_performances_single_taxon <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_single_taxon_ab_", args$most_abundant, ".rds"))
  output_performances_avg_over_taxa_n <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_ab_", args$most_abundant, ".rds"))
  output_performances_avg_over_n <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_n_ab_", args$most_abundant, ".rds"))
  output_performances_avg_over_taxa_n_mn <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_mn_ab_", args$most_abundant, ".rds"))
  output_performances_avg_over_n_mn <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_n_mn_ab_", args$most_abundant, ".rds"))
  output_performances_avg_over_taxa_n_var <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_var_ab_", args$most_abundant, ".rds"))
  output_performances_avg_over_n_var <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_n_var_ab_", args$most_abundant, ".rds"))
  output_performances_avg_over_taxa_n_nofilter <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_nofilter_ab_", args$most_abundant, ".rds"))
  output_performances_avg_over_n_nofilter <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_n_nofilter_ab_", args$most_abundant, ".rds"))
}

## transform to lists with a matrix for each q_obs (so that it can be plotted on x-axis)
q_obs_performance <- q_to_q_obs(output_performances)
q_obs_performance_single_taxon <- q_to_q_obs(output_performances_single_taxon)

## transform into a giant matrix for the average over n and taxa
performance_avg_over_taxa_n <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n)
performance_avg_over_taxa_n$grouping <- paste(performance_avg_over_taxa_n$q, performance_avg_over_taxa_n$q_obs, sep = "_")

performance_avg_over_taxa_n_mn <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n_mn)
performance_avg_over_taxa_n_mn$grouping <- paste(performance_avg_over_taxa_n_mn$q, performance_avg_over_taxa_n_mn$q_obs, sep = "_")

performance_avg_over_taxa_n_var <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n_var)
performance_avg_over_taxa_n_var$grouping <- paste(performance_avg_over_taxa_n_var$q, performance_avg_over_taxa_n_var$q_obs, sep = "_")

performance_avg_over_taxa_n_nofilter <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n_nofilter)
performance_avg_over_taxa_n_nofilter$grouping <- paste(performance_avg_over_taxa_n_nofilter$q, performance_avg_over_taxa_n_nofilter$q_obs, sep = "_")

## based on colSums(W)
performance_by_type <- performance_avg_over_taxa_n %>%
  filter(!is.na(q)) %>%
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_pred_cover = pred_cover,
         no_ve_wald_cover = wald_cover,
         no_ve_width = width, no_ve_bias = bias) %>%
  select(-mse, -mspe, -cover, -pred_cover, -width, -bias, -wald_cover) %>%
  gather(key, value, -q, -q_obs, -mc_id, -grouping) %>%
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>%
  spread(measure, value) %>%
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))

## based on colMeans(W)
performance_by_type_mn <- performance_avg_over_taxa_n_mn %>%
  filter(!is.na(q)) %>%
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_pred_cover = pred_cover,
         no_ve_wald_cover = wald_cover,
         no_ve_width = width, no_ve_bias = bias) %>%
  select(-mse, -mspe, -cover, -pred_cover, -width, -bias, -wald_cover) %>%
  gather(key, value, -q, -q_obs, -mc_id, -grouping) %>%
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>%
  spread(measure, value) %>%
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))

## based on colVars(W)
performance_by_type_var <- performance_avg_over_taxa_n_var %>%
  filter(!is.na(q)) %>%
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_pred_cover = pred_cover,
         no_ve_wald_cover = wald_cover,
         no_ve_width = width, no_ve_bias = bias) %>%
  select(-mse, -mspe, -cover, -pred_cover, -width, -bias, -wald_cover) %>%
  gather(key, value, -q, -q_obs, -mc_id, -grouping) %>%
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>%
  spread(measure, value) %>%
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))

## no filter
performance_by_type_nofilter <- performance_avg_over_taxa_n_nofilter %>%
  filter(!is.na(q)) %>%
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_pred_cover = pred_cover,
         no_ve_wald_cover = wald_cover,
         no_ve_width = width, no_ve_bias = bias) %>%
  select(-mse, -mspe, -cover, -pred_cover, -width, -bias, -wald_cover) %>%
  gather(key, value, -q, -q_obs, -mc_id, -grouping) %>%
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>%
  spread(measure, value) %>%
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))

## do the same thing for the average only over n
performance_avg_over_n <- do.call(rbind.data.frame, output_performances_avg_over_n)
performance_avg_over_n$grouping <- paste(performance_avg_over_n$q, performance_avg_over_n$q_obs, performance_avg_over_n$taxon_id, sep = "_")

performance_avg_over_n_mn <- do.call(rbind.data.frame, output_performances_avg_over_n_mn)
performance_avg_over_n_mn$grouping <- paste(performance_avg_over_n_mn$q, performance_avg_over_n_mn$q_obs, performance_avg_over_n_mn$taxon_id, sep = "_")

performance_avg_over_n_var <- do.call(rbind.data.frame, output_performances_avg_over_n_var)
performance_avg_over_n_var$grouping <- paste(performance_avg_over_n_var$q, performance_avg_over_n_var$q_obs, performance_avg_over_n_var$taxon_id, sep = "_")

performance_avg_over_n_nofilter <- do.call(rbind.data.frame, output_performances_avg_over_n_nofilter)
performance_avg_over_n_nofilter$grouping <- paste(performance_avg_over_n_nofilter$q, performance_avg_over_n_nofilter$q_obs, performance_avg_over_n_nofilter$taxon_id, sep = "_")

## colSums(W)
performance_by_type_over_n <- performance_avg_over_n %>%
  filter(!is.na(q)) %>%
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_pred_cover = pred_cover,
         no_ve_wald_cover = wald_cover,
         no_ve_width = width, no_ve_bias = bias) %>%
  select(-mse, -mspe, -cover, -pred_cover, -width, -bias, -wald_cover) %>%
  gather(key, value, -q, -q_obs, -mc_id, -taxon_id, -grouping) %>%
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>%
  spread(measure, value) %>%
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))

## colMeans(W)
performance_by_type_over_n_mn <- performance_avg_over_n_mn %>%
  filter(!is.na(q)) %>%
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_pred_cover = pred_cover,
         no_ve_wald_cover = wald_cover,
         no_ve_width = width, no_ve_bias = bias) %>%
  select(-mse, -mspe, -cover, -pred_cover, -width, -bias, -wald_cover) %>%
  gather(key, value, -q, -q_obs, -mc_id, -taxon_id, -grouping) %>%
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>%
  spread(measure, value) %>%
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))

## colVars(W)
performance_by_type_over_n_var <- performance_avg_over_n_var %>%
  filter(!is.na(q)) %>%
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_pred_cover = pred_cover,
         no_ve_wald_cover = wald_cover,
         no_ve_width = width, no_ve_bias = bias) %>%
  select(-mse, -mspe, -cover, -pred_cover, -width, -bias, -wald_cover) %>%
  gather(key, value, -q, -q_obs, -mc_id, -taxon_id, -grouping) %>%
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>%
  spread(measure, value) %>%
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))

## no filter
performance_by_type_over_n_nofilter <- performance_avg_over_n_nofilter %>%
  filter(!is.na(q)) %>%
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_pred_cover = pred_cover,
         no_ve_wald_cover = wald_cover,
         no_ve_width = width, no_ve_bias = bias) %>%
  select(-mse, -mspe, -cover, -pred_cover, -width, -bias, -wald_cover) %>%
  gather(key, value, -q, -q_obs, -mc_id, -taxon_id, -grouping) %>%
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>%
  spread(measure, value) %>%
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))

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
  group_by(q, q_obs, estimator, grouping) %>%
  summarize(mn_pred_cover = mean(pred_cover)) %>%
  ggplot(aes(x = q_obs, y = mn_pred_cover, group = factor(paste(grouping, "_", estimator, sep = "")),
             color = factor(q), shape = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  # fill = factor(q), alpha = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression(Coverage)) +
  ggtitle("Prediction interval coverage for V") +
  # geom_boxplot(width = 1, position = position_dodge(0.8), outlier.shape = NA) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  # scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  # scale_alpha_discrete(range = c(0.5, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.90, 1)) +
  guides(fill = FALSE, shape = FALSE, color = FALSE)
## next, get credible interval coverage
cred_cover <- performance_by_type %>%
  group_by(q, q_obs, estimator, grouping) %>%
  summarize(mn_cover = mean(cover)) %>%
  ggplot(aes(x = q_obs, y = mn_cover, group = factor(paste(grouping, "_", estimator, sep = "")),
             color = factor(q), shape = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression(Coverage)) +
  labs(color = "q") +
  labs(shape = "Estimator") +
  ggtitle(expression(bold(paste("Interval coverage for ", mu, sep = "")))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.90, 1)) +
  guides(alpha = guide_legend(override.aes = list(fill = "blue"))) +
  theme(legend.position = c(0.05, 0.15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing = unit(.1, "cm"),
        legend.box = "horizontal")

## get rmspe
rmspe <- performance_by_type %>%
  group_by(q, q_obs, estimator, grouping) %>%
  summarize(mn_rmspe = mean(rmspe)) %>%
  ggplot(aes(x = q_obs, y = mn_rmspe, group = factor(paste(grouping, "_", estimator, sep = "")),
             color = factor(q), shape = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression(RMSPE)) +
  ggtitle("Root mean squared prediction error") +
  # geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 350)) +
  guides(color = FALSE, shape = FALSE)

## get rmspe
rmse <- performance_by_type %>%
  group_by(q, q_obs, estimator, grouping) %>%
  summarize(mn_rmse = mean(rmse)) %>%
  ggplot(aes(x = q_obs, y = mn_rmse, group = factor(paste(grouping, "_", estimator, sep = "")),
             color = factor(q), shape = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression(RMSE)) +
  ggtitle("Root mean squared error") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 350)) +
  guides(color = FALSE, shape = FALSE)


## actually plot them
ggsave(paste0(plots_dir, "vary_q_n_", args$N, ".png"), plot = plot_grid(cred_cover, pred_cover, rmse, rmspe),
       device = "png", width = 30, height = 25, units = "cm", dpi = 300)

## get widths, for supplement
width <- performance_by_type %>%
  filter(estimator == "no_ve") %>%
  group_by(q, q_obs, grouping) %>%
  summarize(mn_width = mean(width)) %>%
  ggplot(aes(x = q_obs, y = mn_width, group = grouping,
             color = factor(q))) +
  xlab(expression(q^obs)) +
  ylab(expression(Width)) +
  ggtitle("Prediction interval width") +
  # geom_boxplot(width = 0.7, position = position_identity(), outlier.shape = NA) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  labs(color = "q") +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 450))

ggsave(paste0(plots_dir, "vary_q_width_n_", args$N, ".png"), plot = width,
       width = 10, height = 10, units = "cm", dpi = 300)


## -------------------------------------------------------------------------------------------------
## debug: colMeans, colVars of W
## -------------------------------------------------------------------------------------------------
cred_cover_mn <- performance_by_type_mn %>%
  group_by(q, q_obs, estimator, grouping) %>%
  summarize(mn_cover = mean(cover)) %>%
  ggplot(aes(x = q_obs, y = mn_cover, group = factor(paste(grouping, "_", estimator, sep = "")),
             color = factor(q), shape = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression(Coverage)) +
  labs(color = "q") +
  labs(shape = "Estimator") +
  ggtitle("Credible interval coverage: filter by mean(W) > 0.5") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.9, 1)) +
  guides(alpha = guide_legend(override.aes = list(fill = "blue"))) +
  theme(legend.position = c(0.05, 0.15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing = unit(.1, "cm"),
        legend.box = "horizontal")

cred_cover_var <- performance_by_type_var %>%
  group_by(q, q_obs, estimator, grouping) %>%
  summarize(mn_cover = mean(cover)) %>%
  ggplot(aes(x = q_obs, y = mn_cover, group = factor(paste(grouping, "_", estimator, sep = "")),
             color = factor(q), shape = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression(Coverage)) +
  labs(color = "q") +
  labs(shape = "Estimator") +
  ggtitle(expression(bold(paste("Credible interval coverage: filter by var(", mu, ") > 1")))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.9, 1)) +
  guides(alpha = guide_legend(override.aes = list(fill = "blue"))) +
  theme(legend.position = c(0.05, 0.15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing = unit(.1, "cm"),
        legend.box = "horizontal")

cred_cover_nofilter <- performance_by_type_nofilter %>%
  group_by(q, q_obs, estimator, grouping) %>%
  summarize(mn_cover = mean(cover)) %>%
  ggplot(aes(x = q_obs, y = mn_cover, group = factor(paste(grouping, "_", estimator, sep = "")),
             color = factor(q), shape = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression(Coverage)) +
  labs(color = "q") +
  labs(shape = "Estimator") +
  ggtitle("Credible interval coverage: no filter") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.9, 1)) +
  guides(alpha = guide_legend(override.aes = list(fill = "blue"))) +
  theme(legend.position = c(0.05, 0.15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing = unit(.1, "cm"),
        legend.box = "horizontal")

ggsave(paste0(plots_dir, "vary_q_n_", args$N, "_mn_var_nofilter.png"),
       plot = plot_grid(cred_cover_mn, cred_cover_var, cred_cover_nofilter),
       device = "png", width = 30, height = 25, units = "cm", dpi = 300)

## also, coverage of wald-type intervals for all three
cred_wald_cover <- performance_by_type %>%
  filter(estimator == "no_ve") %>%
  group_by(q, q_obs, estimator, grouping) %>%
  summarize(mn_wald_cover = mean(wald_cover)) %>%
  ggplot(aes(x = q_obs, y = mn_wald_cover, group = grouping,
             color = factor(q))) +
  xlab(expression(q^obs)) +
  ylab(expression(Coverage)) +
  labs(color = "q") +
  ggtitle("Wald-type 'credible' interval coverage") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.9, 1)) +
  theme(legend.position = c(0.6, 0.3),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing = unit(.1, "cm"))
cred_wald_cover_mn <- performance_by_type_mn %>%
  filter(estimator == "no_ve") %>%
  group_by(q, q_obs, estimator, grouping) %>%
  summarize(mn_wald_cover = mean(wald_cover)) %>%
  ggplot(aes(x = q_obs, y = mn_wald_cover, group = grouping,
             color = factor(q))) +
  xlab(expression(q^obs)) +
  ylab(expression(Coverage)) +
  labs(color = "q") +
  ggtitle("Wald-type 'credible' interval coverage: filter by mean(W) > 0.5") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.9, 1)) +
  theme(legend.position = c(0.6, 0.3),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing = unit(.1, "cm"))
cred_wald_cover_var <- performance_by_type_var %>%
  filter(estimator == "no_ve") %>%
  group_by(q, q_obs, estimator, grouping) %>%
  summarize(mn_wald_cover = mean(wald_cover)) %>%
  ggplot(aes(x = q_obs, y = mn_wald_cover, group = grouping,
             color = factor(q))) +
  xlab(expression(q^obs)) +
  ylab(expression(Coverage)) +
  labs(color = "q") +
  ggtitle(expression(bold(paste("Wald-type 'credible interval' coverage: filter by var(", mu, ") > 1")))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.9, 1)) +
  theme(legend.position = c(0.6, 0.3),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing = unit(.1, "cm"))
ggsave(paste0(plots_dir, "vary_q_n_", args$N, "_wald_cover.png"), plot = plot_grid(cred_wald_cover, cred_wald_cover_mn, cred_wald_cover_var),
       device = "png", width = 30, height = 25, units = "cm", dpi = 300)

## -------------------------------------------------------------------------------------------------
## look into contributions to mse
## -------------------------------------------------------------------------------------------------
sq_bias <- performance_by_type %>%
  mutate(sq_bias = bias ^ 2) %>%
  group_by(q, q_obs, estimator, grouping) %>%
  summarize(mn_sq_bias = mean(sq_bias)) %>%
  ggplot(aes(x = q_obs, y = mn_sq_bias, group = factor(paste(grouping, "_", estimator, sep = "")),
             color = factor(q), shape = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression("Squared bias")) +
  ggtitle("Squared estimation bias") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  labs(color = "q") +
  labs(shape = "Estimator") +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 80))
variance <- performance_by_type %>%
  mutate(variance = mse - bias ^ 2) %>%
  group_by(q, q_obs, estimator, grouping) %>%
  summarize(mn_var = mean(variance)) %>%
  ggplot(aes(x = q_obs, y = mn_var, group = factor(paste(grouping, "_", estimator, sep = "")),
             color = factor(q), shape = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression("Estimated variance")) +
  ggtitle("Estimated variance") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  labs(color = "q") +
  labs(shape = "Estimator") +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 80000))

ggsave(paste0(plots_dir, "vary_q_bias_var_n_", args$N, ".png"), plot = plot_grid(sq_bias, variance),
       width = 30, height = 15, units = "cm", dpi = 300)

## -------------------------------------------------------------------------------------------------
## bias vs rank order with no filter; this gets further into MSE
## -------------------------------------------------------------------------------------------------
bias_vs_rank_order_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 60, estimator == "no_ve") %>%
  group_by(taxon_id, estimator) %>%
  summarize(mn_bias = mean(bias)) %>%
  # ggplot(aes(x = taxon_id, y = mn_bias, group = taxon_id, color = factor(taxon_id), shape = factor(estimator, levels = c("naive", "no_ve"), labels = c("Naive", "Proposed, no ve")))) +
  ggplot(aes(x = taxon_id, y = mn_bias, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(-3, 3)) +
  labs(color = "Taxon") +
  labs(y = "Bias") +
  labs(x = "Rank order of W") +
  ggtitle("Bias vs. rank order, q = 60")

bias_vs_rank_order_40_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 40, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_bias = mean(bias)) %>%
  ggplot(aes(x = taxon_id, y = mn_bias, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(-3, 3)) +
  labs(color = "Taxon") +
  labs(y = "Bias") +
  labs(x = "Rank order of W") +
  ggtitle("Bias vs. rank order, q = 40")

bias_vs_rank_order_20_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 20, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_bias = mean(bias)) %>%
  ggplot(aes(x = taxon_id, y = mn_bias, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(-3, 3)) +
  labs(color = "Taxon") +
  labs(y = "Bias") +
  labs(x = "Rank order of W") +
  ggtitle("Bias vs. rank order, q = 20")

bias_vs_rank_order_10_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 10, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_bias = mean(bias)) %>%
  ggplot(aes(x = taxon_id, y = mn_bias, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(-3, 3)) +
  labs(color = "Taxon") +
  labs(y = "Bias") +
  labs(x = "Rank order of W") +
  ggtitle("Bias vs. rank order, q = 10")

ggsave(paste0(plots_dir, "vary_q_bias_vs_rank_order_n_", args$N, "_nofilter.png"),
       plot = plot_grid(bias_vs_rank_order_10_nofilter, bias_vs_rank_order_20_nofilter,
                        bias_vs_rank_order_40_nofilter, bias_vs_rank_order_nofilter),
       width = 30, height = 25, units = "cm", dpi = 300)

## -----------------------------------------------------------------------------------------------------
## plot of coverage vs rank order
## NB: since I'm using the most abundant as qPCR, and the rest are ordered already,
## taxon_id is equivalent to rank order
## -----------------------------------------------------------------------------------------------------

cover_vs_rank_order <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 60, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.6, 1)) +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of credible intervals vs. rank order, q = 60")

pred_cover_vs_rank_order <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 60, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(pred_cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.6, 1)) +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of prediction intervals vs. rank order of W")

width_vs_rank_order <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 60, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_width = mean(width)) %>%
  ggplot(aes(x = taxon_id, y = mn_width, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  labs(color = "Taxon") +
  labs(y = "Width") +
  labs(x = "Rank order of W") +
  ggtitle("Prediction interval width vs. rank order of W")

cover_vs_rank_order_40 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 40, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.6, 1)) +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of credible intervals vs. rank order, q = 40")

pred_cover_vs_rank_order_40 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 40, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(pred_cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  ylim(c(0.6, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of prediction intervals vs. rank order of W")

width_vs_rank_order_40 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 40, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_width = mean(width)) %>%
  ggplot(aes(x = taxon_id, y = mn_width, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  labs(color = "Taxon") +
  labs(y = "Width") +
  labs(x = "Rank order of W") +
  ggtitle("Prediction interval width vs. rank order of W")

cover_vs_rank_order_20 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 20, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  ylim(c(0.6, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of credible intervals vs. rank order, q = 20")

pred_cover_vs_rank_order_20 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 20, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(pred_cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.6, 1)) +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of prediction intervals vs. rank order of W")

width_vs_rank_order_20 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 20, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_width = mean(width)) %>%
  ggplot(aes(x = taxon_id, y = mn_width, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  labs(color = "Taxon") +
  labs(y = "Width") +
  labs(x = "Rank order of W") +
  ggtitle("Prediction interval width vs. rank order of W")

cover_vs_rank_order_10 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 10, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.6, 1)) +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of credible intervals vs. rank order, q = 10")

pred_cover_vs_rank_order_10 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 10, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(pred_cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  ylim(c(0.6, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of prediction intervals vs. rank order of W")

width_vs_rank_order_10 <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 10, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_width = mean(width)) %>%
  ggplot(aes(x = taxon_id, y = mn_width, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  labs(color = "Taxon") +
  labs(y = "Width") +
  labs(x = "Rank order of W") +
  ggtitle("Prediction interval width vs. rank order of W")


ggsave(paste0(plots_dir, "vary_q_cover_vs_rank_order_40_n_", args$N, ".png"),
       plot = cover_vs_rank_order_40,
       width = 15, height = 15, units = "cm", dpi = 300)
ggsave(paste0(plots_dir, "vary_q_pred_cover_vs_rank_order_40_n_", args$N, ".png"), plot = plot_grid(pred_cover_vs_rank_order_40,
                                                                                                    width_vs_rank_order_40),
       width = 30, height = 15, units = "cm", dpi = 300)
ggsave(paste0(plots_dir, "vary_q_cover_vs_rank_order_20_n_", args$N, ".png"), cover_vs_rank_order_20,
       width = 15, height = 15, units = "cm", dpi = 300)
ggsave(paste0(plots_dir, "vary_q_pred_cover_vs_rank_order_20_n_", args$N, ".png"), plot = plot_grid(pred_cover_vs_rank_order_20,
                                                                                                    width_vs_rank_order_20),
       width = 30, height = 15, units = "cm", dpi = 300)
ggsave(paste0(plots_dir, "vary_q_cover_vs_rank_order_10_n_", args$N, ".png"), cover_vs_rank_order_10,
       width = 15, height = 15, units = "cm", dpi = 300)
ggsave(paste0(plots_dir, "vary_q_pred_cover_vs_rank_order_10_n_", args$N, ".png"), plot = plot_grid(pred_cover_vs_rank_order_10,
                                                                                                    width_vs_rank_order_10),
       width = 30, height = 15, units = "cm", dpi = 300)

ggsave(paste0(plots_dir, "vary_q_cover_vs_rank_order_60_n_", args$N, ".png"), cover_vs_rank_order,
       width = 15, height = 15, units = "cm", dpi = 300)
ggsave(paste0(plots_dir, "vary_q_pred_cover_vs_rank_order_60_n_", args$N, ".png"),
       plot = plot_grid(pred_cover_vs_rank_order,
                        width_vs_rank_order),
       width = 30, height = 15, units = "cm", dpi = 300)

## full plot
ggsave(paste0(plots_dir, "vary_q_cover_vs_rank_order_n_", args$N, ".png"),
       plot = plot_grid(cover_vs_rank_order_10, cover_vs_rank_order_20,
                        cover_vs_rank_order_40, cover_vs_rank_order),
       width = 30, height = 25, units = "cm", dpi = 300)

## -------------------------------------------------------------------------------------------------
## coverage vs rank order with no filter
## -------------------------------------------------------------------------------------------------
cover_vs_rank_order_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 60, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.6, 1)) +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of credible intervals vs. rank order, q = 60")

pred_cover_vs_rank_order_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 60, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(pred_cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.6, 1)) +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of prediction intervals vs. rank order of W")

width_vs_rank_order_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 60, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_width = mean(width)) %>%
  ggplot(aes(x = taxon_id, y = mn_width, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  labs(color = "Taxon") +
  labs(y = "Width") +
  labs(x = "Rank order of W") +
  ggtitle("Prediction interval width vs. rank order of W")

cover_vs_rank_order_40_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 40, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.6, 1)) +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of credible intervals vs. rank order, q = 40")

pred_cover_vs_rank_order_40_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 40, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(pred_cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  ylim(c(0.6, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of prediction intervals vs. rank order of W")

width_vs_rank_order_40_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 40, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_width = mean(width)) %>%
  ggplot(aes(x = taxon_id, y = mn_width, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  labs(color = "Taxon") +
  labs(y = "Width") +
  labs(x = "Rank order of W") +
  ggtitle("Prediction interval width vs. rank order of W")

cover_vs_rank_order_20_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 20, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  ylim(c(0.6, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of credible intervals vs. rank order, q = 20")

pred_cover_vs_rank_order_20_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 20, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(pred_cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.6, 1)) +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of prediction intervals vs. rank order of W")

width_vs_rank_order_20_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 20, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_width = mean(width)) %>%
  ggplot(aes(x = taxon_id, y = mn_width, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  labs(color = "Taxon") +
  labs(y = "Width") +
  labs(x = "Rank order of W") +
  ggtitle("Prediction interval width vs. rank order of W")

cover_vs_rank_order_10_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 10, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.6, 1)) +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of credible intervals vs. rank order, q = 10")

pred_cover_vs_rank_order_10_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 10, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(pred_cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  ylim(c(0.6, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(color = "Taxon") +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of prediction intervals vs. rank order of W")

width_vs_rank_order_10_nofilter <- performance_by_type_over_n_nofilter %>%
  filter(q_obs == 7, q == 10, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_width = mean(width)) %>%
  ggplot(aes(x = taxon_id, y = mn_width, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  labs(color = "Taxon") +
  labs(y = "Width") +
  labs(x = "Rank order of W") +
  ggtitle("Prediction interval width vs. rank order of W")


ggsave(paste0(plots_dir, "vary_q_cover_vs_rank_order_40_n_", args$N, "_nofilter.png"),
       plot = cover_vs_rank_order_40_nofilter,
       width = 15, height = 15, units = "cm", dpi = 300)
ggsave(paste0(plots_dir, "vary_q_pred_cover_vs_rank_order_40_n_", args$N, "_nofilter.png"), plot = plot_grid(pred_cover_vs_rank_order_40_nofilter,
                                                                                                             width_vs_rank_order_40_nofilter),
       width = 30, height = 15, units = "cm", dpi = 300)
ggsave(paste0(plots_dir, "vary_q_cover_vs_rank_order_20_n_", args$N, "_nofilter.png"), cover_vs_rank_order_20_nofilter,
       width = 15, height = 15, units = "cm", dpi = 300)
ggsave(paste0(plots_dir, "vary_q_pred_cover_vs_rank_order_20_n_", args$N, "_nofilter.png"), plot = plot_grid(pred_cover_vs_rank_order_20_nofilter,
                                                                                                             width_vs_rank_order_20_nofilter),
       width = 30, height = 15, units = "cm", dpi = 300)
ggsave(paste0(plots_dir, "vary_q_cover_vs_rank_order_10_n_", args$N, "_nofilter.png"), cover_vs_rank_order_10_nofilter,
       width = 15, height = 15, units = "cm", dpi = 300)
ggsave(paste0(plots_dir, "vary_q_pred_cover_vs_rank_order_10_n_", args$N, "_nofilter.png"), plot = plot_grid(pred_cover_vs_rank_order_10_nofilter,
                                                                                                             width_vs_rank_order_10_nofilter),
       width = 30, height = 15, units = "cm", dpi = 300)

ggsave(paste0(plots_dir, "vary_q_cover_vs_rank_order_60_n_", args$N, "_nofilter.png"), cover_vs_rank_order_nofilter,
       width = 16, height = 15, units = "cm", dpi = 300)
ggsave(paste0(plots_dir, "vary_q_pred_cover_vs_rank_order_60_n_", args$N, "_nofilter.png"),
       plot = plot_grid(pred_cover_vs_rank_order_nofilter,
                        width_vs_rank_order_nofilter),
       width = 30, height = 15, units = "cm", dpi = 300)

## full plot
ggsave(paste0(plots_dir, "vary_q_cover_vs_rank_order_n_", args$N, "_nofilter.png"),
       plot = plot_grid(cover_vs_rank_order_10_nofilter, cover_vs_rank_order_20_nofilter,
                        cover_vs_rank_order_40_nofilter, cover_vs_rank_order_nofilter),
       width = 30, height = 25, units = "cm", dpi = 300)

## -------------------------------------------------------------------------------------------------
## debugging:
## -------------------------------------------------------------------------------------------------
## look at bias vs rank order
performance_by_type_over_n %>%
  filter(estimator == "no_ve", q == 40) %>%
  mutate(sq_bias = bias ^ 2) %>%
  ggplot(aes(x = taxon_id, y = sq_bias, group = taxon_id, fill = factor(taxon_id))) +
  xlab(expression(q^obs)) +
  ylab(expression("Squared bias")) +
  ggtitle("Squared estimation bias vs rank order") +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  labs(fill = "Taxon") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 20)) +
  guides(alpha = guide_legend(override.aes = list(fill = "blue")))

## look at number of MC reps for each taxon
num_reps <- performance_by_type_over_n %>%
  group_by(taxon_id, q, q_obs) %>%
  summarize(num_reps = length(unique(mc_id)))

num_reps_mn <- performance_by_type_over_n_mn %>%
  group_by(taxon_id, q, q_obs) %>%
  summarize(num_reps = length(unique(mc_id)))

num_reps_var <- performance_by_type_over_n_var %>%
  group_by(taxon_id, q, q_obs) %>%
  summarize(num_reps = length(unique(mc_id)))

## check for each taxon id
num_reps %>% filter(taxon_id == 46) %>% print(n = Inf)
num_reps_mn %>% filter(taxon_id == 34) %>% print(n = Inf)
num_reps_var %>% filter(taxon_id == 34) %>% print(n = Inf)

## wald coverage vs rank order
wald_cover_vs_rank_order <- performance_by_type_over_n %>%
  filter(q_obs == 7, q == 60, estimator == "no_ve") %>%
  group_by(taxon_id) %>%
  summarize(mn_cover = mean(wald_cover)) %>%
  ggplot(aes(x = taxon_id, y = mn_cover, group = taxon_id, color = factor(taxon_id))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.6, 1)) +
  labs(color = "Taxon") +
  labs(y = "Wald-type 'credible' coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of credible intervals vs. rank order of W")
wald_cover_vs_rank_order
ggsave(paste0(plots_dir, "vary_q_wald_cover_vs_rank_order_n_", args$N, ".png"), wald_cover_vs_rank_order,
       device = "png", width = 15, height = 15, units = "cm", dpi = 300)

## ---------------------------------------------------------------------------------
## convergence diagnostics
## ---------------------------------------------------------------------------------
rhats <- list(length(args$q))
for (j in 1:length(args$q)) {
  results_dir <- paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/q_", args$q[i], "/q_obs_", args$q_obs, "/")
  
  dir_mat <- expand.grid(job = 1:args$num_jobs, dir = results_dir)
  mod_nms_lst <- as.list(paste0(dir_mat$dir, args$stan_model, "_mod_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
  data_nms_lst <- as.list(paste0(dir_mat$dir, args$stan_model, "_data_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
  samp_nms_lst <- as.list(paste0(dir_mat$dir, args$stan_model, "_samps_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
  
  mod_lst <- lapply(mod_nms_lst, read_func)
  data_lst <- lapply(data_nms_lst, read_func)
  beta_summ_lst <- lapply(mod_lst, function(x) summary(x[grepl("beta", rownames(x)), 7]))
  beta_mat <- do.call(rbind.data.frame, beta_summ_lst)
  colnames(beta_mat) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  beta_mat$q_obs <- rep(2:7, each = 50)
  Sigma_summ_lst <- lapply(mod_lst, function(x) summary(x[grepl("Sigma", rownames(x)), 7]))
  Sigma_mat <- do.call(rbind.data.frame, Sigma_summ_lst)
  colnames(Sigma_mat) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  Sigma_mat$q_obs <- rep(2:7, each = 50)
  mu_summ_lst <- lapply(mod_lst, function(x) summary(x[grepl("mu", rownames(x)), 7])[1:6])
  mu_mat <- do.call(rbind.data.frame, mu_summ_lst)
  colnames(mu_mat) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  mu_mat$q_obs <- rep(2:7, each = 50)
  
  beta_means <- beta_mat %>%
    group_by(q_obs) %>%
    summarize(min = mean(Min.), first_q = mean(`1st Qu.`), median = mean(Median), mean = mean(Mean), third_q = mean(`3rd Qu.`), max = mean(Max.))
  Sigma_means <- Sigma_mat %>%
    group_by(q_obs) %>%
    summarize(min = mean(Min.), first_q = mean(`1st Qu.`), median = mean(Median), mean = mean(Mean), third_q = mean(`3rd Qu.`), max = mean(Max.))
  mu_means <- mu_mat %>%
    group_by(q_obs) %>%
    summarize(min = mean(Min.), first_q = mean(`1st Qu.`), median = mean(Median), mean = mean(Mean), third_q = mean(`3rd Qu.`), max = mean(Max.))
  
  rhats[[j]] <- tibble(q = args$q[j], q_obs = rep(args$q_obs, each = 3), parameter = rep(c("beta", "Sigma", "mu"), length(args$q_obs)))
  rhats[[j]] <- tibble::add_column(rhats[[j]], median = c(beta_means$median, Sigma_means$median, mu_means$median))
  rhats[[j]] <- tibble::add_column(rhats[[j]], iqr_low = c(beta_means$first_q, Sigma_means$first_q, mu_means$first_q))
  rhats[[j]] <- tibble::add_column(rhats[[j]], iqr_high = c(beta_means$third_q, Sigma_means$third_q, mu_means$third_q))
}
rhats_all <- do.call(rbind.data.frame, rhats)
write.csv(rhats_all, "plots/vary_q/rhats.csv")

## ------------------------------------------------------------------------------------------------------
## get it without using most abundant, for supplement (this was only the results with no samples)
## ------------------------------------------------------------------------------------------------------
# output_performances_avg_over_taxa_n_ab_0 <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_ab_", 0, ".rds"))
#
# ## transform into a giant matrix for the average over n and taxa
# performance_avg_over_taxa_n_ab_0 <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n_ab_0)
# performance_avg_over_taxa_n_ab_0$grouping <- paste(performance_avg_over_taxa_n_ab_0$q, performance_avg_over_taxa_n_ab_0$q_obs, sep = "_")
#
# performance_by_type_ab_0 <- performance_avg_over_taxa_n_ab_0 %>%
#   mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_pred_cover = pred_cover,
#          no_ve_width = width) %>%
#   select(-mse, -mspe, -cover, -pred_cover, -width) %>%
#   gather(key, value, -q, -q_obs, -mc_id, -grouping) %>%
#   tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>%
#   spread(measure, value) %>%
#   mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))
#
# pred_cover_ab_0 <- performance_by_type_ab_0 %>%
#   ggplot(aes(x = q_obs, y = pred_cover, group = grouping, fill = factor(q))) +
#   xlab(expression(q^obs)) +
#   ylab(expression(Coverage)) +
#   labs(fill = "q") +
#   ggtitle("Prediction interval coverage") +
#   geom_boxplot(size = 0.5, position = position_dodge(0.5), outlier.shape = NA) +
#   scale_fill_manual(values = cols) +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
#   ylim(c(0.7, 1)) +
#   # guides(fill = FALSE) +
#   theme_bw() +
#   theme(legend.position = c(0.9, 0.3))
# ## next, get credible interval coverage
# cred_cover_ab_0 <- performance_by_type_ab_0 %>%
#   ggplot(aes(x = q_obs, y = cover, group = grouping, fill = factor(q))) +
#   xlab(expression(q^obs)) +
#   ylab(expression(Coverage)) +
#   labs(fill = "q") +
#   ggtitle("Credible interval coverage") +
#   geom_boxplot(size = 0.5, position = position_dodge(0.5), outlier.shape = NA) +
#   scale_fill_manual(values = cols) +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
#   ylim(c(0.7, 1)) +
#   guides(fill = FALSE) +
#   theme_bw()
#
# ## get rmspe
# rmspe_ab_0 <- performance_by_type_ab_0 %>%
#   ggplot(aes(x = q_obs, y = rmspe, group = factor(paste(grouping, "_", estimator, sep = "")),
#              fill = factor(q), alpha = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
#   xlab(expression(q^obs)) +
#   ylab(expression(RMSPE)) +
#   ggtitle("Root mean squared prediction error") +
#   geom_point(position = position_dodge(point_dodge), size = point_cex) +
#   labs(alpha = "Estimator") +
#   scale_fill_manual(values = cols) +
#   scale_alpha_discrete(range = c(0.5, 1)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#   # ylim(c(0, 80)) +
#   guides(fill = FALSE, alpha = guide_legend(override.aes = list(fill = "blue"))) +
#   theme_bw() +
#   theme(legend.position = c(0.75, 0.8))
#
# ## get rmspe
# rmse_ab_0 <- performance_by_type_ab_0 %>%
#   ggplot(aes(x = q_obs, y = rmse, group = factor(paste(grouping, "_", estimator, sep = "")),
#              fill = factor(q), alpha = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
#   xlab(expression(q^obs)) +
#   ylab(expression(RMSE)) +
#   ggtitle("Root mean squared error") +
#   geom_point(position = position_dodge(point_dodge), size = point_cex) +
#   labs(alpha = "Estimator") +
#   scale_fill_manual(values = cols) +
#   scale_alpha_discrete(range = c(0.5, 1)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#   # ylim(c(0, 80)) +
#   guides(fill = FALSE, alpha = FALSE) +
#   theme_bw()
#
#
# ## actually plot them
# png(paste0(plots_dir, "vary_q_ab_0.png"), width = fig_width, height = fig_height, units = "px", res = 300)
# gridExtra::grid.arrange(pred_cover_ab_0, cred_cover_ab_0,
#                         rmspe_ab_0, rmse_ab_0,
#                         nrow = 2, ncol = 2)
# dev.off()
#
# ## get widths, for supplement
# width_ab_0 <- performance_by_type_ab_0 %>%
#   filter(estimator == "no_ve") %>%
#   ggplot(aes(x = q_obs, y = width, group = grouping, 
#              fill = factor(q))) +
#   xlab(expression(q^obs)) +
#   ylab(expression(Width)) +
#   ggtitle("Prediction interval width") +
#   geom_boxplot(width = 0.7, position = position_identity(), outlier.shape = NA) +
#   labs(fill = "q") +
#   scale_fill_manual(values = cols) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#   # ylim(c(0, 80)) +
#   theme_bw()
#
# png(paste0(plots_dir, "vary_q_width_ab_0.png"), width = fig_width, height = fig_height, units = "px", res = 300)
# width_ab_0
# dev.off()
