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
library("paramedic")
source("load_qpcr_sim_helpers.R") # provides function get_summaries


## grab command-line arguments
parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "vary_q", help = "name of the simulation")
parser$add_argument("--stan-model", default = "predict_qpcr",
                    help = "Which Stan model file to use.")
parser$add_argument("--N", type = "double", default = 50, help = "sample size")
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

# args$most_abundant <- 1

plots_dir <- paste0("plots/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

## --------------------------------------------------------------------------------------------------
## CREATE LIST OF OUTPUT FOR EASY SUMMARIES; GET THOSE SUMMARIES
## --------------------------------------------------------------------------------------------------
if (args$read_data) {
output_performances <- vector(mode = "list", length = length(args$q))
output_performances_single_taxon <- vector(mode = "list", length = length(args$q))

output_performances_avg_over_taxa_n <- vector(mode = "list", length = length(args$q))
output_performances_avg_over_n <- vector(mode = "list", length = length(args$q))

for (i in 1:length(args$q)) {
  ## load in results and data
  results_dir <- paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/q_", args$q[i], "/q_obs_", args$q_obs, "/")
  
  dir_mat <- expand.grid(job = 1:args$num_jobs, dir = results_dir)
  mod_nms_lst <- as.list(paste0(dir_mat$dir, args$stan_model, "_mod_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
  data_nms_lst <- as.list(paste0(dir_mat$dir, args$stan_model, "_data_jobid_", dir_mat$job, "_ad_", args$ad, "0000_mt_", args$mt, "_ab_", args$most_abundant, ".rds"))
  
  mod_lst <- lapply(mod_nms_lst, read_func)
  data_lst <- lapply(data_nms_lst, read_func)
  
  if (args$q[i] == 10 & args$most_abundant) { # need to realign the mus
    for (j in 1:length(data_lst)) {
      data_lst[[j]]$mu <- data_lst[[j]]$mu[, order(data_lst[[1]]$beta, decreasing = TRUE)]
    }
  }
  
  ## (1) pair model summaries with relevant data, for all taxa, all q_obs
  summary_df <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type) get_summaries(w, x, y, z, type), mod_lst, data_lst, NA, MoreArgs = list(z = 1:args$q[i], type = "no_ve"), SIMPLIFY = FALSE))
  ## pair with the monte-carlo id
  summary_df$mc_id <- rep(rep(1:50, each = args$N*args$q[i]), length(args$q_obs))
  
  ## (2) compute performance for each row
  performance_df <- summary_df %>%
    mutate(mse = (mu - est)^2, naive_mse = (mu - naive_est)^2,
           mspe = (qpcr - est)^2, naive_mspe = (qpcr - naive_est)^2,
           cover = cil <= mu & ciu >= mu,
           pred_cover = wald_pred_cil <= qpcr & wald_pred_ciu >= qpcr, # choose wald
           width = abs(wald_pred_ciu - wald_pred_cil),
           bias = mu - est, naive_bias = mu - naive_est)
  
  # performance_df$include_taxa_in_avg <- apply(performance_df, 1, function(x) x[4] %in% tail(1:x[1], x[1] - x[2]))
  performance_df$include_taxa_in_avg <- performance_df$taxon_id > performance_df$q_obs
  
  ## (3) average over MC reps for each taxon
  mc_averaged_performance <- performance_df %>%
    select(q, q_obs, subj_id, taxon_id, mse, naive_mse, mspe, naive_mspe, cover, pred_cover, width, bias, naive_bias) %>%
    group_by(q, q_obs, subj_id, taxon_id) %>%
    summarize(mse = mean(mse), naive_mse = mean(naive_mse), 
              mspe = mean(mspe), naive_mspe = mean(naive_mspe), 
              cover = mean(cover), pred_cover = mean(pred_cover),
              width = mean(width), bias = mean(bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE)) %>%
    ungroup()
  
  ## (4) set flag for taxa of interest
  mc_averaged_performance$include_taxa_in_avg <- apply(mc_averaged_performance, 1, 
                                                       function(x) x[4] %in% tail(1:x[1], x[1] - x[2]))
  ## (5) average over taxa of interest for each q_obs
  performance_across_taxa <- mc_averaged_performance %>%
    group_by(q, q_obs, subj_id) %>%
    filter(include_taxa_in_avg) %>%
    summarize(mse = mean(mse), naive_mse = mean(naive_mse),
              mspe = mean(mspe), naive_mspe = mean(naive_mspe),
              cover = mean(cover), pred_cover = mean(pred_cover),
              width = mean(width), bias = mean(bias), naive_bias = mean(naive_bias))
  
  ## (6) average over n, for each q_obs
  average_over_n <- performance_across_taxa %>%
    group_by(q_obs) %>%
    summarize(mse = mean(mse), naive_mse = mean(naive_mse),
              mspe = mean(mspe), naive_mspe = mean(naive_mspe),
              cover = mean(cover), pred_cover = mean(pred_cover),
              width = mean(width), bias = mean(bias), naive_bias = mean(naive_bias))
  average_over_n_single_taxon <- mc_averaged_performance %>%
    filter(taxon_id == 10) %>%
    group_by(q_obs) %>%
    summarize(mse = mean(mse), naive_mse = mean(naive_mse),
              mspe = mean(mspe), naive_mspe = mean(naive_mspe),
              cover = mean(cover), pred_cover = mean(pred_cover),
              width = mean(width), bias = mean(bias), naive_bias = mean(naive_bias))
  
  ## (5) add rmse, transpose
  performance_matrix <- average_over_n %>%
    mutate(rmse = sqrt(mse), naive_rmse = sqrt(naive_mse),
           rmspe = sqrt(mspe), naive_rmspe = sqrt(naive_mspe))
  performance_matrix_single_taxon <- average_over_n_single_taxon %>%
    mutate(rmse = sqrt(mse), naive_rmse = sqrt(naive_mse),
           rmspe = sqrt(mspe), naive_rmspe = sqrt(naive_mspe))
  
  ## (6) average only over n, taxa
  average_over_n_taxa <- performance_df %>%
    filter(include_taxa_in_avg) %>% 
    select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe, 
           cover, pred_cover, width, bias, naive_bias) %>%
    group_by(q, q_obs, mc_id) %>%
    summarize(mse = mean(mse), naive_mse = mean(naive_mse), 
              mspe = mean(mspe), naive_mspe = mean(naive_mspe), 
              cover = mean(cover), pred_cover = mean(pred_cover),
              width = mean(width), bias = mean(bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE)) %>% 
    ungroup()
  
  ## (7) average only over n
  average_over_n <- performance_df %>% 
    filter(include_taxa_in_avg) %>% 
    select(q, q_obs, subj_id, taxon_id, mc_id, mse, naive_mse, mspe, naive_mspe,
           cover, pred_cover, width, bias, naive_bias) %>% 
    group_by(q, q_obs, taxon_id, mc_id) %>% 
    summarize(mse = mean(mse), naive_mse = mean(naive_mse), 
              mspe = mean(mspe), naive_mspe = mean(naive_mspe), 
              cover = mean(cover), pred_cover = mean(pred_cover),
              width = mean(width), bias = mean(bias, na.rm = TRUE), naive_bias = mean(naive_bias, na.rm = TRUE)) %>% 
    ungroup()
  
  output_performances[[i]] <- performance_matrix
  output_performances_single_taxon[[i]] <- performance_matrix_single_taxon
  output_performances_avg_over_taxa_n[[i]] <- average_over_n_taxa
  output_performances_avg_over_n[[i]] <- average_over_n
}
saveRDS(output_performances, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_ab_", args$most_abundant, ".rds"))
saveRDS(output_performances_single_taxon, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_single_taxon_ab_", args$most_abundant, ".rds"))
saveRDS(output_performances_avg_over_taxa_n, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_ab_", args$most_abundant, ".rds"))
saveRDS(output_performances_avg_over_n, paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_n_ab_", args$most_abundant, ".rds"))
} else {
  output_performances <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_ab_", args$most_abundant, ".rds"))
  output_performances_single_taxon <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_single_taxon_ab_", args$most_abundant, ".rds"))
  output_performances_avg_over_taxa_n <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_ab_", args$most_abundant, ".rds"))
  output_performances_avg_over_n <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_n_ab_", args$most_abundant, ".rds"))
}

## transform to lists with a matrix for each q_obs (so that it can be plotted on x-axis)
q_obs_performance <- q_to_q_obs(output_performances)
q_obs_performance_single_taxon <- q_to_q_obs(output_performances_single_taxon)

## transform into a giant matrix for the average over n and taxa
performance_avg_over_taxa_n <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n)
performance_avg_over_taxa_n$grouping <- paste(performance_avg_over_taxa_n$q, performance_avg_over_taxa_n$q_obs, sep = "_")

performance_by_type <- performance_avg_over_taxa_n %>% 
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_pred_cover = pred_cover,
         no_ve_width = width, no_ve_bias = bias) %>% 
  select(-mse, -mspe, -cover, -pred_cover, -width, -bias) %>% 
  gather(key, value, -q, -q_obs, -mc_id, -grouping) %>% 
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>% 
  spread(measure, value) %>% 
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))

## do the same thing for the average only over n
performance_avg_over_n <- do.call(rbind.data.frame, output_performances_avg_over_n)
performance_avg_over_n$grouping <- paste(performance_avg_over_n$q, performance_avg_over_n$q_obs, performance_avg_over_n$taxon_id, sep = "_")

performance_by_type_over_n <- performance_avg_over_n %>% 
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_pred_cover = pred_cover,
         no_ve_width = width, no_ve_bias = bias) %>% 
  select(-mse, -mspe, -cover, -pred_cover, -width, -bias) %>% 
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
point_cex <- 2
axis_cex <- 1.75

## --------------------------------------------------------------------------------------------------
## FINAL PLOT: 4 panels (prediction coverage, coverage, rmspe, mspe)
## add widths, non-abundance for supplement
## --------------------------------------------------------------------------------------------------
## create the plots
## first, get prediction interval coverage (only the proposed estimator)
pred_cover <- performance_by_type %>%
  ggplot(aes(x = q_obs, y = pred_cover, group = grouping, fill = factor(q))) +
  xlab(expression(q^obs)) +
  ylab(expression(Coverage)) + 
  labs(fill = "q") +
  ggtitle("Prediction interval coverage") +
  geom_boxplot(size = 0.5, position = position_dodge(0.5), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.7, 1)) +
  # guides(fill = FALSE) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.3))
## next, get credible interval coverage
cred_cover <- performance_by_type %>%
  ggplot(aes(x = q_obs, y = cover, group = grouping, fill = factor(q))) +
  xlab(expression(q^obs)) +
  ylab(expression(Coverage)) + 
  labs(fill = "q") +
  ggtitle("Credible interval coverage") +
  geom_boxplot(size = 0.5, position = position_dodge(0.5), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.7, 1)) +
  guides(fill = FALSE) +
  theme_bw()

## get rmspe
rmspe <- performance_by_type %>% 
  ggplot(aes(x = q_obs, y = rmspe, group = factor(paste(grouping, "_", estimator, sep = "")), 
             fill = factor(q), alpha = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression(RMSPE)) + 
  ggtitle("Root mean squared prediction error") +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  labs(alpha = "Estimator") +
  scale_fill_manual(values = cols) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 80)) +
  guides(fill = FALSE, alpha = guide_legend(override.aes = list(fill = "blue"))) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.8))
  
## get rmspe
rmse <- performance_by_type %>% 
  ggplot(aes(x = q_obs, y = rmse, group = factor(paste(grouping, "_", estimator, sep = "")), 
             fill = factor(q), alpha = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression(RMSE)) + 
  ggtitle("Root mean squared error") +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  labs(alpha = "Estimator") +
  scale_fill_manual(values = cols) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 80)) +
  guides(fill = FALSE, alpha = FALSE) +
  theme_bw() 


## actually plot them
png(paste0(plots_dir, "vary_q.png"), width = fig_width, height = fig_height, units = "px", res = 300)
gridExtra::grid.arrange(pred_cover, cred_cover,
                        rmspe, rmse,
                        nrow = 2, ncol = 2)
dev.off()

## get widths, for supplement
width <- performance_by_type %>% 
  filter(estimator == "no_ve") %>% 
  ggplot(aes(x = q_obs, y = width, group = grouping, 
             fill = factor(q))) +
  xlab(expression(q^obs)) +
  ylab(expression(Width)) + 
  ggtitle("Prediction interval width") +
  geom_boxplot(width = 0.7, position = position_identity(), outlier.shape = NA) +
  labs(fill = "q") +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 80)) +
  theme_bw() 

png(paste0(plots_dir, "vary_q_width.png"), width = fig_width, height = fig_height, units = "px", res = 300)
width
dev.off()

## -------------------------------------------------------------------------------------------------
## look into contributions to mse
## -------------------------------------------------------------------------------------------------
sq_bias <- performance_by_type %>% 
  mutate(sq_bias = bias ^ 2) %>% 
  ggplot(aes(x = q_obs, y = sq_bias, group = factor(paste(grouping, "_", estimator, sep = "")), 
             fill = factor(q), alpha = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression("Squared bias")) + 
  ggtitle("Squared estimation bias") +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  labs(fill = "q") +
  labs(alpha = "Estimator") +
  scale_fill_manual(values = cols) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylim(c(0, 5)) +
  guides(alpha = guide_legend(override.aes = list(fill = "blue"))) +
  theme_bw()


png(paste0(plots_dir, "vary_q_bias_2.png"), width = fig_width, height = fig_height, units = "px", res = 300)
sq_bias
dev.off()

## -----------------------------------------------------------------------------------------------------
## plot of coverage vs rank order
## NB: since I'm using the most abundant as qPCR, and the rest are ordered already,
## taxon_id is equivalent to rank order
## -----------------------------------------------------------------------------------------------------

cover_vs_rank_order <- performance_by_type_over_n %>% 
  filter(q_obs == 7, q == 60, estimator == "no_ve") %>%
  ggplot(aes(x = taxon_id, y = cover, group = taxon_id, fill = factor(taxon_id))) +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(fill = "Taxon") +
  # scale_fill_manual(values = cols) +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of credible intervals vs. rank order of W") +
  theme_bw()

pred_cover_vs_rank_order <- performance_by_type_over_n %>% 
  filter(q_obs == 7, q == 60, estimator == "no_ve") %>%
  ggplot(aes(x = taxon_id, y = pred_cover, group = taxon_id, fill = factor(taxon_id))) +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(fill = "Taxon") +
  # scale_fill_manual(values = cols) +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of prediction intervals vs. rank order of W") +
  theme_bw()

width_vs_rank_order <- performance_by_type_over_n %>% 
  filter(q_obs == 7, q == 60, estimator == "no_ve") %>%
  ggplot(aes(x = taxon_id, y = width, group = taxon_id, fill = factor(taxon_id))) +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  labs(fill = "Taxon") +
  # scale_fill_manual(values = cols) +
  labs(y = "Width") +
  labs(x = "Rank order of W") +
  ggtitle("Prediction interval width vs. rank order of W") +
  theme_bw()

pred_cover_vs_rank_order_40 <- performance_by_type_over_n %>% 
  filter(q_obs == 7, q == 40, estimator == "no_ve") %>%
  ggplot(aes(x = taxon_id, y = pred_cover, group = taxon_id, fill = factor(taxon_id))) +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(fill = "Taxon") +
  # scale_fill_manual(values = cols) +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of prediction intervals vs. rank order of W") +
  theme_bw()

width_vs_rank_order_40 <- performance_by_type_over_n %>% 
  filter(q_obs == 7, q == 40, estimator == "no_ve") %>%
  ggplot(aes(x = taxon_id, y = width, group = taxon_id, fill = factor(taxon_id))) +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  labs(fill = "Taxon") +
  # scale_fill_manual(values = cols) +
  labs(y = "Width") +
  labs(x = "Rank order of W") +
  ggtitle("Prediction interval width vs. rank order of W") +
  theme_bw()

pred_cover_vs_rank_order_20 <- performance_by_type_over_n %>% 
  filter(q_obs == 7, q == 20, estimator == "no_ve") %>%
  ggplot(aes(x = taxon_id, y = pred_cover, group = taxon_id, fill = factor(taxon_id))) +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(fill = "Taxon") +
  # scale_fill_manual(values = cols) +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of prediction intervals vs. rank order of W") +
  theme_bw()

width_vs_rank_order_20 <- performance_by_type_over_n %>% 
  filter(q_obs == 7, q == 20, estimator == "no_ve") %>%
  ggplot(aes(x = taxon_id, y = width, group = taxon_id, fill = factor(taxon_id))) +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  labs(fill = "Taxon") +
  # scale_fill_manual(values = cols) +
  labs(y = "Width") +
  labs(x = "Rank order of W") +
  ggtitle("Prediction interval width vs. rank order of W") +
  theme_bw()

pred_cover_vs_rank_order_10 <- performance_by_type_over_n %>% 
  filter(q_obs == 7, q == 10, estimator == "no_ve") %>%
  ggplot(aes(x = taxon_id, y = pred_cover, group = taxon_id, fill = factor(taxon_id))) +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(fill = "Taxon") +
  # scale_fill_manual(values = cols) +
  labs(y = "Coverage") +
  labs(x = "Rank order of W") +
  ggtitle("Coverage of prediction intervals vs. rank order of W") +
  theme_bw()

width_vs_rank_order_10 <- performance_by_type_over_n %>% 
  filter(q_obs == 7, q == 10, estimator == "no_ve") %>%
  ggplot(aes(x = taxon_id, y = width, group = taxon_id, fill = factor(taxon_id))) +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  labs(fill = "Taxon") +
  # scale_fill_manual(values = cols) +
  labs(y = "Width") +
  labs(x = "Rank order of W") +
  ggtitle("Prediction interval width vs. rank order of W") +
  theme_bw()

  
png(paste0(plots_dir, "vary_q_cover_vs_rank_order.png"), width = fig_width, height = fig_height, units = "px", res = 300)
cover_vs_rank_order
dev.off()

png(paste0(plots_dir, "vary_q_pred_cover_vs_rank_order_60.png"), width = fig_width, height = fig_height, units = "px", res = 300)
gridExtra::grid.arrange(pred_cover_vs_rank_order, 
                        width_vs_rank_order)
dev.off()
png(paste0(plots_dir, "vary_q_pred_cover_vs_rank_order_40.png"), width = fig_width, height = fig_height, units = "px", res = 300)
gridExtra::grid.arrange(pred_cover_vs_rank_order_40, 
                        width_vs_rank_order_40)
dev.off()
png(paste0(plots_dir, "vary_q_pred_cover_vs_rank_order_20.png"), width = fig_width, height = fig_height, units = "px", res = 300)
gridExtra::grid.arrange(pred_cover_vs_rank_order_20, 
                        width_vs_rank_order_20)
dev.off()
png(paste0(plots_dir, "vary_q_pred_cover_vs_rank_order_10.png"), width = fig_width, height = fig_height, units = "px", res = 300)
gridExtra::grid.arrange(pred_cover_vs_rank_order_10, 
                        width_vs_rank_order_10)
dev.off()
## ------------------------------------------------------------------------------------------------------
## get it without using most abundant, for supplement
## ------------------------------------------------------------------------------------------------------
output_performances_avg_over_taxa_n_ab_0 <- readRDS(paste0("results/", args$sim_name, "/ve_", args$sigma, "/cov_", args$corr, "/n_", args$N, "/output_performances_avg_over_taxa_n_ab_", 0, ".rds"))

## transform into a giant matrix for the average over n and taxa
performance_avg_over_taxa_n_ab_0 <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n_ab_0)
performance_avg_over_taxa_n_ab_0$grouping <- paste(performance_avg_over_taxa_n_ab_0$q, performance_avg_over_taxa_n_ab_0$q_obs, sep = "_")

performance_by_type_ab_0 <- performance_avg_over_taxa_n_ab_0 %>% 
  mutate(no_ve_mse = mse, no_ve_mspe = mspe, no_ve_cover = cover, no_ve_pred_cover = pred_cover,
         no_ve_width = width) %>% 
  select(-mse, -mspe, -cover, -pred_cover, -width) %>% 
  gather(key, value, -q, -q_obs, -mc_id, -grouping) %>% 
  tidyr::extract(key, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)") %>% 
  spread(measure, value) %>% 
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))

pred_cover_ab_0 <- performance_by_type_ab_0 %>%
  ggplot(aes(x = q_obs, y = pred_cover, group = grouping, fill = factor(q))) +
  xlab(expression(q^obs)) +
  ylab(expression(Coverage)) + 
  labs(fill = "q") +
  ggtitle("Prediction interval coverage") +
  geom_boxplot(size = 0.5, position = position_dodge(0.5), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.7, 1)) +
  # guides(fill = FALSE) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.3))
## next, get credible interval coverage
cred_cover_ab_0 <- performance_by_type_ab_0 %>%
  ggplot(aes(x = q_obs, y = cover, group = grouping, fill = factor(q))) +
  xlab(expression(q^obs)) +
  ylab(expression(Coverage)) + 
  labs(fill = "q") +
  ggtitle("Credible interval coverage") +
  geom_boxplot(size = 0.5, position = position_dodge(0.5), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0.7, 1)) +
  guides(fill = FALSE) +
  theme_bw()

## get rmspe
rmspe_ab_0 <- performance_by_type_ab_0 %>% 
  ggplot(aes(x = q_obs, y = rmspe, group = factor(paste(grouping, "_", estimator, sep = "")), 
             fill = factor(q), alpha = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression(RMSPE)) + 
  ggtitle("Root mean squared prediction error") +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  labs(alpha = "Estimator") +
  scale_fill_manual(values = cols) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  # ylim(c(0, 80)) +
  guides(fill = FALSE, alpha = guide_legend(override.aes = list(fill = "blue"))) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.8))

## get rmspe
rmse_ab_0 <- performance_by_type_ab_0 %>% 
  ggplot(aes(x = q_obs, y = rmse, group = factor(paste(grouping, "_", estimator, sep = "")), 
             fill = factor(q), alpha = factor(estimator, levels = c("naive", "no_ve"), ordered = FALSE, labels = c("Naive", "Proposed, no ve")))) +
  xlab(expression(q^obs)) +
  ylab(expression(RMSE)) + 
  ggtitle("Root mean squared error") +
  geom_boxplot(width = 1, position = position_dodge(width = 0.7, preserve = "total"), outlier.shape = NA) +
  labs(alpha = "Estimator") +
  scale_fill_manual(values = cols) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  # ylim(c(0, 80)) +
  guides(fill = FALSE, alpha = FALSE) +
  theme_bw() 


## actually plot them
png(paste0(plots_dir, "vary_q_ab_0.png"), width = fig_width, height = fig_height, units = "px", res = 300)
gridExtra::grid.arrange(pred_cover_ab_0, cred_cover_ab_0,
                        rmspe_ab_0, rmse_ab_0,
                        nrow = 2, ncol = 2)
dev.off()

## get widths, for supplement
width_ab_0 <- performance_by_type_ab_0 %>% 
  filter(estimator == "no_ve") %>% 
  ggplot(aes(x = q_obs, y = width, group = grouping, 
             fill = factor(q))) +
  xlab(expression(q^obs)) +
  ylab(expression(Width)) + 
  ggtitle("Prediction interval width") +
  geom_boxplot(width = 0.7, position = position_identity(), outlier.shape = NA) +
  labs(fill = "q") +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  # ylim(c(0, 80)) +
  theme_bw() 

png(paste0(plots_dir, "vary_q_width_ab_0.png"), width = fig_width, height = fig_height, units = "px", res = 300)
width_ab_0
dev.off()