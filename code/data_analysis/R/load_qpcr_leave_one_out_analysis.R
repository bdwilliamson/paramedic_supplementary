######################################################################################
##
## FILE: load_qpcr_leave_one_out_analysis.R
##
## CREATED: 27 January 2019 by Brian Williamson
##
## PURPOSE: data analysis of the qPCR + br16s data from Hutch collaborators
##          to showcase qPCR estimation method; leave-one-out
##
## INPUTS: ../../data/p2_brngs_qPCR_merge_20180314.csv
##         qpcr_data_analysis_est_*.rds 
## 
## OUTPUTS: qPCR/tables/data_analysis/sample_*_most_abundant_*/est_efficiencies.tex
##          qPCR/tables/data_analysis/sample_*_most_abundant_*/point_ests_plus_intervals.tex
######################################################################################
## -----------------------------------------------------------------------------------
## load required functions and libraries
## -----------------------------------------------------------------------------------
code_dir <- "code/R/"
data_dir <- "data/"
stan_dir <- "stan/"
results_dir <- "results/"
library("tidyr")
library("dplyr")
library("ggplot2")
library("rstan")
library("paramedic")
source(paste0(code_dir, "naive_qpcr_estimator.R"))
source(paste0(code_dir, "analyze_data/get_most_abundant_taxa.R"))

## FUNCTION: interval_long_to_wide
## ARGS: intervals - long matrix of intervals
##          digits - number of digits to round to
## RETURNS: intervals as a wide matrix (one column for each taxon, one row for each woman)
interval_long_to_wide <- function(intervals, digits, n_taxa, n_women) {
  chr_intervals <- apply(intervals[, c(2, 3)], 1, function(x) paste0("[", round(x[1], digits), ", ", round(x[2], digits), "]"))
  chr_intervals_df <- data.frame(taxon = 1:n_taxa, interval = chr_intervals, stringsAsFactors = FALSE)
  wide_intervals_mat <- matrix(NA, nrow = n_women, ncol = n_taxa)
  for (i in 1:n_women) {
    wide_intervals_mat[i, ]  <- as.character(tidyr::spread(chr_intervals_df[1:n_taxa + (i-1)*n_taxa, ], taxon, interval))
  }
  return(wide_intervals_mat)
}

## set up arguments
q <- 7
q_obs <- 7
samp_num <- 100
div_num <- 1000
analysis_name <- paste0("data_analysis/sample_", samp_num, "_most_abundant_", q_obs, "/")

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
data_lst <- process_data(full_data, br_inds, qpcr_inds, pcr_plus_br_inds, llod, m_min, div_num)
qpcr <- data_lst$qpcr
br <- data_lst$br

## set up matrices without ids
qpcr_mat <- as.matrix(qpcr[, 2:dim(qpcr)[2]])
br16_mat <- as.matrix(br[, 2:dim(br)[2]])
## calculate the read numbers
m <- rowSums(br16_mat)

set.seed(4747)
samp <- sample(1:dim(br16_mat)[1], samp_num)  

## if q is smaller than the total number, select the most abundant taxa (it always is, for the leave-one-out analysis)
if (q < dim(br16_mat)[2]) {
  ## get the order
  ordered_by_abundance <- get_most_abundant_taxa(br16_mat, m)
  ## remove the ones corresponding to observed qPCR
  taxa_to_estimate <- ordered_by_abundance[!(ordered_by_abundance %in% 1:7)]
  ## select the most abundant taxa (always select the first 7, corresponding to qPCR)
  most_abundant_16S <- br16_mat[, 1:7]
  ## rename br16_mat
  br16_mat <- most_abundant_16S
  ## re-normalize
  m <- rowSums(br16_mat)
}

## --------------------------------------------------------------------------------------------------------
##
## Leave-one-out analysis:
## Compute point estimates, Wald-type PIs, quantile PIs, credible intervals (for range of levels)
## Compute coverage, MSE, MSPE for the left-out one
##
## --------------------------------------------------------------------------------------------------------
leave_one_out_performance <- vector("list", length = 7)
mult_num <- 1 ## report per 1000; setting mult_num <- div_num would report on original scale
for (i in 1:7) {
  ## read in the datasets
  naive_i <- readRDS(paste0(results_dir, analysis_name, "qpcr_data_analysis_est_naive_q_", q, "_q_obs_", q_obs - 1, "_sample_", samp_num, "_loo_", i, ".rds", sep = ""))
  no_ves_i <- readRDS(paste0(results_dir, analysis_name, "qpcr_data_analysis_est_no_ve_q_", q, "_q_obs_", q_obs - 1, "_sample_", samp_num, "_loo_", i, ".rds", sep = ""))
  ves_i <- readRDS(paste0(results_dir, analysis_name, "qpcr_data_analysis_est_ve_q_", q, "_q_obs_", q_obs - 1, "_sample_", samp_num, "_loo_", i, ".rds", sep = ""))
  
  ## -----------------------------------------------------------------------------------
  ## extract the estimates of qPCR and prediction intervals for qPCR
  ## -----------------------------------------------------------------------------------
  ## taxa of interest
  taxa_of_interest <- 7 # always the final one in the estimates
  ## naive estimator
  naive_qpcr_i <- naive_i$mod[, taxa_of_interest]*mult_num
  ## proposed estimator without and with varying efficiency
  no_ve_samps_i <- rstan::extract(no_ves_i$stan_out)
  ve_samps_i <- rstan::extract(ves_i$stan_out)
  no_ve_qpcr_lst_i <- extract_posterior_summaries(no_ves_i$mod, no_ve_samps_i, q = q, taxa_of_interest, mult_num = mult_num, level = 0.95, interval_type = "quantile")
  ve_qpcr_lst_i <- extract_posterior_summaries(ves_i$mod, ve_samps_i, q = q, taxa_of_interest, mult_num = mult_num, level = 0.95, interval_type = "quantile")
  
  ## get mspe
  ests <- cbind(naive_qpcr_i, no_ve_qpcr_lst_i$estimates, ve_qpcr_lst_i$estimates)
  spe_mat <- (ests - qpcr_mat[samp, i]*mult_num)^2
  mspe <- colMeans(spe_mat)
  
  ## get prediction intervals, coverage
  no_ve_pis_quantile_cred <- gen_quantile_interval(mu_quantiles = no_ve_qpcr_lst_i$cred_intervals[,,1], alpha = 0.05, type = "credible_quantiles")
  ve_pis_quantile_cred <- gen_quantile_interval(mu_quantiles = ve_qpcr_lst_i$cred_intervals[,,1], alpha = 0.05, type = "credible_quantiles")
  set.seed(4747)
  no_ve_pis_quantile_samp <- gen_quantile_interval(mu_samps = no_ve_samps_i$mu[, , taxa_of_interest], alpha = 0.05, type = "sample_quantiles")
  set.seed(4747)
  ve_pis_quantile_samp <- gen_quantile_interval(mu_samps = ve_samps_i$mu[, , taxa_of_interest], alpha = 0.05, type = "sample_quantiles")
  no_ve_pis_wald <- gen_wald_interval(no_ve_qpcr_lst_i$estimates, no_ve_qpcr_lst_i$sd, alpha = 0.05)
  ve_pis_wald <- gen_wald_interval(ve_qpcr_lst_i$estimates, ve_qpcr_lst_i$sd, alpha = 0.05)
  
  no_ve_quantile_cred_covers <- no_ve_pis_quantile_cred[,1] <= qpcr_mat[samp, i]*mult_num & no_ve_pis_quantile_cred[,2] >= qpcr_mat[samp, i]*mult_num
  no_ve_quantile_cred_cover <- mean(no_ve_quantile_cred_covers)
  ve_quantile_cred_covers <- ve_pis_quantile_cred[,1] <= qpcr_mat[samp, i]*mult_num & ve_pis_quantile_cred[,2] >= qpcr_mat[samp, i]*mult_num
  ve_quantile_cred_cover <- mean(ve_quantile_cred_covers)
  no_ve_quantile_samp_covers <- no_ve_pis_quantile_samp[,1] <= qpcr_mat[samp, i]*mult_num & no_ve_pis_quantile_samp[,2] >= qpcr_mat[samp, i]*mult_num
  no_ve_quantile_samp_cover <- mean(no_ve_quantile_samp_covers)
  ve_quantile_samp_covers <- ve_pis_quantile_samp[,1] <= qpcr_mat[samp, i]*mult_num & ve_pis_quantile_samp[,2] >= qpcr_mat[samp, i]*mult_num
  ve_quantile_samp_cover <- mean(ve_quantile_samp_covers)
  no_ve_wald_covers <- no_ve_pis_wald[,1] <= qpcr_mat[samp, i]*mult_num & no_ve_pis_wald[,2] >= qpcr_mat[samp, i]*mult_num
  no_ve_wald_cover <- mean(no_ve_wald_covers)
  ve_wald_covers <- ve_pis_wald[,1] <= qpcr_mat[samp, i]*mult_num & ve_pis_wald[,2] >= qpcr_mat[samp, i]*mult_num
  ve_wald_cover <- mean(ve_wald_covers)
  ## get credible intervals
  no_ve_cis_95 <- no_ve_qpcr_lst_i$cred_intervals
  ve_cis_95 <- ve_qpcr_lst_i$cred_intervals
  ## credible interval coverage
  no_ve_cred_cover <- mean(no_ve_cis_95[, 1, 1] <= qpcr_mat[samp, i]*mult_num & no_ve_cis_95[, 2, 1] >= qpcr_mat[samp, i])
  ve_cred_cover <- mean(ve_cis_95[, 1, 1] <= qpcr_mat[samp, i]*mult_num & ve_cis_95[, 2, 1] >= qpcr_mat[samp, i])
  
  ## -----------------------------------------------------------------------------------
  ## extract the estimates of efficiency
  ## -----------------------------------------------------------------------------------
  est_efficiency <- ve_qpcr_lst_i$est_efficiency
  
  ## save them all
  leave_one_out_performance[[i]] <- list(mspe = mspe, spe_mat = spe_mat, no_ve_pis = cbind(no_ve_pis_quantile_cred, no_ve_pis_quantile_samp, no_ve_pis_wald),
                                         ve_pis = cbind(ve_pis_quantile_cred, ve_pis_quantile_samp, ve_pis_wald), 
                                         no_ve_cis = no_ve_cis_95[,,1], ve_cis = ve_cis_95[,,1],
                                         no_ve_quantile_cred_cover = no_ve_quantile_cred_cover, ve_quantile_cred_cover = ve_quantile_cred_cover,
                                         no_ve_quantile_samp_cover = no_ve_quantile_samp_cover, ve_quantile_samp_cover = ve_quantile_samp_cover,
                                         no_ve_wald_cover = no_ve_wald_cover, ve_wald_cover = ve_wald_cover,
                                         obs_qpcr = qpcr_mat[samp, i]*mult_num,
                                         est_efficiency = est_efficiency,
                                         no_ve_cred_cover = no_ve_cred_cover,
                                         ve_cred_cover = ve_cred_cover,
                                         efficiency_samps = ve_samps_i$e,
                                         no_ve_wald_width = abs(no_ve_pis_wald[, 2] - no_ve_pis_wald[, 1]),
                                         ve_wald_width = abs(ve_pis_wald[, 2] - ve_pis_wald[, 1]),
                                         no_ve_quantile_cred_width = abs(no_ve_pis_quantile_cred[, 2] - no_ve_pis_quantile_cred[, 1]),
                                         ve_quantile_cred_width = abs(ve_pis_quantile_cred[, 2] - ve_pis_quantile_cred[, 1]),
                                         no_ve_quantile_samp_width = abs(no_ve_pis_quantile_samp[, 2] - no_ve_pis_quantile_samp[, 1]),
                                         ve_quantile_samp_width = abs(ve_pis_quantile_samp[, 2] - ve_pis_quantile_samp[, 1]),
                                         no_ve_quantile_cred_covers = no_ve_quantile_cred_covers, ve_quantile_cred_covers = ve_quantile_cred_covers,
                                         no_ve_quantile_samp_covers = no_ve_quantile_samp_covers, ve_quantile_samp_covers = ve_quantile_samp_covers,
                                         no_ve_wald_covers = no_ve_wald_covers, ve_wald_covers = ve_wald_covers)
}
## make a tibble with all results for qPCR
names(leave_one_out_performance) <- paste0("LOO", 1:7)
loo_perf_init <- tibble::tibble(naive_mspe = unlist(lapply(leave_one_out_performance, function(x) x$mspe[1])), 
                                no_ve_mspe = unlist(lapply(leave_one_out_performance, function(x) x$mspe[2])), 
                                ve_mspe = unlist(lapply(leave_one_out_performance, function(x) x$mspe[3])),
                                naive_quantile_cred_cover = rep(NA, 7), 
                                no_ve_quantile_cred_cover = unlist(lapply(leave_one_out_performance, function(x) x$no_ve_quantile_cred_cover)),
                                ve_quantile_cred_cover = unlist(lapply(leave_one_out_performance, function(x) x$ve_quantile_cred_cover)),
                                naive_quantile_samp_cover = rep(NA, 7), 
                                no_ve_quantile_samp_cover = unlist(lapply(leave_one_out_performance, function(x) x$no_ve_quantile_samp_cover)),
                                ve_quantile_samp_cover = unlist(lapply(leave_one_out_performance, function(x) x$ve_quantile_samp_cover)),
                                naive_wald_cover = rep(NA, 7), 
                                no_ve_wald_cover = unlist(lapply(leave_one_out_performance, function(x) x$no_ve_wald_cover)), 
                                ve_wald_cover = unlist(lapply(leave_one_out_performance, function(x) x$ve_wald_cover)))
loo_perf_init$loo_index <- 1:7

## get credible interval "coverage"
unlist(lapply(leave_one_out_performance, function(x) x$no_ve_cred_cover))
unlist(lapply(leave_one_out_performance, function(x) x$ve_cred_cover))

loo_perf <- loo_perf_init %>% 
  gather(key, value, -loo_index) %>% 
  tidyr::extract(key, c("estimator", "measure"), regex = "([nv].*_?e)(._?[qmw].*)") %>% 
  spread(measure, value) %>% 
  mutate(mspe = `_mspe`, quantile_cred_cover = `_quantile_cred_cover`, quantile_samp_cover = `_quantile_samp_cover`, wald_cover = `_wald_cover`) %>% 
  select(-`_mspe`, -`_quantile_cred_cover`, `_quantile_samp_cover`, -`_wald_cover`) %>% 
  gather(cover_type, cover, quantile_cred_cover, quantile_samp_cover, wald_cover) %>% 
  mutate(rmspe = sqrt(mspe))

loo_perf$mspe_groups <- paste(loo_perf$estimator, loo_perf$loo_index, sep = "_")
loo_perf$cover_groups <- paste(loo_perf$estimator, loo_perf$cover_type, loo_perf$loo_index, sep = "_")

## read in the estimates on the full data (with 7 observed, 7 total) to compare efficiencies
no_ve_full_results <- readRDS(paste0(results_dir, analysis_name, "qpcr_data_analysis_est_no_ve_q_", q, "_q_obs_", q_obs, "_sample_", samp_num, "_loo_", 999, ".rds", sep = ""))
ve_full_results <- readRDS(paste0(results_dir, analysis_name, "qpcr_data_analysis_est_ve_q_", q, "_q_obs_", q_obs, "_sample_", samp_num, "_loo_", 999, ".rds", sep = ""))
no_ve_full_samps <- rstan::extract(no_ve_full_results$stan_out)
ve_full_samps <- rstan::extract(ve_full_results$stan_out)
no_ve_qpcr_lst_full <- extract_posterior_summaries(no_ve_full_results$mod, no_ve_full_samps, q = q, 1:7, mult_num = mult_num, level = 0.95, interval_type = "quantile")
ve_qpcr_lst_full <- extract_posterior_summaries(ve_full_results$mod, ve_full_samps, q = q, 1:7, mult_num = mult_num, level = 0.95, interval_type = "quantile")
full_est_efficiency <- ve_qpcr_lst_full$est_efficiency
full_efficiency_samps <- ve_full_samps$e

## make a tibble for estimated efficiencies
eff_samps_lst <- lapply(leave_one_out_performance, function(x) x$efficiency_samps)

## make a list of widths
widths_lst <- tibble::tibble(no_ve_wald_width = unlist(lapply(leave_one_out_performance, function(x) x$no_ve_wald_width)),
                             ve_wald_width = unlist(lapply(leave_one_out_performance, function(x) x$ve_wald_width)),
                             no_ve_quantile_cred_width = unlist(lapply(leave_one_out_performance, function(x) x$no_ve_quantile_cred_width)),
                             ve_quantile_cred_width = unlist(lapply(leave_one_out_performance, function(x) x$ve_quantile_cred_width)),
                             no_ve_quantile_samp_width = unlist(lapply(leave_one_out_performance, function(x) x$no_ve_quantile_samp_width)),
                             ve_quantile_samp_width = unlist(lapply(leave_one_out_performance, function(x) x$ve_quantile_samp_width)))
widths_lst$loo_index <- rep(1:7, each = samp_num)
## make a table with width, width_type
widths <- widths_lst %>% 
  gather(key = "width_type", value = "width", -loo_index) %>% 
  tidyr::extract(width_type, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)")


## -------------------------------------------------------------------------------------------------
## PLOTS
## -------------------------------------------------------------------------------------------------
## set up colors
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

## create directory for plots if it doesn't already exist
plots_dir <- paste0("plots/", analysis_name)
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}
## plot coverage against left-out taxon (both quantile and wald-type intervals, no ve and ve)
pred_interval_plot <- ggplot(loo_perf, aes(x = loo_index, y = cover, group = cover_groups, color = factor(estimator, levels = c("naive", "no_ve", "ve"), labels = c("Naive", "Proposed, no ve", "Proposed, ve")), shape = factor(cover_type, levels = c("quantile_cred_cover", "quantile_samp_cover", "wald_cover"), labels = c("Quantile - Cred.", "Quantile - Samp.", "Wald")))) +
  xlab("Left-out taxon index") +
  ylab("Coverage") +
  ggtitle("Prediction interval coverage") +
  labs(color = "Estimator type", shape = "Interval type") +
  geom_point(size = 4, position = position_dodge(width = 1, preserve = "total")) +
  scale_shape_manual(values = pchs) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.35),
        text = element_text(size = 20))
pred_interval_plot

## plot mspe against left-out taxon (naive, no ve, ve)
mspe_plot <- ggplot(loo_perf, aes(x = loo_index, y = rmspe, group = mspe_groups, color = factor(estimator, levels = c("naive", "no_ve", "ve"), labels = c("Naive", "Proposed, no ve", "Proposed, ve")))) +
  xlab("Left-out taxon index") +
  ylab("RMSPE") +
  ggtitle("Root mean squared prediction error") +
  labs(color = "Estimator type") +
  geom_point(size = 4, position = position_dodge(width = 0.5, preserve = "total")) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  guides(color = FALSE) +
  theme_bw() +
  theme(text = element_text(size = 20))
mspe_plot

## plot width against left-out index for each type of interval
width_plot <- widths %>% 
  ggplot(aes(x = loo_index, y = log10(width), group = paste(estimator, measure, loo_index, sep = "_"), 
             fill = factor(estimator, levels = c("naive", "no_ve", "ve"), labels = c("Naive", "Proposed, no ve", "Proposed, ve")),
             color = factor(measure, levels = c("wald_width", "quantile_cred_width", "quantile_samp_width"), labels = c("Wald", "Quantile-cred.", "Quantile-samp.")))) +
  xlab("Left-out taxon index") +
  ylab(expression(paste(Log[10], " width", sep = ""))) +
  ggtitle("Prediction interval width for the left-out taxon") +
  labs(fill = "Estimator") +
  labs(color = "Type of interval") +
  geom_boxplot(width = 0.75, position = position_dodge(), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols[c(3, 5, 7)]) +
  theme_bw() +
  theme(text = element_text(size = 20))
width_plot

## plot estimated efficiencies
for (i in 1:7) {
  ## plot the left-out taxon
  left_out_eff_tibble_i <- tibble(eff = c(eff_samps_lst[[i]][, 7], full_efficiency_samps[, i]), loo = c(rep(1, dim(eff_samps_lst[[i]])[1]), rep(0, dim(full_efficiency_samps)[1])))
  eval(parse(text = paste0("eff_plot_", i, " <- ggplot(left_out_eff_tibble_i, aes(x = eff, group = loo, fill = factor(loo, labels = c('Full data', 'Leave-one-out')))) + 
  geom_density() + 
  xlab('Estimated efficiency') +
  ylab('Density') +
  ggtitle('Density plots of estimated efficiencies (-index ", i,")') +
  labs(fill = 'Data source') +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  xlim(c(0, 15))")))
  eval(parse(text = paste0("eff_plot_", i)))
  
  ## plot all other taxa
  eff_tibble_i <- tibble(eff = c(as.vector(eff_samps_lst[[i]][, -7]), as.vector(full_efficiency_samps[, -i])), 
                         index = rep(rep(c((1:7)[-i]), each = dim(eff_samps_lst[[i]])[1]), 2),
                         loo = c(rep(1, dim(eff_samps_lst[[i]])[1]*6), rep(0, dim(full_efficiency_samps)[1]*6))) %>% 
    mutate(grp = paste0(index, "_", loo))
  eval(parse(text = paste0("eff_all_plot_", i, " <- ggplot(eff_tibble_i, aes(x = eff, group = grp, fill = factor(index), alpha = factor(loo, labels = c('Full data', 'Leave-one-out')))) + 
  geom_density() + 
  xlab('Estimated efficiency') +
  ylab('Density') +
  ggtitle('Density plots of estimated efficiencies (-index ", i,")') +
  labs(fill = 'Taxon index') +
  labs(alpha = 'Data source') +
  scale_fill_manual(values = cols) +
  scale_alpha_discrete(range = c(0.5, 0.1)) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  guides(alpha = guide_legend(override.aes = list(fill = 'blue'))) +
  xlim(c(0, 15))")))
  eval(parse(text = paste0("eff_all_plot_", i)))
  
}

## save all of the plots
png(paste0(plots_dir, "/loo_perf.png"), width = 1.5*fig_width, height = fig_height, units = "px", res = 300)
gridExtra::grid.arrange(pred_interval_plot, mspe_plot, nrow = 1, ncol = 2)
dev.off()

png(paste0(plots_dir, "/loo_widths.png"), width = fig_width, height = fig_height, units = "px", res = 300)
width_plot
dev.off()


for (i in 1:7) {
  png(paste0(plots_dir, "/efficiencies_vs_left_out_index_", i, ".png"), width = fig_width, height = fig_height, units = 'px', res = 300)
  eval(parse(text = paste0("gridExtra::grid.arrange(eff_plot_", i, ", eff_all_plot_", i, ")")))
  dev.off()
}
