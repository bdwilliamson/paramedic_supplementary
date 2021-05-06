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
library("rstan")
library("tidyverse")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("paramedic")
source(paste0(code_dir, "analyze_data/qpcr_analysis_helper_functions.R"))
source(paste0(code_dir, "naive_qpcr_estimator.R"))
source(paste0(code_dir, "analyze_data/get_most_abundant_taxa.R"))
source(paste0(code_dir, "gen_naive_interval.R"))

## set up arguments
q <- 13
q_obs <- 13
div_num <- 1000
analysis_name <- paste0("data_analysis/leave_one_out/")
# valid values here are 1 (adjusted, main ms fig 6) and 0 (unadjusted, fig s14)
adjust <- 0

## -----------------------------------------------------------------------------------
## load in the data, clean it up
## -----------------------------------------------------------------------------------
analysis_data <- readRDS(paste0(data_dir, "analysis_data.rds"))
qpcr <- analysis_data$qpcr %>%
  mutate_all(~ . / div_num)
br <- analysis_data$br16S
case_control <- (analysis_data$case_control %>%
                   mutate(case_bin = as.numeric(case == "Case")) %>%
                   select(-case))$case_bin
data_lst <- list(qpcr = qpcr, br = br, cc = case_control)

qpcr_mat <- as.matrix(qpcr)
br16_mat <- as.matrix(br)
br16_mat_nms <- colnames(br16_mat)
## calculate the read numbers
m <- rowSums(br16_mat)
set.seed(4747)
samp_num <- 55
## sample women
samp <- if (samp_num == nrow(br16_mat)) {1:nrow(br16_mat)} else {sample(1:nrow(br16_mat), samp_num)}

## --------------------------------------------------------------------------------------------------------
##
## Leave-one-out analysis:
## Compute point estimates, Wald-type PIs, quantile PIs, credible intervals (for range of levels)
## Compute coverage, MSE, MSPE for the left-out one
##
## --------------------------------------------------------------------------------------------------------
leave_one_out_performance <- vector("list", length = q)
mult_num <- 1 ## report per 1000; setting mult_num <- div_num would report on original scale
for (i in 1:length(leave_one_out_performance)) {
  ## read in the datasets
  naive_i <- readRDS(paste0(results_dir, analysis_name, "qpcr_data_analysis_est_naive_q_", q, "_q_obs_", q_obs - 1, "_sample_", samp_num, "_loo_", i, "_adjust_", adjust, ".rds", sep = ""))
  no_ves_i <- readRDS(paste0(results_dir, analysis_name, "qpcr_data_analysis_est_no_ve_q_", q, "_q_obs_", q_obs - 1, "_sample_", samp_num, "_loo_", i, "_adjust_", adjust, ".rds", sep = ""))
  ves_i <- readRDS(paste0(results_dir, analysis_name, "qpcr_data_analysis_est_ve_q_", q, "_q_obs_", q_obs - 1, "_sample_", samp_num, "_loo_", i, "_adjust_", adjust, ".rds", sep = ""))
  
  ## extract the data inputted to stan for estimators
  taxa_of_interest <- q # always the final one in the estimates
  W <- no_ves_i$stan_data_lst$W
  V <- no_ves_i$stan_data_lst$V
  m <- rowSums(W)
  samp_i <- 1:nrow(W)
  full_v <- qpcr_mat[samp, ]
  ## -----------------------------------------------------------------------------------
  ## extract the estimates of qPCR and prediction intervals for qPCR
  ## -----------------------------------------------------------------------------------
  ## naive estimator
  naive_qpcr_i <- naive_estimator(idx = taxa_of_interest, brs = W, qpcrs = V, known_qpcr = 1:q_obs)
  ## proposed estimator without and with varying efficiency
  no_ve_samps_i <- tryCatch(rstan::extract(no_ves_i$stan_out), error = function(e) NA) # purposeful NA if stan model wasn't saved
  ve_samps_i <- tryCatch(rstan::extract(ves_i$stan_out), error = function(e) NA)
  no_ve_qpcr_lst_i <- paramedic::extract_posterior_summaries(no_ves_i$mod, no_ve_samps_i, taxa_of_interest, mult_num = mult_num, level = 0.95, interval_type = "wald")
  ve_qpcr_lst_i <- paramedic::extract_posterior_summaries(ves_i$mod, ve_samps_i, taxa_of_interest, mult_num = mult_num, level = 0.95, interval_type = "wald")
  
  ## get mspe
  ests <- cbind(naive_qpcr_i, no_ve_qpcr_lst_i$estimates, ve_qpcr_lst_i$estimates)
  spe_mat <- (ests - full_v[, i]*mult_num)^2
  mspe <- colMeans(spe_mat)
  
  ## get prediction intervals, coverage
  naive_intervals <- generate_naive_intervals_single_subject(v = V, w = W, m = m, indx = 7, q_obs = 6, level = 0.05)
  no_ve_pis_wald <- paramedic::gen_wald_interval(no_ve_qpcr_lst_i$estimates, no_ve_qpcr_lst_i$sd, alpha = 0.05)
  ve_pis_wald <- paramedic::gen_wald_interval(ve_qpcr_lst_i$estimates, ve_qpcr_lst_i$sd, alpha = 0.05)
  
  naive_wald_covers <- naive_intervals$pi[, 1] <= full_v[, i]*mult_num & naive_intervals$pi[, 2] >= full_v[, i]*mult_num
  naive_wald_cover <- mean(naive_wald_covers, na.rm = TRUE)
  no_ve_wald_covers <- no_ve_pis_wald[,1] <= full_v[, i]*mult_num & no_ve_pis_wald[,2] >= full_v[, i]*mult_num
  no_ve_wald_cover <- mean(no_ve_wald_covers)
  ve_wald_covers <- ve_pis_wald[,1] <= full_v[, i]*mult_num & ve_pis_wald[,2] >= full_v[, i]*mult_num
  ve_wald_cover <- mean(ve_wald_covers)
  
  ## -----------------------------------------------------------------------------------
  ## extract the estimates of efficiency
  ## -----------------------------------------------------------------------------------
  est_efficiency <- ve_qpcr_lst_i$est_efficiency
  
  ## save them all
  leave_one_out_performance[[i]] <- list(mspe = mspe, spe_mat = spe_mat, no_ve_pis = no_ve_pis_wald,
                                         ve_pis = ve_pis_wald, 
                                         naive_wald_cover = naive_wald_cover, no_ve_wald_cover = no_ve_wald_cover, ve_wald_cover = ve_wald_cover,
                                         obs_qpcr = full_v[, i]*mult_num,
                                         est_efficiency = est_efficiency,
                                         efficiency_samps = tryCatch(ve_samps_i$e, error = function(e) NA),
                                         naive_wald_width = abs(naive_intervals$pi[, 2] - naive_intervals$pi[, 1]),
                                         no_ve_wald_width = abs(no_ve_pis_wald[, 2] - no_ve_pis_wald[, 1]),
                                         ve_wald_width = abs(ve_pis_wald[, 2] - ve_pis_wald[, 1]),
                                         no_ve_wald_covers = no_ve_wald_covers, ve_wald_covers = ve_wald_covers)
}
## make a tibble with all results for qPCR
names(leave_one_out_performance) <- paste0("LOO_", colnames(qpcr_mat))
loo_perf_init <- tibble::tibble(naive_mspe = unlist(lapply(leave_one_out_performance, function(x) x$mspe[1])), 
                                no_ve_mspe = unlist(lapply(leave_one_out_performance, function(x) x$mspe[2])), 
                                ve_mspe = unlist(lapply(leave_one_out_performance, function(x) x$mspe[3])),
                                naive_wald_cover = unlist(lapply(leave_one_out_performance, function(x) x$naive_wald_cover)), 
                                no_ve_wald_cover = unlist(lapply(leave_one_out_performance, function(x) x$no_ve_wald_cover)), 
                                ve_wald_cover = unlist(lapply(leave_one_out_performance, function(x) x$ve_wald_cover)))
loo_perf_init$loo_index <- 1:q

loo_perf <- loo_perf_init %>% 
  gather(key, value, -loo_index) %>% 
  tidyr::extract(key, c("estimator", "measure"), regex = "([nv].*_?e)(._?[qmw].*)") %>% 
  spread(measure, value) %>% 
  mutate(mspe = `_mspe`, wald_cover = `_wald_cover`) %>% 
  select(-`_mspe`, -`_wald_cover`) %>% 
  gather(cover_type, cover, wald_cover) %>% 
  mutate(rmspe = sqrt(mspe))

loo_perf$mspe_groups <- paste(loo_perf$estimator, loo_perf$loo_index, sep = "_")
loo_perf$cover_groups <- paste(loo_perf$estimator, loo_perf$cover_type, loo_perf$loo_index, sep = "_")

## read in the estimates on the full data (with 7 observed, 7 total) to compare efficiencies
no_ve_full_results <- readRDS(paste0(results_dir, analysis_name, "qpcr_data_analysis_est_no_ve_q_", q, "_q_obs_", q_obs, "_sample_", samp_num, "_loo_", 999, "_adjust_", adjust, ".rds", sep = ""))
ve_full_results <- readRDS(paste0(results_dir, analysis_name, "qpcr_data_analysis_est_ve_q_", q, "_q_obs_", q_obs, "_sample_", samp_num, "_loo_", 999, "_adjust_", adjust, ".rds", sep = ""))
# no_ve_full_samps <- rstan::extract(no_ve_full_results$stan_out)
# ve_full_samps <- rstan::extract(ve_full_results$stan_out)
no_ve_full_samps <- no_ve_full_results$samps
ve_full_samps <- ve_full_results$samps
no_ve_qpcr_lst_full <- paramedic::extract_posterior_summaries(no_ve_full_results$mod, no_ve_full_samps, 1:13, mult_num = mult_num, level = 0.95, interval_type = "quantile")
ve_qpcr_lst_full <- paramedic::extract_posterior_summaries(ve_full_results$mod, ve_full_samps, 1:13, mult_num = mult_num, level = 0.95, interval_type = "quantile")

## convergence diagnostics:
## rhats for mu 
summary(no_ve_full_results$mod[grepl("mu", rownames(no_ve_full_results$mod)), 7])
summary(no_ve_full_results$mod[grepl("beta", rownames(no_ve_full_results$mod)), 7])
summary(no_ve_full_results$mod[grepl("Sigma", rownames(no_ve_full_results$mod)), 7])

summary(ve_full_results$mod[grepl("mu", rownames(ve_full_results$mod)), 7])
summary(ve_full_results$mod[grepl("beta", rownames(ve_full_results$mod)), 7])
summary(ve_full_results$mod[grepl("Sigma", rownames(ve_full_results$mod)), 7])
summary(ve_full_results$mod[grepl("e", rownames(ve_full_results$mod)) & !grepl("beta", rownames(ve_full_results$mod)), 7])

## make a list of widths
widths_lst <- tibble::tibble(naive_wald_width = unlist(lapply(leave_one_out_performance, function(x) x$naive_wald_width)),
                             no_ve_wald_width = unlist(lapply(leave_one_out_performance, function(x) x$no_ve_wald_width)),
                             ve_wald_width = unlist(lapply(leave_one_out_performance, function(x) x$ve_wald_width)))
widths_lst$loo_index <- rep(1:q, each = samp_num)
## make a table with width, width_type
widths <- widths_lst %>% 
  gather(key = "width_type", value = "width", -loo_index) %>% 
  tidyr::extract(width_type, c("estimator", "measure"), regex = "(n*.*ve)._?([qmwcbp].*)")


## -------------------------------------------------------------------------------------------------
## PLOTS
## -------------------------------------------------------------------------------------------------
loo_perf_plot_tib <- loo_perf %>% 
  mutate(shape_fct = factor(estimator,
                            levels = c("naive", "no_ve", "ve"),
                            labels = c("Naive", 
                                       "Efficiency-naive Bayes",
                                       "Varying-efficiency Bayes")))
width_plot_tib <- widths %>% 
  mutate(shape_fct = factor(estimator,
                            levels = c("naive", "no_ve", "ve"),
                            labels = c("Naive", 
                                       "Efficiency-naive Bayes",
                                       "Varying-efficiency Bayes")))
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
fig_width <- 18
fig_height <- fig_width/1.5
point_cex <- 2
axis_cex <- 1.75
text_size <- 20

## create directory for plots if it doesn't already exist
plots_dir <- paste0("plots/", analysis_name)
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}
## plot coverage against left-out taxon (both quantile and wald-type intervals, no ve and ve)
taxon_nms <- c("A. christensenii", "A. vaginae",
               "BVAB2 spp.", "D. micraerophilus",
               "E. spp. type 1",
               "G. vaginalis", "L. crispatus",
               "L. iners", "L. jensenii",
               "M. hominis", "P. spp. type 1",
               "P. bennonis", "P. micra")
pred_interval_plot <- loo_perf_plot_tib %>% 
  filter(cover_type == "wald_cover") %>% 
  ggplot(aes(x = factor(loo_index, labels = taxon_nms), 
             y = cover, group = cover_groups, 
             shape = shape_fct)) +
  xlab("Left-out taxon") +
  ylab("Coverage") +
  ggtitle("Prediction interval coverage") +
  labs(shape = "Estimator type") +
  geom_point(size = 4, position = position_dodge(width = 0.3, preserve = "total")) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0, 1)) +
  theme(legend.position = c(0.6, 0.25),
        legend.box.background = element_rect(color = "black"),
        text = element_text(size = text_size),
        plot.title = element_text(size = text_size*4/3),
        axis.text = element_text(size = text_size*3/2),
        axis.title = element_text(size = text_size*3/2),
        axis.text.x = element_text(face = "italic", angle = 280, hjust = 0),
        legend.text = element_text(size = text_size),
        legend.title = element_text(size = text_size),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
pred_interval_plot

## plot mspe against left-out taxon (naive, no ve, ve)
mspe_plot <- loo_perf_plot_tib %>% 
  ggplot(aes(x = factor(loo_index, labels = taxon_nms), y = rmspe, 
             group = mspe_groups, 
             shape = shape_fct)) +
  xlab("Left-out taxon") +
  ylab("RMSPE") +
  scale_y_log10() +
  ggtitle("Root mean squared prediction error") +
  labs(color = "Estimator type") +
  geom_point(size = 4, position = position_dodge(width = 0.5, preserve = "total")) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  guides(shape = FALSE) +
  theme(text = element_text(size = text_size),
        plot.title = element_text(size = text_size*4/3),
        axis.text = element_text(size = text_size*3/2),
        axis.title = element_text(size = text_size*3/2),
        axis.text.x = element_text(face = "italic", angle = 280, hjust = 0),
        legend.text = element_text(size = text_size),
        plot.margin = unit(c(0, 1.75, 0, 0.75), "cm"))
mspe_plot

## save the plot
ggsave(paste0(plots_dir, "/loo_perf_adjust_", adjust, ".png"),
       plot = plot_grid(pred_interval_plot, mspe_plot),
       width = fig_width, height = 10)

## plot width against left-out index for each type of interval
width_plot <- width_plot_tib %>% 
  ggplot(aes(x = loo_index, y = log10(width), 
             group = paste(estimator, measure, loo_index, sep = "_"), 
             fill = shape_fct)) +
  xlab("Left-out taxon index") +
  ylab(expression(paste(Log[10], " width", sep = ""))) +
  ggtitle("Prediction interval width for the left-out taxon") +
  labs(fill = "Estimator") +
  geom_boxplot(width = 0.75, position = position_dodge(), outlier.shape = NA) +
  scale_fill_manual(values = cols) +
  guides(color = FALSE) +
  theme(text = element_text(size = text_size),
        plot.title = element_text(size = text_size),
        axis.text = element_text(size = text_size*3/4))
width_plot

## save the plot
ggsave(paste0(plots_dir, "/loo_widths_adjust_", adjust, ".png"),
       plot = width_plot,
       width = fig_width, height = fig_height)

## plot estimated efficiencies
## violin plot with two panels
main_font_size <- 30
main_font_size_lab <- 9
title_font_size <- 26
str_width <- 20

## make a tibble for estimated efficiencies
eff_samps_lst <- lapply(leave_one_out_performance, function(x) x$efficiency_samps)
full_est_efficiency <- ve_qpcr_lst_full$est_efficiency
full_efficiency_samps <- ve_full_samps$e
colnames(full_efficiency_samps) <- colnames(qpcr_mat)
sd(full_est_efficiency)

## (1) create tibbles with samples and the correct taxon index, latin name
## want the full samples for each taxon, both when leaving out the named one and not
g_vaginalis_tibble <- tibble(est_eff = c(full_efficiency_samps[, "aero"], 
                                         full_efficiency_samps[, "atop"],
                                         full_efficiency_samps[, "bvab2"], 
                                         full_efficiency_samps[, "dial1"],
                                         full_efficiency_samps[, "egg"], 
                                         full_efficiency_samps[, "gvag"],
                                         full_efficiency_samps[, "lcrisp"],
                                         full_efficiency_samps[, "liners"],
                                         full_efficiency_samps[, "ljens"],
                                         full_efficiency_samps[, "myco"],
                                         full_efficiency_samps[, "porphtype1"],
                                         full_efficiency_samps[, "porphbenn"],
                                         full_efficiency_samps[, "parv"],
                                         eff_samps_lst$LOO_gvag[, 1], 
                                         eff_samps_lst$LOO_gvag[, 2],
                                         eff_samps_lst$LOO_gvag[, 3], 
                                         eff_samps_lst$LOO_gvag[, 4],
                                         eff_samps_lst$LOO_gvag[, 5], 
                                         eff_samps_lst$LOO_gvag[, 13],
                                         eff_samps_lst$LOO_gvag[, 6],
                                         eff_samps_lst$LOO_gvag[, 7],
                                         eff_samps_lst$LOO_gvag[, 8],
                                         eff_samps_lst$LOO_gvag[, 9],
                                         eff_samps_lst$LOO_gvag[, 10],
                                         eff_samps_lst$LOO_gvag[, 11],
                                         eff_samps_lst$LOO_gvag[, 12]),
                             taxon_name = c(rep(taxon_nms, each = nrow(full_efficiency_samps)), 
                                            rep(taxon_nms, each = nrow(eff_samps_lst$LOO_gvag))),
                             type = c(rep("full", dim(full_efficiency_samps)[1]*q_obs), 
                                      rep("loo", dim(eff_samps_lst$LOO_gvag)[1]*q_obs)),
                             iter = c(rep(1:dim(full_efficiency_samps)[1], q_obs),
                                      rep(1:nrow(eff_samps_lst$LOO_gvag), q_obs)))
bvab2_tibble <- tibble(est_eff = c(full_efficiency_samps[, "aero"], 
                                   full_efficiency_samps[, "atop"],
                                   full_efficiency_samps[, "bvab2"], 
                                   full_efficiency_samps[, "dial1"],
                                   full_efficiency_samps[, "egg"], 
                                   full_efficiency_samps[, "gvag"],
                                   full_efficiency_samps[, "lcrisp"],
                                   full_efficiency_samps[, "liners"],
                                   full_efficiency_samps[, "ljens"],
                                   full_efficiency_samps[, "myco"],
                                   full_efficiency_samps[, "porphtype1"],
                                   full_efficiency_samps[, "porphbenn"],
                                   full_efficiency_samps[, "parv"],
                                   eff_samps_lst$LOO_bvab2[, 1], 
                                   eff_samps_lst$LOO_bvab2[, 2],
                                   eff_samps_lst$LOO_bvab2[, 13], 
                                   eff_samps_lst$LOO_bvab2[, 3],
                                   eff_samps_lst$LOO_bvab2[, 4], 
                                   eff_samps_lst$LOO_bvab2[, 5],
                                   eff_samps_lst$LOO_bvab2[, 6],
                                   eff_samps_lst$LOO_bvab2[, 7],
                                   eff_samps_lst$LOO_bvab2[, 8],
                                   eff_samps_lst$LOO_bvab2[, 9],
                                   eff_samps_lst$LOO_bvab2[, 10],
                                   eff_samps_lst$LOO_bvab2[, 11],
                                   eff_samps_lst$LOO_bvab2[, 12]),
                       taxon_name = c(rep(taxon_nms, each = nrow(full_efficiency_samps)), 
                                      rep(taxon_nms, each = nrow(eff_samps_lst$LOO_bvab2))),
                       type = c(rep("full", dim(full_efficiency_samps)[1]*q_obs), 
                                rep("loo", dim(eff_samps_lst$LOO_bvab2)[1]*q_obs)),
                       iter = c(rep(1:dim(full_efficiency_samps)[1], q_obs),
                                rep(1:nrow(eff_samps_lst$LOO_bvab2), q_obs)))
get_adjacent_indices <- function(index_vec, q_obs) {
  full_vec <- index_vec[1:q_obs]
  loo_vec <- index_vec[(q_obs + 1):length(index_vec)]
  return(c(rbind(full_vec, loo_vec)))
}
## (2) create first panel: violin plots of full, leave-out when excluding G. vaginalis
g_vaginalis_plot_tibble <- g_vaginalis_tibble %>% 
  mutate(plot_levels_str = stringr::str_c(taxon_name, "_", type),
         plot_factor_levels = factor(plot_levels_str, 
                                     levels = unique(plot_levels_str)[get_adjacent_indices(1:(2*q_obs), q_obs)], # make full and loo next to each other
                                     ordered = TRUE),
         nice_type = ifelse(type == "full", "Full", "Leave out"),
         wide_taxon_name = stringr::str_pad(taxon_name, width = str_width - str_length("Full") + str_length(taxon_name), side = "right"),
         plot_labels_str = stringr::str_pad(stringr::str_c(wide_taxon_name, nice_type), width = str_width, side = "left"),
         plot_labels = factor(plot_labels_str,
                              levels = unique(plot_labels_str)[get_adjacent_indices(1:(2*q_obs), q_obs)],
                              ordered = TRUE)) %>% 
  select(-plot_levels_str, -plot_labels_str, -wide_taxon_name)
g_vaginalis_plot <- g_vaginalis_plot_tibble %>%  # make levels correspond to taxon number
  ggplot(aes(y = est_eff, x = factor(plot_labels, levels = rev(levels(plot_labels))), 
             group = plot_factor_levels)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  scale_x_discrete(labels = ifelse(grepl("Leave out", rev(levels(g_vaginalis_plot_tibble$plot_labels))),
                                   stringr::str_pad("Leave out", width = str_width, side = "left"), 
                                   rev(levels(g_vaginalis_plot_tibble$plot_labels)))) +
  ggtitle("Distribution of estimated efficiencies with G. vaginalis excluded") +
  ylab("Efficiency") + 
  xlab("") +
  labs(fill = "Taxon") +
  labs(alpha = "Analysis") +
  geom_hline(yintercept = 1, linetype = "dashed", col = "red") +
  guides(alpha = guide_legend(override.aes = list(fill = "blue"))) +
  coord_flip() +
  theme(axis.text.y = element_text(hjust = 1),
        legend.position = c(0.5, 0.7),
        text = element_text(size = main_font_size),
        axis.text = element_text(size = main_font_size),
        plot.title = element_text(size = main_font_size - 8),
        plot.margin = unit(c(0, 1, 0, 0), "cm")) # remove y-axis text
## the full plot
g_vaginalis_full_plot <- ggdraw(g_vaginalis_plot) +
  draw_label("Taxon", size = title_font_size, x = 0.18, y = 0.975, fontface = "bold") +
  draw_label("Dataset", size = title_font_size, x = 0.31, y = 0.975, fontface = "bold")
ggsave(filename = paste0(plots_dir, "distribution_of_efficiencies_gvag_adjust_", adjust, ".png"),
       plot = g_vaginalis_full_plot,
       width = 18, height = 10)


## (2) create second panel: violin plots of full, leave-out when excluding BVAB2
bvab2_plot_tibble <- bvab2_tibble %>% 
  mutate(plot_levels_str = stringr::str_c(taxon_name, "_", type),
         plot_factor_levels = factor(plot_levels_str, 
                                     levels = unique(plot_levels_str)[get_adjacent_indices(1:(2*q_obs), q_obs)], # make full and loo next to each other
                                     ordered = TRUE),
         nice_type = ifelse(type == "full", "Full", "Leave out"),
         wide_taxon_name = stringr::str_pad(taxon_name, width = str_width - str_length("Full") + str_length(taxon_name), side = "right"),
         plot_labels_str = stringr::str_pad(stringr::str_c(wide_taxon_name, nice_type), width = str_width, side = "left"),
         plot_labels = factor(plot_labels_str,
                              levels = unique(plot_labels_str)[get_adjacent_indices(1:(2*q_obs), q_obs)],
                              ordered = TRUE)) %>% 
  select(-plot_levels_str, -plot_labels_str, -wide_taxon_name)
bvab2_plot <- bvab2_plot_tibble %>%  # make levels correspond to taxon number
  ggplot(aes(y = est_eff, x = factor(plot_labels, levels = rev(levels(plot_labels))), 
             group = plot_factor_levels)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  scale_x_discrete(labels = ifelse(grepl("Leave out", rev(levels(bvab2_plot_tibble$plot_labels))),
                                   stringr::str_pad("Leave out", width = str_width, side = "left"), 
                                   rev(levels(bvab2_plot_tibble$plot_labels)))) +
  ggtitle("Distribution of estimated efficiencies with BVAB2 spp. excluded") +
  ylab("Efficiency") + 
  xlab("") +
  labs(fill = "Taxon") +
  labs(alpha = "Analysis") +
  geom_hline(yintercept = 1, linetype = "dashed", col = "red") +
  guides(alpha = guide_legend(override.aes = list(fill = "blue"))) +
  coord_flip() +
  theme(axis.text.y = element_text(hjust = 1),
        legend.position = c(0.5, 0.7),
        text = element_text(size = main_font_size),
        axis.text = element_text(size = main_font_size),
        plot.title = element_text(size = main_font_size - 8),
        plot.margin = unit(c(0, 1, 0, 0), "cm")) # remove y-axis text
## the full plot
bvab2_full_plot <- ggdraw(bvab2_plot) +
  draw_label("Taxon", size = title_font_size, x = 0.18, y = 0.975, fontface = "bold") +
  draw_label("Dataset", size = title_font_size, x = 0.31, y = 0.975, fontface = "bold")
ggsave(filename = paste0(plots_dir, "distribution_of_efficiencies_bvab2_adjust_", adjust, ".png"),
       plot = bvab2_full_plot,
       width = 18, height = 10)

