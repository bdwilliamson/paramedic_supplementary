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
source(paste0(code_dir, "data_analysis/R/get_most_abundant_taxa.R"))

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

## for backwards compatability with R version 3.4.3 (where sims are run)
RNGkind(sample.kind = "Rounding")
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

  ## extract the data inputted to stan for estimators
  taxa_of_interest <- 7 # always the final one in the estimates
  W <- no_ves_i$stan_data_lst$W
  V <- no_ves_i$stan_data_lst$V
  m <- rowSums(W)
  samp_i <- rownames(W)
  full_v <- qpcr_mat[samp, ]
  ## -----------------------------------------------------------------------------------
  ## extract the estimates of qPCR and prediction intervals for qPCR
  ## -----------------------------------------------------------------------------------

  ## naive estimator
  naive_qpcr_i <- naive_estimator(idx = taxa_of_interest, brs = W, qpcrs = V, known_qpcr = 1:6)
  ## proposed estimator without and with varying efficiency
  no_ve_samps_i <- tryCatch(rstan::extract(no_ves_i$stan_out), error = function(e) NA) # purposeful NA if stan model wasn't saved
  ve_samps_i <- tryCatch(rstan::extract(ves_i$stan_out), error = function(e) NA)
  no_ve_qpcr_lst_i <- extract_posterior_summaries(no_ves_i$mod, no_ve_samps_i, q = q, taxa_of_interest, mult_num = mult_num, level = 0.95, interval_type = "wald")
  ve_qpcr_lst_i <- extract_posterior_summaries(ves_i$mod, ve_samps_i, q = q, taxa_of_interest, mult_num = mult_num, level = 0.95, interval_type = "wald")

  ## get mspe
  ests <- cbind(naive_qpcr_i, no_ve_qpcr_lst_i$estimates, ve_qpcr_lst_i$estimates)
  spe_mat <- (ests - full_v[, i]*mult_num)^2
  mspe <- colMeans(spe_mat)

  ## get prediction intervals, coverage
  naive_intervals <- generate_naive_intervals_single_subject(v = V, w = W, m = m, indx = 7, q_obs = 6, level = 0.05)
  no_ve_pis_wald <- gen_wald_interval(no_ve_qpcr_lst_i$estimates, no_ve_qpcr_lst_i$sd, alpha = 0.05)
  ve_pis_wald <- gen_wald_interval(ve_qpcr_lst_i$estimates, ve_qpcr_lst_i$sd, alpha = 0.05)

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
loo_perf_init$loo_index <- 1:7

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
no_ve_full_results <- readRDS(paste0(results_dir, analysis_name, "qpcr_data_analysis_est_no_ve_q_", q, "_q_obs_", q_obs, "_sample_", samp_num, "_loo_", 999, ".rds", sep = ""))
ve_full_results <- readRDS(paste0(results_dir, analysis_name, "qpcr_data_analysis_est_ve_q_", q, "_q_obs_", q_obs, "_sample_", samp_num, "_loo_", 999, ".rds", sep = ""))
no_ve_full_samps <- rstan::extract(no_ve_full_results$stan_out)
ve_full_samps <- rstan::extract(ve_full_results$stan_out)
no_ve_qpcr_lst_full <- extract_posterior_summaries(no_ve_full_results$mod, no_ve_full_samps, q = q, 1:7, mult_num = mult_num, level = 0.95, interval_type = "quantile")
ve_qpcr_lst_full <- extract_posterior_summaries(ve_full_results$mod, ve_full_samps, q = q, 1:7, mult_num = mult_num, level = 0.95, interval_type = "quantile")

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
pred_interval_plot <- loo_perf %>%
  filter(cover_type == "wald_cover") %>%
  ggplot(aes(x = factor(loo_index, labels = c("G. vaginalis", "L. crispatus",
                                              "L. iners", "L. jensenii",
                                              "M. type 1",
                                              "BVAB2 spp.", "A. vaginae")),
             y = cover, group = cover_groups,
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), labels = c("Naive", "Proposed, no ve", "Proposed, ve")))) +
  xlab("Left-out taxon") +
  ylab("Coverage") +
  ggtitle("Prediction interval coverage") +
  labs(shape = "Estimator type") +
  geom_point(size = 4, position = position_dodge(width = 0.3, preserve = "total")) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0, 1)) +
  theme(legend.position = c(0.6, 0.35),
        text = element_text(size = text_size),
        plot.title = element_text(size = text_size*4/3),
        axis.text = element_text(size = text_size*3/2),
        axis.title = element_text(size = text_size*3/2),
        axis.text.x = element_text(face = "italic", angle = 310, hjust = 0),
        legend.text = element_text(size = text_size*4/3),
        legend.title = element_text(size = text_size*3/2),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
pred_interval_plot

## plot mspe against left-out taxon (naive, no ve, ve)
mspe_plot <- loo_perf %>%
  ggplot(aes(x = factor(loo_index, labels = c("G. vaginalis", "L. crispatus",
                                              "L. iners", "L. jensenii",
                                              "M. type 1",
                                              "BVAB2 spp.", "A. vaginae")), y = rmspe,
             group = mspe_groups,
             shape = factor(estimator, levels = c("naive", "no_ve", "ve"), labels = c("Naive", "Proposed, no ve", "Proposed, ve")))) +
  xlab("Left-out taxon") +
  ylab("RMSPE") +
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
        axis.text.x = element_text(face = "italic", angle = 310, hjust = 0),
        legend.text = element_text(size = text_size),
        plot.margin = unit(c(0, 1.75, 0, 0), "cm"))
mspe_plot

## save the plot
ggsave(paste0(plots_dir, "/loo_perf.png"),
       plot = plot_grid(pred_interval_plot, mspe_plot),
       width = fig_width, height = 10)

## plot width against left-out index for each type of interval
width_plot <- widths %>%
  ggplot(aes(x = loo_index, y = log10(width), group = paste(estimator, measure, loo_index, sep = "_"),
             fill = factor(estimator, levels = c("naive", "no_ve", "ve"), labels = c("Naive", "Proposed, no ve", "Proposed, ve")))) +
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
ggsave(paste0(plots_dir, "/loo_widths.png"),
       plot = width_plot,
       width = fig_width, height = fig_height)

## plot estimated efficiencies
## violin plot with two panels
main_font_size <- 30
main_font_size_lab <- 10
title_font_size <- 26

## make a tibble for estimated efficiencies
eff_samps_lst <- lapply(leave_one_out_performance, function(x) x$efficiency_samps)
full_est_efficiency <- ve_qpcr_lst_full$est_efficiency
full_efficiency_samps <- ve_full_samps$e
colnames(full_efficiency_samps) <- colnames(qpcr_mat)
sd(full_est_efficiency)

## (1) create tibbles with samples and the correct taxon index, latin name
## want the full samples for each taxon, both when leaving out the named one and not
g_vaginalis_tibble <- tibble(est_eff = c(full_efficiency_samps[, "gvag_cps"], full_efficiency_samps[, "lcrisp_cps"],
                                         full_efficiency_samps[, "liners_cps"], full_efficiency_samps[, "ljensii_cps"],
                                         full_efficiency_samps[, "mega_cps"], full_efficiency_samps[, "bvab2_cps"],
                                         full_efficiency_samps[, "avag_cps"],
                                         eff_samps_lst$LOO_gvag_cps[, 7], eff_samps_lst$LOO_gvag_cps[, 1],
                                         eff_samps_lst$LOO_gvag_cps[, 2], eff_samps_lst$LOO_gvag_cps[, 3],
                                         eff_samps_lst$LOO_gvag_cps[, 4], eff_samps_lst$LOO_gvag_cps[, 5],
                                         eff_samps_lst$LOO_gvag_cps[, 6]),
                             taxon_name = rep(c(rep("G. vaginalis", dim(full_efficiency_samps)[1]),
                                                rep("L. crispatus", dim(full_efficiency_samps)[1]),
                                                rep("L. iners", dim(full_efficiency_samps)[1]),
                                                rep("L. jensenii", dim(full_efficiency_samps)[1]),
                                                rep("M. type 1", dim(full_efficiency_samps)[1]),
                                                rep("BVAB2 spp.", dim(full_efficiency_samps)[1]),
                                                rep("A. vaginae", dim(full_efficiency_samps)[1])), 2),
                             type = c(rep("full", dim(full_efficiency_samps)[1]*7),
                                      rep("loo", dim(full_efficiency_samps)[1]*7)),
                             iter = rep(1:dim(full_efficiency_samps)[1], 7*2))
bvab2_tibble <- tibble(est_eff = c(full_efficiency_samps[, "gvag_cps"], full_efficiency_samps[, "lcrisp_cps"],
                                         full_efficiency_samps[, "liners_cps"], full_efficiency_samps[, "ljensii_cps"],
                                         full_efficiency_samps[, "mega_cps"], full_efficiency_samps[, "bvab2_cps"],
                                         full_efficiency_samps[, "avag_cps"],
                                         eff_samps_lst$LOO_bvab2_cps[, 1], eff_samps_lst$LOO_bvab2_cps[, 2],
                                         eff_samps_lst$LOO_bvab2_cps[, 3], eff_samps_lst$LOO_bvab2_cps[, 4],
                                         eff_samps_lst$LOO_bvab2_cps[, 5], eff_samps_lst$LOO_bvab2_cps[, 7],
                                         eff_samps_lst$LOO_bvab2_cps[, 6]),
                             taxon_name = rep(c(rep("G. vaginalis", dim(full_efficiency_samps)[1]),
                                                rep("L. crispatus", dim(full_efficiency_samps)[1]),
                                                rep("L. iners", dim(full_efficiency_samps)[1]),
                                                rep("L. jensenii", dim(full_efficiency_samps)[1]),
                                                rep("M. type 1", dim(full_efficiency_samps)[1]),
                                                rep("BVAB2 spp.", dim(full_efficiency_samps)[1]),
                                                rep("A. vaginae", dim(full_efficiency_samps)[1])), 2),
                             type = c(rep("full", dim(full_efficiency_samps)[1]*7),
                                      rep("loo", dim(full_efficiency_samps)[1]*7)),
                             iter = rep(1:dim(full_efficiency_samps)[1], 7*2))

## (2) create first panel: violin plots of full, leave-out when excluding G. vaginalis
g_vaginalis_plot <- g_vaginalis_tibble %>%
  mutate(plot_factor_levels = factor(paste0(taxon_name, "_", type),
                                     levels = unique(paste0(taxon_name, "_", type))[c(1, 8, 2, 9, 3, 10, 4, 11, 5, 12, 6, 13, 7, 14)],
                                     ordered = TRUE)) %>%  # make levels correspond to taxon number
  ggplot(aes(y = est_eff, x = factor(plot_factor_levels, levels = rev(levels(plot_factor_levels))),
             group = plot_factor_levels)) +
  coord_flip() +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  ggtitle("Distribution of estimated efficiencies with G. vaginalis excluded") +
  ylab("Efficiency") +
  xlab("") +
  labs(fill = "Taxon") +
  labs(alpha = "Analysis") +
  geom_hline(yintercept = 1, linetype = "dashed", col = "red") +
  guides(alpha = guide_legend(override.aes = list(fill = "blue"))) +
  theme(axis.text.y = element_blank(),
        legend.position = c(0.5, 0.7),
        text = element_text(size = main_font_size),
        axis.text = element_text(size = main_font_size),
        plot.title = element_text(size = main_font_size - 8),
        plot.margin = unit(c(0, 1, 0, 0), "cm")) # remove y-axis text
## make the labels nice
g_vaginalis_labels <- g_vaginalis_tibble %>%
  select(taxon_name, type) %>%
  group_by(taxon_name, type) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(type = ifelse(type == "full", "Full", "Leave out"),
         var = 1) %>%
  dplyr::arrange(taxon_name)
g_vaginalis_labs <- reshape2::melt(g_vaginalis_labels, id.var = "var")
g_vaginalis_labs$x_coord <- rep(c(0, 0.025), each = 14)
g_vaginalis_labs$y_coord <- rep(c(2, 1, 4, 3, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5), 2)
g_vaginalis_labs$value[c(2, 4, 6, 8, 10, 12, 14)] <- ""
## plot the labels
g_vaginalis_label_plot <- g_vaginalis_labs %>%
  ggplot(aes(x = x_coord, y = y_coord, label = value)) +
  geom_text(size = main_font_size_lab, hjust = 0, vjust = 0.15) +
  xlim(c(0, 0.04)) +
  ylim(c(1, 14)) +
  theme(legend.position = "",
        axis.line = element_blank(),
        axis.text = element_blank(),
        text = element_text(size = main_font_size_lab),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
## the full plot
g_vaginalis_full_plot <- plot_grid(g_vaginalis_label_plot, g_vaginalis_plot, align = "h",
                                   rel_widths = c(1, 2.1)) +
  draw_label("Taxon", size = title_font_size, x = 0.05, y = 0.97, fontface = "bold") +
  draw_label("Dataset", size = title_font_size, x = 0.24, y = 0.97, fontface = "bold")
ggsave(filename = paste0(plots_dir, "distribution_of_efficiencies_gvag.png"),
       plot = g_vaginalis_full_plot,
       width = 14, height = 10)


## (2) create second panel: violin plots of full, leave-out when excluding BVAB2
bvab2_plot <- bvab2_tibble %>%
  mutate(plot_factor_levels = factor(paste0(taxon_name, "_", type),
                                     levels = unique(paste0(taxon_name, "_", type))[c(1, 8, 2, 9, 3, 10, 4, 11, 5, 12, 6, 13, 7, 14)],
                                     ordered = TRUE)) %>%  # make levels correspond to taxon number
  ggplot(aes(y = est_eff, x = factor(plot_factor_levels, levels = rev(levels(plot_factor_levels))),
             group = plot_factor_levels)) +
  coord_flip() +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  scale_alpha_discrete(range = c(0.5, 1)) +
  ggtitle("Distribution of estimated efficiencies with BVAB2 spp. excluded") +
  ylab("Efficiency") +
  xlab("") +
  labs(fill = "Taxon") +
  labs(alpha = "Analysis") +
  geom_hline(yintercept = 1, linetype = "dashed", col = "red") +
  guides(alpha = guide_legend(override.aes = list(fill = "blue"))) +
  theme(axis.text.y = element_blank(),
        legend.position = c(0.5, 0.7),
        text = element_text(size = main_font_size),
        axis.text = element_text(size = main_font_size),
        plot.title = element_text(size = main_font_size - 8),
        plot.margin = unit(c(0, 1, 0, 0), "cm")) # remove y-axis text
## make the labels nice
bvab2_labels <- bvab2_tibble %>%
  select(taxon_name, type) %>%
  group_by(taxon_name, type) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(type = ifelse(type == "full", "Full", "Leave out"),
         var = 1) %>%
  dplyr::arrange(taxon_name)
bvab2_labs <- reshape2::melt(bvab2_labels, id.var = "var")
bvab2_labs$x_coord <- rep(c(0, 0.025), each = 14)
bvab2_labs$y_coord <- rep(c(2, 1, 4, 3, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5), 2)
bvab2_labs$value[c(2, 4, 6, 8, 10, 12, 14)] <- ""
## plot the labels
bvab2_label_plot <- bvab2_labs %>%
  ggplot(aes(x = x_coord, y = y_coord, label = value)) +
  geom_text(size = main_font_size_lab, hjust = 0, vjust = 0.15) +
  xlim(c(0, 0.04)) +
  ylim(c(1, 14)) +
  theme(legend.position = "",
        axis.line = element_blank(),
        axis.text = element_blank(),
        text = element_text(size = main_font_size_lab),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
## the full plot
bvab2_full_plot <- plot_grid(bvab2_label_plot, bvab2_plot, align = "h",
                             rel_widths = c(1, 2.1)) +
  draw_label("Taxon", size = title_font_size, x = 0.05, y = 0.97, fontface = "bold") +
  draw_label("Dataset", size = title_font_size, x = 0.24, y = 0.97, fontface = "bold")
ggsave(filename = paste0(plots_dir, "distribution_of_efficiencies_bvab2.png"),
       plot = bvab2_full_plot,
       width = 14, height = 10)
