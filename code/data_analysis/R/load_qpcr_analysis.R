######################################################################################
##
## FILE: load_qpcr_analysis.R
##
## CREATED: 27 November 2018 by Brian Williamson
##
## PURPOSE: data analysis of the qPCR + br16s data from Hutch collaborators
##          to showcase qPCR estimation method
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
source(paste0(code_dir, "data_analysis/R/qpcr_analysis_helper_functions.R"))
source(paste0(code_dir, "naive_qpcr_estimator.R"))
source(paste0(code_dir, "data_analysis/R/get_most_abundant_taxa.R"))
library("paramedic")
library("cowplot")

## set up arguments
q <- 17
q_obs <- 7
samp_num <- 1213
do_sample <- TRUE
analysis_name <- "data_analysis/all_women_q_17"

## -----------------------------------------------------------------------------------
## load in the data, clean it up
## -----------------------------------------------------------------------------------
full_data <- read.csv(file=paste0(data_dir, "p2_brngs_qPCR_merge_20180314.csv"))

## set m_min, llod for processing
m_min <- 1000
llod <- 0
div_num <- 1000 # divide qpcr by 1000, setting llod to (essentially) zero during processing
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

## set q, q_obs
q_obs <- dim(qpcr_mat)[2]

## if q is smaller than the total number, select the most abundant taxa
if (q < dim(br16_mat)[2]) {
  ## get the order
  ordered_by_abundance <- get_most_abundant_taxa(br16_mat, m)
  ## remove the ones corresponding to observed qPCR
  taxa_to_estimate <- ordered_by_abundance[!(ordered_by_abundance %in% 1:7)]
  ## select the most abundant taxa (always select the first 7, corresponding to qPCR)
  most_abundant_16S <- br16_mat[, c(1:7, taxa_to_estimate[1:(q - q_obs)])]
  ## rename br16_mat
  br16_mat <- most_abundant_16S
  ## re-normalize
  m <- rowSums(br16_mat)
}

## for backwards compatability with R version 3.4.3 (where sims are run)
RNGkind(sample.kind = "Rounding")
set.seed(4747)
samp <- sample(1:dim(br16_mat)[1], samp_num)

## -----------------------------------------------------------------------------------
## load in the estimators (all women)
## -----------------------------------------------------------------------------------

## read in the datasets
naive <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_centered_est_naive_q_", q, "_q_obs_", q_obs, "_sample_", samp_num, "_loo_999.rds", sep = ""))
no_ves <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_centered_est_no_ve_q_", q, "_q_obs_", q_obs, "_sample_", samp_num, "_loo_999.rds", sep = ""))
ves <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_centered_est_ve_q_", q, "_q_obs_", q_obs, "_sample_", samp_num, "_loo_999.rds", sep = ""))

## -----------------------------------------------------------------------------------
## extract the estimates of qPCR and prediction intervals for qPCR
## report everything in 1000s (don't multiply by div_num)
## -----------------------------------------------------------------------------------
## taxa of interest
taxa_of_interest <- (q_obs + 1):q
## naive estimator
naive_qpcr <- naive$mod[, taxa_of_interest]
## proposed estimator without and with varying efficiency
no_ve_samps <- no_ves$samps
ve_samps <- ves$samps
no_ve_qpcr_lst <- extract_posterior_summaries(no_ves$mod, no_ve_samps, q = q, taxa_of_interest, mult_num = 1, level = 0.95)
ve_qpcr_lst <- extract_posterior_summaries(ves$mod, ve_samps, q = q, taxa_of_interest, mult_num = 1, level = 0.95)

## check convergence diagnostics
summary(no_ves$mod)
summary(ves$mod)

## check for efficiency only
summary(ves$mod[grepl("e", rownames(ves$mod)) & !grepl("beta", rownames(ves$mod)), ])
ves$mod[grepl("e", rownames(ves$mod)) & !grepl("beta", rownames(ves$mod)), ]
ves$mod[grepl("sigma", rownames(ves$mod)), ]

## -----------------------------------------------------------------------------------
## summarize efficiencies
## -----------------------------------------------------------------------------------
## table of point estimates, credible intervals
efficiency_table <- ves$mod[grepl("e", rownames(ves$mod)) & !grepl("beta", rownames(ves$mod)), c(1, 4, 5)]
sigma <- ves$mod[grepl("sigma", rownames(ves$mod)), c(1, 4, 5)]
print(xtable::xtable(efficiency_table, caption = paste0("Posterior means and 95\\% credible intervals for the efficiencies $e_j$, $j = 1, \\dots, 17$, based on all 1213 women. The posterior mean for $\\sigma$ is ", round(sigma[1], 2), " with a 95\\% credible interval of [", round(sigma[2], 2), ", ", round(sigma[3], 2), "]."),
                     label = "tab:efficiencies"),
      file = paste0("tables/", analysis_name, "/efficiencies.tex"))

## -----------------------------------------------------------------------------------
## summarize qpcr, for each taxon
## -----------------------------------------------------------------------------------
naive_taxon_means <- colMeans(naive_qpcr)
no_ve_taxon_means <- colMeans(no_ve_qpcr_lst$estimates)
ve_taxon_means <- colMeans(ve_qpcr_lst$estimates)

naive_taxon_means
no_ve_taxon_means
ve_taxon_means

## differences between the two
taxon_differences <- naive_taxon_means - no_ve_taxon_means
taxon_differences_ve <- naive_taxon_means - ve_taxon_means
taxon_differences_ve_vs_no <- no_ve_taxon_means - ve_taxon_means

summary(taxon_differences)
summary(taxon_differences_ve)
summary(taxon_differences_ve_vs_no)

## -----------------------------------------------------------------------------------
## create output tables/plots from qpcr
## -----------------------------------------------------------------------------------
## table with point estimates, credible intervals, prediction intervals for all women, all taxa
## only for the variable efficiency estimator
all_taxa_ests <- extract_posterior_summaries(ves$mod, ve_samps, q = q, taxa_of_interest = 1:q, mult_num = 1, level = 0.95)
est_vec <- as.vector(all_taxa_ests$estimates)
## add on observed qPCR for those that we have
obs_vec <- as.vector(cbind(qpcr_mat, matrix(NA, nrow = nrow(qpcr_mat), ncol = ncol(all_taxa_ests$estimates) - ncol(qpcr_mat))))
cred_vec <- apply(all_taxa_ests$cred_intervals, 2, c)
pred_vec <- apply(all_taxa_ests$pred_intervals, 2, c)
est_table_init <- tibble::tibble(est = est_vec, obs = obs_vec, cred_int_lower = cred_vec[, 1], cred_int_upper = cred_vec[, 2],
                                 pred_int_lower = pred_vec[, 1], pred_int_upper = pred_vec[, 2]) %>%
  mutate(taxon = rep(1:q, each = samp_num), subj_id = rep(1:samp_num, q), log_est = log10(est)) %>%
  group_by(taxon, subj_id)
est_table_init[est_table_init$taxon < q_obs + 1, c("pred_int_lower", "pred_int_upper")] <- NA

est_table <- est_table_init[, c(7, 8, 1:6, 9)]
readr::write_csv(est_table[, 1:8], path = paste0("tables/", analysis_name, "/est_concentrations.csv"))

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
## correct labels for printing
taxon_labels <- c("G. vaginalis", "L. crispatus", "L. iners", "L. jensenii", "M. genomosp. type 1",
                  "BVAB2 spp.", "A. vaginae", "BVAB1 spp.", "S. amnii", "L. gasseri johnsonii",
                  "S. sanguinegens", "P. timonensis", "L. gasseri", "P. amnii",
                  "Eggerthella spp. type 1", "Dialister spp. type 2",
                  "P. bivia")

## create directory for plots if it doesn't already exist
plots_dir <- paste0("plots/", analysis_name)
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}
## plot estimated mu/qPCR for all women
est_plot <- est_table %>%
  ggplot(aes(x = log_est, group = taxon, fill = factor(taxon, labels = taxon_labels, ordered = FALSE))) +
  geom_histogram(binwidth = .025) +
  scale_fill_manual(values = cols) +
  labs(x = expression(log[10](hat(mu)))) +
  labs(y = "Frequency") +
  labs(fill = "Taxon") +
  ggtitle("Estimated concentrations: each taxon, all women") +
  xlim(c(0, 12)) +
  ylim(c(0, 400)) +
  theme_bw() +
  theme(text = element_text(size = 20), legend.text = element_text(face = "italic")) +
  geom_vline(xintercept = 1, col = "red", linetype = "dashed")
png(paste0(plots_dir, "/est_concentrations.png"), width = fig_width, height = fig_height, res = 300, units = "px")
est_plot
dev.off()

## -----------------------------------------
## heatmap of all taxa for a subset of women
## -----------------------------------------
## sample 40 women
set.seed(4747)
samp_women <- sample(1:dim(qpcr_mat)[1], 40, replace = FALSE)
## get the subset
table_for_heatmap <- est_table %>%
  filter(subj_id %in% samp_women)
table_for_heatmap$subj_id <- rep(1:length(samp_women), length(unique(table_for_heatmap$taxon)))
## plot a heatmap of log mu for these women, all taxa
concentration_heatmap <- table_for_heatmap %>%
  mutate(normalized_log_est = log_est/max(table_for_heatmap$log_est)) %>%
  ggplot(aes(y = subj_id, x = factor(taxon, ordered = FALSE, labels = taxon_labels), fill = normalized_log_est)) +
  geom_tile(color = "white") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(y = "Subject ID") +
  labs(x = "") +
  labs(fill = expression(paste("Normalized ", log[10](hat(mu)), sep = ""))) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 330, hjust = 0, face = "italic"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) # top, right, bottom, left

## ----------------------------------------------------------
## plot recreating Fig 1 but with estimated conc. vs obs qPCR
## ----------------------------------------------------------
## qpcr proportions; need to reorder by samp to match up with estimates
compositional_qpcr_init <- as.data.frame(qpcr_mat) %>%
  mutate(row_sum = rowSums(.)) %>%
  mutate_at(vars(-row_sum), list(~`/`(., row_sum))) %>%
  select(-row_sum)
compositional_qpcr <- compositional_qpcr_init[samp, ]
compositional_qpcr$id <- 1:dim(compositional_qpcr)[1]
compositional_qpcr_long <- compositional_qpcr %>%
  gather(key = taxon, value = proportion, -id)

## concentration proportions
obs_taxa_ests <- extract_posterior_summaries(ves$mod, ve_samps, q = q, taxa_of_interest = 1:q_obs, mult_num = 1, level = 0.95)
colnames(obs_taxa_ests$estimates) <- colnames(qpcr_mat)
compositional_conc <- as.data.frame(obs_taxa_ests$estimates*div_num) %>%
  mutate(row_sum = rowSums(.)) %>%
  mutate_at(vars(-row_sum), list(~`/`(., row_sum))) %>%
  select(-row_sum)
compositional_conc$id <- 1:dim(compositional_conc)[1]
compositional_conc_long <- compositional_conc %>%
  gather(key = taxon, value = proportion, -id)

obs_taxa_ests_no_ve <- extract_posterior_summaries(no_ves$mod, no_ve_samps, q = q, taxa_of_interest = 1:q_obs, mult_num = 1, level = 0.95)
colnames(obs_taxa_ests_no_ve$estimates) <- colnames(qpcr_mat)
compositional_conc_no_ve <- as.data.frame(obs_taxa_ests_no_ve$estimates*div_num) %>%
  mutate(row_sum = rowSums(.)) %>%
  mutate_at(vars(-row_sum), list(~`/`(., row_sum))) %>%
  select(-row_sum)
compositional_conc_no_ve$id <- 1:dim(compositional_conc_no_ve)[1]
compositional_conc_long_no_ve <- compositional_conc_no_ve %>%
  gather(key = taxon, value = proportion, -id)

## create dataset for plotting
compositional_data <- dplyr::left_join(compositional_conc_long, compositional_qpcr_long, by = c("id", "taxon"), suffix = c(".est_conc", ".qPCR")) %>%
  select(id, taxon, proportion.qPCR, proportion.est_conc)

compositional_data_no_ve <- dplyr::left_join(compositional_conc_long_no_ve, compositional_qpcr_long, by = c("id", "taxon"), suffix = c(".est_conc", ".qPCR")) %>%
  select(id, taxon, proportion.qPCR, proportion.est_conc)

## create the plot
compositional_data %$% taxon %>% unique # print the taxon names currently
taxon_labels[1:q_obs] # print what I want them to be
compositional_plot <- compositional_data %>%
  mutate(taxon = stringr::str_replace(taxon, "_cps", "")) %>%
  mutate(taxon = stringr::str_replace(taxon, "gvag", taxon_labels[1])) %>%
  mutate(taxon = stringr::str_replace(taxon, "lcrisp", taxon_labels[2])) %>%
  mutate(taxon = stringr::str_replace(taxon, "liners", taxon_labels[3])) %>%
  mutate(taxon = stringr::str_replace(taxon, "ljensii", taxon_labels[4])) %>%
  mutate(taxon = stringr::str_replace(taxon, "mega", taxon_labels[5])) %>%
  mutate(taxon = stringr::str_replace(taxon, "bvab2", taxon_labels[6])) %>%
  mutate(taxon = stringr::str_replace(taxon, "avag", taxon_labels[7])) %>%
  ggplot(aes(x = proportion.qPCR, y = proportion.est_conc, color = taxon)) +
  geom_point() +
  scale_color_manual(values = cols) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Proportion within subcomposition: qPCR", y = "Proportion within subcomposition:\n estimated concentration", color = "Taxon") +
  theme(text = element_text(size = 15), legend.text=element_text(size = 15))

compositional_data_no_ve %>%
  mutate(taxon = stringr::str_replace(taxon, "_cps", "")) %>%
  mutate(taxon = stringr::str_replace(taxon, "gvag", taxon_labels[1])) %>%
  mutate(taxon = stringr::str_replace(taxon, "lcrisp", taxon_labels[2])) %>%
  mutate(taxon = stringr::str_replace(taxon, "liners", taxon_labels[3])) %>%
  mutate(taxon = stringr::str_replace(taxon, "ljensii", taxon_labels[4])) %>%
  mutate(taxon = stringr::str_replace(taxon, "mega", taxon_labels[5])) %>%
  mutate(taxon = stringr::str_replace(taxon, "bvab2", taxon_labels[6])) %>%
  mutate(taxon = stringr::str_replace(taxon, "avag", taxon_labels[7])) %>%
  ggplot(aes(x = proportion.qPCR, y = proportion.est_conc, color = taxon)) +
  geom_point() +
  scale_color_manual(values = cols) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Proportion within subcomposition: qPCR", y = "Proportion within subcomposition:\n estimated concentration", color = "Taxon") +
  theme(text = element_text(size = 15), legend.text=element_text(size = 15))

compositional_plot

## ----------------------------------------------------------
## save off the plot
## ----------------------------------------------------------
ggsave(paste0(plots_dir, "/concentration_full_fig.png"),
       plot = plot_grid(concentration_heatmap, compositional_plot, labels = c("A", "B")),
       width = 18, height = 8)


## spike in L. crisp?
as.data.frame(obs_taxa_ests$estimates*div_num) %>%
  select(lcrisp_cps) %>%
  ggplot(aes(x = lcrisp_cps)) +
  geom_histogram(bins = 200) +
  geom_vline(xintercept = 100, color = "red", linetype = "dashed")
as.data.frame(qpcr_mat) %>%
  select(lcrisp_cps) %>%
  ggplot(aes(x = lcrisp_cps)) +
  geom_histogram(bins = 200) +
  geom_vline(xintercept = 100, color = "red", linetype = "dashed")
as.data.frame(br16_mat) %>%
  select(Lactobacillus.crispatus) %>%
  ggplot(aes(x = Lactobacillus.crispatus)) +
  geom_histogram(bins = 200) +
  geom_vline(xintercept = 100, color = "red", linetype = "dashed")
## check to see if it's the same ones
zero_lcrisp_w <- as.data.frame(br16_mat)[samp, ] %>%
  mutate(id = 1:dim(br16_mat)[1]) %>%
  select(id, Lactobacillus.crispatus) %>%
  filter(Lactobacillus.crispatus == 0)
zero_lcrisp_w_est <- as.data.frame(obs_taxa_ests$estimates*div_num) %>%
  mutate(id = 1:dim(obs_taxa_ests$estimates)[1]) %>%
  select(id, lcrisp_cps) %>%
  filter(id %in% zero_lcrisp_w$id)
dplyr::left_join(zero_lcrisp_w, zero_lcrisp_w_est, by = "id")

summary(obs_taxa_ests$estimates[, "lcrisp_cps"]*div_num)

zero_lcrisp_est <- as.data.frame(obs_taxa_ests$estimates*div_num) %>%
  mutate(id = 1:dim(obs_taxa_ests$estimates)[1]) %>%
  select(lcrisp_cps) %>%
  filter(lcrisp_cps < 10)
