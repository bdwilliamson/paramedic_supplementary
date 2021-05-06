######################################################################################
##
## FILE: load_qpcr_analysis.R
##
## CREATED: 27 November 2018 by Brian Williamson
##
## PURPOSE: data analysis of the qPCR + br16s data from Hutch collaborators
##          to showcase qPCR estimation method
##
## INPUTS: analysis_data.rds
##
## OUTPUTS: qPCR/tables/data_analysis/sample_*_most_abundant_*/est_efficiencies.tex
##          qPCR/tables/data_analysis/sample_*_most_abundant_*/point_ests_plus_intervals.tex
######################################################################################
## -----------------------------------------------------------------------------------
## load required functions and libraries
## -----------------------------------------------------------------------------------
analysis_name <- "data_analysis/full_analysis"
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) { # if running locally
  code_dir <- "code/R/"
  data_dir <- "data/"
  stan_dir <- "stan/"
  results_dir <- "results/"
  tables_dir <- "tables/"
  plots_dir <- paste0("plots/", analysis_name)
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }
} else {
  code_dir <- "../R/"
  data_dir <- "../../data/"
  stan_dir <- "../../stan/"
  results_dir <- "../../results/"
  tables_dir <- "../../tables/"
  plots_dir <- paste0("../../plots/", analysis_name)
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }
}

library("methods")
library("StanHeaders")
library("rstan")
library("argparse")
library("Cairo")
# remotes::install_github("statdivlab/paramedic", ref = "v0.0.3") # local
# withr::with_libpaths("/home/bwillia2/", remotes::install_github("statdivlab/paramedic", ref = "v0.0.2")) # rhino
library("paramedic")
library("tidyverse")
# source(paste0(code_dir, "analyze_data/qpcr_analysis_helper_functions.R"))
source(paste0(code_dir, "naive_qpcr_estimator.R"))
source(paste0(code_dir, "analyze_data/get_most_abundant_taxa.R"))
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())

parser <- ArgumentParser()
parser$add_argument("--q", type = "double", default = 127, help = "number of total taxa")
parser$add_argument("--q-obs", type = "double", default = 13, help = "number of observed qPCR taxa")
parser$add_argument("--sample-num", type = "double", default = 55, help = "sample sample-num women or not")
parser$add_argument("--div-num", type = "double", default = 1000, help = "number to divide qpcr by")
parser$add_argument("--adjust", type = "double", default = 1, help = "adjusted analysis?")
parser$add_argument("--sens", type = "double", default = 0, help = "sensitivity analysis with different alpha_sigma and kappa_sigma?")
args <- parser$parse_args()

# three sets of valid args:
# adjust = 1, sens = 0
# adjust = 0, sens = 0
# adjust = 0, sens = 1

## set up arguments
q <- args$q
q_obs <- args$q_obs
do_sample <- TRUE
sigma_txt <- ifelse(args$sens, "_as_2_ks_3", "")
n <- 110

## -----------------------------------------------------------------------------------
## load in the data, clean it up
## -----------------------------------------------------------------------------------
analysis_data <- readRDS(paste0(data_dir, "analysis_data.rds"))
qpcr <- analysis_data$qpcr %>%
  mutate_all(~ . / args$div_num)
br <- analysis_data$filtered_br16S
case_control <- (analysis_data$case_control %>%
                   mutate(case_bin = as.numeric(case == "Case")) %>%
                   select(-case))$case_bin
data_lst <- list(qpcr = qpcr, br = br, cc = case_control)

qpcr_mat <- as.matrix(qpcr)
br16_mat <- as.matrix(br)
br16_mat_nms <- colnames(br16_mat)
## calculate the read numbers
m <- rowSums(br16_mat)
observed_taxa <- 1:args$q_obs

set.seed(4747)
samp_num <- args$sample_num

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
} else {
  taxa_to_estimate <- (1:q)[-observed_taxa]
}

cat("\n Taxa to estimate: \n")
print(taxa_to_estimate[1:(q - q_obs)])

## -----------------------------------------------------------------------------------
## load in the estimators (no adjustment for now)
## -----------------------------------------------------------------------------------

## read in the datasets
naive_0 <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_est_naive_q_", q, "_q_obs_", q_obs, "_sample_", args$sample_num, "_loo_999", sigma_txt, "_adjust_", args$adjust, ".rds", sep = ""))
no_ves_0 <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_est_no_ve_q_", q, "_q_obs_", q_obs, "_sample_", args$sample_num, "_loo_999", sigma_txt, "_adjust_", args$adjust, ".rds", sep = ""))
ves_0 <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_est_ve_q_", q, "_q_obs_", q_obs, "_sample_", args$sample_num, "_loo_999", sigma_txt, "_adjust_", args$adjust, ".rds", sep = ""))

## -----------------------------------------------------------------------------------
## extract the estimates of qPCR and prediction intervals for qPCR
## report everything in 1000s (don't multiply by div_num)
## -----------------------------------------------------------------------------------
## taxa of interest
taxa_of_interest <- (q_obs + 1):q
## naive estimator
naive_qpcr <- naive_0$mod[, taxa_of_interest]
## proposed estimator without and with varying efficiency
no_ve_samps <- no_ves_0$samps
ve_samps <- ves_0$samps
no_ve_qpcr_lst <- extract_posterior_summaries(no_ves_0$mod, no_ve_samps, taxa_of_interest, mult_num = 1, level = 0.95)
ve_qpcr_lst <- extract_posterior_summaries(ves_0$mod, ve_samps, taxa_of_interest, mult_num = 1, level = 0.95)

## check convergence diagnostics
summary(no_ves_0$mod)
summary(ves_0$mod)

## check for efficiency only
summary(ves_0$mod[grepl("e", rownames(ves_0$mod)) & !grepl("beta", rownames(ves_0$mod)), ])
ves_0$mod[grepl("e", rownames(ves_0$mod)) & !grepl("beta", rownames(ves_0$mod)), ]
ves_0$mod[grepl("sigma", rownames(ves_0$mod)), ]

## -----------------------------------------------------------------------------------
## summarize efficiencies
## -----------------------------------------------------------------------------------
## table of point estimates, credible intervals
efficiency_table <- ves_0$mod[grepl("e", rownames(ves_0$mod)) & !grepl("beta", rownames(ves_0$mod)), c(1, 4, 5)]
sigma <- ves_0$mod[grepl("sigma", rownames(ves_0$mod)), c(1, 4, 5)]
efficiency_table_cap <- paste0("Posterior means and 95\\% credible intervals for the efficiencies $e_j$, $j = 1, \\dots,",
                        q, "$, based on ",
                        ifelse(args$sample_num < nrow(qpcr), 
                        paste0(args$sample_num, " randomly sampled "),
                        "all "), "women",
                        ifelse(args$adjust == 1, ", adjusted for case-control status", ""),
                        ". The posterior mean for $\\sigma$ is ", 
                        round(sigma[1], 2), " with a 95\\% credible interval of [", 
                        round(sigma[2], 2), ", ", round(sigma[3], 2), "].")
print(xtable::xtable(efficiency_table, caption = efficiency_table_cap,
                     label = "tab:efficiencies"),
      file = paste0(tables_dir, analysis_name, "/efficiencies_n_", args$sample_num, sigma_txt, "_adjust_", args$adjust, ".tex"))

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
all_taxa_ests <- extract_posterior_summaries(ves_0$mod, ve_samps, taxa_of_interest = 1:q, mult_num = 1, level = 0.95)
est_vec <- as.vector(all_taxa_ests$estimates)
## add on observed qPCR for those that we have
ids <- ves_0$stan_data_lst$V[, 1]
qpcr_plots_mat <- qpcr_mat[ids, ]
obs_vec <- as.vector(cbind(qpcr_plots_mat, matrix(NA, nrow = nrow(qpcr_plots_mat), ncol = ncol(all_taxa_ests$estimates) - ncol(qpcr_mat))))
cred_vec <- apply(all_taxa_ests$cred_intervals, 2, c)
pred_vec <- apply(all_taxa_ests$pred_intervals, 2, c)
est_table_init <- tibble::tibble(est = est_vec, obs = obs_vec, cred_int_lower = cred_vec[, 1], cred_int_upper = cred_vec[, 2],
                                 pred_int_lower = pred_vec[, 1], pred_int_upper = pred_vec[, 2]) %>%
  mutate(taxon = rep(1:q, each = args$sample_num), subj_id = rep(ids, q), log_est = log10(est)) %>%
  group_by(taxon, subj_id)
est_table_init[est_table_init$taxon < q_obs + 1, c("pred_int_lower", "pred_int_upper")] <- NA

est_table <- est_table_init %>% 
  select(taxon, subj_id, est, obs, cred_int_lower, cred_int_upper, 
         pred_int_lower, pred_int_upper, log_est)
readr::write_csv(est_table[, 1:8], 
                 path = paste0(tables_dir, analysis_name, "/est_concentrations_n_", 
                               args$sample_num, sigma_txt, "_", args$adjust, ".csv"))

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
# obs_taxon_nms <- c("A. christensenii", "A. vaginae",
#                             "BVAB2 spp.", "D. micraerophilus",
#                             "E. spp. type 1",
#                             "G. vaginalis", "L. crispatus",
#                             "L. iners", "L. jensenii",
#                             "M. hominis", "P. micra",
#                             "P. bennonis", "P. spp. type 1")
taxon_nms <- tibble(taxon = names(br)) %>% 
  mutate(taxon = str_replace(taxon, "\\.", " ")) %>%   
  # the different ways to code spaces
  mutate(taxon = str_replace(taxon, "\\_", ". ")) %>% 
  mutate(taxon = str_replace(taxon, "\\..", ". ")) %>% 
  # species abbreviations
  mutate(taxon = str_replace(taxon, "ssp. t", "s sp. t")) %>%  
  mutate(taxon = str_replace(taxon, "sptype", "sp. type")) %>% 
  mutate(taxon = str_replace(taxon, "erasp", "era sp")) %>% 
  # the different ways that names are pasted together
  mutate(taxon = str_replace(taxon, "coccus{1}", "coccus ")) %>% 
  mutate(taxon = str_replace(taxon, "ella{1}", "ella ")) %>% 
  mutate(taxon = str_replace(taxon, "monas{1}", "monas ")) %>% 
  mutate(taxon = str_replace(taxon, "type{1}", "type ")) %>% 
  mutate(taxon = str_replace(taxon, "chia{1}", "chia ")) %>% 
  mutate(taxon = str_replace(taxon, "thia{1}", "thia ")) %>% 
  mutate(taxon = str_replace(taxon, "illus{1}", "illus ")) %>% 
  mutate(taxon = str_replace(taxon, "lister{1}", "lister ")) %>% 
  mutate(taxon = str_replace(taxon, "asma{1}", "asma ")) %>% 
  mutate(taxon = str_replace(taxon, "ium{1}", "ium ")) %>% 
  mutate(taxon = str_replace(taxon, "ilus{1}", "ilus ")) %>% 
  mutate(taxon = str_replace(taxon, "ldia{1}", "ldia ")) %>% 
  mutate(taxon = str_replace(taxon, "ncus{1}", "ncus ")) %>% 
  mutate(taxon = str_replace(taxon, "culum{1}", "culum ")) %>% 
  mutate(taxon = str_replace(taxon, "bacteru{1}", "bacter u")) %>% 
  mutate(taxon = str_replace(taxon, "bacterh{1}", "bacter h")) %>% 
  mutate(taxon = str_replace(taxon, "myces{1}", "myces ")) %>% 
  mutate(taxon = str_replace(taxon, "itis{1}", "itis ")) %>% 
  mutate(taxon = str_replace(taxon, "genogroup{1}", "genogroup ")) %>% 
  # several require special attention
  mutate(taxon = str_replace(taxon, "Veillonella ceae", "Veillonellaceae")) %>% 
  mutate(taxon = str_replace(taxon, "Clostridiales. Family_XI__Incerta",
                             "Clostridiales Family XI Incerta")) %>% 
  mutate(taxon = str_replace(taxon, "ClostridialesFamilyXI. IncertaeSe",
                             "Clostridiales Family XI Incertae Sedis")) %>% 
  # make first word (family?) single letter (unless there is only one word)
  mutate(first = str_split_fixed(taxon, " ", n = 2)[, 1], 
         second = str_trim(str_split_fixed(taxon, " ", n = 2)[, 2], side = "both")) %>% 
  mutate(first = ifelse(nchar(second) > 0 & !grepl("Family", second), paste0(str_sub(first, 1, 1), "."), first),
         second = ifelse(nchar(second) == 0 & !grepl("bacteria", first) & !grepl("eae", first), "spp.", second)) %>% 
  mutate(abbr_taxon = paste(first, second, sep = " ")) %>% 
  rename(taxon_name = taxon) %>% 
  mutate(taxon = 1:q, fct_taxon = factor(abbr_taxon))

est_table_with_names <- dplyr::left_join(est_table, taxon_nms, by = "taxon")
## plot estimated mu/qPCR for all women
est_plot <- est_table_with_names %>%
  ggplot(aes(x = log_est, group = taxon, fill = factor(taxon, labels = unique(abbr_taxon), ordered = FALSE))) +
  geom_histogram(binwidth = .025) +
  scale_fill_manual(values = rep(cols, ceiling(q / length(cols)))) +
  labs(x = expression(log[10](hat(mu)))) +
  labs(y = "Frequency") +
  labs(fill = "Taxon") +
  ggtitle("Estimated concentrations: each taxon, all women") +
  xlim(c(0, 8)) +
  ylim(c(0, 100)) +
  theme(text = element_text(size = 20), legend.text = element_text(face = "italic")) +
  geom_vline(xintercept = 1, col = "red", linetype = "dashed")
ggsave(filename = paste0(plots_dir, "/est_concentrations_n_", 
                         args$sample_num, sigma_txt, "_", args$adjust, ".png"),
       plot = est_plot,
       width = 75, height = 15, units = "cm")

## -----------------------------------------
## heatmap of all taxa for a subset of women
## -----------------------------------------
samp_women <- ids
## get the subset
set.seed(123456)
qpcr_taxa <- 1:13
other_taxa <- sample(14:127, 7, replace = FALSE)
sub_taxa <- c(qpcr_taxa, other_taxa)
table_for_heatmap <- est_table_with_names %>%
  mutate(fct_taxon_2 = forcats::fct_reorder(fct_taxon, taxon)) %>% 
  filter(subj_id %in% samp_women, taxon %in% sub_taxa) 
  
table_for_heatmap$subj_id <- rep(1:length(samp_women), length(unique(table_for_heatmap$taxon)))
## plot a heatmap of log mu for these women, all taxa
concentration_heatmap <- table_for_heatmap %>%
  mutate(normalized_log_est = log_est/max(table_for_heatmap$log_est)) %>%
  ggplot(aes(y = subj_id, x = fct_taxon, fill = normalized_log_est)) +
  geom_tile(color = "white") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(y = "Subject ID") +
  labs(x = "") +
  labs(fill = expression(paste("Normalized ", log[10](hat(mu)), sep = ""))) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 310, hjust = 0, face = "italic", size = 20),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) # top, right, bottom, left

## ----------------------------------------------------------
## plot recreating Fig 1 but with estimated conc. vs obs qPCR
## ----------------------------------------------------------
## qpcr proportions; need to reorder by samp to match up with estimates
compositional_qpcr_init <- as.data.frame(qpcr_plots_mat) %>%
  mutate(row_sum = rowSums(.)) %>% 
  mutate_at(vars(-row_sum), list(~`/`(., row_sum))) %>% 
  select(-row_sum)
compositional_qpcr <- compositional_qpcr_init
compositional_qpcr$id <- 1:dim(compositional_qpcr)[1]
compositional_qpcr_long <- compositional_qpcr %>% 
  gather(key = taxon, value = proportion, -id)

## concentration proportions
obs_taxa_ests <- extract_posterior_summaries(ves_0$mod, ve_samps, taxa_of_interest = 1:q_obs, mult_num = 1, level = 0.95)
colnames(obs_taxa_ests$estimates) <- colnames(qpcr_plots_mat)
compositional_conc <- as.data.frame(obs_taxa_ests$estimates) %>% 
  mutate(row_sum = rowSums(.)) %>% 
  mutate_at(vars(-row_sum), list(~`/`(., row_sum))) %>% 
  select(-row_sum)
compositional_conc$id <- 1:dim(compositional_conc)[1]
compositional_conc_long <- compositional_conc %>% 
  gather(key = taxon, value = proportion, -id)

obs_taxa_ests_no_ve <- extract_posterior_summaries(no_ves_0$mod, no_ve_samps, taxa_of_interest = 1:q_obs, mult_num = 1, level = 0.95)
colnames(obs_taxa_ests_no_ve$estimates) <- colnames(qpcr_plots_mat)
compositional_conc_no_ve <- as.data.frame(obs_taxa_ests_no_ve$estimates) %>% 
  mutate(row_sum = rowSums(.)) %>% 
  mutate_at(vars(-row_sum), list(~`/`(., row_sum))) %>% 
  select(-row_sum)
compositional_conc_no_ve$id <- 1:dim(compositional_conc_no_ve)[1]
compositional_conc_long_no_ve <- compositional_conc_no_ve %>% 
  gather(key = taxon, value = proportion, -id)

obs_taxa_ests_naive <- naive_0$mod[, 1:q_obs]
colnames(obs_taxa_ests_naive) <- colnames(qpcr_plots_mat)
compositional_conc_naive <- as.data.frame(obs_taxa_ests_naive) %>% 
  mutate(row_sum = rowSums(.)) %>% 
  mutate_at(vars(-row_sum), list(~`/`(., row_sum))) %>% 
  select(-row_sum)
compositional_conc_naive$id <- 1:dim(compositional_conc_naive)[1]
compositional_conc_long_naive <- compositional_conc_naive %>% 
  gather(key = taxon, value = proportion, -id)

## create dataset for plotting
compositional_data <- dplyr::left_join(compositional_conc_long, compositional_qpcr_long, by = c("id", "taxon"), suffix = c(".est_conc", ".qPCR")) %>% 
  select(id, taxon, proportion.qPCR, proportion.est_conc) %>% 
  mutate(taxon_id = rep(1:q_obs, each = length(unique(id)))) %>% 
  left_join(taxon_nms %>% 
              rename(taxon_id = taxon) %>%  
              select(taxon_id, abbr_taxon), by = "taxon_id")

compositional_data_no_ve <- dplyr::left_join(compositional_conc_long_no_ve, compositional_qpcr_long, by = c("id", "taxon"), suffix = c(".est_conc", ".qPCR")) %>% 
  select(id, taxon, proportion.qPCR, proportion.est_conc) %>% 
  mutate(taxon_id = rep(1:q_obs, each = length(unique(id)))) %>% 
  left_join(taxon_nms %>% 
              rename(taxon_id = taxon) %>%  
              select(taxon_id, abbr_taxon), by = "taxon_id")

compositional_data_naive <- dplyr::left_join(compositional_conc_long_naive, compositional_qpcr_long, by = c("id", "taxon"), suffix = c(".est_conc", ".qPCR")) %>% 
  select(id, taxon, proportion.qPCR, proportion.est_conc) %>% 
  mutate(taxon_id = rep(1:q_obs, each = length(unique(id)))) %>% 
  left_join(taxon_nms %>% 
              rename(taxon_id = taxon) %>%  
              select(taxon_id, abbr_taxon), by = "taxon_id")


## create the plot
compositional_data %$% abbr_taxon %>% unique # print the taxon names currently
compositional_plot <- compositional_data %>% 
  ggplot(aes(x = proportion.qPCR, y = proportion.est_conc, color = abbr_taxon)) +
  geom_point() + 
  scale_color_manual(values = cols) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Proportion within subcomposition: qPCR", y = "Proportion within subcomposition:\n estimated concentration", color = "Taxon") +
  theme(text = element_text(size = 15), legend.text=element_text(size = 15))

compositional_plot_no_ve <- compositional_data_no_ve %>% 
  ggplot(aes(x = proportion.qPCR, y = proportion.est_conc, color = abbr_taxon)) +
  geom_point() + 
  scale_color_manual(values = cols) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Proportion within subcomposition: qPCR", y = "Proportion within subcomposition:\n estimated concentration", color = "Taxon") +
  theme(text = element_text(size = 15), legend.text=element_text(size = 15))

compositional_plot_naive <- compositional_data_naive %>% 
  ggplot(aes(x = proportion.qPCR, y = proportion.est_conc, color = abbr_taxon)) +
  geom_point() + 
  scale_color_manual(values = cols) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Proportion within subcomposition: qPCR", y = "Proportion within subcomposition:\n estimated concentration", color = "Taxon") +
  theme(text = element_text(size = 15), legend.text=element_text(size = 15))

compositional_plot
compositional_plot_no_ve
compositional_plot_naive
## ----------------------------------------------------------
## save off the plot
## ----------------------------------------------------------
# width was 40 for all 127 taxa... for 20, bring it down
# full_fig_width <- 40
full_fig_width <- 20
ggsave(paste0(plots_dir, "/concentration_full_fig_n_", args$sample_num, 
              sigma_txt, "_", args$adjust, ".png"),
       plot = plot_grid(concentration_heatmap, 
                        compositional_plot, labels = c("A", "B"),
                        rel_widths = c(2, 1)),
       width = full_fig_width, height = 8)
ggsave(paste0(plots_dir, "/compositional_plot_ve_n_", args$sample_num, 
              sigma_txt, "_", args$adjust, ".png"),
       plot = compositional_plot,
       width = 9, height = 8)
ggsave(paste0(plots_dir, "/compositional_plot_no_ve_n_", args$sample_num, 
              sigma_txt, "_", args$adjust, ".png"),
       plot = compositional_plot_no_ve,
       width = 9, height = 8)
ggsave(paste0(plots_dir, "/compositional_plot_naive_n_", args$sample_num, 
              sigma_txt, "_", args$adjust, ".png"),
       plot = compositional_plot_naive,
       width = 9, height = 8)

# on the log-log scale
ggsave(paste0(plots_dir, "/concentration_full_fig_n_", args$sample_num, 
              sigma_txt, "_", args$adjust, "_loglog.png"),
       plot = plot_grid(concentration_heatmap, 
                        compositional_plot + scale_x_log10() +
                          scale_y_log10(), 
                        labels = c("A", "B"),
                        rel_widths = c(2, 1)),
       width = full_fig_width, height = 8)
ggsave(paste0(plots_dir, "/compositional_plot_ve_n_", args$sample_num, 
              sigma_txt, "_", args$adjust, "_loglog.png"),
       plot = compositional_plot + scale_x_log10() +
         scale_y_log10(),
       width = 9, height = 8)
ggsave(paste0(plots_dir, "/compositional_plot_no_ve_n_", args$sample_num, 
              sigma_txt, "_", args$adjust, "_loglog.png"),
       plot = compositional_plot_no_ve + scale_x_log10() +
         scale_y_log10(),
       width = 9, height = 8)
ggsave(paste0(plots_dir, "/compositional_plot_naive_n_", args$sample_num, 
              sigma_txt, "_", args$adjust, "_loglog.png"),
       plot = compositional_plot_naive + scale_x_log10() +
         scale_y_log10(),
       width = 9, height = 8)

## ----------------------------------------------------------
## plot of beta_1's (only for adjust); also save all as a table
## ----------------------------------------------------------
# pick the top 10 by absolute value
num_beta1s <- 10
main_font_size <- 30
main_font_size_lab <- 9
title_font_size <- 26
str_width <- 20
if (args$adjust) {
  beta_1s <- ves_0$mod[grepl("beta_1", rownames(ves_0$mod)), ]
  beta_1_tib <- as_tibble(beta_1s) %>% 
    mutate(taxon = taxon_nms$abbr_taxon, taxon_num = 1:127)
  readr::write_csv(beta_1_tib, path = paste0(tables_dir, analysis_name, "/est_beta1s_n_", args$sample_num, "_", args$adjust, ".csv"))
  abs_est_beta_1 <- abs(pull(beta_1_tib, mean))
  top_beta_1s <- beta_1_tib %>% 
    filter(taxon_num %in% order(abs_est_beta_1, decreasing = TRUE)[1:num_beta1s])
  
  # plot them
  beta_1_plot <- top_beta_1s %>%  # make levels correspond to taxon number
    mutate(taxon = ifelse(taxon == "G. sp", "Gardnerella sp.", taxon)) %>% 
    ggplot(aes(y = mean, x = forcats::fct_reorder(factor(taxon, ordered = TRUE), 
                                                  mean, .desc = FALSE))) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
    geom_point() +
    ggtitle(bquote("Estimated"~beta[1]~"with 95% credible intervals")) +
    ylab(bquote("Estimated"~beta[1])) + 
    xlab("Taxon") +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    coord_flip() +
    theme(axis.text.y = element_text(hjust = 1),
          legend.position = c(0.5, 0.7),
          text = element_text(size = main_font_size),
          axis.text = element_text(size = main_font_size),
          plot.title = element_text(size = main_font_size - 8),
          plot.margin = unit(c(0, 1, 0, 0), "cm"))
  ggsave(paste0(plots_dir, "/est_beta1_n_", args$sample_num, "_", args$adjust, ".png"),
         plot = beta_1_plot,
         width = 9, height = 8)
  # check W for the ones with top beta_1s
  threshold <- 5
  br16_mat[samp_women, pull(top_beta_1s, taxon_num)] %>% 
    as_tibble() %>% 
    pivot_longer(cols = tidyselect::everything(), names_to = "taxon", values_to = "w") %>%
    mutate(subj_id = rep(samp_women, each = num_beta1s),
           above_thresh = w > threshold) %>% 
    group_by(taxon) %>% 
    summarize(thresh_perc = mean(above_thresh))
}

# ---------------------------------------------------------------
# Histogram of estimated taxon means / sum of estimated taxon means
# for L crispatus and L iners (for review response)
# ---------------------------------------------------------------
l_iners_crisp_compositional_data <- compositional_data %>% 
  filter(abbr_taxon %in% c("L. iners", "L. crispatus"))
l_iners_crisp_subcomposition <- l_iners_crisp_compositional_data %>% 
  ggplot(aes(x = proportion.est_conc, fill = abbr_taxon)) +
  geom_histogram(binwidth = 25 * 1e-3) +
  scale_color_manual(values = cols) +
  labs(x = "Proportion within subcomposition:\n estimated concentration", 
       fill = "Taxon") +
  theme(text = element_text(size = 15), legend.text=element_text(size = 15))
ggsave(filename = paste0(plots_dir, "/l-iners-crisp_subcomposition_n_", 
                         args$sample_num, "_", args$adjust, ".png"),
       l_iners_crisp_subcomposition,
       width = 9, height = 8)
## ----------------------------------------------------------
## plot of mean interval width (unadjusted, adjusted, different sigma_e?)
## ----------------------------------------------------------

## read in both unadjusted and adjusted
ves_0 <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_est_ve_q_", q, "_q_obs_", q_obs, "_sample_", args$sample_num, "_loo_999_adjust_", 0, ".rds", sep = ""))
ves_1 <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_est_ve_q_", q, "_q_obs_", q_obs, "_sample_", args$sample_num, "_loo_999_adjust_", 1, ".rds", sep = ""))
ves_sens <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_est_ve_q_", q, "_q_obs_", q_obs, "_sample_", args$sample_num, "_loo_999_as_2_ks_3_adjust_", 0, ".rds", sep = ""))
ve_samps <- ves_0$samps
ve_samps_1 <- ves_1$samps
ve_samps_sens <- ves_sens$samps
all_taxa_ests <- extract_posterior_summaries(ves_0$mod, ve_samps, taxa_of_interest = 1:q, mult_num = 1, level = 0.95)
all_taxa_ests_1 <- extract_posterior_summaries(ves_1$mod, ve_samps_1, taxa_of_interest = 1:q, mult_num = 1, level = 0.95)
all_taxa_ests_sens <- extract_posterior_summaries(ves_sens$mod, ve_samps_sens, taxa_of_interest = 1:q, mult_num = 1, level = 0.95)

## check mean credible interval width
pred_width_0 <- colMeans(apply(all_taxa_ests$pred_intervals, 3, function(x) x[, 2] - x[, 1]))
pred_width_1 <- colMeans(apply(all_taxa_ests_1$pred_intervals, 3, function(x) x[, 2] - x[, 1]))
pred_width_sens <- colMeans(apply(all_taxa_ests_sens$pred_intervals, 3, function(x) x[, 2] - x[, 1]))

cred_width_0 <- colMeans(apply(all_taxa_ests$cred_intervals, 3, function(x) x[, 2] - x[, 1]))
cred_width_1 <- colMeans(apply(all_taxa_ests_1$cred_intervals, 3, function(x) x[, 2] - x[, 1]))
cred_width_sens <- colMeans(apply(all_taxa_ests_sens$cred_intervals, 3, function(x) x[, 2] - x[, 1]))

width_tib <- tibble(taxon = rep(1:q, 6), 
                    model = rep(rep(c("adjusted", "unadjusted", "sensitivity"), 
                                    each = q), 2), 
       interval_type = rep(c("credible", "prediction"), each = 3 * q),
       width = c(cred_width_1, cred_width_0, cred_width_sens,
                 pred_width_1, pred_width_0, pred_width_sens))
pal <- 3
width_plt_all <- width_tib %>% 
  ggplot(aes(x = taxon, y = width, group = factor(paste0(model, "_", interval_type)),
             color = model, shape = interval_type)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = pal) +
  labs(x = "Taxon", y = "Mean interval width", 
       color = "Model", shape = "Type of interval") +
  ggtitle("Mean interval width")

width_plt_subset <- width_tib %>% 
  filter(width < 1e5) %>% 
  ggplot(aes(x = taxon, y = width, group = factor(paste0(model, "_", interval_type)),
             color = model, shape = interval_type)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = pal) +
  labs(x = "Taxon", y = "Mean interval width", 
       color = "Model", shape = "Type of interval") +
  ggtitle("Mean interval width (subset)")

ggsave(filename = paste0(plots_dir, "/adjusted-unadjusted_widths.png"),
       plot = plot_grid(width_plt_all, width_plt_subset),
       width = 40, height = 20, units = "cm")

# check distribution of qpcrs
formatC(apply(qpcr, 2, summary), format = "e", digits = 2)
formatC(apply(qpcr, 2, var), format = "e", digits = 2)
formatC(apply(qpcr, 2, mean), format = "e", digits = 2)
# check for br16S
formatC(apply(br[, 1:13], 2, summary), format = "e", digits = 2)
formatC(apply(br[, 1:13], 2, var), format = "e", digits = 2)
