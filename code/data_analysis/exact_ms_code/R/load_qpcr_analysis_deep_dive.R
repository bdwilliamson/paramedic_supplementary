######################################################################################
##
## FILE: load_qpcr_analysis_deep_dive.R
##
## CREATED: 18 January 2019 by Brian Williamson
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
library("dplyr")
source(paste0(code_dir, "analyze_data/qpcr_analysis_helper_functions.R"))
source(paste0(code_dir, "naive_qpcr_estimator.R"))
source(paste0(code_dir, "analyze_data/get_most_abundant_taxa.R"))
source(paste0(code_dir, "gen_bootstrap_interval.R"))
source(paste0(code_dir, "gen_wald_interval.R"))
source(paste0(code_dir, "gen_quantile_interval.R"))
source(paste0(code_dir, "extract_posterior_summaries.R"))

## set up arguments
q <- 8
q_obs <- 7
samp_num <- 10
div_num <- 1000
analysis_name <- "data_analysis/sample_10_most_abundant_8/iter_10K"

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

set.seed(4747)
## break up the sample into even chunks
# folds_init <- rep(seq_len(args$num_folds), length = dim(br16_mat)[1])
# folds <- sample(folds_init)
## sample women
samp <- sample(1:dim(br16_mat)[1], samp_num)  

## -----------------------------------------------------------------------------------
## load in the estimators (all observed qPCRs)
## -----------------------------------------------------------------------------------

## read in the datasets
naive <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_est_naive_q_", q, "_q_obs_", q_obs, "_sample_", samp_num, "_loo_999.rds", sep = ""))
no_ves <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_est_no_ve_q_", q, "_q_obs_", q_obs, "_sample_", samp_num, "_loo_999.rds", sep = ""))
ves <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_est_ve_q_", q, "_q_obs_", q_obs, "_sample_", samp_num, "_loo_999.rds", sep = ""))

## -----------------------------------------------------------------------------------
## extract the estimates of qPCR and prediction intervals for qPCR
## -----------------------------------------------------------------------------------
## taxa of interest
taxa_of_interest <- (q_obs + 1):q
## naive estimator
naive_qpcr <- naive$mod[, taxa_of_interest]*div_num
## proposed estimator without and with varying efficiency
no_ve_qpcr_lst <- extract_posterior_summaries(no_ves$mod, q = q, taxa_of_interest, mult_num = div_num, level = 0.95, interval_type = "wald")
ve_qpcr_lst <- extract_posterior_summaries(ves$mod, q = q, taxa_of_interest, mult_num = div_num, level = 0.95, interval_type = "wald")

## check convergence diagnostics
summary(no_ves$mod)
summary(ves$mod)

## check for mu only

## check for efficiency only
summary(ves$mod[grepl("e", rownames(ves$mod)) & !grepl("beta", rownames(ves$mod)), ])
ves$mod[grepl("e", rownames(ves$mod)) & !grepl("beta", rownames(ves$mod)), ]
ves$mod[grepl("sigma", rownames(ves$mod)), ]

## -----------------------------------------------------------------------------------
## summarize, for each taxon
## -----------------------------------------------------------------------------------
if (is.matrix(naive_qpcr)) {
  naive_taxon_means <- colMeans(naive_qpcr)  
} else {
  naive_taxon_means <- mean(naive_qpcr)
}
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
## create output tables
## -----------------------------------------------------------------------------------
efficiency_table <- ve_qpcr_lst$est_efficiency
hist(ve_qpcr_lst$est_efficiency, breaks = 20, xlab = "Estimated efficiency", main = "Histogram of estimated efficiencies")
print(xtable::xtable(matrix(efficiency_table, nrow = 1), caption = "Estimated efficiencies for each taxon.", label = "tab:data_est_efficiencies"),
      file = paste0("tables/", analysis_name, "/est_efficiencies.tex"))

ests_sub <- ve_qpcr_lst$estimates
cis_sub <- plyr::adply(ve_qpcr_lst$cred_intervals, 3)
pis_sub <- plyr::adply(ve_qpcr_lst$pred_intervals, 3)
ests_sub_table <- as.vector(t(apply(ests_sub, 1, function(x) as.character(round(x, 3)))))
cis_sub_table <- as.vector(interval_long_to_wide(cis_sub, digits = 3, n_taxa = 1, n_women = 10))
pis_sub_table <- as.vector(interval_long_to_wide(pis_sub, digits = 3, n_taxa = 1, n_women = 10))
est_table_sub <- cbind(ests_sub_table, cis_sub_table, pis_sub_table)
colnames(est_table_sub) <- c("Point Est.", "Cred. Int.", "Pred. Int.")
print(xtable::xtable(est_table_sub, caption = "Point estimates, posterior credible intervals for $\\mu$, and posterior prediction intervals for $\\V$ for the first five taxa with unobserved qPCR and five randomly sampled women.",
                     label = "tab:data_concentration_ests_sub"), file = paste0("tables/", analysis_name, "/concentration_ests_sub.tex"),
      sanitize.text.function = function(x) {x})

## print everything else out for supplement
ests <- ve_qpcr_lst$estimates
cis <- plyr::adply(ve_qpcr_lst$cred_intervals, 3)
pis <- plyr::adply(ve_qpcr_lst$pred_intervals, 3)
set.seed(4747)
pis_boot <- gen_bootstrap_interval(ve_qpcr_lst$estimates, b = 100000, alpha = 0.05, na.rm = TRUE)
pis_boot
ests_table <- as.vector(t(apply(ests, 1, function(x) as.character(round(x, 3)))))
cis_table <- as.vector(interval_long_to_wide(cis, digits = 3, n_taxa = q - q_obs, n_women = samp_num))
pis_table <- as.vector(interval_long_to_wide(pis, digits = 3, n_taxa = q - q_obs, n_women = samp_num))
est_table <- cbind(ests_table, cis_table, pis_table)
taxon_woman_combos <- expand.grid(1:samp_num, tail(1:q, q - q_obs))
rownames(est_table) <- apply(taxon_woman_combos, 1, function(x) paste0("id_", x[1], "_taxon_", x[2]))
colnames(est_table) <- c("Point Est.", "Cred. Int.", "Pred. Int.")
print(xtable::xtable(est_table, caption = "Point estimates, posterior credible intervals for $\\mu$, and posterior prediction intervals for $\\V$ for all taxa with unobserved qPCR and all women.",
                     label = "tab:data_concentration_ests"), include.rownames = TRUE, file = paste0("tables/", analysis_name, "/concentration_ests.tex"),
      sanitize.text.function = function(x) {x})
write.csv(est_table, file = paste0("tables/", analysis_name, "/concentration_ests.csv"))

## do the same for leave-one-out analysis 1

## --------------------------------------------------------------------------------------------------------
##
## Leave-one-out analysis:
## Compute point estimates, Wald-type PIs, model-based PIs, credible intervals (for range of levels)
## Compute coverage, MSE, MSPE for the left-out one
##
## --------------------------------------------------------------------------------------------------------

leave_one_out_performance <- vector("list", length = 7)
for (i in 1:7) {
  ## read in the datasets
  naive_i <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_est_naive_q_", q, "_q_obs_", q_obs - 1, "_sample_", samp_num, "_loo_", i, ".rds", sep = ""))
  no_ves_i <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_est_no_ve_q_", q, "_q_obs_", q_obs - 1, "_sample_", samp_num, "_loo_", i, ".rds", sep = ""))
  ves_i <- readRDS(paste0(results_dir, analysis_name, "/qpcr_data_analysis_est_ve_q_", q, "_q_obs_", q_obs - 1, "_sample_", samp_num, "_loo_", i, ".rds", sep = ""))
  
  ## -----------------------------------------------------------------------------------
  ## extract the estimates of qPCR and prediction intervals for qPCR
  ## -----------------------------------------------------------------------------------
  ## taxa of interest
  taxa_of_interest <- 7 # always the first one outside of q_obs for the leave-one-out
  ## naive estimator
  naive_qpcr_i <- naive_i$mod[, taxa_of_interest]*div_num
  ## proposed estimator without and with varying efficiency
  set.seed(4747)
  no_ve_qpcr_lst_i <- extract_posterior_summaries(no_ves_i$mod, q = q, taxa_of_interest, mult_num = div_num, level = 0.95, na.rm = TRUE)
  set.seed(4747)
  ve_qpcr_lst_i <- extract_posterior_summaries(ves_i$mod, q = q, taxa_of_interest, mult_num = div_num, level = 0.95, na.rm = TRUE)
  
  ## get mspe
  ests <- cbind(naive_qpcr_i, no_ve_qpcr_lst_i$estimates, ve_qpcr_lst_i$estimates)
  spe_mat <- (ests - qpcr_mat[samp, i]*div_num)^2
  mspe <- colMeans(spe_mat)
  
  ## get prediction intervals, coverage
  no_ve_pis_boot <- no_ve_qpcr_lst_i$pred_intervals
  ve_pis_boot <- ve_qpcr_lst_i$pred_intervals
  no_ve_pis_wald <- gen_wald_interval(no_ve_qpcr_lst_i$estimates, no_ve_qpcr_lst_i$sd[, taxa_of_interest], alpha = 0.05)
  ve_pis_wald <- gen_wald_interval(ve_qpcr_lst_i$estimates, ve_qpcr_lst_i$sd[, taxa_of_interest], alpha = 0.05)
  no_ve_pis_quantile <- gen_quantile_interval(extract(no_ves_i$stan_out)$mu[, , taxa_of_interest], alpha = 0.05)
  
  no_ve_boot_cover <- mean(no_ve_pis_boot[,1,1] <= qpcr_mat[samp, i]*div_num & no_ve_pis_boot[,2,1] >= qpcr_mat[samp, i]*div_num)
  ve_boot_cover <- mean(ve_pis_boot[,1,1] <= qpcr_mat[samp, i]*div_num & ve_pis_boot[,2,1] >= qpcr_mat[samp, i]*div_num, na.rm = TRUE)
  no_ve_wald_cover <- mean(no_ve_pis_wald[,1] <= qpcr_mat[samp, i]*div_num & no_ve_pis_wald[,2] >= qpcr_mat[samp, i]*div_num)
  ve_wald_cover <- mean(ve_pis_wald[,1] <= qpcr_mat[samp, i]*div_num & ve_pis_wald[,2] >= qpcr_mat[samp, i]*div_num)
  ## get credible intervals
  no_ve_cis_95 <- no_ve_qpcr_lst_i$cred_intervals
  ve_cis_95 <- no_ve_qpcr_lst_i$cred_intervals
  
  no_ve_cis_75 <- extract_posterior_summaries(summary(no_ves_i$stan_out, probs = c(0.125, 0.875))$summary, q = q, taxa_of_interest, mult_num = div_num, level = 0.75, na.rm = TRUE)$cred_intervals
  ve_cis_75 <- extract_posterior_summaries(summary(ves_i$stan_out, probs = c(0.125, 0.875))$summary, q = q, taxa_of_interest, mult_num = div_num, level = 0.75, na.rm = TRUE)$cred_intervals
  
  
  no_ve_cis_50 <- extract_posterior_summaries(summary(no_ves_i$stan_out, probs = c((1-0.5)/2, 1 - (1-0.5)/2))$summary, q = q, taxa_of_interest, mult_num = div_num, level = 0.75, na.rm = TRUE)$cred_intervals
  ve_cis_50 <- extract_posterior_summaries(summary(ves_i$stan_out, probs = c((1-0.5)/2, 1 - (1-0.5)/2))$summary, q = q, taxa_of_interest, mult_num = div_num, level = 0.75, na.rm = TRUE)$cred_intervals
  
  ## save them all
  leave_one_out_performance[[i]] <- list(mspe = mspe, spe_mat = spe_mat, no_ve_pis = cbind(no_ve_pis_boot[,,1], no_ve_pis_wald),
                                    ve_pis = cbind(ve_pis_boot[,,1], ve_pis_wald), no_ve_cis = cbind(no_ve_cis_95[,,1], no_ve_cis_75[,,1], no_ve_cis_50[,,1]),
                                    ve_cis = cbind(ve_cis_95[,,1], ve_cis_75[,,1], ve_cis_50[,,1]),
                                    no_ve_boot_cover = no_ve_boot_cover, ve_boot_cover = ve_boot_cover,
                                    no_ve_wald_cover = no_ve_wald_cover, ve_wald_cover = ve_wald_cover,
                                    obs_qpcr = qpcr_mat[samp, i]*div_num)
}