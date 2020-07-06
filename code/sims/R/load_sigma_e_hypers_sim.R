# load_sigma_e_hypers_sim.R
# load results from the sigma e hyperparameters simulation, create plots

# ------------------------------------------------------------------------------
# setup
# ------------------------------------------------------------------------------
## load required functions and packages
library("methods")
library("argparse")
library("Cairo")
library("tidyr")
library("dplyr")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
source("code/R/load_qpcr_sim_helpers.R") # provides get_summaries
source("code/R/naive_qpcr_estimator.R")
source("code/R/gen_wald_interval.R")
source("code/R/gen_quantile_interval.R")
source("code/R/plot_qpcr_sim_helpers.R")
source("code/R/analyze_data/get_most_abundant_taxa.R")
source("code/R/gen_naive_interval.R")

# ---------------------------------
# grab command-line args
# ---------------------------------
parser <- ArgumentParser()
parser$add_argument("--N", type = "double", default = 100, help = "sample size")
parser$add_argument("--q", type = "double", default = 40, help = "how many taxa do we have?")
parser$add_argument("--q_obs", type = "double", default = 7, help = "how many taxa do we have qPCR measured on?")
parser$add_argument("--n-chains", type = "double", default = 1, help = "number of chains for MCMC")
parser$add_argument("--B", type = "double", default = 50, help = "total number of MC reps per q, q^obs")
parser$add_argument("--adapt-delta", type = "double", default = 0.85, help = "adapt_delta, for Stan fitting")
parser$add_argument("--max-treedepth", type = "integer", default = 18, help = "max_treedepth, for Stan fitting")
parser$add_argument("--read-data", type = "double", default = 1, help = "whether or not to read in raw output, or instead do summaries")
args <- parser$parse_args()

# -----------------------------
# a bit more setup
# -----------------------------
read_func <- function(x) tryCatch(readRDS(x), error = function(e) NA)
plots_dir <- paste0("plots/sigma_e_hypers/n_", args$N, "/")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# ------------------------------------------------------------------
# create a list of output for easy summaries and get those summaries
# ------------------------------------------------------------------
all_e_combos <- expand.grid(alpha_sigma = 1:3, kappa_sigma = c(0.5, 1))
e_nms <- paste0("sigma-e-", apply(all_e_combos, 1, function(x) paste(x, collapse = "-")))
get_e_labs <- function(x) {
       gsub("-", "/", gsub("sigma-e-", "", x))
}
e_labels <- get_e_labs(e_nms)
results_dir <- paste0("results/sigma_e_hypers/")
if (args$read_data) {
  output_performances <- vector(mode = "list", length = length(e_nms))
  
  output_performances_avg_over_taxa_n <- output_performances_avg_over_taxa_n_mn <- output_performances_avg_over_taxa_n_var <- output_performances_avg_over_taxa_n_nofilter <- vector(mode = "list", length = length(e_nms))
  output_performances_avg_over_n <- output_performances_avg_over_n_mn <- output_performances_avg_over_n_var <- output_performances_avg_over_n_nofilter <- vector(mode = "list", length = length(e_nms))
  e_output_performances <- e_output_performances_avg_over_n <- e_output_performances_avg_over_n_taxa <- vector(mode = "list", length = length(e_nms))
  for (i in 1:length(e_nms)) {
    dir_mat <- expand.grid(job = 1:args$B, dir = paste0(results_dir, e_nms[i], "/n_", args$N, "/q_", args$q, "/q_obs_", args$q_obs, "/"))
    mod_nms_lst <- as.list(paste0(dir_mat$dir, "run_paramedic_mod_jobid_", dir_mat$job, "_ad_", args$adapt_delta, "0000_mt_", args$max_treedepth, ".rds"))
    data_nms_lst <- as.list(paste0(dir_mat$dir, "run_paramedic_data_jobid_", dir_mat$job, "_ad_", args$adapt_delta, "0000_mt_", args$max_treedepth, ".rds"))
    samp_nms_lst <- as.list(paste0(dir_mat$dir, "run_paramedic_samps_jobid_", dir_mat$job, "_ad_", args$adapt_delta, "0000_mt_", args$max_treedepth, ".rds"))
    
    mod_lst <- lapply(mod_nms_lst, read_func)
    data_lst <- lapply(data_nms_lst, read_func)
    mc_ids <- 1:args$B
    
    ## (1) pair model summaries with relevant data, for all taxa, all q_obs; replace NA with samps_lst if you want
    summary_df <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type, mc_id) get_summaries(mod_summ = w, data = x, samps = y, mc_id = mc_id, indx = z, type), w = mod_lst, x = data_lst, y = NA, mc_id = mc_ids, MoreArgs = list(z = 1:args$q, type = "ve"), SIMPLIFY = FALSE))
    summary_df <- as_tibble(summary_df) %>% 
      filter(!is.na(summary_df$q)) %>% 
      mutate(alpha_sigma = all_e_combos$alpha_sigma[i], kappa_sigma = all_e_combos$kappa_sigma[i])
    
    e_summary_df <- do.call(rbind.data.frame, mapply(function(w, x, y, mc_id) get_e_summaries(mod_summ = w, data = x, mc_id = mc_id, indx = y), w = mod_lst, x = data_lst, mc_id = mc_ids, MoreArgs = list(y = 1:args$q), SIMPLIFY = FALSE)) %>% 
      filter(!is.na(q)) %>% 
      mutate(alpha_sigma = all_e_combos$alpha_sigma[i], kappa_sigma = all_e_combos$kappa_sigma[i]) %>% 
      as_tibble()

    ## (2) compute performance for each row
    performance_df <- summary_df %>%
      mutate(bias = mu - est, mse = (mu - est) ^ 2, mspe = (qpcr - est) ^ 2,
             cover = cil <= mu & ciu >= mu, 
             wald_cover = wald_cred_cil <= mu & wald_cred_ciu >= mu,
             pred_cover = wald_pred_cil <= qpcr & wald_pred_ciu >= qpcr,
             cred_width = abs(ciu - cil), wald_cred_width = abs(wald_cred_ciu - wald_cred_cil),
             pred_width = abs(wald_pred_ciu - wald_pred_cil),
             include_taxa_in_v_avg = taxon_id > q_obs)
    e_performance_df <- e_summary_df %>% 
      mutate(bias = e - est_e, mse = (e - est_e) ^ 2,
             cover = e_cil <= e & e_ciu >= e, 
             cred_width = abs(e_ciu - e_cil), estimator = "ve", # the next things are all to fool helper functions
             subj_id = NA, mspe = NA,
             wald_cover = NA, pred_cover = NA, wald_cred_width = NA,
             pred_width = NA, include_taxa_in_v_avg = NA)
    
    ## (3) average over MC reps for each taxon
    mc_averaged_performance <- get_average_performance(performance_df %>% filter(include_taxa_in_avg),
                                                       group = "drop_mc", var = c("kappa_sigma", "alpha_sigma"))
    mc_averaged_performance_nofilter <- get_average_performance(performance_df,
                                                                group = "drop_mc", var = c("kappa_sigma", "alpha_sigma"))
    e_mc_averaged_performance <- get_average_performance(e_performance_df,
                                                       group = "drop_mc", var = c("kappa_sigma", "alpha_sigma"))
    ## check how many naives got filtered out since the intervals were NA
    num_na_naive_over_mc <- performance_df %>%
      filter(estimator == "naive") %>% 
      group_by(q, q_obs, subj_id, taxon_id, alpha_sigma, kappa_sigma) %>%
      summarize(perc_na = mean(is.na(cover)))
    num_na_naive_by_taxon_id <- num_na_naive_over_mc %>%
      group_by(q, q_obs, taxon_id, alpha_sigma, kappa_sigma) %>%
      summarize(perc_na = mean(perc_na))
    num_na_naive <- num_na_naive_by_taxon_id %>%
      group_by(q, q_obs, alpha_sigma, kappa_sigma) %>%
      summarize(perc_na = mean(perc_na))
    print(num_na_naive)
    ## (4) set flag for taxa of interest: for V, need taxon_id > q^obs
    mc_averaged_performance <- mc_averaged_performance %>% 
      mutate(include_taxa_in_v_avg = taxon_id > q_obs)
    mc_averaged_performance_nofilter <- mc_averaged_performance_nofilter %>% 
      mutate(include_taxa_in_v_avg = taxon_id > q_obs)
    e_mc_averaged_performance <- e_mc_averaged_performance %>% 
      mutate(include_taxa_in_v_avg = NA) # again, fool the helper functions
    
    ## (5) average over taxa of interest for each q_obs
    performance_across_taxa <- get_performance_filter(mc_averaged_performance,
                                                      group = "taxa", 
                                                      type = c("alpha_sigma", "kappa_sigma"))
    e_performance_across_taxa <- get_performance_filter(e_mc_averaged_performance,
                                                        group = "taxa",
                                                        type = c("alpha_sigma", "kappa_sigma"))
    ## without filter
    performance_across_taxa_nofilter <- get_performance_filter(mc_averaged_performance_nofilter,
                                                               group = "taxa", 
                                                               type = c("alpha_sigma", "kappa_sigma"))
    ## (6) average over n, for each q_obs
    average_over_n <- get_average_performance(performance_across_taxa,
                                              group = "drop_n", 
                                              var = c("alpha_sigma", "kappa_sigma"))
    e_average_over_n <- get_average_performance(e_performance_across_taxa,
                                              group = "drop_n", 
                                              var = c("alpha_sigma", "kappa_sigma"))
    ## no filter
    average_over_n_nofilter <- get_average_performance(performance_across_taxa_nofilter,
                                                       group = "drop_n", 
                                                       var = c("alpha_sigma", "kappa_sigma"))
    ## (5) add rmse, transpose
    performance_matrix <- average_over_n %>%
      mutate(rmse = sqrt(mse),
             rmspe = sqrt(mspe))
    performance_matrix_nofilter <- average_over_n_nofilter %>%
      mutate(rmse = sqrt(mse),
             rmspe = sqrt(mspe))
    e_performance_matrix <- e_average_over_n %>% 
      mutate(rmse = sqrt(mse),
             rmspe = sqrt(mspe))
    
    ## (6) average only over n, taxa
    ## no filter
    average_over_n_taxa_nofilter <- get_performance_filter(performance_df, 
                                                           group = "taxa+subj", 
                                                           type = c("alpha_sigma", "kappa_sigma"))
    e_average_over_n_taxa <- get_performance_filter(e_performance_df, 
                                                           group = "taxa+subj", 
                                                           type = c("alpha_sigma", "kappa_sigma"))
    ## (a) filter by colSums(W) > 0
    average_over_n_taxa <- get_performance_filter(performance_df %>% filter(include_taxa_in_avg),
                                                  group = "taxa+subj", 
                                                  type = c("alpha_sigma", "kappa_sigma"))
    ## (b) filter by colMeans(W) > 0.5
    average_over_n_taxa_mn <- get_performance_filter(performance_df %>% filter(include_taxa_in_avg_w_mn),
                                                     group = "taxa+subj", 
                                                     type = c("alpha_sigma", "kappa_sigma"))
    ## (c) filter by colVars(W) > 1
    average_over_n_taxa_var <- get_performance_filter(performance_df %>% filter(include_taxa_in_avg_var),
                                                      group = "taxa+subj", 
                                                      type = c("alpha_sigma", "kappa_sigma"))
    
    ## (7) average only over n
    ## (a) filter by colSums(W) > 0
    average_over_n <- get_performance_filter(performance_df %>% filter(include_taxa_in_avg),
                                             group = "subj", 
                                             type = c("alpha_sigma", "kappa_sigma"))
    e_average_over_n <- get_performance_filter(e_performance_df,
                                             group = "subj", 
                                             type = c("alpha_sigma", "kappa_sigma"))
    ## (b) filter by colMeans(W) > 0
    average_over_n_mn <- get_performance_filter(performance_df %>% filter(include_taxa_in_avg_w_mn),
                                                group = "subj", 
                                                type = c("alpha_sigma", "kappa_sigma"))
    ## (a) filter by colVars(W) > 0
    average_over_n_var <- get_performance_filter(performance_df %>% filter(include_taxa_in_avg_var),
                                                 group = "subj", 
                                                 type = c("alpha_sigma", "kappa_sigma"))
    ## no filter
    average_over_n_nofilter <- get_performance_filter(performance_df,
                                                      group = "subj", 
                                                      type = c("alpha_sigma", "kappa_sigma"))
    
    output_performances[[i]] <- performance_matrix
    e_output_performances[[i]] <- e_performance_matrix
    output_performances_avg_over_taxa_n[[i]] <- average_over_n_taxa
    output_performances_avg_over_n[[i]] <- average_over_n
    e_output_performances_avg_over_n[[i]] <- e_average_over_n
    e_output_performances_avg_over_n_taxa[[i]] <- e_average_over_n_taxa
    output_performances_avg_over_taxa_n_mn[[i]] <- average_over_n_taxa_mn
    output_performances_avg_over_n_mn[[i]] <- average_over_n_mn
    output_performances_avg_over_taxa_n_var[[i]] <- average_over_n_taxa_var
    output_performances_avg_over_n_var[[i]] <- average_over_n_var
    output_performances_avg_over_taxa_n_nofilter[[i]] <- average_over_n_taxa_nofilter
    output_performances_avg_over_n_nofilter[[i]] <- average_over_n_nofilter
  }
  if (!dir.exists(paste0(results_dir, "n_", args$N))) {
    dir.create(paste0(results_dir, "n_", args$N))
  }
  saveRDS(output_performances, paste0(results_dir, "n_", args$N, "/output_performances.rds"))
  saveRDS(output_performances_avg_over_taxa_n, paste0(results_dir, "n_", args$N, "/output_performances_avg_over_taxa_n.rds"))
  saveRDS(output_performances_avg_over_n, paste0(results_dir, "n_", args$N, "/output_performances_avg_over_n.rds"))
  saveRDS(output_performances_avg_over_taxa_n_mn, paste0(results_dir, "n_", args$N, "/output_performances_avg_over_taxa_n_mn.rds"))
  saveRDS(output_performances_avg_over_n_mn, paste0(results_dir, "n_", args$N, "/output_performances_avg_over_n_mn.rds"))
  saveRDS(output_performances_avg_over_taxa_n_var, paste0(results_dir, "n_", args$N, "/output_performances_avg_over_taxa_n_var.rds"))
  saveRDS(output_performances_avg_over_n_var, paste0(results_dir, "n_", args$N, "/output_performances_avg_over_n_var.rds"))
  saveRDS(output_performances_avg_over_taxa_n_nofilter, paste0(results_dir, "n_", args$N, "/output_performances_avg_over_taxa_n_nofilter.rds"))
  saveRDS(output_performances_avg_over_n_nofilter, paste0(results_dir, "n_", args$N, "/output_performances_avg_over_n_nofilter.rds"))
  saveRDS(e_output_performances, paste0(results_dir, "n_", args$N, "/e_output_performances.rds"))
  saveRDS(e_output_performances_avg_over_n, paste0(results_dir, "n_", args$N, "/e_output_performances_avg_over_n.rds"))
  saveRDS(e_output_performances_avg_over_n_taxa, paste0(results_dir, "n_", args$N, "/e_output_performances_avg_over_n_taxa.rds"))
} else {
  output_performances <- readRDS(paste0(results_dir, "n_", args$N, "/output_performances.rds"))
  output_performances_avg_over_taxa_n <- readRDS(paste0(results_dir, "n_", args$N, "/output_performances_avg_over_taxa_n.rds"))
  output_performances_avg_over_n <- readRDS(paste0(results_dir, "n_", args$N, "/output_performances_avg_over_n.rds"))
  output_performances_avg_over_taxa_n_mn <- readRDS(paste0(results_dir, "n_", args$N, "/output_performances_avg_over_taxa_n_mn.rds"))
  output_performances_avg_over_n_mn <- readRDS(paste0(results_dir, "n_", args$N, "/output_performances_avg_over_n_mn.rds"))
  output_performances_avg_over_taxa_n_var <- readRDS(paste0(results_dir, "n_", args$N, "/output_performances_avg_over_taxa_n_var.rds"))
  output_performances_avg_over_n_var <- readRDS(paste0(results_dir, "n_", args$N, "/output_performances_avg_over_n_var.rds"))
  output_performances_avg_over_taxa_n_nofilter <- readRDS(paste0(results_dir, "n_", args$N, "/output_performances_avg_over_taxa_n_nofilter.rds"))
  output_performances_avg_over_n_nofilter <- readRDS(paste0(results_dir, "n_", args$N, "/output_performances_avg_over_n_nofilter.rds"))
  e_output_performances <- readRDS(paste0(results_dir, "n_", args$N, "/e_output_performances.rds"))
  e_output_performances_avg_over_n <- readRDS(paste0(results_dir, "n_", args$N, "/e_output_performances_avg_over_n.rds"))
  e_output_performances_avg_over_n_taxa <- readRDS(paste0(results_dir, "n_", args$N, "/e_output_performances_avg_over_n_taxa.rds"))
}

# -----------------------------------------------------------
# get the output ready for plots
# -----------------------------------------------------------

# bind together
overall_performance <- do.call(rbind.data.frame, output_performances)

# averaging over different things; plus filters
perf_avg_over_all <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n)
perf_avg_over_all_mn <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n_mn)
perf_avg_over_all_var <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n_var)
perf_avg_over_all_nofilter <- do.call(rbind.data.frame, output_performances_avg_over_taxa_n_nofilter)

perf_avg_over_n <- do.call(rbind.data.frame, output_performances_avg_over_n)
perf_avg_over_n_mn <- do.call(rbind.data.frame, output_performances_avg_over_n_mn)
perf_avg_over_n_var <- do.call(rbind.data.frame, output_performances_avg_over_n_var)
perf_avg_over_n_nofilter <- do.call(rbind.data.frame, output_performances_avg_over_n_nofilter)

e_perf_avg_over_all <- do.call(rbind.data.frame, e_output_performances_avg_over_n_taxa)
e_perf_avg_over_n <- do.call(rbind.data.frame, e_output_performances_avg_over_n)

performance_by_type <- perf_avg_over_all_mn %>% 
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))
e_performance_by_type <- e_perf_avg_over_all %>% 
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))
# ----------------------------------------------------------------------------------
# set up plot options
# ----------------------------------------------------------------------------------
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
fig_width <- 75
fig_height <- 15
point_cex <- 3
point_dodge <- 0.5
axis_cex <- 1.75

# ----------------------------------------------------------------------------------
# create plots
# ----------------------------------------------------------------------------------
make_fct_e <- function(.x, .y) {
  covers <- .x[.y == "cover"]
  mean(covers[order(covers, na.last = FALSE)])
}
e_main_plot_tib <- e_performance_by_type %>%
  filter(estimator == "ve") %>% 
  select(-q, -q_obs, -mse, -mspe, -wald_cover, -wald_cred_width,
         -pred_cover, -pred_width, -rmspe) %>% 
  pivot_longer(cols = tidyselect::all_of(c("cover", "bias",
                                           "cred_width", 
                                           "rmse")), 
               names_to = "summary", values_to = "value") %>% 
  mutate(object = factor("e",
                         levels = c("e", "V", "mu"),
                         labels = c("e", "V", "mu"))) 
main_plot_tib <- performance_by_type %>%
  filter(estimator == "ve") %>% 
  select(-q, -q_obs, -mse, -mspe) %>% 
  pivot_longer(cols = tidyselect::all_of(c("cover", "wald_cover", "bias",
                                           "cred_width", "wald_cred_width",
                                           "pred_cover", "pred_width", 
                                           "rmse", "rmspe")), 
               names_to = "summary", values_to = "value") %>% 
  mutate(object = factor(ifelse(grepl("p", summary), "V", "mu"),
                         levels = c("mu", "V", "e"),
                         labels = c(bquote(mu), "V", "e"))) %>% 
  bind_rows(e_main_plot_tib) %>% 
  mutate(fct_e = forcats::fct_reorder2(factor(paste0("sigma-e-", alpha_sigma, "-", kappa_sigma),
                                             levels = e_nms, 
                                             labels = e_labels,
                                             ordered = TRUE), value, summary, .fun = make_fct_e, .desc = TRUE))
x_lab <- bquote("Combination of "*alpha[sigma]*"/"*kappa[sigma])
cover_plt <- main_plot_tib %>%
  group_by(object, fct_e) %>% 
  filter(summary %in% c("cover", "pred_cover")) %>% 
  summarize(mn_cover = mean(value), se_cover = sqrt(var(value) * sqrt(100) / sqrt(50))) %>% 
  ggplot(aes(x = fct_e, y = mn_cover, shape = object)) +
  xlab(x_lab) +
  ylab(expression(Coverage)) + 
  ylim(c(0, 1)) +
  ggtitle("Interval coverage") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  # geom_errorbar(aes(ymin = mn_cover - 1.96*se_cover, ymax = mn_cover + 1.96*se_cover), 
  #               position = position_dodge(point_dodge)) +
  scale_shape_manual(values = pchs[c(2,1,3)], labels = c(bquote(mu), "V", "e")) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(shape = "Object to estimate") +
  theme(legend.position = c(0, 0.45),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing = unit(.1, "cm"),
        legend.box.background = element_rect(colour = "black"))

mse_plt <- main_plot_tib %>%
  group_by(object, fct_e) %>% 
  filter(summary %in% c("rmse", "rmspe")) %>% 
  summarize(mn_rmse = mean(value)) %>% 
  ggplot(aes(x = fct_e, y = mn_rmse, shape = object)) +
  xlab(x_lab) +
  ylab("RMS(P)E") + 
  ggtitle("Root mean squared (prediction) error") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_shape_manual(values = pchs[c(2,1,3)]) +
  guides(shape = FALSE)

width_plt <- main_plot_tib %>% 
  group_by(object, fct_e) %>% 
  filter(summary %in% c("cred_width", "pred_width")) %>% 
  summarize(mn_width = mean(value)) %>% 
  ggplot(aes(x = fct_e, y = mn_width, shape = object)) +
  xlab(x_lab) +
  ylab("Interval width") + 
  ggtitle("Mean interval width") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_shape_manual(values = pchs[c(2,1,3)]) +
    guides(shape = FALSE)

plot_grid(cover_plt, mse_plt, width_plt,
          nrow = 1, ncol = 3,
          align = "hv", rel_heights = c(1.1, 1, 1),
          labels = "AUTO")
ggsave(filename = paste0(plots_dir, "e_hypers.png"),
       plot = plot_grid(cover_plt, mse_plt, width_plt,
                        nrow = 1, ncol = 3,
                        align = "hv", rel_heights = c(1.1, 1, 1),
                        labels = "AUTO"),
       width = fig_width, height = fig_height, units = "cm", dpi = 300)


# --------------------------------------
# check coverage by taxon id
# --------------------------------------
performance_by_type_taxa <- perf_avg_over_n %>% 
  filter(!is.na(q)) %>% 
  mutate(rmse = sqrt(mse), rmspe = sqrt(mspe))

get_indx_in_main <- function(model_label, lev) {
  which(model_label == lev)
}
indices <- apply(matrix(e_labels), 1, function(x) get_indx_in_main(x, levels(main_plot_tib$fct_e)))
num_each_mod <- performance_by_type_taxa %>% filter(estimator == "ve") %>% group_by(alpha_sigma, kappa_sigma) %>% summarize(n = n())
all_indices <- vector("numeric", sum(num_each_mod$n))
for (i in 1:length(e_nms)) {
  all_indices[grepl(e_nms[i], (performance_by_type_taxa %>% filter(estimator == "ve") %>% mutate(e_nm = paste0("sigma-e-", alpha_sigma, "-", kappa_sigma)) %>% pull(e_nm)))] <- indices[i]  
}
taxon_check_plot_tib <- performance_by_type_taxa %>%
  filter(estimator == "ve") %>% 
  mutate(mod_lab_indx = all_indices,
         fct_e = forcats::fct_reorder(factor(paste0("sigma-e-", alpha_sigma, "-", kappa_sigma), levels = e_nms, 
                                               labels = e_labels,
                                               ordered = FALSE), mod_lab_indx)) %>% 
  select(-q, -q_obs, -mse, -mspe) %>% 
  pivot_longer(cols = tidyselect::all_of(c("cover", "wald_cover", "bias",
                                           "cred_width", "wald_cred_width",
                                           "pred_cover", "pred_width", 
                                           "rmse", "rmspe")), 
               names_to = "summary", values_to = "value") %>% 
  mutate(object = factor(ifelse(grepl("p", summary), "V", "mu"),
                         levels = c("mu", "V"),
                         labels = c(bquote(mu), "V"), ordered = TRUE))

cover_by_taxon_plt <- taxon_check_plot_tib %>%
  group_by(object, fct_e, taxon_id) %>% 
  filter(summary %in% c("cover")) %>% 
  summarize(mn_cover = mean(value)) %>% 
  ggplot(aes(x = fct_e, y = mn_cover, shape = object, color = taxon_id)) +
  xlab(x_lab) +
  ylab(expression(Coverage)) + 
  ggtitle(expression(bold(paste("Interval coverage for ", mu, sep = "")))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_shape_manual(values = pchs[c(2,1)]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  guides(shape = FALSE, color = FALSE)
rmse_by_taxon_plt <- taxon_check_plot_tib %>%
  group_by(object, fct_e, taxon_id) %>% 
  filter(summary %in% c("rmse")) %>% 
  summarize(mn_rmse = mean(value)) %>% 
  ggplot(aes(x = fct_e, y = mn_rmse, shape = object, color = taxon_id)) +
  xlab(x_lab) +
  ylab(expression(RMSE)) + 
  ggtitle(expression(bold(paste("Root mean squared error for ", mu, sep = "")))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_shape_manual(values = pchs[c(2,1)]) +
  guides(shape = FALSE) + 
  theme(legend.position = c(0, 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing = unit(.1, "cm"),
        legend.box.background = element_rect(color = "black"))
plot_grid(cover_by_taxon_plt, rmse_by_taxon_plt, labels = "AUTO",
          align = "hv",
          rel_heights = c(1, 1.1))
ggsave(filename = paste0(plots_dir, "e_hypers_cover_rmse_by_taxon.png"),
       plot = plot_grid(cover_by_taxon_plt, rmse_by_taxon_plt, labels = "AUTO",
                        align = "hv",
                        rel_heights = c(1, 1.1)),
       width = fig_width, height = fig_height, units = "cm", dpi = 300)

# ---------------------------------------------------------------
# check some specific MC ids
# ---------------------------------------------------------------
perf_avg_over_all %>% 
  mutate(fct_e = factor(paste0("sigma-e-", alpha_sigma, "-", kappa_sigma), levels = e_nms, 
                          labels = e_labels,
                          ordered = TRUE)) %>% 
  filter(estimator == "ve") %>% 
  ggplot(aes(x = fct_e, y = cover, color = factor(mc_id))) + 
  geom_point() + 
  geom_hline(yintercept = 0.95, color = "red") + 
  labs(color = "MC id") +
  xlab(x_lab)
perf_avg_over_all %>% 
  mutate(fct_e = factor(paste0("sigma-e-", alpha_sigma, "-", kappa_sigma), levels = e_nms, 
                          labels = e_labels,
                          ordered = TRUE)) %>% 
  filter(mc_id %in% c(9, 10, 11, 17), estimator == "ve") %>% 
  ggplot(aes(x = fct_e, y = cover, color = factor(mc_id))) + 
  geom_point() + 
  scale_color_manual(values = cols) + 
  geom_hline(yintercept = 0.95, color = "red") + 
  labs(color = "MC id") +
  xlab(x_lab)

