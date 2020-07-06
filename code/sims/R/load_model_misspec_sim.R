# load_model-misspec_sim.R
# load results from the model-misspecification simulation, create plots

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
plots_dir <- paste0("plots/misspecification/n_", args$N, "/")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# ------------------------------------------------------------------
# create a list of output for easy summaries and get those summaries
# ------------------------------------------------------------------
all_sim_combos <- expand.grid(mu_dist = c("gamma", "halft", "normal"), e_dist = c("gamma", "halft", "normal"), stringsAsFactors = FALSE) 
# scramble the rows so that they match the original ones
all_sim_combos <- all_sim_combos[c(1, 2, 4, 5, 3, 6, 7, 8, 9), ]
all_sim_combos$v_dist <- "poisson"
all_sim_combos$w_dist <- "mult"
last_combos <- data.frame(mu_dist = "normal", e_dist = "normal", 
                          v_dist = c("poisson", "negbin", "negbin"), 
                          w_dist = c("dirimult", "mult", "dirimult"))
final_sim_combos <- rbind.data.frame(all_sim_combos, last_combos)
# for now, don't have all results yet
model_nms_init <- c(paste0("misspec-", apply(final_sim_combos[-9, ], 1, function(x) paste(x, collapse = "-"))),
               paste0("spec-", paste(final_sim_combos[9, ], collapse = "-")))[c(1:8, 12, 9:11)]
model_nms <- model_nms_init
get_model_labs <- function(x) {
  gsub("diriM", "CDM", 
       gsub("negbin", "NB", 
            gsub("mult", "M",
                 gsub("normal", "N", 
                      gsub("poisson", "P", 
                           gsub("halft", "H-t", 
                                gsub("gamma", "G",     
                                     gsub("-", "/", 
                                          gsub("spec-", "",
                                               gsub("misspec-", "", x))))))))))
}
model_labels <- get_model_labs(model_nms)
results_dir <- paste0("results/misspecification/")
if (args$read_data) {
    output_performances <- vector(mode = "list", length = length(model_nms))

    output_performances_avg_over_taxa_n <- output_performances_avg_over_taxa_n_mn <- output_performances_avg_over_taxa_n_var <- output_performances_avg_over_taxa_n_nofilter <- vector(mode = "list", length = length(model_nms))
    output_performances_avg_over_n <- output_performances_avg_over_n_mn <- output_performances_avg_over_n_var <- output_performances_avg_over_n_nofilter <- vector(mode = "list", length = length(model_nms))
    
    data_tib <- NULL
    for (i in 1:length(model_nms)) {
        dir_mat <- expand.grid(job = 1:args$B, dir = paste0(results_dir, model_nms[i], "/n_", args$N, "/q_", args$q, "/q_obs_", args$q_obs, "/"))
        mod_nms_lst <- as.list(paste0(dir_mat$dir, "run_paramedic_mod_jobid_", dir_mat$job, "_ad_", args$adapt_delta, "0000_mt_", args$max_treedepth, ".rds"))
        data_nms_lst <- as.list(paste0(dir_mat$dir, "run_paramedic_data_jobid_", dir_mat$job, "_ad_", args$adapt_delta, "0000_mt_", args$max_treedepth, ".rds"))
        samp_nms_lst <- as.list(paste0(dir_mat$dir, "run_paramedic_samps_jobid_", dir_mat$job, "_ad_", args$adapt_delta, "0000_mt_", args$max_treedepth, ".rds"))

        mod_lst <- lapply(mod_nms_lst, read_func)
        data_lst <- lapply(data_nms_lst, read_func)
        mu_tib <- do.call(rbind.data.frame, sapply(1:length(data_lst[!is.na(data_lst)]), function(x) make_data_tibble(data_lst[!is.na(data_lst)], x, obj = "mu", model_nm = model_nms[i]), simplify = FALSE))
        V_tib <- do.call(rbind.data.frame, sapply(1:length(data_lst[!is.na(data_lst)]), function(x) make_data_tibble(data_lst[!is.na(data_lst)], x, obj = "V", model_nm = model_nms[i]), simplify = FALSE))
        W_tib <- do.call(rbind.data.frame, sapply(1:length(data_lst[!is.na(data_lst)]), function(x) make_data_tibble(data_lst[!is.na(data_lst)], x, obj = "W", model_nm = model_nms[i]), simplify = FALSE))
        e_tib <- do.call(rbind.data.frame, sapply(1:length(data_lst[!is.na(data_lst)]), function(x) make_data_tibble(data_lst[!is.na(data_lst)], x, obj = "e", model_nm = model_nms[i]), simplify = FALSE))
        this_data_tib <- mu_tib %>% 
          left_join(V_tib, by = c("taxon_id", "mc_id", "subj_id", "model", "observed", "unobserved")) %>% 
          left_join(W_tib, by = c("taxon_id", "mc_id", "subj_id", "model", "observed", "unobserved")) %>% 
          left_join(e_tib, by = c("taxon_id", "mc_id", "model", "observed", "unobserved"))
        data_tib <- data_tib %>% 
          bind_rows(this_data_tib)
        mc_ids <- 1:args$B
        
        ## (1) pair model summaries with relevant data, for all taxa, all q_obs; replace NA with samps_lst if you want
        summary_df <- do.call(rbind.data.frame, mapply(function(w, x, y, z, type, mc_id) get_summaries(mod_summ = w, data = x, samps = y, mc_id = mc_id, indx = z, type), w = mod_lst, x = data_lst, y = NA, mc_id = mc_ids, MoreArgs = list(z = 1:args$q, type = "ve"), SIMPLIFY = FALSE))
        summary_df <- as_tibble(summary_df) %>% 
          filter(!is.na(summary_df$q)) %>% 
          mutate(model = model_nms[i])
        
        e_summary_df <- do.call(rbind.data.frame, mapply(function(w, x, y, mc_id) get_e_summaries(mod_summ = w, data = x, mc_id = mc_id, indx = y), w = mod_lst, x = data_lst, mc_id = mc_ids, MoreArgs = list(y = 1:args$q), SIMPLIFY = FALSE)) %>% 
          filter(!is.na(q)) %>% 
          mutate(model = model_nms[i]) %>% 
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
          
        ## (3) average over MC reps for each taxon
        mc_averaged_performance <- get_average_performance(performance_df %>% filter(include_taxa_in_avg),
                                                           group = "drop_mc")
        mc_averaged_performance_nofilter <- get_average_performance(performance_df,
                                                                    group = "drop_mc")
        ## check how many naives got filtered out since the intervals were NA
        num_na_naive_over_mc <- performance_df %>%
          filter(estimator == "naive") %>% 
          group_by(q, q_obs, subj_id, taxon_id, model) %>%
          summarize(perc_na = mean(is.na(cover)))
        num_na_naive_by_taxon_id <- num_na_naive_over_mc %>%
          group_by(q, q_obs, taxon_id, model) %>%
          summarize(perc_na = mean(perc_na))
        num_na_naive <- num_na_naive_by_taxon_id %>%
          group_by(q, q_obs, model) %>%
          summarize(perc_na = mean(perc_na))
        print(num_na_naive)
        ## (4) set flag for taxa of interest: for V, need taxon_id > q^obs
        mc_averaged_performance <- mc_averaged_performance %>% 
          mutate(include_taxa_in_v_avg = taxon_id > q_obs)
        mc_averaged_performance_nofilter <- mc_averaged_performance_nofilter %>% 
          mutate(include_taxa_in_v_avg = taxon_id > q_obs)
      
        ## (5) average over taxa of interest for each q_obs
        performance_across_taxa <- get_performance_filter(mc_averaged_performance,
                                                          group = "taxa")
        ## without filter
        performance_across_taxa_nofilter <- get_performance_filter(mc_averaged_performance_nofilter,
                                                                   group = "taxa")
        ## (6) average over n, for each q_obs
        average_over_n <- get_average_performance(performance_across_taxa,
                                                 group = "drop_n")
        ## no filter
        average_over_n_nofilter <- get_average_performance(performance_across_taxa_nofilter,
                                                           group = "drop_n")
        ## (5) add rmse, transpose
        performance_matrix <- average_over_n %>%
          mutate(rmse = sqrt(mse),
                 rmspe = sqrt(mspe))
        performance_matrix_nofilter <- average_over_n_nofilter %>%
          mutate(rmse = sqrt(mse),
             rmspe = sqrt(mspe))

        ## (6) average only over n, taxa
        ## no filter
        average_over_n_taxa_nofilter <- get_performance_filter(performance_df, 
                                                      group = "taxa+subj")
        ## (a) filter by colSums(W) > 0
        average_over_n_taxa <- get_performance_filter(performance_df %>% filter(include_taxa_in_avg),
                                                      group = "taxa+subj")
        ## (b) filter by colMeans(W) > 0.5
        average_over_n_taxa_mn <- get_performance_filter(performance_df %>% filter(include_taxa_in_avg_w_mn),
                                                      group = "taxa+subj")
        ## (c) filter by colVars(W) > 1
        average_over_n_taxa_var <- get_performance_filter(performance_df %>% filter(include_taxa_in_avg_var),
                                                         group = "taxa+subj")
        
        ## (7) average only over n
        ## (a) filter by colSums(W) > 0
        average_over_n <- get_performance_filter(performance_df %>% filter(include_taxa_in_avg),
                                                 group = "subj")
        ## (b) filter by colMeans(W) > 0
        average_over_n_mn <- get_performance_filter(performance_df %>% filter(include_taxa_in_avg_w_mn),
                                                 group = "subj")
        ## (a) filter by colVars(W) > 0
        average_over_n_var <- get_performance_filter(performance_df %>% filter(include_taxa_in_avg_var),
                                                     group = "subj")
        ## no filter
        average_over_n_nofilter <- get_performance_filter(performance_df,
                                                          group = "subj")
        
        output_performances[[i]] <- performance_matrix
        output_performances_avg_over_taxa_n[[i]] <- average_over_n_taxa
        output_performances_avg_over_n[[i]] <- average_over_n
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
  saveRDS(data_tib, paste0(results_dir, "n_", args$N, "/data.rds"))
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
  data_tib <- readRDS(paste0(results_dir, "n_", args$N, "/data.rds"))
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

performance_by_type <- perf_avg_over_all_mn %>% 
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
main_plot_tib <- performance_by_type %>%
  filter(estimator == "ve") %>% 
  mutate(fct_mod = forcats::fct_reorder(factor(model, levels = model_nms, 
                          labels = model_labels,
                          ordered = TRUE), cover, .desc = TRUE)) %>% 
  select(-q, -q_obs, -mse, -mspe) %>% 
  pivot_longer(cols = tidyselect::all_of(c("cover", "wald_cover", "bias",
                                           "cred_width", "wald_cred_width",
                                           "pred_cover", "pred_width", 
                                           "rmse", "rmspe")), 
               names_to = "summary", values_to = "value") %>% 
  mutate(object = factor(ifelse(grepl("p", summary), "V", "mu"),
                         levels = c("mu", "V"),
                         labels = c(bquote(mu), "V"), ordered = TRUE))
cover_plt <- main_plot_tib %>%
    group_by(object, fct_mod) %>% 
    filter(summary %in% c("cover", "pred_cover")) %>% 
    summarize(mn_cover = mean(value)) %>% 
    ggplot(aes(x = fct_mod, y = mn_cover, shape = object)) +
    xlab(expression(paste("True data-generating mechanism combination (", mu, "/e/V/W)", sep = ""))) +
    ylab(expression(Coverage)) + 
    ylim(c(0, 1)) +
    ggtitle("Interval coverage") +
    geom_point(position = position_dodge(point_dodge), size = point_cex) +
    scale_shape_manual(values = pchs[c(2,1)], labels = c(bquote(mu), "V")) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    guides(x = guide_axis(n.dodge = 2)) +
    labs(shape = "Object to estimate",
         caption = "N = normal, P = poisson, M = mult., CDM = compound dirichlet mult., H-t = half-t, G = gamma, NB = neg. bin.") +
    theme(legend.position = c(0, 0.45),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.spacing = unit(.1, "cm"),
          legend.box.background = element_rect(colour = "black"))

mse_plt <- main_plot_tib %>%
    group_by(object, fct_mod) %>% 
    filter(summary %in% c("rmse", "rmspe")) %>% 
    summarize(mn_rmse = mean(value)) %>% 
    ggplot(aes(x = fct_mod, y = mn_rmse, shape = object)) +
    xlab(expression(paste("True data-generating mechanism combination (", mu, "/e/V/W)", sep = ""))) +
    ylab("RMS(P)E") + 
    ggtitle("Root mean squared (prediction) error") +
    geom_point(position = position_dodge(point_dodge), size = point_cex) +
    scale_shape_manual(values = pchs[c(2,1)]) +
    guides(shape = FALSE, x = guide_axis(n.dodge = 2))

width_plt <- main_plot_tib %>% 
  group_by(object, fct_mod) %>% 
  filter(summary %in% c("cred_width", "pred_width")) %>% 
  summarize(mn_width = mean(value)) %>% 
  ggplot(aes(x = fct_mod, y = mn_width, shape = object)) +
  xlab(expression(paste("True data-generating mechanism combination (", mu, "/e/V/W)", sep = ""))) +
  ylab("Interval width") + 
  ggtitle("Mean interval width") +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_shape_manual(values = pchs[c(2,1)]) +
  guides(shape = FALSE, x = guide_axis(n.dodge = 2))

plot_grid(cover_plt, mse_plt, width_plt,
          nrow = 1, ncol = 3,
          align = "hv", rel_heights = c(1.1, 1, 1),
          labels = "AUTO")
ggsave(filename = paste0(plots_dir, "misspec.png"),
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
indices <- apply(matrix(model_labels), 1, function(x) get_indx_in_main(x, levels(main_plot_tib$fct_mod)))
num_each_mod <- performance_by_type_taxa %>% filter(estimator == "ve") %>% group_by(model) %>% summarize(n = n())
all_indices <- vector("numeric", sum(num_each_mod$n))
for (i in 1:length(model_nms)) {
  all_indices[grepl(model_nms[i], (performance_by_type_taxa %>% filter(estimator == "ve"))$model)] <- indices[i]  
}
taxon_check_plot_tib <- performance_by_type_taxa %>%
  filter(estimator == "ve") %>% 
  mutate(mod_lab_indx = all_indices,
    fct_mod = forcats::fct_reorder(factor(model, levels = model_nms, 
                                               labels = model_labels,
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
  group_by(object, fct_mod, taxon_id) %>% 
  filter(summary %in% c("cover")) %>% 
  summarize(mn_cover = mean(value)) %>% 
  ggplot(aes(x = fct_mod, y = mn_cover, shape = object, color = taxon_id)) +
  xlab(expression(paste("True data-generating mechanism combination (", mu, "/e/V/W)", sep = ""))) +
  ylab(expression(Coverage)) + 
  ggtitle(expression(bold(paste("Interval coverage for ", mu, sep = "")))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_shape_manual(values = pchs[c(2,1)]) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  guides(shape = FALSE, color = FALSE, x = guide_axis(n.dodge = 2))
rmse_by_taxon_plt <- taxon_check_plot_tib %>%
  group_by(object, fct_mod, taxon_id) %>% 
  filter(summary %in% c("rmse")) %>% 
  summarize(mn_rmse = mean(value)) %>% 
  ggplot(aes(x = fct_mod, y = mn_rmse, shape = object, color = taxon_id)) +
  xlab(expression(paste("True data-generating mechanism combination (", mu, "/e/V/W)", sep = ""))) +
  ylab(expression(RMSE)) + 
  ggtitle(expression(bold(paste("Root mean squared error for ", mu, sep = "")))) +
  geom_point(position = position_dodge(point_dodge), size = point_cex) +
  scale_shape_manual(values = pchs[c(2,1)]) +
  labs(caption = "N = normal, P = poisson, M = mult., CDM = compound dirichlet mult., H-t = half-t, G = gamma, NB = neg. bin.") +
  guides(shape = FALSE, x = guide_axis(n.dodge = 2)) + 
  theme(legend.position = c(0, 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing = unit(.1, "cm"),
        legend.box.background = element_rect(color = "black"))
plot_grid(cover_by_taxon_plt, rmse_by_taxon_plt, labels = "AUTO",
          align = "hv",
          rel_heights = c(1, 1.1))
ggsave(filename = paste0(plots_dir, "misspec_cover_rmse_by_taxon.png"),
       plot = plot_grid(cover_by_taxon_plt, rmse_by_taxon_plt, labels = "AUTO",
                        align = "hv",
                        rel_heights = c(1, 1.1)),
       width = fig_width, height = fig_height, units = "cm", dpi = 300)


# ---------------------------------------------------------------
# check some specific MC ids
# ---------------------------------------------------------------
perf_avg_over_all %>% 
  mutate(fct_mod = factor(model, levels = model_nms, 
                          labels = model_labels,
                          ordered = TRUE)) %>% 
  filter(estimator == "ve") %>% 
  ggplot(aes(x = fct_mod, y = cover, color = factor(mc_id))) + 
  geom_point() + 
  geom_hline(yintercept = 0.95, color = "red") + 
  labs(color = "MC id") +
  xlab("Model")
perf_avg_over_all %>% 
  mutate(fct_mod = factor(model, levels = model_nms, 
                          labels = model_labels,
                          ordered = TRUE)) %>% 
  filter(mc_id %in% c(9, 10, 11, 17), estimator == "ve") %>% 
  ggplot(aes(x = fct_mod, y = cover, color = factor(mc_id))) + 
  geom_point() + 
  scale_color_manual(values = cols) + 
  geom_hline(yintercept = 0.95, color = "red") + 
  labs(color = "MC id")

# ---------------------------------------------------------------
# diagnostic plots: distribution summaries of mu, e, V, W
# for each model
# ---------------------------------------------------------------
num_each_mod_data <- data_tib %>% group_by(model) %>% summarize(n = n())
all_data_indices <- vector("numeric", sum(num_each_mod_data$n))
for (i in 1:length(model_nms)) {
  all_data_indices[grepl(model_nms[i], data_tib$model)] <- indices[i]  
}
distn_plot_tib <- data_tib %>% 
  mutate(mod_lab_indx = all_data_indices,
         fct_mod = forcats::fct_reorder(factor(model, levels = model_nms, 
                          labels = model_labels,
                          ordered = TRUE), mod_lab_indx),
         fct_obs = factor(observed, levels = c("FALSE", "TRUE"), labels = c("Unobserved", "Observed")))
mu_boxplot <- distn_plot_tib %>% 
  # ggplot(aes(x = fct_mod, y = mu + 1, color = factor(mc_id))) +
  ggplot(aes(x = fct_mod, y = mu + 1, color = fct_obs)) +
  geom_boxplot() +
  labs(color = "V observed?",
       caption = "N = normal, P = poisson, M = mult., CDM = compound dirichlet mult., H-t = half-t, G = gamma, NB = neg. bin.") +
  xlab("Model") +
  ylab(expression(mu)) +
  scale_y_log10() +
  theme(legend.position = c(0, 0.9)) + 
  guides(x = guide_axis(n.dodge = 2))
v_boxplot <- distn_plot_tib %>% 
  ggplot(aes(x = fct_mod, y = V + 1, color = fct_obs)) +
  geom_boxplot() +
  labs(color = "V observed?") +
  xlab("Model") +
  ylab("V") +
  scale_y_log10() + 
  guides(color = FALSE, x = guide_axis(n.dodge = 2))
w_boxplot <- distn_plot_tib %>% 
  ggplot(aes(x = fct_mod, y = W + 1, color = fct_obs)) +
  geom_boxplot() +
  labs(color = "V observed?") +
  xlab("Model") +
  ylab("W") +
  scale_y_log10() + 
  guides(color = FALSE, x = guide_axis(n.dodge = 2))
e_boxplot <- distn_plot_tib %>% 
  mutate(e = value) %>% 
  filter(subj_id == 1) %>% # all the same
  ggplot(aes(x = fct_mod, y = e, color = fct_obs)) +
  geom_boxplot() +
  labs(color = "V observed?") +
  xlab("Model") +
  ylab("e") + 
  scale_y_sqrt(breaks = c(1, 10, 20, 40, 60, 80)) +
  guides(color = FALSE, x = guide_axis(n.dodge = 2))

plot_grid(mu_boxplot, v_boxplot, w_boxplot, e_boxplot, nrow = 2, ncol = 2, labels = "AUTO",
          rel_heights = c(1.1, 1, 1, 1), align = "hv")
ggsave(filename = paste0(plots_dir, "misspec_data_by_obs.png"),
       plot = plot_grid(mu_boxplot, v_boxplot, w_boxplot, e_boxplot, nrow = 2, ncol = 2, labels = "AUTO",
                        align = "hv",
                        rel_heights = c(1.1, 1, 1, 1)),
       width = fig_width, height = fig_height * 3, units = "cm", dpi = 300)
