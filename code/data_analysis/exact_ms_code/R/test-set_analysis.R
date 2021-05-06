# test-set analysis

# ------------------------------------------------------------------------------
# load required functions and packages
# ------------------------------------------------------------------------------
library("tidyverse")
library("cowplot")
theme_set(theme_cowplot())
library("rstan")
library("rstantools")
library("paramedic")
library("here")

# ------------------------------------------------------------------------------
# load the dataset and results of the full data analysis
# ------------------------------------------------------------------------------
div_num <- 1e3
analysis_data <- readRDS(here("data", "analysis_data.rds"))
qpcr <- analysis_data$qpcr %>%
  mutate_all(~ . / div_num)
br <- analysis_data$filtered_br16S
case_control <- (analysis_data$case_control %>%
                   mutate(case_bin = as.numeric(case == "Case")) %>%
                   select(-case))$case_bin

# static args
q <- ncol(br)
q_obs <- ncol(qpcr)

# read in the results; need to make these class "paramedic"
adjust <- 1
naive_results <- readRDS(
  here("results", "data_analysis", "full_analysis",
       paste0("qpcr_data_analysis_est_naive_q_", q, "_q_obs_", q_obs,
              "_sample_55_loo_999_adjust_", adjust, ".rds"))
)
no_ve_results <- readRDS(
  here("results", "data_analysis", "full_analysis",
       paste0("qpcr_data_analysis_est_no_ve_q_", q, "_q_obs_", q_obs,
              "_sample_55_loo_999_adjust_", adjust, ".rds"))
)
names(no_ve_results)[5] <- "stan_fit"
names(no_ve_results)[2] <- "summary"
names(no_ve_results)[3] <- "stan_data"
no_ve_results <- structure(no_ve_results, class = "paramedic")

ve_results <- readRDS(
  here("results", "data_analysis", "full_analysis",
       paste0("qpcr_data_analysis_est_ve_q_", q, "_q_obs_", q_obs,
              "_sample_55_loo_999_adjust_", adjust, ".rds"))
)
names(ve_results)[5] <- "stan_fit"
names(ve_results)[2] <- "summary"
names(ve_results)[3] <- "stan_data"
ve_results <- structure(ve_results, class = "paramedic")
sampled_women <- ve_results$stan_data$V[, 1]
test_women <- (1:110)[-sampled_women]

# ------------------------------------------------------------------------------
# fit the prediction model based on the original model output
# ------------------------------------------------------------------------------
# test data
test_W <- cbind.data.frame(subj_id = test_women,
                           br[test_women, , drop = FALSE])
test_V <- cbind.data.frame(subj_id = test_women,
                           qpcr[test_women, , drop = FALSE])
if (adjust) {
  test_X <- cbind.data.frame(subj_id = test_women,
                             case_control[test_women])  
} else {
  test_X <- data.frame(subj_id = test_women)
}
names(test_V) <- names(test_W)[1:(q_obs + 1)]
# get predictions: note that
# alpha_sigma = 4 & kappa_sigma = 3 -> varying-efficiency estimator
# alpha_sigma = 0 & kappa_sigma = 0 -> efficiency-naive estimator
# alpha_sigma = 2 & kappa_sigma = 3 -> sensitivity analysis
set.seed(4747)
seeds <- round(runif(n = 2, 1e3, 1e4))
set.seed(seeds[1])
pred_mod_ve <- posterior_predict(
  ve_results, W = test_W, V = test_V, X = test_X,
  alpha_sigma = 4, kappa_sigma = 3, 
  alpha_phi = 0, beta_phi = 0
)
set.seed(seeds[2])
pred_mod_no_ve <- posterior_predict(
  no_ve_results, W = test_W, V = test_V, X = test_X,
  alpha_sigma = 0, kappa_sigma = 0, 
  alpha_phi = 0, beta_phi = 0
)

# extract predictions of V
V_pred_means <- apply(pred_mod_ve$V, c(1, 2), mean)
V_pred_intervals <- aperm(apply(pred_mod_ve$V, c(1, 2), 
                          function(x) quantile(x, c(0.025, 0.975))),
                          c(2, 3, 1))
V_pred_means_no_ve <- apply(pred_mod_no_ve$V, c(1, 2), mean)
V_pred_intervals_no_ve <- aperm(apply(pred_mod_no_ve$V, c(1, 2), 
                                function(x) quantile(x, c(0.025, 0.975))),
                          c(2, 3, 1))
# prediction MSE
V_pred_means_obs <- V_pred_means[, 1:q_obs]
obs_mse_vec <- colMeans((V_pred_means_obs - as.matrix(test_V[, -1])) ^ 2)
obs_mse <- mean(obs_mse_vec)

V_pred_means_obs_no_ve <- V_pred_means_no_ve[, 1:q_obs]
obs_mse_vec_no_ve <- colMeans((V_pred_means_obs_no_ve - 
                                 as.matrix(test_V[, -1])) ^ 2)
obs_mse_no_ve <- mean(obs_mse_vec_no_ve)

# prediction coverage
V_pred_intervals_obs <- V_pred_intervals[, 1:q_obs, ]
obs_cover_vec <- colMeans(sapply(1:q_obs, function(j) {
  V_pred_intervals_obs[, j, 1] <= test_V[, -1][, j] & 
    test_V[, -1][, j] <= V_pred_intervals_obs[, j, 2]
}))
obs_cover <- mean(obs_cover_vec)
obs_cover

V_pred_intervals_obs_no_ve <- V_pred_intervals_no_ve[, 1:q_obs, ]
obs_cover_vec_no_ve <- colMeans(sapply(1:q_obs, function(j) {
  V_pred_intervals_obs_no_ve[, j, 1] <= test_V[, -1][, j] & 
    test_V[, -1][, j] <= V_pred_intervals_obs_no_ve[, j, 2]
}))
obs_cover_no_ve <- mean(obs_cover_vec_no_ve)
obs_cover_no_ve

taxon_nms <- c("A. christensenii", "A. vaginae",
               "BVAB2 spp.", "D. micraerophilus",
               "E. spp. type 1",
               "G. vaginalis", "L. crispatus",
               "L. iners", "L. jensenii",
               "M. hominis", "P. spp. type 1",
               "P. bennonis", "P. micra")

output_tib_ve <- tibble::tibble(
  taxon = taxon_nms,
  mse = obs_mse_vec,
  cover = obs_cover_vec,
  estimator = "ve"
)
output_tib_no_ve <- tibble::tibble(
  taxon = taxon_nms,
  mse = obs_mse_vec_no_ve,
  cover = obs_cover_vec_no_ve,
  estimator = "no_ve"
)
output_tib <- dplyr::bind_rows(
  output_tib_ve, output_tib_no_ve
) %>% 
  mutate(est_fct = factor(estimator,
                          levels = c("no_ve", "ve"),
                          labels = c("Efficiency-naive Bayes",
                                     "Varying-efficiency Bayes"))) %>% 
  arrange(taxon)
# ------------------------------------------------------------------------------
# make a similar plot to the one for the leave-one-out analysis
# ------------------------------------------------------------------------------
# set up colors
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

# create directory for plots if it doesn't already exist
plots_dir <- here("plots", "data_analysis", "test-set")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}
# plot coverage against taxon 
pred_interval_plot <- output_tib %>% 
  ggplot(aes(x = factor(taxon), 
             y = cover, shape = est_fct)) +
  xlab("Taxon") +
  ylab("Test-set coverage") +
  ggtitle("Prediction interval coverage") +
  labs(shape = "Estimator type") +
  geom_point(size = 4, position = position_dodge(width = 0.5, preserve = "total")) +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = c(17, 15)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0, 1)) +
  theme(legend.position = c(0.5, 0.15),
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

# plot mspe against left-out taxon (naive, no ve, ve)
mspe_plot <- output_tib %>% 
  ggplot(aes(x = factor(taxon), 
             y = mse, shape = est_fct)) +
  xlab("Taxon") +
  ylab("Test-set RMSPE") +
  scale_y_log10() +
  ggtitle("Root mean squared prediction error") +
  labs(color = "Estimator type") +
  geom_point(size = 4, position = position_dodge(width = 0.5, preserve = "total")) +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = c(17, 15)) +
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

# save the plot
ggsave(filename = here("plots", "data_analysis", "test-set", 
                       paste0("test-set_perf_adjust_", adjust, ".png")),
       plot = plot_grid(pred_interval_plot, mspe_plot),
       width = fig_width, height = 10)

