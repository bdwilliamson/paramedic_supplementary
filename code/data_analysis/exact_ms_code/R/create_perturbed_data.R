## create perturbed data from UW STICRC data

code_dir <- "code/R/"
data_dir <- "data/"
stan_dir <- "stan/"
library("methods")
library("rstan")
library("argparse")
library("Cairo")
library("dplyr")
library("tibble")
source(paste0(code_dir, "analyze_data/qpcr_analysis_helper_functions.R"))
source(paste0(code_dir, "naive_qpcr_estimator.R"))
source(paste0(code_dir, "analyze_data/get_most_abundant_taxa.R"))

## -----------------------------------------------------------------------------------
## load in the data, clean it up
## -----------------------------------------------------------------------------------
full_data <- read.csv(file=paste0(data_dir, "p2_brngs_qPCR_merge_20180314.csv"),
                      stringsAsFactors = FALSE)
split_ids <- get_ids(full_data$SampleID)
full_data$ptid <- split_ids$woman_ids
full_data$obs <- split_ids$obs_ids
## set m_min, llod for processing
m_min <- 1000
llod <- 0
div_num <- 1000
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
br16_mat_nms <- colnames(br16_mat)
## calculate the read numbers
m <- rowSums(br16_mat)

## sample a subset of women
set.seed(4747)
samp <- sample(1:dim(br16_mat)[1], 20)
## set up q, q_obs
q <- dim(br16_mat)[2]
q_obs <- dim(qpcr_mat)[2]
qpcr_samp <- qpcr_mat[samp, ]
br_samp <- br16_mat[samp, ]

naive_est <- cbind(qpcr_mat[samp, 1:q_obs], apply(matrix((q_obs + 1):q), 1, naive_estimator, br16_mat[samp, ], qpcr_mat[samp, ], 1:q_obs))

## add random noise (Poisson(1) to qPCR, Mult(1, m[j], prob = naive_j/sum(naive_j)))
set.seed(4747)
qpcr_perturbed <- qpcr_samp + replicate(q_obs, rpois(dim(qpcr_samp)[1], 1))
br_perturbed <- br_samp 
for (i in 1:dim(br_perturbed)[1]) {
  br_perturbed[i, ] <- br_perturbed[i, ] + rmultinom(n = 1, size = m[i], prob = naive_est[i, ]/rowSums(naive_est)[i])
}

## add on new row names
rownames(br_perturbed) <- 1:dim(br_perturbed)[1]
rownames(qpcr_perturbed) <- 1:dim(qpcr_perturbed)[1]
colnames(br_perturbed) <- colnames(br_samp)
colnames(qpcr_perturbed) <- colnames(br_samp)[1:7]
mode(qpcr_perturbed) <- "integer"

example_data <- list(br = br_perturbed, qpcr = qpcr_perturbed)
saveRDS(example_data, paste0(data_dir, "full_example_data.rds"))
small_example_data <- list(br = br_perturbed[, 1:17], qpcr = qpcr_perturbed)
saveRDS(small_example_data, paste0(data_dir, "small_example_data.rds"))

## create tibbles
br_tibble <- tibble::as_tibble(br_perturbed) %>% 
  mutate(sample_id = 1:dim(br_perturbed)[1]) %>% 
  select(sample_id, colnames(br_perturbed))
qpcr_tibble <- tibble::as_tibble(qpcr_perturbed) %>% 
  mutate(sample_id = 1:dim(qpcr_perturbed)[1]) %>% 
  select(sample_id, colnames(qpcr_perturbed))

saveRDS(br_tibble, paste0(data_dir, "example_16S_data.rds"))
saveRDS(qpcr_tibble, paste0(data_dir, "example_qPCR_data.rds"))

## copy and paste this code if I need to update it in package:
example_16S_data <- readRDS("~/Projects/HPTN/qPCR/data/example_16S_data.rds")
example_qPCR_data <- readRDS("~/Projects/HPTN/qPCR/data/example_qPCR_data.rds")
usethis::use_data(example_16S_data, overwrite = TRUE)
usethis::use_data(example_qPCR_data, overwrite = TRUE)
