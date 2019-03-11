##################################################################################
## FILE: process_data.R
##
## CREATED: 14 August 2018 by Brian Williamson
##
## PURPOSE: read in data from Hutch  collaborators; create .rds files
##          with qPCR and br16s in final form. Also, create .rds files
##          with filled-in qPCR information.
##
## INPUTS: ~/Projects/HPTN/qPCR/data/qPCR/p2_brngs_qPCR_merge_20180314.csv
##         
## OUTPUTS: ~/Projects/HPTN/qPCR/data/final_qpcr.rds - dataset with final qpcr matrix
##          ~/Projects/HPTN/qPCR/data/final_br.rds - dataset with final br16s matrix
##          ~/Projects/HPTN/qPCR/data/imputed_qpcr.rds - dataset with imputed qpcr matrix
##          ~/Projects/HPTN/qPCR/data/person-specific_scaling_estimated_mean_qpcr.rds - mean for sims
##          ~/Projects/HPTN/qPCR/data/person-specific_scaling_estimated_cov_qpcr.rds - covariance for sims
##################################################################################

## set up data directory
project_dir <- "~/Projects/HPTN/qPCR/"
## read in the data
dat <- read.csv(file = paste0(project_dir, "data/p2_brngs_qPCR_merge_20180314.csv"))

## check how many people we have
dim(dat)[1]

## -------------------------------------------------------------------------------------------
## The actual study data
## -------------------------------------------------------------------------------------------
## Deal with qPCR data
## First column is id
tmp_qpcr <- dat[ , c(1, 481:494)]

## set qPCR values at detection limit to 0 (for now)
tmp_qpcr[, 2] <- ifelse(tmp_qpcr[, 2] <= tmp_qpcr[, 3], 0, tmp_qpcr[, 2])
tmp_qpcr[, 4] <- ifelse(tmp_qpcr[, 4] <= tmp_qpcr[, 5], 0, tmp_qpcr[, 4])
tmp_qpcr[, 6] <- ifelse(tmp_qpcr[, 6] <= tmp_qpcr[, 7], 0, tmp_qpcr[, 6])
tmp_qpcr[, 8] <- ifelse(tmp_qpcr[, 8] <= tmp_qpcr[, 9], 0, tmp_qpcr[, 8])
tmp_qpcr[, 10] <- ifelse(tmp_qpcr[, 10] <= tmp_qpcr[, 11], 0, tmp_qpcr[, 10])
tmp_qpcr[, 12] <- ifelse(tmp_qpcr[, 12] <= tmp_qpcr[, 13], 0, tmp_qpcr[, 12])
tmp_qpcr[, 14] <- ifelse(tmp_qpcr[, 14] <= tmp_qpcr[, 15], 0, tmp_qpcr[, 14])

## take only the qpcr information (not the *_t info)
qpcr_0 <- tmp_qpcr[, c(1,2,4,6,8,10,12,14)]

## Deal with br16s data
## First column is id
allbr_0 <- dat[, 1:476]
## reorder br so first 7 bugs correspond to qpcr
allbr_1 <- cbind(allbr_0[, c(1, 198, 476, 3, 2, 27, 64, 217)], allbr_0[, -c(1, 198, 476, 3, 2, 27, 64, 217)])
#drop br bugs with no counts
keep <- c(1, (2:476)[apply(allbr_1[, 2:476] > 0, 2, sum) > 0])
allbr <- allbr_1[, keep]

## compute the number of reads for each person from the br16s data
m = apply(allbr[, -1], 1, sum)

## delete samples with fewer than 1000 reads
qpcr <- qpcr_0[m >= 1000,]
br <- allbr[m >= 1000,]

## check how many people we have
dim(qpcr)
dim(br)

## We now have 7 qPCR bugs and 433 br16s bugs on 1325 data points

## check it out
head(qpcr)
head(br)

## create new datasets with removed nas for qpcr reads
na_qpcr <- apply(qpcr[, 2:8], 1, function(x) any(is.na(x)))
qpcr_final <- qpcr[!na_qpcr, ]
br_final <- br[!na_qpcr, ]

# get true means for 7 bugs
true_mean_qpcr <- colMeans(qpcr_final[, 2:8], na.rm = TRUE) 
names(true_mean_qpcr) <- c("Gardnerella vaginalis", "Lactobactillus crispatus",
                           "Lactobactillus iners", "Lactobactillus jensenii",
                           "Megasphaera genomosp.", "BVAB2 species", "Atopobium vaginae")

## ------------------------------------------
## get average scaling factor for each person
## ------------------------------------------
br_equals_zero <- (br_final == 0)[, 2:8]
qpcr_over_br16s <- qpcr_final[, 2:8]/br_final[, 2:8]

# on a row-by-row basis, do the following:
# remove infinite values
# take the mean
# append to avg vector
avg_scaling_per_person <- vector("numeric", dim(qpcr_over_br16s)[1])
for (i in 1:length(avg_scaling_per_person)) {
  sub_i <- as.numeric(qpcr_over_br16s[i, ])
  # get the row without infinite values
  avg_scaling_per_person[i] <- mean(sub_i[!is.infinite(sub_i)], na.rm = TRUE)
}

## --------------------
## get the imputed qPCR
## --------------------
## set up matrices
qpcr_mat <- as.matrix(qpcr_final[, 2:dim(qpcr_final)[2]])
br16s_mat <- as.matrix(br_final[, 2:dim(br_final)[2]])

## get the qPCR values where br16s == 0 for each person
sample_qpcr <- list()
for (i in 1:dim(qpcr_mat)[1]) {
  sample_qpcr[[i]] <- qpcr_mat[i, br16s_mat[i, 1:7] == 0]
}

## function to get imputed qpcrs
impute_qpcr <- function(br16s, scaling_factor, qpcr_sample) {
  # pre-allocate matrix
  ret <- matrix(NA, nrow = nrow(br16s), ncol = ncol(br16s))
  # loop through
  for (i in 1:dim(br16s)[1]) {
    # pick out the current participant
    sub_i <- br16s[i, ]
    # if any of the ones with actual qPCR information are zero, sample; otherwise, don't
    if (length(qpcr_sample[[i]]) > 0) {
      replace_vec_i <- sample(qpcr_sample[[i]], sum(sub_i == 0), replace = TRUE)
      ret[i, ] <- ifelse(sub_i == 0, replace_vec_i, sub_i*scaling_factor[i])
    } else {
      ret[i, ] <- sub_i*scaling_factor[i]
    }
  }
  return(ret)
}
# impute the qpcr
# first, sample a random value
# sample_qpcr <- qpcr_mat[br16s_mat[, 1:7] == 0]
set.seed(4747)
imputed_qpcr <- impute_qpcr(br16s_mat, avg_scaling_per_person, sample_qpcr)
imputed_qpcr[, 1:dim(qpcr_mat)[2]] <- qpcr_mat
colnames(imputed_qpcr) <- names(br_final)[2:dim(br_final)[2]]
# imputed_qpcr <- cbind(qpcr_final, br_final[, 9:dim(br_final)[2]]*avg_scaling)

## log mus
set.seed(1234)
log_mu <- log(imputed_qpcr)
# log_mu[is.infinite(log_mu)] <- rnorm(sum(is.infinite(log_mu)), log(0.5), 1)

# taxa-specific means
est_beta <- vector("numeric", dim(log_mu)[2]) 
for (j in 1:length(est_beta)) {
  sub_j <- log_mu[, j]
  est_beta[j] <- mean(sub_j[!is.infinite(sub_j)])
}
# est_beta <- colMeans(log_mu)

# covariance matrix -- for now, only do diagonal and use column-wise variances (removing infinite values)
est_sigma_vec <- vector("numeric", dim(log_mu)[2]) 
for (j in 1:length(est_sigma_vec)) {
  sub_j <- log_mu[, j]
  est_sigma_vec[j] <- var(sub_j[!is.infinite(sub_j)])
}
est_sigma <- diag(est_sigma_vec, ncol = length(est_sigma_vec))

## save final qpcr, br16s matrices
saveRDS(qpcr_final, paste0(project_dir, "data/final_qpcr.rds"))
saveRDS(br_final, paste0(project_dir, "data/final_br.rds"))
## save imputed qpcr matrix
saveRDS(imputed_qpcr, paste0(project_dir, "data/imputed_qpcr.rds"))
## save beta, sigma
saveRDS(est_beta, paste0(project_dir, "data/person-specific_scaling_estimated_mean_qpcr.rds"))
saveRDS(est_sigma, paste0(project_dir, "data/person-specific_scaling_estimated_cov_qpcr.rds"))

## -------------------------------------------------------------------------------------------
## Data for copy-number correction
## -------------------------------------------------------------------------------------------
## read in the copy number data 
copy_nums_full <- readr::read_csv(paste0(project_dir, "data/copy_numbers.csv"))

## get the ones for each of the 7 bugs
which(grepl("Lactobacillus crispatus", copy_nums_full$tax_name))

print(copy_nums_full[which(grepl("Gardnerella vaginalis", copy_nums_full$tax_name)), ], n = Inf)
print(copy_nums_full[which(grepl("Lactobacillus crispatus", copy_nums_full$tax_name)), ], n = 26)
print(copy_nums_full[which(grepl("Lactobacillus iners", copy_nums_full$tax_name)), ], n = Inf)
print(copy_nums_full[which(grepl("Lactobacillus jensenii", copy_nums_full$tax_name)), ], n = Inf)
print(copy_nums_full[which(grepl("Megasphaera genomosp. type_1", copy_nums_full$tax_name)), ], n = Inf)
print(copy_nums_full[which(grepl("BVAB2", copy_nums_full$tax_name)), ], n = Inf)
print(copy_nums_full[which(grepl("Atopobium vaginae", copy_nums_full$tax_name)), ], n = Inf)

## trying to find BVAB2
print(copy_nums_full[which(grepl("BVA", copy_nums_full$tax_name)), ], n = Inf)
print(copy_nums_full[which(grepl("BVAB", copy_nums_full$tax_name)), ], n = Inf)
## closest thing is BVAB3... use this?

## get the copy numbers for each one; the first row in the tables above
copy_nums <- c(copy_nums_full[which(grepl("Gardnerella vaginalis", copy_nums_full$tax_name)), ]$median[1], 
               copy_nums_full[which(grepl("Lactobacillus crispatus", copy_nums_full$tax_name)), ]$median[1],
               copy_nums_full[which(grepl("Lactobacillus iners", copy_nums_full$tax_name)), ]$median[1],
               copy_nums_full[which(grepl("Lactobacillus jensenii", copy_nums_full$tax_name)), ]$median[1],
               copy_nums_full[which(grepl("Megasphaera genomosp. type_1", copy_nums_full$tax_name)), ]$median[1],
               copy_nums_full[which(grepl("BVAB", copy_nums_full$tax_name)), ]$median[1],
               copy_nums_full[which(grepl("Atopobium vaginae", copy_nums_full$tax_name)), ]$median[1])
saveRDS(copy_nums, paste0(project_dir, "data/copy_numbers.rds"))
