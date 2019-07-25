######################################################################################
##
## FILE: qpcr_analysis_helper_functions.R
##
## CREATED: 15 October 2018 by Brian Williamson
##
## PURPOSE: data analysis of the qPCR + br16s data from Hutch collaborators
##          to showcase qPCR estimation method; helper functions
##
######################################################################################

## FUNCTION: get_ids
## @param sample_id - the sample ids
## @return ids split into woman id (first five digits) and observation id (digits after hyphen)
get_ids <- function(sample_ids) {
  full_ids_lst <- strsplit(sample_ids, "-", fixed = TRUE)
  woman_ids <- unlist(lapply(full_ids_lst, function(x) x[1]))
  obs_ids <- unlist(lapply(full_ids_lst, function(x) x[2]))
  return(list(woman_ids = woman_ids, obs_ids = obs_ids))
}

## FUNCTION: process_data
## ARGS:  full_data - the full dataset
##          br_inds - indices of the br16s measurements
##        qpcr_inds - indices of the qPCR measurements
## pcr_plus_br_inds - indices in br16s that also have qpcr
##             llod - the lower limit of detection for qPCR
##            m_min - the minimum number of reads to consider
##          div_num - the number to divide qpcr by (should be approx lower limit of detection, e.g., 100)
## RETURNS: the qPCR and br16s data, cleaned up; only includes those with >= min number of reads
process_data <- function(full_data, br_inds, qpcr_inds, pcr_plus_br_inds, llod = 0, m_min = 1000, div_num = 100) {
  ## subset the qpcr data; set counts at lower limit of detection to llod
  sub_qpcr <- full_data[, qpcr_inds]
  clean_qpcr <- process_qpcr(sub_qpcr, grepl("_t", names(sub_qpcr), fixed = TRUE), grepl("_cps", names(sub_qpcr), fixed = TRUE), llod, div_num)
  ## add on ids
  tmp_qpcr <- data.frame(full_data[, 1], clean_qpcr)
  names(tmp_qpcr)[1] <- names(full_data)[1]

  ## subset the br16s data
  subset_br <- full_data[, br_inds]
  ## re-order so that the bugs with qpcr are first
  all_br <- cbind(subset_br[, pcr_plus_br_inds], subset_br[, -c(pcr_plus_br_inds)])
  ## check whether or not any have a nonzero count
  nonzero_counts <- apply(all_br > 0, 2, sum) > 0
  tmp_br <- data.frame(full_data[, 1], all_br[, nonzero_counts])
  names(tmp_br)[1] <- names(full_data)[1]

  ## remove any rows that have less than m_min reads
  m <- apply(tmp_br[, -1], 1, sum)
  qpcr_m <- tmp_qpcr[m >= m_min, ]
  br_m <- tmp_br[m >= m_min, ]
  ## remove any rows that have missing qpcr
  na_qpcr <- apply(qpcr_m[, -1], 1, function(x) any(is.na(x)))
  qpcr <- qpcr_m[!na_qpcr, ]
  br <- br_m[!na_qpcr, ]
  return(list(qpcr = qpcr, br = br))
}

## FUNCTION: process_qpcr
## ARGS: qpcr - the qpcr matrix
##  llod_inds - the indices giving the lower limit of detection
##  qpcr_inds - the indices giving actual qpcr values
##       llod - value to set things at lower limit of detection to
## RETURNS: the qpcr matrix, with only actual values, and values at lower limit of detection equal to llod
process_qpcr <- function(qpcr, llod_inds, qpcr_inds, llod = 0, div_num = 100) {
  ## subset with qpcr values
  sub_qpcr <- qpcr[, qpcr_inds]/div_num
  sub_llod <- qpcr[, llod_inds]/div_num
  ## clean up with llod
  clean_qpcr <- mapply(function(x, y) ifelse(x <= y, llod, x), sub_qpcr, sub_llod)
  return(clean_qpcr)
}

## FUNCTION: interval_long_to_wide
## ARGS: intervals - long matrix of intervals
##          digits - number of digits to round to
## RETURNS: intervals as a wide matrix (one column for each taxon, one row for each woman)
interval_long_to_wide <- function(intervals, digits, n_taxa, n_women) {
  chr_intervals <- apply(intervals[, c(2, 3)], 1, function(x) paste0("[", round(x[1], digits), ", ", round(x[2], digits), "]"))
  chr_intervals_df <- data.frame(taxon = 1:n_taxa, interval = chr_intervals, stringsAsFactors = FALSE)
  wide_intervals_mat <- matrix(NA, nrow = n_women, ncol = n_taxa)
  for (i in 1:n_women) {
    wide_intervals_mat[i, ]  <- as.character(tidyr::spread(chr_intervals_df[1:n_taxa + (i-1)*n_taxa, ], taxon, interval))
  }
  return(wide_intervals_mat)
}
