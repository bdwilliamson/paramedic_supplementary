######################################################################################
##
## FILE: get_most_abundant_taxa.R
##
## CREATED: 19 November 2018 by Brian Williamson
##
## PURPOSE: return the indices of the most abundant taxa based on br16S
##
## INPUTS: 
## 
## OUTPUTS: 
######################################################################################

## FUNCTION: get_most_abundant_taxa
## ARGS: br16S_mat - the matrix of br16S reads
##             m - the vector of read counts for each participant
## RETURNS: the ordered column indices, from most abundant to least abundant
get_most_abundant_taxa <- function(br16S_mat, m) {
  ## divide by m, take the average down columns
  col_avg <- colMeans(br16S_mat/m)
  ## order, return
  ret <- order(col_avg, decreasing = TRUE)
  return(ret)
}