##################################################################################
## FILE: naive_qpcr_estimator.R
##
## CREATED: 20 August 2018 by Brian Williamson
##
## PURPOSE: compute the naive estimator of qPCR from br16s
##
## INPUTS: 
##         
## OUTPUTS: 
##################################################################################

## ARGS:      idx - the index of interest (from brs)
##            brs - the br16s information
##          qpcrs - the qpcr information
##     known_qpcr - the indices with known qpcr
## RETURNS: estimated qPCR for bug idx
naive_estimator <- function(idx, brs, qpcrs, known_qpcr) {
  ## get the sum of the br16s's
  sum_br16s <- rowSums(brs[, known_qpcr, drop = FALSE])
  ## get the sum of the qpcr's
  sum_qpcr <- rowSums(qpcrs)
  
  ## naive estimator
  qpcr <- brs[, idx]*sum_qpcr/sum_br16s
  ## if sum_br16s was zero, predict zero (?)
  qpcr <- ifelse (sum_br16s == 0, 0, qpcr)
  return(qpcr)
}