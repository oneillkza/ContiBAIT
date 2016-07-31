#####################################################################################
#' Attempt to order contigs within linkage groups using travelling salesperson algorithm
#' @useDynLib contiBAIT
#' @import Rcpp TSP
#' @import TSP
# 
#' @param linkageGroupReadTable dataframe of strand calls (product of combineZeroDists or preprocessStrandTable)
#' @param reps Number of repetitions for he 2-opt TSP solver algorithm
#' @return list of two members: 1) contig names in order, 2) the original data.frame entered into function correctly ordered 
########################################################################################

orderContigsTSP <- function(linkageGroupReadTable, reps)
{

  dists <- as.matrix(daisy(linkageGroupReadTable))
  
  dists <- cbind(dists, rep(1, nrow(dists)))
  dists <- rbind(dists, rep(1, ncol(dists)))
  rownames(dists)[nrow(dists)] <- 'dummy'
  colnames(dists)[nrow(dists)] <- 'dummy'
  
  contig.tsp <- TSP(dists)
  contigsOrder <- solve_TSP(contig.tsp, method='2-opt', control=list(rep=reps))
  contigsOrder <- labels(contigsOrder)[-c(which(labels(contigsOrder)=='dummy'))]
  
  #order the table entered into function and convert to factor, then add factor levels.
  linkageGroupReadTable[] <- lapply(linkageGroupReadTable[contigsOrder,], factor)
  linkageGroupReadTable <- data.frame(lapply(linkageGroupReadTable, 
  										   function(x){levels(x) <- c(1,2,3); x}) )
  rownames(linkageGroupReadTable) <- contigsOrder

  return(list(orderVector=contigsOrder, orderedMatrix=linkageGroupReadTable))
}
