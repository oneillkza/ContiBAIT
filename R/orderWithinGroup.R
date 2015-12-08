#' Attempt to order contigs within 
#' @useDynLib contiBAIT
#' @import Rcpp TSP
#' @export
#' 
#' @param linkageGroup vector of names of contigs in this linkage group to try to order
#' @param readTable table of W:C read proportions (used for QC)
#' @param allStrands table of strand calls for all contigs
#' @param lg numeric value defining which linkage groups
#' @param method algorithm to order the contigs
#' @param contigWeight optional order in which to place contigs using greedy algorithm
#' @return list of two members: 1) contig names in order, 2) contigs that were combined due to being zero-distance from each other


orderWithinGroup <- function(linkageGroup, readTable, allStrands, method = 'greedy', contigWeight = NULL)
{
  if (method == 'greedy'){
    return( orderContigsGreedy(linkageGroup,readTable,allStrands,contigWeight))
  } else {
  return (orderContigsTSP(linkageGroup,readTable,allStrands))
  }
}