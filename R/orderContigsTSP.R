# Attempt to order contigs within 
# @useDynLib contiBAIT
# @import Rcpp TSP

# 
# @param linkageGroups list of vector of contig names for all linkage groups 
# @param allStrands table of strand calls for all contigs
# @param readTable table of W:C read proportions (used for QC)
# @param lg iterger specifying the element of linkageGroups to be used (i.e. specific linkage group to try to order)
# @return list of two members: 1) contig names in order, 2) contigs that were combined due to being zero-distance from each other


orderContigsTSP <- function(linkageGroups, allStrands, readTable, lg)
{

  linkageGroup <- linkageGroups[[lg]]
  #Combine zero-distance contigs:
  strands <- allStrands[linkageGroup,]
  contigs.combined <- combineZeroDistContigs(strands, readTable, lg)
  
  strands <- contigs.combined[[1]]
  dists <- as.matrix(daisy(strands))
  
  dists <- cbind(dists, rep(1, nrow(dists)))
  dists <- rbind(dists, rep(1, ncol(dists)))
  rownames(dists)[nrow(dists)] <- 'dummy'
  colnames(dists)[nrow(dists)] <- 'dummy'
  
  contig.tsp <- TSP(dists)
  contigsOrder <- solve_TSP(contig.tsp, method='2-opt', control=list(rep=1000))
  contigsOrder <- labels(contigsOrder)[-c(which(labels(contigsOrder)=='dummy'))]
  
  #Re-expand zero-distance contigs:
  expandedOrder <- expandMergedContigs(contigsOrder, contigs.combined[[2]],linkageGroup)
  
  #TODO: find a good way to return the list of combined contigs
  
  return(list(order=expandedOrder, zeroDist=contigs.combined[[2]]))
}