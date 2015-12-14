####################################################################################################
#' Function to order contigs within a single linkage group using a greedy algorithms
#' Attempt to order contigs within 
#' @useDynLib contiBAIT
#' @import Rcpp TSP
# 
#' @param linkageGroups list of vector of contig names for all linkage groups (product of clusterContigs)
#' @param allStrands table of strand calls for all contigs
#' @param readTable table of W:C read proportions (used for QC) (product of strandSeqFreqTable[[1]])
#' @param countTable table of read counts (product of strandSeqFreqTable[[2]])
#' @param lg iterger specifying the element of linkageGroups to be used (i.e. specific linkage group to try to order)
#' @param randomAttempts iterger specifying number of randomized clusterings to identify the best ordering. Default is 75
#' @example inst/examples/orderContigsGreedy.R
#' @export
#' @return list of two members: 1) contig names in order, 2) contigs that were combined due to being zero-distance from each other
####################################################################################################


orderContigsGreedy <- function(linkageGroups, allStrands, readTable, countTable, lg, randomAttempts=75, verbose=TRUE)
{  
  linkageGroup <- linkageGroups[[lg]]
  linkageGroupReadTable <- allStrands[linkageGroup,]

  zeroGroups <- combineZeroDistContigs(linkageGroupReadTable, readTable, lg)
  libWeight <- apply(countTable[which(rownames(countTable) %in% rownames(allStrands)  ),] , 1, median)

  zeroGroups[[2]]$weights <- libWeight[zeroGroups[[2]][,2]]
  linkageGroupReadTable <- zeroGroups[[1]]
  libWeight <- sapply(unique(zeroGroups[[2]][,1]), function(x) sum(zeroGroups[[2]]$weights[which(zeroGroups[[2]][,1] == x)]))

  for (i in 1:ncol(linkageGroupReadTable)){
    linkageGroupReadTable[,i] <- as.numeric(as.character( linkageGroupReadTable[,i]))
  }

  linkageGroupReadTable[is.na(linkageGroupReadTable)] <- 0

  best_order <- list(order = 1:nrow(linkageGroupReadTable),score = 0)
  temp_order <- list(order = 1:nrow(linkageGroupReadTable),score = 0)
  best_table <- linkageGroupReadTable

  if(nrow(linkageGroupReadTable) > 1)
  {
    if (!is.null(libWeight)){

  #    contigWeights <- libWeight[linkageGroup]
  #    contigWeights <- sort(contigWeights,decreasing = TRUE)
      contigWeights <- sort(libWeight, decreasing=TRUE)
      best_table <- as.matrix(linkageGroupReadTable[names(contigWeights),])
      best_order <- .Call('orderContigsGreedy', best_table)
    }
    for (i in 1:randomAttempts){
      #temp_order <- list(order = 1:length(linkageGroup),score = 0)
      temp_table <- as.matrix(linkageGroupReadTable[sample(nrow(linkageGroupReadTable)),])
      temp_order <- .Call('orderContigsGreedy', temp_table)

      if ( temp_order$score > best_order$score){
        if(verbose){message('     -> Found better ordering!')}  
        best_order <- temp_order
        best_table <- temp_table
      }
    }
  }
  order <- row.names(best_table)[best_order$order]

   mergedMatrix <- zeroGroups[[1]][order,]

  mergedGroups <- data.frame(LG=vector(), contig=vector())
  for(gp in 1:length(order)){
    mergedGroups <- rbind(mergedGroups, zeroGroups[[2]][which(zeroGroups[[2]] == order[gp]),1:2] )
  }
  mergedGroups <- new("ContigOrdering", mergedGroups)

  result <- list(order=order, mergedMatrix=mergedMatrix, mergedOrder=mergedGroups)
  result
}