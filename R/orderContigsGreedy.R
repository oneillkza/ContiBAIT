####################################################################################################
#' Function to order contigs within a single linkage group using a greedy algorithms
#' Attempt to order contigs within 
#' @useDynLib contiBAIT
#' @import Rcpp TSP
# 
#' @param linkageGroups list of vector of contig names for all linkage groups (product of clusterContigs)
#' @param allStrands table of strand calls for all contigs
#' @param readTable list with table of W:C read proportions (used for QC) and read counts (product of strandSeqFreqTable)
#' @param lg iterger specifying the element of linkageGroups to be used (i.e. specific linkage group to try to order)
#' @param randomAttempts iterger specifying number of randomized clusterings to identify the best ordering. Default is 75
#' @export
#' @return list of two members: 1) contig names in order, 2) contigs that were combined due to being zero-distance from each other
####################################################################################################


orderContigsGreedy <- function(linkageGroup, allStrands, readTable, lg, randomAttempts=75)
{  
  linkageGroup <- linkageGroups[[lg]]
  linkageGroupReadTable <- allStrands[linkageGroup,]

  libWeight <- apply(readTable[[2]][which(rownames(readTable[[2]]) %in% rownames(allStrands)  ),] , 1, median)

  for (i in 1:ncol(linkageGroupReadTable)){
    linkageGroupReadTable[,i] <- as.numeric(as.character( linkageGroupReadTable[,i]))
  }

  linkageGroupReadTable[is.na(linkageGroupReadTable)] <- 0

  best_order <- list(order = 1:length(linkageGroup),score = 0)
   temp_order <- list(order = 1:length(linkageGroup),score = 0)
  best_table <- linkageGroupReadTable
  if (!is.null(libWeight)){
    contigWeights <- libWeight[linkageGroup]
    contigWeights <- sort(contigWeights,decreasing = TRUE)
    best_table <- as.matrix(linkageGroupReadTable[names(contigWeights),])
    best_order <- .Call('orderContigsGreedy', best_table)
  }
  for (i in 1:randomAttempts){
    #temp_order <- list(order = 1:length(linkageGroup),score = 0)
    temp_table <- as.matrix(linkageGroupReadTable[sample(length(linkageGroup)),])
    temp_order <- .Call('orderContigsGreedy', temp_table)

    print(temp_order$score)
    if ( temp_order$score > best_order$score){
      print(temp_order$score)
      best_order <- temp_order
      best_table <- temp_table
    }
  }
  order <- row.names(best_table)[best_order$order]

  mergedOrder <- combineZeroDistContigs(allStrands[order,], readTable[[1]], lg)

  result <- list(order=order, mergedOrder=mergedOrder)
  result
}