
#' @useDynLib contiBAIT
#' @import Rcpp


orderContigsGreedy <- function(linkageGroups,readTable,allStrands, lg, libWeight=NULL, randomAttempts = 75)
{
  linkageGroup <- linkageGroups[[lg]]
  linkageGroupReadTable <- allStrands[linkageGroup,]
  for (i in 1:ncol(linkageGroupReadTable)){
    linkageGroupReadTable[,i] <- as.numeric(as.character( linkageGroupReadTable[,i]))
  }

  linkageGroupReadTable[is.na(linkageGroupReadTable)] <- 0

  best_order <- list(order = 1:length(linkageGroup),score = 0)
  best_table <- linkageGroupReadTable
  if (!is.null(libWeight)){
  
    contigWeights <- libWeight[linkageGroup]
    contigWeights <- sort(contigWeights,decreasing = TRUE)
    best_table <- as.matrix(linkageGroupReadTable[names(contigWeights),])
    best_order$score <- .Call('orderContigsGreedy', best_table)
  }
  for (i in 1:randomAttempts){
temp_order <- list(order = 1:length(linkageGroup),score = 0)

    temp_table <- as.matrix(linkageGroupReadTable[sample(length(linkageGroup)),])
    temp_order$score <- .Call('orderContigsGreedy', temp_table)
    print(temp_order$score)
    if ( temp_order$score > best_order$score){
      print(temp_order$score)
      best_order <- temp_order
      best_table <- temp_table
    }
  }
  order <- row.names(best_table)[best_order$order]
#  plotContigOrder(row.names(best_table)[best_order$order])
  mergedOrder <- combineZeroDistContigs(allStrands[order,],readTable, lg)
#  plotContigOrder(mergedOrder[[2]])
  result <- list(order,mergedOrder[[1]], mergedOrder[[2]])
  result
}