
#' @useDynLib contiBAIT
#' @import Rcpp



orderContigsGreedy <- function(linkageGroup,readTable,allStrands,libWeight, randomAttempts = 75)
{  
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
    best_order <- .Call('orderContigsGreedy', best_table)
  }
  for (i in 1:randomAttempts){
temp_order <- list(order = 1:length(linkageGroup),score = 0)

    temp_table <- as.matrix(linkageGroupReadTable[sample(length(linkageGroup)),])
    temp_order <- .Call('orderContigsGreedy', temp_table)
    if ( temp_order$score >best_order$score){
      best_order <- temp_order
      best_table <- temp_table
    }
  }
  order <- row.names(best_table)[best_order$order]

  mergedOrder <- combineZeroDistContigs(allStrands[order,],readTable[[1]])

  result <- list(order =order,mergedOrder = mergedOrder)
  result
}