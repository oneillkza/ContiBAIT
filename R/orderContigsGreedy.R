
#' @useDynLib contiBAIT
#' @import Rcpp


orderContigsGreedy <- function(linkageGroups,readTable,allStrands,libWeight, lg)
{
  linkageGroup <- linkageGroups[[lg]]
  linkageGroupReadTable <- allStrands[linkageGroup,]
  for (i in 1:ncol(linkageGroupReadTable)){
    linkageGroupReadTable[,i] <- as.numeric(as.character( linkageGroupReadTable[,i]))
  }

  linkageGroupReadTable[is.na(linkageGroupReadTable)] <- 0
  if(length(libWeight) != 1)
  {
    contigWeights <- libWeight[linkageGroup]
    contigWeights <- sort(contigWeights,decreasing = TRUE)
    weightOrderedTable <- as.matrix(linkageGroupReadTable[names(contigWeights),])
  }else{
    weightOrderedTable <- as.matrix(linkageGroupReadTable)
  }
  order <- row.names(weightOrderedTable)[.Call('orderContigsGreedy', weightOrderedTable)]
  mergedOrder <- combineZeroDistContigs(allStrands[order,],readTable, lg)

  result <- list(order,mergedOrder[[1]], mergedOrder[[2]])
  result
}