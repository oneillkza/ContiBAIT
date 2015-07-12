
# @useDynLib contiBAIT
# @import Rcpp


orderContigsGreedy <- function(linkageGroup,readTable,allStrands,libWeight)
{
  
  linkageGroupReadTable <- allStrands[linkageGroup,]
  for (i in 1:ncol(linkageGroupReadTable)){
    linkageGroupReadTable[,i] <- as.numeric(as.character( linkageGroupReadTable[,i]))
  }
  linkageGroupReadTable[is.na(linkageGroupReadTable)] <- 0
  contigWeights <- libWeight[linkageGroup]
  contigWeights <- sort(contigWeights,decreasing = TRUE)
  weightOrderedTable <- as.matrix(linkageGroupReadTable[names(contigWeights),])
  order <- row.names(weightOrderedTable)[.Call('orderContigsGreedy', weightOrderedTable)]
  mergedOrder <- combineZeroDistContigs(allStrands[order,],readTable)
  result <- list(order,mergedOrder[[2]])
  result
}