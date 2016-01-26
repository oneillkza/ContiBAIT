####################################################################################################
#' Function to order contigs within a single linkage group using a greedy algorithms
#' Attempt to order contigs within 
#' @useDynLib contiBAIT
#' @import Rcpp TSP
# 
#' @param linkageGroupReadTable dataframe of strand calls (product of combineZeroDists or preprocessStrandTable)
#' @param randomAttempts number of times to repeat the greedy algortihm with a random restart
#' @param verbose whether to print verbose messages
#' @return list of two members: 1) contig names in order, 2) the original data.frame entered into function correctly ordered  
####################################################################################################


orderContigsGreedy <- function(linkageGroupReadTable, randomAttempts=75, verbose=TRUE)
{  
  factorizedLinkageGroupReadTable <- linkageGroupReadTable

  for (i in 1:ncol(linkageGroupReadTable)){
    linkageGroupReadTable[,i] <- as.numeric(as.character( linkageGroupReadTable[,i]))
  }
  linkageGroupReadTable[is.na(linkageGroupReadTable)] <- 0

  if(nrow(linkageGroupReadTable) > 1)
  {
    best_order <- .Call('orderContigsGreedy', as.matrix(linkageGroupReadTable))
    best_table <- linkageGroupReadTable

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

  linkageGroupReadTable <- factorizedLinkageGroupReadTable[row.names(best_table)[best_order$order],]
  return(list(orderVector=row.names(best_table)[best_order$order], orderedMatrix=linkageGroupReadTable))

  }else{
    linkageGroupReadTable <- factorizedLinkageGroupReadTable
    return(list(orderVector=row.names(linkageGroupReadTable), orderedMatrix=linkageGroupReadTable))
  }
}