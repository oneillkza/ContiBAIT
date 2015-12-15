####################################################################################################
#' Function to order contigs within a single linkage group using a greedy algorithms
#' Attempt to order contigs within 
#' @useDynLib contiBAIT
#' @import Rcpp TSP
# 
#' @param linkageGroupReadTable dataframe of strand calls (product of combineZeroDists or preprocessStrandTable)
#' @return list of two members: 1) contig names in order, 2) the original data.frame entered into function correctly ordered  
####################################################################################################


orderContigsGreedy <- function(linkageGroupReadTable, randomAttempts=75, verbose=TRUE)
{  
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
  }

  #order the table entered into function and convert to factor, then add factor levels.
  linkageGroupReadTable[] <- lapply(linkageGroupReadTable[row.names(best_table)[best_order$order],], factor)
  linkageGroupReadTable <- data.frame(lapply(linkageGroupReadTable, function(x){levels(x) <- c(1,2,3); x}) )
  rownames(linkageGroupReadTable) <- row.names(best_table)[best_order$order]

  return(list(orderVector=row.names(best_table)[best_order$order], orderedMatrix=linkageGroupReadTable))
}