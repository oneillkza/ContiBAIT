####################################################################################################
#' Function to order contigs within a single linkage group using a greedy algorithms
#' Attempt to order contigs within 
#' @useDynLib contiBAIT
#' @import Rcpp TSP
# 
#' @param linkageGroupReadTable dataframe of strand calls (product of combineZeroDists or preprocessStrandTable)
#' @param randomAttempts number of times to repeat the greedy algortihm with a random restart
#' @param nProcesses number of processes to attempt ordering in parallel
#' @param verbose whether to print verbose messages
#' @return list of two members: 1) contig names in order, 2) the original data.frame entered into function correctly ordered  
####################################################################################################


orderContigsGreedy <- function(linkageGroupReadTable, randomAttempts=75,nProcesses = 1, verbose=TRUE)
{  
  factorizedLinkageGroupReadTable <- linkageGroupReadTable

  for (i in seq_len(ncol(linkageGroupReadTable))) {
    linkageGroupReadTable[,i] <- as.numeric(as.character( linkageGroupReadTable[,i]))
  }
  linkageGroupReadTable[is.na(linkageGroupReadTable)] <- 0

  if(nrow(linkageGroupReadTable) > 1)
  {
    order_contigs <- function(linkageGroupReadTable,randomAttempts,verbose,factorizedLinkageGroupReadTable)
    {
      best_order <- .Call('orderContigsGreedy', as.matrix(linkageGroupReadTable))
      best_table <- linkageGroupReadTable
      
      
      for (i in seq_len(randomAttempts)) {
        #temp_order <- list(order = 1:length(linkageGroup),score = 0)
        temp_table <- as.matrix(linkageGroupReadTable[sample(nrow(linkageGroupReadTable)),])
        temp_order <- .Call('orderContigsGreedy', temp_table)
        
        if ( temp_order$score < best_order$score){
          if(verbose){message('     -> Found better ordering!')}  
          best_order <- temp_order
          best_table <- temp_table
        }
      }
      
      linkageGroupReadTable <- factorizedLinkageGroupReadTable[row.names(best_table)[best_order$order],]
      return(c(list(orderVector=row.names(best_table)[best_order$order], orderedMatrix=linkageGroupReadTable),best_order$score))
    }
    if (nProcesses > 1)
    {
      cl <- makeCluster(getOption("cl.cores",nProcesses))
      
      
      clusterExport(cl,'contiBAIT')
      res <- clusterCall(cl,order_contigs,linkageGroupReadTable,randomAttempts,verbose,factorizedLinkageGroupReadTable)
      
      stopCluster(cl)
      best_order <- res[[1]]
      for (i in res[-1]){
        if (i[[3]] < best_order[[3]])
        {
          best_order <<- i
        }
        
      }
      return(best_order[-3])
    } else
    {
      return(order_contigs(linkageGroupReadTable,randomAttempts,verbose,factorizedLinkageGroupReadTable)[-3])
    }
  }else
  {
    linkageGroupReadTable <- factorizedLinkageGroupReadTable
    return(list(orderVector=row.names(linkageGroupReadTable), orderedMatrix=linkageGroupReadTable))
  }
}
