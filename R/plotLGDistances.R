plotLGDistances.func <- function(object, allStrands)
{
  linkageGroups <- object
  linkageStrands <- data.frame(do.call(rbind, lapply(linkageGroups, computeConsensus, allStrands)))
  rownames(linkageStrands) <- sapply(1:nrow(linkageStrands), function(x){paste('LG', x, ' (', length(linkageGroups[[x]]), ')', sep='') })
  sim <- 1-as.matrix(daisy(data.frame(linkageStrands)))
  rownames(sim) <- rownames(linkageStrands)
  colnames(sim) <- rownames(linkageStrands)
  breaks <- c(0/15, 1/15, 2/15, 3/15, 4/15, 5/15, 6/15, 7/15, 8/15, 9/15, 10/15, 11/15, 12/15, 12/15, 14/15, 15/15)
  cols <- c("cyan","cyan3","blue","blue4","gray22","gray0","gray0","gray0","gray0","gray0","gray22","red4","red3","red","darkorange")
  heatmap.2(sim, trace='none', col=cols, breaks=breaks)	
}
####################################################################################################
#' plotLGDistances -- plots a heatmap of the distances between linkage groups 
#' @param object LinkageGroupList 
#' @param allStrands StrandStateMatrix for all linkageGroups (usually reoriented by reorientStrandTable)
#' @aliases plotLGDistances plotLGDistances,LinkageGroupList,LinkageGroupList-method
#' @export
#' @importFrom gplots heatmap.2
####################################################################################################

setMethod('plotLGDistances',
          signature = signature(object='LinkageGroupList', allStrands = 'StrandStateMatrix'),
          definition = plotLGDistances.func
)
