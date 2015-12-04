plotLGDistances.func <- function(object, allStrands)
{
  linkageGroups <- object
  linkageStrands <- data.frame(do.call(rbind, lapply(linkageGroups, computeConsensus, allStrands)))
  rownames(linkageStrands) <- sapply(1:nrow(linkageStrands), function(x){paste('LG', x, ' (', length(linkageGroups[[x]]), ')', sep='') })
  sim <- 1-as.matrix(daisy(data.frame(linkageStrands)))
  rownames(sim) <- rownames(linkageStrands)
  colnames(sim) <- rownames(linkageStrands)
  breaks <- seq(0, 100, length.out=101)/100 
  cols <- colorRampPalette(c("cyan", "blue", "grey30", "black", "grey30", "red", "orange"))
  heatmap.2(sim, trace='none', col=cols(100), breaks=breaks )
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
