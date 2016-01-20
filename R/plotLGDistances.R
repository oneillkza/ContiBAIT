plotLGDistances.func <- function(object, allStrands, lg='all', labels=TRUE)
{

if(lg[1] == 'all')
{
  linkageStrands <- data.frame(do.call(rbind, lapply(object, computeConsensus, allStrands)))
#  rownames(linkageStrands) <- sapply(1:nrow(linkageStrands), function(x){paste('LG', x, ' (', length(object[[x]]), ')', sep='') })
}else{
  linkageStrands <- allStrands[unlist(object[lg]),]
  #For brevity, only report chromosome name, not split locations.
  rownames(linkageStrands) <- paste(as.character(matrix(unlist(strsplit(rownames(linkageStrands), ':')), ncol=2, byrow=TRUE)[,1]), '_', seq(1,nrow(linkageStrands)), sep='')
}

  sim <- 1-as.matrix(daisy(data.frame(linkageStrands)))
  rownames(sim) <- rownames(linkageStrands)
  colnames(sim) <- rownames(linkageStrands)

  if(labels)
  {
    labRow=NULL
    labCol=NULL
  }else
  {
    labRow=rep('', nrow(sim))
    labCol=rep('', ncol(sim))
  }
  labelScale <- min(1, 32/nrow(sim))

  breaks <- seq(0, 100, length.out=101)/100 
  cols <- colorRampPalette(c("cyan", "blue", "grey30", "black", "grey30", "red", "orange"))
 
  if(length(lg) > 1)
  {
    chrLabels <- melt(object[lg])
    chrCols <- rainbow_hcl(length(unique(chrLabels[,2])), c=90, l=60)
    names(chrCols) <- unique(chrLabels[,2])
    rowCols <- chrCols[chrLabels[,2]]
    colCols <- rowCols
    heatmap.2(sim, 
              trace='none', 
              col=cols(100), 
              labRow=labRow, 
              breaks=breaks, 
              labCol=labCol, 
              RowSideColors=rowCols, 
              ColSideColors=colCols, 
              main=paste('Distances of ', nrow(sim),' linkage groups', sep=''))
    legend("topright",
            legend=unique(names(rowCols)),
            fill=unique(rowCols), 
            border=FALSE, 
            bty="n", 
            cex=0.7)
  }else{
    heatmap.2(sim, 
              trace='none', 
              col=cols(100), 
              breaks=breaks, 
              labRow=labRow, 
              labCol=labCol, 
              cexRow=labelScale, 
              cexCol=labelScale, 
              main=paste('Distances of ', nrow(sim) , ' linkage groups', sep=''))
  }
}

####################################################################################################
#' plotLGDistances -- plots a heatmap of the distances between linkage groups 
#' @param object LinkageGroupList 
#' @param allStrands StrandStateMatrix for all linkageGroups (usually reoriented by reorientStrandTable)
#' @param lg ='all' vector of integers to determine which linkage group(s) to plot. 'all' will calculate consensus
#' strand calls for all linkage groups and plot them side by side (default it 'all')
#' @param labels =TRUE if TRUE, contig names will be plotted on the axes
#' @param ... additional parameters to pass to heatmap.2
#' @aliases plotLGDistances plotLGDistances,LinkageGroupList,LinkageGroupList-method
#' @example inst/examples/plotLGDistances.R
#' @return a heatplot of linkage group calls
#' @export
#' @importFrom gplots heatmap.2
#' @importFrom colorspace rainbow_hcl
#' @importFrom reshape2 melt 
####################################################################################################

setMethod('plotLGDistances',
          signature = signature(object='LinkageGroupList', allStrands = 'StrandStateMatrix'),
          definition = plotLGDistances.func
)
