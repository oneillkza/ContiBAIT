plotLGDistances.func <- function(object, allStrands, lg='all', labels=TRUE, state='all', alreadyOrdered=FALSE)
{

  if(lg[1] == 'all')
  {
    linkageStrands <- data.frame(do.call(rbind, lapply(object, computeConsensus, allStrands)))
  #  rownames(linkageStrands) <- sapply(1:nrow(linkageStrands), function(x){paste('LG', x, ' (', length(object[[x]]), ')', sep='') })
  }else{
    linkageStrands <- allStrands[unlist(object[lg]),]
    #For brevity, only report chromosome name, not split locations.
    chrName <- as.character(matrix(unlist(strsplit(rownames(linkageStrands), ':')), ncol=2, byrow=TRUE)[,1])
    chrName <- paste(chrName, '_', seq(1,nrow(linkageStrands)), sep='')

    if(state == 'homo')
    {
      linkageStrands <- replace(linkageStrands, linkageStrands == 2, NA)
    }else if(state == 'hetero'){
      linkageStrands <- replace(linkageStrands, linkageStrands == 3, 1)
    }

    linkageStrands <- data.frame(linkageStrands)
    linkageStrands <-  data.frame(lapply(linkageStrands, function(x) factor(x, levels=c(1,2,3))))
    rownames(linkageStrands) <- chrName

  }

  sim <- suppressWarnings(1-as.matrix(daisy(data.frame(linkageStrands))))
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

  if(alreadyOrdered)
  {
    plotDen=FALSE
    dend='none'
  }else{
    plotDen=TRUE
    dend='both'
  }

  breaks <- seq(0, 100, length.out=101)/100 
  cols <- colorRampPalette(c("cyan", "blue", "grey30", "black", "grey30", "red", "orange"))

  if(lg[1] == 'all')
  {
    chrLabels <- data.frame(1:length(object), names(object))   
  }else{
   chrLabels <- melt(object[lg])
  }
  chrCols <- rainbow_hcl(length(unique(chrLabels[,2])), c=90, l=60)
  names(chrCols) <- unique(chrLabels[,2])
  rowCols <- chrCols[chrLabels[,2]]
  colCols <- rowCols
  heatmap.2(sim, 
            trace='none', 
            col=cols(100),
            Rowv=plotDen,
            Colv=plotDen,
            dendrogram=dend, 
            labRow=labRow, 
            breaks=breaks, 
            labCol=labCol, 
            cexRow=labelScale, 
            cexCol=labelScale,
            RowSideColors=rowCols, 
            ColSideColors=colCols, 
            main=paste('Distances of ', nrow(sim),' linkage groups', sep=''))
  legend("topright",
          legend=unique(names(rowCols)),
          fill=unique(rowCols), 
          border=FALSE, 
          bty="n", 
          cex=0.7)

}
 
####################################################################################################
#' plotLGDistances -- plots a heatmap of the distances between linkage groups 
#' @param object LinkageGroupList 
#' @param allStrands StrandStateMatrix for all linkageGroups (usually reoriented by reorientStrandTable)
#' @param lg ='all' vector of integers to determine which linkage group(s) to plot. 'all' will calculate consensus
#' strand calls for all linkage groups and plot them side by side (default it 'all')
#' @param labels =TRUE if TRUE, contig names will be plotted on the axes
#' @param state string denoting whether only homozygous states should be used ('homo'; just WW vs CC), if a comparison should be made between
#' homozygous and heterozygous states only ('hetero'; just homo vs hetero), or if all three states ('all'; WW vs WC vs CC). Default is 'all'
#' @param alreadyOrdered if TRUE, the function will assume that the linkageGroupList is already ordered and not 
#' create a dendrogram. Default is FALSE
#' @param ... additional parameters to pass to heatmap.2
#' @aliases plotLGDistances plotLGDistances,LinkageGroupList,LinkageGroupList-method
#' @example inst/examples/plotLGDistances.R
#' @return a heatplot of linkage group calls
#' @export
#' @importFrom cluster daisy
#' @importFrom gplots heatmap.2
#' @importFrom colorspace rainbow_hcl
#' @importFrom reshape2 melt 
####################################################################################################

setMethod('plotLGDistances',
          signature = signature(object='LinkageGroupList', allStrands = 'StrandStateMatrix'),
          definition = plotLGDistances.func
)
