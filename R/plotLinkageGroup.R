####################################################################################################
#' Plot a heatmap of contigs within a single linkage group.
#' @param groupMembers vector with the names of the contigs to plot
#' @param allStrands matrix of all contigs with strand inheritance for all libraries
#' @param lgChr =NULL vector of true calls for chromosomes of contigs, will be plotted as side colours if 
#' not NULL
#' @param labels =TRUE if TRUE, contig names will be plotted on the axes
#' @param ... additional parameters to pass to heatmap.2
#' @example inst/examples/plotLinkageGroup.R
#' @export
#' @importFrom gplots heatmap.2
#' @importFrom colorspace rainbow_hcl
#' @importFrom cluster daisy
####################################################################################################

plotLinkageGroup <- function(groupMembers, allStrands, lgChr=NULL, labels=TRUE, ...)
{
	groupStrands <- allStrands[groupMembers,]
	groupDist <- 1-as.matrix(daisy(data.frame(groupStrands)))
	
  #If chromosome assignments provided, plot them as sidecolours
  if(is.null(lgChr))
  {
    rowCols <- NULL
    colCols <- NULL
  } else
  {
    chrLabels <- unique(lgChr)
    chrCols <- rainbow_hcl(length(chrLabels), c=90, l=60)
    names(chrCols) <- chrLabels
    rowCols <- chrCols[lgChr]
    colCols <- rowCols
  }
  
  #Omit labels if labels is set to false
  if(labels)
  {
    labRow=NULL
    labCol=NULL
  }else
  {
    labRow=rep('', length(groupMembers))
    labCol=rep('', length(groupMembers))
  }
  
  breaks <- seq(0, 100, length.out=101)/100 
  cols <- colorRampPalette(c("cyan", "blue", "grey30", "black", "grey30", "red", "orange"))
  if(is.null(lgChr))
  {
    heatmap.2(groupDist, trace='none', col=cols(100), breaks=breaks, labRow=labRow, labCol=labCol, ...)
  } else
  {
	 heatmap.2(groupDist, trace='none', col=cols(100), labRow=labRow, breaks=breaks, labCol=labCol, RowSideColors=rowCols, ColSideColors=colCols,  ...)
  }
}