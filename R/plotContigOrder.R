####################################################################################################
#' Plot ordering of contigs within a single linkage group.
#' @param contigOrder vector with the names of the contigs to plot
#' @export
####################################################################################################

plotContigOrder <- function(contigOrder)
{

  contigChr <- sub(':.*', '', contigOrder)
  primaryContigChr <- names(sort(table(contigChr), decreasing=T))[1]

  contigStarts <- sub('.*:', '', contigOrder)
  contigStarts <- sub('-.*', '', contigStarts)
  contigStarts <- as.numeric(contigStarts)
  contigStarts[which(!( sub(':.*', '', contigOrder) == primaryContigChr))] <- 0
  plot(contigStarts/10^6, pch='_', 
       ylab='True Contig Start Position (Mb)', xlab='Predicted Contig Order',
       cex=2)
}