####################################################################################################
#' Plot ordering of contigs within a single linkage group.
#' @param contigOrder vector with the names of the contigs to plot
#' @export
####################################################################################################

plotContigOrder <- function(contigOrder)
{
  contigStarts <- sub('.*:', '', contigOrder)
  contigStarts <- sub('-.*', '', contigStarts)
  contigStarts <- as.numeric(contigStarts)
  plot(contigStarts/10^6, pch='_', 
       ylab='True Contig Start Position (Mb)', xlab='Predicted Contig Order',
       cex=2)
}