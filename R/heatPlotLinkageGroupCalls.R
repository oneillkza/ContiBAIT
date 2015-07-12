####################################################################################################
#' Heat plot linkage groups against true chromosomes, by lengths in MB
#' @param linkageGroups list of vectors, each specifying which contigs belong in which linkage group
#' @param assemblyBED table from a BED containing assembly information about the contigs, including length and chromosome
#' Note that the rownames of assemblyBED should be the contig names, as they are used in linkageGroups
#' @param ... any additional parameters to pass down to heatmap.2()
#' @export
#' @importFrom gplots heatmap.2
####################################################################################################
heatPlotLinkageGroupCalls <- function(linkageGroups, assemblyBED, ...)
{
	chr.table <- computeBarPlotMatrix(linkageGroups, assemblyBED)
	
	heatmap.2(chr.table, dendrogram='none', trace='none', 
			  col=c('gray',rev(heat.colors(500))), 
			  xlab='Linkage Groups', 
			  ylab='Chromosomes')
}