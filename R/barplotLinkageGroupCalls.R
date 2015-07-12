####################################################################################################
#' Bar plot all linkage groups, with the true chromosomes of contigs coloured.
#' @param linkageGroups list of vectors, each specifying which contigs belong in which linkage group
#' @param assemblyBED table from a BED containing assembly information about the contigs, including length and chromosome
#' Note that the rownames of assemblyBED should be the contig names, as they are used in linkageGroups. 
#' To use a bam file header, the product of makeChrTable(bamFile) with or without the asBed option is suitable for input
#' contig within that linkage group
#' @param by ='lg' whether to plot by linkage group (if 'lg') or chromosomes ('chr')
#' @param returnTable TRUE to return chromosome length matrix
#' @param ... any additional parameters to pass down to barplot()
#' @return a matrix of lengths of each chromosome (rows) in each linkage group (columns)
#' @export
#' @importFrom colorspace rainbow_hcl
####################################################################################################



barplotLinkageGroupCalls <- function(linkageGroups, assemblyBED, by='lg', returnTable=FALSE,  ...)
{

	linkage.chr <- lapply(linkageGroups, function(x){as.character(assemblyBED[x, colnames(assemblyBED[1])])})
	complete.list <- unique(unlist(linkage.chr))
		
	chr.table <- computeBarPlotMatrix(linkageGroups, assemblyBED)
	#colnames(chr.table) <- seq(1:ncol(chr.table))

	roundUpNice <- function(x, nice=c(1,2,3,4,5,6,7,8,9,10)) 
	{
    	if(length(x) != 1) stop("'x' must be of length 1")
    	10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
	}
	  
	#Plot by linkage group
	if(by=='lg')
	{	
		maxX <- max(sapply(seq(1,ncol(chr.table)), function(x) sum(chr.table[,x])))
		maxX <- roundUpNice(max(maxX))
		chr.cols <- rainbow_hcl(length(complete.list), c=90, l=60)
		barplot(chr.table, col = chr.cols, names.arg=colnames(chr.table), las=2, xlab='Linkage group', ylab='Length in Mb', ylim=c(0,maxX), ...)
	}
	
	#Alternately, plot by chromosome:
	if(by=='chr')
	{
		#change order to chromosome number
		chr.table <- chr.table[order(rownames(chr.table)),]
		maxX <- max(sapply(seq(1,nrow(chr.table)), function(x) sum(chr.table[x,])))
		maxX <- roundUpNice(max(maxX))

		linkage.cols <- rainbow_hcl(ncol(chr.table), c=100, l=60) 
		barplot(t(chr.table), col=linkage.cols, names.arg=rownames(chr.table), xlab='Chromosome', ylab='Length in Mb', ylim=c(0,maxX), ...) 
	}
	if(returnTable == TRUE)
	{
		return(chr.table)
	}
}
