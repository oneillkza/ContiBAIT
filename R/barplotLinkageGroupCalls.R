
####################################################################################################
#' Bar plot all linkage groups, with the true chromosomes of contigs coloured.
#' @param linkageGroups list of vectors, each specifying which contigs belong in which linkage group
#' @param assemblyBED table from a BED containing assembly information about the contigs, including length and chromosome
#' Note that the rownames of assemblyBED should be the contig names, as they are used in linkageGroups. 
#' To use a bam file header, the product of makeChrTable(bamFile) with or without the asBed option is suitable for input
#' @param by whether to plot by linkage group (if 'lg') or chromosomes ('chr')
#' @param returnTable TRUE to return chromosome length matrix.
#' @param saveFile string providing a file name to save the plot as a pdf (with legend). Default is no saved file.
#' @param ... any additional parameters to pass down to barplot()
#' Note to include legend, use legend=rownames(chr.table) for by='lg', and legend=colnames(chr.table) for by='chr'
#' 
#' @return a matrix of lengths of each chromosome (rows) in each linkage group (columns)
#' 
#' @export
#' @importFrom colorspace rainbow_hcl
#' @importFrom gtools mixedsort
####################################################################################################

barplotLinkageGroupCalls <- function(linkageGroups, assemblyBED, by='lg', returnTable=FALSE, saveFile=FALSE,  ...)
{

	linkage.chr <- lapply(linkageGroups, function(x){as.character(assemblyBED[x, colnames(assemblyBED[1])])})
	complete.list <- unique(unlist(linkage.chr))

	# If using chromosome notation, order by chromosome
	if( length(grep('chr', complete.list)) == length(complete.list) )
	{ 
		#if name in format chr:start-end, then find location of colon in string
		locationOfColon <- sapply(lapply(strsplit(complete.list, ''), function(x) which(x == ':')), "[[", 1)-1
		complete.list <- complete.list[order(as.numeric(substring(complete.list,4, locationOfColon)))]
	}

	chr.table <- computeBarPlotMatrix(linkageGroups, assemblyBED)
	chr.table <- chr.table[mixedsort(rownames(chr.table)),]

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
		if(saveFile != FALSE){pdf(paste(saveFile, '_barplot_LG.pdf', sep=""))}
		barplot(chr.table, col = chr.cols, names.arg=colnames(chr.table), las=2, xlab='Linkage group', ylab='Length in Mb', ylim=c(0,maxX), ...)

		if(saveFile != FALSE){
			pie(rep(1,length(complete.list)), labels=complete.list, col=chr.cols, clockwise=T)
			dev.off()
		}

	}
	
	#Alternately, plot by chromosome:
	if(by=='chr')
	{
		maxX <- max(sapply(seq(1,nrow(chr.table)), function(x) sum(chr.table[x,])))
		maxX <- roundUpNice(max(maxX))
		linkage.cols <- rainbow_hcl(ncol(chr.table), c=100, l=60) 
		if(saveFile != FALSE){pdf(paste(saveFile, '_barplot_chr.pdf', sep=""))}
		barplot(t(chr.table), col=linkage.cols, names.arg=rownames(chr.table), xlab='Chromosome', ylab='Length in Mb', ylim=c(0,maxX), ...) 
		if(saveFile != FALSE){
			pie(rep(1,ncol(chr.table)), labels=colnames(chr.table), col=linkage.cols, clockwise=T)
			dev.off()
		}

	}
	if(returnTable == TRUE)
	{
		return(chr.table)
	}
}
