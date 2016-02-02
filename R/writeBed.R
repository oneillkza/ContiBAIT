
writeBed.func <- function(chrTable, 
					 	  orientationData, 
					 	  contigOrder,
					 	  libWeight=NULL,
					 	  file='contiBAIT_assembly.bed')
{

	if(is.null(libWeight))
	{
		libWeight <- rep(0, length(contigOrder[,2]))
		names(libWeight) <- contigOrder[,2]
	}

	chrTable <- as.data.frame(chrTable)
	rownames(chrTable) <- chrTable$name
	chrTable <- chrTable[,1:3]
	toBed <- data.frame(chrTable[contigOrder[,2],], strand=orientationData[contigOrder[,2],2] )
	bedNames <- paste(contigOrder[,1],contigOrder[,2], sep="_")
	bedRange <- GRanges(toBed[,1:3], strand=toBed$strand, score=libWeight[contigOrder[,2]], name=bedNames)
	export.bed(con =file,bedRange)

}

####################################################################################################
#' function to write contig order to BED file
#' @param file character string for bed file name to write
#' @param chrTable a GRanges object with a 'name' meta column matching contig names. Product of makeChrTable
#' @param orientationData data.frame of contig and strand (with rownames matching contig names). Product of reorientLinkageGroups[[2]]
#' @param contigOrder an object of type ContigOrdering with ordered Linkage Groups and contigs. Product of orderAllLinkageGroups 
#' @param libWeight average quality across all libraries for a contig
#' @importFrom rtracklayer export.bed
#' @import GenomicRanges
#' @aliases writeBed writeBed-OrientationFrame-ContigOrdering-method
#' @rdname writeBed
#' @example inst/examples/writeBed.R
#' @return NULL; BED file written to file
#' 
#' @export
####################################################################################################

setMethod('writeBed',
		  signature = signature(orientationData='OrientationFrame', contigOrder='ContigOrdering'),
		  definition = writeBed.func
		  )
