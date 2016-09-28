
writeBed.func <- function(chrTable, 
					 	  orientationData, 
					 	  contigOrder,
					 	  #misorientations=NULL,
					 	  #chimeras=NULL,
					 	  libWeight=NULL,
					 	  file='contiBAIT_assembly.bed')
{

	chrTable <- as.data.frame(chrTable)
	rownames(chrTable) <- chrTable$name
	chrTable <- chrTable[as.character(contigOrder[,2]),1:3]
	if(is.null(libWeight)){libWeight <- rep(0, nrow(chrTable))}else{libWeight <- libWeight[contigOrder[,2]]}
	bedRange <- GRanges(seqnames=chrTable[,1],
						IRanges(start=chrTable[,2], end=chrTable[,3]),
						strand=as.character(strand(orientationData[match(as.character(contigOrder[,2]), orientationData$name)])),
						score=libWeight,
						name=paste(contigOrder[,1],contigOrder[,2], sep="_"))
	export.bed(con =file,bedRange)

}

####################################################################################################
#' function to write contig order to BED file
#' @param file character string for bed file name to write
#' @param chrTable a GRanges object with a 'name' meta column matching contig names. Product of makeChrTable
#' @param orientationData ChrTable of contig and strand (with rownames matching contig names). Product of reorientAndMergeLGs[[2]]
#' @param contigOrder an object of type ContigOrdering with ordered Linkage Groups and contigs. Product of orderAllLinkageGroups 
#' @param libWeight average quality across all libraries for a contig
#' @importFrom rtracklayer export.bed
#' @import GenomicRanges
#' @aliases writeBed writeBed-ChrTable-ContigOrdering-method
#' @rdname writeBed
#' @example inst/examples/writeBed.R
#' @return NULL; BED file written to file
#' 
#' @export
####################################################################################################

setMethod('writeBed',
		  signature = signature(orientationData='ChrTable', 
		  					  contigOrder='ContigOrdering'),
		  definition = writeBed.func
		  )
