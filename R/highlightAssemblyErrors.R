####################################################################################################
#' highlightAssemblyErrors -- Master function to identify misorientations and chimeras in the assembly 
#' 
#' @param path String denoting location of Strand-seq bam files
#' @param splitBy integer determining the average size contigs should be split by
#' @param cluster Number of times to recluster and take the consensus of. If NULL, clustering is 
#' run only once.
#' @param clusterParam optional \code{BiocParallelParam} specifying cluster to use for parallel execution.
#' When \code{NULL}, execution will be serial.
#' @param pairedEnd  Whether the bam files being read are in paired end format. Default is TRUE. Note,
#' @param qual  Mapping quality threshold. Default is 1
#' @param gapFile A GRanges object consisting of start and end locations of assembly gaps (defaul it NULL)
#' @param writeBed Character vector, this option will write the resulting bed file to a specified location with the character as the file name. Defulat is NULL
#' @param verbose prints messages to the terminal (default is TRUE)
#' 
#' @return a directional ChrTable object that can be used in downstream functions (strandSeqFreqTable)
#' @aliases highlightAssemblyErrors highlightAssemblyErrors,highlightAssemblyErrors-GRanges-method
#' @rdname highlightAssemblyErrors
#' @import GenomicRanges
#' @importFrom S4Vectors DataFrame
#' @export
#' @include AllClasses.R
####################################################################################################


highlightAssemblyErrors <- function(path, 
									splitBy=1000000,
									cluster=1,
									clusterParam=NULL,
									pairedEnd=TRUE,
									qual=10,
									gapFile=NULL,
									writeBed=NULL,
									verbose=TRUE)
{
	bamFileList <- list.files(path=path, pattern=".bam$", full.names=TRUE)
	filter <- makeChrTable(bamFileList[1], splitBy=splitBy, verbose=verbose)

	#if(length(which(width(filter) >= splitBy)))
	strandFrequencyList <- strandSeqFreqTable(bamFileList, 
	                                             filter=filter,
	                                             field=FALSE, 
	                                             qual=qual, 
	                                             pairedEnd=pairedEnd,
	                                             verbose=verbose)
	strandFrequencyList[[1]][which(strandFrequencyList[[2]] < 10)] <- NA 

	strandStateMatrixList <- preprocessStrandTable(strandFrequencyList[[1]], 
	                                                filterThreshold=0.95,
	                                                 lowQualThreshold=NULL,
	                                                 verbose=verbose)

	chrGrange <- filter[which(filter$name %in% rownames(strandStateMatrixList[[1]]))]
	chrLength <- unique(seqnames(chrGrange)[which(duplicated(as.character(seqnames(chrGrange))))])
	relatedLibList <- lapply(seq_along(chrLength), function(x) findSimilarLibraries(strandStateMatrixList[[1]], strandFrequencyList[[2]], chrGrange, x, cluster=cluster, clusterParam=clusterParam, verbose=verbose))
	relatedLibList <- LibraryGroupList(relatedLibList[!sapply(relatedLibList, is.null)])

	dirRange <- thoroughBed(bamFileList, relatedLibList, verbose=verbose)

	reorientFile <- locateMisorients(dirRange, gapFile=gapFile, writeBed=writeBed, verbose=verbose)

	filter <- makeChrTable(bamFileList[1], splitFile=reorientFile, splitBy=splitBy, verbose=verbose)
	return(filter)
}
