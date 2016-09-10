# Compute similarity between contig and linkage group
# @useDynLib contiBAIT
# @import Rcpp
# @import BH
# @keywords internal 

computeSim <- function(contigStrand, linkageStrand, minimumLibraryOverlap)
{
	numCommon <- sum(!is.na(contigStrand) & !is.na(linkageStrand))
	numEqual <- sum((contigStrand == 2) == (linkageStrand == 2), na.rm=TRUE)
	if(numCommon < minimumLibraryOverlap) NA else numEqual / numCommon
}