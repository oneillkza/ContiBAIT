# Compute similarity between contig and linkage group
# @useDynLib contiBAIT
# @import Rcpp
# @keywords internal 

computeSim <- function(contigStrand, linkageStrand, minimumLibraryOverlap)
{
	.Call('computeSim', as.integer(contigStrand), as.integer(linkageStrand), minimumLibraryOverlap)
}