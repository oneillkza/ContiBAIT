# Copyright (c) 2014, Mark Hills & Kieran O'Neill
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

####################################################################################################
#' makeChrTable -- Pulls out chromosome and length data from the header of a bam file
#' 
#' @param bamFile string of location of a bam file to extract header data from
#' @param splitFile GRanges object (of type chr, start and end: no strand or meta columns) of locations in which to split
#' the assembly, such as previously determined locations of contig chimerism
#' @param splitBy integer determining the average size contigs should be split by
#' @param verbose if FALSE, no messages appear on terminal
#'  
#' @details makeChrTable creates a table with chromosome name and chromosome length by extracting 
#' header data from the supplied bam file.
#' @example inst/examples/makeChrTable.R
#' @return a GRanges object of class ChrTable, containing information on the organism's chromosomes as extracted from the BAM file header.
#' 
#' @import Rsamtools
#' @import GenomicRanges
#' @importFrom exomeCopy subdivideGRanges
#' @importFrom stats aggregate rbinom runif
#' @include AllClasses.R
#' @export
#
####################################################################################################


#function to create a chromosome table (chromosome, length) from a user-input bam file.
makeChrTable <- function(bamFile,
						 splitFile=NULL,
						 splitBy=NULL, 
					     verbose=TRUE) 
{
	if(verbose){message(paste("-> Creating chromosome table from ", bamFile, sep=""))}
	lengthOfContigs <- scanBamHeader(bamFile)[[1]][["text"]]
	bamChr <- sapply(lengthOfContigs[grep("SN:", lengthOfContigs)], "[",1)
	bamChr <- gsub("SN:", "", as.character(bamChr))

	bamLength <- sapply(lengthOfContigs[grep("LN:", lengthOfContigs)], "[",2)
	bamLength <-  as.numeric(gsub("LN:", "", bamLength))

	chrTable <- GRanges(bamChr, IRanges(start=1, end=bamLength))

	if(!(is.null(splitFile)))
	{
		#If the splitFile contains a name meta (ie is a ChrTable instance), remove it to allow appending.
		if(length(elementMetadata(splitFile)) != 0)
		{
			strippedSplitFile <- splitFile
			strand(strippedSplitFile) <- '*'
			elementMetadata(strippedSplitFile) <- NULL
			appendedChrTable <- append(chrTable, strippedSplitFile)			
		}else{ 
			appendedChrTable <- append(chrTable, splitFile)
		}
		chrTable <- disjoin(appendedChrTable)
	}

	if(!(is.null(splitBy)))
	{
		#Only split those fragments bigger than splitBy; improves speed of subdivideGRanges
		splitBig <- chrTable[which(width(chrTable) >= splitBy)]
		dontSplit <- chrTable[which(width(chrTable) < splitBy)]
		#Alternate strand for fragments. Prevents inherent 'reduce' within subdivideGRanges that will merge fragments previously split with splitFile
		strand(splitBig) <- rep(c('+','-'), (length(splitBig)/2)+1)[1:length(splitBig)]
		splitBig <- subdivideGRanges(splitBig, subsize=splitBy)
		strand(splitBig) <- "*"
		#then merge everything back together
		chrTable <- append(splitBig, dontSplit)
		chrTable <- sort(chrTable)
	}

	#now attempt to reappend metadata. First if there's strand information
	olap <- NULL
	if(length(which(strand(splitFile) != "*")) != 0)
	{
		olap <- findOverlaps(chrTable, splitFile)
		strand(chrTable[queryHits(olap)]) <- strand(splitFile[subjectHits(olap)])
	}

	if(length(splitFile$name) > 0)
	{
		if(is.null(olap)){olap <- findOverlaps(chrTable, splitFile)}
		chrTable$name <- c(1:length(chrTable))
		chrTable$name[queryHits(olap)] <- splitFile$name[subjectHits(olap)]
		unLocation <- grep("chrUn", chrTable$name)
		chrTable$name <- paste(seqnames(chrTable), ":", start(chrTable), "-", end(chrTable), sep="")	
		chrTable$name[unLocation] <- paste("chrUn_", seqnames(chrTable)[unLocation], ":", start(chrTable)[unLocation], "-", end(chrTable)[unLocation], sep="")
	}else{
		chrTable$name <- paste(seqnames(chrTable), ":", start(chrTable), "-", end(chrTable), sep="")	
	}

	#removes bug where occasionally empty or 1 nucleotide splits occur
	chrTable <- chrTable[width(chrTable) > 1]
	chrTable <- ChrTable(chrTable)
	return(chrTable)
}