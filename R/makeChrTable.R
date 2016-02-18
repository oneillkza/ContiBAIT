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
	if(verbose){message(paste("-> Creating chromosome table from", bamFile, sep=""))}
	lengthOfContigs <- scanBamHeader(bamFile)[[1]][["text"]]
	bamChr <- sapply(lengthOfContigs[grep("SN:", lengthOfContigs)], "[",1)
	bamChr <- gsub("SN:", "", as.character(bamChr))

	bamLength <- sapply(lengthOfContigs[grep("LN:", lengthOfContigs)], "[",2)
	bamLength <-  as.numeric(gsub("LN:", "", bamLength))


	chrTable <- GRanges(bamChr, IRanges(start=0, end=bamLength))

	if(!(is.null(splitFile)))
	{
		appendedChrTable <- append(chrTable, splitFile)
		chrTable <- disjoin(appendedChrTable)
	}

	if(!(is.null(splitBy)))
	{
		chrTable <- subdivideGRanges(chrTable, subsize=splitBy)
	}

	chrTable$name <- paste(seqnames(chrTable), ":", start(chrTable), "-", end(chrTable), sep="")
	chrTable <- ChrTable(chrTable)
	return(chrTable)
}