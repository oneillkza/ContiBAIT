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
#' @param asBed if TRUE, a bed format table will be generated (chr, start, end)
#' @param verbose if FALSE, no messages appear on terminal
#' @param asRownames Boolean. If false, rownames will not be returned, otherwise rownames will be equal to chromosome names
#'  
#' @details makeChrTable creates a table with chromosome name and chromosome length by extracting 
#' header data from the supplied bam file.
#'
#' 
#' @import Rsamtools
#' @export
#
####################################################################################################


#function to create a chromosome table (chromosome, length) from a user-input bam file.
makeChrTable <- function(bamFile, verbose=TRUE, asBed=FALSE, asRownames=TRUE)
{
	options(scipen=20)
	if(verbose){message(paste("-> Creating chromosome table from", bamFile, sep=""))}
	lengthOfContigs <- scanBamHeader(bamFile)[[1]][["text"]]
	bamChr <- sapply(lengthOfContigs[grep("SN:", lengthOfContigs)], "[",1)
	bamChr <- gsub("SN:", "", as.character(bamChr))

	bamLength <- sapply(lengthOfContigs[grep("LN:", lengthOfContigs)], "[",2)
	bamLength <-  as.numeric(gsub("LN:", "", bamLength))

	if(asBed)
	{
		chrTable <- data.frame(chr=bamChr, start=0, end=bamLength)
	}else{
		chrTable <- data.frame(chr=bamChr, length=bamLength)
	}
	if(asRownames){rownames(chrTable) <- chrTable[,1] }  

	return(new("ChrTable", chrTable))
}