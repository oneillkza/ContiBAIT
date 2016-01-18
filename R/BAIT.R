# Copyright (c) 2016, Mark Hills
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


####################################################################################################
#' BAIT -- master function to process strand-seq libraries into BAIT ideograms 
#' 
#' @param path  String denoting location of Strand-seq bam files (default is ".")
#' @param readQual Integer dictating the minimal mapping quality required for a read to be accepted. Default is 10. 
#' @param pairedEnd  Whether the bam files being read are in paired end format. Default is TRUE. Note,
#' @param plotBy Whether to plot by library ('lib') or chromosome ('chr')
#' @param plotName character which determines file name to be saved. Default is to open an R plot from the terminal
#' @param chroms vector of chromosome number to prevent contig plotting. eg for humans use 1:24. Default is 'all'
#' @param verbose prints messages to the terminal (default is TRUE)
#' 
#' @return ordered contigs in bed format. Depending on options, intermediate files and plots will also be generated
#' @import snow
#' @importFrom S4Vectors DataFrame
#' @example inst/examples/contiBAIT.R
#' @export
#' @include AllClasses.R
####################################################################################################


BAIT <- function(path=".", 
                splitBy=200000,
                readQual=10,
                pairedEnd=TRUE,
                plotBy='lib',
                plotName=NULL,
                chroms='all',
                verbose=TRUE)
{

	if(!(is.null(plotName)))
	{
		pdf(paste(plotName, "_", plotBy, ".pdf", sep=""))	
	}

	bamFileList <- list.files(path=path, pattern=".bam$", full.names=TRUE)

	chrTable <- makeChrTable(bamFileList[1])
	if(chroms[1] != 'all')
	{
		chrTable <- chrTable[chroms,]
	}


	dividedChr <- divideMyChr(chrTable, splitBy=splitBy)

	strandFrequencyList <- strandSeqFreqTable(bamFileList, filter=dividedChr, qual=readQual, pairedEnd=pairedEnd, BAITtables=TRUE)

	ideogramPlot(strandFrequencyList[[3]], strandFrequencyList[[4]], dividedChr, plotBy=plotBy, verbose=TRUE)

	if(!(is.null(plotName))){dev.off()}

}