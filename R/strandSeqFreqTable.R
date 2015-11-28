# Copyright (c) 2015, Mark Hills
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


####################################################################################################
#' strandSeqFreqTable -- function to process bam files for contiBAIT 
#' 
#' @param bamFileIst  vector containing the location of the bams file to be read
#' @param filed  The field of the bam file name to use as an index (default is 1)
#' @param fieldSep  The field seperator of the bam file to use to define the field. Default is '.'
#' @param qual  Mapping quality threshold. Default is 0
#' @param rmdup  remove duplicates in output file. Default is TRUE 
#' @param filter  additional file to split chromosomes based on locations. If this parameter is blank,
#' a filter table will be automatically generated from the header of the first file in bamFileList
#' @param tileChunk  Number of reads to split bam files into (smaller number requires less RAM). Default is 100000.
#' @param pairedEnd  Whether the bam files being read are in paired end format. Default is TRUE. Note,
#' since paired reads will be the same direction, only first mate read of pair is used in output
#' @param verbose prints messages to the terminal (default is TRUE)
#' 
#' @return all reads in bed format if frequencies=FALSE, or list of frequency of strand calls for each contig and number of reads if frequencies=TRUE
#' @import Rsamtools
#' @import GenomicFiles
#' @importFrom S4Vectors DataFrame
#' @export
#' @include AllClasses.R
####################################################################################################

strandSeqFreqTable <- function(bamFileList, fieldSep='.', field=1, qual=0, rmdup=TRUE, verbose=TRUE, filter=FALSE, tileChunk=100000, pairedEnd=TRUE)
{
	##### DEFINE FUNCTIONS

	#Function to tidy up the matrices, convert first row (index) to the column name
	#And remove convert value from strandTable to NA if number of reads is less than countLimit
	limitStrandCounts <- function(strandTable) {
		colnames(strandTable) <- strandTable[1,]
		strandTable <- strandTable[-1,]
		mode(strandTable) <- "numeric"
   		return(strandTable)
	}

	strandInfo <- function(x, strand=TRUE, ...) {
		# Handle paired end read by looking only at first mate read
		if(pairedEnd){
		    param <- ScanBamParam(flag=scanBamFlag(isFirstMateRead=TRUE, isUnmappedQuery=FALSE, isMinusStrand=strand ), mapqFilter=qual, what=c("rname","pos","strand"))
	    }else{
	    	# Handle single end reads by just looking at strand direction	
		    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, isMinusStrand=strand ), mapqFilter=qual, what=c("rname","pos","strand"))
	    }
	    scanBam(x, param=param)[[1]]
	}

	overlapStrand <- function(x, ..., grfilter) {
	    ## 'x' is the return value from strandInfo
	    grinput <- with(x, GRanges(rname, IRanges(pos, width=1), strand))
	    if(rmdup){grinput <- grinput[!duplicated(grinput)]}
	    countOverlaps(grfilter, grinput)
	}

	loopedChunk <- function(x, ...) {
	    length(x[["rname"]]) == 0L
	}

	##### PROCESS DATA

	bamFileLength <- length(bamFileList)
	indexCounter <- 1
	#create empty matrices with the same number of columns as contigs

	if(verbose){message(paste("-> Counting number of fragments", sep=""))}
	if(length(filter) == 1)
	{
		filter <- makeChrTable(bamFileList[1], asBed=T)
		filter$name <- filter$chr
		lengthOfContigs <- nrow(filter)
		if(verbose){message(paste("-> ", lengthOfContigs," fragments found", sep=""))}
	}else{
		if(ncol(filter) == 3){filter$name <- paste(filter[,1], ':', filter[,2], '-', filter[,3], sep='')}
		lengthOfContigs <- nrow(filter)
		if(verbose){message(paste("-> ", lengthOfContigs," fragments found", sep=""))}
		colnames(filter) <- c('chr', 'start', 'end', 'name')
	}

	strandTable <- matrix(nrow=lengthOfContigs, ncol=bamFileLength)
	colnames(strandTable) <- sapply(bamFileList, function(x) strsplit(basename(x), paste('\\', fieldSep, sep=""))[[1]][field] )
	rownames(strandTable) <- filter[,4]
	countTable <- strandTable

	for(fileName in bamFileList)
	{
		# Find the index
		index <- strsplit(basename(fileName), paste('\\', fieldSep, sep=""))[[1]][field]

		# Make GRanges object from filter
		grfilter <- makeGRangesFromDataFrame(filter)
		# Read bamfile into tileChunk pieces
		bf = BamFile(fileName, yieldSize=tileChunk)
		# Count plus strand reads from first read
		resultPos <- reduceByYield(bf, strandInfo, overlapStrand, DONE=loopedChunk, grfilter=grfilter)
		# Count minus strand reads from first read
		resultNeg <- reduceByYield(bf, strandInfo, overlapStrand, DONE=loopedChunk, grfilter=grfilter, strand=FALSE)
		# Total read number
		absCount <- resultPos + resultNeg
		# Calculate strand call
		strandCall <- (resultPos - resultNeg)/absCount
			
		strandTable[,indexCounter] <- strandCall 
		countTable[,indexCounter] <- absCount

		if(verbose){message(paste('-> Creating contig table for index ', index, " [", indexCounter, "/", bamFileLength, "]      ", sep=""), "\r", appendLF=FALSE)}
		indexCounter <- indexCounter+1
	}
	if(verbose){message(" ")}

  return(list(strandTable, countTable))
}