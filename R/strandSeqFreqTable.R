# Copyright (c) 2015, Mark Hills
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


####################################################################################################
#' strandSeqFreqTable -- function to process bam files for contiBAIT 
#' 
#' @param bamFileList  vector containing the location of the bams file to be read
#' @param field  The field of the bam file name to use as an index (default is 1)
#' @param fieldSep  The field seperator of the bam file to use to define the field. Default is '.'
#' @param qual  Mapping quality threshold. Default is 0
#' @param rmdup  remove duplicates in output file. Default is TRUE 
#' @param filter  additional file of type ChrTable (GRanges with a meta column titled 'name' determining contig name) to split chromosomes based on locations. If this parameter is blank,
#' a filter table will be automatically generated from the header of the first file in bamFileList.
#' @param misorientations additional file of type ChrTable (GRanges with a meta column titled 'name' determining contig name) with putative misorientation locations. These location's
#' reads are reversed in the output
#' @param tileChunk  Number of reads to split bam files into (smaller number requires less RAM). Default is 100000.
#' @param pairedEnd  Whether the bam files being read are in paired end format. Default is TRUE. Note,
#' since paired reads will be the same direction, only first mate read of pair is used in output
#' @param verbose prints messages to the terminal (default is TRUE)
#' @param BAITtables creates additional matrices in the returned list with just Watson and Crick read counts to be used in downstreat BAIT plotting. Default is FALSE
#' 
#' @return a list containing two matrices: a StrandFreqMatrix of W:C read frequencies, and a StrandReadMatrix of read counts
#' @import Rsamtools
#' @import IRanges
#' @import GenomicFiles
#' @import GenomicRanges
#' @importFrom S4Vectors DataFrame
#' @example inst/examples/strandSeqFreqTable.R
#' @export
#' @include AllClasses.R
####################################################################################################


strandSeqFreqTable <- function(bamFileList, 
							   fieldSep='.', 
							   field=1, 
							   qual=0, 
							   rmdup=TRUE, 
							   verbose=TRUE, 
							   filter=NULL, 
							   misorientations=NULL,
							   tileChunk=100000, 
							   pairedEnd=TRUE,
							   BAITtables=FALSE)
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
		    param <- ScanBamParam(flag=scanBamFlag(isFirstMateRead=TRUE, 
		    									   isUnmappedQuery=FALSE, 
		    									   isMinusStrand=strand ), 
		    					  mapqFilter=qual, 
		    					  what=c("rname","pos","strand"))
	    }else{
	    	# Handle single end reads by just looking at strand direction	
		    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, 
		    									   isMinusStrand=strand ), 
		    					  mapqFilter=qual, what=c("rname","pos","strand"))
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
	if(is.null(filter))
	{
		filter <- makeChrTable(bamFileList[1])
		lengthOfContigs <- length(filter)
		if(verbose){message(paste("-> ", lengthOfContigs," fragments found", sep=""))}
	}else{
		lengthOfContigs <- length(filter)
		if(verbose){message(paste("-> ", lengthOfContigs," fragments found", sep=""))}
	}

	strandTable <- matrix(nrow=lengthOfContigs, ncol=bamFileLength)
	if(field != FALSE)
	{
		colnames(strandTable) <- sapply(bamFileList, 
										function(x) 
											strsplit(basename(x), 
													paste('\\', fieldSep, sep=""))[[1]][field] )
	}else{
		colnames(strandTable) <- bamFileList 
	}

	#if library names are numeric, add 'lib_' in front to avoid numeric colnames 
	colnames(strandTable) [grep("^[0-9]", colnames(strandTable) )] <- 
		paste('lib', colnames(strandTable) [grep("[0-9]", colnames(strandTable) )], sep='_')

	rownames(strandTable) <- filter$name
	countTable <- strandTable
	
	if(BAITtables == TRUE){
		WatsonTable <- strandTable
		CrickTable <- strandTable
	}


	#if there is strand information in the filter file, use to flip misoriented fragments
	strand(filter) <- "*" 

	if(!(is.null(misorientations)))
	{
		#First, remove strandedness to prevent weird Rsamtools behaviour
		strand(misorientations) <- "*"
		#remove the 'chrUn' fragments, if present
		tempRange <- filter[-(grep("chrUn", filter$name))]
		#Now, split the misorientations by wherever there are splits in filter. Prevents a misorientation that covers a portion of
		#two fragments to be excluded. This will chop up the misorientations into chunks that cover each fragment.
		tempRange <- GRanges(append(tempRange, misorientations))
		tempRange <- ChrTable(disjoin(tempRange))		
		misorientations <- tempRange[queryHits(findOverlaps(tempRange, misorientations))]
		lengthOfNeg <- length(misorientations)
	}

	for(fileName in bamFileList)
	{
		index <- colnames(strandTable)[indexCounter]
		# Read bamfile into tileChunk pieces
		bf = BamFile(fileName, yieldSize=tileChunk)
		# Count plus strand reads from first read
		resultPos <- reduceByYield(bf, strandInfo, 
								   overlapStrand, DONE=loopedChunk, 
								   grfilter=filter)
		#fixes bug where no reads returns a list rather than an integer
		if(class(resultPos) == 'list'){resultPos <- as.integer(rep(0, lengthOfContigs))}

		if(is.list(resultPos)){
			warning('\n####################\n WARNING! BAM FILE', 
					index, 'APPEARS TO BE SINGLE-END. TRY RERUNNING WITH pairedEnd=FALSE \n####################')
			break
		}
		# Count minus strand reads from first read
		resultNeg <- reduceByYield(bf, 
								   strandInfo, 
								   overlapStrand, 
								   DONE=loopedChunk, 
								   grfilter=filter, 
								   strand=FALSE)

		if(class(resultNeg) == 'list'){resultNeg <- as.integer(rep(0, lengthOfContigs))}
		# Total read number
		absCount <- resultPos + resultNeg
		# Calculate strand call

		if(!(is.null(misorientations)))
		{
			resultPosMinus <- reduceByYield(bf, strandInfo, 
										   overlapStrand, DONE=loopedChunk, 
										   grfilter=misorientations)
			if(class(resultPosMinus) == 'list'){resultPosMinus <- as.integer(rep(0, lengthOfNeg))}

			resultNegMinus <- reduceByYield(bf, 
										   strandInfo, 
										   overlapStrand, 
										   DONE=loopedChunk, 
										   grfilter=misorientations, 
										   strand=FALSE)
			if(class(resultNeg) == 'list'){resultNegMinus <- as.integer(rep(0, lengthOfNeg))}

			#now find the overlaps, and assign each misorientation a name that corresponds to the location of the fragment in filter
			names(resultPosMinus) <- subjectHits(findOverlaps(misorientations,filter, type="within"))
			names(resultNegMinus) <- names(resultPosMinus)

			#condense regions where >1 misorientations occur in the same fragment, by summing the values that have the same name
			toSwitchToNeg <- tapply(resultPosMinus, names(resultPosMinus), sum)
			toSwitchToPos <- tapply(resultNegMinus, names(resultNegMinus), sum)

			#Switch the minus reads to plus for the internal misorients
			resultPos[as.numeric(names(toSwitchToPos))] <- resultPos[as.numeric(names(toSwitchToPos))]+toSwitchToPos
			resultNeg[as.numeric(names(toSwitchToPos))] <- resultNeg[as.numeric(names(toSwitchToPos))]-toSwitchToPos
			#And switch the plus reads to minus for internal misorients
			resultPos[as.numeric(names(toSwitchToNeg))] <- resultPos[as.numeric(names(toSwitchToNeg))]-toSwitchToNeg
			resultNeg[as.numeric(names(toSwitchToNeg))] <- resultNeg[as.numeric(names(toSwitchToNeg))]+toSwitchToPos
		}

		strandCall <- (resultPos - resultNeg)/absCount
			
		strandTable[,indexCounter] <- strandCall 
		countTable[,indexCounter] <- absCount


		if(BAITtables == TRUE){
			WatsonTable[,indexCounter] <- resultPos
			CrickTable[,indexCounter] <- resultNeg
		}

		if(verbose){message('-> Creating contig table for index ',
							index, " [", indexCounter, "/", bamFileLength, "]")}
		indexCounter <- indexCounter+1
	}

	#Get rid of contigs that are entirely empty to prevent low quality calls in downstream functions
	strandTable <- strandTable[which(apply(countTable, 1, sum) > 0),]
	countTable <- countTable[which(apply(countTable, 1, sum) > 0),]

	if(BAITtables == FALSE){
  		return(list(strandTable=StrandFreqMatrix(strandTable), 
  					countTable=StrandReadMatrix(countTable)))
  	}else{
  		return(list(strandTable=StrandFreqMatrix(strandTable), 
  					countTable=StrandReadMatrix(countTable), 
  					WatsonReads=StrandReadMatrix(WatsonTable), 
  					CrickReads=StrandReadMatrix(CrickTable)))

  	}
}	