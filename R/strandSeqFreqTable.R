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
#' @param filter  additional file of type ChrTable (GRanges with a meta column titled 'name' determining contig name) to split chromosomes based on locations. If this parameter is NULL (default),
#' a filter table will be automatically generated from the header of the first file in bamFileList.
#' @param misorientations additional file of type ChrTable (GRanges with a meta column titled 'name' determining contig name) with putative misorientation locations. These location's
#' reads are reversed in the regions of filter that they lie within. Product of locateMisorients[[1]]. Default is NULL
#' @param chimeras additional file of type ChrTable (GRanges with a meta column titled 'name' determining contig name) with putative chimeric locations. These location's
#' reads are masked in the regions of filter that they lie within, and are appended to the end of the output file. Product of locateMisorients[[2]]. Default is NULL
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
							   chimeras=NULL,
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

	ProcessBam <- function(filt, str, bf) {
		# Count plus strand reads from first read
		result <- reduceByYield(bf, strandInfo, 
								   overlapStrand, 
								   DONE=loopedChunk, 
								   grfilter=filt,
								   strand=str)
		#fixes bug where no reads returns a list rather than an integer
		if(class(result) == 'list'){result <- as.integer(rep(0, length(filter)))}
		return(result)
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
		if(verbose){message(paste(" -> ", lengthOfContigs," fragments found", sep=""))}
	}else{
		if(!(is.null(chimeras)))
		{
			lengthOfContigs <- length(filter)+length(chimeras)
		}else{
			lengthOfContigs <- length(filter)
		}
		if(verbose){message(paste(" -> ", lengthOfContigs," fragments found", sep=""))}
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

	if(!(is.null(chimeras))){rownames(strandTable) <- c(filter$name, chimeras$name)}else{rownames(strandTable) <- filter$name}
	countTable <- strandTable
	
	if(BAITtables == TRUE){
		WatsonTable <- strandTable
		CrickTable <- strandTable
	}


	#if there is strand information in the filter file, use to flip misoriented fragments
	strand(filter) <- "*" 

	if(!(is.null(misorientations)))
	{
		if(verbose){message(' -> Processing ', length(misorientations), ' misorientations')}
		#First, remove strandedness to prevent weird Rsamtools behaviour
		strand(misorientations) <- "*"
		#Now, split the misorientations by wherever there are splits in filter. Prevents a misorientation that covers a portion of
		#two fragments to be excluded. This will chop up the misorientations into chunks that cover each fragment.
		tempRange <- GRanges(append(filter, misorientations))
		tempRange <- ChrTable(disjoin(tempRange))		
		misorientations <- tempRange[queryHits(findOverlaps(tempRange, misorientations))]
		lengthOfNeg <- length(misorientations)
	}

	if(!(is.null(chimeras)))
	{
		if(verbose){message(' -> Processing ', length(chimeras), ' chimeras')}
		#First, remove strandedness to prevent weird Rsamtools behaviour
		strand(chimeras) <- "*"
		#Now, split the chimeras by wherever there are splits in filter. Prevents a misorientation that covers a portion of
		#two fragments to be excluded. This will chop up the chimeras into chunks that cover each fragment.
		tempRange <- GRanges(append(filter, chimeras))
		tempRange <- ChrTable(disjoin(tempRange))		
		chimerasSplit <- tempRange[queryHits(findOverlaps(tempRange, chimeras))]
		#Now prevent bug where misorientation regions with chimeras screw things up
		if(!(is.null(misorientations)))
		{
			tempRange <- GRanges(append(ChrTable(chimerasSplit), misorientations))
			tempRange <- ChrTable(disjoin(tempRange))		
			chimerasSplit <- tempRange[queryHits(findOverlaps(tempRange, chimeras))]
		}
		lengthOfChim <- length(chimerasSplit)
	}

	for(fileName in bamFileList)
	{
		index <- colnames(strandTable)[indexCounter]
		# Read bamfile into tileChunk pieces
		bf = BamFile(fileName, yieldSize=tileChunk)
		# Count plus strand reads from first read

		resultPos <- ProcessBam(filter, TRUE, bf)
		resultNeg <- ProcessBam(filter, FALSE, bf)

		if(!(is.null(misorientations)))
		{
			resultPosMinus <- ProcessBam(misorientations, TRUE, bf)
			resultNegMinus <- ProcessBam(misorientations, FALSE, bf)

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
			resultNeg[as.numeric(names(toSwitchToNeg))] <- resultNeg[as.numeric(names(toSwitchToNeg))]+toSwitchToNeg
		}

		if(!(is.null(chimeras)))
		{

			#calculate dijoined chimeras to mask filter
			resultPosChim <- ProcessBam(chimerasSplit, TRUE, bf)
			resultNegChim <- ProcessBam(chimerasSplit, FALSE, bf)

			if(!(is.null(misorientations)))
			{
				resultPosTemp <- resultPosChim
				hits <- subjectHits(findOverlaps(misorientations, chimerasSplit, type='within'))
				resultPosChim[hits] <- resultNegChim[hits]
				resultNegChim[hits] <- resultPosTemp[hits]
			}

			#and calculate undisjoined chimeras to append to results
			chimResultPos <- ProcessBam(chimeras, TRUE, bf)
			chimResultNeg <- ProcessBam(chimeras, FALSE, bf)

			#now find the overlaps, and assign each misorientation a name that corresponds to the location of the fragment in filter
			names(resultPosChim) <- subjectHits(findOverlaps(chimerasSplit,filter, type="within"))
			names(resultNegChim) <- names(resultPosChim)

			#condense regions where >1 misorientations occur in the same fragment, by summing the values that have the same name
			toRemovePosChim <- tapply(resultPosChim, names(resultPosChim), sum)
			toRemoveNegChim <- tapply(resultNegChim, names(resultNegChim), sum)

			#Remove reads from chimeric regions
			resultPos[as.numeric(names(toRemovePosChim))] <- resultPos[as.numeric(names(toRemovePosChim))]-toRemovePosChim
			resultNeg[as.numeric(names(toRemoveNegChim))] <- resultNeg[as.numeric(names(toRemoveNegChim))]-toRemoveNegChim
			resultPos <- c(resultPos, chimResultPos)
			resultNeg <- c(resultNeg, chimResultNeg)

			resultPos[which(resultPos < 0)] <- 0
			resultNeg[which(resultNeg < 0)] <- 0
		}


		# Total read number
		absCount <- resultPos + resultNeg

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