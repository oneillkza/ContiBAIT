preprocessStrandTable.func <- function(strandTable, 
									   strandTableThreshold=0.8, 
									   filterThreshold=0.8, 
									   orderMethod='libsAndConc', 
									   lowQualThreshold=0.9, 
									   verbose=TRUE, 
									   minLib=10)
{
	strandTableLength <- nrow(strandTable)

	# Filter low quality libraries.  Scan "WW" and "CC" background regions.
	lowQualList <- data.frame(library=vector(), quality=vector())
	qualList <- lowQualList

	if(!(is.null(lowQualThreshold)))
	{
		if(verbose){message("-> Checking for high quality libraries")}

		for( col in seq_len(ncol(strandTable)) )
		{
			libName <- colnames(strandTable, do.NULL=FALSE)[col]
			backGroundC <- abs(mean(strandTable[,col][which(strandTable[,col] < -0.6)], na.rm=TRUE)) 
			backGroundW <- abs(mean(strandTable[,col][which(strandTable[,col] > 0.6)], na.rm=TRUE))
			libraryQual <- round(backGroundC + backGroundW / 2, digits=3)
			colQual <- data.frame(library=libName, quality=libraryQual)
			
			if(libraryQual < lowQualThreshold || libraryQual == "NaN")
			{
				if(libraryQual == "NaN" & verbose){message(paste("    -> ", libName, " has insufficient reads. Removing", sep=""))}else{
				if(verbose){message(paste("    -> ",libName, " has high background (", (1-libraryQual)*100, " %). Removing", sep=""))}}
				lowQualList <- rbind(lowQualList, colQual)
			} else {
				qualList <- rbind(qualList, colQual)
			}
		}
		
		if(nrow(lowQualList) == ncol(strandTable))
		{
			warning("-> WARNING! ALL LIBRARIES ARE OF LOW QUALITY!! UNABLE TO REMOVE HIGH BACKGROUND LIBRARIES!")
		} else if(nrow(lowQualList) == 0 & verbose) { 
			message("-> All libraries of good quality" )
		}else{
			if(verbose){
				stCol <- ncol(strandTable)
				lqRow <- nrow(lowQualList)
				message(paste("-> Removed ", lqRow,
						  	" libraries from a total of ", stCol, ". ", 
							stCol-lqRow, " remaining (", 
							round((stCol-lqRow)/stCol*100, digits=1), "%)", sep="") )}
			strandTable <- strandTable[,!(colnames(strandTable) %in% lowQualList[,1])]
		}	
	}
	rawTable <- strandTable
	
	strandTable[strandTable >= strandTableThreshold] <- 1
	strandTable[strandTable <= -strandTableThreshold] <- 3
	strandTable[strandTable < strandTableThreshold & strandTable > -strandTableThreshold] <- 2

	preFilterData <- function(strandTable, filterThreshold, minLib, onlyWC=FALSE)
	{
		##### PRE-FILTER CONTIGS ##### 
		#Ignore contigs present in fewer than 10 libraries
		strandTable <- strandTable[which(apply(strandTable, 1, 
											   function(x){length(which(!is.na(x)))} >= minLib)),]

		#And ignore libraries that are entirely NA (indicating no cell present)	
		strandTable <- strandTable[,which(apply(strandTable, 2, 
												function(x){length(which(is.na(x)))} <= length(x)*filterThreshold ))]

		#Ignore contigs that are entirely WC (likely contain inversions) (NB all <10 contigs have already been excluded)
		WCvaluesCon <- apply(strandTable, 1, function(row) { sum(row == 2, na.rm=TRUE) / sum(is.element(row, c(1,2,3)), na.rm=TRUE) })		

		#And libraries that are entirely WC (whole genome libraries, not strandSeq)
		WCvaluesLib <- apply(strandTable, 2, function(col) { sum(col == 2, na.rm=TRUE) / sum(is.element(col, c(1,2,3)), na.rm=TRUE) })

		if(onlyWC)
		{
			strandTable <- strandTable[which(WCvaluesCon > filterThreshold),]
		}else{
			strandTable <-  strandTable[which(WCvaluesCon <= filterThreshold),]

			##### PRE-FILTER LIBRARIES #####
			#Ignore libraries that are mostly WC (indicating Strand-Seq failure)
			strandTable <- strandTable[,which(WCvaluesLib <= filterThreshold)]
		}
		return(strandTable)
	}

	#Create new data.frame of contigs that are entirely WC to investigate further
	strandTableAWC <- preFilterData(strandTable, 
									filterThreshold, 
									minLib,
									onlyWC=TRUE)

	strandTable <- preFilterData(strandTable, 
								 filterThreshold,
								 minLib)

	rawTable <- rawTable[rownames(strandTable),]
	rawTable <- rawTable[,colnames(strandTable)]

	#Order rows of strandTable by contig quality (best first)
	if(orderMethod=='libsAndConc')
	{
		if(verbose){message("-> Computing QA measures for contigs and sorting by best quality first")}
		contigNAs <- apply(strandTable, 1, 
						   function(x){length(which(!is.na(x)))}) / ncol(strandTable) #number of libraries contig is non-NA 
		#Compute the divergence between the call for a contig and the raw value used to make that call:
		computeOneAgreement <- function(contigName) {
			contigCall <- strandTable[contigName,]			
			contigRaw <- rawTable[contigName,]

			mean(ifelse(contigCall == 2, 
					(strandTableThreshold - abs(contigRaw)) / strandTableThreshold,
					(abs(contigRaw) - strandTableThreshold) / (1 - strandTableThreshold))
				[!is.na(contigCall)])
		}
		
		contigAgreement <- sapply(rownames(strandTable), computeOneAgreement)
		contigQA <- contigAgreement * contigNAs #compute a compound QA measure
		strandTable <- strandTable[names(sort(contigQA, decreasing=TRUE)),] # and sort
	}
		
	# Convert NaNs to NAs
	strandTable[is.na(strandTable)] <- NA
		
#	strandMatrix <- data.frame(lapply(strandTable, function(x){factor(x, levels=c(1,2,3))}))  
#	rownames(strandMatrix) <- rownames(strandTable)
	strandTable <- StrandStateMatrix(strandTable)

	return(list(strandMatrix=strandTable,
				qualList=qualList, 
				lowQualList=lowQualList, 
				AWCcontigs=row.names(strandTableAWC)))
}

# Copyright (c) 2015, Mark Hills & Kieran O'Neill
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

####################################################################################################
#' preprocessStrandTable -- remove low quality libraries and contigs before attempting to build
#' a genome
#' @param strandTable data.frame containing the strand table to use as input
#' @param strandTableThreshold threshold at which to call a contig WW or CC rather than WC
#' @param filterThreshold maximum number of libraries a contig can be NA or WC in
#' @param orderMethod the method to oder contigs. currently libsAndConc only option. Set to FALSE to not order contigs based on library quality
#' @param lowQualThreshold background threshold at which to toss an entire library. If NULL, function will not make an overall assessment of library quality.
#' Very chimeric assemblies can appear low quality across all libraries.
#' @param minLib minimum number of libraries a contig must be present in to be included in the output
#' @param verbose messages written to terminal
#' @aliases preprocessStrandTable preprocessStrandTable,StrandFreqMatrix,StrandFreqMatrix-method
#' 
#' @example inst/examples/preprocessStrandTable.R
#' 
#' @return A list of one matrix and three quality data.frames -- 1: a matrix of WW/WC/WW calls for all contigs; 3: the quality of libraries used (based on frequencies outside expected ranges); 4: A data.frame of libraries that are of low quality and therefore excluded from analysis; 5: contigs that are present as WC in more libraries than expected. These are excluded from the strandStateMatrix, but are potentially worth investigating for chimerism.
#' 
#' @export
#
####################################################################################################

setMethod('preprocessStrandTable',
		  signature = signature(strandTable='StrandFreqMatrix'),
		  definition = preprocessStrandTable.func
		  )
