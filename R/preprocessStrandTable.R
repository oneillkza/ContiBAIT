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
#' @param strandTableThreshold=0.8 threshold at which to call a contig WW or CC rather than WC
#' @param filterThreshold=0.8 maximum number of libraries a contig can be NA or WC in
#' @param orderMethod='libsAndConc' the method to oder contigs. currently libsAndConc only option. Set to FALSE to not order contigs based on library quality
#' @param lowQualThreshold=0.9 background threshold at which to toss an entire library
#' @param ignoreInternalQual=FALSE logical that prevents function for making an overall assessment of library quality. Very chimeric assemblies can appear low quality across all libraries. 
#' @param verbose=TRUE messages written to terminal
#' 
#' @return A list of three matrices -- 1: WW/WC/WW; 2: WW/CC; 3: sex contigs
#' 
#' @export
#
####################################################################################################

preprocessStrandTable <- function(strandTable, strandTableThreshold=0.8, filterThreshold=0.8, orderMethod='libsAndConc', lowQualThreshold=0.9, verbose=TRUE, minLib=10, ignoreInternalQual=FALSE)
{
	strandTableLength <- nrow(strandTable)
	
	if(!is.data.frame(strandTable))
			strandTable <- data.frame(strandTable)
	
	# Filter low quality libraries.  Scan "WW" and "CC" background regions.
	lowQualList <- matrix(ncol=2, nrow=0)
	qualList <- matrix(ncol=2, nrow=0)

if(ignoreInternalQual == FALSE)
{
	if(verbose){message("-> Checking for high quality libraries")}

	for( col in seq(1:ncol(strandTable)) )
	{
		libraryQual <- round((abs(mean(strandTable[,col][which(strandTable[,col] < -0.6)], na.rm=TRUE) ) + abs(mean(strandTable[,col][which(strandTable[,col] > 0.6)], na.rm=TRUE))) / 2, digits=3)
		
		if(libraryQual < lowQualThreshold || libraryQual == "NaN")
		{
			if(libraryQual == "NaN" & verbose){message(paste("    -> ", colnames(strandTable[col], do.NULL=FALSE), " has insufficient reads. Removing", sep=""))}else{
			if(verbose){message(paste("    -> ", colnames(strandTable[col], do.NULL=FALSE), " has high background (", (1-libraryQual)*100, " %). Removing", sep=""))}}
			lowQualList <- rbind(lowQualList, cbind(colnames(strandTable[col], do.NULL=FALSE), libraryQual))
		} else {
			qualList <- rbind(qualList, cbind(colnames(strandTable[col], do.NULL=FALSE), libraryQual))
		}
	}

	
	if(nrow(lowQualList) == ncol(strandTable))
	{
		warning("-> WARNING! ALL LIBRARIES ARE OF LOW QUALITY!! UNABLE TO REMOVE HIGH BACKGROUND LIBRARIES!")
	} else if(nrow(lowQualList) == 0) { 
		if(verbose){message("-> All libraries of good quality" )}
	} else {
		if(verbose){message(paste("-> Removed ", nrow(lowQualList), " libraries from a total of ", ncol(strandTable), ". ", ncol(strandTable)-nrow(lowQualList), " remaining (", round((ncol(strandTable)-nrow(lowQualList))/ncol(strandTable), digits=1), "%)", sep="") )}
		strandTable <- strandTable[,!(names(strandTable) %in% lowQualList[,1])]
	}	
}
	rawTable <- strandTable
	
	strandTable[strandTable >= strandTableThreshold] <- 1
	strandTable[strandTable <= -strandTableThreshold] <- 3
	strandTable[strandTable < strandTableThreshold & strandTable > -strandTableThreshold] <- 2

	preFilterData <- function(strandTable, filterThreshold=0.8, onlyWC=FALSE, minLib=minLib)
	{
		##### PRE-FILTER CONTIGS ##### 
		#Ignore contigs present in fewer than 10 libraries
		strandTable <- strandTable[which(apply(strandTable, 1, function(x){length(which(!is.na(x)))} >= minLib)),]

		#Ignore contigs that are entirely WC (likely contain inversions) (NB all <10 contigs have already been excluded)
		if(onlyWC)
		{
			strandTable <- strandTable[which(apply(strandTable, 1, function(x){length(which(x == 2))} > ncol(strandTable)*filterThreshold )),]
		}else{
			strandTable <- strandTable[which(apply(strandTable, 1, function(x){length(which(x == 2))} <= ncol(strandTable)*filterThreshold )),]
		}
		##### PRE-FILTER LIBRARIES #####
		#Ignore libraries that are mostly WC (indicating Strand-Seq failure)
		strandTable <- strandTable[,which(apply(strandTable, 2, function(x){length(which(x == 2))} <= nrow(strandTable)*filterThreshold ))]
		#And ignore libraries that are entirely NA (indicating no cell present)	
		strandTable <- strandTable[,which(apply(strandTable, 2, function(x){length(which(is.na(x)))} <= nrow(strandTable)*filterThreshold ))]
		return(strandTable)
	}

	#Create new data.frame of contigs that are entirely WC to investigate further
	strandTableAWC <- preFilterData(strandTable, filterThreshold=filterThreshold, onlyWC=TRUE, minLib=minLib)

	strandTable <- preFilterData(strandTable, filterThreshold=filterThreshold, minLib=minLib)
	rawTable <- rawTable[rownames(strandTable),]
	rawTable <- rawTable[,colnames(strandTable)]

	#Order rows of strandTable by contig quality (best first)
	if(orderMethod=='libsAndConc')
	{
		if(verbose){message("-> Computing QA measures for contigs and sorting by best quality first")}
		contigNAs <- apply(strandTable, 1, function(x){length(which(!is.na(x)))}) / ncol(strandTable) #number of libraries contig is non-NA 
		#Compute the divergence between the call for a contig and the raw value used to make that call:
		computeOneAgreement <- function(contigName)
		{
			contigDists <- vector()
			for(lib in 1:ncol(strandTable))
			{
				contigCall <- strandTable[contigName, lib]
				if(!is.na(contigCall))
				{
					contigRaw <- rawTable[contigName, lib]
					if(contigCall == 2)
					{
						contigDist <- (strandTableThreshold - abs(contigRaw)) / strandTableThreshold
					}else
					{
						contigDist <- (abs(contigRaw) - strandTableThreshold) / (1-strandTableThreshold)
					}
					contigDists <- append(contigDists, contigDist)
				}
			}
			mean(contigDists)
		}
		
		contigAgreement <- sapply(rownames(strandTable), computeOneAgreement)
		contigQA <- contigAgreement * contigNAs #compute a compound QA measure
		strandTable <- strandTable[names(sort(contigQA, decreasing=T)),] # and sort
	}
		
	# Convert NaNs to NAs
	strandTable[is.na(strandTable)] <- NA
		
	#Create new table where WC reads are converted to NA, so only WW and CC relationships are considered
	strandTable2 <- replace(strandTable, strandTable == 2, NA)
	#Filter this table to exclude contigs that are now only present in fewer than 10 libraries
	strandTable2 <- strandTable2[which(apply(strandTable2, 1, function(x){length(which(!is.na(x)))} >= minLib)),]
	strandTable <- strandTable[rownames(strandTable2),]

	#Scan for sex chromosomes based on no WC inheritance
	#First, count non-WC's for each line...
	if(verbose){message("-> Searching for contigs on sex chromosomes (no WC in libraries)")}
	
	#First, count non-WC's for each line and divide by total number of libraries (which aren't NA).Only include if above strandTableThreshold
	includeSex <- which(apply(strandTable, 1, function(x){length(which(x != 2)) /length(which(x != "NA")) }) >= strandTableThreshold) 
	
	if( length(includeSex) > 1)
	{
		strandMatrixSex <- strandTable[includeSex,] 
		if(verbose){message(paste("    -> ", nrow(strandMatrixSex), " found!", sep="") )}
	} else {
		if(verbose){message("    -> None found")}
		strandMatrixSex <- matrix(nrow=2, ncol=ncol(strandTable))
		strandMatrixSex <- data.frame(apply(strandMatrixSex, 2, function(x){as.factor(x)}) )
 	}

	strandMatrix <- data.frame(lapply(strandTable, function(x){factor(x, levels=c(1,2,3))}))  
	rownames(strandMatrix) <- rownames(strandTable)
	strandMatrix2 <- data.frame(lapply(strandTable2, function(x){factor(x, levels=c(1,2,3))}))	
	rownames(strandMatrix2) <- rownames(strandTable2)


	if (nrow(strandMatrixSex) > 2)
	{
		#Ignore contigs present in fewer than 10 libraries
		strandMatrixSex <- strandMatrixSex[which(apply(strandMatrixSex, 1, function(x){length(which(!is.na(x)))} >= minLib)),]
		#And ignore libraries that are entirely NA (indicating no cell present)	
		strandMatrixSex <- strandMatrixSex[,which(apply(strandMatrixSex, 2, function(x){length(which(is.na(x)))} <= nrow(strandMatrixSex)*filterThreshold ))]
		
		strandMatrixSex <- data.frame(lapply(strandMatrixSex, function(x){factor(x, levels=c(1,2,3))}))  
		strandMatrixSex <- new('StrandStateMatrix', strandMatrixSex)
	}
	
	#Filter out contigs that look like allosomes:
	strandMatrix <- strandMatrix[setdiff(rownames(strandMatrix),rownames(strandMatrixSex)),]
	strandMatrix2 <- strandMatrix2[setdiff(rownames(strandMatrix2),rownames(strandMatrixSex)),]


	strandMatrix <- new('StrandStateMatrix', strandMatrix)
	strandMatrix2 <- new('StrandStateMatrix', strandMatrix2)

	
	return(list(strandMatrix=strandMatrix, strandMatrixWWCC=strandMatrix2, strandMatrixSex=strandMatrixSex, qualList=qualList, lowQualList=lowQualList, AWCcontigs=row.names(strandTableAWC)))
}