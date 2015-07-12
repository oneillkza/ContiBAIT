####################################################################################################
#' preprocessStrandTable -- remove low quality libraries and contigs before attempting to build
#' a genome
#' @param strandTable data.frame containing the strand table to use as input
#' @param strandTableThreshold =0.8 threshold at which to call a contig WW or CC rather than WC
#' @param filterThreshold =0.8 maximum number of libraries a contig can be NA or WC in
#' @param orderContigs =TRUE whether to sort contigs by background quality or read as is
#' @param orderMethod ='libsAndConc' the method to oder contigs. currently libsAndConc only option
#' @param lowQualThreshold =0.9 background threshold at which to toss an entire library
#' @param verbose =TRUE messages written to terminal
#' 
#' @return A list of three matrices -- 1: WW/WC/WW; 2: WW/CC; 3: sex contigs
#' 
#' @export
#
####################################################################################################

preprocessStrandTable <- function(strandTable, strandTableThreshold=0.8, filterThreshold=0.8, orderContigs=TRUE, orderMethod='libsAndConc', lowQualThreshold=0.9, verbose=TRUE)
{
	strandTableLength <- nrow(strandTable)
	
	if(!is.data.frame(strandTable))
		strandTable <- data.frame(strandTable)
	
	# Filter low quality libraries.  Scan "WW" and "CC" background regions.
	lowQualList <- matrix(ncol=2, nrow=0)
	qualList <- matrix(ncol=2, nrow=0)
	if(verbose){message("-> Checking for high quality libraries")}

	for( col in seq(1:ncol(strandTable)) )
	{
		libraryQual <- round((abs(mean(strandTable[,col][which(strandTable[,col] < -0.6)], na.rm=TRUE) ) + abs(mean(strandTable[,col][which(strandTable[,col] > 0.6)], na.rm=TRUE))) / 2, digits=3)
		
		if(libraryQual < lowQualThreshold || libraryQual == "NaN")
		{
			if(verbose){message(paste("    -> ", colnames(strandTable[col], do.NULL=FALSE), " has high background (", (1-libraryQual)*100, " %). Removing", sep=""))}
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
		if(verbose){message(paste("-> Removed ", nrow(lowQualList), " libraries from a total of ", ncol(strandTable), sep="") )}
		strandTable <- strandTable[,!(names(strandTable) %in% lowQualList[,1])]
	}	

	rawTable <- strandTable
	
	strandTable[strandTable >= strandTableThreshold] <- 1
	strandTable[strandTable <= -strandTableThreshold] <- 3
	strandTable[strandTable < strandTableThreshold & strandTable > -strandTableThreshold] <- 2

	preFilterData <- function(strandTable, filterThreshold=0.8, onlyWC=FALSE)
	{
		##### PRE-FILTER CONTIGS ##### 
		#Ignore contigs present in fewer than 10 libraries
		strandTable <- strandTable[which(apply(strandTable, 1, function(x){length(which(!is.na(x)))} >= 10)),]

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
	strandTableAWC <- preFilterData(strandTable, filterThreshold=filterThreshold, onlyWC=TRUE)

	strandTable <- preFilterData(strandTable, filterThreshold=filterThreshold)
	rawTable <- rawTable[rownames(strandTable),]
	rawTable <- rawTable[,colnames(strandTable)]

	#Order rows of strandTable by contig quality (best first)
	if(orderContigs)
	{
		if(verbose){message("-> Computing QA measures for contigs and sorting by best quality first")}
		contigNAs <- apply(strandTable, 1, function(x){length(which(!is.na(x)))}) / ncol(strandTable) #number of libraries contig is non-NA 
		if(orderMethod=='libsAndConc')
		{
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
		} else
		{
			contigQA <- contigNAs #compute a compound QA measure
		}
		
		
		strandTable <- strandTable[names(sort(contigQA, decreasing=T)),] # and sort
	}
		
	#Create new table where WC reads are converted to NA, so only WW and CC relationships are considered
	strandTable2 <- replace(strandTable, strandTable == 2, NA)
	#Filter this table to exclude contigs that are now only present in fewer than 10 libraries
	strandTable2 <- strandTable2[which(apply(strandTable2, 1, function(x){length(which(!is.na(x)))} >= 10)),]
	strandTable <- strandTable[rownames(strandTable2),]
#	strandTable <- strandTable[which(apply(strandTable2, 1, function(x){length(which(!is.na(x)))} >= 10)),]

	#Scan for sex chromosomes based on no WC inheritance
	#First, count non-WC's for each line...
	if(verbose){message("-> Searching for contigs on sex chromosomes (no WC in libraries)")}
	numSex <- apply(strandTable, c(1), function(x){length(which(x != 2))})
	
	#Then, count number of NAs in each line
	naNum <- apply(strandTable, c(1), function(x){length(which(x != "NA"))})
	
	includeSex <- which(numSex > naNum*strandTableThreshold )
	
	if( length(includeSex) > 1)
	{
		strandMatrixSex <- strandTable[includeSex,] 
		if(verbose){message(paste("    -> ", nrow(strandMatrixSex), " found!", sep="") )}
	} else {
		if(verbose){message("    -> None found")}
		strandMatrixSex <-matrix(nrow=1, ncol=1)
	}
	#Turn each column into factors (categorical -- we don't want the distance metric to interpret 1,2,3 as numerical)
  strandMatrix <- data.frame(lapply(strandTable, function(x){factor(x, levels=c(1,2,3))}))  
  rownames(strandMatrix) <- rownames(strandTable)
	strandMatrix2 <- data.frame(lapply(strandTable2, function(x){factor(x, levels=c(1,2,3))}))	
  rownames(strandMatrix2) <- rownames(strandTable2)

	if (nrow(strandMatrixSex) > 3)
	{
		#Ignore contigs present in fewer than 10 libraries
		strandMatrixSex <- strandMatrixSex[which(apply(strandMatrixSex, 1, function(x){length(which(!is.na(x)))} >= 10)),]
		#And ignore libraries that are entirely NA (indicating no cell present)	
		strandMatrixSexTemp <- strandMatrixSex[,which(apply(strandMatrixSex, 2, function(x){length(which(is.na(x)))} <= nrow(strandMatrixSex)*filterThreshold ))]
		strandMatrixSex <- data.frame(apply(strandMatrixSex, c(2), function(x){factor(x, levels=c(1,2,3))})) 
		strandMatrixSex <- data.frame(lapply(strandMatrixSex, function(x){factor(x, levels=c(1,2,3))}))  
		rownames(strandMatrixSex) <- rownames(strandMatrixSexTemp)
	}
	
	#Filter out contigs that look like allosomes:
	strandMatrix <- strandMatrix[setdiff(rownames(strandMatrix),rownames(strandMatrixSex)),]
	strandMatrix2 <- strandMatrix2[setdiff(rownames(strandMatrix2),rownames(strandMatrixSex)),]
	
	strandMatrix <- new('StrandStateMatrix', strandMatrix)
	strandMatrix2 <- new('StrandStateMatrix', strandMatrix2)
	strandMatrixSex <- new('StrandStateMatrix', strandMatrixSex)

	
	return(list(strandMatrix=strandMatrix, strandMatrixWWCC=strandMatrix2, strandMatrixSex=strandMatrixSex, qualList=qualList, lowQualList=lowQualList, AWCcontigs=row.names(strandTableAWC)))
}