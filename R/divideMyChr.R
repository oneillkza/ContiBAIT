####################################################################################################
#' Creates a data.frame of chromosomes split into user-derived bin sizes
#' @param chrTable data.frame consisting of either chromosome name and length, or chromosome name, start and end. Can
#' generate from makeChrTable
#' @param splitBy numeric value dictating the size of each bin in nucleotides. Default is 200000. 
#' @example inst/examples/divideMyChr.R
#' @return a ChrTable object, containing the newly-divided chromosomes 
#' @export
####################################################################################################

divideMyChr <- function(chrTable, splitBy=200000)
{
	options(scipen=20)
	if(ncol(chrTable) == 3 )
	{
		chrTable$length <- chrTable[,3] - chrTable[,2]
	}else{
		chrTable$end <- chrTable[,2]
		chrTable$start <- 0
		chrTable <- subset(chrTable, select=c('chr', 'start', 'end', 'length'))
	}

	splitChromosomes <- data.frame(chr=vector(), start=vector(), end=vector())

	largeFragTab <- chrTable[which(chrTable$length >= splitBy),]
	for(chr in 1:nrow(largeFragTab))
	{
		myChr <- data.frame(chr=largeFragTab[chr,1], start=seq(0, largeFragTab$length[chr], by=splitBy), end=c(seq(splitBy, largeFragTab$length[chr], by=splitBy), largeFragTab$length[chr]) ) 
		splitChromosomes <- rbind(splitChromosomes, myChr)
	}
	smallFragTab <- data.frame(chr=chrTable[which(!(chrTable$length >= splitBy)),1], start= chrTable[which(!(chrTable$length >= splitBy)),2], end=chrTable[which(!(chrTable$length >= splitBy)),3] )
	splitChromosomes <- rbind(splitChromosomes, smallFragTab)
	rownames(splitChromosomes) <- paste(splitChromosomes$chr, ":", as.character(splitChromosomes$start), "-", as.character(splitChromosomes$end), sep="")
	splitChromosomes <- new("ChrTable", splitChromosomes)
	return(splitChromosomes)
}

