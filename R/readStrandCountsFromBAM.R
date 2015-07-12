####################################################################################################
#' countDirectionReads -- function to process bam files for contiBAIT 
#' 
#' @param bamDir directory to scan for bam files
#' @param qualLimit minimum quality for an individual read to be accepted. default is 10
#' @param readLimit minimum number of reads to make a state call. contigs with less reads than this limit return NA. Default is 10
#' @param dups  remove duplicates in output file. Default is FALSE 
#' @param field  The region to use as a file identifier eg. 401.sort.bam has 3 fields seperated by periods. default is 1.
#' @param fieldSep  The seperator to denote different fields in a fileName. Default is '.'
#' @param verbose prints messages to the terminal (default is TRUE)
#' @return list of two data.frames. The first displays state calls for every contig (rows) in every library (cols). The second is a matching data.frame with the number of reads used to make the call
#' @import Rsamtools
#' @importFrom S4Vectors DataFrame
#' @export
####################################################################################################



readStrandCountsFromBAM <- function(bamDir, 
                                    qualLimit=10, 
                                    readLimit=10, 
                                    dups=FALSE, 
                                    verbose=TRUE, 
                                    fieldSep=".", 
                                    field=1 
                                    #filter=FALSE, 
                                    #freq=FALSE
                                    )
{
  filter=FALSE 
  freq=FALSE
	findIndex <- function(fileName, fieldSep=".", field=1)
	{
		#removes the ./ as the start of fileName if present
		fileName <- sub("./","", fileName)
		indexer <- as.character(paste("echo ", fileName, "| awk -F", fieldSep, " '{print $", field, "}' ", sep=""))
		return(index=system(indexer, intern=TRUE))
	}
	
	#Function to tidy up the matrices, convert first row (index) to the column name
	#And remove convert value from strandTable to NA if number of reads is less than countLimit
	limitStrandCounts <- function(strandTable)
	{
		colnames(strandTable) <- strandTable[1,]
		strandTable <- strandTable[-1,]
		mode(strandTable) <- "numeric"
   		return(strandTable)
	}
	
	# Create a list of all bam files in the working directory
	bamFileList <- list.files(path=bamDir, pattern=".bam$")
	bamFileLength <- length(bamFileList)
	indexCounter <- 1
	
	#read in first bam file and count the number of fragments in the header
if(length(filter) != 3)
{
	if(verbose){message(paste("-> Counting number of fragments", sep=""))}
	lengthOfContigs <- scanBamHeader(paste(bamDir,bamFileList[1],sep='/'))[[1]][["text"]]
	lengthOfContigs <- length(lengthOfContigs[which(grepl("\\<SN:", lengthOfContigs))])
	if(verbose){message(paste("-> ", lengthOfContigs," fragments found", sep=""))}
}else{
	lengthOfContigs <- nrow(filter)
}
	
	#create empty matrices with the same number of columns as contigs
	if(freq == TRUE){strandTable <- matrix(nrow=lengthOfContigs+1, ncol=0)}else{strandTable <- matrix(nrow=0, ncol=5)}
	countTable <- matrix(nrow=lengthOfContigs+1, ncol=0)
	
	for(fileName in bamFileList)
	{
		
		index <- findIndex(fileName, fieldSep=fieldSep, field=field)
		#processBam <- countDirectionReads(paste(bamDir,fileName,sep='/'), index=index, path=bamDir, qual=qualLimit, rmdup=dups, indexCounter=indexCounter, bamFileLength=bamFileLength, verbose=verbose, filter=filter, frequencies=freq)
		processBam <- countDirectionReads(fileName, index=index, path=bamDir, qual=qualLimit, rmdup=dups, indexCounter=indexCounter, bamFileLength=bamFileLength, verbose=verbose, filter=filter, frequencies=freq)

		if(freq == FALSE)
		{
			if(nrow(processBam) > 0)
			{
				processBam <- cbind(processBam, index)
				strandTable <- rbind(strandTable, processBam)
			}
		}else{
			strandTable <- cbind(strandTable, processBam[[1]])
			countTable <- cbind(countTable, processBam[[2]])
		}		

		indexCounter <- indexCounter + 1
	}	

	if(freq == FALSE)
	{
		return(strandTable)
	}else{

		strandTable <- limitStrandCounts(strandTable)
		countTable <- limitStrandCounts(countTable)
		#For fragments with less than countLimit reads, change strandTable value to NA
		strandTable <- replace(strandTable, countTable <= readLimit, NA)
		#Remove any rows that are entirely NA

	  return(list(strandTable, countTable))
	}
}
