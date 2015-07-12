
####################################################################################################
# countDirectionReads -- function to process bam files for contiBAIT 
# 
# @param fileName  The name of the bam file to be read
# @param index  The index name for the file. Default is to use fileName
# @param qual  Mapping quality threshold. Default is 0
# @param rmdup  remove duplicates in output file. Default is FALSE 
# @param frequencies  create a table of frequencies for strands. Use TRUE when clustering contigs, and FALSE to
# return bed-formatted data.frame of reads
# @param indexCounter  counter to use when looping mutiple bame files (default is unknown)
# @param bamFileLength  complete number of bamFiles being analyzed (default is unknown)
# @param verbose prints messages to the terminal (default is TRUE)
# 
# @return all reads in bed format if frequencies=FALSE, or list of frequency of strand calls for each contig and number of reads if frequencies=TRUE
# @import Rsamtools
# @importFrom S4Vectors DataFrame
# @include AllClasses.R
####################################################################################################

countDirectionReads <- function(fileName, index=1, path='.', qual=0, rmdup=FALSE, frequencies=TRUE, verbose=TRUE, indexCounter='?', bamFileLength='?', filter=FALSE)
{
  if(index == 1){index <- fileName}
  
  if(verbose){message(paste('-> Creating contig table for index ', index, " [", indexCounter, "/", bamFileLength, "]", sep=""))}
  
  #read in first reads & create dataframe, then convert to table
  #
  if(length(filter) == 3 & length(list.files(path=path, pattern=paste(sub("./", "", fileName), '.bai', sep=''))) == 1) 
  { 
    # if we want to filter data we can do it quickly with an index
    temp.1.bam <- scanBam(fileName, param=ScanBamParam(flag=scanBamFlag(isFirstMateRead=TRUE), which=GRanges(seqnames = c(as.character(filter[,1])), ranges = IRanges(c(filter[,2]), c(filter[,3]) )), what=c("rname","pos","strand","mapq")))

    #if filter is a data.frame, merge the lists together before proceeding
    if(nrow(filter) > 1)
    {
      temp.1.dataframe <- data.frame(rname=vector(), pos=vector(), strand=vector(), mapq=vector(), stringsAsFactors=FALSE)
      fragNames <- c(1:nrow(filter))
      for(n in seq(1, nrow(filter)))
      {
        partialDataFrame <- as.data.frame(do.call("DataFrame", temp.1.bam[[n]]))
        if(nrow(partialDataFrame) != 0)
        {
          newContigName <- paste(as.character(filter[n,1]), fragNames[n], sep='~')
          levels(partialDataFrame[,1]) <- c(levels(partialDataFrame[,1]), newContigName)
          partialDataFrame$rname <- newContigName
          temp.1.dataframe <- rbind(temp.1.dataframe, partialDataFrame)
        }
        temp.1.dataframe <- temp.1.dataframe[which(temp.1.dataframe$pos != "NA"), ]
      }
    }else{
      temp.1.dataframe <- do.call("DataFrame", temp.1.bam)
      names(temp.1.dataframe) <- names(temp.1.bam[[1]])
    }
  
  }else{
    temp.1.bam <- scanBam(fileName, param=ScanBamParam(flag=scanBamFlag(isFirstMateRead=TRUE), what=c("rname","pos","strand","mapq")))
    temp.1.dataframe <- do.call("DataFrame", temp.1.bam)
  }
  
  if(qual != 0)
  {
    temp.1.dataframe <- temp.1.dataframe[which(temp.1.dataframe$mapq >= qual),]
  }
  if(rmdup)
  {
    temp.1.dataframe <- temp.1.dataframe[!duplicated(paste(temp.1.dataframe$rname, temp.1.dataframe$pos)),]
  }
  
  #read in second pair reads, create dataframe then convert to table
  if(length(filter) == 3 & length(list.files(path=path, pattern=paste(sub("./", "", fileName), '.bai', sep=''))) == 1) 
  {
    # if we want to filter data we can do it quickly with an index
    temp.2.bam <- scanBam(fileName, param=ScanBamParam(flag=scanBamFlag(isSecondMateRead=TRUE), which=GRanges(seqnames = c(as.character(filter[,1])), ranges = IRanges(c(filter[,2]), c(filter[,3]) )), what=c("rname","pos","strand","mapq")))

    #if filter is a data.frame, merge the lists together before proceeding
    if(nrow(filter) > 1)
    {
      temp.2.dataframe <- data.frame(rname=vector(), pos=vector(), strand=vector(), mapq=vector(), stringsAsFactors=FALSE)
      fragNames <- c(1:nrow(filter))
      for(n in seq(1, nrow(filter)))
      {
        partialDataFrame <- as.data.frame(do.call("DataFrame", temp.2.bam[[n]]))
        if(nrow(partialDataFrame) != 0)
        {         
          newContigName <- paste(as.character(filter[n,1]), fragNames[n], sep='~')
          levels(partialDataFrame[,1]) <- c(levels(partialDataFrame[,1]), newContigName)
          partialDataFrame$rname <- newContigName
          temp.2.dataframe <- rbind(temp.2.dataframe, partialDataFrame)
        }
        temp.2.dataframe <- temp.2.dataframe[which(temp.2.dataframe$pos != "NA"), ]
      }
    }else{
      temp.2.dataframe <- do.call("DataFrame", temp.2.bam)
      names(temp.2.dataframe) <- names(temp.2.bam[[1]])
    }
  }else{
    temp.2.bam <- scanBam(fileName, param=ScanBamParam(flag=scanBamFlag(isSecondMateRead=TRUE), what=c("rname","pos","strand","mapq")))
    temp.2.dataframe <- do.call("DataFrame", temp.2.bam)
  }

  if(qual != 0)
  {
    temp.2.dataframe <- temp.2.dataframe[which(temp.2.dataframe$mapq >= qual),]
  }
  if(rmdup)
  {
    temp.2.dataframe <- temp.2.dataframe[!duplicated(paste(temp.2.dataframe$rname, temp.2.dataframe$pos)),]
  }
  
  #Now, flip all reverse complement reads to the different strand
  #add a new level to the factor
  levels(temp.2.dataframe[,2]) <- c(levels(temp.2.dataframe[,2]), "1")
  temp.2.dataframe$strand <- replace(temp.2.dataframe$strand, temp.2.dataframe$strand == '+', '1')
  temp.2.dataframe$strand <- replace(temp.2.dataframe$strand, temp.2.dataframe$strand == '-', '+')
  temp.2.dataframe$strand <- replace(temp.2.dataframe$strand, temp.2.dataframe$strand == '1', '-')		
  bam.dataframe <- rbind(as.data.frame(temp.1.dataframe), as.data.frame(temp.2.dataframe))
   
  if(frequencies == TRUE)
  {
    if(length(filter) == 3)
    {
      fragments <- paste(as.character(filter[,1]), fragNames, sep='~')    
      missingValues <- fragments[!(fragments %in% bam.dataframe$rname)]
      #if any values not represented, append dummy line to the end of the dataframe
      if(length(missingValues) >= 1 & nrow(bam.dataframe) != 0)
      {
        for(dummy in seq(1,length(missingValues)))
        {
          bam.dataframe <- rbind(bam.dataframe, data.frame(rname=missingValues[dummy], pos=0, strand="*", mapq=0, stringsAsFactors=FALSE))
        }
      }
    }
    table.freq <- table(bam.dataframe$rname, bam.dataframe$strand)

 
#    table.2 <- table(temp.2.dataframe$rname, temp.2.dataframe$strand)
    #Calculate strand frequencies by adding the Watson reads together; plus strand read 1 (table.1, column1) and the minus strand read 2 (table.2, column2)
    #then subtracting the Crick reads; minus strand read 1 (table.1, column1) and the plus strand read 2 (table.2, column1)
    #then divide by total number of reads to get frequency between -1 and +1
    if(length(table.freq) == 0)
    {
      #duplicateLines <- function(x, y, rows=1){x[rep(rows,y),]}
      #strand.frequencies <- data.frame(index, rep(0,2))
      #strand.frequencies <- duplicateLines(strand.frequencies, length(fragments))
      #read.counts <- data.frame(index, 0)
      #read.counts <- duplicateLines(read.counts, length(fragments))
      strand.frequencies <- c(index, rep(0, length(fragments)))
      read.counts <- c(index, rep(0, length(fragments)))

    }else{
      strand.frequencies <- c(index, ((table.freq[,2])-(table.freq[,1]))/(table.freq[,1]+table.freq[,2]))
      strand.frequencies[strand.frequencies=="NaN"] <- 0
      read.counts <- c(index, table.freq[,2]+table.freq[,1])

      #Create quality metric read.counts by adding all + and - reads from 1st and 2nd pairs to get total
    }
  
    #  if(length(strand.frequencies) != nrow(filter)+1)
      return(list(strandFreq=strand.frequencies, readCounts=read.counts))
       
  }else{
    bam.dataframe <- rbind(as.data.frame(temp.1.dataframe), as.data.frame(temp.2.dataframe))
    bam.dataframe <- bam.dataframe[order(bam.dataframe$rname, bam.dataframe$pos),]
    rownames(bam.dataframe) <- NULL
    #bam.dataframe <- new('RawReadStrands', bam.dataframe)

    if(length(filter) == 3 & length(list.files(path=path, pattern=paste(sub("./", "", fileName), '.bai', sep=''))) != 1) 
    {
      if(verbose){message("BAM FILES ARE NOT INDEXED! FILTERING UNINDEXED BAM FILES WILL BE SLOW")}
 
      filteredBam <- data.frame(rname=vector(), pos=vector(), strand=vector(), mapq=vector())
      for(line in seq(1, nrow(filter)))
      {
        bam.bit <- bam.dataframe[which(bam.dataframe$rname == as.character(filter[1, line]) & bam.dataframe$pos >= filter[2, line] & bam.dataframe$pos <= filter[3, line]),]
        filteredBam <- cbind(filteredBam, bam.bit)
      }
      bam.dataframe <- filteredBam
    }

    return(bam.dataframe)
  }
}
