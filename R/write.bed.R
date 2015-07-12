####################################################################################################
#' write.bed -- function to merge and write like-data to a bed file
#' 
#' @param libList vector of library names to be analyzed; created using mergeLibraries 
#' @param chr  chromoosome name to be analyzed
#' @param libCall  Number reflecting the state of the bed files; 1=Crick, 3=Watson, 10=analyze allFiles. Default is 10 
#' @param splitBed  Logical value which when TRUE will place bed headers within file, generating split (unmerged) bedfiles. Default is FALSE
#' @param skip File containing known gap locations in bed format. Either gapFile or binData must be specified
#' @param verbose prints messages to the terminal (default is TRUE)
#' @param byChr  Creates bedFile consisting of reads from only the chromosome being analyzed. Default is TRUE
#' @param file outPut fileName. Default is 'bedfile'
#' @param path the location of the bed files in which to merge. Default is current directory.
#' 
#' @return bed file printed to working directory names by fileName, chromosme and state.
#' 
#' @export
####################################################################################################


write.bed <- function(libList, chr, libCall=10, splitBed=FALSE, skip=3, verbose=TRUE, byChr=TRUE, file='bedfile', path='.')
{
  #combine bed files from bedfilelist and save as crick bedfile
  #	makeBed <- data.frame(chr=vector(), start=vector(), end=vector(), name=vector(), qual=vector(), strand=vector())
  
  #	cutPlotNames <- names(cutPlot[which(cutPlot == cluster[,1])])
  #	BedList <- bedFileList[bedFileList %in% paste(cutPlotNames, ".bed", sep="")]
  
  if(libCall == 1){state <- "Crick"}else if (libCall == 3){state <- "Watson"} else if (libCall == 10){state <- "allFiles"} else { warning('Cannot call fragment WW or CC. Consider changing libNum, skipping')}
  if(verbose){message(paste("  -> Merging ", state, " clustered files [", length(libList), " beds]", sep=""))}
  if(skip !=0)
  {
    #If there's a first header print it for the bedgraph
    firstHead <- data.frame(readLines(paste(path, '/', libList[1], sep=""), n=skip), as.integer(c(rep("", skip))), as.integer(c(rep("", skip))), c(rep("", skip)), as.integer(c(rep("", skip))), c(rep("", skip))  )
  }else{
    #Or create a generic header
    firstHead <- data.frame(chr=vector(), start=vector(), end=vector(), name=vector() , qual=vector(), strand=vector()  )
  }
  
  #then save this onto the harddrive
  write.table(firstHead, file=paste(file, "_", chr, "_", state, "_", length(libList), "merged.txt", sep=""), row.names=FALSE, col.names=FALSE, na="", quote=FALSE, append=TRUE)

  bedCount <- 1
  
  for(bedLine in seq(1:length(libList)))
  {
    bed <- libList[bedLine]
    if(verbose){message(paste("    -> Writing ", bed, " for ", state, " ", chr, " [", bedCount, "/", length(libList), " bed files]", sep=""))}
    fullBedFile <- read.table(paste(path, '/', bed, sep=""), skip=skip)
    if(byChr){bedFile <- fullBedFile[which(fullBedFile[,1] == as.character(chr)),]}else{bedFile <- fullBedFile}
    
    
    if(splitBed & skip != 0)
    {
      header <- data.frame(readLines(paste(path, '/', libList[bedLine], sep=""), n=skip), as.integer(c(rep("", skip))), as.integer(c(rep("", skip))), c(rep("", skip)), as.integer(c(rep("", skip))), c(rep("", skip))  )
      colnames(header) <- colnames(bedFile)
      bedFile <- rbind(header, bedFile)
      write.table(bedFile, file=paste(file, "_", chr, "_", state, "_", length(libList), "merged.txt", sep=""), row.names=FALSE, col.names=FALSE, na="", quote=FALSE, append=TRUE)
    } else {
      #makeBed <- rbind(makeBed, bedFile)
      write.table(bedFile, file=paste(file, "_", chr, "_", state, "_", length(libList), "merged.txt", sep=""), row.names=FALSE, col.names=FALSE, na="", quote=FALSE, append=TRUE)			
    }

#    if(shrinkBed == TRUE)
 #   {
  #    if(state == "Crick")
   #   {
        #remove all crick reads for the current chromosome, then resaves bed file
 #       if(verbose){message()}
  #      truncatedFile <- fullBedFile[which(paste(fullBedFile[,1], fullBedFile[,6], sep="") != paste(as.character(chr), '+', sep="")),]
   #     write.table(truncatedFile, file=paste(path, '/', bed, sep=""), append=FALSE)
#      } 
 #     else if
  #    (state == "Watson")
   #   {
        #remove all watson reads for the current chromosome, then resaves bed file 
#        truncatedFile <- fullBedFile[which(paste(fullBedFile[,1], fullBedFile[,6], sep="") != paste(as.character(chr), '+', sep="")),]
 #       write.table(truncatedFile, file=paste(path, '/', bed, sep=""), append=FALSE)        
  #    }
    #}
    bedCount <- bedCount + 1
  }
}
