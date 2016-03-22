
####################################################################################################
#' thoroughBed -- function to merge chromosomes from libraries that have the same strand states 
#' 
#' @param bamFileList vector containing the location of the bams file to be read
#' @param relatedLibList list where each element contains all library names that show similar strand pattern. The product of findSimilarLibraries 
#' @param qual Mapping quality threshold. Default is 10 
#' @param pairedEnd Whether the bam files being read are in paired end format. Default is TRUE. Note,
#' since paired reads will be the same direction, only first mate read of pair is used in output to reduce file size
#' @param verbose prints messages to the terminal (default is TRUE)
#' 
#' @return a GRanges object comprising merged directional reads from all libraries in relatedLibList. 
#' @import Rsamtools
#' @import IRanges
#' @import GenomicRanges
#' @import GenomicFiles
#' @import GenomicAlignments
#' @importFrom S4Vectors DataFrame
#' @example inst/examples/thoroughBed.R
#' @export
#' @include AllClasses.R
####################################################################################################


thoroughBed <- function(bamFileList, relatedLibList, qual=10, pairedEnd=TRUE, verbose=TRUE)
{   
  directionalRange <- GRanges()
  #set parameters. Parameters are different if paired end data is used.  
  if(pairedEnd){
    # Handle Pair End by only taking the first mate read (the second mate read should be opposite strand in every case, so ignore)
    param <- ScanBamParam(flag=scanBamFlag(isFirstMateRead=TRUE, isUnmappedQuery=FALSE), mapqFilter=qual)
  }else{
    # Handle single end reads by just looking at strand direction 
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE), mapqFilter=qual)
  }

  #Open each bam and allocate chromosomes: ie. chrck library name against findSimilarLibrary Output.
  for(fileName in bamFileList)
  {
    includedInList <- names(which(unlist(relatedLibList) == fileName))

    if(length(includedInList) > 0)
    {
      if(verbose){message(' -> Processing ', fileName, ' [', which(bamFileList == fileName), '/', length(bamFileList), ']')}
      readLocations <- sapply(includedInList, function(x) strsplit(x, '_mostly')[[1]][1])
      readLocationsW  <- sapply(includedInList, function(x) strsplit(x, '_mostly_Watson')[[1]][1])

      readTable <- granges(readGAlignments(fileName, param=param))
      readTable <- readTable[which(seqnames(readTable) %in% readLocations)]
      #get strand data from mostly Watson libraries
      getWatsonStrand <- as.vector(strand(readTable[which(seqnames(readTable) %in% readLocationsW)] ))    
      #flip the read direction
      getWatsonStrand <- as.character(sapply(getWatsonStrand, function(x) if(x == '+') {return('-')}else if(x == '-'){return('+')} ))
      #and replace
      strand(readTable[which(seqnames(readTable) %in% readLocationsW)]) <- getWatsonStrand

      directionalRange <- append(directionalRange, readTable)
      #RETURN A GRANGES TO GO INTO REORIENTR!!!
    }else{ 
      if(verbose){message(' -> Processing ', fileName, ' [', which(bamFileList == fileName), '/', length(bamFileList), ']')}
    }
  }
  return(directionalRange)
}
