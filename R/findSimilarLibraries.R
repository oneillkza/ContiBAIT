findSimilarLibraries.func <- function(strandStateMatrix, 
                                      strandReadMatrix, 
                                      chrGrange, 
                                      chrNum, 
                                      cluster=1, 
                                      clusterParam=NULL, 
                                      verbose=TRUE)
{

  getLibWeight <- function(allStrands, strandNum)
  {
    filtWeight <- strandNum[which(rownames(strandNum) %in% rownames(allStrands)  ),]
    libWeight <- apply(filtWeight, 1, function(x) median(x, na.rm=TRUE))
  }

  proportionHet <- function(groupMembers, allStrands)
  {
    if(length(groupMembers) > 1)
    {
      groupStrands <- allStrands[groupMembers,]
      mean(sapply(1:nrow(groupStrands), function(y) length(grep(2, groupStrands[y,]))/length(grep("1|2|3", groupStrands[y,]))))
    }
  }

  #Only take the top two groups: the 'mostly Watson' and 'mostly Crick'
  findBestMatch <- function(libNum)
  { 
    groupCon <- computeConsensus(linkage.libraries[[libNum]], strandStateLibrary)
    maxOne <- sapply(c(1,3), function(x) length(which(groupCon == x)))
    maxOne <- which(maxOne == max(maxOne))
    maxOne <- if(maxOne == 1){paste("Watson")}else{paste("Crick")}
    return(paste(chrToUse[chrNum], "_mostly_", maxOne, "_(", length(linkage.libraries[[libNum]]), ")", sep="" ) )
  }

  chrToUse <- unique(seqnames(chrGrange)[which(duplicated(as.character(seqnames(chrGrange))))])

  if(verbose){message("  -> Running analysis on ", chrToUse[chrNum], " [", chrNum, "/", length(chrToUse), "]")}

  strandStateMatrix <- strandStateMatrix[which(rownames(strandStateMatrix) %in% chrGrange$name[which(seqnames(chrGrange) == chrToUse[chrNum])]),]
  strandReadMatrix <- strandReadMatrix[which(rownames(strandReadMatrix) %in% chrGrange$name[which(seqnames(chrGrange) == chrToUse[chrNum])]),]
  if(!(is.null(nrow(strandStateMatrix))))
  { 
    if(class(strandStateMatrix) == 'matrix')
    {
      strandStateMatrix <- StrandStateMatrix(strandStateMatrix)

      strandStateLibrary <- t(strandStateMatrix)
      strandReadLibrary <- t(strandReadMatrix)

      libWeight <- getLibWeight(strandStateLibrary, strandReadLibrary)

      linkage.libraries <- clusterContigs(strandStateLibrary, 
                                         randomWeight=libWeight, 
                                         clusterParam=clusterParam, 
                                         similarityCutoff=0.7,
                                         recluster=cluster,
                                         clusterBy='hetero',
                                         randomise=TRUE,
                                         verbose=FALSE)

      tables <- do.call(cbind, lapply(linkage.libraries, proportionHet, strandStateLibrary))

      if(!(is.null(tables)))
      {
        hom.libraries <- linkage.libraries[[which(tables < 0.9)[1]]]
        strandStateHom <- strandStateLibrary[hom.libraries,]
        strandReadHom <- strandReadLibrary[hom.libraries,]
        libWeight <- getLibWeight(strandStateHom, strandReadHom)

        strandStateHom <- StrandStateMatrix(strandStateHom)
        #at least 20% of contigs must overlap for a call to be made
        atLeast20percent <- floor(length(which(seqnames(chrGrange) == chrToUse[chrNum]))/5)

        linkage.libraries <- clusterContigs(strandStateHom, 
                                           randomWeight=libWeight, 
                                           clusterParam=clusterParam, 
                                           similarityCutoff=0.9,
                                           minimumLibraryOverlap=atLeast20percent,
                                           recluster=cluster,
                                           clusterBy='homo',
                                           randomise=TRUE,
                                           verbose=FALSE)
        newNames <- sapply(1:2, function(x) findBestMatch(x))
        topTwo <- LinkageGroupList(linkage.libraries[1:2])
        names(topTwo) <- newNames

        return(topTwo)
      }
    }
  }
}


####################################################################################################
#' findSimilarLibraries -- function to identify libraries that hare similar WC patterns on chromosomes 
#' 
#' @param strandStateMatrix  A strandStateMatrix object for all libraries across split fragments, derived from preprocessStrandTable
#' @param strandReadMatrix  The number of reads present for each strandStateMatrix element. An object of type strandReadMatrix from strandSeqFreqTable
#' @param chrGrange File of type ChrTable (a GRanges object with a meta column titled 'name' determining contig name) to split chromosomes based on locations.
#' name meta should match the rownames of strandStateMatrix and strandReadMatrix
#' @param chrNum The chromosome number to analyse
#' @param cluster Number of times to recluster and take the consensus of. If NULL, clustering is 
#' run only once
#' @param clusterParam optional \code{BiocParallelParam} specifying cluster to use for parallel execution.
#' When \code{NULL}, execution will be serial. 
#' @param verbose prints messages to the terminal (default is TRUE)
#' 
#' @return a list of type LinkageGroupList with two elements; libraries that are mostly Watson, and those that are mostly Crick
#' @import Rsamtools
#' @import IRanges
#' @import GenomicRanges
#' @importFrom S4Vectors DataFrame
#' @example inst/examples/findSimilarLibraries.R
#' @aliases findSimilarLibraries findSimilarLibraries,findSimilarLibraries-StrandStateMatrix-StrandReadMatrix-ChrTable-method
#' @rdname findSimilarLibraries
#' @export
#' @include AllClasses.R
####################################################################################################

setMethod('findSimilarLibraries',
      signature = signature(strandStateMatrix='StrandStateMatrix', strandReadMatrix='StrandReadMatrix', chrGrange='ChrTable'),
      definition = findSimilarLibraries.func
      )

