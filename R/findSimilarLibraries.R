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
        if(!(is.null(hom.libraries)))
        {
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

          if(length(linkage.libraries) > 1)
          {
            groupCon <- computeConsensus(linkage.libraries[[1]], strandStateLibrary)
            maxOne <- sapply(c(1,3), function(x) length(which(groupCon == x)))
            maxOne <- which(maxOne == max(maxOne))[1]
            maxOne <- if(maxOne == 1){c("Watson", "Crick")}else{c("Crick", "Watson")}
            newNames <- c(paste(chrToUse[chrNum], "_mostly_", maxOne[1], "_(", length(linkage.libraries[[1]]), ")", sep="" ),
                           paste(chrToUse[chrNum], "_mostly_", maxOne[2], "_(", length(linkage.libraries[[2]]), ")", sep=""))

            if(maxOne[1] == "Crick")
            {
              topTwo <- LinkageGroupList(linkage.libraries[1:2], names=names(linkage.libraries[1:2]))
            }else{
              topTwo <- LinkageGroupList(list(linkage.libraries[[2]], linkage.libraries[[1]]), names=c(names(linkage.libraries[2]), names(linkage.libraries[1])))
            }
            names(topTwo) <- newNames[order(newNames)]
            return(topTwo)
          }
        }
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

