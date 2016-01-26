orderAllLinkageGroups.func <- function(linkageGroupList, strandStateMatrix, strandFreqMatrix, strandReadCount, whichLG=NULL, saveOrdered=FALSE, orderCall='greedy', randomAttempts=75, verbose=TRUE)
{

  combineZeroDistContigs <- function(linkageStrands, rawStrandTable, lg)
  {
    #Filter out borderline calls:
    
    contigNames <- rownames(linkageStrands)
    
    linkageMat <- as.matrix(linkageStrands)
    linkageMat <- apply(linkageMat, 2, as.numeric)
    
    rawStrandTable <- rawStrandTable[rownames(linkageStrands), ]
    linkageMat[which(abs(rawStrandTable) > 0.2 & abs(rawStrandTable) < 0.1 ) ] <- NA
    
    linkageStrands <- data.frame(linkageMat)

    linkageStrands <- data.frame(lapply(linkageStrands, function(x) factor(x, levels=c(1,2,3))))  
    rownames(linkageStrands) <- contigNames
    linkageStrands <- linkageStrands[apply(linkageStrands, 1, function(x) length(which(is.na(x))) != ncol(linkageStrands)) ,apply(linkageStrands, 2, function(x) length(which(is.na(x))) != nrow(linkageStrands))]
    
    ##Combine zero dist contigs:
    strandDist <- daisy(linkageStrands)
    strandDist <- as.matrix(strandDist)
    
    mergedContigs <- list()
    beenMerged <- vector()
    mergedStrands <- matrix(nrow=0, ncol=ncol(linkageStrands))
    groupCount <- 1
    
    for(contig in rownames(linkageStrands))
    {
      #Only merge contigs not already pulled in by other merges:
      if(!contig %in% beenMerged)
      {
        toMerge <- which(strandDist[contig,] == 0)
        #And don't pull in contigs that are present in toMerge if they've already been asigned
        toMerge <- toMerge[!names(toMerge) %in% beenMerged]
        beenMerged <- append(beenMerged, names(toMerge))
        mergedContigs[[paste('LG', lg, '.', groupCount, sep='')]] <- names(toMerge)
        mergedStrands <- rbind(mergedStrands, linkageStrands[contig,])
        groupCount <- groupCount +1
      }
    }
    
    orderedContigMatrix <- data.frame(LG=unlist(lapply(1:length(mergedContigs), 
                                      function(x) rep(names(mergedContigs[x]), 
                                      length(mergedContigs[[x]]) ))), 
                                      contig=unlist(mergedContigs), 
                                      row.names=NULL, 
                                      stringsAsFactors=FALSE )

    contigsByLG <- sapply(orderedContigMatrix$LG, function(x){strsplit(x, '[.]')[[1]]})
    contigStarts <- sub('-.*', '', sub('.*:', '', orderedContigMatrix$contig))
    orderedContigMatrix <- orderedContigMatrix[order(contigsByLG[1,], as.numeric(contigsByLG[2,] ), as.numeric(contigStarts)),]

    orderedContigMatrix <- new("ContigOrdering", orderedContigMatrix)

    mergedStrands <- data.frame(lapply(mergedStrands, function(x){factor(x, levels=c(1,2,3))}))  
    rownames(mergedStrands) <- names(mergedContigs)
    mergedStrands <- new("StrandStateMatrix", mergedStrands)

    return(list(mergedStrands=mergedStrands, contigKey=orderedContigMatrix))
  }

  if(is.null(whichLG)){whichLG=c(1:length(linkageGroupList))}
 orderedGroups <- data.frame(LG=vector(), name=vector())
 
  for(lg in whichLG)
  {
    if(verbose){message(paste('-> Ordering fragments in LG', lg, sep=""))}
    if(length(linkageGroupList[[lg]]) > 1)
    {

      linkageGroup <- linkageGroupList[[lg]]
      linkageGroupReadTable <- strandStateMatrix[linkageGroup,]
      zeroGroups <- combineZeroDistContigs(linkageGroupReadTable, strandFreqMatrix, lg)
      #Make a contig Weight vector
      zeroGroups[[2]]$weights <- apply(strandReadCount[which(rownames(strandReadCount) %in% zeroGroups[[2]]$contig ),] , 1, median)
      #The make a LG weight by taking the sum of all contigs within that LG, and order the linkageGroupTable based on the deepest LG
      linkageGroupReadTable <- zeroGroups[[1]][names(sort(sapply(unique(zeroGroups[[2]][,1]), function(x) sum(zeroGroups[[2]]$weights[which(zeroGroups[[2]][,1] == x)])), decreasing=TRUE)),]
      
      if(orderCall == 'greedy')
      {
        outOfOrder <- orderContigsGreedy(linkageGroupReadTable, randomAttempts=randomAttempts)
  
      }else if (orderCall == 'TSP')
      {
        outOfOrder <- orderContigsTSP(linkageGroupReadTable)
      }else{
        warning('###### WARNING! orderCall parameter not recognized! No ordering Performed. ######')
        break
      }

      mergedGroups <- data.frame(LG=vector(), name=vector())
      for(gp in 1:length(outOfOrder[[1]])){
        mergedGroups <- rbind(mergedGroups, zeroGroups[[2]][which(zeroGroups[[2]] == outOfOrder[[1]][gp]),1:2] )
      }
      orderedGroups <- rbind(orderedGroups, mergedGroups)

      chromosome <- strsplit(linkageGroupList[[lg]][1],':')[[1]][1]
      if(saveOrdered != FALSE)
      {
        similarLinkageStrands <- as.matrix(1-daisy(outOfOrder[[2]]))
        diag(similarLinkageStrands) <- 1
        if(nrow(similarLinkageStrands) > 1)
        {
          breaks <- seq(0, 100, length.out=101)/100 
          cols <- colorRampPalette(c("cyan", "blue", "grey30", "black", "grey30", "red", "orange"))
          suppressWarnings(heatmap.2(similarLinkageStrands, Rowv=NA, Colv=NA, dendrogram="none", revC=TRUE, col=cols(100), breaks=breaks, trace='none', main=paste('greedy-ordered ', chromosome, sep="")))
        }
      }
    }
  }  
  orderedGroups <- new("ContigOrdering", orderedGroups)
  return(orderedGroups)
}

####################################################################################################
#' Function to call contig ordering algorithms iteratively across each linkage group element
#' @useDynLib ContiBAIT
#' @param linkageGroupList list of vectors, each specifying which contigs belong in which linkage group (product of clusterContigs)
#' @param strandStateMatrix table of strand calls for all contigs (product of preprocessStrandTable)
#' @param strandFreqMatrix table of W:C read proportions (used for QC) (product of strandSeqFreqTable[[1]])
#' @param strandReadCount table of read counts (product of strandSeqFreqTable[[2]])  
#' @param whichLG vector of integers specifying the element(s) of linkageGroupList to be ordered (i.e. which specific linkage groups to try to order). Default is all LGs.
#' @param saveOrdered Will return a pdf of heatmaps for each linkage group; String entered becomes the fileName (default is saveOrderedPDF=FALSE)
#' @param orderCall currently either 'greedy' for greedy algorithm or 'TSP' for travelling salesperson alogrithm (default is 'greedy')
#' @param verbose Pringts messages to the terminal. Default is TRUE
#' @param randomAttempts iterger specifying number of randomized clusterings to identify the best ordering. Default is 75
#' @aliases orderAllLinkageGroups orderAllLinkageGroups,orderAllLinkageGroups-LinkageGroupList-StrandStateMatrix-StrandFreqMatrix-StrandReadMatrix-method
#' @rdname orderAllLinkageGroups
#' 
#' @return a data.frame of ordered contigs with linkage group names
#' 
#' @example inst/examples/orderAllLinkageGroups.R
#' @export
#' @include AllClasses.R
#' @import Rcpp
#' @import BH
#' @importFrom cluster daisy
#' @importFrom gplots heatmap.2
####################################################################################################

setMethod('orderAllLinkageGroups',
      signature = signature(linkageGroupList='LinkageGroupList', strandStateMatrix= 'StrandStateMatrix', strandFreqMatrix='StrandFreqMatrix', strandReadCount='StrandReadMatrix'),
      definition = orderAllLinkageGroups.func
      )
