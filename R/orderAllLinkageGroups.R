orderAllLinkageGroups.func <- function(linkageGroupList, strandStateMatrix, strandFreqMatrix, strandReadCount, whichLG=NULL, saveOrdered=NULL, orderCall='greedy', randomAttempts=75, verbose=TRUE)
{

  combineZeroDistContigs <- function(linkageStrands, rawStrandTable, lg)
  {
    #Filter out borderline calls:
    
    contigNames <- rownames(linkageStrands)
    
    rawStrandTable <- rawStrandTable[contigNames, ]
    linkageStrands[which(abs(rawStrandTable) > 0.2 & abs(rawStrandTable) < 0.1 ) ] <- NA

    linkageMat <- linkageStrands
    linkageStrands <- data.frame(linkageStrands)
    linkageStrands <- data.frame(lapply(linkageStrands, function(x) factor(x, levels=c(1,2,3))))  
    rownames(linkageStrands) <-  contigNames
 
    linkRows <- apply(linkageStrands, 1, function(x) length(which(is.na(x))) != ncol(linkageStrands))
    linkCols <- apply(linkageStrands, 2, function(x) length(which(is.na(x))) != nrow(linkageStrands))
    linkageStrands <- linkageStrands[linkRows , linkCols]
    linkageMat <- linkageMat[linkRows, linkCols]   
    ##Combine zero dist contigs:
    strandDist <- as.matrix(daisy(linkageStrands))
    
    mergedContigs <- list()
    beenMerged <- vector()
    mergedStrands <- matrix(nrow=nrow(linkageStrands), ncol=ncol(linkageStrands))
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
        mergedStrands[groupCount,] <- linkageMat[contig,]
        groupCount <- groupCount +1
      }
    }
    colnames(mergedStrands) <- colnames(linkageMat)
    mergedStrands <- mergedStrands[rowSums(is.na(mergedStrands)) != ncol(mergedStrands),]

    orderedContigMatrix <- matrix(c(unlist(lapply(1:length(mergedContigs), 
                                      function(x) rep(names(mergedContigs[x]), 
                                      length(mergedContigs[[x]]) ))), 
                                      unlist(mergedContigs)), ncol=2)


    contigsByLG <- sapply(orderedContigMatrix[,1], function(x){strsplit(x, '[.]')[[1]]})
    contigStarts <- sub('-.*', '', sub('.*:', '', orderedContigMatrix[,2]))
    orderedContigMatrix <- orderedContigMatrix[order(contigsByLG[1,], as.numeric(contigsByLG[2,] ), as.numeric(contigStarts)),]

    orderedContigMatrix <- ContigOrdering(orderedContigMatrix)

    if(class(mergedStrands) != 'matrix')
    {
      mergedStrands <- t(as.matrix(mergedStrands))
      rownames(mergedStrands) <- names(mergedContigs)
      mergedStrands <- StrandStateMatrix(mergedStrands)
    }else{
      rownames(mergedStrands) <- names(mergedContigs)
      mergedStrands <- StrandStateMatrix(mergedStrands)
    }
    return(list(mergedStrands=mergedStrands, contigKey=orderedContigMatrix))
  }

  if(is.null(whichLG)){whichLG=seq_len(length(linkageGroupList))}
  orderedGroups <- matrix(nrow=0, ncol=2)

  plotGroups <- vector("list", length(whichLG))
  for(lg in whichLG)
  {
    if(verbose){message(paste('-> Ordering fragments in LG', lg, sep=""))}
    if(length(linkageGroupList[[lg]]) > 1)
    {

      linkageGroup <- linkageGroupList[[lg]]
      linkageGroupReadTable <- strandStateMatrix[linkageGroup,]
      zeroGroups <- combineZeroDistContigs(linkageGroupReadTable, strandFreqMatrix, lg)
      #Make a contig Weight vector
      ordMat <- cbind(zeroGroups[[2]], apply(strandReadCount[which(rownames(strandReadCount) %in% zeroGroups[[2]][,2] ),] , 1, median))
      #Then make a LG weight by taking the sum of all contigs within that LG, and order the linkageGroupTable based on the deepest LG
      uniqueZeros <- unique(ordMat[,1])
      orderZeros <-sapply(uniqueZeros, function(x) sum(as.numeric(ordMat[which(ordMat[,1] == x),3])))
      orderZeros <- names(sort(orderZeros, decreasing=TRUE))
      #solves issue that subsetting one value matrix coreces it to a vector and causes crash
      if(nrow(zeroGroups[[1]]) > 1)
      {
        linkageGroupReadTable <- zeroGroups[[1]][orderZeros,]
      }else{
        linkageGroupReadTable <- zeroGroups[[1]]
      }
      #Choose which ordering method to use
      if(orderCall == 'greedy')
      {
        outOfOrder <- orderContigsGreedy(linkageGroupReadTable, randomAttempts=randomAttempts)
      }else if (orderCall == 'TSP')
      {
        outOfOrder <- orderContigsTSP(linkageGroupReadTable, reps=randomAttempts)
      }else{
        warning('###### WARNING! orderCall parameter not recognized! No ordering Performed. ######')
        break
      }

      mergedGroups <- matrix(nrow=0, ncol=2)
      for(gp in seq_len(length(outOfOrder[[1]]))) {
        mergedGroups <- rbind(mergedGroups, ordMat[which(ordMat[,1] == outOfOrder[[1]][gp]),1:2] )
      }
      orderedGroups <- rbind(orderedGroups, mergedGroups)

      plotGroups[[lg]] <- StrandStateMatrix(outOfOrder[[2]])
      if(!(is.null(saveOrdered)))
      {
        #find fragment names from lg
        chromosome <- sort(table(sapply(seq_len(length(linkageGroupList[[lg]])), function(x) strsplit(linkageGroupList[[lg]][x] ,':')[[1]][1])), decreasing=TRUE)
        #find predominant fragment name
        chromosome <- round(chromosome[1]/length(linkageGroupList[[lg]])*100, digits=1)

         plotFrame <- data.frame(outOfOrder[[2]])
         plotFrame <- data.frame(lapply(plotFrame, function(x) factor(x, levels=c(1,2,3))))
         rownames(plotFrame) <- rownames(outOfOrder[[2]])
         similarLinkageStrands <- as.matrix(1-daisy(plotFrame))
         diag(similarLinkageStrands) <- 1

         if(nrow(similarLinkageStrands) > 1)
         {
           breaks <- seq(0, 100, length.out=101)/100 
           cols <- colorRampPalette(c("cyan", "blue", "grey30", "black", "grey30", "red", "orange"))
           suppressWarnings(heatmap.2(similarLinkageStrands, 
                                      Rowv=NA, 
                                      Colv=NA, 
                                      dendrogram="none", 
                                      revC=TRUE, 
                                      col=cols(100), 
                                      breaks=breaks, 
                                      trace='none', 
                                      main=paste(orderCall, '-ordered LG', lg, '\n main fragment: ', names(chromosome), ' (', chromosome, '%)', sep="")))
         }
      }
    }else{
      groupName <-  paste('LG', lg, '.1', sep='')
      linkageGroup <- linkageGroupList[[lg]]
      orderedGroups <- rbind(orderedGroups, c(groupName, linkageGroup))
      linkageGroupReadTable <- t(as.matrix(strandStateMatrix[linkageGroup,]))
      rownames(linkageGroupReadTable) <- groupName
      plotGroups[[which(whichLG == lg)]] <- StrandStateMatrix(linkageGroupReadTable)
    }
  }  
  rownames(orderedGroups) <- NULL
  colnames(orderedGroups) <- c('LG', 'contig')
  rownames(orderedGroups) <- sapply(seq_len(nrow(orderedGroups)), function(x) strsplit(orderedGroups[x,], "[.]")[[1]][1])
  orderedGroups <- ContigOrdering(orderedGroups, reps)

  plotGroups <- StrandStateList(plotGroups, names=paste('LG', whichLG, ' StrandStateMatrix', sep=''))
  return(list(orderedGroups, plotGroups))
}

####################################################################################################
#' Function to call contig ordering algorithms iteratively across each linkage group element
#' @useDynLib contiBAIT
#' @param linkageGroupList list of vectors, each specifying which contigs belong in which linkage group (product of clusterContigs)
#' @param strandStateMatrix table of strand calls for all contigs (product of preprocessStrandTable)
#' @param strandFreqMatrix table of W:C read proportions (used for QC) (product of strandSeqFreqTable[[1]])
#' @param strandReadCount table of read counts (product of strandSeqFreqTable[[2]])  
#' @param whichLG vector of integers specifying the element(s) of linkageGroupList to be ordered (i.e. which specific linkage groups to try to order). Default is all LGs.
#' @param saveOrdered Will return a pdf of heatmaps for each linkage group; String entered becomes the fileName (default is saveOrderedPDF=FALSE)
#' @param orderCall currently either 'greedy' for greedy algorithm or 'TSP' for travelling salesperson alogrithm (default is 'greedy')
#' @param verbose Pringts messages to the terminal. Default is TRUE
#' @param randomAttempts itenger specifying number of iterations of either the Monte Carlo or the 2-opt TSP solver to identify the best ordering. Default is 75, but higher values (eg 1000) make more sense for the TSP.
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
