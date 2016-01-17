orderAllLinkageGroups.func <- function(linkageGroupList, strandStateMatrix, strandFreqMatrix, strandReadCount, whichLG=NULL, saveOrdered=FALSE, orderCall='greedy', randomAttempts=75, verbose=TRUE)
{
  if(is.null(whichLG)){whichLG=c(1:length(linkageGroupList))}
#  if(saveOrderedPDF != FALSE) {pdf(paste(saveOrderedPDF, 'contig_order.pdf', sep='_'))}
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

#      orderedGroups <- lapply(1:length(outOfOrder[[1]]), function(group){zeroGroups[[2]][which(zeroGroups[[2]] == outOfOrder[[1]][group]),1:2]})
#      orderedGroups <- do.call(rbind, orderedGroups)
      
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
#  if(saveOrderedPDF != FALSE){dev.off()}
  orderedGroups <- new("ContigOrdering", orderedGroups)
  return(orderedGroups)
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
#' @param randomAttempts iterger specifying number of randomized clusterings to identify the best ordering. Default is 75
#' @aliases orderAllLinkageGroups orderAllLinkageGroups,LinkageGroupList,LinkageGroupList-method, StrandStateMatrix, StrandStateMatrix-method, StrandFreqMatrix, StrandFreqMatrix-method, StrandReadMatrix, StrandReadMatrix-method
#' 
#' @return a data.frame of ordered contigs with linkage group names
#' 
#' @example inst/examples/orderAllLinkageGroups.R
#' @export
#' @importFrom cluster daisy
#' @importFrom gplots heatmap.2
####################################################################################################

setMethod('orderAllLinkageGroups',
      signature = signature(linkageGroupList='LinkageGroupList', strandStateMatrix= 'StrandStateMatrix', strandFreqMatrix='StrandFreqMatrix', strandReadCount='StrandReadMatrix'),
      definition = orderAllLinkageGroups.func
      )
