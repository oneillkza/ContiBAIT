####################################################################################################
#' Function to call contig ordering algorithms iteratively across each linkage group element
#' @param linkageGroupList list of vectors, each specifying which contigs belong in which linkage group (product of clusterContigs)
#' @param strandStateMatrix table of strand calls for all contigs (product of preprocessStrandTable)
#' @param strandFreqMatrix list with table of W:C read proportions (used for QC) and read counts (product of strandSeqFreqTable)
#' @param saveOrderedPDF Will return a pdf of heatmaps for each linkage group; String entered becomes the fileName (default is saveOrderedPDF=FALSE)
#' @param verbose Pringts messages to the terminal. Default is TRUE
#' 
#' @return a data.frame of ordered contigs with linkage group names
#' 
#' @export
#' @importFrom cluster daisy
#' @importFrom gplots heatmap.2
####################################################################################################


orderAllLinkageGroups <- function(linkageGroupList, strandStateMatrix, strandFreqMatrix, saveOrderedPDF=FALSE, verbose=TRUE)
{
  orderedGroups <- data.frame(LG=vector(), name=vector())
  if(saveOrderedPDF != FALSE) {pdf(paste(saveOrderedPDF, 'contig_order.pdf', sep='_'))}

  for( lg in seq(1, length(linkageGroupList)))
  {
    if(verbose){message(paste('  -> Ordering fragments in LG', lg, sep=""))}
    if(length(linkageGroupList[[lg]]) > 1)
    {
      if(orderCall == 'greedy')
      {
        outOfOrder <- orderContigsGreedy(linkageGroupList, strandStateMatrix, strandFreqMatrix, lg, randomAttempts=1)
      }else{
        outOfOrder <- orderContigsTSP(linkageGroupList, strandStateMatrix, strandFreqMatrix[[1]], lg)
      }
      orderFrame <- outOfOrder[[3]]
      orderedGroups <- rbind(orderedGroups, orderFrame)
      chromosome <- strsplit(linkageGroupList[[lg]][1],':')[[1]][1]
      if(saveOrderedPDF != FALSE)
      {
        similarLinkageStrands <- as.matrix(1-daisy(outOfOrder[[2]]))
        diag(similarLinkageStrands) <- 1
        if(nrow(similarLinkageStrands) > 1)
        {
          breaks <- seq(0, 100, length.out=101)/100 
          cols <- colorRampPalette(c("cyan", "blue", "grey30", "black", "grey30", "red", "orange"))
          suppressWarnings(heatmap.2(similarLinkageStrands, Rowv=NA, Colv=NA, dendrogram="none", col=cols(100), breaks=breaks, trace='none', main=paste('greedy-ordered ', chromosome, sep="")))
        }
      }
    }
  }
  if(saveOrderedPDF != FALSE){dev.off()}
  orderedGroups <- new("ContigOrdering", orderedGroups)

  return(orderedGroups)
}

