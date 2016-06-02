mergeFlankedLGs.func <- function(linkageGroupList, 
                                          strandStateMatrix, 
                                          buildConsensus=1, 
                                          cluster=NULL, 
                                          clusterParam=NULL, 
                                          verbose=TRUE)
{
  switchAroo <- function(dataFrameToSwitch)
  {
    data.frame(
      apply(dataFrameToSwitch, c(1,2),
      function(entry)
      {
        if(!is.na(entry)&&entry=='3') return('1')
        if(!is.na(entry)&&entry=='1') return('3')
        entry
      }
     ))
  }

  #Find those elements greater than twice the build consensus (so flanking regions aren't made from the same contigs)
  buildOver <- sapply(seq_len(length(linkageGroupList)), function(x) length(linkageGroupList[[x]]) > (buildConsensus*2))
  includeCon <- linkageGroupList[buildOver]
  #Use first contig name from LG as an identified. We'll use this to merge groups back later.
  consensusNames <- sapply(includeCon, head, 1)
  #calculate consensus for upstream flank of LG
  upCont <- lapply(seq_len(length(includeCon)), function(x) head(includeCon[[x]], buildConsensus))
  consensusTableUp <- data.frame(do.call(rbind, lapply(upCont, computeConsensus, strandStateMatrix))) 
  #Calculate consensus for downstream flank of LG
  downCont <- lapply(seq_len(length(includeCon)), function(x) tail(includeCon[[x]], buildConsensus))
  consensusTableDown <- data.frame(do.call(rbind, lapply(downCont, computeConsensus, strandStateMatrix))) 
  #merge these together
  conTab <- rbind(consensusTableUp, consensusTableDown)
  colnames(conTab) <- colnames(strandStateMatrix)
 
  #Those contigs with not enough elements to make a flank, just compute consensus
  excludeCon <- linkageGroupList[buildOver == FALSE]
  consensusNamesEx <- sapply(excludeCon, head, 1)
  nameIndex <- c(consensusNames, consensusNamesEx)

  consensusTableExclude <- data.frame(do.call(rbind, lapply(excludeCon, computeConsensus, strandStateMatrix))) 
  rownames(consensusTableExclude) <- NULL
  colnames(consensusTableExclude) <- colnames(strandStateMatrix)
  justFlankMatrix <- StrandStateMatrix(as.matrix(rbind(conTab, consensusTableExclude)))
  #Make a ket so the rownames can be traced back to the contig names/LGs
  flankKey <- data.frame(clusterName=rownames(justFlankMatrix), 
                           contigName=c(consensusNames, consensusNames, consensusNamesEx), 
                           location=c(rep('up', length(consensusNames)), rep('down', length(consensusNames)), rep('all', length(consensusNamesEx)) ))

  #Cluster these consensus regions to see if LGs need to be merged
  linkageflank <- clusterContigs(justFlankMatrix, 
                                    clusterParam=clusterParam, 
                                    recluster=cluster, 
                                    similarityCutoff=0.9,
                                    minimumLibraryOverlap=10,
                                    randomise=TRUE)

  #take only those groups with more than one member (ie these should be merged)
  linkClusters <- sapply(seq_len(length(linkageflank)), function(x) length(linkageflank[[x]]) > 1)
  linkageflank <- linkageflank[linkClusters] 

  flankCount <- seq_len(length(linkageflank))
  if(verbose){message(' -> ', length(flankCount), ' mergable LGs found!')}
  #Now check orientation is ok...
  toReorient <- vector()
  for(listElement in flankCount)
  {
    if(verbose){message(' -> checking orientation for group ', listElement, ' of ', length(flankCount))}
    clus <- linkageflank[[listElement]]
    thisElement <- StrandStateMatrix(justFlankMatrix[clus,])
    orientClusters <- clusterContigs(thisElement, 
                                      clusterParam=clusterParam, 
                                      recluster=cluster, 
                                      similarityCutoff=0.9,
                                      clusterBy='homo',
                                      randomise=TRUE, 
                                      verbose=FALSE)
    #If it falls nicely into two groups, take the smaller group to reorient
    if(length(orientClusters) == 2)
    {
      flipGroup <- as.character(flankKey[match(orientClusters[[2]], flankKey$clusterName), 2])
      flippedIndex <- names(nameIndex[which(nameIndex %in% toReorient)])
      flipNames <- unlist(linkageGroupList[names(linkageGroupList) %in% flippedIndex])
      flippedStateMatrix  <- switchAroo(strandStateMatrix[flipNames,])
      flippedStateMatrix <-  as.matrix(data.frame(lapply(data.frame(flippedStateMatrix), 
                                                          function(x){as.numeric(x)}))) 
      strandStateMatrix[flipNames,] <- flippedStateMatrix
    }
  }

  #extract names of the contigs from the contig key
  linkageflank <- lapply(flankCount, function(x) as.character(flankKey[match(linkageflank[[x]], flankKey$clusterName), 2]))

  #and remove elements where both upstream and downstream merged...
  linkageflank <- lapply(flankCount, function(x) unique(linkageflank[[x]]))

  #and remove any cluster where the same contig is represented more than one LG...
  duplicateElements <- sort(unlist(linkageflank))
  uniqueElements <- duplicateElements[!(duplicated(duplicateElements) | duplicated(duplicateElements, fromLast=TRUE))]
  linkageflank <- lapply(flankCount, function(x) linkageflank[[x]][linkageflank[[x]] %in% uniqueElements])
 
  #take only those groups with more than one member (ie these should be merged)
  linkClusters <- sapply(seq_len(length(linkageflank)), function(x) length(linkageflank[[x]]) > 1)
  linkageflank <- linkageflank[linkClusters] 


  #take only those groups with more than one member (ie these should be merged)
  linkClusters <- sapply(seq_len(length(linkageflank)), function(x) length(linkageflank[[x]]) > 1)
  linkageflank <- linkageflank[linkClusters] 

  groupsToMerge <- lapply(seq_len(length(linkageflank)), function(x) names(nameIndex[which(nameIndex %in% linkageflank[[x]])]))

  revisedList <- linkageGroupList[!(names(linkageGroupList) %in% unlist(groupsToMerge))]
  mergeList <- lapply(seq_len(length(groupsToMerge)), function(x) unname(unlist(linkageGroupList[which(names(linkageGroupList) %in% groupsToMerge[[x]])])))
  revisedList <- c(revisedList, mergeList)
  revisedList <- revisedList[order(sapply(revisedList, length), decreasing=TRUE)]
  revisedList <- LinkageGroupList(
                          revisedList, 
                          names= sapply(1:length(revisedList), 
                                  function(x)
                                    {
                                    paste('LG', x, ' (', length(revisedList[[x]]), ')', sep='')
                                    }))
  return(list(revisedList, strandStateMatrix))
}

####################################################################################################
#' mergeFlankedLGs -- searches for similarities at the ends of ordered linkage groups to chain groups together
#'
#' @param linkageGroupList List of ordered vectors containing names of contigs belonging to each LG, of type LinkageGroupList
#' @param strandStateMatrix Table of type strandStateMatrix encompassing strand state for all contigs. Product of StrandSeqFreqTable.
#' @param buildConsensus number of contigs to take at the end of the linkage group to build a consensus strand state. Default is 5
#' @param cluster Number of times to recluster and take the consensus of. If NULL, clustering is 
#' run only once.
#' @param clusterParam optional \code{BiocParallelParam} specifying cluster to use for parallel execution. When NULL, execution will be serial.
#' @param verbose Outputs information to the terminal. Default is TRUE.
#' @return a list containing a revised LinkageGroupList with merged groups, if appropriate, and 
#' a StrandStateMatrix with contigs reoriented, if newly merged groups were in opposite orientations.
#' @aliases mergeFlankedLGs mergeFlankedLGs,LinkageGroupList,StrandStateMatrix-method
#' @rdname mergeFlankedLGs
#' 
#' @importFrom cluster daisy
#' @export
#' @include AllClasses.R
####################################################################################################

 setMethod('mergeFlankedLGs',
      signature = signature(linkageGroupList='LinkageGroupList', strandStateMatrix='StrandStateMatrix'),
      definition = mergeFlankedLGs.func
      )