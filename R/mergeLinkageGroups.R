
mergeLinkageGroups.func <- function(object, allStrands, clusNum=1, cluster=1, similarityCutoff=0.7)
{

  consensusStrands <- data.frame(do.call(rbind, lapply(object, computeConsensus, allStrands)))
  consensusStrands <- data.frame(lapply(consensusStrands, function(x){factor(x, levels=c(1,2,3))}))  
  colnames(consensusStrands) <- colnames(allStrands)
  rownames(consensusStrands) <- names(object)

  consensusStrands <- new('StrandStateMatrix', consensusStrands)
  slaveNum <- makeCluster(clusNum)

  consensus.groups <- clusterContigs(consensusStrands, snowCluster=slaveNum, recluster=cluster, randomise=TRUE, similarityCutoff=similarityCutoff, verbose=FALSE)
  stopCluster(slaveNum)
 
  mergeThoseGroups <- lapply(consensus.groups, function(x) as.character(melt(object[x])[,1]))

  mergeThoseGroups <- new('LinkageGroupList', 
                          mergeThoseGroups[order(sapply(mergeThoseGroups, length), decreasing=TRUE)], 
                          names= sapply(1:length(mergeThoseGroups), function(x){paste('LG', x, ' (', length(mergeThoseGroups[[x]]), ')', sep='')}) )

  return(mergeThoseGroups)
}
####################################################################################################
#' mergeLinkageGroups -- merge very similar linkage groups, including those in reverse orientation
#' @param object LinkageGroupList 
#' @param allStrands StrandStateMatrix for all linkageGroups (usually reoriented by reorientStrandTable)
#' @param clusNum  Number of parallel processors to use when clustering contigs. Default is 1. 
#' @param cluster  Integer denoting the number of reclusterings to be performed for creating linkage groups (default is 1)
#' @param similarityCutoff merge contigs that are more similar this this
#' @return list of indices within the allStrands matrix indicating linkage group membership,  
#' a list consisting of a strandStateMatrix (a reoriented version of allStrands), 
#' and a data.frame of type OrientationFrame containing contig names and orientations, as '+' or '-'.
#' @aliases mergeLinkageGroups mergeLinkageGroups,LinkageGroupList,LinkageGroupList-method
#' @export
#' @import snow
#' @importFrom reshape2 melt 
#' @include AllClasses.R
####################################################################################################

setMethod('mergeLinkageGroups',
          signature = signature(object='LinkageGroupList', allStrands = 'StrandStateMatrix'),
          definition = mergeLinkageGroups.func
          )
