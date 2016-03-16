
mergeLinkageGroups.func <- function(object, 
									allStrands, 
									clusterParam=NULL, 
									cluster=1, 
									similarityCutoff=0.7)
{

  consensusStrands <- data.matrix(do.call(rbind, 
  										lapply(object, computeConsensus, allStrands)))

  consensusStrands <- StrandStateMatrix(consensusStrands)

  consensus.groups <- clusterContigs(consensusStrands, 
  								   clusterParam=clusterParam, 
  								   recluster=cluster, 
  								   randomise=TRUE, 
  								   similarityCutoff=similarityCutoff, 
  								   verbose=FALSE)

 
  mergeThoseGroups <- lapply(consensus.groups, function(x) as.character(melt(object[x])[,1]))

  mergeThoseGroups <- LinkageGroupList(
                          mergeThoseGroups[order(sapply(mergeThoseGroups, length), decreasing=TRUE)], 
                          names= sapply(1:length(mergeThoseGroups), 
                          			  function(x)
                          			  	{
                          			  	paste('LG', x, ' (', length(mergeThoseGroups[[x]]), ')', sep='')
                          			  	}))

  return(mergeThoseGroups)
}
####################################################################################################
#' mergeLinkageGroups -- merge very similar linkage groups, including those in reverse orientation
#' @param object LinkageGroupList 
#' @param allStrands StrandStateMatrix for all linkageGroups (usually reoriented by reorientStrandTable)
#' @param clusterParam  optional \code{BiocParallelParam} specifying cluster to use for parallel execution.
#' When \code{NULL}, execution will be serial.
#' @param cluster  Integer denoting the number of reclusterings to be performed for creating linkage groups (default is 1)
#' @param similarityCutoff merge contigs that are more similar this this
#' @return list of indices within the allStrands matrix indicating linkage group membership,  
#' a list consisting of a strandStateMatrix (a reoriented version of allStrands), 
#' and a data.frame of type OrientationFrame containing contig names and orientations, as '+' or '-'.
#' @aliases mergeLinkageGroups mergeLinkageGroups-LinkageGroupList-StrandStateMatrix-method
#' @rdname mergeLinkageGroups
#' @example inst/examples/mergeLinkageGroups.R
#' @export
#' @import BiocParallel
#' @importFrom reshape2 melt 
#' @include AllClasses.R
####################################################################################################

setMethod('mergeLinkageGroups',
          signature = signature(object='LinkageGroupList', allStrands = 'StrandStateMatrix'),
          definition = mergeLinkageGroups.func
          )
