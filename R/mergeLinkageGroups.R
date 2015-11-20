
mergeLinkageGroups.func <- function(object, allStrands, similarityCutoff=0.7)
{
  linkageGroups <- object
  linkageStrands <- data.frame(do.call(rbind, lapply(linkageGroups, computeConsensus, allStrands)))
  rownames(linkageStrands) <- 1:nrow(linkageStrands)
  sim <- 1-as.matrix(daisy(data.frame(linkageStrands)))
  rownames(sim) <- rownames(linkageStrands)
  removedGroups <- vector()
  
  #Remove groups with NAs in the distance matrix:
  dist.nas <- apply(sim, 2, function(x){length(which(is.na(x)))})
  while(max(dist.nas) > 0)
  {
    sim <- sim[-c(which.max(dist.nas)), -c(which.max(dist.nas))]
    dist.nas <- apply(sim, 2, function(x){length(which(is.na(x)))})
    removedGroups <- append(removedGroups, which.max(dist.nas))
  }
  #TODO: removal not actually doing anything right now
  
  #Use hclust to merge remaining:
  linkageClust <- hclust(as.dist(1-sim))
  newLabs <- cutree(linkageClust, h=1-similarityCutoff)
  
  
  
  newGroups <- list()
  for (i in names(newLabs))
  {
    groupNum <- as.numeric(i)
    
    clustNum <- newLabs[i]
    
    if(clustNum < length(newGroups))	
    {			
      newGroups[[clustNum]] <- append(newGroups[[clustNum]], linkageGroups[[groupNum]])
    }
    else
    {
      newGroups[[clustNum]] <- linkageGroups[[groupNum]]
    }
  }	
  
  return(newGroups)
  
  
  
}
####################################################################################################
#' mergeLinkageGroups -- merge very similar linkage groups, including those in reverse orientation
#' @param object LinkageGroupList 
#' @param allStrands StrandStateMatrix for all linkageGroups (usually reoriented by reorientStrandTable)
#' @param similarityCutoff merge contigs that are more similar this this
#' @return list of indices within the allStrands matrix indicating linkage group membership, 
#'  vector of orientation changes
#' @aliases mergeLinkageGroups mergeLinkageGroups,LinkageGroupList,LinkageGroupList-method
#' @export
#' @importFrom cluster daisy
#' @include AllClasses.R
####################################################################################################

setMethod('mergeLinkageGroups',
          signature = signature(object='LinkageGroupList', allStrands = 'StrandStateMatrix'),
          definition = mergeLinkageGroups.func
          )