reorientLinkageGroups.func <- function(object, allStrands, verbose=TRUE)
{

	completeOrientation <- data.frame(contig=vector(), orientation=vector())

	for(lg in seq(1:length(object)))
	{
		if(verbose){message('Reorienting fragments from LG', lg, ' [', lg, '/', length(object), ']' )}

		linkageGroup <- object[[lg]]
		if(length(linkageGroup) > 1)
		{
			subsetStrands <- allStrands[which(rownames(allStrands) %in% linkageGroup),]
			subsetStrands <- replace(subsetStrands, subsetStrands == 2, NA)
			sim <- suppressWarnings(1-as.matrix(daisy(data.frame(subsetStrands))))
			sim[is.na(sim)] <- 0

			rownames(sim) <- rownames(subsetStrands)
			findGroups <- cutree(hclust(dist(sim)), h=10)
	  		getMax <- names(sort(table(findGroups), decreasing=TRUE))
	  		forwardStrands <- names(findGroups)[which(findGroups == getMax[1])]
			reverseStrands <- names(findGroups)[which(findGroups == getMax[2])]

			orientVec <- linkageGroup
			orientVec[which(orientVec %in% forwardStrands)] <- '+'
			orientVec[which(orientVec %in% reverseStrands)] <- '-'
			orientVec[which(!((orientVec == '+') | (orientVec == '-')))] <- '*'
			orientationFrame <- data.frame(contig=linkageGroup, orientation=orientVec)
		}else{
			orientationFrame <-  data.frame(contig=linkageGroup, orientation='+')
		}
		completeOrientation <- rbind(completeOrientation, orientationFrame)
	}

	toReorient <- as.character(completeOrientation$contig[which(completeOrientation$orientation == '-')])


	toReorientStrands <- allStrands[toReorient,]
	toReorientStrands <- data.frame(
		apply(toReorientStrands, c(1,2),
			  function(entry)
			  {
			  	if(!is.na(entry)&&entry=='3') return('1')
			  	if(!is.na(entry)&&entry=='1') return('3')
			  	entry
			  }
	))

	allStrands[toReorient,] <- toReorientStrands
	return(list(new('StrandStateMatrix', allStrands), completeOrientation))
}

####################################################################################################
#' reorientLinkageGroups uses a simple dissimilarity to find misoriented fragments within linkage groups.
#' @param object List of vectors containing names of contigs belonging to each LG.
#' @param allStrands Table of strand state for all contigs. Product of StrandSeqFreqTable.
#' @param verbose Outputs information to the terminal. Default is TRUE.
#' @aliases reorientLinkageGroups reorientLinkageGroups,LinkageGroupList,LinkageGroupList-method, strandStateMatrix, strandStateMatrix-method
#' @return a list consisting of a strandStateMatrix (a reoriented version of allStrands), and a data.frame of contig names and orientations, as '+' or '-'.
#' 
#' @export
####################################################################################################

setMethod('reorientLinkageGroups',
		  signature = signature(object='LinkageGroupList', allStrands = 'StrandStateMatrix'),
		  definition = reorientLinkageGroups.func
		  )
