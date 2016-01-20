reorientLinkageGroups.func <- function(object, allStrands, verbose=TRUE)
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

	completeOrientation <- data.frame(contig=vector(), orientation=vector())

	#find consensus
	linkageStrands <- data.frame(do.call(rbind, lapply(object, computeConsensus, allStrands)))
	linkageStrands <- switchAroo(linkageStrands)
	colnames(linkageStrands) <- colnames(allStrands) 
	counter=1

	for(lg in names(object))
	{
		if(verbose){message('Reorienting fragments from ', lg, ' [', counter, '/', length(object), ']' )}

		linkageGroup <- object[[lg]]
		if(length(linkageGroup) > 1)
		{
			subsetStrands <- allStrands[which(rownames(allStrands) %in% linkageGroup),]

			#add dummy opposite line to subsetStrands
			subsetStrands <- rbind(subsetStrands, linkageStrands[which(rownames(linkageStrands) == lg),])

			subsetStrands <- replace(subsetStrands, subsetStrands == 2, NA)
			sim <- suppressWarnings(1-as.matrix(daisy(data.frame(subsetStrands))))
			sim[is.na(sim)] <- 0

			rownames(sim) <- rownames(subsetStrands)
			findGroups <- cutree(hclust(dist(sim)), k=2)
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
		counter <- counter+1
	}

	toReorient <- as.character(completeOrientation$contig[which(completeOrientation$orientation == '-')])

	toReorientStrands <- switchAroo(allStrands[toReorient,])

	allStrands[toReorient,] <- toReorientStrands
	rownames(completeOrientation) <- completeOrientation$contig
	return(list(new('StrandStateMatrix', allStrands), new('OrientationFrame', completeOrientation)))
}

####################################################################################################
#' reorientLinkageGroups uses a simple dissimilarity to find misoriented fragments within linkage groups.
#' @param object List of vectors containing names of contigs belonging to each LG.
#' @param allStrands Table of strand state for all contigs. Product of StrandSeqFreqTable.
#' @param verbose Outputs information to the terminal. Default is TRUE.
#' @aliases reorientLinkageGroups reorientLinkageGroups,LinkageGroupList,LinkageGroupList-method, strandStateMatrix, strandStateMatrix-method
#' @return a list consisting of a strandStateMatrix (a reoriented version of allStrands), and a data.frame of type OrientationFrame containing contig names and orientations, as '+' or '-'.
#' 
#' @export
####################################################################################################

setMethod('reorientLinkageGroups',
		  signature = signature(object='LinkageGroupList', allStrands = 'StrandStateMatrix'),
		  definition = reorientLinkageGroups.func
		  )
