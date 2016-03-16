reorientLinkageGroups.func <- function(object, 
									   allStrands, 
									   previousOrient=NULL, 
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

	completeOrientation <- matrix(nrow=0, ncol=2)
	colnames(completeOrientation) <- c('contig', 'orientation')
	rememberRows <- rownames(allStrands)
	
	allStrands <-  data.frame(lapply(data.frame(allStrands), 
									 function(x){factor(x, levels=c(1,2,3))})) 
	rownames(allStrands) <- rememberRows
  
	#find consensus
	linkageStrands <- data.frame(do.call(rbind, lapply(object, computeConsensus, allStrands)))
	#find the 'opposite' of consensus to append to each LG for downstream cutree
	linkageStrands <- switchAroo(linkageStrands)
	colnames(linkageStrands) <- colnames(allStrands) 
	counter=1

	for(lg in names(object))
	{
		if(verbose){message('Reorienting fragments from ', lg, 
							' [', counter, '/', length(object), ']' )}

		linkageGroup <- object[[lg]]
		if(length(linkageGroup) > 1)
		{
			subsetStrands <- allStrands[which(rownames(allStrands) %in% linkageGroup),]

			#add dummy opposite line to subsetStrands
			subsetStrands <- rbind(subsetStrands, 
								   linkageStrands[which(rownames(linkageStrands) == lg),])
			subsetStrands <- replace(subsetStrands, subsetStrands == 2, NA)
      
			sim <- suppressWarnings(1-as.matrix(daisy(data.frame(subsetStrands))))
			sim[is.na(sim)] <- 0

			findGroups <- cutree(hclust(dist(sim)), k=2)
	  		getMax <- names(sort(table(findGroups), decreasing=TRUE))
	  		forwardStrands <- names(findGroups)[which(findGroups == getMax[1])]
			reverseStrands <- names(findGroups)[which(findGroups == getMax[2])]

			orientVec <- linkageGroup
			orientVec[which(orientVec %in% forwardStrands)] <- '+'
			orientVec[which(orientVec %in% reverseStrands)] <- '-'
			orientVec[which(!((orientVec == '+') | (orientVec == '-')))] <- '*'
			orientationFrame <- matrix(c(linkageGroup, orientVec), ncol=2)
		}else{
			orientationFrame <-  matrix(c(linkageGroup, '+'), ncol=2)
		}
		completeOrientation <- rbind(completeOrientation, orientationFrame)
		counter <- counter+1
	}

	toReorient <- as.character(completeOrientation[which(completeOrientation[,2] == '-'),1])

	toReorientStrands <- switchAroo(allStrands[toReorient,])

	if(!is.null(previousOrient))
	{
	  toInvert <- previousOrient[which(previousOrient[,1] %in% toReorient),2]
	  levels(toInvert) <- rev(levels(toInvert))
	  previousOrient[which(previousOrient[,1] %in% toReorient),2] <- toInvert
	  completeOrientation <- previousOrient
	}

	allStrands[toReorient,] <- toReorientStrands
	rownames(completeOrientation) <- completeOrientation[,1]
	allStrands <- data.matrix(allStrands)


	return(list(StrandStateMatrix(allStrands), OrientationFrame(completeOrientation)))
}

####################################################################################################
#' reorientLinkageGroups uses a simple dissimilarity to find misoriented fragments within linkage groups.
#' @param object List of vectors containing names of contigs belonging to each LG.
#' @param allStrands Table of type \code{strandStateMatrix} encompassing strand state for all contigs. Product of StrandSeqFreqTable.
#' @param previousOrient data.frame of type \code{OrientationFrame} of previous orientation states if performed. Default is NULL
#' @param verbose Outputs information to the terminal. Default is TRUE.
#' @aliases reorientLinkageGroups reorientLinkageGroups-LinkageGroupList-StrandStateMatrix-method
#' @rdname reorientLinkageGroups
#' @example inst/examples/reorientLinkageGroups.R
#' @return a list consisting of a strandStateMatrix (a reoriented version of allStrands), and a data.frame of type OrientationFrame containing contig names and orientations, as '+' or '-'.
#' 
#' @export
####################################################################################################

setMethod('reorientLinkageGroups',
		  signature = signature(object='LinkageGroupList', 
		  					  allStrands = 'StrandStateMatrix'),
		  definition = reorientLinkageGroups.func
		  )
