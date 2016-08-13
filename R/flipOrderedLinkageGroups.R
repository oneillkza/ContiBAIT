flipOrderedLinkageGroups.func <- function(contigOrdering, orderFrame, linkageGroupList, strandStateMatrix, whichLG=NULL, maxiter=100, dissimilarityCutoff=0.3)
{
    exdiag <- function(mat, off){mat[row(mat)+off == col(mat)]}

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


    flipperMeTimbers <- function(groupCheck, its=maxiter, cutOff=dissimilarityCutoff)
    {

		checkForMiso <- clusterContigs(groupCheck, clusterBy='homo', verbose=FALSE)
	   	consensusTable <- data.frame(do.call(rbind, lapply(checkForMiso, computeConsensus, groupCheck)))
		consensusTable <- replace(consensusTable, consensusTable == 2, NA) 
		consensusTable <- consensusTable[colSums(!is.na(consensusTable)) > 0] 	
		orientation <- rep('+', nrow(consensusTable))
			
		#Greedy algorithm to try to reorient linkage groups.
		#Picks the most misoriented-looking of the contis with one or more apparent misorientation
		for(iteration in 1:its)
		{
			#Compute similarity matrix:
			sim <- suppressWarnings(1-as.matrix(daisy(data.frame(consensusTable))))
			#rownames(sim) <- 1:nrow(sim)
			
			#If no LGs need reorientation, then we're done
			if(min(sim[which(!is.na(sim))]) > cutOff)
				break
			
			#Figure out which LGs need reorientation, and score them, then pick one:
			groupNeedsInversion <- apply(sim, 1, function(x){any(x < (cutOff))})
			groupScores <- apply(sim, 1, sum, na.rm=T) / nrow(sim)
			toInvert <- names(which.min(groupScores[which(groupNeedsInversion)]))[1]
				
			#Invert LG and record:
			toInvertStrands <- consensusTable[toInvert, ]
			suppressWarnings(toInvertStrands[which(toInvertStrands==3)] <- 2)
			suppressWarnings(toInvertStrands[which(toInvertStrands==1)] <- 3)
			suppressWarnings(toInvertStrands[which(toInvertStrands==2)] <- 1)
			suppressWarnings(consensusTable[toInvert, ] <- toInvertStrands)

			rowNum <- match(toInvert, rownames(consensusTable))
			if(orientation[rowNum]=='+') 
			{
				suppressWarnings(orientation[rowNum] <- '-')
			}else
			{
				suppressWarnings(orientation[rowNum] <- '+')
			}
			if(iteration==maxiter){
				warning('Maximum iterations reached while reorienting without convergence.')
			}
		}	
		#Now flip the reads that need flipping.
		contigsToFlip <-  unlist(checkForMiso[rownames(consensusTable)[which(orientation == '-')]])
		names(contigsToFlip) <- NULL
		return(contigsToFlip)
	}    

	# CODE STARTS HERE. TAKE ALL LGs WITH >1 ELEMENTS
	if(is.null(whichLG)){whichLG=length(which(sapply(seq_len(length(orderFrame)), function(x) nrow(orderFrame[[x]])) > 1))}

 	rownames(contigOrdering) <- contigOrdering[,1]

	for(lg in seq_len(whichLG))
  	{
  		orderGroup <- orderFrame[[lg]]
  		chrNames <- rownames(orderGroup)
		lgOrdering <- contigOrdering[which(contigOrdering[,1] %in% rownames(orderGroup)),]

		#run flipper here.
		misoContig <- flipperMeTimbers(orderGroup)
		if(!(is.null(misoContig)))
		{
			fetchNames <- lgOrdering[which(lgOrdering[,1] %in% misoContig),2]
			#now convert the reoriented table to be returned
			toReorientMatrix <- switchAroo(strandStateMatrix[fetchNames,,drop=FALSE])
			reorientedMatrixElements <-  as.matrix(data.frame(lapply(toReorientMatrix, 
					 								   function(x){as.numeric(as.character(x))}))) 
			strandStateMatrix[fetchNames,] <- reorientedMatrixElements
		}
	}
		return(strandStateMatrix)
}


####################################################################################################
#' flipOrderedLinkageGroups -- Simple greedy algorithm to find misorientations within previously ordered LG.
#'
#' @param contigOrdering a data.frame of ordered contigs with linkage group names of class ContigOrdering
#' @param orderFrame a list of StrandStateMatrix elements of class StrandStateList
#' @param linkageGroupList List of vectors containing names of contigs belonging to each LG of type LinkageGroupList. Should be ordered as the output of orderAllLinkageGroups[[3]]
#' @param strandStateMatrix a StrandStateMatrix to be reoriented
#' @param whichLG  vector of integers specifying the element(s) of linkageGroupList to be ordered (i.e. which specific linkage groups to try to order). Default is all LGs.
#' @param maxiter The number of iterations used to find the optimal path. Default is 100 
#' @param dissimilarityCutoff The minimum amount of dissimilarity needed to start reorienting contigs. Default is 0.3.
#' @aliases flipOrderedLinkageGroups flipOrderedLinkageGroups-ContigOrdering-StrandStateList-method
#' @rdname flipOrderedLinkageGroups
#' @return a StrandStateMatrix with misorientations resolved
#' 
#' @importFrom cluster daisy
#' @export
#' @include AllClasses.R
####################################################################################################

 setMethod('flipOrderedLinkageGroups',
 		  signature = signature(contigOrdering='ContigOrdering',
 		  						orderFrame='StrandStateList',
 		  						linkageGroupList='LinkageGroupList',
 		  						strandStateMatrix='StrandStateMatrix'),
 		  definition = flipOrderedLinkageGroups.func
 		  )
