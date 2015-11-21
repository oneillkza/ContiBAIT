####################################################################################################
#' reorientLinkageGroups uses a greedy algorithm to try to find the correct orientations of linkage groups
#' @param linkageGroups List of vectors containing names of contigs belonging to each LG.
#' @param allStrands Table of strand state for all contigs.
#' @param dissimilarityCutoff Try to reorient all LGs more dissimilar than this to any other LG. Default is 0.3
#' @param maxiter Maximum number of iterations to try reoriented for (to prevent possible infinite loops). Default is 100.
#' @param verbose Outputs information to the terminal. Default is TRUE.
#' @return a vector of orientations, as '+' or '-', in the order of linkageGroups.
#' 
#' @export
####################################################################################################

reorientLinkageGroups <- function(linkageGroups, allStrands, dissimilarityCutoff=0.3, maxiter=100, verbose=TRUE)
{
	linkageStrands <- data.frame(do.call(rbind, lapply(linkageGroups, computeConsensus, allStrands)))
	orientation <- rep('+', nrow(linkageStrands))
	rownames(linkageStrands) <- 1:nrow(linkageStrands)
	sim <- 1-as.matrix(daisy(data.frame(linkageStrands)))
	rownames(sim) <- rownames(linkageStrands)
	
	#Greedy algorithm to try to reorient linkage groups.
	#Picks the most misoriented-looking of the contis with one or more apparent misorientation
	for(iteration in 1:maxiter)
	{
		#Compute similarity matrix:
		sim <- 1-as.matrix(daisy(data.frame(linkageStrands)))
		rownames(sim) <- rownames(linkageStrands)
		
		#If no LGs need reorientation, then we're done
		if(min(sim[which(!is.na(sim))]) > dissimilarityCutoff)
			break
		
		#Figure out which LGs need reorientation, and score them, then pick one:
		groupNeedsInversion <- apply(sim, 1, function(x){any(x < (dissimilarityCutoff))})
		groupScores <- apply(sim, 1, sum, na.rm=T) / nrow(sim)
		toInvert <- names(which.min(groupScores[which(groupNeedsInversion)]))[1]
		
		#Invert LG and record:
		toInvertStrands <- linkageStrands[toInvert, ]
		suppressWarnings(toInvertStrands[which(toInvertStrands==3)] <- 2)
		suppressWarnings(toInvertStrands[which(toInvertStrands==1)] <- 3)
		suppressWarnings(toInvertStrands[which(toInvertStrands==2)] <- 1)
		suppressWarnings(linkageStrands[toInvert, ] <- toInvertStrands)
		
		if(verbose){message('Reorienting LG ', toInvert)}
		
		if(orientation[as.numeric(toInvert)]=='+') 
		{
			suppressWarnings(orientation[as.numeric(toInvert)] <- '-')
		}else
		{
			suppressWarnings(orientation[as.numeric(toInvert)] <- '+')
		}
		if(iteration==maxiter)
			warning('Maximum iterations reached while reorienting without convergence.')
	}	
	return(orientation)
}