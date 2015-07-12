####################################################################################################
# computeConsensus -- Function to compute the consensus strand inheritance of all members of
#   a linkage group
# @param groupMembers members within a single linkage group list element
# @param allStrands data.frame of all contigs
# @param misSupport the minimal value that supports the call (default is 0.05)
# @return A vector containing the contig to be added to the best matched LG
#  
#
####################################################################################################



computeConsensus <- function(groupMembers, allStrands, minSupport=0.05)
{
	#counter <- 1
	if(!is.data.frame(allStrands))
		stop('the allStrands variable passed to computeConsensus function must be a data frame.')
	
	
	if(length(groupMembers) > 1)
	{
		groupStrands <- allStrands[groupMembers,]
		tables <- table(groupStrands[,1])
		
		for (i in 2:ncol(groupStrands))
		{
			newTable <- table(groupStrands[,i])
			if(length(newTable) == 0) newTable <- c(0)
			tables <- cbind(tables, table(groupStrands[,i]))
		}
		
		strandVec <- apply(tables, 2, function(x){names(which.max(x))})
		qcScores <- apply(groupStrands, 2, function(x){length(which(!is.na(x)))}) / nrow(groupStrands)
		strandVec[which(qcScores<minSupport)] <- NA
		#scores <- apply(tables, 2, function(x){})
		#counter <- counter + 1
	}
	else
	{
		strandVec <- allStrands[groupMembers,]
	}
	strandVec
}
