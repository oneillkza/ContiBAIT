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
#	if(!is.data.frame(allStrands))
#		stop('the allStrands variable passed to computeConsensus function must be a data frame.')

	strandVec <- allStrands[groupMembers,]

	if(length(groupMembers) > 1) {
		# score support
		qcScores <- apply(strandVec, 2, function(col) { sum(!is.na(col)) / length(groupMembers) })

		# make vector which, for each column, holds that column's most occuring state (1, 2, or 3)
		strandVec <- apply(strandVec, 2, function(col) { 
			which.max(sapply(1:3, function(i) { sum(col == i, na.rm=TRUE) })) })			
			
		# replace with NA those columns that were predominantly NA
		strandVec[qcScores < minSupport] <- NA
	}

	strandVec
}
