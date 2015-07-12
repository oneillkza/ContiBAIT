# Function to combine contigs that are right next to each other.
# @param linkageStrands correctly oriented strandTable, but only containing the rows for this linkage group
# @return a two-member list, first=new strand table, second=list mapping new merged contigs to old

combineZeroDistContigs <- function(linkageStrands, rawStrandTable)
{
	#Filter out borderline calls:
	
	contigNames <- rownames(linkageStrands)
	
	linkageMat <- as.matrix(linkageStrands)
	linkageMat <- apply(linkageMat, 2, as.numeric)
	
	rawStrandTable <- rawStrandTable[rownames(linkageStrands), colnames(linkageStrands)]
	linkageMat[which(abs(rawStrandTable) > 0.2 & abs(rawStrandTable) < 0.9 ) ] <- NA
	
	linkageStrands <- data.frame(linkageMat)
	
	linkageStrands <- data.frame(lapply(linkageStrands, as.factor))
	rownames(linkageStrands) <- contigNames
	
	
	##Combine zero dist contigs:
	strandDist <- daisy(linkageStrands)
	strandDist <- as.matrix(strandDist)
	
	mergedContigs <- list()
	beenMerged <- vector()
	mergedStrands <- matrix(nrow=0, ncol=ncol(linkageStrands))
	
	for(contig in rownames(linkageStrands))
	{
		#Only merge contigs not already pulled in by other merges:
		if(!contig %in% beenMerged)
		{
			toMerge <- which(strandDist[contig,] == 0)
			beenMerged <- append(beenMerged, names(toMerge))
			mergedContigs[[paste(names(toMerge), collapse='+')]] <- names(toMerge)
			mergedStrands <- rbind(mergedStrands, linkageStrands[contig,])
		}
	}
	rownames(mergedStrands) <- names(mergedContigs)
	
	return(list(mergedStrands=mergedStrands, contigKey=mergedContigs))
}