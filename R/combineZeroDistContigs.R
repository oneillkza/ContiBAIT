#' Function to combine contigs that are right next to each other.
#' @param linkageStrands correctly oriented strandTable, but only containing the rows for this linkage group
#' @return a two-member list, first=new strand table, second=list mapping new merged contigs to old
#' @include AllClasses.R
#' @importFrom cluster daisy

combineZeroDistContigs <- function(linkageStrands, rawStrandTable, lg)
{
	#Filter out borderline calls:
	
	contigNames <- rownames(linkageStrands)
	
	linkageMat <- as.matrix(linkageStrands)
	linkageMat <- apply(linkageMat, 2, as.numeric)
	
	rawStrandTable <- rawStrandTable[rownames(linkageStrands), ]
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
	groupCount <- 1
	
	for(contig in rownames(linkageStrands))
	{
		#Only merge contigs not already pulled in by other merges:
		if(!contig %in% beenMerged)
		{
			toMerge <- which(strandDist[contig,] == 0)
			#And don't pull in contigs that are present in toMerge if they've already been asigned
			toMerge <- toMerge[!names(toMerge) %in% beenMerged]
			beenMerged <- append(beenMerged, names(toMerge))
#			mergedContigs[[paste(names(toMerge), collapse='+')]] <- names(toMerge)
			mergedContigs[[paste('LG', lg, '.', groupCount, sep='')]] <- names(toMerge)
			mergedStrands <- rbind(mergedStrands, linkageStrands[contig,])
			groupCount <- groupCount +1

		}
	}
	
	orderedContigMatrix <- data.frame(LG=unlist(sapply(1:length(mergedContigs), function(x) rep(names(mergedContigs[x]), length(mergedContigs[[x]]) ))), contig=unlist(mergedContigs), row.names=NULL, stringsAsFactors=FALSE )
	orderedContigMatrix <- new("ContigOrdering", orderedContigMatrix)

	mergedStrands <- data.frame(lapply(mergedStrands, function(x){factor(x, levels=c(1,2,3))}))  
	rownames(mergedStrands) <- names(mergedContigs)
	mergedStrands <- new("StrandStateMatrix", mergedStrands)

	return(list(mergedStrands=mergedStrands, contigKey=orderedContigMatrix))
}