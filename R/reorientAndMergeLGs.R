reorientAndMergeLGs.func <- function(object, 
									   allStrands,
									   cluster=NULL, 
									   clusterParam=NULL,
									   similarityCutoff=0.9,
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

 	findConsensus <- function(linkList, strandMat, flip=FALSE)
 	{
		#find consensus
		consensusTable <- data.frame(do.call(rbind, lapply(linkList, computeConsensus, strandMat)))
		#find the 'opposite' of consensus to append to each LG 
		if(flip)
		{
			consensusTable <- switchAroo(consensusTable)
		}
		rememberRows <- rownames(consensusTable)
		consensusTable <-  as.matrix(data.frame(lapply(data.frame(consensusTable), 
					 								   function(x){as.numeric(x)}))) 
		rownames(consensusTable) <- rememberRows
		return(consensusTable)
	}

	completeOrientation <- matrix(nrow=0, ncol=2)
	colnames(completeOrientation) <- c('contig', 'orientation')
	
	counter=1

	linkageStrands <- findConsensus(object, allStrands, flip=TRUE)

	for(lg in names(object))
	{
		if(verbose){message('Reorienting fragments from ', lg, 
							' [', counter, '/', length(object), ']' )}

		linkageGroup <- object[[lg]]
		if(length(linkageGroup) > 1)
		{
			subsetStrands <- allStrands[linkageGroup,]
			#add dummy opposite line to subsetStrands
			tempNames <- rownames(subsetStrands)
			subsetStrands <- rbind(subsetStrands, linkageStrands[lg,])
			rownames(subsetStrands) <- c(tempNames, "OPPOSITE")
			subsetStrands <- StrandStateMatrix(subsetStrands)
      
		    linkedContigs <- clusterContigs(subsetStrands, 
		                                        clusterParam=clusterParam, 
		                                        similarityCutoff=0.7,
		                                        recluster=cluster,
		                                        clusterBy='homo',
		                                        randomise=TRUE,
		                                        verbose=FALSE)

			#find location of list element that contains the "OPPOSITE" strands. These are our misorientations
			oppoLoco <- which(sapply(seq_along(linkedContigs), function(x) which("OPPOSITE" %in% unlist(linkedContigs[[x]])) ) == 1)
   			reverseStrands <- linkedContigs[[oppoLoco]]
			forwardStrands <- as.character(unlist(linkedContigs))
			forwardStrands <- forwardStrands[which(!(forwardStrands %in% reverseStrands))]

			orientVec <- linkageGroup
			orientVec[which(orientVec %in% forwardStrands)] <- '+'
			orientVec[which(orientVec %in% reverseStrands)] <- '-'
			orientationFrame <- matrix(c(linkageGroup, orientVec), ncol=2)
		}else{
			orientationFrame <-  matrix(c(linkageGroup, '+'), ncol=2)
		}
		completeOrientation <- rbind(completeOrientation, orientationFrame)
		counter <- counter+1

	if(verbose){message(' -> ', nrow(orientationFrame[which(orientationFrame[,2] == '-'),]), 
							' fragments reoriented from ', length(linkageGroup), ' within ', lg )}

	}

	toReorient <- as.character(completeOrientation[which(completeOrientation[,2] == '-'),1])

	toReorientStrands <- switchAroo(allStrands[toReorient,])
	toReorientStrands <-  as.matrix(data.frame(lapply(data.frame(toReorientStrands), 
				 								   function(x){as.numeric(x)}))) 

	allStrands[toReorient,] <- toReorientStrands
	rownames(completeOrientation) <- completeOrientation[,1]

	#Now the reorientation is done, MERGE reoriented linkage groups by checking similarity (note need to redo computeConsensus as allStrands has changed)

 	consensusStrands <- findConsensus(object, allStrands)
	consensusStrandsOp <- findConsensus(object, allStrands, flip=TRUE)

  	consensusStrands <- StrandStateMatrix(consensusStrands)

  	#See if any LGs now merge together during regular clustering
  	consensus.groups <- clusterContigs(consensusStrands, 
	  								   clusterParam=clusterParam, 
	  								   recluster=cluster, 
	  								   randomise=TRUE, 
	  								   similarityCutoff=similarityCutoff, 
	  								   verbose=FALSE)

  	#All those that do, check to see if there are misorientations
  	subsetGroup <- consensus.groups[which(as.numeric(sapply(consensus.groups, length)) > 1)]

	for(group in subsetGroup )
	{
	    checkOrient <- consensusStrands[group,]
		checkOrientOp <- consensusStrandsOp[group,]
		rownames(checkOrientOp) <- paste( rownames(checkOrientOp),"_OPPOSITE", sep="")
 		checkOrient <- StrandStateMatrix(rbind(checkOrient, checkOrientOp))

	    sub.groups <- clusterContigs(checkOrient, 
			                        clusterParam=clusterParam, 
			                        recluster=cluster,
			                        clusterBy='homo', 
			                        randomise=TRUE, 
			                        similarityCutoff=0.9, 
			                        verbose=FALSE)  

	 	if(length(sub.groups) > 1)
	    {
	      	#If groups have been merged, make sure they are oriented in the same direction
			findGroupWithOp <- unlist(sapply(seq_along(sub.groups), function(x) if(length(grep("OPPOSITE", sub.groups[x])) > 0){return(names(sub.groups[x]))}))

			for(opposite in findGroupWithOp)
			{
				getGroup <- unlist(sub.groups[opposite])
				getElements <- getGroup[grep("OPPOSITE", getGroup)]
				getElements <- sapply(seq_along(getElements), function(x) strsplit(getElements[x], "_")[[1]][1] )
				reverseStrands <- unlist(linkedContigs[getElements])
		
				#Now flip the opposite elements for allStrands and orientationTable
				toReorientStrands <- switchAroo(allStrands[reverseStrands,])
				toReorientStrands <-  as.matrix(data.frame(lapply(data.frame(toReorientStrands), 
				 												   function(x){as.numeric(x)}))) 
				allStrands[reverseStrands,] <- toReorientStrands
				flipPlus <- reverseStrands[which(completeOrientation[reverseStrands,2] == "+")]
				flipMinus <- reverseStrands[which(completeOrientation[reverseStrands,2] == "-")]
				completeOrientation[flipPlus,2] <- '-'
				completeOrientation[flipMinus,2] <- '+'
			}
	    }
	}

	  mergeThoseGroups <- lapply(consensus.groups, function(x) as.character(melt(object[x])[,1]))
	  mergeThoseGroups <- mergeThoseGroups[order(sapply(mergeThoseGroups, length), decreasing=TRUE)]
	  mergeThoseGroups <- LinkageGroupList(
                          mergeThoseGroups, 
                          names= sapply(1:length(mergeThoseGroups), 
                          			  function(x)
                          			  	{
                          			  	paste('LG', x, ' (', length(mergeThoseGroups[[x]]), ')', sep='')
                          			  	}))

	   flipLengths <- sub('.*:', '', completeOrientation[,1])
	   completeOrientation <- ChrTable(GRanges(seqnames=sub(':.*', '', completeOrientation[,1]), 
	   										   IRanges(start=as.numeric(sub('-.*', '', flipLengths)), 
	   										   		   end=as.integer(sub('.*-', '', flipLengths))),
	   										   strand=as.character(completeOrientation[,2]),
	  								 		   name=as.character(completeOrientation[,1])))

	return(list(StrandStateMatrix(allStrands), completeOrientation, mergeThoseGroups))
}

####################################################################################################
#' reorientAndMergeLGs uses a simple dissimilarity to find misoriented fragments within linkage groups.
#' @param object List of vectors containing names of contigs belonging to each LG.
#' @param allStrands Table of type \code{strandStateMatrix} encompassing strand state for all contigs. Product of StrandSeqFreqTable.
#' @param cluster Number of times to recluster and take the consensus of. If NULL, clustering is 
#' run only once.
#' @param clusterParam optional \code{BiocParallelParam} specifying cluster to use for parallel execution.
#' When \code{NULL}, execution will be serial.
#' @param similarityCutoff merge contigs that are more similar this this
#' @param verbose Outputs information to the terminal. Default is TRUE.
#' @aliases reorientAndMergeLGs reorientAndMergeLGs-LinkageGroupList-StrandStateMatrix-method
#' @rdname reorientAndMergeLGs
#' @example inst/examples/reorientAndMergeLGs.R
#' @return a list consisting of a strandStateMatrix (a reoriented version of allStrands), a ChrTable 
#' containing contig names and orientations, as '+' or '-' and a merged LinkageGroupList.
#' 
#' @import GenomicRanges
#' @importFrom reshape2 melt 
#' @export
####################################################################################################

setMethod('reorientAndMergeLGs',
		  signature = signature(object='LinkageGroupList', 
		  					  allStrands = 'StrandStateMatrix'),
		  definition = reorientAndMergeLGs.func
		  )
