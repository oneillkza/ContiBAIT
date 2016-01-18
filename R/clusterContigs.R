clusterContigs.func <- function(object, #heatFile from contiBAIT; a data frame containing strand inheritance for every contig
						   similarityCutoff=0.7, #
						   recluster=NULL, #whether to 
						   minimumLibraryOverlap=5,
                           randomise=TRUE,
						   randomSeed=NULL,
						   randomWeight=NULL,
						   snowCluster=NULL,
						   clusterBy='hetero',
						   verbose=TRUE)
{
	
	#similarityCutoff <- 0.7
	if(!is.null(randomSeed))
		set.seed(randomSeed)
	
	
	runOnce <- function(object, randomise, randomWeight)
	{
		randomOrder <- rownames(object) # by default, no reordering
		names(randomOrder) <- rownames(object)
		#Randomise object contig order:
		if(randomise)
			if(!is.null(randomWeight))
			{
				#dequantise by adding some tiny random noise in case there are many zeroes:
				randomWeight <- randomWeight + runif(length(randomWeight)) / (10^5*max(randomWeight))
				names(randomOrder) <- sample(rownames(object), prob=randomWeight)
				
			}else
			{
				names(randomOrder) <- sample(rownames(object))
			}
		#Use as data.frame in case only one column is present
		object <- as.data.frame(object[names(randomOrder),])
		
	
		#Set up a list for storing linkage groups (clusters):
		linkageGroups <- list()
		linkageStrands <- object[1,]
		
		contigs <- rownames(object)
		linkageGroups[[1]] <- c(1)
		
		#Function to assign a new contig to an existing linkage group, or create a new one:
		if(verbose){message(paste('Inicializing contig ', contigs[1], ' [1/', nrow(object), '] as LG1', sep=""))}		

		for (contig.num in 2:nrow(object))
		{
			if(verbose){message('Clustering contig ', contigs[contig.num], ' [', contig.num, '/', nrow(object), ']   \r', appendLF=FALSE )}
			computePairwiseSim <- function(linkage.num, contig.num)
			{
				contigStrand <- object[contig.num,]
				linkageStrand <- linkageStrands[linkage.num,]
				#numCommon <- length(which(!is.na(linkageStrand)&!is.na(contigStrand)))
				#if(numCommon < minimumLibraryOverlap)
				#	return(NA)
				#suppressWarnings(1 - daisy(rbind(contigStrand, linkageStrand) )[1] )
				contigStrand[which(contigStrand==3)] <- 2
				linkageStrand[which(linkageStrand==3)] <- 2	
				computeSim(contigStrand, linkageStrand, minimumLibraryOverlap)
			}
			similarities <- sapply(1:nrow(linkageStrands), computePairwiseSim, contig.num)
			best.match <- which.max(similarities)
			#If no good match, make this contig the founder of a new linkage group:
			if (length(best.match) ==0 || similarities[best.match] < similarityCutoff)
			{
				linkageGroups[[length(linkageGroups)+1]] <- c(contig.num)
				linkageStrands <- rbind(linkageStrands, object[contig.num,])
			}else
				#Otherwise, add this to the best matched group, and recompute the strand state for that group:
			{
				if(verbose){message(paste('\n  -> Adding ', contigs[contig.num],' to LG', best.match, ' for a cluster of ',length(linkageGroups[[best.match]])+1 , sep=""))}
				linkageGroups[[best.match]] <- append(linkageGroups[[best.match]], contig.num)
				strandVec <- computeConsensus(linkageGroups[[best.match]], object)
				linkageStrands[best.match,] <- strandVec
			}
		}
		linkageGroups <- lapply(linkageGroups, function(lg){names(randomOrder)[lg]})

		return(new('LinkageGroupList',linkageGroups))
	}	
	
        
	linkageGroups <- list()
    
    #If no reclustering, just run once:
    #Code starts here!
    if(clusterBy == 'homo')
    {
		object <- replace(object, object == 2, NA)
	}else if (clusterBy == 'hetero'){
		object <- replace(object, object == 3, 1)
	}else{
		warning('### Unrecognized clusterBy parameter! ###')
		break
	}
	

	if(is.null(recluster))
	{
		linkageGroups <- runOnce(object, randomise, randomWeight)
	}else
	{
		runOneForTheEnsemble <- function(dump, object, randomise, randomWeight)
		{
			clusterVec <- rep(0, nrow(object))
			names(clusterVec) <- rownames(object)
			linkageGroups <- runOnce(object, randomise, randomWeight)
			#linkageGroups <- mergeLinkageGroups(linkageGroups, object)
			for (lg in 1:length(linkageGroups))
				clusterVec[linkageGroups[[lg]]] <- lg
			clusterVec
			#linkageGroups
		}
		
		#Then get consensus:
		if(!is.null(snowCluster))
		{
			if(verbose){message(paste('-> Running ', recluster, ' clusterings in parallel on ', length(snowCluster), ' processors', sep=""))}
			multiClust <- parLapply(snowCluster, 1:recluster, runOneForTheEnsemble, object, randomise, randomWeight)
		}else
		{
			multiClust <- lapply(1:recluster, runOneForTheEnsemble, object, randomise, randomWeight)
		}
		
		
		multiClue <- lapply(multiClust, as.cl_partition)
		multiEnsemble <- cl_ensemble(list=multiClue)
		multiConsensus <- cl_consensus(multiEnsemble)
		multiLabels <- apply(multiConsensus$.Data, 1, which.max) #take most likely cluster for each contig
		
		#Roll back into a linkageGroups list:
		linkageGroups <- list()
		for (label in unique(multiLabels))
		{
			lgContigs <- names(multiLabels)[which(multiLabels==label)]
			if(length(lgContigs) > 0)
				linkageGroups[[as.character(label)]] <- lgContigs
		}
	}

	#order linkage groups by biggest first
	#names(linkageGroups) <- sapply(1:length(linkageGroups), function(x){paste('LG', x, ' (', length(linkageGroups[[x]]), ')', sep='') })
  	linkageGroups <- new('LinkageGroupList',linkageGroups[order(sapply(linkageGroups, length), decreasing=TRUE)], names= sapply(1:length(linkageGroups), function(x){paste('LG', x, ' (', length(linkageGroups[[x]]), ')', sep='') }))


	return(linkageGroups)
}

####################################################################################################
#' clusterContigs -- agglomeratively clusters contigs into linkage groups based on strand inheritance
#' 
#' @param object \code{data.frame} containing strand inheritance information for every contig (rows)
#' in every library (columns). This should be the product of strandSeqFreqTable
#' @param similarityCutoff place contigs in a cluster when their strand state is at least this similar
#' @param recluster =NULL Number of times to recluster and take the consensus of. If NULL, clustering is 
#' run only once.
#' @param minimumLibraryOverlap for two contigs to be clustered together, the strand inheritance must 
#' be present for both contigs in at least this many libraries (in addition to their similarity being at least 
#' similarityCutoff)
#' @param randomise whether to reorder contigs before clustering
#' @param randomSeed random seed to initialize clustering
#' @param randomWeight vector of weights for contigs for resampling. If NULL, uniform resampling is used.
#' @param clusterBy Method for performing clustering. Default is 'hetero' (for comparing heterozygous calls to homozygous). Alternative is 'homo' (for compairson between the two homozygous calls)
#' Typically this should be a measure of contig quality, such as library coverage, so that clustering tends to
#' start from the better quality contigs.
#' @param snowCluster optional snowCluster for parallel execution
#' @param verbose = TRUE prints function progress
#' @details Note that a more stringent similarity cutoff will result in more clusters, and a longer run time,
#' since at every iteration a distance is computed to the existing clusters. However, in lower-quality data, a
#' more stringent cutoff may be necessary to reduce the number of contigs that are erroneously grouped.
#' @return \code{LinkageGroupList} of vectors containing labels of contigs belonging to each linkage
#' group
#' 
#' @aliases clusterContigs clusterContigs,StrandStateMatrix,StrandStateMatrix-method
#' 
#' @example inst/examples/clusterContigs.R
#' 
#' @importFrom cluster daisy
#' @import clue
#' @export
#' @include AllClasses.R
#
####################################################################################################

setMethod('clusterContigs',
		  signature = signature(object='StrandStateMatrix'),
		  definition = clusterContigs.func
		  )
