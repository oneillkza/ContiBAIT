clusterContigs.func <- function(object, #heatFile from contiBAIT; a data frame containing strand inheritance for every contig
						   similarityCutoff=0.7, #
						   recluster=NULL, #whether to 
						   minimumLibraryOverlap=5,
                           randomise=TRUE,
						   randomSeed=NULL,
						   randomWeight=NULL,
						   clusterParam=NULL,
						   clusterBy='hetero',
						   verbose=TRUE)
{	
	if(!is.null(randomSeed))
		set.seed(randomSeed)
	
	
	runOnce <- function(object, randomise, randomWeight, similarityCutoff)
	{
		#Randomise object contig order:
		if(randomise) 
		{			
			if(!is.null(randomWeight))
			{
				#dequantise by adding some tiny random noise in case there are many zeroes:				
				randomWeight <- randomWeight + runif(length(randomWeight)) / (10^5*max(randomWeight))
			}

			object <- object[sample(1:nrow(object), prob=randomWeight),]
		}
		
		#Use as data.frame in case only one column is present
		object <- as.data.frame(object)

		linkageGroups <- .Call('buildLinkageGroups', data.matrix(object), similarityCutoff, minimumLibraryOverlap, 
			verbose, if(verbose) rownames(object) else vector("character", 0))
			
		linkageGroups <- lapply(linkageGroups, function(lg){rownames(object)[lg]})
		return(linkageGroups)
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
		stop('### Unrecognized clusterBy parameter! ###')
	}
	

	if(is.null(recluster))
	{
		linkageGroups <- runOnce(object, randomise, randomWeight, similarityCutoff)
	}else
	{
		runOneForTheEnsemble <- function(dump, object, randomise, randomWeight, similarityCutoff)
		{
			clusterVec <- rep(0, nrow(object))
			names(clusterVec) <- rownames(object)
			linkageGroups <- runOnce(object, randomise, randomWeight, similarityCutoff)
			for (lg in seq_len(length(linkageGroups)))
				clusterVec[linkageGroups[[lg]]] <- lg
				clusterVec
		}
		
		#Then get consensus:
		if(!is.null(clusterParam))
		{
			if(verbose){message(paste('-> Running ', recluster, 
									  ' clusterings in parallel on ', 
									  clusterParam$workers, ' processors', sep=""))}
			#Prevent spew of multiple core verbose messages by turning verbose off
			verbose=FALSE
			multiClust <-bplapply(  seq_len(recluster), 
						  			runOneForTheEnsemble, 
									object, 
									randomise, 
									randomWeight, 
									similarityCutoff,
									BPPARAM=clusterParam)
		}else
		{
			multiClust <- lapply(seq_len(recluster), 
								 runOneForTheEnsemble, 
								 object, 
								 randomise, 
								 randomWeight, 
								 similarityCutoff)
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

	linkageGroups <- linkageGroups[order(sapply(linkageGroups, length), decreasing=TRUE)]
	#order linkage groups by biggest first
	linkageGroups <- LinkageGroupList(
  						 linkageGroups, 
  						  names= sapply(1:length(linkageGroups), 
  						  			  function(x){
  						  			  	paste('LG', x, ' (', length(linkageGroups[[x]]), ')', sep='') 
  						  			  	}))


	return(linkageGroups)
}

####################################################################################################
#' clusterContigs -- agglomeratively clusters contigs into linkage groups based on strand inheritance
#' 
#' @param object \code{data.frame} containing strand inheritance information for every contig (rows)
#' in every library (columns). This should be the product of strandSeqFreqTable
#' @param similarityCutoff place contigs in a cluster when their strand state is at least this similar
#' @param recluster Number of times to recluster and take the consensus of. If NULL, clustering is 
#' run only once.
#' @param minimumLibraryOverlap for two contigs to be clustered together, the strand inheritance must 
#' be present for both contigs in at least this many libraries (in addition to their similarity being at least 
#' similarityCutoff)
#' @param randomise whether to reorder contigs before clustering
#' @param randomSeed random seed to initialize clustering
#' @param randomWeight vector of weights for contigs for resampling. If NULL, uniform resampling is used.
#' Typically this should be a measure of contig quality, such as library coverage, so that clustering tends to
#' start from the better quality contigs.
#' @param clusterBy Method for performing clustering. Default is 'hetero' (for comparing heterozygous calls to homozygous). 
#' Alternative is 'homo' (for compairson between the two homozygous calls)
#' @param clusterParam optional \code{BiocParallelParam} specifying cluster to use for parallel execution.
#' When \code{NULL}, execution will be serial.
#' @param verbose prints function progress
#' @details Note that a more stringent similarity cutoff will result in more clusters, and a longer run time,
#' since at every iteration a distance is computed to the existing clusters. However, in lower-quality data, a
#' more stringent cutoff may be necessary to reduce the number of contigs that are erroneously grouped.
#' @return \code{LinkageGroupList} of vectors containing labels of contigs belonging to each linkage
#' group
#' 
#' @details Note that \code{clusterParam} requires \code{BiocParallel} to be installed.
#' 
#' @aliases clusterContigs clusterContigs,StrandStateMatrix,StrandStateMatrix-method
#' 
#' @example inst/examples/clusterContigs.R
#' 
#' @importFrom cluster daisy
#' @import clue
#' @import BiocParallel
#' @export
#' @include AllClasses.R
#
####################################################################################################

setMethod('clusterContigs',
		  signature = signature(object='StrandStateMatrix'),
		  definition = clusterContigs.func
		  )
