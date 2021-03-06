findSexGroups.func <- function(linkageGroupList, strandStateMatrix, callThreshold=0.2)
{

	createTable <- function(group, allStrands)												
	{	
		groupStrands <- allStrands[group,, drop=FALSE]

		tables <- sapply(1:ncol(groupStrands), 
					 function(y) sapply(c(paste(c(1,2,3), collapse="|"), 2), 
					 function(x) length(grep(x, groupStrands[,y]))))
		tabTots <- apply(tables, 1, sum)
		tabTots <- tabTots[2]/tabTots[1]

	}

	LGconsensus <- data.frame(do.call(rbind, lapply(linkageGroupList, createTable, strandStateMatrix))) 
	LGconsensus <- LGconsensus[order(LGconsensus),1, drop=FALSE]
	LGconsensus <- rownames(LGconsensus[which(LGconsensus[,1] <= callThreshold), ,drop=FALSE])

	if(length(LGconsensus) == 0)
	{
		warning('NO SEX LINKAGE GROUPS FOUND.')
		return(linkageGroupList)
	}else{
		LGNames <- replace(names(linkageGroupList), names(linkageGroupList) == LGconsensus, paste("sex", LGconsensus, sep='_')) 
		names(linkageGroupList) <- LGNames
		return(linkageGroupList)
	}
}

####################################################################################################
#' Bar plot all linkage groups, with the true chromosomes of contigs coloured.
#' @param linkageGroupList a LinkageGroupList object as generated by clusterContigs
#' @param strandStateMatrix a StrandStateMatrix object as generated by preprocessStrandTable 
#' @param callThreshold a value between 0 and 1 to threshold the proportion of libraries in which scaffolds are heterozygous (ie
#' unlikely to be a monosome sex chromosome). Default value is 0.2
#' 
#' @return the same LinkageGroupList entered into the function, but with potential sex chromosomes highlighted in the name attribute
#' 
#' @aliases findSexGroups 
#' @export
####################################################################################################

setMethod('findSexGroups',
          signature = signature(linkageGroupList = 'LinkageGroupList', strandStateMatrix='StrandStateMatrix'),
          definition = findSexGroups.func
)
