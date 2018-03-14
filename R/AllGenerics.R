## =========================================================================
## Generic for preprocessStrandTable.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export preprocessStrandTable
setGeneric("preprocessStrandTable", 
		   function(strandTable, 
		   			strandTableThreshold=NULL, 
		   			filterThreshold=NULL, 
		   			orderMethod=NULL, 
		   			lowQualThreshold=NULL, 
		   			verbose=NULL, 
		   			minLib=NULL, 
		   			ignoreInternalQual=NULL) standardGeneric("preprocessStrandTable"),
		   signature='strandTable')


## =========================================================================
## Generic for clusterContigs.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export clusterContigs
setGeneric("clusterContigs", 
		   	function(object, 
		   			 similarityCutoff=NULL,
			   		 recluster=NULL, 
			   		 minimumLibraryOverlap=NULL,
			   		 randomise=NULL,
			   		 randomSeed=NULL,
			   		 randomWeight=NULL,
			   		 clusterParam=NULL,
			   		 clusterBy=NULL,
			   		 verbose=NULL) standardGeneric("clusterContigs"),
			   		 signature='object')

## =========================================================================
## Generic for findSimilarLibraries.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export findSimilarLibraries
setGeneric("findSimilarLibraries", 
		   function(strandStateMatrix, 
		   		 	strandReadMatrix,
		   		 	chrGrange,
		   		 	chrNum,
			   		cluster=NULL, 
			   		clusterParam=NULL,
			   		verbose=NULL) standardGeneric("findSimilarLibraries"),
			   		signature=c('strandStateMatrix', 'strandReadMatrix', 'chrGrange', 'chrNum'))

## =========================================================================
## Generic for thoroughBed.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export thoroughBed
setGeneric("thoroughBed", 
		   function(bamFileList,
		   			relatedLibList, 
		   		 	qual=NULL,
		   		 	pairedEnd=NULL,
		   		 	rmdup=NULL,
		   			verbose=NULL) standardGeneric("thoroughBed"),
			   		signature=c('bamFileList', 'relatedLibList'))

## =========================================================================
## Generic for locateMisorients.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export locateMisorients
setGeneric("locateMisorients", 
		   function(compiledGrange, 
		   		 	gapFile=NULL,
		   		 	stateNum=NULL,
		   		 	readCutOff=NULL,
			   		verbose=NULL) standardGeneric("locateMisorients"),
			   		signature='compiledGrange')


## =========================================================================
## Generic for reorientAndMergeLGs.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export reorientAndMergeLGs
setGeneric("reorientAndMergeLGs", 
		   function(object, 
		   			allStrands,
			   		cluster=NULL, 
			   		clusterParam=NULL,
			   		similarityCutoff=NULL,
		   		 	verbose=NULL) standardGeneric("reorientAndMergeLGs"),
		   signature=c('object', 'allStrands'))

## =========================================================================
## Generic for orderAllLinkageGroups.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export orderAllLinkageGroups
setGeneric("orderAllLinkageGroups", 
           function(linkageGroupList,
            		strandStateMatrix, 
            		strandFreqMatrix, 
            		strandReadCount, 
            		whichLG=NULL, 
            		saveOrdered=NULL, 
            		orderCall=NULL, 
            		randomAttempts=NULL,
            		nProcesses = 1,
            		verbose=NULL) standardGeneric("orderAllLinkageGroups"),
		   signature=c('linkageGroupList', 'strandStateMatrix','strandFreqMatrix', 'strandReadCount' ))

## =========================================================================
## Generic for plotLGDistances.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export plotLGDistances
setGeneric("plotLGDistances", 
           function(object, 
                    allStrands,
                    lg=NULL,
                    labels=NULL,
                    state=NULL,
                    alreadyOrdered=NULL) standardGeneric("plotLGDistances"),
		   signature=c('object', 'allStrands'))

## =========================================================================
## Generic for barplotLinkageGroupCalls
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export barplotLinkageGroupCalls
setGeneric("barplotLinkageGroupCalls", 
           function(object, 
                    chrTable, 
                    by=NULL, 
                    bySize=NULL,
                    returnTable=NULL,
                    whichGroup=NULL,
                    percentage=NULL) standardGeneric("barplotLinkageGroupCalls"),
		   signature=c('object', 'chrTable'))

## =========================================================================
## Generic for plotWCdistribution
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export plotWCdistribution
setGeneric("plotWCdistribution", 
           function(object, 
                    allStrands,
                    filterThreshold=NULL) standardGeneric("plotWCdistribution"),
		   signature=c('object', 'allStrands'))

## =========================================================================
## Generic for makeBoxPlot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export makeBoxPlot
setGeneric("makeBoxPlot", 
           function(chrTable, 
                    linkage.contigs) standardGeneric("makeBoxPlot"),
		   signature=c('chrTable', 'linkage.contigs'))

## =========================================================================
## Generic for plotContigOrder
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export plotContigOrder
setGeneric("plotContigOrder", 
           function(contigOrder, 
                    lg=NULL,
                    verbose=NULL) standardGeneric("plotContigOrder"),
		   signature='contigOrder')

## =========================================================================
## Generic for ideogramPlot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export ideogramPlot
setGeneric("ideogramPlot", 
           function(WatsonFreqList, 
                    CrickFreqList,
                    chrTable,
                    plotBy=NULL,
                    sizeProportional=NULL,
                    showPage=NULL,
                    orderFrame=NULL,
   		            orientationData=NULL,
                    verbose=NULL) standardGeneric("ideogramPlot"),
		   signature=c('WatsonFreqList', 'CrickFreqList', 'chrTable'))

## =========================================================================
## Generic for writeBed
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export writeBed
setGeneric("writeBed",
		   function(chrTable, 
					orientationData, 
					contigOrder,
					libWeight=NULL,
					file=NULL) standardGeneric("writeBed"),
		   signature=c('chrTable', 'orientationData', 'contigOrder'))

## =========================================================================
## Generic for fixLinkageGroups
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export fixLinkageGroups
setGeneric("fixLinkageGroups",
		   function(contigOrdering, 
					orderFrame, 
					linkageGroupList,
					whichLG=NULL,
					relatedCutOff=NULL,
					verbose=NULL) standardGeneric("fixLinkageGroups"),
		   signature=c('contigOrdering', 'orderFrame', 'linkageGroupList'))

## =========================================================================
## Generic for flipOrderedLinkageGroups
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export flipOrderedLinkageGroups
setGeneric("flipOrderedLinkageGroups",
		   function(contigOrdering, 
					orderFrame, 
					linkageGroupList,
					strandStateMatrix,
					whichLG=NULL,
					maxiter=NULL,
					dissimilarityCutoff=NULL) standardGeneric("flipOrderedLinkageGroups"),
		   signature=c('contigOrdering', 'orderFrame', 'linkageGroupList', 'strandStateMatrix'))

## =========================================================================
## Generic for mergeFlankedLGs
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export mergeFlankedLGs
setGeneric("mergeFlankedLGs",
		   function(linkageGroupList, 
					strandStateMatrix, 
					buildConsensus=NULL,
					cluster=NULL,
					clusterParam=NULL,
					similarityCutoff=NULL,
					verbose=NULL) standardGeneric("mergeFlankedLGs"),
		   signature=c('linkageGroupList', 'strandStateMatrix'))


## =========================================================================
## Generic for findSexGroups
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export findSexGroups
setGeneric("findSexGroups",
		   function(linkageGroupList, 
					strandStateMatrix, 
					callThreshold=0.2) standardGeneric("findSexGroups"),
		   signature=c('linkageGroupList', 'strandStateMatrix'))



## =========================================================================
## show Methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## show StrandStateMatrix
#' @name show,StrandStateMatrix-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a StrandStateMatrix
#' @return nothing
#' @description Shows a StrandStateMatrix
setMethod("show",
		  signature=signature(object="StrandStateMatrix"),
		  definition=function(object)
		  {
        	d <- dim(object)
  
		  	cat('A strand state matrix for', d[[1]], 'contigs over',d[[2]],'libraries.\n')
		  }
)


## show StrandFreqMatrix
#' @name show,StrandFreqMatrix-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a StrandFreqMatrix
#' @return nothing
#' @description Shows a StrandFreqMatrix
setMethod("show",
		  signature=signature(object="StrandFreqMatrix"),
		  definition=function(object)
		  {
		  	d <- dim(object)
		  	
		  	cat('A matrix of strand frequencies for', d[[1]], 'contigs over',d[[2]],'libraries.\n')
		  }
)

## show StrandReadMatrix
#' @name show,StrandReadMatrix-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a StrandReadMatrix
#' @return nothing
#' @description Shows a StrandReadMatrix
setMethod("show",
		  signature=signature(object="StrandReadMatrix"),
		  definition=function(object)
		  {
		  	d <- dim(object)
		  	
		  	cat('A matrix of read counts for', d[[1]], 'contigs over',d[[2]],'libraries.\n')
		  }
)


## show LinkageGroupList
#' @name show,LinkageGroupList-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a LinkageGrouplist
#' @return nothing
#' @description Shows a LinkageGroupList
setMethod("show",
		  signature=signature(object="LinkageGroupList"),
		  definition=function(object)
		  {
		  	cat('A linkage group list containing', length(object), 'linkage groups.\n\n')
		  	if(length(object) < 10)	{	  	
		  	show(data.frame(NumberOfContigs=sapply(object, length), row.names=NULL))
		  	}else{
 		  	show(data.frame(NumberOfContigs=head(sapply(object, length)), row.names=NULL))
            show(data.frame("...           "=tail(sapply(object, length)), 
            				row.names=seq(length(object)-5, length(object) )))
        	}
		  }
)

## show LibraryGroupList
#' @name show,LibraryGroupList-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a LibraryGroupList
#' @return nothing
#' @description Shows a LibraryGroupList
setMethod("show",
		  signature=signature(object="LibraryGroupList"),
		  definition=function(object)
		  {
		  	cat('A library group list containing', length(object), 'linkage groups.\n\n')
		  	if(length(object) < 10)	{	  	
		  	show(data.frame(Contig=names(object),
		  					MostlyCrickLibraries=sapply(seq_along(object), function(x) length(object[[x]][[1]])),
		  					MostlyWatsonLibraries=sapply(seq_along(object), function(x) length(object[[x]][[2]])), 
		  					row.names=NULL))
		  	}else{
		  	show(data.frame(Contig=c(head(names(object)),
		  							"...",
		  							tail(names(object))),
		  					MostlyCrickLibraries=c(head(sapply(seq_along(object), function(x) length(object[[x]][[1]]))), 
		  										 "...", 
		  										 tail(sapply(seq_along(object), function(x) length(object[[x]][[1]])))),
		  					MostlyWatsonLibraries=c(head(sapply(seq_along(object), function(x) length(object[[x]][[2]]))),
		  										 "...", 
		  										 tail(sapply(seq_along(object), function(x) length(object[[x]][[2]])))),
		  					row.names=c(1:6, "...", seq(length(object)-5, length(object) ) )))
		  	}
		  }
)

## show ContigOrdering
#' @name show,ContigOrdering-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a ContigOrdering
#' @return nothing
#' @description Shows a ContigOrdering
setMethod("show",
		  signature=signature(object="ContigOrdering"),
		  definition=function(object)
		  {

    		lg <- strsplit(as.character(object[,1]), "\\.") # why not character() already?
    		len1 <- length(unique(sapply(lg, '[', 1)))     # first element of object$LG
    		len2 <- length(unique(sapply(lg, '[', 2)))   # second element of object$LG
 
		  	cat('A matrix of', len1, 'LGs split into', len2, 
		  		'sub-groups from', nrow(object), 'ordered fragments.\n')
		  	if(length(table(object[,1])) > 10)
		  	{
			  	show(head(table(object[,1])))
			  	cat('...')
			  	show(tail(table(object[,1])))
		  	}else{
		  		show(table(object[,1]))
		  	}
		  }
)

## show StrandStateList
#' @name show,StrandStateList-method
#' @export
#' @importFrom utils head tail
#' @docType methods
#' @title show-methods
#' @param object a StrandStateList
#' @return nothing
#' @description Shows a StrandStateList
setMethod("show",
		  signature=signature(object="StrandStateList"),
		  definition=function(object)
		  {
	  		cat('A strandStateList containing', length(object), 'StrandStateMatrices.\n\n')
		  	if(length(object) < 10)	{	  	
		  		show(data.frame(Name=names(object), NumberOfContigs= sapply(1:length(object), function(x) nrow(object[[x]])), row.names=NULL))
		  	}else{
				endLength <- seq(length(object)-4, length(object))
	 		  	show(data.frame(Name=head(names(object),5), NumberOfContigs= sapply(1:5, function(x) nrow(object[[x]])), row.names=NULL))
	            show(data.frame("...           "=tail(names(object),5), 
	            				NumberOfContigs=sapply(endLength, function(x) nrow(object[[x]])), row.names=endLength ))
		  	}
		  }
)