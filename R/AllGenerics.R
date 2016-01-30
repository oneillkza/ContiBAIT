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
		   		 snowCluster=NULL,
		   		 clusterBy=NULL,
		   		 verbose=NULL) standardGeneric("clusterContigs"),
		   signature='object')


## =========================================================================
## Generic for reorientLinkageGroups.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export reorientLinkageGroups
setGeneric("reorientLinkageGroups", 
		   function(object, 
		   			allStrands,
		   			previousOrient=NULL,
		   		 	verbose=NULL) standardGeneric("reorientLinkageGroups"),
		   signature=c('object', 'allStrands'))

## =========================================================================
## Generic for mergeLinkageGroups.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export mergeLinkageGroups
setGeneric("mergeLinkageGroups", 
           function(object, 
                    allStrands,
                    clusNum=NULL, 
                    cluster=NULL,
                    similarityCutoff=NULL) standardGeneric("mergeLinkageGroups"),
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
                    labels=NULL) standardGeneric("plotLGDistances"),
		   signature=c('object', 'allStrands'))

## =========================================================================
## Generic for barplotLinkageGroupCalls
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export barplotLinkageGroupCalls
setGeneric("barplotLinkageGroupCalls", 
           function(object, 
                    chrTable, 
                    by=NULL, 
                    returnTable=NULL) standardGeneric("barplotLinkageGroupCalls"),
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
## Generic for ideogramPlot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export ideogramPlot
setGeneric("ideogramPlot", 
           function(WatsonFreqList, 
                    CrickFreqList,
                    chrTable,
                    plotBy=NULL,
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
  
		  	cat('A strand state matrix for ', d[[1]], ' contigs over ',d[[2]],' libraries.\n')
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
		  	
		  	cat('A matrix of strand frequencies for ', d[[1]], ' contigs over ',d[[2]],' libraries.\n')
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
		  	
		  	cat('A matrix of read counts for ', d[[1]], ' contigs over ',d[[2]],' libraries.\n')
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
		  	cat('A linkage group list containing ', length(object), ' linkage groups.\n\n')
		  	show(data.frame(NumberOfContigs=head(sapply(object, length)), row.names=NULL))
		  	if(length(object) > 5)
            show(data.frame("...           "=tail(sapply(object, length)), row.names=seq(length(object)-5, length(object) )))

		  }
)

## show RawReadStrands
#' @name show,RawReadStrands-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a RawReadStrands
#' @return nothing
#' @description Shows a RawReadStrands
setMethod("show",
		  signature=signature(object="RawReadStrands"),
		  definition=function(object)
		  {
		  	cat('A data.frame containing ', nrow(object), ' reads.\n\n')
		  	#show(data.frame(NumberOfContigs=sapply(object, length)))
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

    		lg <- strsplit(as.character(object$LG), "\\.") # why not character() already?
    		len1 <- length(unique(sapply(lg, '[', 1)))     # first element of object$LG
    		len2 <- length(unique(sapply(lg, '[', 2)))   # second element of object$LG
 
		  	cat('A data.frame of', len1, 'LGs split into', len2, 'sub-groups from', nrow(object), 'ordered fragments.\n')
		  	if(length(unique(sapply(1:nrow(object), function(x) strsplit(as.character(object$LG), "\\.")[[x]][2]))) > 25)
		  	{
			  	show(head(table(object[,1])))
			  	cat('...')
			  	show(tail(table(object[,1])))
		  	}else{
		  		show(table(object[,1]))
		  	}
		  }
)

## show OrientationFrame
#' @name show,OrientationFrame-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a OrientationFrame
#' @return nothing
#' @description Shows a OrientationFrame
setMethod("show",
		  signature=signature(object="OrientationFrame"),
		  definition=function(object)
		  {
		  	elements <- nrow(object)
		  	misorientations <- nrow(object[which(object[,2] == '-'),])
		  	cat('A data.frame of ', elements, ' contigs with ',misorientations,' identified misorientations.\n')
		  }
)
